// Übergabe von AT und BPXT zusätzlich für adjointSolver -> Operator statt OperatorWithTranspose?

#include "PDP.hh"

#include <Spacy/Spaces/ProductSpace.h>
#include <Spacy/ZeroVectorCreator.h>
#include <Spacy/Util/Mixins/Get.h>
#include <Spacy/Util/Logger.h>
#include <Spacy/Util/Log.h>
#include <iostream>

#define LAPACK_COL_MAJOR 102

namespace Spacy
{
    namespace
    {
        Logger< unsigned > logIterationsState( "stateIt.log" );
        Logger< unsigned > logIterationsAdjoint( "adjointIt.log" );
        Logger< unsigned > logIterationsIcg( "icgIt.log" );
        Logger< double > logTheta( "Theta.log" );
    }
    namespace PDP
    {
        
        DEFINE_LOG_TAG(static const char* log_tag = "PDP");
        
        
        Solver::Solver(Operator M, OperatorWithTranspose A, OperatorWithTranspose AT, OperatorWithTranspose minusB, CallableOperator BPX, CallableOperator BPXT, CallableOperator controlSolver, ::Spacy::CG::Regularization R,const VectorSpace& totalSpace):
        M_((M)), A_((A)), AT_((AT)), minusB_((minusB)), BPX_(BPX), BPXT_(BPXT), R_(R),
        totalSpace_(totalSpace), primalSpace_(M_.domain()),adjointSpace_(minusB_.range()), controlSpace_(minusB_.domain()), stateSpace_(A_.domain()),
        ChebPrec_(A,BPX), controlSolver_(std::move(controlSolver))
        {
        }
        
        Vector Solver::C_(const Vector& input) const
        {
            auto result = A_(stateSpace_.project(input));
            result += minusB_(controlSpace_.project(input));
            return result;
            
        }
        
        Vector Solver::C_transpose(const Vector& input) const
        {
            auto result = zero(primalSpace_);
            result += primalSpace_.embed( A_.transposed(input) );
            
            
            result += primalSpace_.embed( minusB_.transposed(input) );
            return result; 
        }
        
        void Solver::setInternalAccuracies(double stateAcc, double adjointAcc, double modifiedCGAcc, double chebyshevAcc)
        {
            stateAcc_=stateAcc;
            adjointAcc_=adjointAcc;
            modifiedCGAcc_=modifiedCGAcc;
            chebyshevAcc_=chebyshevAcc;
        }
        
        
        Vector Solver::operator()(const Vector &b) const
        {
            std::cout << std::endl;
            ::Spacy::Vector x = zero(primalSpace_);
            ::Spacy::Vector result = zero(b.space());
            theta_ = 0.0;
            convex_flag_ = false;
            bool logConvex = true;
            result_ = Result::Failed;
            definiteness_ = DefiniteNess::PositiveDefinite;
            
            while(!convex_flag_)
            {
                result = PDPLoop(x,b);
                convex_flag_ = isPositiveDefinite();
                
                LOG(log_tag, "PPCG is Convex: ", convex_flag_);
                LOG(log_tag, "Energy in PPCG is Convex: ", signalConvex_(false,false));
                
                if(!signalConvex_(false,false))
                {
                    return result;
                }
                
                if(theta_ == 0.0)
                    theta_ = 1e-8;
                else
                    theta_ *= 100.0;
                
                if(!convex_flag_)
                    logConvex = false;
            }
            
            if(!logConvex)
                definiteness_ = DefiniteNess::Indefinite;
            
            return result;
        }
        bool Solver::isPositiveDefinite() const
        {
            return definiteness_ == DefiniteNess::PositiveDefinite;
        }
        
        Vector Solver::PDPLoop(Vector x,const Vector& b) const
        {
            LOG_SEPARATOR(log_tag);
            
            LOG(log_tag, "Regularization Parameter: ", theta_);
            LOG(log_tag, "Relative Accuracy: ", getRelativeAccuracy());
            LOG(log_tag, "Eps: ", eps());
            
            cheb_created_ = false;
            convex_flag_ = true;
            energyIsConvex_ = true;
            definiteness_ = DefiniteNess::PositiveDefinite;
            normDxOld_ = 0;
            
            Real stateAcc = stateAcc_;
            Real adjointAcc = adjointAcc_;
            Real modifiedCGAcc = modifiedCGAcc_;
            
            auto absAccIncrease = 1.0;
            
            /// Initialize Adjoint Solver
            Spacy::CG::NoRegularization NR;
            auto adjointSolver = Spacy::CG::Solver(AT_, BPXT_, NR, true);
            adjointSolver.setEps(eps());
            adjointSolver.setAbsoluteAccuracy(eps());
            adjointSolver.setMaxSteps(getMaxSteps());
            adjointSolver.setTerminationCriterion( CG::Termination::AdaptiveRelativeEnergyError());
            
            /// Initialize Primal Solver
            auto stateSolver = Spacy::CG::Solver(A_, BPX_, NR, true);
            stateSolver.setEps(eps());
            stateSolver.setAbsoluteAccuracy(eps());
            stateSolver.setMaxSteps(getMaxSteps());
            stateSolver.setTerminationCriterion( CG::Termination::AdaptiveRelativeEnergyError());
            
            std::vector<Spacy::Real> alpha_vec;
            std::vector<Spacy::Real> beta_vec;
            adjointSolver.callback = [&alpha_vec,&beta_vec]( const std::vector<Spacy::Real>& vec1, const std::vector<Spacy::Real>& vec2)
            {
                alpha_vec = vec1;
                beta_vec = vec2;
            };
            
            
            stateSolver.callback = [&alpha_vec,&beta_vec]( const std::vector<Spacy::Real>& vec1, const std::vector<Spacy::Real>& vec2)
            {
                alpha_vec = vec1;
                beta_vec = vec2;
            };
            
            
            
            auto condBPXA = 0.0;
            auto condQH = 0.0;
            
            sumNormSquaredUpdate_ = 0.0;
            
            
            auto bx = primalSpace_.project(b);
            
            /// x==0 
            /// !!!!!!!!!!! Does not work if theta != 0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            auto rx = primalSpace_.project( - bx);
            
            
            //if(theta_ != 0) rx+= theta_*R_(x_);
            
            
            Real omega = 0.0;
            
            double eigMin = std::numeric_limits<double>::max();
            double eigMax = std::numeric_limits<double>::min();
            
            auto dx=zero(primalSpace_);
            auto Mdx=zero(primalSpace_);
            auto p = zero(adjointSpace_);
            auto dp = zero(adjointSpace_);
            auto rp = zero(adjointSpace_);
            auto bp = zero(adjointSpace_);
            auto ddy=zero(stateSpace_);
            
            
//             ChebPrec_.setRelativeAccuracy(chebyshevAcc_);
            
            Spacy::ModifiedPPCG::TriangularConstraintPreconditioner triH(BPX_,controlSolver_,BPX_,minusB_,totalSpace_,stateSpace_);
            Spacy::ModifiedPPCG::Solver modPPCG_(triH, M_, minusB_, totalSpace_, stateSpace_);
            
            
            modPPCG_.setTerminationCriterion(  CG::Termination::AdaptiveRelativeEnergyError() );
            
            adjointAcc = std::max(1e-12,Spacy::Mixin::get(adjointAcc));
            stateAcc = std::max(1e-12,Spacy::Mixin::get(stateAcc));
            
            LOG(log_tag, "AdjointAcc: ", adjointAcc);
            LOG(log_tag, "StateAcc: ", stateAcc);
            LOG(log_tag, "ModifiedCGAcc: ", modifiedCGAcc);
            
            if(cheb_created_ == false)
            {
                LOG(log_tag, "ChebyshevAcc: ", chebyshevAcc_);
            }
            
            adjointSolver.setRelativeAccuracy(adjointAcc);
            stateSolver.setRelativeAccuracy(stateAcc);
            modPPCG_.setRelativeAccuracy(modifiedCGAcc);
            modPPCG_.setMaxSteps(5000);
            
            for(unsigned int step = 1; step < getMaxSteps(); step++)
            {
                LOG_SEPARATOR(log_tag);
                // Check regularization part
                dp = adjointSpace_.embed(adjointSolver(stateSpace_.project(-rx),adjointSpace_)); //Vermerk: b ^= -c in BA

//                 Spacy::Preconditioner::computeEigsFromCGCoefficients(alpha_vec, beta_vec);

                logIterationsAdjoint(adjointSolver.getIterations());
                LOG(log_tag, "Adjoint Solver Number of Iterations: ", adjointSolver.getIterations());
                dual_iterations_ += adjointSolver.getIterations();
                if(adjointSolver.indefiniteOperator())
                {
                    signalConvex_(true,false);
                    return x;
                }

                p += dp;
                
                
                rx += C_transpose(dp);
                
                triH.setSpectralBounds(eigMin, eigMax);
                dx = primalSpace_.project(modPPCG_(totalSpace_.embed(-rx)));
                
                LOG(log_tag, "PPCG Number of Iterations: ", modPPCG_.getIterations());
                

                convex_flag_ = modPPCG_.isPositiveDefinite();
                
                //         if(modPPCG_->isSingular())
                //         {
                //             LOG_INFO(log_tag, "inexactCG is singular");
                //             signalConvex_(true,false);
                //             return x;
                //         }
                
                if(!convex_flag_)
                {
                    definiteness_ = DefiniteNess::Indefinite;
                    result_ = Result::TruncatedAtNonConvexity;
                    return x;
                }
               
                bp = -( rp + C_(dx) );
               
                auto ddy = stateSolver(bp, stateSpace_ ); 
              
                dx += primalSpace_.embed(ddy);
               
                LOG(log_tag, "State Solver Number of Iterations: ", stateSolver.getIterations());
               
                proj_iterations_ += stateSolver.getIterations();
                
//                 Spacy::Preconditioner::computeEigsFromCGCoefficients(alpha_vec, beta_vec);
                
                if(stateSolver.indefiniteOperator())
                {
                    signalConvex_(true,false);
                    return x;
                }
                Mdx=M_(dx);//+theta_*R_(dx);
                
                Real dxMdx = Mdx(dx);
                
                omega = -rx(dx)/dxMdx;
                LOG(log_tag, "omega: ", omega);
                
                
                auto normX = sqrt(M_(x)(x));
                normDx_ =std::fabs(get(omega)) *sqrt(dxMdx);
                
                //LOG(log_tag, "normX (old): ", normX, "normDx_: ", normDx_);
                
                energyIsConvex_ = signalConvex_(false,false);
                bool converged = convergenceTest(convex_flag_, energyIsConvex_, M_, x ,  omega*dx);
                
                if(converged)// && !projected_)
                {
                    
                    LOG(log_tag, "Number of Iterations  PDP: ", step);
                    LOG(log_tag, "Number of Dual Iterations: ", dual_iterations_);
                    LOG(log_tag, "Number of PPCG Iterations: ", modPPCG_.getIterationsInLifeTime());
                    LOG(log_tag, "Number of Proj Iterations: ", proj_iterations_);
                    LOG(log_tag, "Number of MG calls in  CG: ", dual_iterations_ + proj_iterations_);
                    LOG(log_tag, "Number of MG calls   PPCG: ", triH.getBPXCallsInLifeTime());
                    LOG(log_tag, "Number of MG calls  total: ", triH.getBPXCallsInLifeTime()+dual_iterations_ + proj_iterations_);
                    LOG(log_tag, "Average no. MG calls /  A: ", triH.getBPXCallsInLifeTime()/(2.0*modPPCG_.getIterationsInLifeTime()));
                    
                    result_ = Result::Converged;
                    iterations_ = step;
                    std::cout << std::endl;
                    return x;
                }
                
                x +=get(omega)*dx;
                
                rx += omega * Mdx;
                rp += omega * C_(dx);
                
                
                normDxOld_ = normDx_;
            }
            LOG_INFO(log_tag, "pcgStep > getMaxSteps()");
            result_ = Result::NotConverged;
            iterations_ = getMaxSteps();
            return x;
        }
        
        bool Solver::convergenceTest(bool isConvex, bool energyIsConvex, CallableOperator H, const Vector &x, const Vector &dx) const
        {
            
            if(!isConvex)
            {
                LOG_INFO(log_tag, "Not Converged: !isConvex");
                return false;
            }
            if(!energyIsConvex)
            {
                LOG_INFO(log_tag, "Not Converged: !energyIsConvex");
                return false;
            }
            
            
            sumNormSquaredUpdate_ += get(normDx_*normDx_);
            
            if(normDxOld_ == 0) //normDxOld_ = 0 -> noch im ersten Schritt, erzwingt so mind. 2 Schritte (es ei denn dx ist sehr sehr klein)
                logTheta(0.0);
            
            else
                logTheta(get(normDx_/normDxOld_));
            
            if(normDx_ < eps())
            {
                LOG_INFO(log_tag, "Converged: normDx_ < eps()");
                return true;
            }
            
            if(normDxOld_ == 0) //normDxOld_ = 0 -> noch im ersten Schritt, erzwingt so mind. 2 Schritte (es ei denn dx ist sehr sehr klein)
            {
                // LOG_INFO(log_tag, "Not Converged: normDxOld = 0");
                return false;
            }
            else
            {
                auto normX = sqrt(M_(x)(x));
                
                auto contraction = normDx_/normDxOld_; //wenn contraction > desiredContraction -> erhöhe lookAhead
                if(contraction >= 1)
                {
                    LOG_INFO(log_tag, "Not Converged: contraction >= 1");
                    return false;
                }
                
                
                
                auto normSolutionEstimate = sqrt(sumNormSquaredUpdate_);
                auto normErrorEstimate = (contraction/sqrt(1-contraction*contraction))*normDx_;
                
                if(normErrorEstimate < getRelativeAccuracy() * normSolutionEstimate)
                {
                    std::cout << " Convergence Test: "  << normErrorEstimate << " < " << getRelativeAccuracy() << " * " << normSolutionEstimate << " theta " << contraction << std::endl;
                    std::cout << " Exceeded desired accuracy by factor: " << normErrorEstimate/(getRelativeAccuracy() * normSolutionEstimate) << std::endl;
                    return true;
                }
                else
                {
                std::cout << " Convergence Test: "  << normErrorEstimate << " > " << getRelativeAccuracy() << " * " << normSolutionEstimate << " theta " << contraction << std::endl;
                }
            }
            return false;
        }//*/
    }
}
