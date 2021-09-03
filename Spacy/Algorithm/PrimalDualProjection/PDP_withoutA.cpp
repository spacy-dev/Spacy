#include "PDP_withoutA.hh"

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
//         Logger< unsigned > logIterationsIcg( "icgIt.log" );
        Logger< double > logTheta( "Theta.log" );
    }
    namespace PDP_withoutA
    {
        
        DEFINE_LOG_TAG(static const char* log_tag = "PDP");
        
        
        Solver::Solver(Operator M, Spacy::SolverBase stateSolver, Spacy::SolverBase adjointSolver, OperatorWithTranspose minusB, CallableOperator surrogateStateSolver, CallableOperator surrogateAdjointSolver, CallableOperator controlSolver, const VectorSpace& totalSpace):
        M_((M)), stateSolver_(std::move(stateSolver)), adjointSolver_(std::move(adjointSolver)), minusB_((minusB)), surrogateStateSolver_(std::move(surrogateStateSolver)), surrogateAdjointSolver_(std::move(surrogateAdjointSolver)), totalSpace_(totalSpace), primalSpace_(M_.domain()),adjointSpace_(minusB_.range()), controlSpace_(minusB_.domain()), stateSpace_(stateSolver_.domain()), controlSolver_(std::move(controlSolver))
        {
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
            
            convex_flag_ = true;
            energyIsConvex_ = true;
            definiteness_ = DefiniteNess::PositiveDefinite;
            normDxOld_ = 0;
            
            Real stateAcc = stateAcc_;
            Real adjointAcc = adjointAcc_;
            Real modifiedCGAcc = modifiedCGAcc_;
            
            sumNormSquaredUpdate_ = 0.0;
            

            auto rx = primalSpace_.project(-b);
            
            
            Real omega = 0.0;
            
            
            auto dx=zero(primalSpace_);
            auto Mdx=zero(primalSpace_);
            auto p = zero(adjointSpace_);
            auto dp = zero(adjointSpace_);
            auto rp = zero(adjointSpace_);
            auto bp = zero(adjointSpace_);
            auto ddy=zero(stateSpace_);
         
            
            
            Spacy::ModifiedPPCG::TriangularConstraintPreconditioner triH(surrogateStateSolver_,controlSolver_,surrogateAdjointSolver_,minusB_,totalSpace_,stateSpace_);
            Spacy::ModifiedPPCG::Solver modPPCG_(triH, M_, minusB_, totalSpace_, stateSpace_);
            
            
            modPPCG_.setTerminationCriterion(  CG::Termination::AdaptiveRelativeEnergyError() );
            
            adjointAcc = std::max(1e-12,Spacy::Mixin::get(adjointAcc));
            stateAcc = std::max(1e-12,Spacy::Mixin::get(stateAcc));
            
            LOG(log_tag, "AdjointAcc: ", adjointAcc);
            LOG(log_tag, "StateAcc: ", stateAcc);
            LOG(log_tag, "ModifiedCGAcc: ", modifiedCGAcc);
            
//             adjointSolver_.setRelativeAccuracy(adjointAcc);
//             stateSolver_.setRelativeAccuracy(stateAcc);
            modPPCG_.setRelativeAccuracy(modifiedCGAcc);
            modPPCG_.setMaxSteps(5000);
            
            for(unsigned int step = 1; step < getMaxSteps(); step++)
            {
                LOG_SEPARATOR(log_tag);
                
                // Check regularization part
                dp = adjointSolver_(stateSpace_.project(-rx)); //Vermerk: b ^= -c in BA

                logIterationsAdjoint(adjointSolver_.getIterations());
                LOG(log_tag, "Adjoint Solver Number of Iterations: ", adjointSolver_.getIterations());
                dual_iterations_ += adjointSolver_.getIterations();
                
//                 if(adjointSolver.indefiniteOperator())
//                 {
//                     signalConvex_(true,false);
//                     return x;
//                 }

                p += dp;
                
                rx = primalSpace_.embed( controlSpace_.project(rx) ); // ry -> 0
                rx += primalSpace_.embed( minusB_.transposed(dp) ); // ru -= B^T(dp)
                

                dx = primalSpace_.project(modPPCG_(totalSpace_.embed(-rx)));
                
                LOG(log_tag, "PPCG Number of Iterations: ", modPPCG_.getIterations());
                

                convex_flag_ = modPPCG_.isPositiveDefinite();
                
//                         if(modPPCG_.isSingular())
//                         {
//                             LOG_INFO(log_tag, "inexactCG is singular");
//                             signalConvex_(true,false);
//                             return x;
//                         }
                
                if(!convex_flag_)
                {
                    definiteness_ = DefiniteNess::Indefinite;
                    result_ = Result::TruncatedAtNonConvexity;
                    return x;
                }
                
                bp = -( minusB_( controlSpace_.project(dx) ) );
                
                auto ddy = stateSolver_(bp);
                
                
                dx = primalSpace_.embed(ddy) + primalSpace_.embed( controlSpace_.project(dx) );
               
                logIterationsState(stateSolver_.getIterations());
                LOG(log_tag, "State Solver Number of Iterations: ", stateSolver_.getIterations());
                proj_iterations_ += stateSolver_.getIterations();
                
//                 if(stateSolver.indefiniteOperator())
//                 {
//                     signalConvex_(true,false);
//                     return x;
//                 }
                
                Mdx=M_(dx);
                
                Real dxMdx = Mdx(dx);
                
                omega = -rx(dx)/dxMdx;
                LOG(log_tag, "omega: ", omega);
                
                
                auto normX = sqrt(M_(x)(x));
                normDx_ =std::fabs(get(omega)) *sqrt(dxMdx);
                
                LOG(log_tag, "normX (old): ", normX, "normDx_: ", normDx_);
                
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
                
                if (normDxOld_.get() != 0) std::cout << "Contraction: " << normDx_/normDxOld_ << std::endl;
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
//                 LOG_INFO(log_tag, "Not Converged: normDxOld = 0");
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
