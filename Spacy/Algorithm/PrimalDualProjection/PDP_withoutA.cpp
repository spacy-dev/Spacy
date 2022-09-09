#include "PDP_withoutA.hh"

#include <Spacy/Spaces/ProductSpace.h>
#include <Spacy/ZeroVectorCreator.h>
#include <Spacy/Util/Mixins/Get.h>
#include <Spacy/Util/Logger.h>
#include <Spacy/Util/Log.h>
#include <iostream>
#include <chrono>

#define LAPACK_COL_MAJOR 102

using namespace std::chrono;

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
        
        
        Solver::Solver(Operator M, Operator Mys, Operator Mus, Spacy::SolverBase stateSolver, Spacy::SolverBase adjointSolver, OperatorWithTranspose minusB, OperatorWithTranspose minusBs, Spacy::SolverBase surrogateStateSolver, Spacy::SolverBase surrogateAdjointSolver, Operator controlSolver, Operator controlSolverS, const VectorSpace& totalSpace, const VectorSpace& totalSpaceS) :
        M_(std::move(M)), My_s(std::move(Mys)), Mu_s(std::move(Mus)), stateSolver_(std::move(stateSolver)), adjointSolver_(std::move(adjointSolver)), minusB_(std::move(minusB)), minusB_s(std::move(minusBs)), surrogateStateSolver_(std::move(surrogateStateSolver)), surrogateAdjointSolver_(std::move(surrogateAdjointSolver)), totalSpace_(totalSpace), totalSpace_s(totalSpaceS), adjointSpace_(stateSolver_.domain()), controlSpace_(controlSolver.domain()), stateSpace_(stateSolver_.range()), controlSolver_((controlSolver)), controlSolver_s((controlSolverS))
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
            ::Spacy::Vector x = zero(totalSpace_);
            ::Spacy::Vector result = zero(b.space());
            theta_ = 0.0;
            convex_flag_ = false;
            bool logConvex = true;
            result_ = Result::Failed;
            definiteness_ = DefiniteNess::PositiveDefinite;
            
//             while(!convex_flag_)
//             {
                result = PDPLoop(x,b);
                convex_flag_ = isPositiveDefinite();
                
                if(verbose()) LOG(log_tag, "PPCG is Convex: ", convex_flag_);
                if(verbose()) LOG(log_tag, "Energy in PPCG is Convex: ", signalConvex_(false,false));
                
                if(!signalConvex_(false,false))
                {
                    return result;
                }
                
                if(!convex_flag_)
                    logConvex = false;
//             }
            
            if(!logConvex)
                definiteness_ = DefiniteNess::Indefinite;
            
            return result;
        }
        
        bool Solver::isPositiveDefinite() const
        {
            return definiteness_ == DefiniteNess::PositiveDefinite;
        }
        
        Vector Solver::PDPLoop(Vector x, const Vector& b) const
        { 
            if(verbose()) { LOG_SEPARATOR(log_tag); }
            
            if(verbose()) LOG(log_tag, "Regularization Parameter: ", theta_);
            if(verbose()) LOG(log_tag, "Relative Accuracy: ", getRelativeAccuracy());
            if(verbose()) LOG(log_tag, "Eps: ", eps());
            
            convex_flag_ = true;
            energyIsConvex_ = true;
            definiteness_ = DefiniteNess::PositiveDefinite;
            normDxOld_ = 0;
            
            Real stateAcc = stateAcc_;
            Real adjointAcc = adjointAcc_;
            Real modifiedCGAcc = modifiedCGAcc_;
            
            sumNormSquaredUpdate_ = 0.0;
            
            Real omega = 0.0;
            
            auto y      = stateSpace_.project(x);
            auto u      = controlSpace_.project(x);
            auto p      = adjointSpace_.project(x);
            
            auto ry     = - stateSpace_.project(b);
            auto ru     = - controlSpace_.project(b);
            
            auto d      = zero(totalSpace_);
            auto dy     = zero(stateSpace_);
            auto du     = zero(controlSpace_);
            auto dp     = zero(adjointSpace_);
            
            auto Mdx    = zero(M_.domain());
            auto bp     = zero(adjointSpace_);
         
            
            Spacy::ModifiedPPCG::PPCG_Test testPPCG(surrogateStateSolver_,controlSolver_s,surrogateAdjointSolver_,My_s,Mu_s,minusB_s);
            Spacy::ModifiedPPCG::TriangularConstraintPreconditioner triH(surrogateStateSolver_,controlSolver_s,surrogateAdjointSolver_,minusB_s,totalSpace_);
            Spacy::ModifiedPPCG::Solver modPPCG_(triH, My_s, Mu_s, minusB_s, totalSpace_, testPPCG);
            
            
            modPPCG_.setTerminationCriterion(  CG::Termination::AdaptiveRelativeEnergyError() );
            modPPCG_.setVerbosity(verbose());
            
            adjointAcc = std::max(1e-12,Spacy::Mixin::get(adjointAcc));
            stateAcc = std::max(1e-12,Spacy::Mixin::get(stateAcc));
            
            if(verbose()) LOG(log_tag, "AdjointAcc: ", adjointAcc);
            if(verbose()) LOG(log_tag, "StateAcc: ", stateAcc);
            if(verbose()) LOG(log_tag, "ModifiedCGAcc: ", modifiedCGAcc);
            
            modPPCG_.setRelativeAccuracy(modifiedCGAcc);
            modPPCG_.setMaxSteps(5000);
            
            iterations_ = 0;
            PPCGiterations_ = 0;
            MGiterations_ = 0;
            
            auto startTime = high_resolution_clock::now();
            auto endTime = high_resolution_clock::now() - startTime;
            
            for(unsigned int step = 1; step < getMaxSteps(); step++)
            {
                if(verbose()) { LOG_SEPARATOR(log_tag); }
                
                // Check regularization part
                dp = adjointSolver_(-ry);
                
                if(verbose()) logIterationsAdjoint(adjointSolver_.getIterations());
                if(verbose()) LOG(log_tag, "Adjoint Solver Number of Iterations: ", adjointSolver_.getIterations());
                dual_iterations_ += adjointSolver_.getIterations();
                
//                 if(adjointSolver.indefiniteOperator())
//                 {
//                     signalConvex_(true,false);
//                     return x;
//                 }

                p += dp;
                
                ry = zero(stateSpace_); // ry -> 0
                ru += minusB_.transposed(dp); // ru -= B^T(dp)
                
                startTime = high_resolution_clock::now();
                d = totalSpace_.embed(modPPCG_(-ru)); 
                endTime += high_resolution_clock::now() - startTime;
                
                if(verbose())
                {
                    LOG_INFO(log_tag, "nach Gitterwechsel: ");
                    
                    testPPCG(surrogateStateSolver_.range().project(d),
                            Mu_s.domain().project(d),
                            surrogateAdjointSolver_.range().project(d),
                            Mu_s.domain().project(ru));
                }
                
                du = controlSpace_.project(totalSpace_.embed(d));
                
                if(verbose()) LOG(log_tag, "PPCG Number of Iterations: ", modPPCG_.getIterations());
                

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
                    iterations_ = getMaxSteps();
                    PPCGiterations_ = modPPCG_.getIterationsInLifeTime();
                    MGiterations_ = triH.getBPXCallsInLifeTime() + dual_iterations_ + proj_iterations_;
                    
                    x = x.space().embed( y ) + x.space().embed( u ) + x.space().embed(p);
                    return x;
                }
                
                bp = -minusB_( du );
                
                
                dy = stateSolver_(bp);
                
               
                if(verbose()) logIterationsState(stateSolver_.getIterations());
                if(verbose()) LOG(log_tag, "State Solver Number of Iterations: ", stateSolver_.getIterations());
                proj_iterations_ += stateSolver_.getIterations();
                
//                 if(stateSolver.indefiniteOperator())
//                 {
//                     signalConvex_(true,false);
//                     return x;
//                 }
                
                Mdx=M_( M_.domain().embed(dy) + M_.domain().embed(du) );
                
                Real dxMdx = Mdx( M_.domain().embed(dy) + M_.domain().embed(du) );
                
                omega = -ru(du)/dxMdx;
                if(verbose()) LOG(log_tag, "omega: ", omega);
                
                auto Mx = M_( M_.domain().embed(y) + M_.domain().embed(u) );
                auto normX = sqrt( Mx( M_.domain().embed(y) + M_.domain().embed(u) ) );
                normDx_ =std::fabs(get(omega)) * sqrt(dxMdx);
                
                if(verbose()) LOG(log_tag, "normX (old): ", normX, "normDx_: ", normDx_);
                
                energyIsConvex_ = signalConvex_(false,false);
                bool converged = convergenceTest(convex_flag_, energyIsConvex_);
                
                if(converged)// && !projected_)
                {
                    result_ = Result::Converged;
                    iterations_ = step;
                    PPCGiterations_ = modPPCG_.getIterationsInLifeTime();
                    MGiterations_ = triH.getBPXCallsInLifeTime() + dual_iterations_ + proj_iterations_;
                    
                    if(verbose()) 
                    {
                        LOG(log_tag, "Number of Iterations  PDP: ", iterations_);
                        LOG(log_tag, "Number of Dual Iterations: ", dual_iterations_);
                        LOG(log_tag, "Number of PPCG Iterations: ", PPCGiterations_);
                        LOG(log_tag, "Number of Proj Iterations: ", proj_iterations_);
                        LOG(log_tag, "Number of MG calls in  CG: ", dual_iterations_ + proj_iterations_);
                        LOG(log_tag, "Number of MG calls   PPCG: ", triH.getBPXCallsInLifeTime());
                        LOG(log_tag, "Number of MG calls  total: ", MGiterations_);
                        LOG(log_tag, "Average no. MG calls /  A: ", MGiterations_/(2.0*PPCGiterations_));
                        std::cout << std::endl;
                    }
                    
                    MillisecondsPPCG = duration_cast< milliseconds >( endTime ).count();
//                     std::cout << "Time spent in PPCG: " << duration_cast< seconds >( endTime ).count() << "s, " 
//                     << duration_cast< milliseconds >( endTime ).count() % 1000 << "ms." << std::endl;
                    
                    x = x.space().embed( y ) + x.space().embed( u ) + x.space().embed(p);
                    return x;
                }
                
                y += get(omega)*dy;
                u += get(omega)*du;
                ry = omega * stateSpace_.project(Mdx); // kein += weil ry bereits 0 ist
                ru += omega * controlSpace_.project(Mdx);
                
                auto vSE   = y + stateSolver_(minusB_(u));
                auto errSE = sqrt(vSE(vSE));
                if(verbose()) LOG(log_tag, "Error State Equation: ", errSE);
                
                normDxOld_ = normDx_; 
            }
            if(verbose()) LOG_INFO(log_tag, "pcgStep > getMaxSteps()");
            result_ = Result::NotConverged;
            iterations_ = getMaxSteps();
            PPCGiterations_ = modPPCG_.getIterationsInLifeTime();
            MGiterations_ = triH.getBPXCallsInLifeTime() + dual_iterations_ + proj_iterations_;
            
            x = x.space().embed( y ) + x.space().embed( u ) + x.space().embed(p);
            return x;
        }
        
        bool Solver::convergenceTest(bool isConvex, bool energyIsConvex) const
        {
            
            if(!isConvex)
            {
                if(verbose()) LOG_INFO(log_tag, "Not Converged: !isConvex");
                return false;
            }
            if(!energyIsConvex)
            {
                if(verbose()) LOG_INFO(log_tag, "Not Converged: !energyIsConvex");
                return false;
            }
            
            
            sumNormSquaredUpdate_ += get(normDx_*normDx_);
            
            if(normDxOld_ == 0) //normDxOld_ = 0 -> noch im ersten Schritt, erzwingt so mind. 2 Schritte (es ei denn dx ist sehr sehr klein)
                if(verbose()) logTheta(0.0);
            
            else
                if(verbose()) logTheta(get(normDx_/normDxOld_));
            
            if(normDx_ < eps())
            {
                if(verbose()) LOG_INFO(log_tag, "Converged: normDx_ < eps()");
                return true;
            }
            
            if(normDxOld_ == 0) //normDxOld_ = 0 -> noch im ersten Schritt, erzwingt so mind. 2 Schritte (es ei denn dx ist sehr sehr klein)
            {
//                 LOG_INFO(log_tag, "Not Converged: normDxOld = 0");
                return false;
            }
            else
            {
                auto contraction = normDx_/normDxOld_; //wenn contraction > desiredContraction -> erhöhe lookAhead
                if(contraction >= 1)
                {
                    if(verbose()) LOG_INFO(log_tag, "Not Converged: contraction >= 1");
                    return false;
                }
                
                if(ThetaMin > contraction) ThetaMin = contraction; 
                if(ThetaMax < contraction) ThetaMax = contraction;
//                 Theta = contraction;
                if(Theta < contraction) Theta = contraction;
//                 std::cout << "Contraction: " << contraction << std::endl;
                
                auto normSolutionEstimate = sqrt(sumNormSquaredUpdate_);
                auto normErrorEstimate = (Theta/sqrt(1-Theta*Theta))*normDx_;
                
//                 std::cout << " Convergence Test: " << normErrorEstimate << " < " << getRelativeAccuracy() << " * " << normSolutionEstimate << " Contraction: " << contraction << std::endl;
                
                if(normErrorEstimate < getRelativeAccuracy() * normSolutionEstimate)
                {std::cout << "Verhältnis ThetaMin/ThetaMax: " << ThetaMin/ThetaMax << std::endl;
                    if(verbose()) std::cout << " Convergence Test: " << normErrorEstimate << " < " << getRelativeAccuracy() << " * " << normSolutionEstimate << " Contraction: " << contraction << std::endl;
                    if(verbose()) LOG(log_tag, " Exceeded desired accuracy by factor: ", normErrorEstimate/(getRelativeAccuracy() * normSolutionEstimate) );
                    return true;
                }
                else
                {
                    if(verbose()) std::cout << " Convergence Test: " << normErrorEstimate << " < " << getRelativeAccuracy() << " * " << normSolutionEstimate << " Contraction: " << contraction << std::endl;
                }
            }
            return false;
        }//*/
    }
}
