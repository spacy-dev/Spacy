#include "ModifiedPPCG.h"

#include <Spacy/ZeroVectorCreator.h>
#include <Spacy/Util/Exceptions.h>
#include <Spacy/Util/Cast.h>
#include <iostream>
#include <Spacy/Util/Logger.h>
#include <Spacy/Util/Log.h>

#include <limits>
#include <typeinfo>

//#include "linalg/lapack/lapacke.h"

#include "TriangConstraintPrec.h"


#define LAPACK_COL_MAJOR 102

namespace Spacy
{
    
    namespace ModifiedPPCG
    {
        using namespace Spacy;

        
        
        DEFINE_LOG_TAG(static const char* log_tag = "ModPPCG");
        
        Solver::Solver(const TriangularConstraintPreconditioner& P, const Operator& My, const Operator& Mu, const OperatorWithTranspose& minusB, const VectorSpace& domainIn, const PPCG_Test& test, Real theta )
        : OperatorBase( domainIn, domainIn ), P_(P), My_(My), Mu_(Mu), minusB_(minusB), stateSpace_(My.domain()), controlSpace_(Mu.domain()), adjointSpace_(minusB_.range()), test_(test), stee_(), theta_(theta)
        {
        }
        
        
        bool Solver::isPositiveDefinite() const noexcept
        {
            return definiteness_ == DefiniteNess::PositiveDefinite;
        }
        
        
        void Solver::setRegularizationOperators(const CallableOperator& R)
        {
            R_=R;
        }
        
        Vector Solver::operator()( const Vector& b) const
        {
            return solve(b);
        }
        
        unsigned int Solver::getIterations() const
        {
            return iterations_;
        }

        unsigned int Solver::getIterationsInLifeTime() const
        {
            return iterationsInLifeTime_;
        }

        bool Solver::isSingular() const
        {
            return result_ == Result::TooSingular;
        }
        
        Vector Solver::solve( const Vector& b ) const
        {
            
            terminate_.clear();
            stee_.clear();
            stee_.setLookAhead(2);
            
            iterations_ = 0;
            
            auto x = zero(domain());
            auto y = zero(stateSpace_);
            auto u = zero(controlSpace_);
            auto p = zero(adjointSpace_);
            
            
            // Todo check if correct for simplified normal step
            terminate_.setEps( get( eps() ) );
            stee_.setEps( get( eps() ) );
                        
            terminate_.setRelativeAccuracy( get( getRelativeAccuracy() ) );
            terminate_.setAbsoluteAccuracy( get( getAbsoluteAccuracy() ) );
                        
            stee_.setRelativeAccuracy( get( getRelativeAccuracy() ) );
            stee_.setAbsoluteAccuracy( get( getAbsoluteAccuracy() ) );
            
            auto H_py = zero(stateSpace_);
            auto H_pu = zero(controlSpace_);
            auto ry = zero(stateSpace_);
            auto ru = -controlSpace_.project(b);
            
            Real beta = 0.0;
            Real alpha = 0.0;
            
            std::vector<Spacy::Real> alpha_vec;
            std::vector<Spacy::Real> beta_vec;
            
            
            auto gy=zero(stateSpace_);
            auto gu=zero(controlSpace_);
            auto gp=zero(adjointSpace_);
            //P_.apply(g,range().embed(-r));
            P_(-ry,-ru,gy,gu,gp);
            
            if(verbose()) { std::cout << "o" << std::flush; }
            iterations_++;
            iterationsInLifeTime_++;
            
            auto py=gy;
            auto pu=gu;
            auto pp=gp;
            
            Real sigma = P_.getSigma();
            
            
            if(sigma < 0)
            {
                if(verbose()) LOG_INFO(log_tag, "Fail Preconditioner not positive definite in PPCG !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
                if(verbose()) LOG(log_tag, "sigma: ", sigma);
                if(verbose()) LOG_INFO(log_tag, "In iteration: 1");
                definiteness_ = DefiniteNess::Indefinite;
                result_ = Result::TruncatedAtNonConvexity;
                return x;
            }

            // P_p  = P_(p)
            auto P_py = ry;  // required for termination criteria
            auto P_pu = ru;  // required for termination criteria
            
            stee_.setRelativeAccuracy( get( getRelativeAccuracy() ) );
//             stee_.setVerbosityLevel(2);
            terminate_.setRelativeAccuracy( get( getRelativeAccuracy() ) );
            
            double ratioY = 1;
            double ratioU = 1;
            
            for(unsigned step = 1; step <= getMaxSteps(); step++)
            {
                H_py = My_(py);
                H_pu = Mu_(pu);
                
                // p_H_p = p_M_p because p is in the nullspace of tilde C
                Real p_H_p = H_py(py) * ratioY + H_pu(pu) * ratioU;
                
                // remaining parts of H_p, which include the update of the Lagrangian multiplier
                H_pu += minusB_.transposed(pp);
                
                // note that \tilde A cannot be applied, so we use the first component of -tilde P_ p
                H_py -= P_py;
                
                Real p_P_p = P_py(py) * ratioY + P_pu(pu) * ratioU;
                
                
                
                
                if(p_H_p < 0)
                {
                    if(verbose()) LOG(log_tag, "Problem not positive definite: pHp: ", p_H_p);
                    if(verbose()) LOG(log_tag, "In iteration: ", step);
                    definiteness_ = DefiniteNess::Indefinite;
                    result_ = Result::TruncatedAtNonConvexity;
                    return x;
                }
                
                alpha = sigma/p_H_p;
                alpha_vec.push_back(get(alpha));
                
                stee_.update(get( alpha ), get( p_H_p ), get( p_P_p ), get( sigma ), x );
                terminate_.update( get( alpha ), get( p_H_p ), get( p_P_p ), get( sigma ), x );
                
                if ( vanishingStep( step, alpha_vec, beta_vec ) )
                {
                    if(verbose()) std::cout << std::endl;
                    if ( step == 1 )
                    {
                        return x.space().embed(py) + x.space().embed(pu) + x.space().embed(pp);
                    }
                    result_ = Result::Converged;
                    return x;
                }
                
                y += alpha * py;
                u += alpha * pu; 
                p += alpha * pp;
                x  = x.space().embed(y) + x.space().embed(u) + x.space().embed(p);
                ry += alpha * H_py;
                ru += alpha * H_pu;
                
                
                if ( stee_() )
                {
                    if(verbose())
                    {
                        std::cout << std::endl;
//                         test_(y,u,p,-controlSpace_.project(b));
//                         test_.diffVec(stateSpace_.project(x),controlSpace_.project(x),adjointSpace_.project(x));
                    }
                    return x;
                }
                
                P_(-ry,-ru,gy,gu,gp);
                //P_.apply(g,range().embed(-r));
                if(verbose()) { std::cout << "o" << std::flush; }
                iterations_++;
                iterationsInLifeTime_++;
                
                Real sigmaNew = P_.getSigma();
                                
                beta = sigmaNew/sigma;
                beta_vec.push_back(get(beta));
                
                sigma = sigmaNew;
                
                
                if(sigma < 0 )
                {
                    if(verbose()) LOG_INFO(log_tag, "Fail: Preconditioner not positive definite in PPCG !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
                    if(verbose()) LOG(log_tag, "sigma: ", sigmaNew);
                    if(verbose()) LOG(log_tag, "In iteration: ", step);
                    definiteness_ = DefiniteNess::Indefinite;
                    result_ = Result::TruncatedAtNonConvexity;
                    return x;
                }
                
                py *= get(beta);
                pu *= get(beta);
                pp *= get(beta);
                py += gy;
                pu += gu;
                pp += gp;
                
                P_py *= get( beta );
                P_pu *= get( beta );
                P_py += ry;
                P_pu += ru;
            }
            result_ = Result::NotConverged;
            return x;
        }
        
        bool Solver::vanishingStep(unsigned step , const std::vector<Real> &alpha_vec, const std::vector<Real> &beta_vec) const
        {
            if ( terminate_.vanishingStep() )
            {
                result_ = Result::Converged;
                
                return true;
            }
            return false;
        }
        
        double Solver::getConditionNumber()
        {
            return condPrecondMatrix;
        }
        
    }
}
