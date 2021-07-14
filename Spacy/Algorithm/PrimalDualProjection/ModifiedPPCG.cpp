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
        
        Solver::Solver(const TriangularConstraintPreconditioner& P, Operator M, OperatorWithTranspose minusB,
                       const VectorSpace& domainIn, const VectorSpace& stateSpace, Real theta )
        : OperatorBase( domainIn, domainIn ), P_(P), M_(M), minusB_((minusB)), stateSpace_(stateSpace), 
        primalSpace_(M_.domain()), adjointSpace_(minusB_.range()), theta_(theta)
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
        
        Vector Solver::solve( const Vector& b) const
        {
            
            terminate_.clear();
            
            iterations_ = 0;
            
            auto x = zero(domain());
            
            
            // Todo check if correct for simplified normal step
            terminate_.setEps( get( eps() ) );
                        
            terminate_.setRelativeAccuracy( get( getRelativeAccuracy() ) );
            terminate_.setAbsoluteAccuracy( get( getAbsoluteAccuracy() ) );
            
            auto H_p = zero(primalSpace_);
            auto r = -primalSpace_.project(b); 
            
            Real beta = 0.0;
            Real alpha = 0.0;
            
            std::vector<Spacy::Real> alpha_vec;
            std::vector<Spacy::Real> beta_vec;
            
            
            auto g=zero(domain());
            //P_.apply(g,range().embed(-r));
            g=-P_(range().embed(r));
            std::cout << "o" << std::flush;
            iterations_++;
            iterationsInLifeTime_++;
            
            auto p = g;
            
            Real sigma = P_.getSigma();
            
            
            if(sigma < 0)
            {
                LOG_INFO(log_tag, "Fail Preconditionernot positive definite in PPCG !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
                LOG(log_tag, "sigma: ", sigma);
                LOG_INFO(log_tag, "In iteration: 1");
                definiteness_ = DefiniteNess::Indefinite;
                result_ = Result::TruncatedAtNonConvexity;
                return x;
            }

            // P_p  = P_(p)
            auto P_p = r;  // required for termination criteria
            
            terminate_.setRelativeAccuracy( get( getRelativeAccuracy() ) );
            
            for(unsigned step = 1; step <= getMaxSteps(); step++)
            {
                
                H_p = M_(primalSpace_.project(p));
                if(get(theta_) != 0.0) H_p+=theta_*R_(primalSpace_.project(p));
                
                // p_H_p = p_M_p because p is in the nullspace of tilde C
                Real p_H_p = H_p(primalSpace_.project(p));
                
                // remaining parts of H_p, which include the update of the Lagrangian multiplier
                H_p += primalSpace_.embed( minusB_.transposed(adjointSpace_.project(p)) );
                
                // note that \tilde A cannot be applied, so we use the first component of -tilde P_ p
                H_p -= primalSpace_.embed( stateSpace_.project(P_p) );
                
                Real p_P_p = P_p(primalSpace_.project(p));
                
                
                
                
                if(p_H_p < 0)
                {
                    LOG(log_tag, "Problem not positive definite: pHp: ", p_H_p);
                    LOG(log_tag, "In iteration: ", step);
                    definiteness_ = DefiniteNess::Indefinite;
                    result_ = Result::TruncatedAtNonConvexity;
                    return x;
                }
                
                alpha = sigma/p_H_p;
                alpha_vec.push_back(get(alpha));
                
                terminate_.update( get( alpha ), get( p_H_p ), get( p_P_p ), get( sigma ), x );
                
                if ( vanishingStep( step, alpha_vec, beta_vec ) )
                {
                    std::cout << std::endl;
                    if ( step == 1 )
                    {
                        x +=p;
                    }
                    result_ = Result::Converged;
                    return x;
                }
                
                x += alpha * p;
                r += alpha * H_p;
                
                if ( terminate_() )
                {
                    std::cout << std::endl;
                    return x;
                }
                
                g = -P_(range().embed(r));
                //P_.apply(g,range().embed(-r));
                std::cout << "o" << std::flush;
                iterations_++;
                iterationsInLifeTime_++;
                
                Real sigmaNew = P_.getSigma();
                                
                beta = sigmaNew/sigma;
                beta_vec.push_back(get(beta));
                
                sigma = sigmaNew;
                
                
                if(sigma < 0 )
                {
                    LOG_INFO(log_tag, "Fail: Preconditioner not positive definite in PPCG !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
                    LOG(log_tag, "sigma: ", sigmaNew);
                    LOG(log_tag, "In iteration: ", step);
                    definiteness_ = DefiniteNess::Indefinite;
                    result_ = Result::TruncatedAtNonConvexity;
                    return x;
                }
                
                p *= get(beta);
                p += g;
                
                P_p *= get( beta );
                P_p += r;
                
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
