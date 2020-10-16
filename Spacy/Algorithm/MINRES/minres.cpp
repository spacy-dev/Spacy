#include "minres.hh"

#include <Spacy/Spaces/ProductSpace.h>
#include <Spacy/ZeroVectorCreator.h>
#include <Spacy/Util/Mixins/Get.h>
#include <Spacy/Util/Logger.h>
#include <Spacy/Util/Log.h>
#include <iostream>
#include "blockDiagPrecond.hh"

#define LAPACK_COL_MAJOR 102

namespace Spacy
{
    namespace
    {
        
    }
    namespace MINRES
    {
        
        DEFINE_LOG_TAG(static const char* log_tag = "MINRES");
        
        Solver::Solver(CallableOperator H,CallableOperator P) :
        H_(H), P_(P)
        {
        }
        
        Vector Solver::operator()(const Vector &b) const
        {
            auto x = 0*P_(b);
            auto result = minresLoop(x,b);
            return result;
        }
        
        Vector Solver::operator()(const Vector &b, const VectorSpace &domain) const
        {
            auto x = zero(domain);
            auto result = minresLoop(x,b);
            return result;
        }
        
        
        bool Solver::convergenceTest(Real norm_r0, Real norm_r, const Vector &x  ) const
        {
            
            if(norm_r/norm_r0 < getRelativeAccuracy() || norm_r < eps())
            {
                
                LOG(log_tag, "Norm Solution: ", sqrt(x(x)) );
                return true;
            }
            else
                return false;
        }//*/
        
        
        Vector Solver::minresLoop(Vector x,const Vector& b) const
        {
            LOG_INFO(log_tag, "Starting MINRES-Solver.");
            LOG_SEPARATOR(log_tag,"");
            
            LOG(log_tag, "Relative Accuracy: ", getRelativeAccuracy());
            LOG(log_tag, "Eps: ", eps());
            
            iterations_=0;
            
            auto v_old = b*0.0;
            LOG(log_tag, "Iteration: ", 1);
            
            auto v = b - H_(x);
            
            auto v_new = v_old;
            auto w = x*0.0;
            auto w_old = w;
            auto w_new = w;
            
            
            auto z = P_(v);
            iterations_++;
            
            auto z_new = z;
            auto gamma = sqrt(v(z));
            
            z *= 1.0/get(gamma);
            v *= 1.0/get(gamma);
            
            auto eta_old = gamma;
            Real s_old = 0.0;
            Real s = 0.0;
            Real c_old = 1.0;
            Real c = 1.0;
            auto norm_r0 = norm(eta_old);
            
            for(unsigned step = 1; step < getMaxSteps(); step++)
            {
                
                LOG_SEPARATOR(log_tag,step);
                
                auto Hzmv = H_(z)- gamma * v_old;
                auto delta = z(Hzmv);
                
                v_new = Hzmv - delta * v;
                
                z_new = P_(v_new);
                iterations_++;
                auto gamma_new = sqrt(v_new(z_new));
                
                
                z_new *= 1.0/get(gamma_new);
                v_new *= 1.0/get(gamma_new);
                
                auto alpha_0 = c * delta  - c_old * s * gamma;
                auto alpha_1 = sqrt(alpha_0 * alpha_0 + gamma_new * gamma_new);
                auto alpha_2 = (s * delta + c_old * c * gamma);
                auto alpha_3 = s_old * gamma;
                
                
                auto c_new = alpha_0/alpha_1;
                auto s_new = gamma_new/alpha_1;
                
                
                w_new = (1.0/alpha_1) * (z - alpha_3 * w_old - alpha_2 * w);
                x = x + c_new * eta_old *w_new;
                eta_old = -s_new * eta_old;
                
                auto norm_r = norm(eta_old);
                auto converged = convergenceTest(norm_r0, norm_r,x);
                
                v_old = v;
                v = v_new;
                w_old = w;
                w = w_new;
                z = z_new;
                gamma = gamma_new;
                s_old = s;
                s = s_new;
                c_old = c;
                c = c_new;
                
                LOG(log_tag, "Norm residuum: ",  norm_r);
                
                if(converged)
                {
                    
                    LOG(log_tag, "MINRES converged in step: ",  iterations_);
                    return x;
                }
            }
            
            LOG_INFO(log_tag, "Number of MINRES iterations > Maximum Number of Iterations");
            iterations_ = getMaxSteps();
            return x;
        }
    }
}
