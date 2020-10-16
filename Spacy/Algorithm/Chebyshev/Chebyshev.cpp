#include "Chebyshev.h"

#include <Spacy/VectorSpace.h>
#include <Spacy/ZeroVectorCreator.h>
#include <Spacy/ScalarProduct.h>
#include <Spacy/LinearOperator.h>

#include <iostream>
#include "linalg/lapack/lapacke.h"

namespace Spacy
{
    namespace Chebyshev
    {
    
        
      void computeEigsFromCGCoefficients(const std::vector<Spacy::Real>& alpha,const std::vector<Spacy::Real>& beta, double& eigMin, double& eigMax)
      {
        int a_size = alpha.size();
        std::vector<double> eig(a_size);
        std::vector<double> diagonal(a_size);
        std::vector<double> subDiagonal(a_size-1);
        std::vector<double> z_vec(a_size);
        std::vector<int> isuppz_vec(2*a_size);
        std::vector<int> m(a_size);

        for(auto i=0u; i<a_size; ++i)
        {
            eig[i] = 0.0;
            if(i==0)
                diagonal[i] = 1./alpha[i].get();
            else
                diagonal[i] = 1./alpha[i].get() + beta[i-1].get()/alpha[i-1].get();
            if(i!=0)
                subDiagonal[i-1] = std::sqrt(beta[i-1].get())/alpha[i-1].get();
            z_vec[i] = 0.0;
            isuppz_vec[i] = 0;
            isuppz_vec[i+a_size] = 0;
            m[i] = 0;
        }

        LAPACKE_dstevr(LAPACK_COL_MAJOR,'N','A',a_size,&diagonal[0],&subDiagonal[0],0.0,0.0,0,0,1e-8,&m[0],&eig[0],&z_vec[0],a_size,&isuppz_vec[0]);


        eigMin = std::min(eig[0],eigMin);
        eigMax = std::max(eig[a_size-1],eigMax);
     }

     void ChebyshevPreconditioner::setSpectralBounds(Real gamma_min, Real gamma_max)
     {
         gamma_max_= gamma_max;
         gamma_min_ = gamma_min;
     }

    unsigned ChebyshevPreconditioner::getIterations() const
    {
            Spacy::Real condition = gamma_max_/gamma_min_;
            Spacy::Real theta = (sqrt(condition)-1)/(sqrt(condition)+1);
           return std::max(floor(log(get(getRelativeAccuracy()))/log(get(theta))),0.0) + 1; 
    }

     
        ChebyshevPreconditioner::ChebyshevPreconditioner(Spacy::CallableOperator A,Spacy::CallableOperator P, ::Spacy::Real gamma_max, ::Spacy::Real gamma_min)

              : A_(A), P_(P), gamma_max_(gamma_max), gamma_min_(gamma_min)
        {
        }

        Vector ChebyshevPreconditioner::operator()(const Vector& b) const
        {

            Spacy::Real condition = gamma_max_/gamma_min_;
            Spacy::Real theta = (sqrt(condition)-1)/(sqrt(condition)+1);
            int N = std::max(floor(log(get(getRelativeAccuracy()))/log(get(theta))),0.0) + 1;

            Spacy::Real c = 0.5*(gamma_max_-gamma_min_);
            Spacy::Real alpha = 0.5*(gamma_max_+gamma_min_);
            Spacy::Real beta_const = -0.5*c*c/alpha;
            Spacy::Real beta = 0;
            Spacy::Real gamma = -alpha;

            
            Spacy::Vector dx = P_(b);
            Spacy::Vector x_nm1 = zero(dx.space());
            Spacy::Vector x_n = zero(dx.space());

            for(int i=0; i<N ; i++)
            {

                if(i == 1)
                    beta = beta_const;
                if(i >= 2)
                    beta = (1./gamma)*0.25*c*c;
                if(i >= 1)
                    gamma = -(alpha+beta);

                
                 //x_np1 = -(1./gamma)*(dx+alpha*x_n+beta*x_nm1);


                 dx+=alpha*x_n;
                 dx+=beta*x_nm1;
                 dx *= -1.0/get(gamma);
                 x_nm1=x_n;
                 x_n=dx;

                if(i < N-1)
                 {
                    dx = P_(b-A_(x_n));
                 }
  //Explicit computation of residual if preferrable
                
            }

            if(verbose())
            {
                std::cout << "Cheb terminating after " << N << " iterations " << std::endl;
            }

            iterations_ = N;
  //          std::cout << "Used BPX-calls:" << N << " Eigs: [" << gamma_min_ << "," << gamma_max_ << "]" << std::endl;
            return x_n;
        }
    }
}
