#include "chebyshev.hh"

//#include <Spacy/Spaces/productSpace.h>
//#include <Spacy/vectorSpace.h>
//#include <Spacy/zeroVectorCreator.h>
#include <Spacy/ScalarProduct.h>
#include <Spacy/LinearOperator.h>

//#include <utility>
#include <iostream>
//using namespace std;

namespace Spacy
{
    namespace PPCG
    {
        ChebyshevPreconditioner::ChebyshevPreconditioner(Spacy::CallableOperator A, Spacy::CallableOperator P,
                                                         Spacy::Real gamma_max, Spacy::Real gamma_min
                                                         //,const VectorSpace& domain, const VectorSpace& range
                                                         , std::function<::Spacy::Vector(::Spacy::Vector,::Spacy::Vector)> transfer)
            : //OperatorBase(domain,range),
              A_(A), P_(P), gamma_max_(gamma_max), gamma_min_(gamma_min), transfer_(std::move(transfer))
        {
        }

        Vector ChebyshevPreconditioner::operator()(const Vector& b) const
        {

            Spacy::Real cnd = gamma_max_/gamma_min_;
            Spacy::Real theta = (sqrt(cnd)-1)/(sqrt(cnd)+1);
            int N = floor(log(get(acc_))/log(get(theta))) + 1;
           // std::cout << "Cheb: N  = " << N << std::endl;

            Spacy::Real c = 0.5*(gamma_max_-gamma_min_);
            Spacy::Real alpha = 0.5*(gamma_max_+gamma_min_);
            Spacy::Real beta_const = -0.5*c*c/alpha;
            Spacy::Real beta = 0;
            Spacy::Real gamma = -alpha;

            Spacy::Vector dx = P_(b);
            Spacy::Vector x_alt = 0*dx;
            Spacy::Vector x = x_alt;
            Spacy::Vector x_neu = x_alt;

            //auto dxCopy = transfer_(dx,b);
            //auto test = b(dxCopy);

            //if(test < 0)
            //   std::cout << "Fail:  BPX in Chebyshev indefinit: " << test << std::endl;


            for(int i=0; i<N ; i++)
            {

                //std::cout << "i: " << i << std::endl;
                if(i == 1)
                    beta = beta_const;
                if(i >= 2)
                    beta = (1./gamma)*0.25*c*c;
                if(i >= 1)
                    gamma = -(alpha+beta);

                x_neu = -(1./gamma)*(dx+alpha*x+beta*x_alt);

                x_alt = x;
                x = x_neu;

                //auto Ax = A_(x);


                //std::cout <<  "(b-A_(x))((b-A_(x))): " << (b-A_(x))((b-A_(x))) << std::endl;

                if(i < N-1)
                dx = P_(b-A_(x));

               // std::cout <<  "dx" << dx(dx) << std::endl;

               //dxCopy = transfer_(dx,b-A_(x)); //unnecessary call of A
               //dxCopy = transfer_(dx,b);
               //auto test = (b-Ax)(dxCopy);

               //if(test < 0)
               //     std::cout << "Fail:  BPX in Chebyshev indefinit: " << test << std::endl;
            }

            if(verbose())
            {
                std::cout << "Cheb terminating after " << N << " iterations " << std::endl;
            }

            iterations_ = N;
            return x;
        }
    }
}
