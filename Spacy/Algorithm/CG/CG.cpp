#include "CG.h"

#include <Spacy/Util/Cast.h>
#include <Spacy/Util/Exceptions.h>
#include <Spacy/Vector.h>

#include "RegularizeViaPreconditioner.h"
#include "TerminationCriteria.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <utility>
//#include "/home/bt304332/svn-repos/Kaskade7.3/linalg/lapack/lapacke.h"

#define LAPACK_COL_MAJOR 102

namespace Spacy
{
    namespace CG
    {
        Solver::Solver( CallableOperator A, CallableOperator P, Regularization regularization,
                        bool truncated )
            : A_( std::move( A ) ), P_( std::move( P ) ),
              terminate_( CG::Termination::StrakosTichyEnergyError{} ), truncated_( truncated ),
              regularization_( std::move( regularization ) )
        {
            if ( !::Spacy::is< NoRegularization >( regularization_ ) )
                regularized_ = true;
        }

        Vector Solver::solve( const Vector& x, const Vector& b ) const
        {
            regularization_.init();
            definiteness_ = DefiniteNess::PositiveDefinite;
            result = Result::Failed;

            terminate_.set_eps( get( eps() ) );
            terminate_.setRelativeAccuracy( get( getRelativeAccuracy() ) );
            terminate_.setAbsoluteAccuracy( get( getAbsoluteAccuracy() ) );

            if ( !regularization_ )
                return cgLoop( x, b );
            else
            {
                Vector y = x;
                while ( result != Result::Converged && result != Result::TruncatedAtNonConvexity )
                    y = cgLoop( x, b );
                return y;
            }
        }

	Vector Solver::operator()( const Vector& b ) const
        {
            auto x = 0 * P_(b);
            return solve(x,b);
        }

        CG::TerminationCriterion& Solver::terminationCriterion() noexcept
        {
            return terminate_;
        }

        bool Solver::indefiniteOperator() const noexcept
        {
            return definiteness_ == DefiniteNess::Indefinite;
        }

        const CallableOperator& Solver::P() const
        {
            return P_;
        }

        const CallableOperator& Solver::A() const
        {
            return A_;
        }

        const Regularization& Solver::R() const
        {
            return regularization_;
        }

        Vector Solver::cgLoop( Vector x, Vector r ) const
        {
            std::vector<Spacy::Real> alpha_vec;
            std::vector<Spacy::Real> beta_vec;

            terminate_.clear();
            result = Result::Failed;

            if ( is< RegularizeViaCallableOperator >( regularization_ ) )
            {
                auto& regimpl = cast_ref< RegularizeViaCallableOperator >( regularization_ );

                regimpl.resetMemory();
            }
            // initialization phase for conjugate gradients
            auto Ax = A_( x );
            r -= Ax;
            auto Qr = Q( r );

            auto q = Qr;
            auto Rq = r; // required only for regularized or hybrid conjugate gradient methods with
                         // precond as regularization
            auto Qr_copy = transfer(Qr,r);
            auto sigma = abs( r( Qr_copy ) );

            if(r( Qr_copy ) < 0.)
            {
                definiteness_ = DefiniteNess::Indefinite;
                result = Result::TruncatedAtNonConvexity;
                std::cout << "Fail: Preconditioner Indefinit" << std::endl;
                std::cout << r( Qr_copy ) << std::endl;
                return x;
            }
            //auto sigma = abs( r( Qr ) ); // preconditioned residual norm squared


            auto qLast = q;
            // the conjugate gradient iteration
            for ( unsigned step = 1; true; step++ )
            {
                // if( getVerbosityLevel() > 1 ) std::cout << "Iteration: " << step << std::endl;
                const auto Aq = A_( q );
                auto Aq_copy = transfer(Aq,q);
                Real qAq = Aq_copy( q );


              //  std::cout << "Cg qAq: " << qAq << std::endl << std::endl;

               // std::cout << "qLast(Aq):  " << Aq_copy(qLast) << std::endl;
                qLast = q;
                //Real qAq = Aq( q );

                // in the case of arbitrary Regularization, we have to compute Rq
                if ( is< RegularizeViaCallableOperator >( regularization_ ) )
                {
                    auto& regimpl = cast_ref< RegularizeViaCallableOperator >( regularization_ );
                    Rq = regimpl.getRegularization()( q );
                }
                auto Rq_copy = transfer(Rq,q);
                Real qRq = Rq_copy( q );
                //Real qRq = Rq( q );

                regularization_.apply( qAq, qRq );

                const auto alpha = sigma / qAq;
                alpha_vec.push_back(alpha);

                terminate_.update( get( alpha ), get( qAq ), get( qRq ), get( sigma ), x );
                //  don't trust small numbers
                if ( vanishingStep( step ) )
                {
                    if ( step == 1 )
                        x += q;
                    result = Result::Converged;
                    iterations_ = step;
                    break;
                }

              //  std::cout << "Cg alpha " << alpha << std::endl;

                if ( terminateOnNonconvexity( qAq, qRq, x, q, step ) )
                {
                    iterations_ = step;
                    break;
                }
                x += ( get( alpha ) * q );

                // convergence test
                if ( terminate_() )
                {
                    if ( verbose() )
                        std::cout << "    "
                                  << ": Terminating in iteration " << step << ".\n";
                    result = ( step == getMaxSteps() ) ? Result::Failed : Result::Converged;
                    iterations_ = step;
                    break;
                }

                r -= alpha * Aq;

                regularization_.adjustResidual( alpha, Rq, r );

                Qr = Q( r );

                // determine new search direction
                Qr_copy = transfer(Qr,r);
                const auto sigmaNew = abs( r( Qr_copy ) );
                if(r( Qr_copy ) < 0.)
                {
                    std::cout << "Fail: Preconditioner Indefinit" << std::endl;
                    std::cout << r( Qr_copy )  << std::endl;
                   // result = Result::Failed;
                    definiteness_ = DefiniteNess::Indefinite;
                    result = Result::TruncatedAtNonConvexity;
                    return x;
                }

                //const auto sigmaNew = abs( r( Qr ) ); // sigma = <Qr,r>
                const auto beta = sigmaNew / sigma;
                beta_vec.push_back(beta);
                sigma = sigmaNew;

            //    std::cout << "sigmaNew: " << sigmaNew << std::endl << std::endl;
                q *= get( beta );
                q += Qr; //  q = Qr + beta*q

                // if regularization is the preconditioner, we update it, as Pq is not accessible
                if ( is< RegularizeViaCallableOperator >( regularization_ ) )
                    continue;

                Rq *= get( beta );
                Rq += r; // Pq = r + beta*Pq
            }

            callback(alpha_vec,beta_vec);


//            auto a_size = alpha_vec.size();
//            double * eig = new double[a_size];
//            double * diagonal = new double[a_size];
//            double * subDiagonal = new double[a_size-1];
//            double * z_vec = new double[a_size];
//            int * isuppz_vec = new int[2*a_size];
//            int * m = new int[a_size];

//            for(auto i=0u; i<a_size; ++i)
//            {
//                eig[i] = 0.0;
//                if(i==0)
//                    diagonal[i] = 1./alpha_vec[i].get();
//                else
//                    diagonal[i] = 1./alpha_vec[i].get() + beta_vec[i-1].get()/alpha_vec[i-1].get();
//                if(i!=0)
//                    subDiagonal[i-1] = sqrt(beta_vec[i-1].get()).get()/alpha_vec[i-1].get();
//                z_vec[i] = 0.0;
//                isuppz_vec[i] = 0;
//                isuppz_vec[i+a_size] = 0;
//                m[i] = 0;
//            }

//            LAPACKE_dstevr(LAPACK_COL_MAJOR,'N','A',a_size,diagonal,subDiagonal,0.0,0.0,0,0,1e-8,m,eig,z_vec,a_size,isuppz_vec);

//            std::cout << "EigMin: " <<  eig[0] << "   EigMax: " << eig[a_size-1] << std::endl;
        //    delete[] eig;
         //   delete[] diagonal;
         //   delete[] subDiagonal;
         //   delete[] z_vec;
         //   delete[] isuppz_vec;
         //   delete[] m;

            return x;
        }

        Vector Solver::Q( const Vector& r ) const
        {
            auto Qr = P_( r );
            for ( auto i = 0u; i < getIterativeRefinements(); ++i )
                Qr += P_( r - A_( Qr ) );
            return Qr;
        }

        bool Solver::vanishingStep( unsigned step ) const
        {
            if ( terminate_.vanishingStep() )
            {
                if ( verbose() )
                    std::cout
                        << "    "
                        << ": Terminating due to numerically almost vanishing step in iteration "
                        << step << "." << std::endl;
                result = Result::Converged;
                return true;
            }
            return false;
        }

        bool Solver::terminateOnNonconvexity( Real qAq, Real qRq, Vector& x, const Vector& q,
                                              unsigned step ) const
        {
            if ( qAq > 0 )
                return false;
           // if ( verbose() )
                std::cout << "    "
                          << ": Negative curvature: " << qAq << std::endl;

            if ( !truncated_ && !regularized_ )
            {
                if ( verbose() )
                {
                    std::cout << "    "
                              << ": Direction of negative curvature encountered in standard CG "
                                 "Implementation!"
                              << std::endl;
                    std::cout << "    "
                              << ": Either something is wrong with your operator or you should use "
                                 "TCG, RCG or TRCG. Terminating CG!"
                              << std::endl;
                }

                throw Exception::SingularOperator( "CG::terminateOnNonconvexity" );
            }

            if ( truncated_ && ( !regularized_ || terminate_.minimalDecreaseAchieved() ) )
            {
                // At least do something to retain a little chance to get out of the nonconvexity.
                // If a nonconvexity is encountered in the first step something probably went wrong
                // elsewhere. Chances that a way out of the nonconvexity can be found are small in
                // this case.
                if ( step == 1 )
                    x += q;
                if ( verbose() )
                    std::cout << "    "
                              << ": Truncating at nonconvexity in iteration " << step << ": " << qAq
                              << std::endl;
                definiteness_ = DefiniteNess::Indefinite;
                result = Result::TruncatedAtNonConvexity;
                return true;
            }

            regularization_.update( qAq, qRq );

            if ( verbose() )
                std::cout << "    "
                          << ": Regularizing at nonconvexity in iteration " << step << "."
                          << std::endl;
            definiteness_ = DefiniteNess::Indefinite;
            result = Result::EncounteredNonConvexity;
            return true;
        }
    }
}
