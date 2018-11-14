#include "CG.h"

#include <Spacy/Util/Cast.h>
#include <Spacy/Util/Exceptions.h>
#include <Spacy/vector.hh>

#include "RegularizeViaPreconditioner.h"
#include "TerminationCriteria.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <utility>

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

            auto sigma = abs( r( Qr ) ); // preconditioned residual norm squared

            // the conjugate gradient iteration
            for ( unsigned step = 1; true; step++ )
            {
                // if( getVerbosityLevel() > 1 ) std::cout << "Iteration: " << step << std::endl;
                const auto Aq = A_( q );
                Real qAq = Aq( q );

                // in the case of arbitrary Regularization, we have to compute Rq
                if ( is< RegularizeViaCallableOperator >( regularization_ ) )
                {
                    auto& regimpl = cast_ref< RegularizeViaCallableOperator >( regularization_ );
                    Rq = regimpl.getRegularization()( q );
                }
                Real qRq = Rq( q );

                regularization_.apply( qAq, qRq );

                const auto alpha = sigma / qAq;

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
                const auto sigmaNew = abs( r( Qr ) ); // sigma = <Qr,r>
                const auto beta = sigmaNew / sigma;
                sigma = sigmaNew;

                q *= get( beta );
                q += Qr; //  q = Qr + beta*q

                // if regularization is the preconditioner, we update it, as Pq is not accessible
                if ( is< RegularizeViaCallableOperator >( regularization_ ) )
                    continue;

                Rq *= get( beta );
                Rq += r; // Pq = r + beta*Pq
            }

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
            if ( verbose() )
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
