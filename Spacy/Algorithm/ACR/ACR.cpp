#include "ACR.h"

#include <Spacy/Algorithm/CG/CG.h>
#include <Spacy/Algorithm/CG/LinearSolver.h>
#include <Spacy/Algorithm/CompositeStep/QuadraticModel.h>
#include <Spacy/Algorithm/DampingFactor.h>
#include <Spacy/Algorithm/Scalar/FindGlobalMinimizer.h>
#include <Spacy/Operator.h>
#include <Spacy/Util/Cast.h>
#include <Spacy/Util/Log.h>
#include <Spacy/ZeroVectorCreator.h>

#include <cmath>
#include <utility>

namespace
{
    class TrivialPreconditioner
    {
    public:
        // compute A(x)
        const Spacy::Vector& operator()( const Spacy::Vector& x ) const
        {
            return x;
        }
    };
}

namespace Spacy
{
    namespace ACR
    {
        DEFINE_LOG_TAG( static const char* log_tag = "ACR" )

        CompositeStep::CubicModel makeCubicModel( const Vector& dx, const C2Functional& f,
                                                  const Vector& x, Spacy::Real omega )
        {
            return CompositeStep::CubicModel(
                CompositeStep::makeQuadraticModel( Spacy::DampingFactor( 0.0 ), dx, dx, f, x ),
                CompositeStep::makeQuadraticNormModel( Spacy::DampingFactor( 0.0 ), dx, dx ),
                omega );
        }

        ACRSolver::ACRSolver( C2Functional f, double eta1, double eta2, double epsilon,
                              double relativeAccuracy, double omegaMax, double lambdaMax )
            : Mixin::RelativeAccuracy( relativeAccuracy ), f_( std::move( f ) ),
              domain_( f_.domain() ), eta1_( eta1 ), eta2_( eta2 ), epsilon_( epsilon ),
              omegaMax_( omegaMax ), lambdaMax_( lambdaMax )
        {
        }

        Vector ACRSolver::operator()()
        {
            return operator()( zero( domain_ ) );
        }

        Vector ACRSolver::operator()( const Vector& x0 )
        {
            LOG_SEPARATOR( log_tag );

            auto x = x0;

            for ( unsigned step = 1; step <= getMaxSteps(); ++step )
            {

                LOG( log_tag, "Iteration", step )

                stepMonitor = StepMonitor::Accepted;
                // TODO domain_.setScalarProduct( );
                Real lambda = 1.0;

                do
                {
                    auto dx = computeStep( x );
                    LOG( log_tag, "|dx|", norm( dx ) );
                    auto cubicModel = makeCubicModel( dx, f_, x, omega_ );
                    LOG_INFO( log_tag, "Computing damping factor" )

                    lambda = Scalar::findGlobalMinimizer( cubicModel, 0, 1, 1e-2 );

                    dx = lambda * dx;

                    LOG( log_tag, "lambda: ", lambda, " omega: ", omega_, " cubicModel: ",
                         cubicModel( lambda ) )

                    if ( ( stepMonitor = acceptanceTest( x, dx, lambda, cubicModel ) ) ==
                         StepMonitor::Accepted )
                    {
                        x += dx;
                        //  if( convergenceTest(/*TODO*/) ) return x;
                        TrivialPreconditioner trivialP;
                        // Todo: Change Norm
                        if ( f_.d1( x )( trivialP( f_.d1( x ) ) ) < epsilon_ * epsilon_ )
                        {
                            LOG_INFO( log_tag, "Converged" )
                            return x;
                        }
                    }

                    else
                        LOG_INFO( log_tag, "Rejected" )

                    omega_ = weightChange( omega_ );

                    if ( omega_ > omegaMax_ )
                        return x;
                }
                // end iteration
                while ( stepMonitor == StepMonitor::Rejected || lambda < lambdaMax_ );
                LOG_SEPARATOR( log_tag );
            }

            LOG_INFO( log_tag, "Reached maximum number of steps" );
            return x;
        }

        ACRSolver::StepMonitor
        ACRSolver::acceptanceTest( const Vector& x, const Vector& dx, const Real& lambda,
                                   const CompositeStep::CubicModel& cubicModel )
        {
            auto diff = cubicModel( lambda ) - cubicModel( 0.0 );

            // Prevent dividing by zero e.g. in the case x_0 = optimal solution
            // Consider the Case numerator = nominator = 0
            if ( abs( diff ) < eps() )
            {
                if ( diff < 0 )
                    diff = -eps();

                else
                    diff = eps();
            }

            rho_ = ( f_( x + dx ) - f_( x ) ) / diff;

            LOG( log_tag, "f_(x+dx): ", f_( x + dx ), "f_(x): ", f_( x ), "CubicModel(lambda): ",
                 cubicModel( lambda ), "CubicModel(0): ", cubicModel( 0 ), "rho: ", rho_, "eta1 ",
                 eta1_, "eta2 ", eta2_ )

            if ( rho_ >= eta1_ )
                return StepMonitor::Accepted;

            return StepMonitor::Rejected;
        }

        Spacy::Real ACRSolver::weightChange( Spacy::Real omega ) const
        {

            LOG( log_tag, "rho: ", rho_, "eta1 ", eta1_, "eta2 ", eta2_ )

            if ( rho_ > eta2_ )
                omega *= 0.5;
            else if ( rho_ < eta1_ )
                omega *= 2;

            return omega;
        }

        Vector ACRSolver::computeStep( const Spacy::Vector& x ) const
        {
            auto tcg = makeTCGSolver( f_.hessian( x ), []( const Vector& x ) { return x; } );
            tcg.setRelativeAccuracy( getRelativeAccuracy() );

            LOG( log_tag, "rhs", norm( f_.d1( x ) ) );

            return tcg( -f_.d1( x ) );
        }
    }
}
