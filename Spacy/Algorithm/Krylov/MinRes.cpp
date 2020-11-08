#include "MinRes.h"

#include <Spacy/Util/Log.h>
#include <Spacy/Util/Logger.h>
#include <Spacy/ZeroVectorCreator.h>

namespace Spacy::Algorithm
{
    DEFINE_LOG_TAG( static const char* log_tag = "MINRES" );

    MinRes::MinRes( Operator H, CallableOperator P ) : H_{ std::move( H ) }, P_{ std::move( P ) }
    {
    }

    Vector MinRes::operator()( const Vector& b ) const
    {
        return minresLoop( zero( H_.domain() ), b );
    }

    int MinRes::getIterations() const noexcept
    {
        return iterations_;
    }

    bool MinRes::convergenceTest( Real norm_r0, Real norm_r, const Vector& x ) const
    {
        if ( norm_r / norm_r0 < getRelativeAccuracy() || norm_r < eps() )
        {
            LOG( log_tag, "Norm Solution: ", sqrt( x( x ) ) );
            return true;
        }

        return false;
    }

    Vector MinRes::minresLoop( Vector x, const Vector& b ) const
    {
        LOG_INFO( log_tag, "Starting MINRES-Solver." );
        LOG_SEPARATOR( log_tag );

        iterations_ = 0;

        auto v_old = zero( b.space() );
        LOG( log_tag, "Iteration: ", 1 );

        auto v = b - H_( x );

        auto v_new = v_old;
        auto w = zero( x.space() );
        auto w_old = w;
        auto w_new = w;

        auto z = P_( v );
        iterations_++;

        auto z_new = z;
        auto gamma = sqrt( v( z ) );

        z *= 1.0 / get( gamma );
        v *= 1.0 / get( gamma );

        auto eta_old = gamma;
        Real s_old = 0.0;
        Real s = 0.0;
        Real c_old = 1.0;
        Real c = 1.0;
        const auto norm_r0 = norm( eta_old );

        for ( unsigned step = 1; step < getMaxSteps(); step++ )
        {
            LOG_SEPARATOR( log_tag );

            auto Hzmv = H_( z ) - gamma * v_old;
            auto delta = z( Hzmv );

            v_new = Hzmv - delta * v;

            z_new = P_( v_new );
            iterations_++;
            const auto gamma_new = sqrt( v_new( z_new ) );

            z_new *= 1.0 / get( gamma_new );
            v_new *= 1.0 / get( gamma_new );

            const auto alpha_0 = c * delta - c_old * s * gamma;
            const auto alpha_1 = sqrt( alpha_0 * alpha_0 + gamma_new * gamma_new );
            const auto alpha_2 = ( s * delta + c_old * c * gamma );
            const auto alpha_3 = s_old * gamma;

            const auto c_new = alpha_0 / alpha_1;
            const auto s_new = gamma_new / alpha_1;

            w_new = ( 1.0 / alpha_1 ) * ( z - alpha_3 * w_old - alpha_2 * w );
            x = x + c_new * eta_old * w_new;
            eta_old = -s_new * eta_old;
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

            const auto norm_r = norm( eta_old );
            LOG( log_tag, "Norm residuum: ", norm_r );
            if ( convergenceTest( norm_r0, norm_r, x ) )
            {
                LOG( log_tag, "MINRES converged in step: ", iterations_ );
                return x;
            }
        }

        LOG_INFO( log_tag, "Number of MINRES iterations > Maximum Number of Iterations" );
        return x;
    }
} // namespace Spacy::Algorithm
