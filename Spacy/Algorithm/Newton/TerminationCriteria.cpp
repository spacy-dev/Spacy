#include "TerminationCriteria.h"

#include <Spacy/Algorithm/DampingFactor.h>
#include <Spacy/C1Operator.h>
#include <Spacy/Operator.h>
#include <Spacy/Vector.h>

namespace Spacy::Newton::Termination
{
    AffineCovariant::AffineCovariant( const C1Operator& /*unused*/, Real relativeAccuracy, Real eps )
        : Mixin::RelativeAccuracy( relativeAccuracy ), Mixin::Eps( eps )
    {
    }

    bool AffineCovariant::operator()( DampingFactor nu, const Vector& x, const Vector& dx ) const
    {
        if ( beforeFirstIteration_ )
        {
            beforeFirstIteration_ = false;
            return false;
        }
        if ( abs( get( nu - 1 ) ) > eps() )
            return false;

        const auto norm_x = norm( x );
        const auto norm_dx = norm( dx );
        if ( norm_x == 0 && norm_dx == 0 )
            return true;

        if ( norm_dx < getRelativeAccuracy() * norm_x )
        {
            LOG( covariant_log_tag, "Terminating -- relative accuracy: ", norm_dx / norm_x );
            return true;
        }

        return false;
    }

    AffineContravariant::AffineContravariant( const C1Operator& F, Real relativeAccuracy, Real eps )
        : Mixin::RelativeAccuracy( relativeAccuracy ), Mixin::Eps( eps ), F_( F )
    {
    }

    bool AffineContravariant::operator()( DampingFactor nu, const Vector& x, const Vector& /*unused*/ ) const
    {
        if ( initialResidual < 0 )
        {
            // if( initialResidual = norm(F_(x)) )
            initialResidual = norm( F_( x ) );
            if ( initialResidual == 0 )
                return true;

            LOG( contravariant_log_tag, "Initial residual:", initialResidual );
            return false;
        }

        if ( abs( get( nu - 1 ) ) > eps() )
            return false;

        if ( norm( F_( x ) ) < getRelativeAccuracy() * initialResidual )
        {
            LOG( contravariant_log_tag, "Terminating -- relative accuracy: ", norm( F_( x ) ) / initialResidual );
            return true;
        }

        return false;
    }
} // namespace Spacy::Newton::Termination
