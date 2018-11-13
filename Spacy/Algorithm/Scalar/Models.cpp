#include "Models.h"

namespace Spacy
{
    namespace Scalar
    {
        Quadratic::Quadratic( Real a, Real b, Real c ) noexcept : a_( a ), b_( b ), c_( c )
        {
        }

        Real Quadratic::operator()( Real t ) const noexcept
        {
            return a_ + t * b_ + t * t * c_;
        }

        Cubic::Cubic( Real a, Real b, Real c, Real d ) noexcept : a_( a ), b_( b ), c_( c ), d_( d )
        {
        }

        Real Cubic::operator()( Real t ) const noexcept
        {
            return a_ + t * ( b_ + t * ( c_ + d_ * t ) );
        }
    }
}
