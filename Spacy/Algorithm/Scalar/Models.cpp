#include "Models.h"

#include <utility>

namespace Spacy
{
    namespace Scalar
    {
        Quadratic::Quadratic( Real a, Real b, Real c ) noexcept : a_( std::move( a ) ), b_( std::move( b ) ), c_( std::move( c ) )
        {
        }

        Real Quadratic::operator()( Real t ) const noexcept
        {
            return a_ + t * b_ + t * t * c_;
        }

        Cubic::Cubic( Real a, Real b, Real c, Real d ) noexcept
            : a_( std::move( a ) ), b_( std::move( b ) ), c_( std::move( c ) ), d_( std::move( d ) )
        {
        }

        Real Cubic::operator()( Real t ) const noexcept
        {
            return a_ + t * ( b_ + t * ( c_ + d_ * t ) );
        }
    } // namespace Scalar
} // namespace Spacy
