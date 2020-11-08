#include "Norm.h"

#include "Vector.h"
#include "VectorSpace.h"

#include <Spacy/Norm.h>
#include <Spacy/Spaces/ScalarSpace/Real.h>
#include <Spacy/Util/Cast.h>
#include <Spacy/Vector.h>

#include <cmath>

namespace Spacy::ProductSpace
{
    Real Norm::operator()( const ::Spacy::Vector& x ) const
    {
        const auto& x_ = cast_ref< Vector >( x );

        auto result = Real{ 0. };
        for ( auto i = 0u; i < x_.numberOfVariables(); ++i )
        {
            auto nx = norm( x_.component( i ) );
            result += nx * nx;
        }
        return sqrt( result );
    }
} // namespace Spacy::ProductSpace
