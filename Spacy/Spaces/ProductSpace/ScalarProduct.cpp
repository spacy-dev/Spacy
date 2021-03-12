#include "ScalarProduct.h"

#include "Vector.h"
#include "VectorSpace.h"

#include <Spacy/Operator.h>
#include <Spacy/Spaces/ScalarSpace/Real.h>
#include <Spacy/Spaces/ScalarSpace/RealSpace.h>
#include <Spacy/Util/Cast.h>
#include <Spacy/Vector.h>
#include <Spacy/VectorSpace.h>

namespace Spacy::ProductSpace
{
    Real ScalarProduct::operator()( const ::Spacy::Vector& x, const ::Spacy::Vector& y ) const
    {
        checkSpaceCompatibility( x.space(), y.space() );
        const auto& x_ = cast_ref< Vector >( x );
        const auto& y_ = cast_ref< Vector >( y );

        auto result = Real{ 0. };
        for ( auto i = 0u; i < x_.numberOfVariables(); ++i )
            result += x_.component( i ) * y_.component( i );
        return result;
    }
} // namespace Spacy::ProductSpace
