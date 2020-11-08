#include "ScalarProduct.h"

#include <Spacy/Spaces/ScalarSpace/Real.h>
#include <Spacy/Vector.h>

namespace Spacy::Generic
{
    Spacy::Real EuclideanScalarProduct::operator()( const ::Spacy::Vector& x, const ::Spacy::Vector& y ) const
    {
        checkSpaceCompatibility( x.space(), y.space() );
        return x( y );
    }
} // namespace Spacy::Generic
