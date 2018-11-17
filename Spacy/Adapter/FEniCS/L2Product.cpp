#include "L2Product.h"

#include <Spacy/Spaces/ScalarSpace/Real.h>
#include <Spacy/Util/Cast.h>
#include <Spacy/Vector.h>
#include <Spacy/VectorSpace.h>

#include "Vector.h"

namespace Spacy
{
    namespace FEniCS
    {
        Real l2Product::operator()( const ::Spacy::Vector& x, const ::Spacy::Vector& y ) const
        {
            checkSpaceCompatibility( x.space(), y.space() );
            return cast_ref< Vector >( x ).get().inner( cast_ref< Vector >( y ).get() );
        }
    }
}
