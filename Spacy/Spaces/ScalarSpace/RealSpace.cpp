#include "RealSpace.h"

#include "Real.h"

#include <Spacy/Operator.h>
#include <Spacy/Util/Cast.h>
#include <Spacy/Vector.h>
#include <Spacy/ZeroVectorCreator.h>

namespace Spacy
{
    VectorSpace makeRealSpace( bool defaultIndex )
    {
        return makeHilbertSpace(
            []( const VectorSpace* space ) { return Real{ *space }; },
            []( const Vector& x, const Vector& y ) { return Real( get( cast_ref< Real >( x ) * cast_ref< Real >( y ) ) ); }, "R",
            defaultIndex );
    }
} // namespace Spacy
