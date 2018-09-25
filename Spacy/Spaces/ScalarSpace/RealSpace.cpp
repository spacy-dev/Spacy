#include "RealSpace.h"

#include "Real.h"
#include <Spacy/Util/cast.hh>
#include <Spacy/vector.hh>
#include <Spacy/zeroVectorCreator.hh>

namespace Spacy
{
    VectorSpace make_real_space( bool defaultIndex )
    {
        return makeHilbertSpace( []( const VectorSpace* space ) { return Real{*space}; },
                                 []( const Vector& x, const Vector& y ) {
                                     return Real(
                                         get( cast_ref< Real >( x ) * cast_ref< Real >( y ) ) );
                                 },
                                 defaultIndex );
    }
}
