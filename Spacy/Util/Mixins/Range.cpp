#include "Range.h"

#include <Spacy/Operator.h>
#include <Spacy/VectorSpace.h>

namespace Spacy::Mixin
{
    Range::Range( const VectorSpace& range ) : range_( range )
    {
    }

    Range& Range::operator=( Range&& other )
    {
        checkSpaceCompatibility( range(), other.range() );
        return *this;
    }

    Range& Range::operator=( const Range& other )
    {
        checkSpaceCompatibility( range(), other.range() );
        return *this;
    }

    const VectorSpace& Range::range() const
    {
        return range_;
    }
} // namespace Spacy::Mixin
