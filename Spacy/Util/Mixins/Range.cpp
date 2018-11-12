#include "Range.h"

#include <Spacy/vectorSpace.hh>

namespace Spacy
{
    namespace Mixin
    {
        Range::Range(const VectorSpace& range)
            : range_(range)
        {}

        Range& Range::operator=(Range&& other)
        {
            checkSpaceCompatibility(range(), other.range());
            return *this;
        }

        Range& Range::operator=(const Range& other)
        {
            checkSpaceCompatibility(range(), other.range());
            return *this;
        }

        const VectorSpace& Range::range() const
        {
            return range_;
        }
    }
}
