#include "OperatorBase.h"

namespace Spacy
{
    OperatorBase::OperatorBase( const VectorSpace& domain, const VectorSpace& range )
        : Mixin::Domain( domain ), Mixin::Range( range )
    {
    }
}
