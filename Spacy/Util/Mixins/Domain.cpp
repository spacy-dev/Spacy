#include "Domain.h"

#include <Spacy/VectorSpace.h>

namespace Spacy::Mixin
{
    Domain::Domain( const VectorSpace& domain ) : domain_( domain )
    {
    }

    Domain& Domain::operator=( Domain&& other )
    {
        checkSpaceCompatibility( domain(), other.domain() );
        return *this;
    }

    Domain& Domain::operator=( const Domain& other )
    {
        checkSpaceCompatibility( domain(), other.domain() );
        return *this;
    }
    const VectorSpace& Domain::domain() const
    {
        return domain_;
    }
} // namespace Spacy::Mixin
