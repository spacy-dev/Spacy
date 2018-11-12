#include "Domain.h"

#include <Spacy/vectorSpace.hh>

namespace Spacy
{
    namespace Mixin
    {
        Domain::Domain(const VectorSpace& domain)
            : domain_(domain)
        {}

        Domain& Domain::operator=(Domain&& other)
        {
            checkSpaceCompatibility(domain(),other.domain());
            return *this;
        }

        Domain& Domain::operator=(const Domain& other)
        {
            checkSpaceCompatibility(domain(),other.domain());
            return *this;
        }
        const VectorSpace& Domain::domain() const
        {
            return domain_;
        }
    }
}
