#pragma once

#include <Spacy/Util/Mixins/Domain.h>
#include <Spacy/Util/Mixins/Range.h>

namespace Spacy
{
    /// Base class for operators \f$A:\ X\rightarrow Y\f$, between function spaces \f$X\f$ and
    /// \f$Y\f$.
    
    class OperatorBase : public Mixin::Domain, public Mixin::Range
    {
    public:
        /**
         * @brief Constructor.
         * @param domain domain space \f$X\f$.
         * @param range range space \f$Y\f$.
         */
        OperatorBase( const VectorSpace& domain, const VectorSpace& range );

    };
}
