#pragma once

#include "Real.h"

#include <Spacy/Util/Base/OperatorBase.h>

#include <functional>

namespace Spacy
{
    /// @cond
    class Vector;
    /// @endcond

    /** @addtogroup ScalarGroup @{ */
    namespace Scalar
    {
        /// Operator \f$A:\mathbb{R}\rightarrow\mathbb{R}\f$.
        class DynamicOperator : public Spacy::OperatorBase
        {
        public:
            /// Constructor for operators between the global real space (index = 0).
            explicit DynamicOperator( std::function< double( double, double ) > value );

            /// Constructor for general operators between real spaces.
            DynamicOperator( std::function< double( double, double ) > value, const VectorSpace& domain,
                      const VectorSpace& range );

            /// Compute \f$A(x)\in\mathbb{R}\f$.
            Spacy::Vector operator()( double t,  const ::Spacy::Vector& x ) const;

        private:
            std::function< double( double, double ) > value_;
        };
    }
    /** @} */
}
