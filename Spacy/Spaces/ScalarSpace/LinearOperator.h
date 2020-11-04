#pragma once

#include <functional>

#include <Spacy/Util/Base/AddArithmeticOperators.h>
#include <Spacy/Util/Base/OperatorBase.h>
#include <Spacy/Util/Base/VectorBase.h>
#include <Spacy/Util/Mixins/Get.h>

namespace Spacy
{
    /// @cond
    class Vector;
    /// @endcond

    /** @addtogroup ScalarGroup  @{ */
    namespace Scalar
    {
        /// Linear operator for scalar problems.
        struct LinearOperator : VectorBase, OperatorBase, Mixin::Get< Real >, AddArithmeticOperators< LinearOperator >
        {
            /// Create linear operator.
            LinearOperator( const VectorSpace& space, Real value );

            /// Apply linear operator.
            ::Spacy::Vector operator()( const ::Spacy::Vector& dx ) const;

            /// Apply as dual space element.
            Real operator()( const LinearOperator& dx ) const;

            /// Get additive inverse.
            LinearOperator operator-() const;

            /// Get solver for the computation of the inverse of this linear operator.
            [[nodiscard]] std::function< ::Spacy::Vector( const ::Spacy::Vector& ) > solver() const;
        };
    } // namespace Scalar
    /** @} */
} // namespace Spacy
