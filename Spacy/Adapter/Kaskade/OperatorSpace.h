#pragma once

#include "Vector.h"

namespace Spacy
{
    /// @cond
    class VectorSpace;
    /// @endcond

    namespace Kaskade
    {
        /// @cond
        template < class, class, class, bool, bool >
        class LinearOperator;
        /// @endcond

        template < class AnsatzVariableSetDescription, class TestVariableSetDescription, bool transposed = false, bool shared = false >
        class LinearOperatorCreator
        {
            using Spaces = typename AnsatzVariableSetDescription::Spaces;
            using Scalar = typename AnsatzVariableSetDescription::Scalar;
            using Domain = typename AnsatzVariableSetDescription::template CoefficientVectorRepresentation<>::type;
            using Range = typename TestVariableSetDescription::template CoefficientVectorRepresentation<>::type;
            using Matrix = ::Kaskade::MatrixAsTriplet< Scalar >;
            using OperatorImpl = ::Kaskade::MatrixRepresentedOperator< Matrix, Domain, Range >;

        public:
            LinearOperatorCreator( const VectorSpace& domain, const VectorSpace& range ) : domain_( domain ), range_( range )
            {
            }

            LinearOperator< AnsatzVariableSetDescription, TestVariableSetDescription, Matrix, transposed, shared >
            operator()( const VectorSpace* space ) const
            {
                return { OperatorImpl{ Matrix{} }, *space };
            }

            [[nodiscard]] const VectorSpace& domain() const
            {
                return domain_;
            }

            [[nodiscard]] const VectorSpace& range() const
            {
                return range_;
            }

        private:
            const VectorSpace &domain_, &range_;
        };
    } // namespace Kaskade
} // namespace Spacy
