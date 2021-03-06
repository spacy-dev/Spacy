#pragma once

#include <functional>

#include <Spacy/Util/Base/FunctionalBase.h>

#include <memory>

namespace Spacy
{
    /// @cond
    class Real;
    class Vector;
    /// @endcond

    /** @addtogroup ScalarGroup @{ */
    namespace Scalar
    {
        /// @cond
        class LinearOperator;
        /// @endcond

        class C2Functional : public ::Spacy::FunctionalBase
        {
        public:
            C2Functional( std::function< double( double ) > value, std::function< double( double ) > derivative,
                          std::function< double( double ) > secDerivative );

            // Compute f(x).
            ::Spacy::Real operator()( const ::Spacy::Vector& x ) const;

            // Compute f'(x) as element of X*.
            [[nodiscard]] ::Spacy::Vector d1( const ::Spacy::Vector& x ) const;

            // Compute f''(x)dx as element of X*.
            [[nodiscard]] ::Spacy::Vector d2( const ::Spacy::Vector& x, const ::Spacy::Vector& dx ) const;

            // Access f''(x) as mapping f''(x): X->X*.
            [[nodiscard]] ::Spacy::Scalar::LinearOperator hessian( const ::Spacy::Vector& x ) const;

        private:
            std::function< double( const double& ) > value_;
            std::function< double( const double& ) > derivative_;
            std::function< double( const double& ) > secDerivative_;
            std::shared_ptr< VectorSpace > operatorSpace_;
        };
    } // namespace Scalar
    /** @} */
} // namespace Spacy
