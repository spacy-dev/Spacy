#pragma once

#include <Spacy/Spaces/ScalarSpace/Real.h>

namespace Spacy::Scalar
{
    /// A one-dimensional quadratic function \f$q(t) = a + bt + ct^2\f$.
    class Quadratic
    {
    public:
        Quadratic( Real a, Real b, Real c ) noexcept;

        /// Compute \f$q(t) = a + bt + ct^2 \f$.
        Real operator()( Real t ) const noexcept;

    private:
        Real a_, b_, c_;
    };

    /// A one-dimensional cubic function \f$q(t) = a + bt + ct^2 + dt^3\f$.
    class Cubic
    {
    public:
        Cubic( Real a, Real b, Real c, Real d ) noexcept;

        /// Compute \f$q(t) = a + bt + ct^2 + dt^3 \f$.
        Real operator()( Real t ) const noexcept;

    private:
        Real a_, b_, c_, d_;
    };
} // namespace Spacy::Scalar
