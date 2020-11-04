#pragma once

#include <functional>

#include <Spacy/Vector.h>
#include <Spacy/VectorSpace.h>

namespace Spacy
{
    /**
     * An operator that does not know about domain and range spaces.
     * Good enough for some algorithms (i.e. for Krylov-methods).
     */
    using Callable = std::function< Vector( const Vector& ) >;

    /// Type-erased operator \f$A:\ X \to Y \f$.
    class Operator
    {
    public:
        /// Apply operator.
        Vector operator()( const Vector& x ) const;

        /// Access domain space \f$X\f$.
        const VectorSpace& domain() const;

        /// Access range space \f$Y\f$.
        const VectorSpace& range() const;
    };
} // namespace Spacy
