#pragma once

#include <string>
#include <dolfin.h>

#include <Spacy/Util/Base/AddArithmeticOperators.h>
#include <Spacy/Util/Base/VectorBase.h>

namespace Spacy
{
    /// @cond
    class Real;
    class VectorSpace;
    /// @endcond

    namespace FEniCS
    {
        /**
         * @ingroup FenicsGroup VectorSpaceGroup
         * @brief %Vector implementation for %FEniCS (single space).
         */
        class Vector : public VectorBase, public AddArithmeticOperators< Vector >
        {
        public:
            using Iterator = VectorIterator< dolfin::GenericVector, dolfin::la_index >;
            using ConstIterator = VectorIterator< const dolfin::GenericVector, dolfin::la_index >;

            /// Construct zero vector \f$x=0\f$ from underlying vector space.
            explicit Vector( const VectorSpace& V );

            Vector( const dolfin::Function& v, const VectorSpace& V );

            /// Assign from dolfin::Function.
            Vector& operator=( const dolfin::Function& v );

            //      Vector& axpy(double a, const AbstractVector& y);

            /// Access reference to underlying coefficient vector.
            dolfin::GenericVector& get();

            /// Access reference to underlying coefficient vector.
            const dolfin::GenericVector& get() const;

            /// Apply as dual element.
            Real operator()( const Vector& y ) const;

            Iterator begin();

            Iterator end();

            ConstIterator begin() const;

            ConstIterator end() const;

        private:
            friend void writeVTK( const Vector&, const std::string& );
            dolfin::Function v_;
        };
    }
}
