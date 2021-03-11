#include "Vector.h"

#include <Spacy/Spaces/ScalarSpace/Real.h>
#include <Spacy/Util/Cast.h>
#include <Spacy/Util/Exceptions.h>
#include <Spacy/Operator.h>
#include <Spacy/VectorSpace.h>
#include <Spacy/ZeroVectorCreator.h>

#include "VectorSpace.h"

#include <utility>

namespace Spacy
{
    const auto& get_creator( const VectorSpace& space )
    {
        auto&& c = creator< FEniCS::VectorCreator >( space );
        return c.get();
    }

    namespace FEniCS
    {
        Vector::Vector( const VectorSpace& V ) : VectorBase( V ), v_( get_creator( V ) )
        {
        }

        Vector::Vector( const dolfin::Function& v, const VectorSpace& V ) : VectorBase( V ), v_( v )
        {
        }

        Vector& Vector::operator=( const dolfin::Function& v )
        {
            v_ = v;
            return *this;
        }

        //    Vector& Vector::axpy(double a, const Vector& y)
        //    {
        //      get().vector()->axpy(a,*castTo<Vector>(y).get().vector());
        //      return *this;
        //    }

        dolfin::GenericVector& Vector::get()
        {
            return *v_.vector();
        }

        const dolfin::GenericVector& Vector::get() const
        {
            return *v_.vector();
        }

        Real Vector::operator()( const Vector& y ) const
        {
            return get().inner( y.get() );
        }

        ContiguousIterator< double > Vector::begin()
        {
            throw Exception::NotImplemented( __func__, "FEniCS::LinearOperator" );
            return ContiguousIterator< double >();
        }

        ContiguousIterator< double > Vector::end()
        {
            throw Exception::NotImplemented( __func__, "FEniCS::LinearOperator" );
            return ContiguousIterator< double >();
        }

        ContiguousIterator< const double > Vector::begin() const
        {
            throw Exception::NotImplemented( __func__, "FEniCS::LinearOperator" );
            return ContiguousIterator< const double >();
        }

        ContiguousIterator< const double > Vector::end() const
        {
            throw Exception::NotImplemented( __func__, "FEniCS::LinearOperator" );
            return ContiguousIterator< const double >();
        }
    }
}
