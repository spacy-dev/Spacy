#include "LinearOperator.h"

#include <Spacy/Spaces/RealSpace.h>
#include <Spacy/Util/Exceptions.h>
#include <Spacy/Operator.h>
#include <Spacy/Vector.h>
#include <Spacy/VectorSpace.h>
#include <Spacy/ZeroVectorCreator.h>

#include "Copy.h"
#include "OperatorSpace.h"
#include "VectorSpace.h"

#include <cassert>

namespace Spacy
{
    namespace FEniCS
    {
        LinearOperator::LinearOperator(
            std::shared_ptr< dolfin::GenericMatrix > A, const VectorSpace& space,
            std::shared_ptr< const dolfin::FunctionSpace > dolfinSpace,
            std::function< LinearSolver( const LinearOperator& ) > solverCreator )
            : OperatorBase( creator< LinearOperatorCreator >( space ).domain(),
                            creator< LinearOperatorCreator >( space ).range() ),
              VectorBase( space ), A_( A ), solverCreator_( std::move( solverCreator ) ),
              space_( dolfinSpace )
        {
        }

        ::Spacy::Vector LinearOperator::operator()( const ::Spacy::Vector& x ) const
        {
            if ( x.space().index() == domain().index() )
                return applyAsOperator( x );

            if ( x.space().isPrimalWRT( space() ) )
                return applyAsDualElement( x );

            assert( false );
        }

        Real LinearOperator::operator()( const LinearOperator& A ) const
        {
            return applyAsDualElement( A );
        }

        LinearSolver LinearOperator::solver() const
        {
            return solverCreator_( *this );
        }

        dolfin::GenericMatrix& LinearOperator::get()
        {
            assert( A_ != nullptr );
            return *A_;
        }

        const dolfin::GenericMatrix& LinearOperator::get() const
        {
            assert( A_ != nullptr );
            return *A_;
        }

        bool LinearOperator::symmetric() const
        {
            return symmetric_;
        }

        std::shared_ptr< const dolfin::FunctionSpace > LinearOperator::dolfinSpace() const
        {
            return space_;
        }

        ::Spacy::Vector LinearOperator::applyAsOperator( const ::Spacy::Vector& x ) const
        {
            auto x_ = dolfin::Function( space_ );
            copy( x, x_ );
            auto y_ = x_.vector()->copy();
            get().mult( *x_.vector(), *y_ );

            auto y = zero( range() );
            copy( *y_, y );

            return y;
        }

        ::Spacy::Real LinearOperator::applyAsDualElement( const ::Spacy::Vector& x ) const
        {
            throw Exception::NotImplemented( __func__, "FEniCS::LinearOperator" );
            return ::Spacy::Real( 0 );
        }

        ContiguousIterator< double > LinearOperator::begin()
        {
            throw Exception::NotImplemented( __func__, "FEniCS::LinearOperator" );
            return ContiguousIterator< double >();
        }

        ContiguousIterator< double > LinearOperator::end()
        {
            throw Exception::NotImplemented( __func__, "FEniCS::LinearOperator" );
            return ContiguousIterator< double >();
        }

        ContiguousIterator< const double > LinearOperator::begin() const
        {
            throw Exception::NotImplemented( __func__, "FEniCS::LinearOperator" );
            return ContiguousIterator< const double >();
        }

        ContiguousIterator< const double > LinearOperator::end() const
        {
            throw Exception::NotImplemented( __func__, "FEniCS::LinearOperator" );
            return ContiguousIterator< const double >();
        }
    }
}
