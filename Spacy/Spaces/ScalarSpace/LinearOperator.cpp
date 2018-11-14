#include "LinearOperator.h"

#include <Spacy/Spaces/RealSpace.h>
#include <Spacy/Util/Cast.h>
#include <Spacy/linearSolver.hh>
#include <Spacy/vector.hh>

#include "LinearSolver.h"

namespace Spacy
{
    namespace Scalar
    {
        LinearOperator::LinearOperator( const VectorSpace& space, Real value )
            : VectorBase( space ), OperatorBase( Space::R, Space::R ), Mixin::Get< Real >( value )
        {
        }

        ::Spacy::Vector LinearOperator::operator()( const ::Spacy::Vector& dx ) const
        {
            return Real( get() * Mixin::get( cast_ref< Real >( dx ) ) );
        }

        ::Spacy::Real LinearOperator::operator()( const LinearOperator& dx ) const
        {
            return Real( get() * dx.get() );
        }

        ::Spacy::LinearSolver LinearOperator::solver() const
        {
            return LinearSolver( get() );
        }

        LinearOperator LinearOperator::operator-() const
        {
            return LinearOperator( space(), -get() );
        }
    }
}
