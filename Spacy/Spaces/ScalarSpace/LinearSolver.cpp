#include "LinearSolver.h"

#include "Real.h"

#include <Spacy/Util/Cast.h>
#include <Spacy/Vector.h>

#include <utility>

namespace Spacy::Scalar
{
    LinearSolver::LinearSolver( Real y ) : y_( std::move( y ) )
    {
    }

    ::Spacy::Vector LinearSolver::operator()( const ::Spacy::Vector& x ) const
    {
        return cast_ref< Real >( x ) / y_;
    }
} // namespace Spacy::Scalar
