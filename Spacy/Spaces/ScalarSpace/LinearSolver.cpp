#include "LinearSolver.h"

#include <Spacy/Spaces/ScalarSpace/Real.h>
#include <Spacy/Util/cast.hh>
#include <Spacy/vector.hh>

namespace Spacy
{
    namespace Scalar
    {
        LinearSolver::LinearSolver( Real y ) : y_( y )
        {
        }

        ::Spacy::Vector LinearSolver::operator()( const ::Spacy::Vector& x ) const
        {
            return cast_ref< Real >( x ) / y_;
        }
    }
}
