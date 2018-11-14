#include "LinearSolver.h"

#include <Spacy/Util/Cast.h>
#include <Spacy/vector.hh>

#include "Real.h"

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
