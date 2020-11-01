#include "linearSolver.hh"

#include <Spacy/Vector.h>

namespace Mock
{
    Spacy::Vector IndefiniteLinearSolver::operator()( const Spacy::Vector& x ) const
    {
        return x;
    }
} // namespace Mock
