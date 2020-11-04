#include "linearOperator.hh"

#include "linearSolver.hh"

#include <Spacy/Spaces/RealSpace.h>

namespace Mock
{
    LinearOperator::LinearOperator( double value )
        : Spacy::Mixin::Get< double >( value ), Spacy::OperatorBase( Spacy::Space::R, Spacy::Space::R ), Spacy::VectorBase(
                                                                                                             Spacy::Space::R )
    {
    }

    ::Spacy::Vector LinearOperator::operator()( const ::Spacy::Vector& dx ) const
    {
        return Spacy::Real( get() );
    }

    ::Spacy::Real LinearOperator::operator()( const LinearOperator& y ) const
    {
        return { get() * y.get() };
    }

    /// Access solver representing \f$A^{-1}\f$
    ::Spacy::LinearSolver LinearOperator::solver() const
    {
        return IndefiniteLinearSolver();
    }
} // namespace Mock
