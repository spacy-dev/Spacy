#include "vectorCreator.hh"

#include <Spacy/Vector.h>
#include <Spacy/VectorSpace.h>

#include "vector.hh"

Mock::Vector Mock::VectorCreator::operator()( const Spacy::VectorSpace* space ) const
{
    return Mock::Vector{*space};
}
