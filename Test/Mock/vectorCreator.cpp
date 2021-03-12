#include "vectorCreator.hh"

#include "vector.hh"

#include <Spacy/Operator.h>
#include <Spacy/Vector.h>
#include <Spacy/VectorSpace.h>

Mock::Vector Mock::VectorCreator::operator()( const Spacy::VectorSpace* space ) const
{
    return Mock::Vector{ *space };
}
