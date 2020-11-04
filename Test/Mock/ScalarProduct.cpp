#include "ScalarProduct.h"

Spacy::Real Mock::ScalarProduct::operator()( const ::Spacy::Vector&, const ::Spacy::Vector& ) const
{
    return { testValue };
}
