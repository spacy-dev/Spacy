#include "norm.hh"

#include <Spacy/Operator.h>

Spacy::Real Mock::Norm::operator()( const ::Spacy::Vector& ) const
{
    return { testValue };
}

Spacy::Real Mock::Norm10::operator()( const ::Spacy::Vector& ) const
{
    return { testValue };
}
