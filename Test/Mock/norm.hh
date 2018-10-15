#pragma once

#include <Spacy/Spaces/RealSpace.h>
#include <Spacy/vector.hh>

namespace Mock
{
    struct Norm
    {
        static constexpr int testValue = 5;
        Spacy::Real operator()( const ::Spacy::Vector& ) const;
    };

    struct Norm10
    {
        static constexpr int testValue = 10;
        Spacy::Real operator()( const ::Spacy::Vector& ) const;
    };
}
