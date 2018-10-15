#pragma once

#include <Spacy/Spaces/RealSpace.h>
#include <Spacy/vector.hh>

#include "norm.hh"

namespace Mock
{
    struct ScalarProduct
    {
        static constexpr int testValue = Norm::testValue * Norm::testValue;
        Spacy::Real operator()( const ::Spacy::Vector&, const ::Spacy::Vector& ) const;
    };
}
