#pragma once

#include <Spacy/Adaptivity/SpatialAdaptivity.h>
#include <Spacy/Operator.h>
#include <Spacy/Vector.h>

namespace Spacy
{
    class ErrorEstimator
    {
    public:
        Vector estimateError( Operator A, const Vector& x );
        ErrorIndicator getErrorIndicator( const Vector& errorEstimate );
    };
}
