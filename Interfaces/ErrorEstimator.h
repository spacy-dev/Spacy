#pragma once

#include <Spacy/Adaptivity/SpatialAdaptivity.h>
#include <Spacy/Operator.h>
#include <Spacy/Vector.h>

namespace Spacy
{
    class ErrorEstimator
    {
    public:
        bool estimateError( const Vector& x, const Vector& dx );
        ErrorIndicator getErrorIndicator( const Vector& errorEstimate );
    };
}
