#pragma once

#include <Spacy/C1Operator.h>
#include <Spacy/C2Functional.h>

namespace Spacy
{
    /**
     * @brief Compute the derivative of a functional \f$ f: X\to \mathbb{R} \f$ as (nonlinear)
     * operator \f$ f':X\to X^* \f$.
     * @param f twice differentiable functional
     * @return \f$f'\f$
     */
    C1Operator derivative( C2Functional f );
}
