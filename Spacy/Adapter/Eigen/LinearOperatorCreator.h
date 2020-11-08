#pragma once

#include "LinearOperator.h"

#include <Spacy/Adapter/Generic/LinearOperatorCreator.h>

namespace Spacy::Rn
{
    /**
     * @ingroup EigenGroup
     * @brief %VectorCreator for linear operators \f$X\rightarrow Y\f$, based on the %Eigen
     * library.
     */
    struct LinearOperatorCreator : Generic::LinearOperatorCreator< LinearOperator >
    {
        LinearOperatorCreator( const VectorSpace& X, const VectorSpace& Y );
    };
} // namespace Spacy::Rn
