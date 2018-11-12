#pragma once

#include <Eigen/Dense>

#include "Spacy/Adapter/Generic/ScalarProduct.h"

namespace Spacy
{
  namespace Rn
  {
    /**
     * @ingroup EigenGroup
     * @brief Euclidean scalar product for %Rn, based on the %Eigen library.
     */
    using EuclideanScalarProduct = Generic::EuclideanScalarProduct;
  }
}
