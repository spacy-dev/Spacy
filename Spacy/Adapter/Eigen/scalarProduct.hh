#pragma once

#include <eigen3/Eigen/Dense>

#include "Spacy/Adapter/Generic/scalarProduct.hh"

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
