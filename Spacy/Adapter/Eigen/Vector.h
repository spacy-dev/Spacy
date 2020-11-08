#pragma once

#include <Eigen/Dense>

#include <Spacy/Adapter/Generic/Vector.h>

namespace Spacy::Rn
{
    /**
     * @ingroup EigenGroup, VectorSpaceGroup
     * @brief %Vector for %Rn, based on the %Eigen library.
     */
    using Vector = Generic::Vector< ::Eigen::VectorXd >;

    unsigned getSize( const ::Spacy::Vector& y );
} // namespace Spacy::Rn
