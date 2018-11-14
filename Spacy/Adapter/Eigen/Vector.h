#pragma once

#include <Spacy/Adapter/Generic/Vector.h>
#include <Eigen/Dense>

namespace Spacy
{
    namespace Rn
    {
        /**
         * @ingroup EigenGroup, VectorSpaceGroup
         * @brief %Vector for %Rn, based on the %Eigen library.
         */
        using Vector = Generic::Vector<::Eigen::VectorXd >;

        unsigned getSize( const ::Spacy::Vector& y );
    }
}
