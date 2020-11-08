#pragma once

#include <Eigen/Dense>

#include <Spacy/Vector.h>

/** @addtogroup EigenGroup @{ */
namespace Spacy::Rn
{
    ///  Copy ::Spacy::Vector to flat coefficient vector of %Eigen .
    void copy( const ::Spacy::Vector& x, Eigen::VectorXd& y );

    ///  Copy flat coefficient vector of %Eigen to ::Spacy::Vector.
    void copy( const Eigen::VectorXd& x, ::Spacy::Vector& y );
} // namespace Spacy::Rn
/** @} */
