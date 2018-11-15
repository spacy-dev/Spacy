#pragma once

#include <Spacy/Adapter/Generic/Vector.h>
#include <deal.II/lac/vector.h>

namespace Spacy
{
    namespace dealII
    {
        /**
         * @ingroup dealIIGroup, VectorSpaceGroup
         * @brief %Vector adapter for deal.II
         */
        using Vector = Generic::Vector< dealii::Vector< double > >;
    }
}
