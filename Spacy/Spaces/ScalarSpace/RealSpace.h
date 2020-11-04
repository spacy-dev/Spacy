#pragma once

#include <Spacy/VectorSpace.h>

namespace Spacy
{
    /**
     * @brief Construct space of real numbers.
     */
    VectorSpace makeRealSpace( bool defaultIndex = false );

    namespace Space
    {
        /**
         * @brief Global space of real numbers with space index 0.
         *
         * This space is used by class Real if the underlying space is left unspecified.
         */
        static VectorSpace R = makeRealSpace( true );
    } // namespace Space
} // namespace Spacy
