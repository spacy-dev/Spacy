#pragma once

#include "MixinConnection.h"

#include <Spacy/Algorithm/DampingFactor.h>

namespace Spacy::Mixin
{
    /// %Mixin class for maximal number of steps/iterations.
    class RegularityTest : public MixinConnection< RegularityTest >
    {
    public:
        RegularityTest() = default;

        /**
         * @brief Constructor.
         * @param lowerBound lower bound for regularity test
         */
        explicit RegularityTest( DampingFactor lowerBound ) noexcept;

        /// Set lower bound of regularity test for termination criteria.
        void setLowerBound( DampingFactor lowerBound );

        [[nodiscard]] DampingFactor getLowerBound() const noexcept;

        /**
         * @brief Apply regularity test.
         * @param nu damping factor
         * @return \f$nu > lowerBound_\f$
         */
        [[nodiscard]] bool regularityTestPassed( DampingFactor nu ) const noexcept;

        /**
         * @brief Apply regularity test.
         * @param nu damping factor
         * @return \f$nu <= lowerBound_\f$
         */
        [[nodiscard]] bool regularityTestFailed( DampingFactor nu ) const noexcept;

        /// update function for observer pattern.
        void update( RegularityTest* changedSubject );

    private:
        DampingFactor lowerBound_{ 1e-12 };
    };
} // namespace Spacy::Mixin
