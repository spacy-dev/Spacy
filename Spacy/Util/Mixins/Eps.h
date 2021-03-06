#pragma once

#include "MixinConnection.h"

#include <Spacy/Spaces/ScalarSpace/Real.h>

namespace Spacy::Mixin
{
    /// %Mixin class for maximal attainable accuracy \f$\varepsilon\f$.
    class Eps : public MixinConnection< Eps >
    {
    public:
        constexpr Eps() = default;
        /**
         * @param eps maximal attainable accuracy \f$\varepsilon\f$
         */
        explicit Eps( Real eps ) noexcept;

        void setEps( Real eps );

        /// Access \f$\varepsilon\f$.
        [[nodiscard]] Real eps() const noexcept;

        /// Access \f$\sqrt\varepsilon\f$.
        [[nodiscard]] Real sqrtEps() const noexcept;

        /// Access \f$\varepsilon^{1/3}\f$.
        [[nodiscard]] Real cbrtEps() const noexcept;

        /// update function for observer pattern.
        void update( Eps* changedSubject );

    private:
        Real eps_{ 1e-15 };
        Real sqrtEps_{ 0 };
        Real cbrtEps_{ 0 };
    };
} // namespace Spacy::Mixin
