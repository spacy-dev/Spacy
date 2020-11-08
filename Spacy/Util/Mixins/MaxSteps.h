#pragma once

#include <Spacy/Util/Mixins/MixinConnection.h>

namespace Spacy::Mixin
{
    class MaxSteps : public MixinConnection< MaxSteps >
    {
        friend class MixinConnection< MaxSteps >;

    public:
        constexpr MaxSteps() = default;

        explicit MaxSteps( unsigned value ) noexcept;

        void setMaxSteps( unsigned value );
        [[nodiscard]] unsigned getMaxSteps() const noexcept;
        void update( MaxSteps* changedSubject );

    private:
        unsigned value_{ 1000 };
    };
} // namespace Spacy::Mixin