#pragma once

#include <Spacy/Util/Mixins/MixinConnection.h>

namespace Spacy
{
    namespace Mixin
    {
        class MaxSteps : public MixinConnection< MaxSteps >
        {
            friend class MixinConnection< MaxSteps >;

        public:
            explicit MaxSteps( unsigned value = 1000 ) noexcept;
            void setMaxSteps( unsigned value );
            [[nodiscard]] unsigned getMaxSteps() const noexcept;
            void update( MaxSteps* changedSubject );

        private:
            unsigned value_;
        };
    } // namespace Mixin
} // namespace Spacy