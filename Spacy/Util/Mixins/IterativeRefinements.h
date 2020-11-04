#pragma once

#include <Spacy/Util/Mixins/MixinConnection.h>

namespace Spacy
{
    namespace Mixin
    {
        class IterativeRefinements : public MixinConnection< IterativeRefinements >
        {
            friend class MixinConnection< IterativeRefinements >;

        public:
            explicit IterativeRefinements( unsigned value = 1000 ) noexcept;
            void setIterativeRefinements( unsigned value );
            [[nodiscard]] unsigned getIterativeRefinements() const noexcept;
            void update( IterativeRefinements* changedSubject );

        private:
            unsigned value_;
        };
    } // namespace Mixin
} // namespace Spacy
