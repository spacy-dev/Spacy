#pragma once

#include <Spacy/Util/Mixins/MixinConnection.h>

namespace Spacy
{
    namespace Mixin
    {
        class StateIndex : public MixinConnection< StateIndex >
        {
            friend class MixinConnection< StateIndex >;

        public:
            explicit StateIndex( unsigned value = 0 ) noexcept;
            void setStateIndex( unsigned value );
            [[nodiscard]] unsigned getStateIndex() const noexcept;
            void update( StateIndex* changedSubject );

        private:
            unsigned value_;
        };

        class ControlIndex : public MixinConnection< ControlIndex >
        {
            friend class MixinConnection< ControlIndex >;

        public:
            explicit ControlIndex( unsigned value = 1 ) noexcept;
            void setControlIndex( unsigned value );
            [[nodiscard]] unsigned getControlIndex() const noexcept;
            void update( ControlIndex* changedSubject );

        private:
            unsigned value_;
        };

        class AdjointIndex : public MixinConnection< AdjointIndex >
        {
            friend class MixinConnection< AdjointIndex >;

        public:
            explicit AdjointIndex( unsigned value = 2 ) noexcept;
            void setAdjointIndex( unsigned value );
            [[nodiscard]] unsigned getAdjointIndex() const noexcept;
            void update( AdjointIndex* changedSubject );

        private:
            unsigned value_;
        };
    } // namespace Mixin
} // namespace Spacy
