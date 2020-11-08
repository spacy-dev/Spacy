#pragma once

#include <Spacy/Util/Mixins/MixinConnection.h>

namespace Spacy::Mixin
{
    class NumberOfThreads : public MixinConnection< NumberOfThreads >
    {
        friend class MixinConnection< NumberOfThreads >;

    public:
        explicit NumberOfThreads( unsigned value = 1 ) noexcept;
        void setNumberOfThreads( unsigned value );
        [[nodiscard]] unsigned getNumberOfThreads() const noexcept;
        void update( NumberOfThreads* changedSubject );

    private:
        unsigned value_;
    };
} // namespace Spacy::Mixin