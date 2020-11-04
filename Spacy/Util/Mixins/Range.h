#pragma once

namespace Spacy
{
    class VectorSpace;

    namespace Mixin
    {
        class Range
        {
        public:
            explicit Range( const VectorSpace& range );

            Range( Range&& ) = default;
            Range( const Range& ) = default;

            Range& operator=( Range&& other );
            Range& operator=( const Range& other );

            /// Access range space \f$X\f$.
            [[nodiscard]] const VectorSpace& range() const;

        private:
            const VectorSpace& range_;
        };
    } // namespace Mixin
} // namespace Spacy
