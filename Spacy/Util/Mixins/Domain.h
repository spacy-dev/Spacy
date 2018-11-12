#pragma once

namespace Spacy
{
    class VectorSpace;

    namespace Mixin
    {
        class Domain
        {
        public:
            explicit Domain(const VectorSpace& domain);

            Domain(Domain&&) = default;
            Domain(const Domain&) = default;

            Domain& operator=(Domain&& other);
            Domain& operator=(const Domain& other);

            /// Access domain space \f$X\f$.
            const VectorSpace& domain() const;

        private:
            const VectorSpace& domain_;
        };
    }
}
