#pragma once

namespace Spacy
{
    namespace Util
    {
        template <class Vector>
        void registerVector(Vector* v)
        {
            if(*v)
                v->space().add(v);
        }

        template <class Vector>
        void deregisterVector(Vector* v)
        {
            if(*v)
                v->space().remove(v);
        }

        template <class Vector>
        class CopyMoveConstructRegistration
        {
        public:
            CopyMoveConstructRegistration() = default;

            CopyMoveConstructRegistration(const CopyMoveConstructRegistration& other) noexcept
            {
                registerVector(other);
            }

            CopyMoveConstructRegistration(CopyMoveConstructRegistration&& other) noexcept
            {
                registerVector(other);
            }

            CopyMoveConstructRegistration& operator=(const CopyMoveConstructRegistration& other)
            {
                deregisterVector();
                registerVector(other);
                return *this;
            }

            CopyMoveConstructRegistration& operator=(CopyMoveConstructRegistration&& other)
            {
                deregisterVector();
                registerVector(other);
                return *this;
            }

        private:
            template <class T>
            void registerVector(T&& other)
            {
                const auto& otherAsVector = static_cast<const Vector&>(other);
                if(otherAsVector)
                    otherAsVector.space().add(static_cast<Vector*>(this));
            }

            void deregisterVector()
            {
                if(static_cast<const Vector&>(*this))
                    static_cast<const Vector*>(this)->space().remove(static_cast<Vector*>(this));
            }
        };
    }
}
