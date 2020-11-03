#pragma once

#include <Spacy/Spaces/ScalarSpace/Real.h>
#include <Spacy/Util/Mixins/MixinConnection.h>

namespace Spacy
{
    namespace Mixin
    {
        class AbsoluteAccuracy : public MixinConnection< AbsoluteAccuracy >
        {
            friend class MixinConnection< AbsoluteAccuracy >;

        public:
            explicit AbsoluteAccuracy( Real value = 1e-15 ) noexcept;
            void setAbsoluteAccuracy( Real value );
            Real getAbsoluteAccuracy() const noexcept;
            void update( AbsoluteAccuracy* changedSubject );

        private:
            Real value_;
        };

        class RelativeAccuracy : public MixinConnection< RelativeAccuracy >
        {
            friend class MixinConnection< RelativeAccuracy >;

        public:
            explicit RelativeAccuracy( Real value = 1e-15 ) noexcept;
            void setRelativeAccuracy( Real value );
            Real getRelativeAccuracy() const noexcept;
            void update( RelativeAccuracy* changedSubject );

        private:
            Real value_;
        };

        class MinimalAccuracy : public MixinConnection< MinimalAccuracy >
        {
            friend class MixinConnection< MinimalAccuracy >;

        public:
            explicit MinimalAccuracy( Real value = 1e-15 ) noexcept;
            void setMinimalAccuracy( Real value );
            Real getMinimalAccuracy() const noexcept;
            void update( MinimalAccuracy* changedSubject );

        private:
            Real value_;
        };

        class DampingAccuracy : public MixinConnection< DampingAccuracy >
        {
            friend class MixinConnection< DampingAccuracy >;

        public:
            explicit DampingAccuracy( Real value = 1e-15 ) noexcept;
            void setDampingAccuracy( Real value );
            Real getDampingAccuracy() const noexcept;
            void update( DampingAccuracy* changedSubject );

        private:
            Real value_;
        };

        class NormalAccuracy : public MixinConnection< NormalAccuracy >
        {
            friend class MixinConnection< NormalAccuracy >;

        public:
            explicit NormalAccuracy( Real value = 1e-15 ) noexcept;
            void setNormalAccuracy( Real value );
            Real getNormalAccuracy() const noexcept;
            void update( NormalAccuracy* changedSubject );

        private:
            Real value_;
        };

        class TangentialAccuracy : public MixinConnection< TangentialAccuracy >
        {
            friend class MixinConnection< TangentialAccuracy >;

        public:
            explicit TangentialAccuracy( Real value = 1e-15 ) noexcept;
            void setTangentialAccuracy( Real value );
            Real getTangentialAccuracy() const noexcept;
            void update( TangentialAccuracy* changedSubject );

        private:
            Real value_;
        };

        class FallBackTangentialAccuracy : public MixinConnection< FallBackTangentialAccuracy >
        {
            friend class MixinConnection< FallBackTangentialAccuracy >;

        public:
            explicit FallBackTangentialAccuracy( Real value = 1e-15 ) noexcept;
            void setFallBackTangentialAccuracy( Real value );
            Real getFallBackTangentialAccuracy() const noexcept;
            void update( FallBackTangentialAccuracy* changedSubject );

        private:
            Real value_;
        };
    } // namespace Mixin
} // namespace Spacy
