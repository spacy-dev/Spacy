#pragma once

#include <Spacy/Spaces/ScalarSpace/Real.h>
#include <Spacy/Util/Mixins/MixinConnection.h>

namespace Spacy::Mixin
{
    class AbsoluteAccuracy : public MixinConnection< AbsoluteAccuracy >
    {
        friend class MixinConnection< AbsoluteAccuracy >;

    public:
        AbsoluteAccuracy() = default;
        explicit AbsoluteAccuracy( Real value ) noexcept;
        void setAbsoluteAccuracy( Real value );
        [[nodiscard]] Real getAbsoluteAccuracy() const noexcept;
        void update( AbsoluteAccuracy* changedSubject );

    private:
        Real value_{ 1e-15 };
    };

    class RelativeAccuracy : public MixinConnection< RelativeAccuracy >
    {
        friend class MixinConnection< RelativeAccuracy >;

    public:
        RelativeAccuracy() = default;
        explicit RelativeAccuracy( Real value ) noexcept;
        void setRelativeAccuracy( Real value );
        [[nodiscard]] Real getRelativeAccuracy() const noexcept;
        void update( RelativeAccuracy* changedSubject );

    private:
        Real value_{ 1e-15 };
    };

    class MinimalAccuracy : public MixinConnection< MinimalAccuracy >
    {
        friend class MixinConnection< MinimalAccuracy >;

    public:
        MinimalAccuracy() = default;
        explicit MinimalAccuracy( Real value ) noexcept;
        void setMinimalAccuracy( Real value );
        [[nodiscard]] Real getMinimalAccuracy() const noexcept;
        void update( MinimalAccuracy* changedSubject );

    private:
        Real value_{ 1e-15 };
    };

    class DampingAccuracy : public MixinConnection< DampingAccuracy >
    {
        friend class MixinConnection< DampingAccuracy >;

    public:
        DampingAccuracy() = default;
        explicit DampingAccuracy( Real value ) noexcept;
        void setDampingAccuracy( Real value );
        [[nodiscard]] Real getDampingAccuracy() const noexcept;
        void update( DampingAccuracy* changedSubject );

    private:
        Real value_{ 1e-15 };
    };

    class NormalAccuracy : public MixinConnection< NormalAccuracy >
    {
        friend class MixinConnection< NormalAccuracy >;

    public:
        NormalAccuracy() = default;
        explicit NormalAccuracy( Real value ) noexcept;
        void setNormalAccuracy( Real value );
        [[nodiscard]] Real getNormalAccuracy() const noexcept;
        void update( NormalAccuracy* changedSubject );

    private:
        Real value_{ 1e-15 };
    };

    class TangentialAccuracy : public MixinConnection< TangentialAccuracy >
    {
        friend class MixinConnection< TangentialAccuracy >;

    public:
        TangentialAccuracy() = default;
        explicit TangentialAccuracy( Real value ) noexcept;
        void setTangentialAccuracy( Real value );
        [[nodiscard]] Real getTangentialAccuracy() const noexcept;
        void update( TangentialAccuracy* changedSubject );

    private:
        Real value_{ 1e-15 };
    };

    class FallBackTangentialAccuracy : public MixinConnection< FallBackTangentialAccuracy >
    {
        friend class MixinConnection< FallBackTangentialAccuracy >;

    public:
        FallBackTangentialAccuracy() = default;
        explicit FallBackTangentialAccuracy( Real value ) noexcept;
        void setFallBackTangentialAccuracy( Real value );
        [[nodiscard]] Real getFallBackTangentialAccuracy() const noexcept;
        void update( FallBackTangentialAccuracy* changedSubject );

    private:
        Real value_{ 1e-15 };
    };
} // namespace Spacy::Mixin
