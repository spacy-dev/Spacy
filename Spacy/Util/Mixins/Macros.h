#include "MixinConnection.h"

#define GENERATE_HEADER_START( NAME )                                                                                                      \
    namespace Spacy                                                                                                                        \
    {                                                                                                                                      \
        namespace Mixin                                                                                                                    \
        {                                                                                                                                  \
            class NAME : public MixinConnection< NAME >                                                                                    \
            {                                                                                                                              \
                friend class MixinConnection< NAME >;                                                                                      \
                                                                                                                                           \
            public:

#define GENERATE_HEADER_BASIC_CONTENT( TYPE, NAME, INITIAL )                                                                               \
    explicit NAME( TYPE value = INITIAL ) noexcept;                                                                                        \
    void set##NAME( TYPE value );                                                                                                          \
    TYPE get##NAME() const noexcept;

#define GENERATE_HEADER_END( TYPE, NAME )                                                                                                  \
private:                                                                                                                                   \
    void update( NAME* changedSubject );                                                                                                   \
    TYPE value_;                                                                                                                           \
    }                                                                                                                                      \
    ;                                                                                                                                      \
    }                                                                                                                                      \
    }

// clang-format off
#define GENERATE_MIXIN_HEADER( TYPE, NAME, INITIAL )                                                                                       \
    namespace Spacy                                                                                                                        \
    {                                                                                                                                      \
        namespace Mixin                                                                                                                    \
        {                                                                                                                                  \
            class NAME : public MixinConnection< NAME >                                                                                    \
            {                                                                                                                              \
                friend class MixinConnection< NAME >;                                                                                      \
                                                                                                                                           \
            public: \
                explicit NAME( TYPE value = INITIAL ) noexcept;                                                                                        \
                void set##NAME( TYPE value );                                                                                                          \
                TYPE get##NAME() const noexcept; \
                void update( NAME* changedSubject );                                                                                                   \
            private:                                                                                                                                   \
                TYPE value_;                                                                                                                           \
            };                                                                                                                                      \
        }                                                                                                                                      \
    }
// clang-format on

#define GENERATE_MIXIN_SOURCE( TYPE, NAME )                                                                                                \
    Spacy::Mixin::NAME::NAME( TYPE value ) noexcept : value_( value )                                                                      \
    {                                                                                                                                      \
    }                                                                                                                                      \
                                                                                                                                           \
    void Spacy::Mixin::NAME::set##NAME( TYPE value )                                                                                       \
    {                                                                                                                                      \
        value_ = value;                                                                                                                    \
        notify();                                                                                                                          \
    }                                                                                                                                      \
                                                                                                                                           \
    TYPE Spacy::Mixin::NAME::get##NAME() const noexcept                                                                                    \
    {                                                                                                                                      \
        return value_;                                                                                                                     \
    }                                                                                                                                      \
                                                                                                                                           \
    void Spacy::Mixin::NAME::update( NAME* changedSubject )                                                                                \
    {                                                                                                                                      \
        set##NAME( changedSubject->get##NAME() );                                                                                          \
    }
