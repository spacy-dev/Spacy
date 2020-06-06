// This file was automatically generated using clang-type-erase.
// Please do not modify.

#pragma once

#include <Spacy/Spaces/ScalarSpace/Real.h>
#include <Spacy/Util/SmartPointerStorage.h>
#include <Spacy/Vector.h>
#include <Spacy/VectorSpace.h>
#include <memory>
#include <type_traits>

namespace Spacy
{
    /// Type-erased functional \f$f:\ X \to \mathbb{R} \f$.
    class Functional
    {
        struct Interface
        {
            virtual ~Interface() = default;
            virtual std::shared_ptr< Interface > clone() const = 0;
            virtual Real call_const_Vector_ref( const Vector& x ) const = 0;
            virtual const VectorSpace& domain() const = 0;
        };

        template < class Impl >
        struct Wrapper : Interface
        {
            template < class T >
            Wrapper( T&& t ) : impl( std::forward< T >( t ) )
            {
            }

            std::shared_ptr< Interface > clone() const override
            {
                return std::make_shared< Wrapper< Impl > >( impl );
            }

            Real call_const_Vector_ref( const Vector& x ) const override
            {
                return impl.operator()( x );
            }

            const VectorSpace& domain() const override
            {
                return impl.domain();
            }

            Impl impl;
        };

        template < class Impl >
        struct Wrapper< std::reference_wrapper< Impl > > : Wrapper< Impl& >
        {
            template < class T >
            Wrapper( T&& t ) : Wrapper< Impl& >( std::forward< T >( t ) )
            {
            }
        };

    public:
        Functional() noexcept = default;

        template <
            class T,
            typename std::enable_if<
                !std::is_same< typename std::decay< T >::type, Functional >::value &&
                !std::is_base_of< Interface, typename std::decay< T >::type >::value >::type* =
                nullptr >
        Functional( T&& value ) : impl_( std::forward< T >( value ) )
        {
        }

        /// Apply functional.
        Real operator()( const Vector& x ) const
        {
            assert( impl_ );
            return impl_->call_const_Vector_ref( x );
        }

        /// Access domain space \f$X\f$.
        const VectorSpace& domain() const
        {
            assert( impl_ );
            return impl_->domain();
        }

        template <
            class T,
            typename std::enable_if<
                !std::is_same< typename std::decay< T >::type, Functional >::value &&
                !std::is_base_of< Interface, typename std::decay< T >::type >::value >::type* =
                nullptr >
        Functional& operator=( T&& value )
        {
            return *this = Functional( std::forward< T >( value ) );
        }

        explicit operator bool() const noexcept
        {
            return bool( impl_ );
        }

        template < class T >
        T* target() noexcept
        {
            return impl_.template target< T >();
        }

        template < class T >
        const T* target() const noexcept
        {
            return impl_.template target< T >();
        }

    private:
        clang::type_erasure::polymorphic::SBOCOWStorage< Interface, Wrapper, 16 > impl_;
    };
} // namespace Spacy
