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
    /// Type-erased differentiable functional \f$f:\ X \to \mathbb{R} \f$.
    class C1Functional
    {
        struct Interface
        {
            virtual ~Interface() = default;
            [[nodiscard]] virtual std::shared_ptr< Interface > clone() const = 0;
            [[nodiscard]] virtual Real call_const_Vector_ref( const Vector& x ) const = 0;
            [[nodiscard]] virtual Vector d1( const Vector& x ) const = 0;
            [[nodiscard]] virtual const VectorSpace& domain() const = 0;
        };

        template < class Impl >
        struct Wrapper : Interface
        {
            template < class T >
            Wrapper( T&& t ) : impl( std::forward< T >( t ) )
            {
            }

            [[nodiscard]] std::shared_ptr< Interface > clone() const override
            {
                return std::make_shared< Wrapper< Impl > >( impl );
            }

            Real call_const_Vector_ref( const Vector& x ) const override
            {
                return impl.operator()( x );
            }

            Vector d1( const Vector& x ) const override
            {
                return impl.d1( x );
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
        C1Functional() noexcept = default;

        template < class T,
                   typename std::enable_if< !std::is_same< typename std::decay< T >::type, C1Functional >::value &&
                                            !std::is_base_of< Interface, typename std::decay< T >::type >::value >::type* = nullptr >
        C1Functional( T&& value ) : impl_( std::forward< T >( value ) )
        {
        }

        /// Apply functional.
        [[nodiscard]] Real operator()( const Vector& x ) const
        {
            assert( impl_ );
            return impl_->call_const_Vector_ref( x );
        }

        /// Compute derivative as function space element in \f$X^*\f$, where \f$x\in X\f$.
        [[nodiscard]] Vector d1( const Vector& x ) const
        {
            assert( impl_ );
            return impl_->d1( x );
        }

        /// Access domain space \f$X\f$.
        [[nodiscard]] const VectorSpace& domain() const
        {
            assert( impl_ );
            return impl_->domain();
        }

        template < class T,
                   typename std::enable_if< !std::is_same< typename std::decay< T >::type, C1Functional >::value &&
                                            !std::is_base_of< Interface, typename std::decay< T >::type >::value >::type* = nullptr >
        C1Functional& operator=( T&& value )
        {
            return *this = C1Functional( std::forward< T >( value ) );
        }

        [[nodiscard]] explicit operator bool() const noexcept
        {
            return bool( impl_ );
        }

        template < class T >
        [[nodiscard]] T* target() noexcept
        {
            return impl_.template target< T >();
        }

        template < class T >
        [[nodiscard]] const T* target() const noexcept
        {
            return impl_.template target< T >();
        }

    private:
        clang::type_erasure::polymorphic::SBOCOWStorage< Interface, Wrapper, 16 > impl_;
    };
} // namespace Spacy
