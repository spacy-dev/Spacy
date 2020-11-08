// This file was automatically generated using clang-type-erase.
// Please do not modify.

#pragma once

#include <functional>

#include <Spacy/Util/SmartPointerStorage.h>
#include <Spacy/Vector.h>
#include <Spacy/VectorSpace.h>

#include <memory>
#include <type_traits>

namespace Spacy
{
    /// An operator that does not know about domain and range spaces.
    /// Good enough for some algorithms (i.e. for Krylov-methods).
    using CallableOperator = std::function< Vector( const Vector& ) >;
    /// Type-erased operator \f$A:\ X \to Y \f$.
    class Operator
    {
        struct Interface
        {
            virtual ~Interface() = default;
            [[nodiscard]] virtual std::shared_ptr< Interface > clone() const = 0;
            [[nodiscard]] virtual Vector call_const_Vector_ref( const Vector& x ) const = 0;
            [[nodiscard]] virtual const VectorSpace& domain() const = 0;
            [[nodiscard]] virtual const VectorSpace& range() const = 0;
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

            Vector call_const_Vector_ref( const Vector& x ) const override
            {
                return impl.operator()( x );
            }

            const VectorSpace& domain() const override
            {
                return impl.domain();
            }

            const VectorSpace& range() const override
            {
                return impl.range();
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
        Operator() noexcept = default;

        template < class T,
                   typename std::enable_if< !std::is_same< typename std::decay< T >::type, Operator >::value &&
                                            !std::is_base_of< Interface, typename std::decay< T >::type >::value >::type* = nullptr >
        Operator( T&& value ) : impl_( std::forward< T >( value ) )
        {
        }

        /// Apply operator.
        [[nodiscard]] Vector operator()( const Vector& x ) const
        {
            assert( impl_ );
            return impl_->call_const_Vector_ref( x );
        }

        /// Access domain space \f$X\f$.
        [[nodiscard]] const VectorSpace& domain() const
        {
            assert( impl_ );
            return impl_->domain();
        }

        /// Access range space \f$Y\f$.
        [[nodiscard]] const VectorSpace& range() const
        {
            assert( impl_ );
            return impl_->range();
        }

        template < class T,
                   typename std::enable_if< !std::is_same< typename std::decay< T >::type, Operator >::value &&
                                            !std::is_base_of< Interface, typename std::decay< T >::type >::value >::type* = nullptr >
        Operator& operator=( T&& value )
        {
            return *this = Operator( std::forward< T >( value ) );
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
