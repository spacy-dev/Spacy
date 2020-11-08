// This file was automatically generated using clang-type-erase.
// Please do not modify.

#pragma once

#include <functional>

#include <Spacy/Util/SmartPointerStorage.h>
#include <Spacy/Vector.h>

#include <memory>
#include <type_traits>

namespace Spacy
{
    /// Type-erased linear solver.
    using LinearSolver = std::function< Vector( const Vector& ) >;
    /// Type-erased indefinite linear solver. Additionally monitors if the underlying operator is
    /// positive definite.
    class IndefiniteLinearSolver
    {
        struct Interface
        {
            virtual ~Interface() = default;
            [[nodiscard]] virtual std::shared_ptr< Interface > clone() const = 0;
            [[nodiscard]] virtual Vector call_const_Vector_ref( const Vector& x ) const = 0;
            [[nodiscard]] virtual bool isPositiveDefinite() const = 0;
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

            bool isPositiveDefinite() const override
            {
                return impl.isPositiveDefinite();
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
        IndefiniteLinearSolver() noexcept = default;

        template < class T,
                   typename std::enable_if< !std::is_same< typename std::decay< T >::type, IndefiniteLinearSolver >::value &&
                                            !std::is_base_of< Interface, typename std::decay< T >::type >::value >::type* = nullptr >
        IndefiniteLinearSolver( T&& value ) : impl_( std::forward< T >( value ) )
        {
        }

        [[nodiscard]] Vector operator()( const Vector& x ) const
        {
            assert( impl_ );
            return impl_->call_const_Vector_ref( x );
        }

        [[nodiscard]] bool isPositiveDefinite() const
        {
            assert( impl_ );
            return impl_->isPositiveDefinite();
        }

        template < class T,
                   typename std::enable_if< !std::is_same< typename std::decay< T >::type, IndefiniteLinearSolver >::value &&
                                            !std::is_base_of< Interface, typename std::decay< T >::type >::value >::type* = nullptr >
        IndefiniteLinearSolver& operator=( T&& value )
        {
            return *this = IndefiniteLinearSolver( std::forward< T >( value ) );
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
        clang::type_erasure::polymorphic::COWStorage< Interface, Wrapper > impl_;
    };
} // namespace Spacy
