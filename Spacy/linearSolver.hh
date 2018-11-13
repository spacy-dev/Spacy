// This file was automatically generated using clang-type-erase.
// Please do not modify.

#pragma once

#include <Spacy/Util/SmartPointerStorage.h>
#include <Spacy/vector.hh>
#include <memory>
#include <type_traits>
#include <functional>

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
            virtual std::unique_ptr< Interface > clone() const = 0;
            virtual Vector call_const_Vector_ref( const Vector& x ) const = 0;
            virtual bool isPositiveDefinite() const = 0;
        };

        template < class Impl >
        struct Wrapper : Interface
        {
            template < class T >
            Wrapper( T&& t ) : impl( std::forward< T >( t ) )
            {
            }

            std::unique_ptr< Interface > clone() const override
            {
                return std::make_unique< Wrapper< Impl > >( impl );
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

        template <
            class T,
            typename std::enable_if<
                !std::is_same< typename std::decay< T >::type, IndefiniteLinearSolver >::value &&
                !std::is_base_of< Interface, typename std::decay< T >::type >::value >::type* =
                nullptr >
        IndefiniteLinearSolver( T&& value ) : impl_( std::forward< T >( value ) )
        {
        }

        Vector operator()( const Vector& x ) const
        {
            return impl_->call_const_Vector_ref( x );
        }

        bool isPositiveDefinite() const
        {
            return impl_->isPositiveDefinite();
        }

        template <
            class T,
            typename std::enable_if<
                !std::is_same< typename std::decay< T >::type, IndefiniteLinearSolver >::value &&
                !std::is_base_of< Interface, typename std::decay< T >::type >::value >::type* =
                nullptr >
        IndefiniteLinearSolver& operator=( T&& value )
        {
            return *this = IndefiniteLinearSolver( std::forward< T >( value ) );
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
        clang::type_erasure::polymorphic::Storage< Interface, Wrapper > impl_;
    };
}
