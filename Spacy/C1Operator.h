// This file was automatically generated using clang-type-erase.
// Please do not modify.

#pragma once

#include <Spacy/LinearOperator.h>
#include <Spacy/Util/SmartPointerStorage.h>
#include <Spacy/Vector.h>
#include <Spacy/VectorSpace.h>
#include <memory>
#include <type_traits>

namespace Spacy
{
    /// Type-erased differentiable operator \f$A:\ X \to Y \f$.
    class C1Operator
    {
        struct Interface
        {
            virtual ~Interface() = default;
            virtual std::shared_ptr< Interface > clone() const = 0;
            virtual Vector call_const_Vector_ref( const Vector& x ) const = 0;
            virtual Vector d1( const Vector& x, const Vector& dx ) const = 0;
            virtual LinearOperator linearization( const Vector& x ) const = 0;
            virtual const VectorSpace& domain() const = 0;
            virtual const VectorSpace& range() const = 0;
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

            Vector call_const_Vector_ref( const Vector& x ) const override
            {
                return impl.operator()( x );
            }

            Vector d1( const Vector& x, const Vector& dx ) const override
            {
                return impl.d1( x, dx );
            }

            LinearOperator linearization( const Vector& x ) const override
            {
                return impl.linearization( x );
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
        C1Operator() noexcept = default;

        template <
            class T,
            typename std::enable_if<
                !std::is_same< typename std::decay< T >::type, C1Operator >::value &&
                !std::is_base_of< Interface, typename std::decay< T >::type >::value >::type* =
                nullptr >
        C1Operator( T&& value ) : impl_( std::forward< T >( value ) )
        {
        }

        /// Apply operator.
        Vector operator()( const Vector& x ) const
        {
            assert( impl_ );
            return impl_->call_const_Vector_ref( x );
        }

        /// Compute directional derivative \f$A'(x)\delta x\f$.
        Vector d1( const Vector& x, const Vector& dx ) const
        {
            assert( impl_ );
            return impl_->d1( x, dx );
        }

        /// Get linearization \f$A
        LinearOperator linearization( const Vector& x ) const
        {
            assert( impl_ );
            return impl_->linearization( x );
        }

        /// Access domain space \f$X\f$.
        const VectorSpace& domain() const
        {
            assert( impl_ );
            return impl_->domain();
        }

        /// Access range space \f$Y\f$.
        const VectorSpace& range() const
        {
            assert( impl_ );
            return impl_->range();
        }

        template <
            class T,
            typename std::enable_if<
                !std::is_same< typename std::decay< T >::type, C1Operator >::value &&
                !std::is_base_of< Interface, typename std::decay< T >::type >::value >::type* =
                nullptr >
        C1Operator& operator=( T&& value )
        {
            return *this = C1Operator( std::forward< T >( value ) );
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
        clang::type_erasure::polymorphic::SBOStorage< Interface, Wrapper, 16 > impl_;
    };
    /// @brief For an operator \f$ A: X\to Y \f$, compute \f$A'\f$ at \f$x\in X\f$ as linear
    /// operator
    ///\f$ A'(x): X \to Y \f$.
    ///     *
    /// Equivalent to calling A.linearization(x).
    ///     *
    ///@param A differentiable operator
    ///@param x point of linearization
    ///@return \f$A'(x)\f$, i.e. A.linearization(x).
    ///@see @ref C1OperatorAnchor "C1Operator", @ref LinearOperatorAnchor "LinearOperator"
    ///
    inline LinearOperator d1( const C1Operator& A, const Vector& x )
    {
        return A.linearization( x );
    }
}
