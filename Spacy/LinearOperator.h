// This file was automatically generated using clang-type-erase.
// Please do not modify.

#pragma once

#include <Spacy/LinearSolver.h>
#include <Spacy/Spaces/ScalarSpace/Real.h>
#include <Spacy/Util/Exceptions.h>
#include <Spacy/Util/SmartPointerStorage.h>
#include <Spacy/Vector.h>
#include <Spacy/VectorSpace.h>
#include <memory>
#include <type_traits>
#include <functional>

namespace Spacy
{
    /// Type-erased linear operator \f$A:\ X \to Y \f$.
    class LinearOperator
    {
        struct Interface
        {
            virtual ~Interface() = default;
            virtual std::shared_ptr< Interface > clone() const = 0;
            virtual Vector call_const_Vector_ref( const Vector& x ) const = 0;
            virtual Real call_const_LinearOperator_ref( const LinearOperator& x ) const = 0;
            virtual void add_const_LinearOperator_ref( const LinearOperator& y ) = 0;
            virtual void subtract_const_LinearOperator_ref( const LinearOperator& y ) = 0;
            virtual void multiply_double( double a ) = 0;
            virtual LinearOperator negate() const = 0;
            virtual bool compare_const_LinearOperator_ref( const LinearOperator& y ) const = 0;
            virtual std::function< Vector( const Vector& ) > solver() const = 0;
            virtual const VectorSpace& domain() const = 0;
            virtual const VectorSpace& range() const = 0;
            virtual const VectorSpace& space() const = 0;
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

            Real call_const_LinearOperator_ref( const LinearOperator& x ) const override
            {
                return impl.operator()( *x.template target< typename std::decay< Impl >::type >() );
            }

            void add_const_LinearOperator_ref( const LinearOperator& y ) override
            {
                impl.operator+=( *y.template target< typename std::decay< Impl >::type >() );
            }

            void subtract_const_LinearOperator_ref( const LinearOperator& y ) override
            {
                impl.operator-=( *y.template target< typename std::decay< Impl >::type >() );
            }

            void multiply_double( double a ) override
            {
                impl.operator*=( std::move( a ) );
            }

            LinearOperator negate() const override
            {
                return impl.operator-();
            }

            bool compare_const_LinearOperator_ref( const LinearOperator& y ) const override
            {
                return impl.operator==( *y.template target< typename std::decay< Impl >::type >() );
            }

            std::function< Vector( const Vector& ) > solver() const override
            {
                return impl.solver();
            }

            const VectorSpace& domain() const override
            {
                return impl.domain();
            }

            const VectorSpace& range() const override
            {
                return impl.range();
            }

            const VectorSpace& space() const override
            {
                return impl.space();
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
        LinearOperator() noexcept = default;

        template <
            class T,
            typename std::enable_if<
                !std::is_same< typename std::decay< T >::type, LinearOperator >::value &&
                !std::is_base_of< Interface, typename std::decay< T >::type >::value >::type* =
                nullptr >
        LinearOperator( T&& value ) : impl_( std::forward< T >( value ) )
        {
        }

        /// Apply operator.
        Vector operator()( const Vector& x ) const
        {
            assert( impl_ );
            return impl_->call_const_Vector_ref( x );
        }

        Real operator()( const LinearOperator& x ) const
        {
            assert( impl_ );
            return impl_->call_const_LinearOperator_ref( x );
        }

        LinearOperator& operator+=( const LinearOperator& y )
        {
            assert( impl_ );
            impl_->add_const_LinearOperator_ref( y );
            return *this;
        }

        LinearOperator& operator-=( const LinearOperator& y )
        {
            assert( impl_ );
            impl_->subtract_const_LinearOperator_ref( y );
            return *this;
        }

        LinearOperator& operator*=( double a )
        {
            assert( impl_ );
            impl_->multiply_double( std::move( a ) );
            return *this;
        }

        LinearOperator operator-() const
        {
            assert( impl_ );
            return impl_->negate();
        }

        bool operator==( const LinearOperator& y ) const
        {
            assert( impl_ );
            return impl_->compare_const_LinearOperator_ref( y );
        }

        std::function< Vector( const Vector& ) > solver() const
        {
            assert( impl_ );
            return impl_->solver();
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

        /// Access underlying space of linear operators.
        const VectorSpace& space() const
        {
            assert( impl_ );
            return impl_->space();
        }

        template <
            class T,
            typename std::enable_if<
                !std::is_same< typename std::decay< T >::type, LinearOperator >::value &&
                !std::is_base_of< Interface, typename std::decay< T >::type >::value >::type* =
                nullptr >
        LinearOperator& operator=( T&& value )
        {
            return *this = LinearOperator( std::forward< T >( value ) );
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
    /// Access solver via A^-1. Throws for k!=-1.
    inline LinearSolver operator^( const LinearOperator& A, int k )
    {
        if ( k == -1 )
            return A.solver();
        throw Exception::InvalidArgument(
            "operator^ for LinearOperator only defined for exponent: k = -1." );
    }

    /// Access solver via A^-1. Throws for k!=-1.
    inline LinearSolver operator^( LinearOperator&& A, int k )
    {
        if ( k == -1 )
            return std::move( A.solver() );
        throw Exception::InvalidArgument(
            "operator^ for LinearOperator only defined for exponent: k = -1." );
    }

    inline LinearOperator& axpy( LinearOperator& A, double a, LinearOperator B )
    {
        return A += ( B *= a );
    }

} // namespace Spacy
