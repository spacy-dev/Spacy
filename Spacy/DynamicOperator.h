// This file was automatically generated using clang-type-erase.
// Please do not modify.

#pragma once

#include <Spacy/LinearOperator.h>
#include <Spacy/Spaces/ScalarSpace/Real.h>
#include <Spacy/Util/SmartPointerStorage.h>
#include <Spacy/Vector.h>
#include <Spacy/VectorSpace.h>
#include <memory>
#include <type_traits>
#include <functional>

namespace Spacy
{
    /// A time-dependent operator that does not know about domain and range spaces.
    using DynamicSimpleOperator = std::function< Vector( double, const Vector& ) >;
    /// Type-erased time-dependent operator \f$A:\ [0,T] \times X \to Y \f$.
    class DynamicOperator
    {
        struct Interface
        {
            virtual ~Interface() = default;
            virtual std::shared_ptr< Interface > clone() const = 0;
            virtual Vector call_double_const_Vector_ref( double t, const Vector& x ) const = 0;
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

            Vector call_double_const_Vector_ref( double t, const Vector& x ) const override
            {
                return impl.operator()( std::move( t ), x );
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
        DynamicOperator() noexcept = default;

        template <
            class T,
            typename std::enable_if<
                !std::is_same< typename std::decay< T >::type, DynamicOperator >::value &&
                !std::is_base_of< Interface, typename std::decay< T >::type >::value >::type* =
                nullptr >
        DynamicOperator( T&& value ) : impl_( std::forward< T >( value ) )
        {
        }

        /// Apply operator.
        Vector operator()( double t, const Vector& x ) const
        {
            assert( impl_ );
            return impl_->call_double_const_Vector_ref( std::move( t ), x );
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
                !std::is_same< typename std::decay< T >::type, DynamicOperator >::value &&
                !std::is_base_of< Interface, typename std::decay< T >::type >::value >::type* =
                nullptr >
        DynamicOperator& operator=( T&& value )
        {
            return *this = DynamicOperator( std::forward< T >( value ) );
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
    /// Type-erased time-dependent linear operator \f$A:\ [0,T] \times X \to Y \f$.
    class DynamicLinearOperator
    {
        struct Interface
        {
            virtual ~Interface() = default;
            virtual std::shared_ptr< Interface > clone() const = 0;
            virtual Vector call_double_const_Vector_ref( double t, const Vector& x ) const = 0;
            virtual void add_const_DynamicLinearOperator_ref( const DynamicLinearOperator& y ) = 0;
            virtual void
            subtract_const_DynamicLinearOperator_ref( const DynamicLinearOperator& y ) = 0;
            virtual void multiply_double( double a ) = 0;
            virtual DynamicLinearOperator negate() const = 0;
            virtual bool
            compare_const_DynamicLinearOperator_ref( const DynamicLinearOperator& y ) const = 0;
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

            Vector call_double_const_Vector_ref( double t, const Vector& x ) const override
            {
                return impl.operator()( std::move( t ), x );
            }

            void add_const_DynamicLinearOperator_ref( const DynamicLinearOperator& y ) override
            {
                impl.operator+=( *y.template target< typename std::decay< Impl >::type >() );
            }

            void subtract_const_DynamicLinearOperator_ref( const DynamicLinearOperator& y ) override
            {
                impl.operator-=( *y.template target< typename std::decay< Impl >::type >() );
            }

            void multiply_double( double a ) override
            {
                impl.operator*=( std::move( a ) );
            }

            DynamicLinearOperator negate() const override
            {
                return impl.operator-();
            }

            bool
            compare_const_DynamicLinearOperator_ref( const DynamicLinearOperator& y ) const override
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
        DynamicLinearOperator() noexcept = default;

        template <
            class T,
            typename std::enable_if<
                !std::is_same< typename std::decay< T >::type, DynamicLinearOperator >::value &&
                !std::is_base_of< Interface, typename std::decay< T >::type >::value >::type* =
                nullptr >
        DynamicLinearOperator( T&& value ) : impl_( std::forward< T >( value ) )
        {
        }

        /// Apply operator.
        Vector operator()( double t, const Vector& x ) const
        {
            assert( impl_ );
            return impl_->call_double_const_Vector_ref( std::move( t ), x );
        }

        DynamicLinearOperator& operator+=( const DynamicLinearOperator& y )
        {
            assert( impl_ );
            impl_->add_const_DynamicLinearOperator_ref( y );
            return *this;
        }

        DynamicLinearOperator& operator-=( const DynamicLinearOperator& y )
        {
            assert( impl_ );
            impl_->subtract_const_DynamicLinearOperator_ref( y );
            return *this;
        }

        DynamicLinearOperator& operator*=( double a )
        {
            assert( impl_ );
            impl_->multiply_double( std::move( a ) );
            return *this;
        }

        DynamicLinearOperator operator-() const
        {
            assert( impl_ );
            return impl_->negate();
        }

        bool operator==( const DynamicLinearOperator& y ) const
        {
            assert( impl_ );
            return impl_->compare_const_DynamicLinearOperator_ref( y );
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
                !std::is_same< typename std::decay< T >::type, DynamicLinearOperator >::value &&
                !std::is_base_of< Interface, typename std::decay< T >::type >::value >::type* =
                nullptr >
        DynamicLinearOperator& operator=( T&& value )
        {
            return *this = DynamicLinearOperator( std::forward< T >( value ) );
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
    /// Type-erased time-dependent differentiable operator \f$A:\ [0,T] \times X \to Y \f$.
    class DynamicC1Operator
    {
        struct Interface
        {
            virtual ~Interface() = default;
            virtual std::shared_ptr< Interface > clone() const = 0;
            virtual Vector call_double_const_Vector_ref( double t, const Vector& x ) const = 0;
            virtual Vector d1( double t, const Vector& x, const Vector& dx ) const = 0;
            virtual LinearOperator linearization( double t, const Vector& x ) const = 0;
            virtual LinearOperator M() const = 0;
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

            Vector call_double_const_Vector_ref( double t, const Vector& x ) const override
            {
                return impl.operator()( std::move( t ), x );
            }

            Vector d1( double t, const Vector& x, const Vector& dx ) const override
            {
                return impl.d1( std::move( t ), x, dx );
            }

            LinearOperator linearization( double t, const Vector& x ) const override
            {
                return impl.linearization( std::move( t ), x );
            }

            LinearOperator M() const override
            {
                return impl.M();
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
        DynamicC1Operator() noexcept = default;

        template <
            class T,
            typename std::enable_if<
                !std::is_same< typename std::decay< T >::type, DynamicC1Operator >::value &&
                !std::is_base_of< Interface, typename std::decay< T >::type >::value >::type* =
                nullptr >
        DynamicC1Operator( T&& value ) : impl_( std::forward< T >( value ) )
        {
        }

        /// Apply operator.
        Vector operator()( double t, const Vector& x ) const
        {
            assert( impl_ );
            return impl_->call_double_const_Vector_ref( std::move( t ), x );
        }

        /// Compute directional derivative \f$A'(x)\delta x\f$.
        Vector d1( double t, const Vector& x, const Vector& dx ) const
        {
            assert( impl_ );
            return impl_->d1( std::move( t ), x, dx );
        }

        /// Get linearization \f$A
        LinearOperator linearization( double t, const Vector& x ) const
        {
            assert( impl_ );
            return impl_->linearization( std::move( t ), x );
        }

        LinearOperator M() const
        {
            assert( impl_ );
            return impl_->M();
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
                !std::is_same< typename std::decay< T >::type, DynamicC1Operator >::value &&
                !std::is_base_of< Interface, typename std::decay< T >::type >::value >::type* =
                nullptr >
        DynamicC1Operator& operator=( T&& value )
        {
            return *this = DynamicC1Operator( std::forward< T >( value ) );
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
}
