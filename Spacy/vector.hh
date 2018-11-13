// This file was automatically generated using clang-type-erase.
// Please do not modify.

#pragma once

#include <Spacy/Spaces/ScalarSpace/Real.h>
#include <Spacy/Util/Exceptions.h>
#include <Spacy/Util/Mixins/Get.h>
#include <Spacy/Util/SmartPointerStorage.h>
#include <Spacy/vectorSpace.hh>
#include <memory>
#include <type_traits>
#include <functional>

namespace Spacy
{
    /// Type-erased vector.
    class Vector
    {
        struct Interface
        {
            virtual ~Interface() = default;
            virtual std::shared_ptr< Interface > clone() const = 0;
            virtual Real call_const_Vector_ref( const Vector& x ) const = 0;
            virtual void add_const_Vector_ref( const Vector& y ) = 0;
            virtual void subtract_const_Vector_ref( const Vector& y ) = 0;
            virtual void multiply_double( double a ) = 0;
            virtual Vector negate() const = 0;
            virtual bool compare_const_Vector_ref( const Vector& y ) const = 0;
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

            Real call_const_Vector_ref( const Vector& x ) const override
            {
                return impl.operator()( *x.template target< typename std::decay< Impl >::type >() );
            }

            void add_const_Vector_ref( const Vector& y ) override
            {
                impl.operator+=( *y.template target< typename std::decay< Impl >::type >() );
            }

            void subtract_const_Vector_ref( const Vector& y ) override
            {
                impl.operator-=( *y.template target< typename std::decay< Impl >::type >() );
            }

            void multiply_double( double a ) override
            {
                impl.operator*=( std::move( a ) );
            }

            Vector negate() const override
            {
                return impl.operator-();
            }

            bool compare_const_Vector_ref( const Vector& y ) const override
            {
                return impl.operator==( *y.template target< typename std::decay< Impl >::type >() );
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
        Vector() noexcept = default;

        template <
            class T,
            typename std::enable_if<
                !std::is_same< typename std::decay< T >::type, Vector >::value &&
                !std::is_base_of< Interface, typename std::decay< T >::type >::value >::type* =
                nullptr >
        Vector( T&& value ) : impl_( std::forward< T >( value ) )
        {
        }

        /// Apply as dual space element.
        Real operator()( const Vector& x ) const
        {
            return impl_->call_const_Vector_ref( x );
        }

        Vector& operator+=( const Vector& y )
        {
            impl_->add_const_Vector_ref( y );
            return *this;
        }

        Vector& operator-=( const Vector& y )
        {
            impl_->subtract_const_Vector_ref( y );
            return *this;
        }

        Vector& operator*=( double a )
        {
            impl_->multiply_double( std::move( a ) );
            return *this;
        }

        Vector operator-() const
        {
            return impl_->negate();
        }

        bool operator==( const Vector& y ) const
        {
            return impl_->compare_const_Vector_ref( y );
        }

        /// Access underlying space.
        const VectorSpace& space() const
        {
            return impl_->space();
        }

        template <
            class T,
            typename std::enable_if<
                !std::is_same< typename std::decay< T >::type, Vector >::value &&
                !std::is_base_of< Interface, typename std::decay< T >::type >::value >::type* =
                nullptr >
        Vector& operator=( T&& value )
        {
            return *this = Vector( std::forward< T >( value ) );
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
    /// Multiplication with arithmetic types (double,float,int,...).
    template < class Arithmetic,
               class = std::enable_if_t< std::is_arithmetic< Arithmetic >::value > >
    Vector operator*( Arithmetic a, Vector x )
    {
        return x *= a;
    }

    /// Multiplication with arithmetic types (double,float,int,...).
    template < class Arithmetic,
               class = std::enable_if_t< std::is_arithmetic< Arithmetic >::value > >
    Vector operator*( Vector x, Arithmetic a )
    {
        return x *= a;
    }

    /// Sum of vectors \f$z=x+y\f$.
    inline Vector operator+( Vector x, const Vector& y )
    {
        return x += y;
    }

    /// Subtract vectors \f$z=x-y\f$.
    inline Vector operator-( Vector x, const Vector& y )
    {
        return x -= y;
    }

    /// Compute scalar product \f$z=x*y=(x,y)\f$.
    inline Real operator*( const Vector& x, const Vector& y )
    {
        return x.space().scalarProduct()( x, y );
    }

    /// Compute norm, where the norm associated with the underlying function space is used \f$ z =
    /// \|x\| \f$.
    inline Real norm( const Vector& x )
    {
        return x.space().norm()( x );
    }

    inline void checkDualPairing( const Vector& x, const Vector& y )
    {
        if ( !y.space().isPrimalWRT( x.space() ) )
            throw Exception::IncompatibleSpace( x.space().index(), y.space().index() );
    }

    template < class T, typename std::enable_if< std::is_arithmetic< T >::value >::type* = nullptr >
    Vector operator*( const Mixin::Get< T >& x, Vector y )
    {
        return y *= get( x );
    }

    template < class T, typename std::enable_if< std::is_arithmetic< T >::value >::type* = nullptr >
    Vector operator*( const Vector& x, const Mixin::Get< T >& y )
    {
        return y * x;
    }

    inline Vector operator*( const Mixin::Get< Real >& x, Vector y )
    {
        return y *= get( get( x ) );
    }

    inline Vector operator*( const Vector& x, const Mixin::Get< Real >& y )
    {
        return y * x;
    }
}
