// This file was automatically generated using clang-type-erase.
// Please do not modify.
// TODO: This file WAS modified by hand. Interface has to be re-written to accomodate for these changes
// Also added a Vector.cpp file

#pragma once

#include <Spacy/Adaptivity/SpaceManager.h>
#include <Spacy/ForwardIterator.h>
#include <Spacy/Spaces/ScalarSpace/Real.h>
#include <Spacy/Util/Exceptions.h>
#include <Spacy/Util/Mixins/Get.h>
#include <Spacy/Util/SmartPointerStorage.h>
#include <Spacy/VectorSpace.h>
#include <memory>
#include <type_traits>
#include <functional>
#include <iostream>

namespace Spacy
{
    class IdEmbedding;
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
            virtual ForwardIterator begin() = 0;
            virtual ForwardIterator end() = 0;
            virtual ConstForwardIterator begin() const = 0;
            virtual ConstForwardIterator end() const = 0;
            virtual void changeSpace(const VectorSpace& V) = 0;
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

            ///TODO: Was modified by hand, how to do this in the "Interface" definition?
            void add_const_Vector_ref( const Vector& y ) override
            {
                /// Check, if an embedding for the added vector is available, if yes, apply 
                if(this->space().index() != y.space().index())
                {
                  impl.operator+=( *(y.space().getEmbedding(this->space()).apply(y)).template target< typename std::decay< Impl >::type >() );
                } else
                {
                  impl.operator+=( *y.template target< typename std::decay< Impl >::type >() );
                }
            }

            ///TODO: Was modified by hand, how to do this in the "Interface" definition?
            void subtract_const_Vector_ref( const Vector& y ) override
            {
                /// Check, if an embedding for the subtracted vector is available, if yes, apply 
                if(this->space().index() != y.space().index())
                {
                  impl.operator-=( *(y.space().getEmbedding(this->space()).apply(y)).template target< typename std::decay< Impl >::type >() );
                } else
                {
                  impl.operator-=( *y.template target< typename std::decay< Impl >::type >() );
                }
          }

            void multiply_double( double a ) override
            {
                impl.operator*=( std::move( a ) );
            }

            Vector negate() const override
            {
                return impl.operator-();
            }

            ///TODO: Was modified by hand, how to do this in the "Interface" definition?
            bool compare_const_Vector_ref( const Vector& y ) const override
            {
                /// Check, if an embedding for the to be compared vector is available, if yes, apply 
                if(this->space().index() == y.space().index())
                  return impl.operator==( *y.template target< typename std::decay< Impl >::type >() );
                return impl.operator==( *(y.space().getEmbedding(this->space()).apply(y)).template target< typename std::decay< Impl >::type >() );
                return false;
	    }

            const VectorSpace& space() const override
            {
                return impl.space();
            }

            ForwardIterator begin() override
            {
                return impl.begin();
            }

            ForwardIterator end() override
            {
                return impl.end();
            }

            ConstForwardIterator begin() const override
            {
                return impl.begin();
            }

            ConstForwardIterator end() const override
            {
                return impl.end();
            }

            void changeSpace(const VectorSpace& V) override
            {
            impl.changeSpace(V);
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

        template < class T, typename std::enable_if<
                                !std::is_same< typename std::decay< T >::type, Vector >::value &&
                                !std::is_base_of< Interface, typename std::decay< T >::type >::
                                    value >::type* = nullptr >
        Vector( T&& value ) : impl_( std::forward< T >( value ) )
        {
            globalSpaceManager().subscribe( this );
        }

        Vector( const Vector& v ) : impl_( v.impl_ )
        {
            if ( impl_ )
                globalSpaceManager().subscribe( this );
        }

        Vector( Vector&& v ) : impl_( std::move( v ).impl_ )
        {
            if ( impl_ )
                globalSpaceManager().subscribe( this );
        }

        Vector& operator=( const Vector& v );

        Vector& operator=( Vector&& v );

        ~Vector()
        {
            if ( impl_ )
                globalSpaceManager().unsubscribe( this );
        }

        /// Apply as dual space element.
        Real operator()( const Vector& x ) const
        {
            assert( impl_ );
            return impl_->call_const_Vector_ref( x );
        }

        Vector& operator+=( const Vector& y )
        {
            assert( impl_ );
            impl_->add_const_Vector_ref( y );
            return *this;
        }

        Vector& operator-=( const Vector& y )
        {
            assert( impl_ );
            impl_->subtract_const_Vector_ref( y );
            return *this;
        }

        Vector& operator*=( double a )
        {
            assert( impl_ );
            impl_->multiply_double( std::move( a ) );
            return *this;
        }

        Vector operator-() const
        {
            assert( impl_ );
            return impl_->negate();
        }

        bool operator==( const Vector& y ) const
        {
            assert( impl_ );
            if(y.space().isEmbeddedIn(this->space()))
              return impl_->compare_const_Vector_ref( y );
             else 
            if(this->space().isEmbeddedIn(y.space()))
               return (y==*this);
            return false;
        }

        /// Access underlying space.
        const VectorSpace& space() const
        {
            assert( impl_ );
            return impl_->space();
        }

        ForwardIterator begin()
        {
            assert( impl_ );
            return impl_->begin();
        }

        ForwardIterator end()
        {
            assert( impl_ );
            return impl_->end();
        }

        ConstForwardIterator begin() const
        {
            assert( impl_ );
            return impl_->begin();
        }

        ConstForwardIterator end() const
        {
            assert( impl_ );
            return impl_->end();
        }

        template < class T, typename std::enable_if<
                                !std::is_same< typename std::decay< T >::type, Vector >::value &&
                                !std::is_base_of< Interface, typename std::decay< T >::type >::
                                    value >::type* = nullptr >
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
        friend IdEmbedding;
        
        void changeSpace(const VectorSpace& V)
        {
            impl_->changeSpace(V);
        }
        
        clang::type_erasure::polymorphic::SBOCOWStorage< Interface, Wrapper, 16 > impl_;
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

} // namespace Spacy
