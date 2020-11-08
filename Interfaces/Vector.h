#pragma once

#include <functional>

#include <Spacy/Adaptivity/SpaceManager.h>
#include <Spacy/ForwardIterator.h>
#include <Spacy/Spaces/ScalarSpace/Real.h>
#include <Spacy/Util/Exceptions.h>
#include <Spacy/Util/Mixins/Get.h>
#include <Spacy/VectorSpace.h>

#include <type_traits>

namespace Spacy
{
    /// Type-erased vector.
    class Vector
    {
    public:
        /// Apply as dual space element.
        Real operator()( const Vector& x ) const;

        Vector& operator+=( const Vector& y );

        Vector& operator-=( const Vector& y );

        Vector& operator*=( double a );

        Vector operator-() const;

        bool operator==( const Vector& y ) const;

        /// Access underlying space.
        const VectorSpace& space() const;

        ForwardIterator begin();
        ForwardIterator end();

        ConstForwardIterator begin() const;
        ConstForwardIterator end() const;
    };

    /// Multiplication with arithmetic types (double,float,int,...).
    template < class Arithmetic, class = std::enable_if_t< std::is_arithmetic< Arithmetic >::value > >
    [[nodiscard]] Vector operator*( Arithmetic a, Vector x )
    {
        return x *= a;
    }

    /// Multiplication with arithmetic types (double,float,int,...).
    template < class Arithmetic, class = std::enable_if_t< std::is_arithmetic< Arithmetic >::value > >
    [[nodiscard]] Vector operator*( Vector x, Arithmetic a )
    {
        return x *= a;
    }

    /// Sum of vectors \f$z=x+y\f$.
    [[nodiscard]] inline Vector operator+( Vector x, const Vector& y )
    {
        return x += y;
    }

    /// Subtract vectors \f$z=x-y\f$.
    [[nodiscard]] inline Vector operator-( Vector x, const Vector& y )
    {
        return x -= y;
    }

    /// Compute scalar product \f$z=x*y=(x,y)\f$.
    [[nodiscard]] inline Real operator*( const Vector& x, const Vector& y )
    {
        return x.space().scalarProduct()( x, y );
    }

    /// Compute norm, where the norm associated with the underlying function space is used \f$ z =
    /// \|x\| \f$.
    [[nodiscard]] inline Real norm( const Vector& x )
    {
        return x.space().norm()( x );
    }

    [[nodiscard]] inline void checkDualPairing( const Vector& x, const Vector& y )
    {
        if ( !y.space().isPrimalWRT( x.space() ) )
            throw Exception::IncompatibleSpace( x.space().index(), x.space().name(), y.space().index(), y.space().name() );
    }

    template < class T, typename std::enable_if< std::is_arithmetic< T >::value >::type* = nullptr >
    [[nodiscard]] Vector operator*( const Mixin::Get< T >& x, Vector y )
    {
        return y *= get( x );
    }

    template < class T, typename std::enable_if< std::is_arithmetic< T >::value >::type* = nullptr >
    [[nodiscard]] Vector operator*( const Vector& x, const Mixin::Get< T >& y )
    {
        return y * x;
    }

    [[nodiscard]] inline Vector operator*( const Mixin::Get< Real >& x, Vector y )
    {
        return y *= get( get( x ) );
    }

    [[nodiscard]] inline Vector operator*( const Vector& x, const Mixin::Get< Real >& y )
    {
        return y * x;
    }
} // namespace Spacy
