#include <Spacy/Adapter/Eigen.h>
#include <Spacy/Algorithm/ODEInt.h>
#include <Spacy/Spaces/ScalarSpace/DynamicOperator.h>
#include <Spacy/Spaces/ScalarSpace/Real.h>
#include <Spacy/Util/Cast.h>
#include <Spacy/ZeroVectorCreator.h>

#include <Test/gtest.hh>

#include <cmath>
#include <iostream>
#include <iomanip>

using testing::DoubleNear;

using Spacy::cast_ref;
using Spacy::Real;
using Spacy::Vector;
using Spacy::Scalar::DynamicOperator;
using Spacy::Algorithm::odeint::integrate;

const auto pi = 4 * std::atan( 1 );
const auto tol = 1e-6;

namespace
{
    struct ExplicitEuler
    {
        using stepper_category = boost::numeric::odeint::stepper_tag;

        template < class System, class State, class Time >
        void do_step( System sys, State& inout, Time t, Time dt )
        {
            State tmp( inout );
            sys( inout, tmp, t );
            inout += dt * tmp;
        }

        template < class System, class State, class Time >
        void do_step( System sys, const State& in, Time t, State& out, Time dt )
        {
            State tmp( in );
            sys( in, tmp, t );
            out += dt * tmp;
        }
    };
}

TEST( ODEIntTest, CosinusWithObserver )
{
    const auto A = DynamicOperator( []( double t, double /*x*/ ) { return -std::sin( t ); } );

    Vector x = Real( 1 );

    int steps = -1;
    x = integrate( A, x, 0, pi, 0.1, [&steps]( const Vector& y, double t ) {
        auto& x_ = cast_ref< Real >( y );
        std::cout << t << "\t" << x_ << std::endl;
        ++steps;
    } );

    EXPECT_THAT( get( cast_ref< Real >( x ) ), DoubleNear( -1, steps * tol ) );
}

TEST( ODEIntTest, Exp )
{
    const auto A = DynamicOperator( []( double t, double /*x*/ ) { return std::exp( t ); } );

    Vector x = Real( 1 );

    int steps = -1;
    x = integrate( A, x, 0, 5, 0.1, [&steps]( const auto&, double ) { ++steps; } );

    EXPECT_THAT( get( cast_ref< Real >( x ) ), DoubleNear( 148.41316, 148.41316 * steps * tol ) );
}

TEST( ODEIntTest, CosinusWithObserverAndExplicitStepper )
{
    const auto A = DynamicOperator( []( double t, double /*x*/ ) { return -std::sin( t ); } );

    Vector x = Real( 1 );

    int steps = -1;
    x = integrate( ExplicitEuler(), A, x, 0, pi, 0.01,
                   [&steps]( const Vector&, double ) { ++steps; } );

    EXPECT_THAT( get( cast_ref< Real >( x ) ), DoubleNear( -1, steps * tol ) );
}

TEST( ODEIntTest, SinusCosinus )
{
    auto V = Spacy::Rn::makeHilbertSpace( 2 );
    const auto A = [&V]( double, Vector x ) {
        auto& x_ = cast_ref< Spacy::Rn::Vector >( x ).get();
        auto tmp = x_[ 1 ];
        x_[ 1 ] = -x_[ 0 ];
        x_[ 0 ] = tmp;
        return x;
    };

    Vector x = zero( V );
    get( cast_ref< Spacy::Rn::Vector >( x ) )[ 0 ] = 0;
    get( cast_ref< Spacy::Rn::Vector >( x ) )[ 1 ] = 1;

    int steps = -1;
    x = integrate( A, x, 0, pi, 0.1, [&steps]( const auto&, double ) { ++steps; } );

    EXPECT_THAT( get( cast_ref< Spacy::Rn::Vector >( x ) )[ 0 ], DoubleNear( 0, steps * tol ) );
    EXPECT_THAT( get( cast_ref< Spacy::Rn::Vector >( x ) )[ 1 ], DoubleNear( -1, steps * tol ) );
}
