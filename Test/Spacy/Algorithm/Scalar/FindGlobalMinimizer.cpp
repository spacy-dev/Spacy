#include <Spacy/Algorithm/Scalar/FindGlobalMinimizer.h>

#include <gtest/gtest.h>

using namespace Spacy;

TEST( FindGlobalMinimizer, ScalarQuadraticInside )
{
    const auto a = Real( -10 );
    const auto b = Real( 10 );
    const auto eps = Real( 1e-2 );
    auto f = []( Real t ) { return ( t - 1 ) * ( t - 1 ); };

    auto t0 = Scalar::findGlobalMinimizer( f, a, b, eps );
    EXPECT_NEAR( get( t0 ), 1, get( ( b - a ) * eps ) );
}

TEST( FindGlobalMinimizer, ScalarQuadraticOnBoundary )
{
    const auto a = Real( -1 );
    const auto b = Real( 2 );
    const auto eps = Real( 1e-2 );
    auto f = []( Real t ) { return -( t - 1 ) * ( t - 1 ); };

    auto t0 = Scalar::findGlobalMinimizer( f, a, b, eps );
    EXPECT_NEAR( get( t0 ), get( a ), get( ( b - a ) * eps ) );
}

TEST( FindGlobalMinimizer, DiscontinuousFourthOrderTwoMinima )
{
    const auto a = Real( -2 );
    const auto b = Real( 2 );
    const auto eps = Real( 1e-2 );

    auto f = []( Real t ) { return ( t - 1 ) * ( t - 1 ) * ( t + 1 ) * ( t + 1 ) - ( t > 0 ? 1 : 0 ); };

    auto t0 = Scalar::findGlobalMinimizer( f, a, b, eps );
    EXPECT_NEAR( get( t0 ), 1, get( ( b - a ) * eps ) );

    auto g = []( Real t ) { return ( t - 1 ) * ( t - 1 ) * ( t + 1 ) * ( t + 1 ) + ( t > 0 ? 1 : 0 ); };

    t0 = Scalar::findGlobalMinimizer( g, a, b, eps );
    EXPECT_NEAR( get( t0 ), -1, get( ( b - a ) * eps ) );
}
