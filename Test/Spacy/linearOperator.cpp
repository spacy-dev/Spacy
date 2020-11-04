#include <Test/Mock/linearOperator.hh>
#include <Test/Mock/linearSolver.hh>
#include <Test/mockSetup.hh>

#include <Spacy/LinearOperator.h>
#include <Spacy/Spaces/RealSpace.h>
#include <Spacy/Util/Cast.h>
#include <Spacy/Vector.h>

#include <gtest/gtest.h>

using namespace Spacy;

namespace
{
    void test( const LinearOperator& f, double value )
    {
        EXPECT_DOUBLE_EQ( toDouble( f( zero( Space::R ) ) ), value );
        EXPECT_EQ( 0u, f.domain().index() );
        EXPECT_EQ( 0u, f.range().index() );

        auto L = f.solver();
        auto expected = is< Mock::IndefiniteLinearSolver >( L );
        EXPECT_TRUE( expected );
    }
} // namespace

#ifndef NDEBUG
TEST( LinearOperator, Assert )
{
    auto X = createMockBanachSpace();
    LinearOperator f;
    EXPECT_DEATH( f( zero( X ) ), "" );
    EXPECT_DEATH( f.solver(), "" );
    EXPECT_DEATH( f.domain(), "" );
    EXPECT_DEATH( f.range(), "" );
    EXPECT_DEATH( -f, "" );

    LinearOperator g;
    EXPECT_DEATH( f *= 1, "" );
    EXPECT_DEATH( f += g, "" );
    EXPECT_DEATH( f -= g, "" );
    EXPECT_DEATH( f( g ), "" );

    f = Mock::LinearOperator( 1 );
    EXPECT_DEATH( f += g, "" );
    EXPECT_DEATH( f -= g, "" );
    EXPECT_DEATH( f( g ), "" );

    EXPECT_DEATH( g += f, "" );
    EXPECT_DEATH( g -= f, "" );
    EXPECT_DEATH( g( f ), "" );
}
#endif

TEST( LinearOperator, IsEmpty )
{
    LinearOperator f;
    LinearOperator g = Mock::LinearOperator( 1 );

    bool f_is_empty = !f;
    bool g_is_empty = !g;
    ASSERT_TRUE( f_is_empty );
    ASSERT_FALSE( g_is_empty );
}

TEST( LinearOperator, Cast )
{
    LinearOperator f = Mock::LinearOperator();

    ASSERT_TRUE( f.target< Mock::LinearOperator >() != nullptr );
}

TEST( LinearOperator, StoreRValue )
{
    LinearOperator f = Mock::LinearOperator( 1 );
    test( f, 1 );
}

TEST( LinearOperator, Copy )
{
    LinearOperator g = Mock::LinearOperator( 1 );
    LinearOperator f = g;
    test( f, 1 );
}

TEST( LinearOperator, StoreCopy )
{
    auto g = Mock::LinearOperator( 1 );
    LinearOperator f = g;
    test( f, 1 );
}

TEST( LinearOperator, Add )
{
    LinearOperator f = Mock::LinearOperator( 1 );
    LinearOperator g = Mock::LinearOperator( 2 );
    f += g;
    test( f, 3 );
}

TEST( LinearOperator, Subtract )
{
    LinearOperator f = Mock::LinearOperator( 1 );
    LinearOperator g = Mock::LinearOperator( 2 );
    f -= g;
    test( f, -1 );
}

TEST( LinearOperator, Multiply )
{
    LinearOperator f = Mock::LinearOperator( 1 );
    f *= 2;
    test( f, 2 );
}

TEST( LinearOperator, Negate )
{
    LinearOperator g = Mock::LinearOperator( 1 );
    auto f = -g;
    test( f, -1 );
}

TEST( LinearOperator, Move )
{
    LinearOperator g = Mock::LinearOperator( 1 );
    bool is_empty_before_move = !g;
    LinearOperator f = std::move( g );

    EXPECT_FALSE( is_empty_before_move );

    test( f, 1 );
}
