#include <Spacy/Operator.h>
#include <Spacy/Spaces/RealSpace.h>
#include <Spacy/Util/Cast.h>
#include <Spacy/Vector.h>
#include <Spacy/ZeroVectorCreator.h>

#include <gtest/gtest.h>

TEST( ScalarSpaceTest, DefaultIndex )
{
    auto r = Spacy::Real{ 2. };
    ASSERT_EQ( r.space().index(), 0 );
}

TEST( ScalarSpaceTest, ElementTest )
{
    using namespace Spacy;
    auto R = makeRealSpace();
    auto x = zero( R );
    EXPECT_DOUBLE_EQ( get( cast_ref< Real >( x ) ), 0. );
}

TEST( ScalarSpaceTest, ScalarProductTest )
{
    using namespace Spacy;
    auto R = makeRealSpace();
    auto x = zero( R );
    auto y = zero( R );
    get( cast_ref< Real >( x ) ) = 1;
    get( cast_ref< Real >( y ) ) = -2;
    EXPECT_DOUBLE_EQ( get( cast_ref< Real >( x ) ), 1. );
    EXPECT_DOUBLE_EQ( get( cast_ref< Real >( y ) ), -2. );
    EXPECT_DOUBLE_EQ( get( x * y ), -2. );
    EXPECT_DOUBLE_EQ( get( x * y ), get( R.scalarProduct()( x, y ) ) );
    EXPECT_TRUE( R.isHilbertSpace() );
}

TEST( ScalarSpaceTest, NormTest )
{
    using namespace Spacy;
    auto R = makeRealSpace();
    auto x = zero( R );
    auto y = zero( R );
    get( cast_ref< Real >( x ) ) = 1;
    get( cast_ref< Real >( y ) ) = -2;
    EXPECT_DOUBLE_EQ( get( cast_ref< Real >( x ) ), 1. );
    EXPECT_DOUBLE_EQ( get( cast_ref< Real >( y ) ), -2. );
    EXPECT_DOUBLE_EQ( get( R.norm()( x ) ), 1. );
    EXPECT_DOUBLE_EQ( get( R.norm()( y ) ), 2. );
}
