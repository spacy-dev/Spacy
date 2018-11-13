#include <Test/gtest.hh>

#include <Spacy/Algorithm/Scalar/Models.h>

TEST( Algorithm_Scalar_Quadratic, CreateAndApply )
{
    auto f = Spacy::Scalar::Quadratic( 3, 2, 1 );

    EXPECT_EQ( f( 2 ), 11 );
}

TEST( Algorithm_Scalar_Cubic, CreateAndApply )
{
    auto f = Spacy::Scalar::Cubic( 4, 3, 2, 1 );

    EXPECT_EQ( f( 2 ), 26 );
}
