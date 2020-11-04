#include <Test/mockSetup.hh>

#include <Spacy/Spaces/ScalarSpace/LinearOperator.h>

#include <gtest/gtest.h>

using namespace Spacy;

TEST( ScalarAdapter_LinearOperator, Apply )
{
    auto space = createMockBanachSpace();
    auto A = Scalar::LinearOperator( space, 2 );

    ASSERT_EQ( toDouble( A( Real( 3 ) ) ), 6. );
}

TEST( ScalarAdapter_LinearOperator, ApplySolver )
{
    auto space = createMockBanachSpace();
    auto A = Scalar::LinearOperator( space, 2 );

    ASSERT_EQ( toDouble( A.solver()( Real( 3 ) ) ), 1.5 );
}
