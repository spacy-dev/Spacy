#include <Test/gtest.hh>

#include <Spacy/Spaces/ScalarSpace/LinearSolver.h>
#include <Spacy/Spaces/ScalarSpace/Real.h>
#include <Spacy/Util/Cast.h>
#include <Spacy/vector.hh>

using namespace Spacy;

TEST( ScalarAdapter_LinearSolver, Apply )
{
    auto solver = Scalar::LinearSolver( 2 );
    auto real = Real{3.};
    ASSERT_EQ( get( cast_ref< Real >( solver( real ) ) ), 1.5 );
}
