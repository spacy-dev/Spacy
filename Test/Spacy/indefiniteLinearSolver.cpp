#include <Test/gtest.hh>

#include <Test/Mock/linearSolver.hh>
#include <Test/Mock/norm.hh>
#include <Test/Mock/vector.hh>
#include <Test/Mock/vectorCreator.hh>
#include <Spacy/Spaces/ScalarSpace/Real.h>
#include <Spacy/linearSolver.hh>
#include <Spacy/vectorSpace.hh>
#include <Spacy/zeroVectorCreator.hh>

using namespace Spacy;

TEST( IndefiniteLinearSolver, Construction )
{
    auto solver = Mock::IndefiniteLinearSolver{};
    IndefiniteLinearSolver typeErasedSolver( solver );
    ASSERT_FALSE( !typeErasedSolver );
}

TEST( IndefiniteLinearSolver, Assignment )
{
    auto solver = Mock::IndefiniteLinearSolver{};
    IndefiniteLinearSolver typeErasedSolver;
    ASSERT_TRUE( !typeErasedSolver );
    typeErasedSolver = solver;
    ASSERT_FALSE( !typeErasedSolver );
}

TEST( IndefiniteLinearSolver, Apply )
{
    auto solver = Mock::IndefiniteLinearSolver{};
    IndefiniteLinearSolver typeErasedSolver( solver );
    const auto V = VectorSpace(Mock::VectorCreator(), Mock::Norm());
    ASSERT_NO_THROW( typeErasedSolver( zero(V) ) );
}

TEST( IndefiniteLinearSolver, IsPositiveDefinite )
{
    auto solver = Mock::IndefiniteLinearSolver{};
    IndefiniteLinearSolver typeErasedSolver( solver );
    ASSERT_NO_THROW( typeErasedSolver.isPositiveDefinite() );
    ASSERT_TRUE( typeErasedSolver.isPositiveDefinite() );
}
