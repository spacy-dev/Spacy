#include <Test/Spacy/Algorithm/Newton/newtonTestSetup.hh>
#include <Test/Spacy/Algorithm/scalarTests.hh>

#include <Spacy/Spaces/RealSpace.h>

#include <gtest/gtest.h>

using namespace Spacy;

namespace LocalNewton
{
    inline auto createSolver( C1Operator F, double relativeAccuracy, double eps, unsigned maxSteps )
    {
        return createLocalNewton( F, relativeAccuracy, eps, maxSteps );
    }

    GENERATE_SCALAR_TEST( LocalNewton, Linear, 2 )
    GENERATE_SCALAR_TEST( LocalNewton, Quadratic, 7 )
} // namespace LocalNewton

namespace CovariantNewton
{
    inline auto createSolver( C1Operator F, double relativeAccuracy, double eps, unsigned maxSteps )
    {
        return createCovariantNewton( F, relativeAccuracy, eps, maxSteps );
    }

    GENERATE_SCALAR_TEST( CovariantNewton, Linear, 2 )
    GENERATE_SCALAR_TEST( CovariantNewton, Quadratic, 7 )
} // namespace CovariantNewton

namespace ContravariantNewton
{
    inline auto createSolver( C1Operator F, double relativeAccuracy, double eps, unsigned maxSteps )
    {
        return createContravariantNewton( F, relativeAccuracy, eps, maxSteps );
    }

    GENERATE_SCALAR_TEST( ContravariantNewton, Linear, 1 )
    GENERATE_SCALAR_TEST( ContravariantNewton, Quadratic, 6 )
} // namespace ContravariantNewton

#include "Test/Spacy/Algorithm/undefScalarTests.hh"
