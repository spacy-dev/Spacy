#include <Spacy/Algorithm/ODEInt.h>
#include <Spacy/Spaces/ScalarSpace/DynamicOperator.h>
#include <Spacy/Spaces/ScalarSpace/Real.h>

#include <Test/gtest.hh>

#include <cmath>
#include <iostream>

using testing::DoubleNear;

using Spacy::cast_ref;
using Spacy::Real;
using Spacy::Vector;
using Spacy::Scalar::DynamicOperator;
using Spacy::Algorithm::odeint::integrate;
using Spacy::Algorithm::odeint::ExplicitEuler;

const auto pi = 4*std::atan(1);
const auto tol = 1e-6;

TEST(ODEIntTest, CosinusWithObserver)
{
    const auto A = DynamicOperator([](double t, double /*x*/){ return -std::sin(t); });

    Vector x = Real(1);

    int steps = -1;
    x = integrate(A, x, 0, pi, 0.1,
                  [&steps](const auto& x, double t)
    {
        std::cout << t << "\t" << x[0] << std::endl;
        ++steps;
    });

    EXPECT_THAT(get(cast_ref<Real>(x)), DoubleNear(-1, steps*tol));
}

TEST(ODEIntTest, Exp)
{
    const auto A = DynamicOperator([](double t, double /*x*/){ return std::exp(t); });

    Vector x = Real(1);

    int steps = -1;
    x = integrate(A, x, 0, 5, 0.1, [&steps](const auto&, double){ ++steps; });

    EXPECT_THAT(get(cast_ref<Real>(x)), DoubleNear(148.41316, 148.41316*steps*tol));
}

TEST(ODEIntTest, CosinusWithObserverAndExplicitStepper)
{
    const auto A = DynamicOperator([](double t, double /*x*/){ return -std::sin(t); });

    Vector x = Real(1);

    int steps = -1;
    x = integrate(ExplicitEuler(), A, x, 0, pi, 0.01,
                  [&steps](const Vector&, double){ ++steps; });

    EXPECT_THAT(get(cast_ref<Real>(x)), DoubleNear(-1, steps*tol));
}


