// Copyright (C) 2017 by Stoecklein Matthias. All rights reserved.
// Released under the terms of the GNU General Public License version 3 or later.

#pragma once

#include <memory>
#include <string>
#include <chrono>

#include "Spacy/C2Functional.h"
#include "Spacy/Operator.h"
#include "Spacy/Util/Mixins/Accuracy.h"
#include "Spacy/Util/Mixins/Eps.h"
#include "Spacy/Util/Mixins/MaxSteps.h"
#include "Spacy/Util/Mixins/Timer.h"
#include "Spacy/Util/Mixins/Verbosity.h"
#include <Spacy/Algorithm/DampingFactor.h>
#include <Spacy/Util/Exceptions.h>

namespace Spacy
{
/// @cond
class Vector;
/// @endcond

namespace PathFollowing
{






/**
 * @ingroup PathFollowingGroup
 * @brief Simple Path-Following method.
 *
 * This implements a simple Path-Following method to solve a parameter-dependent system \f$ F(x,\lambda) = 0 \f$, for scalar parameters.
 */
class SimplePath:
        public Mixin::AbsoluteAccuracy,
        public Mixin::RelativeAccuracy,
        public Mixin::Eps,
        public Mixin::Verbosity,
        public Mixin::MaxSteps,
        public Mixin::Timer< std::chrono::milliseconds >
{
    enum class Result
    {
        Converged, Failed
    };

    using InnerSolver = std::function<
    std::tuple<bool, Vector, unsigned int>(const Vector & x, const Real & lambda)
    >;

public:

    /**
         * \brief Set up ClassicalContinuation method.
         *
     * \param lambdaInit initial path parameter value
     * \param lambdaMax  maximal path parameter value
         * \param initialStepSize  initial stepsize to update parameter \f$ \lambda \f$
         * \param minStepSize  min stepsize to update parameter \f$ \lambda \f$
         * \param thetaMin
         */
    SimplePath( Real lambdaInit, Real lambdaMax, Real stepsizeUpdate_);


    /// Set the solution function for the Newton system
    void setInnerSolver(InnerSolver solver);

    /// Set the plot function
    void setPlot(std::function<void(const Vector & x, Real lambda, unsigned int step)>);

    /**
     * \brief Solve problem via path-following method.
     * @param x0 initial guess
     */
    Vector solve(const Vector & x0) const;


private:

    /// Test feasibility of initial parameters and if the first solution on the path can be computed
    void testSetup(Vector & x) const;

    /// Test if result is reasonable
    bool testResult(bool converged, int step, Real lambda, const Vector & x) const;




    /// Solve the parameter depend system for some parameter \f$x(\lambda_{k})\f$
    std::function<
    std::tuple<bool, Vector, int>(const Vector & x, const Real & lambda)> innerSolver_ =[](const Vector & x, const Real & lambda)
    {
        throw Exception::CallOfUndefinedFunction("Inner solver has not been set for pathfollowing");
        return std::make_tuple(false,x,0);
    };

    /// Plot the solution \f$x(\lambda_{k})\f$ in each step
    std::function<void(const Vector & x, Real lambda, unsigned int step)> plot_ =
            [](const Vector & x, Real lambda, unsigned int step){ };

    mutable Result result_ = Result::Failed;

    // parameters for the classical continuation method
    const Real lambdaInit_ = 1.0;
    const Real lambdaMax_ = 10000000;


    const Real stepsizeUpdate_ = 9.0;

    mutable unsigned int sumOfIterations = 0;
};




}
}

