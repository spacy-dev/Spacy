// Copyright (C) 2015 by Lars Lubkoll. All rights reserved.
// Released under the terms of the GNU General Public License version 3 or later.

#include "pathfollowing.hh"

#include <Spacy/Util/Exceptions.h>
#include <Spacy/Spaces/ProductSpace/Vector.h>
#include <Spacy/Vector.h>
#include <Spacy/Util/Log.h>
#include <Spacy/Util/Logger.h>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <utility>
#include <limits>
#include <fstream>

namespace Spacy
{




namespace
{

auto primalProjection( const Spacy::Vector& v )
{
    auto w = v;
    auto& w_ = cast_ref< ProductSpace::Vector >( w );
    w_.component( DUAL ) *= 0;
    return w;
}

    Logger<double> logLambda("pathParameter.log");
    Logger<double> logNormDx("pathNormDx.log");
    Logger<unsigned int> logIterations("pathIterations.log");

}




namespace PathFollowing
{

DEFINE_LOG_TAG( static const char* log_tag = "PathFollowing" );




SimplePath::SimplePath( Real  lambdaInit,
                        Real  lambdaMax,  Real stepSizeUpdate) :
    lambdaInit_(lambdaInit), lambdaMax_(lambdaMax), stepsizeUpdate_(stepSizeUpdate)
{

}

void SimplePath::setInnerSolver( InnerSolver solver)
{
    innerSolver_ = std::move(solver);
}


void SimplePath::setPlot(
        std::function<void(const Vector & x, Real lambda,unsigned int step)> f)
{
    plot_ = std::move(f);
}

Vector SimplePath::solve(const Vector & x0) const
{
    LOG_INFO(log_tag, "Starting iteration.");

    startTimer();
    Real lambda = lambdaInit_;
    Real normDx = 0.0;

    auto x = x0;
    auto x_next = x;

    bool converged = false;
    const int maxIt = getMaxSteps();

    unsigned int iterations = 0;
     sumOfIterations = 0;

    testSetup(x);

    for (int step = 1; step <= maxIt; step++)
    {

        lambda += lambda*stepsizeUpdate_;
        lambda = std::min(lambdaMax_, lambda);

        std::tie(converged, x_next, iterations) =
                innerSolver_(x, lambda);

        if(converged)
        {
            auto dx = x_next -x;
            normDx = norm(primalProjection(dx));

            x = x_next;
            logLambda(get(lambda));
            logIterations(iterations);
            logNormDx(get(normDx));

            plot_(x,lambda, step);
            sumOfIterations += iterations;

            LOG_SEPARATOR(log_tag);
            LOG(log_tag, "Iteration", step);
            LOG(log_tag, "lambda", lambda);
            LOG(log_tag, "|dx|", normDx);
            LOG(log_tag, "|x|", norm(x));
            LOG_SEPARATOR(log_tag);

        }

        else throw Exception::NotConverged("Inner Solver Did Not Converge");


        converged = testResult(converged, step, lambda,  x);

        LOG_SEPARATOR(log_tag);

        if(converged)
        {
            printElapsedTime();

            std::cout << "Time in hours: " << ((elapsedTime()/1000.0)/3600.0) << std::endl;
            LOG_SEPARATOR(log_tag);
            LOG(log_tag, "PathFollowing converged:" , "True");
            LOG(log_tag, "Number of inner iterations:", sumOfIterations);
            return x;
        }
    }

    throw Exception::NotConverged("Max number of steps reached");
    return x;
}


void SimplePath::testSetup( Vector & x) const
{

    auto xTemp = x;

    if(lambdaInit_ <= 0.0)
        Exception::NotConverged("Negative initial path parameter");

    bool converged = false;

    Real normDx = 0.0;
    unsigned int iterations = 0;

    std::tie(converged, x, iterations) =
            innerSolver_(x, lambdaInit_);

    sumOfIterations += iterations;

    normDx = norm(primalProjection(x-xTemp));
    //std::cout << "Called Composite Step Method: " << std::endl;

    if(!converged)
    {
        throw Exception::NotConverged("Computation of initial iterate for path following");
    }

    LOG_SEPARATOR(log_tag);
    LOG(log_tag, "Inner solver converged for initial guess", "True");
    LOG(log_tag, "lambdaInit", lambdaInit_);
    LOG(log_tag, "|x|", norm(x));
    LOG_SEPARATOR(log_tag);

    logLambda(get(lambdaInit_));
    logIterations(iterations);
    logNormDx(get(normDx));
    plot_(x, lambdaInit_, 0);
}


bool SimplePath::testResult(bool converged, int step, Real lambda, const Vector & x) const
{

    if (converged && lambda == lambdaMax_ )
    {
        std::cout << "Path-Following Converged: " << std::endl;
        //plot_(x, step);
        result_ = Result::Converged;

        return true;
    }

    else if (!converged)
        throw Exception::NotConverged("Computation of initial iterate for path following");

    else
        return false;

}


}
}
