#include "ACR.h"

#include <Spacy/Operator.h>
#include <Spacy/Algorithm/DampingFactor.h>
#include <Spacy/Algorithm/CG/LinearSolver.h>
#include <Spacy/Algorithm/CompositeStep/QuadraticModel.h>

#include <Spacy/Algorithm/CG/CG.h>
#include <Spacy/Algorithm/CG/Regularization.h>
#include <Spacy/Algorithm/CG/TerminationCriteria.h>

#include <Spacy/Algorithm/Scalar/FindGlobalMinimizer.h>
#include <Spacy/Algorithm/CG/LinearSolver.h>
#include <Spacy/Algorithm/CG/RegularizeViaPreconditioner.h>

#include <Spacy/ZeroVectorCreator.h>
#include <Spacy/Util/Cast.h>
#include <Spacy/Util/Log.h>
#include "Spacy/Util/Exceptions.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <utility>
#include <Spacy/InducedScalarProduct.h>
#include <Spacy/Util/Logger.h>



namespace Spacy
{


Logger< double > logLambda( "acrLambda.log" );
Logger< double > logOmega( "acrOmega.log" );
Logger< int > logRejected( "acrRejected.log" );
Logger< double > logCostFunctional( "acrFunctional.log" );
Logger< double > logFunctionalDecrease( "acrFunctionalDecrease.log" );
Logger< bool > logConvexity( "acrConvex.log" );




namespace ACR
{
DEFINE_LOG_TAG(static const char* log_tag = "ACR");



MyCubicModel::MyCubicModel(Real a0, Real a1, Real a2, Real a3): a0_(a0), a1_(a1), a2_(a2), a3_(a3)
{

}

Real MyCubicModel::operator ()(Real x)
{
    update(x);
    return this->operator()();
}

Real MyCubicModel::operator()() const
{
    return   a0_ + a1_*x_ + a2_*x_*x_ +a3_*x_*x_*x_;
}

void MyCubicModel::update(Real x)
{
    x_ = x;
}

// use routines that are already there
MyCubicModel makeMyCubicModel(C2Functional  f,const Vector & x,const Vector & dx, Real normDx, Real omega)
{
    return MyCubicModel(0.0,f.d1(x)(dx),0.5*f.d2(x,dx)(dx),(1.0/6.0)*omega*normDx*normDx*normDx);
}

MyCubicModel makeMyCubicModelFunction(C2Functional  f,const Vector & x,const Vector & dx, Real normDx, Real omega)
{
    return MyCubicModel(f(x),f.d1(x)(dx),0.5*f.d2(x,dx)(dx),(1.0/6.0)*omega*normDx*normDx*normDx);
}




CompositeStep::CubicModel makeCubicModel(const Vector& dx,
                                         const C2Functional& f, const Vector& x, Spacy::Real omega)
{
    return CompositeStep::CubicModel(
                CompositeStep::makeQuadraticModel(Spacy::DampingFactor(0.0), dx, dx,
                                                  f, x),
                CompositeStep::makeQuadraticNormModel(Spacy::DampingFactor(0.0), dx,
                                                      dx), omega);
}

ACRSolver::ACRSolver(C2Functional f, VectorSpace &domain, VectorSpace &tangentSpace, double omega,
                     double relativeAccuracy) :
    Mixin::RelativeAccuracy(relativeAccuracy), f_(std::move(f)), domain_(domain), tangentSpace_(tangentSpace), omega_(omega)
{
}

Vector ACRSolver::operator()()
{
    return operator()(zero(domain_));
}

Vector ACRSolver::operator()(const Vector& x0)
{
    LOG_INFO(log_tag, "Starting ACR-Solver.");

    auto x = x0;
    auto dx = zero(tangentSpace_);
    Real dxQNorm = 10.0;
    Real lambda = 1.0;

    plot_(x0);
    // Not always a good choice
    domain_.setScalarProduct(::Spacy::PrimalInducedScalarProduct(f_.hessian(zero(domain_))));

    // reset after regularization
    auto f_temp_ = f_;

    for (unsigned step = 1; step <= getMaxSteps(); ++step)
    {

        LOG_SEPARATOR(log_tag,"");
        LOG(log_tag, "Iteration", step);
        stepMonitor = StepMonitor::Accepted;

        energyRegularization_ = 0.0;


        
        std::tie(dx,dxQNorm) = computeStep(x,dxQNorm);

        

        LOG(log_tag, "f(x)", f_(x));
        LOG(log_tag, "TestDirection of Descent: f_.d1(x)(dx): ", f_.d1(x)(dx));


        LOG(log_tag, "||dx|| ", dxQNorm);
        //  std::cout << "DxQNormSquared: " << dxQNorm << std::endl;

        if(isConvex_ && convexProjection_ && step > 1 )
        {
            //Todo test at the end and make at least one energy decreasing step
            std::cout << "Convex projection successful" << std::endl;
            return x;
        }

        if (dxQNorm < epsilon_ && isConvex_)
        {
            LOG_INFO(log_tag, "ACR converged");
            return x +dx;
        }

        isConvex_ = true;

        do
        {
            //auto cubicModel = makeCubicModel(dx, f_, x, omega_);
            auto cubicModel = makeMyCubicModelFunction( f_, x, dx,  norm(dx), omega_);

            LOG_INFO(log_tag, "Computing damping factor");

            lambda = Scalar::findLogGlobalMinimizer(cubicModel, 1e-12, 100, 0.01);

            if(lambda >= 0.95 && lambda <= 1.05 )
                lambda = 1.0;



            auto adxlin = get(lambda) * dx;

            auto flin = f_(x + adxlin);

            LOG(log_tag, "FunctionValueBeforeNonlinearUpdate: f_(x+dx): ", flin);

            
            auto adx = get(lambda)*dx;
            
            nonlinUpdate_(x, adx, step);

            auto fnonlin = f_(x + adx);

            LOG(log_tag, "FunctionValueAfterNonlinearUpdate: f_(x+dx): ", fnonlin);


            if(flin < fnonlin) adx = adxlin;
            LOG(log_tag, "lambda: ", lambda, " omega: ", omega_, "lambda |dx|: " , norm(adx), " cubicModel: ", cubicModel(lambda));

            if ((stepMonitor = acceptanceTest(x, adx, lambda, cubicModel))
                    == StepMonitor::Accepted)
            {
                LOG_INFO(log_tag, "Accepted");
                logCostFunctional( get(f_(x+adx)));
                logOmega( get(omega_) );
                logLambda( get(lambda ));
                logFunctionalDecrease(get(f_(x) -f_(x+adx)) );
              std::cout << x.space().name() << " " << adx.space().name() << std::endl;
                x += adx;
               std::cout << x.space().name() << std::endl;
               
               auto y=zero(x.space());
               std::cout << y.space().name() << std::endl;
               y=adx;
               std::cout << y.space().name() << std::endl;
                output_(x,adx);
                plot_(x);
            }

            else
                LOG_INFO(log_tag, "Rejected");

            // Modifikation von omega
            omega_ = weightChange(omega_);


        }
        while (stepMonitor == StepMonitor::Rejected && omega_ <= omegaMax_ );
        if(stepMonitor == StepMonitor::Rejected  )
        {

           // throw Exception::NotConverged("Maximum value for the regularization parameter omega reached");
            std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Max Omega reached" << std::endl;
            return x;
        }
        LOG_SEPARATOR(log_tag,"");

        f_ = f_temp_;
    }
    // throw Exception::NotConverged("Maximum number of iterations reached");
    return x;
}







void ACRSolver::setNonlinUpdate(
        std::function<void(const Vector & x, Vector& dx, unsigned int step)> nonlinUpdate)
{
    nonlinUpdate_ = std::move(nonlinUpdate);
}

void ACRSolver::setOutput(std::function<void(const Vector &, const Vector &)> outputUpdate)
{
    output_ = std::move(outputUpdate);
}

C2Functional &  ACRSolver::getFunctional()
{
    return f_;
}

void ACRSolver::setPositivePrecond(CallableOperator p)
{
    positivePrecond_ = std::move(p);
}




ACRSolver::StepMonitor ACRSolver::acceptanceTest(const Vector &x,
                                                 const Vector &dx, const Real & lambda,
                                                 MyCubicModel &cubicModel)
{

    // Prevent dividing by zero e.g. in the case x_0 = opimal solution
    // Consider the Case numerator = nominator = 0

    const auto diffModel = cubicModel(lambda) - cubicModel(0.0);
    const auto diffFunction = (f_(x + dx) - f_(x));

    //std::cout << "m(dx)-m(0): " << diffModel << std::endl;
    //std::cout << "f_(x+dx)-f_(x): " << diffFunction << std::endl;

    if ( !std::signbit(get(diffFunction)) ||  std::isnan(get(diffFunction))  )
    {
        rho_ = -1.0;
    }

    // better use std::max

    else {
        //Change that Condition
        rho_ = diffFunction / diffModel;
        if (std::isnan(get(rho_)))
            rho_ = -1.0;

    }
    LOG(log_tag, "f_(x+dx): ", f_(x+dx), "f_(x): ", f_(x), "CubicModel(lambda): ", cubicModel(lambda), "CubicModel(0): ", cubicModel(0), "rho: ", rho_, "eta1 ", eta1_, "eta2 ", eta2_ );

    if (rho_ >= eta1_)
        return StepMonitor::Accepted;

    return StepMonitor::Rejected;
}

Spacy::Real ACRSolver::weightChange(Spacy::Real omega) const
{
    if (rho_ > eta2_)
        omega *= 0.5;
    else if (rho_ < eta1_)
        omega *= 1.2;

    return omega;
}

std::tuple<Vector,Real> ACRSolver::computeStep(const Spacy::Vector &x, Real dxNorm) const
{
    LinearSolver solver = { };
    solver = f_.hessian(x) ^ -1;

    if (is<CG::LinearSolver>(solver))
    {
        const auto& cgSolver = cast_ref<CG::LinearSolver>(solver);
        auto accuracy = getRelativeAccuracy();

        if(dxNorm <= 0.05)
            accuracy *= 1e-2;

        auto cg = makeTCGSolver(f_.hessian(x), cgSolver.P(), eps(), verbose());
        cg.setRelativeAccuracy(accuracy);
        cg.setAbsoluteAccuracy(getAbsoluteAccuracy());
        cg.set_eps(eps());
        cg.setVerbosity(2);
        cg.setTerminationCriterion(CG::Termination::AdaptiveRelativeEnergyError());

        auto result = cg(-f_.d1(x));


        // temporarily replace f_
        if(!(cg.isPositiveDefinite()))
        {
            isConvex_ = false;
            if(updateFunctional_)
            {

                std::cout << "Update Functional: " << std::endl;

                bool positiveDef = false;
                do{

                    if(energyRegularization_ == 0.0)
                        energyRegularization_ = 1e-8;
                    else
                        energyRegularization_ *= 5.0;

                    std::cout << "EnergyRegularization: " <<  energyRegularization_ << std::endl;

                    f_ = updateFunctional(energyRegularization_);

                    LinearSolver solverTemp = { };
                    solverTemp = f_.hessian(x) ^ -1;

                    if (is<CG::LinearSolver>(solverTemp) == false)
                    {
                        throw;
                    }

                    const auto& cgSolverTemp = cast_ref<CG::LinearSolver>(solverTemp);

                    auto cgTemp = makeTCGSolver(f_.hessian(x), cgSolverTemp.P(), eps(), verbose());
                    cgTemp.setRelativeAccuracy(accuracy);
                    cgTemp.setAbsoluteAccuracy(getAbsoluteAccuracy());
                    cgTemp.set_eps(eps());
                    cgTemp.setVerbosity(2);
                    cgTemp.setTerminationCriterion(CG::Termination::AdaptiveRelativeEnergyError());

                    result = cgTemp(-f_.d1(x));

                    positiveDef = cgTemp.isPositiveDefinite();


                }while(!positiveDef);

                energyRegularization_ = 0.0;
              // f_ = updateFunctional(0.0);
            }

            else
            {
                std::cout << "Start Regularized CG: " << std::endl;
                isConvex_ = false;
                auto newCg =  CG::LinearSolver( f_.hessian(x),positivePrecond_, false, ::Spacy::CG::RegularizeViaCallableOperator(f_.hessian(zero(domain_)), 0.0  ) );
                newCg.setRelativeAccuracy(getRelativeAccuracy());
                newCg.setAbsoluteAccuracy(getAbsoluteAccuracy());
                newCg.set_eps(eps());
                newCg.setVerbosity(2);
                newCg.setTerminationCriterion(CG::Termination::AdaptiveRelativeEnergyError());
                result = newCg(-f_.d1(x));
                return std::make_tuple(result,norm(result));
            }
        }

        return std::make_tuple(result,norm(result));
    }

    std::cout << "Fail: Call false CG !!!!!!!!!!!" << std::endl;
    auto tcg = makeTCGSolver(f_.hessian(x), solver);
    auto result = tcg(-f_.d1(x));

    // Check if direct solver can detect nonconvexities


    return std::make_tuple(result,norm(result));
}


//std::tuple<bool,Real,Real,Real,Real> ACRSolver::checkContraction(int step, Real normDx,Real normDx1_, Real normDx2_, Real thetaZero_, Real thetaGlobal) const
//{
//    Real normDx1 = normDx1_;
//    Real normDx2 = normDx2_;
//    Real theta = thetaZero_;
//    Real thetaZero = thetaZero_;
//    bool contracted = true;

//    if (step > 2)
//    {
//        normDx1 = normDx2;
//        normDx2 = normDx;

//        theta = normDx2 / normDx1;

//        if (std::isnan(get(theta)) || std::isinf(get(theta)))
//            throw Exception::InvalidArgument("Contraction Factor is NaN");

//        if (theta > thetaGlobal )
//        {
//            std::cout << "No convergence for inner solver in Iterate: " << step  << " "<< theta << " " << normDx1
//                      << " " << normDx2 << '\n';
//            contracted = false;
//        }
//    }

//    else if(step ==1)
//    {
//        theta = thetaZero = 0.1;
//        normDx1 = normDx;
//    }
//    else if (step == 2)
//    {
//        normDx2 = normDx;
//        theta = normDx2 / normDx1;
//        thetaZero = theta;

//        if (std::isnan(get(theta)) || std::isinf(get(theta)))
//            throw Exception::InvalidArgument("Contraction Factor is NaN");


//        if (theta > thetaGlobal )
//        {
//            std::cout << "No convergence for inner solver in Iterate: " << step  << " " << theta <<" " << normDx1
//                      << " " << normDx2 << '\n';
//            contracted = false;
//        }
//    }

//    return std::make_tuple(contracted, normDx1, normDx2, thetaZero, theta);

//}


bool ACRSolver::convergenceTest( Real dxQNorm, const Vector & x) const
{

    if (dxQNorm < epsilon_)
    {
        // std::cout << "ACR Converged: " << '\n';
        //std::cout << "DxQNorm: " << dxQNorm << '\n';
        //std::cout << "f_.d1(x)*f_.d1(x): " <<  f_.d1(x)*f_.d1(x) << '\n';
        return true;
    }

    else
        return false;
}




ElasticityACRSolver::ElasticityACRSolver(C2Functional f, VectorSpace& domain): domain_(domain), f_(std::move(f))
{

}

// ToDo clean up
Vector ElasticityACRSolver::operator()(const Vector & x0)
{
    auto x = x0;
    auto lambda = DampingFactor{1.0};
    Real normDx = 10000000000000000.0;

    LOG_INFO(log_tag, "Starting ACR-Solver.");
    logCostFunctional( get(f_(x)));

    domain_.setScalarProduct(::Spacy::PrimalInducedScalarProduct(f_.hessian(zero(domain_))));

    for (unsigned step = 1; step <= getMaxSteps(); ++step)
    {
        LOG_SEPARATOR(log_tag,"");
        LOG(log_tag, "Iteration", step);
        
        auto dx = computeStep(x,normDx);
        normDx = norm(dx);
        auto normX = norm(x);

        LOG(log_tag, "Is convex: ", isConvex_);

        if(isConvex_ && convexProjection_)
        {
            //Todo test at the end and make at least one energy decreasing step
            std::cout << "Convex projection successful" << std::endl;
            return x;
        }
        // todo set scalar product with hessian

        do
        {
            // Todo
            // inefficient to create model every time this way
            auto cubicModel = makeMyCubicModel(f_, x , dx, normDx, omega_);

            lambda = Scalar::findLogGlobalMinimizer(cubicModel, 1e-10, 1000, 0.001);

            if(lambda >= 0.95 && lambda <=1.0 )
                lambda = 1.0;

            auto lambdaDx = get(lambda) * dx;

            auto f_lin = f_(x+lambdaDx);
            auto dxLin = lambdaDx;

            // choose the best update
            nonlinUpdate_(x,lambdaDx,step);

            auto f_nonlin = f_(x+lambdaDx);

            if(f_lin < f_nonlin)
                lambdaDx = dxLin;

            LOG(log_tag, "f_lin: ", f_lin, "f_nonlin ", f_nonlin);
            LOG(log_tag, "lambda: ", lambda, " omega: ", omega_, "|dx|: " , norm(lambdaDx), " cubicModel(lambda): ", cubicModel(Real(get(lambda)))," cubicModel(0): ", cubicModel(0.0));
            LOG(log_tag, "f_(x+dx): ", f_(x+lambdaDx), "f_(x): ", f_(x));
            //  LOG(log_tag, "f_(x+dx) - f_(x): ", f_(x+lambdaDx) -f_(x));

            stepMonitor_ = acceptanceTest(x,dxLin, f_(x+lambdaDx));

            // actual update is always better or equal to the linear update
            omega_ = weightUpdate(x, dxLin);

            if(stepMonitor_  == StepMonitor::Accepted)
            {
                dx = lambdaDx;
                normDx = norm(dx);
                logCostFunctional( get(f_(x+dx)));
                logOmega( get(omega_) );
                logLambda( get(get(lambda )));
                output_(x,dx);
                plot_(x);
            }
        }
        while (stepMonitor_ == StepMonitor::Rejected );

        logFunctionalDecrease(get(f_(x) -f_(x+dx)) );
        x += dx;

        LOG(log_tag,  " omega: ", omega_, "|dx|: " , normDx, "lambda: ", lambda, "|x|: ", normX );
        LOG_SEPARATOR(log_tag,"");

        if(terminate(normX,normDx,lambda))
        {
            LOG_INFO(log_tag,  "Converged ");
            return x;
        }
    }

    return x;
}

Vector ElasticityACRSolver::operator()()
{
    return operator()(zero(domain_));
}

void ElasticityACRSolver::setNonlinUpdate(
        std::function<void(const Vector & x, Vector& dx, unsigned int step)> nonlinUpdate)
{
    nonlinUpdate_ = std::move(nonlinUpdate);
}

void ElasticityACRSolver::setPositivePrecond(CallableOperator p)
{
    positivePrecond_ = std::move(p);
}


void ElasticityACRSolver::setOutput(std::function<void(const Vector& x, const Vector & dx)> outputUpdate)
{
    output_ = std::move(outputUpdate);
}

Vector ElasticityACRSolver::computeStep(const Spacy::Vector &x, Real dxNorm)
{
    LinearSolver solver = { };
    solver = f_.hessian(x) ^ -1;

    if (is<CG::LinearSolver>(solver))
    {
        const auto& cgSolver = cast_ref<CG::LinearSolver>(solver);

        auto accuracy = getRelativeAccuracy();
        std::cout << "ACR-CG Relative Accuracy: " << getRelativeAccuracy() << std::endl;

        if(dxNorm <= 0.05)
            accuracy *= 1e-2;

        auto cg = makeTCGSolver(f_.hessian(x), cgSolver.P(), eps(), verbose());
        cg.setRelativeAccuracy(accuracy);
        cg.setVerbosity(2);
        cg.setTerminationCriterion(CG::Termination::AdaptiveRelativeEnergyError());

        auto result = cg(-f_.d1(x));

        if(!(cg.isPositiveDefinite()))
        {
            isConvex_ = false;
            std::cout << "getRelativeAccuracy(): " << getRelativeAccuracy() << std::endl;
            auto newCg =  CG::LinearSolver(f_.hessian(x),positivePrecond_, false, ::Spacy::CG::RegularizeViaCallableOperator(f_.hessian(zero(domain_)), 0.0  ) );
            newCg.setRelativeAccuracy(getRelativeAccuracy());
            newCg.setVerbosity(2);
            newCg.setTerminationCriterion(CG::Termination::AdaptiveRelativeEnergyError());
            result = newCg(-f_.d1(x));
            std::cout << "Finished CG: " << std::endl;
            return result;
        }

        isConvex_ = true;

        return result;
    }

    isConvex_ = true;
    std::cout << "Fail: Call false CG" << std::endl;
    auto tcg = makeTCGSolver(f_.hessian(x), solver);
    auto result = tcg(-f_.d1(x));

    // Check if direct solver can detect nonconvexities


    return result;
}



Real ElasticityACRSolver::weightUpdate(const Vector &x ,const Vector & dx)
{
    auto f_x_next =  f_(x+dx);
    auto normDx = norm(dx);

    if(std::isnan(get(f_x_next)) || std::isinf(get(f_x_next)))
        return omega_*2.0;

    else if(normDx < stabilityAcc_)
    {
        //        std::cout << std::endl << std::endl;
        //        std::cout << "TestWeightUpdate: " << std::endl;
        //        std::cout << "-f_.d1(x)(dx) -f_.d2(x,dx)(dx): " << -(f_.d1(x)(dx) +f_.d2(x,dx)(dx)) << std::endl;
        //        std::cout << "Directional Minimizer: " << f_.d1(x)(dx) + f_.d2(x,dx)(dx) + (omega_/2.0)*normDx*normDx*normDx << std::endl;
        //        std::cout << "(omega_/2.0)*normDx*normDx*normDx: " << (omega_/2.0)*normDx*normDx*normDx << std::endl;
        //        std::cout << std::endl << std::endl;

        return (2.0/(normDx*normDx*normDx))*fabs(get(f_.d1(x+dx)(dx)-f_.d1(x)(dx) -f_.d2(x,dx)(dx)));
    }

    else
    {
        //   std::cout << std::endl;
        //        std::cout << "Directional Minimizer: " << f_.d1(x)(dx) + f_.d2(x,dx)(dx) + (omega_/2.0)*normDx*normDx*normDx << std::endl;
        //        std::cout << "TestUpdate: " << "normDx: " << normDx << " f_x_next: " << f_x_next << " f_(x): " << f_(x) << " f_.d1(x)(dx): " << f_.d1(x)(dx) << " f_.d2(x,dx)(dx): " << f_.d2(x,dx)(dx) << std::endl;
        //        std::cout << "(6.0/(normDx*normDx*normDx): " << (6.0/(normDx*normDx*normDx)) << " -f_(x): " << -f_(x) << " -f_.d1(x)(dx): " <<  -f_.d1(x)(dx) <<  " -0.5*f_.d2(x,dx)(dx):"  <<  -0.5*f_.d2(x,dx)(dx) << std::endl;
        //        std::cout << "f_x_next -f_(x) -f_.d1(x)(dx) -0.5*f_.d2(x,dx)(dx): " << f_x_next -f_(x) -f_.d1(x)(dx) -0.5*f_.d2(x,dx)(dx) << std::endl << std::endl;
        //        std::cout << "(6.0/(normDx*normDx*normDx))*fabs(get(f_x_next -f_(x) -f_.d1(x)(dx) -0.5*f_.d2(x,dx)(dx))): " << (6.0/(normDx*normDx*normDx))*fabs(get(f_x_next -f_(x) -f_.d1(x)(dx) -0.5*f_.d2(x,dx)(dx))) << std::endl << std::endl;
        //        std::cout << "(6.0/(normDx*normDx*normDx))*(0.5*f_.d1(x)(dx) - (omega_/36.0)*normDx*normDx*normDx -f_.d1(x)(dx) - 0.5*f_.d2(x,dx)(dx)): " << (6.0/(normDx*normDx*normDx))*(0.5*f_.d1(x)(dx) - (omega_/36.0)*normDx*normDx*normDx -f_.d1(x)(dx) - 0.5*f_.d2(x,dx)(dx)) << std::endl << std::endl;
        //        std::cout << "(6.0/(normDx*normDx*normDx))*(-(omega_/36.0)*normDx*normDx*normDx + (omega_/4.0)*normDx *normDx *normDx): " << (6.0/(normDx*normDx*normDx))*(-(omega_/36.0)*normDx*normDx*normDx + (omega_/4.0)*normDx *normDx *normDx  ) << std::endl << std::endl;
        // std::cout << std::endl;
        return std::max(rhoMin_*omega_,std::min(omega_*rhoMax_,(6.0/(normDx*normDx*normDx))*fabs(get(f_x_next -f_(x) -f_.d1(x)(dx) -0.5*f_.d2(x,dx)(dx)))));
    }
}

bool ElasticityACRSolver::terminate(Real normX, Real normDx, DampingFactor lambda) const
{
    if(lambda < 0.95)
        return false;



    else
    {
        auto relativeAccuracy = getRelativeAccuracy();

        //        relativeAccuracy = 1e-3;
        //        std::cout << "IsConvex: " << isConvex_ << std::endl;
        //        std::cout << "normDx: " << normDx << std::endl;
        //        std::cout << "relativeAccuracy: " << relativeAccuracy << std::endl;
        //        std::cout << "normX: " << normX << std::endl;

        return (isConvex_ && (normDx < relativeAccuracy * normX
                              || (normX < eps() && normDx < eps())));
    }
}


ElasticityACRSolver::StepMonitor ElasticityACRSolver::acceptanceTest(const Vector& x, const Vector& dx, Real f_next)
{

    if(std::isnan(get(f_next)) || std::isinf(get(f_next)))
        return StepMonitor::Rejected;

    auto normDx = norm(dx);

    if( normDx < stabilityAcc_)
    {
        std::cout << std::endl;
        std::cout << "Acceptance test: " << std::endl;
        std::cout << "f_.d1(x+dx)(dx) " << f_.d1(x+dx)(dx) << "(omega_/6.0)*normDx*normDx*normDx " << (omega_/6.0)*normDx*normDx*normDx << std::endl;
        std::cout << std::endl;

        if(f_.d1(x+dx)(dx)< (omega_/6.0)*normDx*normDx*normDx)
        {
            return StepMonitor::Accepted;
        }
        else
            return StepMonitor::Rejected;
    }

    else{

        //std::cout << std::endl;
        //std::cout << "Acceptance test: " << std::endl;
        //std::cout << "f_next: " << f_next << " f_(x): " << f_(x) << " f_.d1(x)(dx): " << f_.d1(x)(dx) << " omega_: " << omega_ << " normDx: " << normDx << std::endl;
        //std::cout << std::endl;

        //std::cout <<  "f_(x) +0.5* f_.d1(x)(dx)- (omega_/36.0)*normDx*normDx*normDx: " << f_(x) +0.5* f_.d1(x)(dx)- (omega_/36.0)*normDx*normDx*normDx << std::endl << std::endl;

        if( f_next < f_(x) +0.5* f_.d1(x)(dx)- (omega_/36.0)*normDx*normDx*normDx )
        {
            return StepMonitor::Accepted;
        }
        else
            return StepMonitor::Rejected;
    }

    return StepMonitor::Rejected;
}




void ElasticityACRSolver::regularityTest(DampingFactor lambda) const
{
    if (!regularityTestPassed(lambda))
        throw Exception::RegularityTestFailed( "ACR::regularityTest (lambda,...)", get(get(lambda)));
}

void ElasticityACRSolver::setLogName(const std::string & logName)
{
    logName_ = logName;
}










//ElasticityACRSolver::ElasticityACRSolver(C2Functional f): domain_(f.domain()), f_(std::move(f))
//{

//}


//Vector ElasticityACRSolver::operator()(const Vector & x0)
//{
//    auto x = x0;
//    auto lambda = DampingFactor{1.0};

//    LOG_INFO(log_tag, "Starting ACR-Solver.");

//    for (unsigned step = 1; step <= getMaxSteps(); ++step)
//    {
//        LOG_SEPARATOR(log_tag);
//        LOG(log_tag, "Iteration", step);

//        auto dx = computeStep(x);
//        auto normDx = norm(dx);
//        auto normX = norm(x);

//        LOG(log_tag, "Is convex: ", isConvex_);

//        // todo set scalar product with hessian

//        do
//        {
//            // Todo
//            // inefficient to create model every time this way
//            auto cubicModel = makeMyCubicModel(f_, x , dx, normDx, omega_);

//            lambda = Scalar::findLogGlobalMinimizer(cubicModel, 1e-10, 2.0, 0.001);

//            auto lambdaDx = get(lambda) * dx;

//            //            for(double i = 1.0; i > 1e-6; i*=0.95)
//            //            {
//            //                auto dxDummy = i*dx;
//            //                auto normDxDummy = norm(dxDummy);
//            //                std::cout << "Directional Minimizer: " << "i: " << i << " " << f_.d1(x)(dxDummy) + f_.d2(x,dxDummy)(dxDummy) + (omega_/2.0)*normDxDummy*normDxDummy*normDxDummy << std::endl;
//            //            }


//            // choose the best update
//            //nonlinUpdate_(x,lambdaDx,step);

//            LOG(log_tag, "lambda: ", lambda, " omega: ", omega_, "|dx|: " , norm(lambdaDx), " cubicModel(lambda): ", cubicModel(Real(get(lambda)))," cubicModel(0): ", cubicModel(0.0));
//            LOG(log_tag, "f_(x+dx): ", f_(x+lambdaDx), "f_(x): ", f_(x));
//            LOG(log_tag, "f_(x+dx) - f_(x): ", f_(x+lambdaDx) -f_(x));

//            // todo order weightUpdate and acceptanceTest


//            stepMonitor_ = acceptanceTest(x,lambdaDx);
//            omega_ = weightUpdate(x, lambdaDx);

//            if(stepMonitor_  == StepMonitor::Accepted)
//            {
//                dx = lambdaDx;
//                normDx = norm(dx);
//            }
//        }
//        while (stepMonitor_ == StepMonitor::Rejected );

//        x += dx;

//        LOG(log_tag,  " omega: ", omega_, "|dx|: " , normDx, "lambda: ", lambda, "|x|: ", normX );
//        LOG_SEPARATOR(log_tag);

//        if(terminate(normX,normDx,lambda))
//        {
//            LOG_INFO(log_tag,  "Converged ");
//            return x;
//        }
//    }

//    return x;
//}

//Vector ElasticityACRSolver::operator()()
//{
//    return operator()(zero(domain_));
//}

//void ElasticityACRSolver::setNonlinUpdate(
//        std::function<void(const Vector & x, Vector& dx, unsigned int step)> nonlinUpdate)
//{
//    nonlinUpdate_ = std::move(nonlinUpdate);
//}



//void ElasticityACRSolver::setOutput(std::function<void(const Vector& x, const Vector & dx)> outputUpdate)
//{
//    output_ = std::move(outputUpdate);
//}



//Vector ElasticityACRSolver::computeStep(const Vector & x)
//{
//    LinearSolver solver = { };
//    solver = f_.hessian(x) ^ -1;

//    if (is<CG::LinearSolver>(solver))
//    {
//        const auto& cgSolver = cast_ref<CG::LinearSolver>(solver);

//        // makeTRCGSolver( Operator A, CallableOperator P, CallableOperator R, Real theta_sugg
//        //auto cg =  makeTRCGSolver( f_.hessian(x), cgSolver.P(), f_.hessian(zero(domain_)), 0.0);

//        auto cg = makeTCGSolver(f_.hessian(x), cgSolver.P(), eps(), verbose());
//        cg.setRelativeAccuracy(getRelativeAccuracy());
//        cg.setVerbosity(2);

//        auto result = cg(-f_.d1(x));

//        isConvex_ = cg.isPositiveDefinite();

//        return result;
//    }

//    auto tcg = makeTCGSolver(f_.hessian(x), solver);
//    auto result = tcg(-f_.d1(x));

//    // Check if direct solver can detect nonconvexities
//    isConvex_ = tcg.isPositiveDefinite();

//    return result;
//}

//Real ElasticityACRSolver::weightUpdate(const Vector &x ,const Vector & dx)
//{
//    auto f_x_next =  f_(x+dx);
//    auto normDx = norm(dx);

//    if(std::isnan(get(f_x_next)) || std::isinf(get(f_x_next)))
//        return omega_*2.0;

//    else if(normDx < stabilityAcc_)
//    {
//        return (2.0/(normDx*normDx*normDx))*fabs(get(f_.d1(x+dx)(dx)-f_.d1(x)(dx) -f_.d2(x,dx)(dx)));
//    }

//    else
//    {
//        //std::cout << "Directionla Minimizer: " << f_.d1(x)(dx) + f_.d2(x,dx)(dx) + (omega_/2.0)*normDx*normDx*normDx;
//        //std::cout << "TestUpdate: " << "normDx: " << normDx << " f_x_next: " << f_x_next << " f_(x): " << f_(x) << " f_.d1(x)(dx): " << f_.d1(x)(dx) << " f_.d2(x,dx)(dx): " << f_.d2(x,dx)(dx) << std::endl;
//        //std::cout << "(6.0/(normDx*normDx*normDx): " << (6.0/(normDx*normDx*normDx)) << " -f_(x): " << -f_(x) << " -f_.d1(x)(dx): " <<  -f_.d1(x)(dx) <<  " -0.5*f_.d2(x,dx)(dx):"  <<  -0.5*f_.d2(x,dx)(dx) << std::endl;
//        //std::cout << "f_x_next -f_(x) -f_.d1(x)(dx) -0.5*f_.d2(x,dx)(dx): " << f_x_next -f_(x) -f_.d1(x)(dx) -0.5*f_.d2(x,dx)(dx) << std::endl;
//        //std::cout << "(6.0/(normDx*normDx*normDx))*fabs(get(f_x_next -f_(x) -f_.d1(x)(dx) -0.5*f_.d2(x,dx)(dx))): " << (6.0/(normDx*normDx*normDx))*fabs(get(f_x_next -f_(x) -f_.d1(x)(dx) -0.5*f_.d2(x,dx)(dx))) << std::endl;

//        return std::max(rhoMin_*omega_,std::min(omega_*rhoMax_,(6.0/(normDx*normDx*normDx))*fabs(get(f_x_next -f_(x) -f_.d1(x)(dx) -0.5*f_.d2(x,dx)(dx)))));
//    }
//}

//bool ElasticityACRSolver::terminate(Real normX, Real normDx, DampingFactor lambda) const
//{
//    if(lambda < 0.95)
//        return false;

//    else
//    {
//        return (isConvex_ && (normDx < getRelativeAccuracy() * normX
//                              || (normX < eps() && normDx < eps())));
//    }
//}


//ElasticityACRSolver::StepMonitor ElasticityACRSolver::acceptanceTest(const Vector& x, const Vector& dx)
//{

//    auto f_next = f_(x+dx);
//    // std::cout << (std::isnan(get(f_next)) || std::isinf(get(f_next))) << std::endl;

//    if(std::isnan(get(f_next)) || std::isinf(get(f_next)))
//        return StepMonitor::Rejected;

//    auto normDx = norm(dx);

//    if( normDx < stabilityAcc_)
//    {
//        if(f_.d1(x+dx)(dx)< (omega_/6.0)*normDx*normDx*normDx)
//        {
//            return StepMonitor::Accepted;
//        }
//        else
//            return StepMonitor::Rejected;
//    }

//    else{

//        std::cout << "Acceptance test: " << std::endl;
//        std::cout << "f_next: " << f_next << " f_(x): " << f_(x) << " f_.d1(x)(dx): " << f_.d1(x)(dx) << " omega_: " << omega_ << " normDx: " << normDx << std::endl;


//        if(f_next < f_(x) +0.5* f_.d1(x)(dx)- (omega_/36.0)*normDx*normDx*normDx)
//        {
//            return StepMonitor::Accepted;
//        }
//        else
//            return StepMonitor::Rejected;
//    }

//    return StepMonitor::Rejected;
//}


//ElasticityACRSolver::StepMonitor ElasticityACRSolver::myAcceptanceTest(const Vector& x, const Vector& dx)
//{
//    auto f_next = f_(x+dx);
//    // std::cout << (std::isnan(get(f_next)) || std::isinf(get(f_next))) << std::endl;

//    if(std::isnan(get(f_next)) || std::isinf(get(f_next)))
//        return StepMonitor::Rejected;

//    auto normDx = norm(dx);

//    auto cubicModel = makeMyCubicModel(f_, x , dx, normDx, omega_);

//    // if dx is already damped, then lambda = 1.0.
//    const auto diffModel = cubicModel(1.0) - cubicModel(0.0);
//    const auto diffFunction = (f_(x + dx) - f_(x));


//    std::cout << "diffFunction: " << diffFunction << std::endl;
//    std::cout << "diffModel: " <<  diffModel << std::endl << std::endl;

//    if( (diffFunction / diffModel) >= 0.25)
//        return StepMonitor::Accepted;

//    return StepMonitor::Rejected;

//}




//void ElasticityACRSolver::regularityTest(DampingFactor lambda) const
//{
//    if (!regularityTestPassed(lambda))
//        throw Exception::RegularityTestFailed( "ACR::regularityTest (lambda,...)", get(get(lambda)));
//}

//void ElasticityACRSolver::setLogName(const std::string & logName)
//{
//    logName_ = logName;
//}


//Vector ElasticityACRSolver::directionalMinimizer(const Vector & x_, const Vector & dx_)
//{
//    auto x = x_;
//    auto dx = dx_;
//    auto normDx = norm(dx);

//    omega_ = 1e-2;


//    // todo make temporäre abstiegs bedingung

//    stepMonitor_ = StepMonitor::Rejected;

//    //    std::cout << "Test wether direction of descent: " << std::endl;

//    //    for(double i = 1.0; i > 1e-8; i*=0.1)
//    //        std::cout << f_(x+i*dx) -f_(x) << std::endl;

//    for(auto i = 0u; i < getMaxSteps() && (stepMonitor_ == StepMonitor::Rejected); i++)
//    {

//        //        for(double i = 1.0; i > 1e-6; i*=0.99)
//        //        {
//        //            auto dxDummy = i*dx;
//        //            auto normDxDummy = norm(dxDummy);
//        //            //std::cout << "Directional Minimizer: " << "i: " << i << " " << f_.d1(x)(dxDummy) + f_.d2(x,dxDummy)(dxDummy) + (omega_/2.0)*normDxDummy*normDxDummy*normDxDummy << std::endl;
//        //            //std::cout << "Model: " << "i: " << i << " " << f_.d1(x)(dxDummy) + 0.5*f_.d2(x,dxDummy)(dxDummy) + (omega_/6.0)*normDxDummy*normDxDummy*normDxDummy << std::endl;
//        //            std::cout << "f_(x+i*dx) -f_(x): " << f_(x+i*dx) -f_(x) << std::endl;
//        //        }

//        std::cout << "I: " << i << " MaxSteps: " << getMaxSteps() << std::endl << std::endl;
//        auto cubicModel = makeMyCubicModel(f_, x , dx, normDx, omega_);
//        // auto lambda = Scalar::findLogGlobalMinimizer(cubicModel, 1e-10, 2.0, 0.05);

//        auto lambda = Scalar::findLogGlobalMinimizer(cubicModel, 1e-10, 2.0, 0.001);
//        // todo why does omega not change
//        std::cout << "lambda: " << lambda << std::endl << std::endl;

//        //        for(double i = 1.0; i > 1e-6; i *= 0.99 )
//        //        {
//        //            auto dxDummy = i*dx;
//        //            auto normDxDummy = norm(dxDummy);
//        //            std::cout << " i: " << i << std::endl;
//        //            std::cout << "Derivative: " << f_.d1(x)(dxDummy) + f_.d2(x,dxDummy)(dxDummy) + (omega_/2.0)*normDxDummy*normDxDummy*normDxDummy << std::endl;
//        //            std::cout << "NewModel: " << f_.d1(x)(dxDummy) + 0.5*f_.d2(x,dxDummy)(dxDummy) + (omega_/6.0)*normDxDummy*normDxDummy*normDxDummy << std::endl;
//        //            std::cout << "Model: " << cubicModel(i) << std::endl;
//        //            //std::cout << "f_(x+i*dx) -f_(x): " << f_(x+i*dx) -f_(x) << std::endl;
//        //        }

//        std::cout << std::endl;

//        //        for(double i = 1.0; i > 1e-8; i*=0.8)
//        //            std::cout << "cubicModel(i): "  << cubicModel(i)  << std::endl;

//        auto lambdaDx = get(lambda) * dx;

//        std::cout << "model(lambda): " << cubicModel(lambda) <<  " model(0): " << cubicModel(0) << std::endl;
//        std::cout << "f_(x+ lambda dx): " << f_(x+lambdaDx) <<  "  f_(x): " << f_(x) << std::endl;

//        // todo before or after omega update ???????????????!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



//        //   stepMonitor_ = myAcceptanceTest(x,lambdaDx);
//        stepMonitor_ = acceptanceTest(x,lambdaDx);

//        omega_ = weightUpdate(x, lambdaDx);

//        std::cout << "Omega Updated: " << omega_ << std::endl << std::endl;

//        if(stepMonitor_ == StepMonitor::Accepted)
//        {
//            std::cout << "Accepted: " << std::endl;
//            dx = lambdaDx;
//        }
//    }

//    if(stepMonitor_ == StepMonitor::Rejected )
//    {
//        throw Exception::RegularityTestFailed( "ACR", 1.0);
//    }
//    return dx;
//}



//Vector ElasticityACRSolver::convexProjection(const Vector & x_)
//{

//    auto x = x_;
//    auto lambda = DampingFactor{1.0};

//    LOG_INFO(log_tag, "Starting ACR-Solver for Convex Projection.");

//    for (unsigned step = 1; step <= getMaxSteps(); ++step)
//    {
//        LOG_SEPARATOR(log_tag);
//        LOG(log_tag, "Iteration", step);

//        auto dx = computeStep(x);
//        auto normDx = norm(dx);
//        auto normX = norm(x);

//        LOG(log_tag, "Is convex: ", isConvex_);

//        if(isConvex_)
//            return x;

//        do
//        {
//            // Todo
//            // inefficient to create model every time this way
//            auto cubicModel = makeMyCubicModel(f_, x , dx, normDx, omega_);

//            lambda = Scalar::findLogGlobalMinimizer(cubicModel, 1e-10, 2.0, 0.001);

//            auto lambdaDx = get(lambda) * dx;

//            LOG(log_tag, "lambda: ", lambda, " omega: ", omega_, "|dx|: " , norm(lambdaDx), " cubicModel(lambda): ", cubicModel(Real(get(lambda)))," cubicModel(0): ", cubicModel(0.0));
//            LOG(log_tag, "f_(x+dx): ", f_(x+lambdaDx), "f_(x): ", f_(x));
//            LOG(log_tag, "f_(x+dx) - f_(x): ", f_(x+lambdaDx) -f_(x));

//            // todo order correct
//            // has to be this way
//            stepMonitor_ = acceptanceTest(x,lambdaDx);
//            omega_ = weightUpdate(x, lambdaDx);

//            if(stepMonitor_ == StepMonitor::Accepted)
//            {
//                dx = lambdaDx;
//                normDx = norm(dx);
//            }
//        }

//        while (stepMonitor_ == StepMonitor::Rejected );

//        x += dx;

//        LOG(log_tag,  " omega: ", omega_, "|dx|: " , normDx, "lambda: ", lambda, "|x|: ", normX );
//        LOG_SEPARATOR(log_tag);

//    }

//    return x;
//}

//bool ElasticityACRSolver::isConvex(const Vector & x_)
//{
//    std::cout << "Test if problem is convex: " << std::endl;
//    computeStep(x_);
//    return isConvex_;
//}



}
}
