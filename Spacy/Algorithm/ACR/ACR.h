#pragma once

#include <chrono>

#include <Spacy/LinearSolver.h>
#include <Spacy/C2Functional.h>
#include <Spacy/Vector.h>
#include <Spacy/VectorSpace.h>
#include <Spacy/Util/Mixins.h>
#include <Spacy/Spaces/ScalarSpace/Real.h>

#include <Spacy/Operator.h>

namespace Spacy
{
/// @cond
namespace CompositeStep { class  CubicModel; }
/// @endcond


namespace ACR
{
class MyCubicModel
{
public:

    MyCubicModel(Real a0, Real a1, Real a2, Real a3);

    Real operator ()(Real x);

    Real operator()() const;

    void update(Real x);

private:

    const Real a0_ = 0.0;
    const Real a1_ = 0.0;
    const Real a2_ = 0.0;
    const Real a3_ = 0.0;

    Real x_ = 1.0;
};

MyCubicModel makeMyCubicModel(C2Functional f, const Vector &x, const Vector & dx, Real normDx, Real omega);
MyCubicModel makeMyCubicModelFunction(C2Functional  f,const Vector & x,const Vector & dx, Real normDx, Real omega);
/**
         * @brief Adaptive cubic regularization approach for solving unconstrained minimization problems
         * Compute direction of descent with one Newton step
         * Use a cubic error model to obtain acceptable corrections
         * Solve problems of the form
         * \f[\min f(x).\f]
         */
class ACRSolver :
        public Mixin::Eps,
        public Mixin::MaxSteps,
        public Mixin::RelativeAccuracy,
        public Mixin::AbsoluteAccuracy,
        public Mixin::Verbosity
{
    enum class StepMonitor { Rejected , Accepted };

public:
    /**
             * @brief Constructor.
             * @param f functional to minimize
             */
    ACRSolver(C2Functional f, VectorSpace & domain, double omega=1, double eta1 = 0.25, double eta2 = 0.5, double epsilon = 1e-8, double relativeAccuracy = 1e-4, double omegaMax = 1e25,  double lambdaMin = 1e-10);

    /// Compute solution starting at \f$x_0=0\f$.
    Vector operator()();

    /**
             * @brief Compute solution.
             * @param x0 initial iterate
             */
    Vector operator()(const Vector& x0);


    /**
             * @brief Compute solution for a parameter-dependent system.
             * @param x0 initial iterate
             */
    std::tuple<bool, Vector, Real, Real, unsigned int> solveParam(
            const Vector & x0, Real lambda, Real thetaGlobal);

    /**
             * @brief Set nonlinear update function
             */
    void setNonlinUpdate(std::function<void(const Vector & x,  Vector& dx, unsigned int step)> nonlinUpdate);


    /**
            * @brief Set nonlinear update function
                       */
    void setOutput(std::function<void(const Vector& x, const Vector & dx)> outputUpdate);

    /**
             * @brief Get C2Functional
                        */
    C2Functional & getFunctional();


    // Just for Testing:

    std::function<void(const Vector &x, const Vector &dx )> testNonlinearity = [](const Vector &x, const Vector & dx){};

    void setPositivePrecond(CallableOperator p);

    std::function<void(const Vector& x)> plot_ = []( const Vector& x){return;};

    /* @brief Flag variable: Only solve the optimization problem until the iterates reach a region of convexity
               */

    std::function<C2Functional(double lambda)> updateFunctional = []( double lambda ){ throw;  C2Functional f_temp;  return f_temp;};
    double energy_ = 0.0;
    bool convexProjection_ = false;
    bool updateFunctional_ = false;
   mutable double energyRegularization_ = 0.0;

private:
    /**
             * @brief Test if dx is an acceptable correction.
             * @param x current iterate
             * @param dx correction
             * @param lambda damping factor
             * @param cubicModel cubic model
             */
    StepMonitor acceptanceTest(const Vector& x, const Vector& dx, const Real& lambda,  MyCubicModel& cubicModel);

    /**
             * @brief Update the weight Parameter of the cubic model.
             * @param omega weight Parameter of the cubic model
             */
    Real weightChange(::Spacy::Real omega) const;

    /**
             * @brief Compute correction dx.
             * @param x current iterate
             */
    std::tuple<Vector,Real> computeStep(const Vector& x, Real dxNorm) const;


    /**
             * @brief Check if the contraction of the update dx is sufficiently small.
             */
    std::tuple<bool,Real,Real,Real,Real> checkContraction(int step, Real normDx,Real normDx1_, Real normDx2_, Real thetaZero_, Real thetaGlobal) const;

    /**
             * @brief Check if the convergence criterion is satisfied
             */
    bool convergenceTest( Real dxQNorm, const Vector & x) const;


    mutable bool isConvex_ = true;

    C2Functional f_temp;
    mutable C2Functional f_;
    VectorSpace& domain_;
    StepMonitor stepMonitor = StepMonitor::Accepted;

    // epsilon_ termination criterion
    const double eta1_, eta2_, epsilon_, omegaMax_,  lambdaMin_;
    Real rho_ = 1.0;
    Real omega_ = 5.0;
    std::function<void(Vector & x,  Vector& dx, unsigned int step)>  nonlinUpdate_ = [](Vector & x,  Vector& dx, unsigned int step){return;};
    std::function<void(const Vector& x, const Vector & dx)> output_ = []( const Vector& x, const Vector & dx){return;};

    ::Spacy::CallableOperator positivePrecond_  = [] ( const Vector& x ) { return x; };
};





class ElasticityACRSolver :
        public Mixin::Eps,
        public Mixin::MaxSteps,
        public Mixin::RelativeAccuracy,
        public Mixin::Verbosity,
        public Mixin::DampingAccuracy,
        public Mixin::RegularityTest
{
    enum class StepMonitor { Rejected , Accepted };

public:


    /**
            * @brief Default constructor is deleted
            */
    ElasticityACRSolver() = delete;

    /**
            * @brief Constructor.
            * @param f functional to minimize
            */
    ElasticityACRSolver(C2Functional f, VectorSpace &domain);

    Vector operator()(const Vector & x0);
    Vector operator()();


    /**
             * @brief Set nonlinear update function
             */
    void setNonlinUpdate(std::function<void(const Vector & x,  Vector& dx, unsigned int step)> nonlinUpdate);


    /**
            * @brief Set output function
                       */
    void setOutput(std::function<void(const Vector& x, const Vector & dx)> outputUpdate);

    /**
     * @brief Set name for different log files
     * */
    void setLogName(const std::string & logName);


    /**
             * @brief Compute directional minimizer along dx
                        */
    Vector directionalMinimizer(const Vector & x_, const Vector & dx_);

    /**
             * @brief Checks if the problem is convex
                        */
    bool isConvex(const Vector & x_);

    /**
    * @brief Set positive preconditioner in case of nonconvex
               */
    void setPositivePrecond(CallableOperator p);


    std::function<void(const Vector& x)> plot_ = []( const Vector& x){return;};


    /**
 * @brief Flag variable: Only solve the optimization problem until the iterates reach a region of convexity
            */
    bool convexProjection_ = false;
private:

    /**
* @brief Compute Newton direction
     **/
    Vector computeStep(const Vector & x, Real normDx);

    /**
* @brief Update the affine conjugate Lipschitz constant \omega
     **/
    Real weightUpdate(const Vector & x, const Vector &dx);


    /**
* @brief Test convergence of the algorithm
     **/
    bool terminate(Real normX, Real normDx, DampingFactor lambda) const;

    /**
* @brief Test is step gets accepted
     **/
    StepMonitor acceptanceTest(const Vector& x, const Vector& dx, Real f_next);


    /**
* @brief Test is step gets accepted
     **/
    StepMonitor myAcceptanceTest(const Vector& x, const Vector& dx);

    /**
* @brief Test regularity of the damping factor \lambda
     **/
    void regularityTest(DampingFactor lambda) const;

    VectorSpace & domain_;
    C2Functional f_ = {};

    mutable Real omega_ = 1.0;
    bool isConvex_ = false;
    StepMonitor stepMonitor_ = StepMonitor::Accepted;

    const Real stabilityAcc_ = {1e-3};

    // Maximum and minimum increase and decrease factors for updating omega
    const Real rhoMin_ = {0.3};
    const Real rhoMax_ = {3.0};

    std::function<void(Vector & x,  Vector& dx, unsigned int step)>  nonlinUpdate_ = [](Vector & x,  Vector& dx, unsigned int step){return;};
    std::function<void(const Vector& x, const Vector & dx)> output_ = []( const Vector& x,  const Vector& dx){return;};

    std::string logName_ = { };
    ::Spacy::CallableOperator positivePrecond_  = [] ( const Vector& x ) { return x; };

};

}
}

