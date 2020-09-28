#include "solver.hh"

#include <Spacy/Spaces/ProductSpace.h>
#include <Spacy/ZeroVectorCreator.h>
#include <Spacy/Util/Mixins/Get.h>
#include <Spacy/Util/Logger.h>
#include <Spacy/Util/Log.h>
#include <iostream>

#define LAPACK_COL_MAJOR 102

namespace Spacy
{
namespace
{
    Logger< unsigned > logIterationsState( "stateIt.log" );
    Logger< unsigned > logIterationsAdjoint( "adjointIt.log" );
    Logger< unsigned > logIterationsIcg( "icgIt.log" );
}
namespace ICG
{

DEFINE_LOG_TAG(static const char* log_tag = "NewICG");


Solver::Solver(const Spacy::ProductSpace::Operator_V2& H, const Spacy::CallableOperator& BPX, const Spacy::CallableOperator& BPXT,
               const Spacy::CallableOperator& controlSolver, const ::Spacy::CallableOperator& diagMuu)
    : OperatorBase( H.domain(), H.range() ), H_(H), BPX_(BPX), BPXT_(BPXT), controlSolver_(controlSolver),  diagMuu_(diagMuu)
{
}

Solver::Solver(const ::Spacy::ProductSpace::Operator_V2& H, const ::Spacy::CallableOperator& BPX, const ::Spacy::CallableOperator& BPXT,
               const ::Spacy::CallableOperator& controlSolver,
               const ::Spacy::CallableOperator& Ry, const ::Spacy::CallableOperator& Ru, const CallableOperator &diagMuu)
    : OperatorBase( H.domain(), H.range() ), H_(H), BPX_(BPX), BPXT_(BPXT), controlSolver_(controlSolver),
      Ry_(Ry), Ru_(Ru),  diagMuu_(diagMuu)
{
}

Vector Solver::operator()(const Vector &b) const
{
    std::cout << std::endl;
    ::Spacy::Vector x = zero(domain());
    ::Spacy::Vector result = zero(domain());
    theta_ = 0.0;
    convex_flag_ = false;
    bool logConvex = true;
    result_ = Result::Failed;
    definiteness_ = DefiniteNess::PositiveDefinite;

    while(!convex_flag_)
    {
        result = outer_pcgLoop(x,b);
        convex_flag_ = isPositiveDefinite();

        LOG(log_tag, "ICG is Convex: ", convex_flag_);
        LOG(log_tag, "Energy in ICG is Convex: ", signalConvex_(false,false));

        if(!signalConvex_(false,false))
        {
            return result;
        }

        if(theta_ == 0.0)
            theta_ = 1e-8;
        else
            theta_ *= 100.0;

        if(!convex_flag_)
            logConvex = false;
    }

    if(!logConvex)
        definiteness_ = DefiniteNess::Indefinite;

    return result;
}

const Spacy::ProductSpace::Operator_V2& Solver::H() const
{
    return H_;
}

const CallableOperator& Solver::BPX() const
{
    return BPX_;
}

const CallableOperator& Solver::BPXT() const
{
    return BPXT_;
}

void Solver::setNorm(CallableOperator Myy , CallableOperator Muu)
{
    Myy_ = std::move(Myy);
    Muu_ = std::move(Muu);
}

const ::Spacy::CallableOperator& Solver::controlSolver() const
{
    return controlSolver_;
}

bool Solver::isPositiveDefinite() const
{
    return definiteness_ == DefiniteNess::PositiveDefinite;
}

Vector Solver::outer_pcgLoop(Vector x,const Vector& b) const
{
    LOG_INFO(log_tag, "Starting ICG-Solver.");
    LOG_SEPARATOR(log_tag);

    LOG(log_tag, "Regularization Parameter: ", theta_);
    LOG(log_tag, "Relative Accuracy: ", getRelativeAccuracy());
    LOG(log_tag, "Eps: ", eps());

    auto A = H_(2,0);
    auto AT = H_(0,2);
    auto B = H_(2,1); //=^ -B
    auto BT = H_(1,2); //=^ -BT

    convex_flag_ = true;
    energyIsConvex_ = true;
    definiteness_ = DefiniteNess::PositiveDefinite;
    icg_iterations_ = 0;

    auto stateAcc = 1e-12;//getRelativeAccuracy();
    auto adjointAcc = 1e-12;//getRelativeAccuracy();
    auto inexactCGAcc = getRelativeAccuracy();
    auto absAccIncrease = 1.0;

    ::Spacy::CG::NoRegularization NR;
    auto adjointSolver = Spacy::CG::Solver(AT, BPXT_, NR, true);
    adjointSolver.set_eps(eps());
    adjointSolver.setAbsoluteAccuracy(eps());
    adjointSolver.setRelativeAccuracy(adjointAcc);
    adjointSolver.setMaxSteps(getMaxSteps());
    adjointSolver.setTerminationCriterion( CG::Termination::AdaptiveRelativeEnergyError() );
    //adjointSolver.setTerminationCriterion( CG::Termination::StrakosTichyEnergyError() );
    adjointSolver.transfer = transfer;
    adjointSolver.setVerbosityLevel(0);

    auto stateSolver = Spacy::CG::Solver(A, BPX_, NR, true);
    stateSolver.set_eps(eps());
    stateSolver.setAbsoluteAccuracy(eps());
    stateSolver.setRelativeAccuracy(stateAcc);
    stateSolver.setMaxSteps(getMaxSteps());
    stateSolver.setTerminationCriterion( CG::Termination::AdaptiveRelativeEnergyError() );
    //stateSolver.setTerminationCriterion( CG::Termination::StrakosTichyEnergyError() );
    stateSolver.transfer = transfer;
    stateSolver.setVerbosityLevel(0);

    ::Spacy::PPCG::TriangularStateConstraintPreconditioner P(stateSolver,controlSolver_,adjointSolver,B,BT,domain(),range());
    P.output = output2;
    P.signalConvex_ = signalConvex_;

    ::Spacy::InexactSolver inexactCg(P,H_,CG::Termination::AdaptiveRelativeEnergyError(),domain(),range(), diagMuu_, theta_);
    //::Spacy::InexactSolver inexactCg(P,BT,Hyy,Hyu,Huu,Huy,CG::Termination::AdaptiveRelativeEnergyError(),domain(),range(), Muu_, diagMuu_, theta_);
    //::Spacy::InexactSolver inexactCg(P,BT,Hyy,Hyu,Huu,Huy,CG::Termination::StrakosTichyEnergyError(),domain(),range(), Muu_, diagMuu_, theta_);
    inexactCg.set_eps(eps()*absAccIncrease);
    inexactCg.setAbsoluteAccuracy(eps()*absAccIncrease);
    inexactCg.setRelativeAccuracy(inexactCGAcc);
    inexactCg.setMaxSteps(getMaxSteps());
    inexactCg.output = output2;
    //inexactCg.setVerbosity(2);
    inexactCg.setRegularizationOperators(Ry_,Ru_);


    LOG(log_tag, "AdjointAcc: ", adjointAcc);
    LOG(log_tag, "StateAcc: ", stateAcc);

    x = inexactCg(b);

    logIterationsIcg(inexactCg.getIterations());
    icg_iterations_ += inexactCg.getIterations();

    LOG(log_tag, "ICG Number of Steps: ", inexactCg.getIterations());
    LOG(log_tag, "ICG is Positive Definite: ", inexactCg.isPositiveDefinite());

    convex_flag_ = inexactCg.isPositiveDefinite();
    if(!convex_flag_)
    {
        definiteness_ = DefiniteNess::Indefinite;
        result_ = Result::TruncatedAtNonConvexity;
        return x;
    }

    energyIsConvex_ = signalConvex_(false,false);
    if(!energyIsConvex_){
        result_ = Result::NotConverged;
        return x;
    }

    result_ = Result::Converged;
    std::cout << std::endl;
    return x;
}

/*
bool Solver::convergenceTest(bool isConvex, bool energyIsConvex, CallableOperator Hyy, CallableOperator Huu , const Vector &x_y, const Vector &x_u, const Vector &dx_y, const Vector &dx_u) const
{
    if(!isConvex || !energyIsConvex)
        return false;
    else
    {
        auto normX = sqrt(x_y(Myy_(x_y)) + x_u(Muu_(x_u)));
        auto normDx = sqrt(dx_y(Myy_(dx_y)) + dx_u(Muu_(dx_u)));

        LOG(log_tag, "normX: ", normX, "normDx: ", normDx);

        if ( (normDx < getRelativeAccuracy() * normX) || ( normX < eps() && normDx < eps() ) ||  normDx < eps()  )
        {
            // if ( verbose() )
            //   std::cout  << "PPCG Terminating (convergent)." << std::endl;
            return true;
        }
    }
    return false;
}*/
}
}
