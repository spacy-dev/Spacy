#define FUSION_MAX_VECTOR_SIZE 15

#include <chrono>
#include <iostream>
#include <functional>

#include <dune/grid/config.h>
#include <dune/grid/uggrid.hh>

#include <Spacy/Adapter/kaskadeParabolic.hh>
#include <Spacy/Algorithm/CompositeStep/affineCovariantSolver.hh>
#include <Spacy/Spaces/ScalarSpace/Real.hh>

#include <fem/forEach.hh>
#include <fem/gridmanager.hh>
#include <fem/lagrangespace.hh>
#include <io/vtk.hh>
#include <utilities/gridGeneration.hh>
#include <utilities/kaskopt.hh>

#define NCOMPONENTS 1
#include "Spacy/Adapter/KaskadeParabolic/Constraints/nonlinOcpHeader.hh"

int main( int argc, char* argv[] )
{
    using namespace Kaskade;
    static constexpr int stateId = 0;
    static constexpr int controlId = 1;
    static constexpr int adjointId = 2;

    constexpr int dim = 2;
    int silence = 0;
    std::unique_ptr< boost::property_tree::ptree > pt =
        getKaskadeOptions( argc, argv, silence, false );

    double desiredAccuracy = getParameter( pt, "desiredAccuracy", 1e-6 );
    double eps = getParameter( pt, "eps", 1e-12 );
    double alpha = getParameter( pt, "alpha", 1e-1 );
    int maxSteps = getParameter( pt, "maxSteps", 500 );
    int initialRefinements = getParameter( pt, "initialRefinements", 5 );
    int iterativeRefinements = getParameter( pt, "iterativeRefinements", 0 );
    int FEorder = getParameter( pt, "FEorder", 1 );
    int verbose = getParameter( pt, "verbose", 2 );
    double c = getParameter( pt, "cPara", 0. );
    double d = getParameter( pt, "dPara", 1. );
    double e = getParameter( pt, "ePara", 0. );
    double mu = getParameter( pt, "muPara", 0. );
    double desContr = getParameter( pt, "desiredContraction", 0.5 );
    double relDesContr = getParameter( pt, "relaxedContraction", desContr + 0.1 );
    double maxContr = getParameter( pt, "maxContraction", 0.75 );

    using std::cout;
    using std::endl;
    cout << "cPara: " << c << " dPara: " << d << " muPara: " << mu << " alpha:" << alpha << endl;
    cout << "dContr: " << desContr << " rdContr: " << relDesContr << " mContr: " << maxContr
         << endl;

    using Grid = Dune::UGGrid< dim >;
    using H1Space = FEFunctionSpace< ContinuousLagrangeMapper< double, Grid::LeafGridView > >;
    using Spaces = boost::fusion::vector< H1Space const*, H1Space const* >;
    using VY = VariableDescription< 0, 1, stateId >;
    using VU = VariableDescription< 1, 1, controlId >;
    using VP = VariableDescription< 0, 1, adjointId >;
    using VariableDescriptions = boost::fusion::vector< VY, VU, VP >;
    using Descriptions = VariableSetDescription< Spaces, VariableDescriptions >;

    std::cout << "alpha " << alpha << " dPara " << d << " cPara " << c << std::endl;

    using NormalFunctionalDefinition =
        NormalStepFunctional< stateId, controlId, adjointId, double, Descriptions >;
    std::function< NormalFunctionalDefinition( const typename Descriptions::VariableSet ) >
        normalFuncGenerator =
            [&c, &d, &e, &mu, &alpha]( const typename Descriptions::VariableSet y_ref ) {
                return NormalFunctionalDefinition( alpha, y_ref, c, d, e, mu );
            };

    using TangentialFunctionalDefinition =
        TangentialStepFunctional< stateId, controlId, adjointId, double, Descriptions >;
    std::function< TangentialFunctionalDefinition( const typename Descriptions::VariableSet ) >
        tangentialFuncGenerator =
            [&c, &d, &e, &mu, &alpha]( const typename Descriptions::VariableSet y_ref ) {
                return TangentialFunctionalDefinition( alpha, y_ref, c, d, e, mu );
            };

    using HeatFunctionalDefinition = NonlinearModelPDE<
        double, ::Spacy::KaskadeParabolic::Detail::ExtractDescription_t< Descriptions, 0 >,
        ::Spacy::KaskadeParabolic::Detail::ExtractDescription_t< Descriptions, 1 > >;
    std::function< HeatFunctionalDefinition(
        const typename ::Spacy::KaskadeParabolic::Detail::ExtractDescription_t< Descriptions,
                                                                                1 >::VariableSet ) >
        forwardFunctionalGenerator = [&c, &d, &e, &mu](
            const typename ::Spacy::KaskadeParabolic::Detail::ExtractDescription_t<
                Descriptions, 1 >::VariableSet control ) {
            return HeatFunctionalDefinition( c, d, e, mu, control );
        };

    ::Spacy::Real t_end = 3;
    ::Spacy::Real tau = 0.5;
    unsigned mpcsteps = 4;
    auto expparam = -1.;
    ::std::vector< unsigned > NVec = {3, 6, 11, 21, 41};

    /// SAVE OBJECTIVE FUNCTION VALUE
    std::ofstream obj;
    obj.open( "data/Objective_Function.txt", std::ofstream::out | std::ofstream::app );
    obj << "alpha: " << alpha << ", c: " << c << ", d: " << d << ", mu: " << mu
        << ", t_end: " << t_end << ", Mpc Steps: " << mpcsteps << ", tau: " << tau << std::endl;
    obj.close();

    for ( auto N : NVec )
    {

        // Solve with an uniform grid
        std::cout << "##########################" << std::endl;
        std::cout << "####    UNIFORM   " << N << " ####" << std::endl;
        std::cout << "##########################" << std::endl;
        ::Spacy::KaskadeParabolic::ModelPredictiveController<
            NormalFunctionalDefinition, TangentialFunctionalDefinition, HeatFunctionalDefinition >
            mpc_uni( normalFuncGenerator, tangentialFuncGenerator, forwardFunctionalGenerator,
                     t_end, N, mpcsteps, tau, 50 );
        mpc_uni.MPCloop();

        return 0;

        // Solve with an exponential grid
        std::cout << "##########################" << std::endl;
        std::cout << "####  Exponential " << N << " ####" << std::endl;
        std::cout << "##########################" << std::endl;
        ::Spacy::KaskadeParabolic::ModelPredictiveController<
            NormalFunctionalDefinition, TangentialFunctionalDefinition, HeatFunctionalDefinition >
            mpc_exp( normalFuncGenerator, tangentialFuncGenerator, forwardFunctionalGenerator,
                     t_end, N, mpcsteps, tau, 50, "exponential", expparam );
        //    mpc_exp(normalFuncGenerator,tangentialFuncGenerator,forwardFunctionalGenerator,t_end,N,mpcsteps,tau,"exponential",-0.5);
        mpc_exp.MPCloop();
    }

    /// Solve with extremely fine grid for reference
    std::cout << "##########################" << std::endl;
    std::cout << "####  Exact sol. " << 401 << " ####" << std::endl;
    std::cout << "##########################" << std::endl;
    {
        ::Spacy::KaskadeParabolic::ModelPredictiveController<
            NormalFunctionalDefinition, TangentialFunctionalDefinition, HeatFunctionalDefinition >
            mpc_uni( normalFuncGenerator, tangentialFuncGenerator, forwardFunctionalGenerator,
                     t_end, 401, mpcsteps, tau, 100 /*,"exponential",-0.5*/ );
        mpc_uni.MPCloop();
    }

    std::cout << "returning from main" << std::endl;
    return 0;
}
