#define FUSION_MAX_VECTOR_SIZE 15

#include <chrono>
#include <iostream>
#include <functional>

#include <dune/grid/config.h>
#include <dune/grid/uggrid.hh>

#include <Spacy/Adapter/kaskadeParabolic.hh>
#include <Spacy/Algorithm/CompositeStep/affineCovariantSolver.hh>
#include <Spacy/Spaces/ScalarSpace/Real.h>

#include <fem/forEach.hh>
#include <fem/gridmanager.hh>
#include <fem/lagrangespace.hh>
#include <io/vtk.hh>
#include <utilities/gridGeneration.hh>
#include <utilities/kaskopt.hh>

#define NCOMPONENTS 1
#include "optimal_control_withoutfung.hh"

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
    double alpha = getParameter( pt, "alpha", 0.5 );
    int maxSteps = getParameter( pt, "maxSteps", 500 );
    int initialRefinements = getParameter( pt, "initialRefinements", 4 );
    int iterativeRefinements = getParameter( pt, "iterativeRefinements", 0 );
    int FEorder = getParameter( pt, "FEorder", 1 );
    int verbose = getParameter( pt, "verbose", 2 );
    double c = getParameter( pt, "cPara", 1 );
    double d = getParameter( pt, "dPara", 1 );
    double e = getParameter( pt, "ePara", 0.0 );
    double desContr = getParameter( pt, "desiredContraction", 0.5 );
    double relDesContr = getParameter( pt, "relaxedContraction", desContr + 0.1 );
    double maxContr = getParameter( pt, "maxContraction", 0.75 );
    double expgridparam = getParameter( pt, "expGrid", -0.345 );

    using std::cout;
    using std::endl;
    cout << "cPara: " << c << " dPara: " << d << " ePara: " << e << " alpha:" << alpha << endl;
    cout << "dContr: " << desContr << " rdContr: " << relDesContr << " mContr: " << maxContr
         << endl;

    using Grid = Dune::UGGrid< dim >;
    using H1Space = FEFunctionSpace< ContinuousLagrangeMapper< double, Grid::LeafGridView > >;
    using Spaces = boost::fusion::vector< H1Space const* >;
    using VY = VariableDescription< 0, 1, stateId >;
    using VU = VariableDescription< 0, 1, controlId >;
    using VP = VariableDescription< 0, 1, adjointId >;
    using VariableDescriptions = boost::fusion::vector< VY, VU, VP >;
    using Descriptions = VariableSetDescription< Spaces, VariableDescriptions >;

    using NormalFunctionalDefinition =
        NormalStepFunctional< stateId, controlId, adjointId, double, Descriptions >;
    std::function< NormalFunctionalDefinition( const typename Descriptions::VariableSet ) >
        normalFuncGenerator =
            [&c, &d, &e, &alpha]( const typename Descriptions::VariableSet y_ref ) {
                return NormalFunctionalDefinition( alpha, y_ref, c, d, e );
            };

    using TangentialFunctionalDefinition =
        TangentialStepFunctional< stateId, controlId, adjointId, double, Descriptions >;

    std::function< TangentialFunctionalDefinition( const typename Descriptions::VariableSet ) >
        tangentialFuncGenerator =
            [&c, &d, &e, &alpha]( const typename Descriptions::VariableSet y_ref ) {
                return TangentialFunctionalDefinition( alpha, y_ref, c, d, e );
            };

    //####################### Exact SOLUTION #######################

    ::Spacy::KaskadeParabolic::GridManager< Spaces > gm1( 101, ::Spacy::Real{30.},
                                                          initialRefinements, FEorder );
    auto domain1 = Spacy::KaskadeParabolic::makeHilbertSpace( gm1, {0u, 1u}, {2u} );

    std::cout << "Creating Functionals" << std::endl;
    auto n_func1 = Spacy::KaskadeParabolic::makeC2Functional( normalFuncGenerator, gm1, domain1 );
    auto t_func1 =
        Spacy::KaskadeParabolic::makeC2Functional( tangentialFuncGenerator, gm1, domain1 );

    cout << "set up solver" << endl;
    // algorithm and parameters
    auto cs1 = Spacy::CompositeStep::AffineCovariantSolver( n_func1, t_func1, domain1 );
    cs1.setRelativeAccuracy( desiredAccuracy );
    cs1.set_eps( eps );
    cs1.setVerbosityLevel( 2 );
    cs1.setMaxSteps( maxSteps );
    cs1.setIterativeRefinements( iterativeRefinements );
    cs1.setDesiredContraction( desContr );
    cs1.setRelaxedDesiredContraction( relDesContr );
    cs1.setMaximalContraction( maxContr );

    cout << "start solver" << endl;
    using namespace std::chrono;
    auto startTime = high_resolution_clock::now();
    auto exact_solution = cs1();

    std::cout << "computation time: "
              << duration_cast< seconds >( high_resolution_clock::now() - startTime ).count()
              << "s." << std::endl;

    std::cout << "Exact grid solution functional value " << t_func1( exact_solution );

    //####################### UNIFORM GRID SOLUTION #######################

    ::Spacy::KaskadeParabolic::GridManager< Spaces > gm2( 11, ::Spacy::Real{30.},
                                                          initialRefinements, FEorder );
    auto domain2 = Spacy::KaskadeParabolic::makeHilbertSpace( gm2, {0u, 1u}, {2u} );

    std::cout << "Creating Functionals" << std::endl;
    auto n_func2 = Spacy::KaskadeParabolic::makeC2Functional( normalFuncGenerator, gm2, domain2 );
    auto t_func2 =
        Spacy::KaskadeParabolic::makeC2Functional( tangentialFuncGenerator, gm2, domain2 );

    cout << "set up solver" << endl;
    // algorithm and parameters
    auto cs2 = Spacy::CompositeStep::AffineCovariantSolver( n_func2, t_func2, domain2 );
    cs2.setRelativeAccuracy( desiredAccuracy );
    cs2.set_eps( eps );
    cs2.setVerbosityLevel( 2 );
    cs2.setMaxSteps( maxSteps );
    cs2.setIterativeRefinements( iterativeRefinements );
    cs2.setDesiredContraction( desContr );
    cs2.setRelaxedDesiredContraction( relDesContr );
    cs2.setMaximalContraction( maxContr );
    auto uniform_solution = cs2();

    std::cout << "Uniform grid solution functional value " << t_func2( uniform_solution )
              << std::endl;

    //####################### EXP GRID SOLUTION #######################

    ::Spacy::KaskadeParabolic::GridManager< Spaces > gm3( 11, ::Spacy::Real{30}, initialRefinements,
                                                          FEorder, "exponential", expgridparam );
    auto domain3 = Spacy::KaskadeParabolic::makeHilbertSpace( gm3, {0u, 1u}, {2u} );

    std::cout << "Creating Functionals" << std::endl;
    auto n_func3 = Spacy::KaskadeParabolic::makeC2Functional( normalFuncGenerator, gm3, domain3 );
    auto t_func3 =
        Spacy::KaskadeParabolic::makeC2Functional( tangentialFuncGenerator, gm3, domain3 );

    cout << "set up solver" << endl;
    // algorithm and parameters
    auto cs3 = Spacy::CompositeStep::AffineCovariantSolver( n_func3, t_func3, domain3 );
    cs3.setRelativeAccuracy( desiredAccuracy );
    cs3.set_eps( eps );
    cs3.setVerbosityLevel( 2 );
    cs3.setMaxSteps( maxSteps );
    cs3.setIterativeRefinements( iterativeRefinements );
    cs3.setDesiredContraction( desContr );
    cs3.setRelaxedDesiredContraction( relDesContr );
    cs3.setMaximalContraction( maxContr );
    auto exp_solution = cs3();

    std::cout << "exponential grid solution functional value " << t_func3( exp_solution )
              << std::endl;

    //####################### PRINT RELATIVE ERROR AND SOLUTION NORM #######################
    auto NormFunc = n_func1.hessian( exact_solution );

    Spacy::KaskadeParabolic::OCP::printNormSolution( exact_solution, NormFunc, gm1, "ex" );
    Spacy::KaskadeParabolic::OCP::printNormSolution(
        uniform_solution, n_func2.hessian( uniform_solution ), gm2, "uni" );
    Spacy::KaskadeParabolic::OCP::printNormSolution( exp_solution, n_func3.hessian( exp_solution ),
                                                     gm3, "exp" );

    std::cout << "distance to uni sol " << std::endl;
    Spacy::KaskadeParabolic::OCP::printErrorDistribution( exact_solution, uniform_solution,
                                                          NormFunc, gm1, "uni" );
    std::cout << "distance to exp sol " << std::endl;
    Spacy::KaskadeParabolic::OCP::printErrorDistribution( exact_solution, exp_solution, NormFunc,
                                                          gm1, "exp" );

    Spacy::KaskadeParabolic::OCP::printNormSolution( exact_solution, NormFunc, gm1, "ex" );
    Spacy::KaskadeParabolic::OCP::printNormSolution(
        uniform_solution, n_func2.hessian( uniform_solution ), gm2, "uni" );
    Spacy::KaskadeParabolic::OCP::printNormSolution( exp_solution, n_func3.hessian( exp_solution ),
                                                     gm3, "exp" );

    ::Spacy::KaskadeParabolic::OCP::writeVTK< Descriptions >( exact_solution, "ex_sol" );

    std::cout << "returning from main" << std::endl;
    return 0;
}
