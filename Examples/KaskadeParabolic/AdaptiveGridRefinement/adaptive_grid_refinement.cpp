#define FUSION_MAX_VECTOR_SIZE 15

#include <chrono>
#include <iostream>
#include <functional>

#include <dune/grid/config.h>
#include <dune/grid/uggrid.hh>

#include <Spacy/Adapter/kaskadeParabolic.hh>
#include <Spacy/Algorithm/CompositeStep/affineCovariantSolver.hh>
#include <Spacy/Spaces/ScalarSpace/Real.h>

#include <fem/lagrangespace.hh>
#include <utilities/kaskopt.hh>

#include "Spacy/Adapter/KaskadeParabolic/Constraints/lqOcpHeader.hh"
#include <Spacy/inducedScalarProduct.hh>

using namespace Kaskade;
int main( int argc, char* argv[] )
{
    static constexpr int stateId = 0;
    static constexpr int controlId = 1;
    static constexpr int adjointId = 2;

    constexpr int dim = 2;
    int silence = 0;
    std::unique_ptr< boost::property_tree::ptree > pt =
        getKaskadeOptions( argc, argv, silence, false );

    double alpha = getParameter( pt, "alpha", 1. );
    double d = getParameter( pt, "dPara", 0.1 );
    double mu = getParameter( pt, "muPara", 0.0 );
    double tolerance = getParameter( pt, "tolerance", 1e-4 );

    using std::cout;
    using std::endl;

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
        normalFuncGenerator = [&d, &mu, &alpha]( const typename Descriptions::VariableSet y_ref ) {
            return NormalFunctionalDefinition( alpha, y_ref, d, mu );
        };

    using TangentialFunctionalDefinition =
        TangentialStepFunctional< stateId, controlId, adjointId, double, Descriptions >;
    std::function< TangentialFunctionalDefinition( const typename Descriptions::VariableSet ) >
        tangentialFuncGenerator =
            [&d, &mu, &alpha]( const typename Descriptions::VariableSet y_ref ) {
                return TangentialFunctionalDefinition( alpha, y_ref, d, mu );
            };

    ::Spacy::KaskadeParabolic::GridManager< Spaces > gm( 2, Spacy::Real{30.} );
    auto domain = Spacy::KaskadeParabolic::makeHilbertSpace( gm, {0u, 1u}, {2u} );

    std::cout << "Creating Functionals" << std::endl;
    ::Spacy::C2Functional t_func =
        Spacy::KaskadeParabolic::makeC2Functional( tangentialFuncGenerator, gm, domain );

    gm.getTempGrid().print();
    using namespace ::Spacy::KaskadeParabolic;

    ::Spacy::C2Functional n_func =
        Spacy::KaskadeParabolic::makeC2Functional( normalFuncGenerator, gm, domain );
    auto& n_func_impl =
        ::Spacy::cast_ref<::Spacy::KaskadeParabolic::C2Functional< NormalFunctionalDefinition > >(
            n_func );

    ::Spacy::CompositeStep::AffineCovariantSolver cs( n_func, t_func, domain );
    cs.setVerbosity( true );

    auto dx = cs();
    OCP::printNormSolution< Descriptions >( dx, n_func_impl.hessian( zero( domain ) ), gm,
                                            "sol" + std::to_string( 0 ) );

    OCP::ErrorEstimator errest( tolerance, 1. );

    // errorestimation with cost functional
    errest.estimateErrCostFunc< TangentialFunctionalDefinition >( dx, t_func );

    // error estimation with quantity of interest
    //  errest.estimateErr<TangentialFunctionalDefinition>(dx,t_func);

    auto ctr = 1u;
    do
    {
        errest.refine< TangentialFunctionalDefinition, NormalFunctionalDefinition >( domain, t_func,
                                                                                     n_func );
        ::Spacy::CompositeStep::AffineCovariantSolver cs2( n_func, t_func, domain );
        dx = cs2();
        OCP::printNormSolution< Descriptions >( dx, n_func_impl.hessian( zero( domain ) ), gm,
                                                "sol" + std::to_string( ctr ) );
        ctr++;
    } while ( errest.estimateErrCostFunc< TangentialFunctionalDefinition >( dx, t_func ) );
    //    while(errest.estimateErr<TangentialFunctionalDefinition>(dx,t_func));

    std::cout << "returning from main" << std::endl;
    return 0;
}
