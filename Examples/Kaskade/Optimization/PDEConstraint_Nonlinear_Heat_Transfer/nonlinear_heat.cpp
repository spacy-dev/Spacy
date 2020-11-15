#define FUSION_MAX_VECTOR_SIZE 15 // NOLINT(cppcoreguidelines-macro-usage)

#include "../../interpolate.h"
#include "nonlinear_control.hh"
#include <functional>
#include <fung/examples/nonlinear_heat.hh>
#include <fung/fung.hh>

#include <Spacy/Adapter/kaskade.hh>
#include <Spacy/Spacy.h>

#include <fem/forEach.hh>
#include <fem/gridmanager.hh>
#include <fem/lagrangespace.hh>
#include <io/vtk.hh>
#include <utilities/gridGeneration.hh>
#include <utilities/kaskopt.hh>

#include <dune/grid/config.h>
#include <dune/grid/uggrid.hh>

#include <chrono>
#include <iostream>

template < class StateVector, class Reference, class ControlVector >
auto trackingTypeCostFunctional( double alpha, const StateVector& y, const Reference& y_ref, const ControlVector& u )
{
    using namespace FunG;
    return finalize( squared( variable< Ids::state >( y ) - variable< Ids::reference >( y_ref ) ) +
                     alpha * squared( variable< Ids::control >( u ) ) );
}

int main( int argc, char* argv[] )
{
    using namespace Kaskade;

    constexpr int dim = 2;
    int silence = 0;
    std::unique_ptr< boost::property_tree::ptree > pt = getKaskadeOptions( argc, argv, silence, false );

    double desiredAccuracy = getParameter( pt, "desiredAccuracy", 1e-6 );
    double eps = getParameter( pt, "eps", 1e-12 );
    double alpha = getParameter( pt, "alpha", 1e-6 );
    int maxSteps = getParameter( pt, "maxSteps", 500 );
    int initialRefinements = getParameter( pt, "initialRefinements", 5 );
    int iterativeRefinements = getParameter( pt, "iterativeRefinements", 0 );
    int FEorder = getParameter( pt, "FEorder", 1 );
    int verbose = getParameter( pt, "verbose", 2 );
    double c = getParameter( pt, "cPara", 1e1 );
    double d = getParameter( pt, "dPara", 1e0 );
    double e = getParameter( pt, "ePara", 0.0 );
    double desContr = getParameter( pt, "desiredContraction", 0.5 );
    double relDesContr = getParameter( pt, "relaxedContraction", desContr + 0.1 );
    double maxContr = getParameter( pt, "maxContraction", 0.75 );

    using std::cout;
    using std::endl;

    cout << "cPara: " << c << " dPara: " << d << " ePara: " << e << " alpha:" << alpha << endl;
    cout << "dContr: " << desContr << " rdContr: " << relDesContr << " mContr: " << maxContr << endl;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
    /* * * * * * * * * * * * * * * * * * * * * * grid generation * * * * * * * * * * * * * * * * * * * * * */
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
    using Grid = Dune::UGGrid< dim >;
    using GridView = Grid::LeafGridView;

    GridManager< Grid > gm( createUnitSquare< Grid >( 1., false ) );
    gm.enforceConcurrentReads( true );
    gm.globalRefine( initialRefinements );

    std::cout << "vertices: " << gm.grid().size( dim ) << std::endl;

    std::cout << std::setprecision( 10 );
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
    /* * * * * * * * * * * * * * * * * * * * * * function spaces * * * * * * * * * * * * * * * * * * * * * */
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
    using H1Space = FEFunctionSpace< ContinuousLagrangeMapper< double, GridView > >;
    using Spaces = boost::fusion::vector< H1Space const* >;
    using PrimalVariables = boost::fusion::vector< Variable< SpaceIndex< 0 >, Components< 1 >, VariableId< Ids::state > >,
                                                   Variable< SpaceIndex< 0 >, Components< 1 >, VariableId< Ids::control > > >;
    using DualVariables = boost::fusion::vector< Variable< SpaceIndex< 0 >, Components< 1 >, VariableId< 0 > > >;
    using VariableDescriptions = boost::fusion::vector< Variable< SpaceIndex< 0 >, Components< 1 >, VariableId< Ids::state > >,
                                                        Variable< SpaceIndex< 0 >, Components< 1 >, VariableId< Ids::control > >,
                                                        Variable< SpaceIndex< 0 >, Components< 1 >, VariableId< Ids::adjoint > > >;
    using Descriptions = VariableSetDescription< Spaces, VariableDescriptions >;
    using PrimalDescriptions = VariableSetDescription< Spaces, PrimalVariables >;
    using DualDescriptions = VariableSetDescription< Spaces, DualVariables >;
    using VarSet = Descriptions::VariableSet;

    // construct variable list.
    std::string names[] = { "y", "u", "p" }; // NOLINT(cppcoreguidelines-avoid-c-arrays)

    auto const& leafView = gm.grid().leafGridView();
    // construct involved spaces.
    H1Space h1Space( gm, leafView, FEorder );
    Spaces spaces( &h1Space );

    Descriptions desc( spaces, names );                           // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay)
    PrimalDescriptions primalDescription( spaces, { "y", "u" } ); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay)
    DualDescriptions dualDescription( spaces, { "p" } );          // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay)

    // Reference solution
    cout << "interpolate" << endl;
    VarSet x_ref( desc );
    interpolateGloballyFromFunctor< PlainAverage >(
        boost::fusion::at_c< Ids::state >( x_ref.data ), []( auto const& cell, auto const& xLocal ) -> Dune::FieldVector< double, 1 > {
            const auto x = cell.geometry().global( xLocal );
            return Dune::FieldVector< double, 1 >( 12 * ( 1 - x[ 1 ] ) * x[ 1 ] * ( 1 - x[ 0 ] ) * x[ 0 ] );
        } );

    cout << "create domain" << endl;
    auto X = Spacy::Kaskade::makeHilbertSpace< Descriptions >( desc, "X" );
    auto Z = Spacy::Kaskade::makeHilbertSpace< PrimalDescriptions >( primalDescription, "Z" );
    auto P = Spacy::Kaskade::makeHilbertSpace< DualDescriptions >( dualDescription, "P" );

    PrimalDescriptions::CoefficientVectorRepresentation<>::type xyz(
        PrimalDescriptions::CoefficientVectorRepresentation<>::init( Spacy::Kaskade::extractSpaces< PrimalDescriptions >( Z ) ) );
    using ZXIdxPairs =
        boost::fusion::vector< Spacy::Kaskade::IdxPair< Ids::state, Ids::state >, Spacy::Kaskade::IdxPair< Ids::control, Ids::control > >;
    Z.setEmbedding( Spacy::Kaskade::ProductSpaceRelation< PrimalDescriptions, Descriptions, ZXIdxPairs >( Z, X ) );
    using XZIdxPairs = ZXIdxPairs;
    Z.setProjection( Spacy::Kaskade::ProductSpaceRelation< Descriptions, PrimalDescriptions, XZIdxPairs >( X, Z ) );
    using PXIdxPairs = boost::fusion::vector< Spacy::Kaskade::IdxPair< 0, Ids::adjoint > >;
    P.setEmbedding( Spacy::Kaskade::ProductSpaceRelation< DualDescriptions, Descriptions, PXIdxPairs >( P, X ) );
    using XPIdxPairs = boost::fusion::vector< Spacy::Kaskade::IdxPair< Ids::adjoint, 0 > >;
    P.setProjection( Spacy::Kaskade::ProductSpaceRelation< Descriptions, DualDescriptions, XPIdxPairs >( X, P ) );
    auto domain = Spacy::Kaskade::makeHilbertSpace< Descriptions >( spaces, { 0u, 1u }, { 2u } );
    // Normal step functional with cg solver
    // auto fn = Spacy::Kaskade::makeLagrangeCGFunctional<stateId,controlId,adjointId>(
    // NormalStepFunctional<stateId,controlId,adjointId,double,Descriptions>(alpha,x_ref,c,d) , domain );

    cout << "create functional" << endl;
    Dune::FieldVector< double, 1 > y0{ 0 };
    Dune::FieldVector< double, 1 > u0{ 0 };
    Dune::FieldVector< double, 1 > y_ref{ 0.5 };
    Dune::FieldMatrix< double, 1, dim > dy0{ 0 };
    auto constraint = FunG::heatModel( c, d, y0, dy0 );
    auto costFunctional = trackingTypeCostFunctional( alpha, y0, y_ref, u0 );

    // Normal step functional with direct solver
    auto normalStepFunctional =
        NormalStepFunctional< Ids, decltype( constraint ), decltype( costFunctional ), Descriptions >( constraint, costFunctional, x_ref );
    auto fn = Spacy::Kaskade::makeC2Functional( std::move( normalStepFunctional ), domain );

    // Lagrange functional
    cout << "make tangential functional " << endl;
    auto tangentialStepFunctional = TangentialStepFunctional< Ids, decltype( constraint ), decltype( costFunctional ), Descriptions >(
        constraint, costFunctional, x_ref );
    auto ft = Spacy::Kaskade::makeC2Functional( std::move( tangentialStepFunctional ), domain );

    cout << "set up solver" << endl;
    // algorithm and parameters
    auto cs = Spacy::CompositeStep::AffineCovariantSolver( fn, ft, domain );
    cs.setRelativeAccuracy( desiredAccuracy );
    cs.setEps( eps );
    cs.setVerbosityLevel( 2 );
    cs.setMaxSteps( maxSteps );
    cs.setIterativeRefinements( iterativeRefinements );
    cs.setDesiredContraction( desContr );
    cs.setRelaxedDesiredContraction( relDesContr );
    cs.setMaximalContraction( maxContr );

    cout << "start solver" << endl;
    using namespace std::chrono;
    auto startTime = high_resolution_clock::now();
    auto result = cs();
    std::cout << "computation time: " << duration_cast< seconds >( high_resolution_clock::now() - startTime ).count() << "s." << std::endl;

    VarSet x( desc );
    Spacy::Kaskade::copy( result, x );

    auto x3 = zero( X );
    Spacy::Kaskade::copy( x, x3 );
    VarSet x2( desc );
    Spacy::Kaskade::copy( x3, x2 );
    // Spacy::Kaskade::writeVTK< VariableDescriptions >( x2, "x2" );

    auto spacy_z = Z.project( x3 );
    PrimalDescriptions::VariableSet z( primalDescription );
    Spacy::Kaskade::copy( spacy_z, z );

    auto spacy_p = P.project( x3 );
    DualDescriptions::VariableSet p( dualDescription );
    Spacy::Kaskade::copy( spacy_p, p );

    auto spacy_z2 = X.embed( spacy_p );
    VarSet x4( desc );
    Spacy::Kaskade::copy( spacy_z2, x4 );

    IoOptions options;
    options.outputType = IoOptions::ascii;
    std::string outfilename( "nonlinear_control" );
    writeVTKFile( gm.grid().leafGridView(), x, outfilename, options, FEorder );
    writeVTKFile( gm.grid().leafGridView(), x_ref, "reference", options, FEorder );
    writeVTKFile( gm.grid().leafGridView(), x4, "x4", options, FEorder );
    writeVTKFile( gm.grid().leafGridView(), z, "primal_variables", options, FEorder );
    writeVTKFile( gm.grid().leafGridView(), p, "dual_variables", options, FEorder );
    // writeVTKFile( gm.grid().leafGridView(), p, "dual_variables", options, FEorder );
}
