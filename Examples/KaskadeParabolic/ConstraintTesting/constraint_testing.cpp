#define FUSION_MAX_VECTOR_SIZE 15

#include <chrono>
#include <iostream>
#include <functional>

#include <dune/grid/config.h>
#include <dune/grid/uggrid.hh>

#include <Spacy/Adapter/kaskadeParabolic.hh>
#include <Spacy/Algorithm/CompositeStep/AffineCovariantSolver.h>
#include <Spacy/Spaces/ProductSpace.h>
#include <Spacy/Spaces/ScalarSpace/Real.h>

#include <fem/forEach.hh>
#include <fem/gridmanager.hh>
#include <fem/lagrangespace.hh>
#include <io/vtk.hh>
#include <utilities/gridGeneration.hh>
#include <utilities/kaskopt.hh>

#define NCOMPONENTS 1
// DISTRIBUTED CONTROL
//#include "Spacy/Adapter/KaskadeParabolic/Constraints/nonlinOcpHeader.hh"

// BOUNDARY CONTROL
#include "Spacy/Adapter/KaskadeParabolic/Constraints/nonlinOcpHeader_boundaryControl.hh"
#include <fem/boundaryspace.hh>

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
    double alpha = getParameter( pt, "alpha", 0.01 );
    int maxSteps = getParameter( pt, "maxSteps", 500 );
    int initialRefinements = getParameter( pt, "initialRefinements", 4 );
    int iterativeRefinements = getParameter( pt, "iterativeRefinements", 0 );
    int FEorder = getParameter( pt, "FEorder", 1 );
    int verbose = getParameter( pt, "verbose", 2 );
    double c = getParameter( pt, "cPara", 0 );
    double d = getParameter( pt, "dPara", 1 );
    double e = getParameter( pt, "ePara", 0 );
    double mu = getParameter( pt, "muPara", 0 );
    double desContr = getParameter( pt, "desiredContraction", 0.5 );
    double relDesContr = getParameter( pt, "relaxedContraction", desContr + 0.1 );
    double maxContr = getParameter( pt, "maxContraction", 0.75 );
    double expgridparam = getParameter( pt, "expGrid", -0.345 );

    using std::cout;
    using std::endl;
    cout << "cPara: " << c << " dPara: " << d << " ePara: " << e << " muPara: " << mu
         << " alpha:" << alpha << endl;

    using Grid = Dune::UGGrid< dim >;

    // if boundary control
    using L2Space =
        FEFunctionSpace< BoundaryMapper< ContinuousLagrangeMapper, double, Grid::LeafGridView > >;
    // if distributed control
    //  using L2Space = FEFunctionSpace<ContinuousLagrangeMapper<double,Grid::LeafGridView> >;

    using H1Space = FEFunctionSpace< ContinuousLagrangeMapper< double, Grid::LeafGridView > >;
    using Spaces = boost::fusion::vector< H1Space const*, L2Space const* >;
    using VY = VariableDescription< 0, 1, stateId >;
    using VU = VariableDescription< 1, 1, controlId >;
    using VP = VariableDescription< 0, 1, adjointId >;
    using VariableDescriptions = boost::fusion::vector< VY, VU, VP >;
    using Descriptions = VariableSetDescription< Spaces, VariableDescriptions >;

    ::Spacy::KaskadeParabolic::GridManager< Spaces > gm( 11, ::Spacy::Real{3.}, initialRefinements,
                                                         FEorder );
    auto domain = Spacy::KaskadeParabolic::makeHilbertSpace( gm, {0u, 1u}, {2u} );

    // if boundary control
    using NormalFunctionalDefinition =
        NormalBoundaryControlStepFunctional< stateId, controlId, adjointId, double, Descriptions >;
    // if distributed control
    //  using NormalFunctionalDefinition =
    //  NormalStepFunctional<stateId,controlId,adjointId,double,Descriptions>;
    std::function< NormalFunctionalDefinition( const typename Descriptions::VariableSet ) >
        normalFuncGenerator =
            [&c, &d, &e, &mu, &alpha]( const typename Descriptions::VariableSet y_ref ) {
                return NormalFunctionalDefinition( alpha, y_ref, c, d, e, mu );
            };

    // if boundary control
    using TangentialFunctionalDefinition =
        TangentialBoundaryControlStepFunctional< stateId, controlId, adjointId, double,
                                                 Descriptions >;
    // if distributed control
    //  using TangentialFunctionalDefinition =
    //  TangentialStepFunctional<stateId,controlId,adjointId,double,Descriptions>;
    std::function< TangentialFunctionalDefinition( const typename Descriptions::VariableSet ) >
        tangentialFuncGenerator =
            [&c, &d, &e, &mu, &alpha]( const typename Descriptions::VariableSet y_ref ) {
                return TangentialFunctionalDefinition( alpha, y_ref, c, d, e, mu );
            };

    //####################### OCP SOLUTION #######################
    std::cout << "Creating Functionals" << std::endl;
    auto n_func2 = Spacy::KaskadeParabolic::makeC2Functional( normalFuncGenerator, gm, domain );
    auto t_func2 = Spacy::KaskadeParabolic::makeC2Functional( tangentialFuncGenerator, gm, domain );

    typename Descriptions::VariableSet x0(
        Descriptions( *gm.getSpacesVec().at( 0 ), {"y", "u", "p"} ) );
    interpolateGloballyWeak< PlainAverage >(
        boost::fusion::at_c< 0 >( x0.data ),
        makeWeakFunctionView(
            []( auto const& cell, auto const& xLocal ) -> Dune::FieldVector< double, 1 > {
                auto x = cell.geometry().global( xLocal );

                auto returnval = Dune::FieldVector< double, 1 >( +12 * ( 1 - x[ 1 ] ) * x[ 1 ] *
                                                                 ( 1 - x[ 0 ] ) * x[ 0 ] );
                if ( x[ 0 ] >= 1 && x[ 0 ] <= 3 )
                    return 0 * returnval;
                return returnval;
            } ) );

    t_func2.setInitialCondition( boost::fusion::at_c< 0 >( x0.data ).coefficients() );

    cout << "set up solver" << endl;
    // algorithm and parameters
    auto cs2 = Spacy::CompositeStep::AffineCovariantSolver( n_func2, t_func2, domain );
    cs2.setRelativeAccuracy( desiredAccuracy );
    cs2.set_eps( eps );
    cs2.setVerbosityLevel( 2 );
    cs2.setMaxSteps( maxSteps );
    cs2.setIterativeRefinements( iterativeRefinements );
    cs2.setDesiredContraction( desContr );
    cs2.setRelaxedDesiredContraction( relDesContr );
    cs2.setMaximalContraction( maxContr );
    auto uniform_solution = cs2();
    ::Spacy::KaskadeParabolic::OCP::writeVTK< Descriptions >( uniform_solution, "sol" );

    //####################### PDE SOLUTION #######################

    // if boundary control
    using HeatFunctionalDefinition = BoundaryControlledNonlinearModelPDE<
        double, ::Spacy::KaskadeParabolic::Detail::ExtractDescription_t< Descriptions, 0 >,
        ::Spacy::KaskadeParabolic::Detail::ExtractDescription_t< Descriptions, 1 > >;
    // if distributed control
    //  using HeatFunctionalDefinition =
    //  NonlinearModelPDE<double,::Spacy::KaskadeParabolic::Detail::ExtractDescription_t<Descriptions,0>
    //  , ::Spacy::KaskadeParabolic::Detail::ExtractDescription_t<Descriptions,1>  >;
    std::function< HeatFunctionalDefinition(
        const typename ::Spacy::KaskadeParabolic::Detail::ExtractDescription_t< Descriptions,
                                                                                1 >::VariableSet ) >
        forwardFunctionalGenerator = [&c, &d, &e, &mu](
            const typename ::Spacy::KaskadeParabolic::Detail::ExtractDescription_t<
                Descriptions, 1 >::VariableSet control ) {
            return HeatFunctionalDefinition( c, d, e, mu, control );
        };

    auto domain_forward_ = Spacy::KaskadeParabolic::makeHilbertSpace( gm );
    auto z = zero( domain );
    auto z_ps = ::Spacy::cast_ref<::Spacy::ProductSpace::Vector >( z );
    ::Spacy::C1Operator A = ::Spacy::KaskadeParabolic::makeC1Operator< HeatFunctionalDefinition >(
        forwardFunctionalGenerator, gm, domain_forward_, domain_forward_.dualSpace(),
        ::Spacy::cast_ref<::Spacy::ProductSpace::Vector >( z_ps.component( 0 ) ).component( 1 ) );
    domain_forward_.setScalarProduct( Spacy::InducedScalarProduct(
        ::Spacy::cast_ref<::Spacy::KaskadeParabolic::C1Operator< HeatFunctionalDefinition > >( A )
            .massMatrix() ) );
    ::Spacy::cast_ref<::Spacy::KaskadeParabolic::C1Operator< HeatFunctionalDefinition > >( A )
        .setInitialCondition( boost::fusion::at_c< 0 >( x0.data ).coefficients() );

    auto result_ps = ::Spacy::cast_ref<::Spacy::ProductSpace::Vector >( uniform_solution );
    ::Spacy::cast_ref<::Spacy::KaskadeParabolic::C1Operator< HeatFunctionalDefinition > >( A )
        .setSource( (::Spacy::cast_ref<::Spacy::ProductSpace::Vector >( result_ps.component( 0 ) ) )
                        .component( 1 ) );

    auto p = Spacy::Newton::Parameter{};
    p.setRelativeAccuracy( 1e-12 );
    auto y = covariantNewton( A, p );

    using DescriptionsY =
        ::Spacy::KaskadeParabolic::Detail::ExtractDescription_t< Descriptions, 0 >;
    using CoefficientVectorY =
        typename DescriptionsY::template CoefficientVectorRepresentation<>::type;

    //####################### Compute Difference #######################

    ::Spacy::KaskadeParabolic::PDE::writeVTK< DescriptionsY >( y, "sol_forward" );

    auto y1 = (::Spacy::cast_ref<::Spacy::ProductSpace::Vector >( result_ps.component( 0 ) ) )
                  .component( 0 );
    auto y1impl = ::Spacy::cast_ref<::Spacy::KaskadeParabolic::Vector< DescriptionsY > >( y1 );
    auto y2impl = ::Spacy::cast_ref<::Spacy::KaskadeParabolic::Vector< DescriptionsY > >( y );

    for ( auto i = 0u; i < gm.getSpacesVec().size(); i++ )
    {
        CoefficientVectorY y1coeff( Descriptions::template CoefficientVectorRepresentation<>::init(
            *gm.getSpacesVec().at( i ) ) );
        CoefficientVectorY y2coeff( Descriptions::template CoefficientVectorRepresentation<>::init(
            *gm.getSpacesVec().at( i ) ) );

        boost::fusion::at_c< 0 >( y1coeff.data ) = y1impl.getCoeffVec( i );
        boost::fusion::at_c< 0 >( y2coeff.data ) = y2impl.getCoeffVec( i );

        y1coeff -= y2coeff;

        auto My = y1coeff;
        ::Spacy::cast_ref<::Spacy::KaskadeParabolic::C1Operator< HeatFunctionalDefinition > >( A )
            .linearization( zero( domain_forward_ ) )
            .getKaskOp( "MassY", i )
            .apply( y1coeff, My );
        std::cout << "norm " << i << y1coeff * My << std::endl;
    }

    std::cout << "returning from main" << std::endl;
    return 0;
}
