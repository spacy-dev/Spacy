#define FUSION_MAX_VECTOR_SIZE 15 // NOLINT(cppcoreguidelines-macro-usage)

#include <chrono>
#include <iostream>
#include <io/vtk.hh>
#include <utilities/kaskopt.hh>
 
#include <Spacy/Spacy.h>
#include <Spacy/VectorSpace.h>
#include <Spacy/Adapter/kaskadeParabolic.hh>
#include <Spacy/Algorithm/PrimalDualProjection/PDP_withoutA.hh>

#include "../optimal_control_withoutfung.hh"

using std::cout;
using std::endl;
using namespace std::chrono;
using namespace Kaskade;
using namespace Spacy::KaskadeParabolic;

constexpr int dim = 2;
constexpr int cp = 1;

int main( int argc, char* argv[] )
{
    cout << std::setprecision( 10 );
    cout << endl;
    
    const auto silence = 0;
    std::unique_ptr< boost::property_tree::ptree > pt = getKaskadeOptions( argc, argv, silence, false );
    
    const auto N                    = getParameter( pt, "stepsInTime", 30 );
    const auto T                    = getParameter( pt, "endTime", 1. );
    
    const auto initialRefinements   = getParameter( pt, "initialRefinements", 5 );
    const auto iterativeRefinements = getParameter( pt, "iterativeRefinements", 0 );
    
    const auto maxSteps             = getParameter( pt, "maxSteps", 1000 );
    const auto FEorder              = getParameter( pt, "FEorder", 1 );
    const auto verbose              = getParameter( pt, "verbose", 1 );
    
    double alpha                    = getParameter( pt, "alpha", 1e-3 );
    double c                        = getParameter( pt, "cPara", 0 );
    double d                        = getParameter( pt, "dPara", 1 );
    double e                        = getParameter( pt, "ePara", 0 );

    const auto stateAcc             = getParameter (pt, "stateAccuracy", 1e-1);
    const auto adjointAcc           = getParameter (pt, "adjointAccuracy", 1e-1);
    const auto modifiedCGAcc        = getParameter (pt, "modifiedCGAccuracy", 1e-2);
    const auto chebyshevAcc         = getParameter (pt, "chebyshevAccuracy", 1e-2);
    const auto desiredAccuracy      = getParameter (pt, "desiredAccuracy", 1e-3);
    
    ///////////////////////////////////////////////// Kaskade FEM-Spaces
    cout << "create FEM-Spaces" << endl;
    
    using Grid = Dune::UGGrid< dim >;
    using GridView = Grid::LeafGridView;
    
    using H1Space = FEFunctionSpace< ContinuousLagrangeMapper< double, GridView > >;
    using Spaces = boost::fusion::vector< FEFunctionSpace< ContinuousLagrangeMapper< double, GridView > > const* >;

    using components = Components<cp>;
    using VariableDescriptions = boost::fusion::vector< Variable< SpaceIndex< 0 >, components, VariableId< 0 > >,
                                                        Variable< SpaceIndex< 0 >, components, VariableId< 1 > >,
                                                        Variable< SpaceIndex< 0 >, components, VariableId< 2 > > >;
    using PrimalVariables = boost::fusion::vector< Variable< SpaceIndex< 0 >, components, VariableId< 0 > >,
                                                Variable< SpaceIndex< 0 >, components, VariableId< 1 > > >;
    using StateVariables = boost::fusion::vector< Variable< SpaceIndex< 0 >, components, VariableId< 0 > > >;
    using ControlVariables = boost::fusion::vector< Variable< SpaceIndex< 0 >, components, VariableId< 0 > > >;
    using DualVariables = boost::fusion::vector< Variable< SpaceIndex< 0 >, components, VariableId< 0 > > >;

    using Descriptions = VariableSetDescription< Spaces, VariableDescriptions >;
    using PrimalDescriptions = VariableSetDescription< Spaces, PrimalVariables >;
    using DualDescriptions = VariableSetDescription< Spaces, DualVariables >;
    using StateDescriptions = VariableSetDescription< Spaces, StateVariables >;
    using ControlDescriptions = VariableSetDescription< Spaces, ControlVariables >;
    using VarSet = Descriptions::VariableSet;
    
            //////////////////////////////// Grid
    Spacy::KaskadeParabolic::GridManager< Spaces > gm( N, ::Spacy::Real{ T }, initialRefinements, FEorder );
    
    ///////////////////////////////////////////////// Spacy Vector-Spaces
    cout << "create Vector-Spaces" << endl;
    
    auto totalSpace_ = makeHilbertSpace< Descriptions >( gm, {"y","u","p"}, "X" );
    
    auto primalSpace_ = makeHilbertSpace< PrimalDescriptions >( gm, {"y","u"}, "Z" );
    using ZXIdxMap = boost::fusion::vector< IdxPair<0,0>, IdxPair<1,1> >;
    setSubSpaceRelation< PrimalDescriptions, Descriptions, ZXIdxMap >( primalSpace_, totalSpace_ );
    
    auto adjointSpace_ = makeHilbertSpace< DualDescriptions >( gm, {"p"}, "P" );
    using PXIdxMap = boost::fusion::vector< IdxPair<0,2> >;
    setSubSpaceRelation< DualDescriptions, Descriptions, PXIdxMap >( adjointSpace_, totalSpace_ );
    
    auto stateSpace_ = makeHilbertSpace< StateDescriptions >( gm, {"y"}, "Y" );
    using YXIdxMap = boost::fusion::vector< IdxPair<0,0> >;
    setSubSpaceRelation< StateDescriptions, Descriptions, YXIdxMap >( stateSpace_, totalSpace_ );
    using YZIdxMap = boost::fusion::vector< IdxPair<0,0> >;
    setSubSpaceRelation< StateDescriptions, PrimalDescriptions, YZIdxMap >( stateSpace_, primalSpace_ );
    
    auto controlSpace_ = makeHilbertSpace< ControlDescriptions >( gm, {"u"}, "U" );
    using UXIdxMap = boost::fusion::vector< IdxPair<0,1> >;
    setSubSpaceRelation< ControlDescriptions, Descriptions, UXIdxMap >( controlSpace_, totalSpace_ );
    using UZIdxMap = boost::fusion::vector< IdxPair<0,1> >;
    setSubSpaceRelation< ControlDescriptions, PrimalDescriptions, UZIdxMap >( controlSpace_, primalSpace_ );
    
    ///////////////////////////////////////////////// Functional 
    cout << "create Functional" << endl;
    
    using FunctionalDefinition = TangentialStepFunctional< 0, 1, 2, double, Descriptions >;
    std::function< FunctionalDefinition( const typename Descriptions::VariableSet ) > FuncGenerator =
        [ &c, &d, &e, &alpha ]( const typename Descriptions::VariableSet y_ref ) {
            return FunctionalDefinition( alpha, y_ref, c, d, e );
        };

    auto fL = makeC2Functional<FunctionalDefinition>( FuncGenerator, gm, totalSpace_ ); fL.setVerbosity( verbose > 1 ? true : false );
    
    cout << "Assembling" << endl;
    auto x0 = zero(totalSpace_);
    fL.hessian(x0);
    auto b = fL.d1(x0);
    
    ///////////////////////////////////////////////// Blocks, Solvers
    cout << "create Blocks and Solvers" << endl;
    
    using KMAT = ::Kaskade::MatrixAsTriplet<double>;
    
    auto A                  = makeCallableBlock<FunctionalDefinition>(fL,stateSpace_,adjointSpace_);
    auto AT                 = makeCallableBlock<FunctionalDefinition>(fL,adjointSpace_,stateSpace_);
    auto minusB             = makeCallableBlock<FunctionalDefinition>(fL,controlSpace_,adjointSpace_);
    auto M                  = makeSymmetricCallableBlock<FunctionalDefinition>(fL,primalSpace_);
    auto Mu                 = makeSymmetricCallableBlock<FunctionalDefinition>(fL,controlSpace_);
    
    auto D                  = makeMassMatrixBlock<StateDescriptions,DualDescriptions>(fL,stateSpace_,adjointSpace_);
    auto DT                 = makeMassMatrixBlock<DualDescriptions,StateDescriptions>(fL,adjointSpace_,stateSpace_);

    auto Prec                = makeSingleBlockPreconditioner<DualDescriptions,StateDescriptions,FunctionalDefinition>(adjointSpace_,stateSpace_,fL );
    auto PrecT               = makeSingleBlockPreconditioner<StateDescriptions,DualDescriptions,FunctionalDefinition>(stateSpace_,adjointSpace_,fL,true );
                                
    auto controlSolver_     = makeDirectSolver<ControlDescriptions,ControlDescriptions>(Mu,gm,controlSpace_, controlSpace_);
    auto stateSolver_       = makePDESolver(D,A,Prec,gm.getTempGrid().getDtVec(),gm,stateAcc);
    auto adjointSolver_     = makePDESolver(DT,AT,PrecT,gm.getTempGrid().getDtVec(),gm,adjointAcc,false);
    auto surrStateSolver_   = makeSurrogatePDESolver<StateDescriptions,DualDescriptions>(A,Prec,gm.getTempGrid().getDtVec());
    auto surrAdjointSolver_ = makeSurrogatePDESolver<DualDescriptions,StateDescriptions>(AT,PrecT,gm.getTempGrid().getDtVec(),false);

    ///////////////////////////////////////////////// Execute
    cout << "Execute" << endl;

    ::Spacy::PDP_withoutA::Solver pdp(M,stateSolver_,adjointSolver_,minusB,surrStateSolver_,surrAdjointSolver_,controlSolver_,totalSpace_);
    
    pdp.setMaxSteps(maxSteps);
    pdp.setRelativeAccuracy(desiredAccuracy);
    pdp.setInternalAccuracies(stateAcc,adjointAcc,modifiedCGAcc,chebyshevAcc);
    
    
    auto startTime = high_resolution_clock::now();
    
    auto res = totalSpace_.embed( pdp(b) );                         // Start PDP
    
    auto endTime = high_resolution_clock::now() - startTime;
    cout << "computation time: " << duration_cast< seconds >( endTime ).count() << "s, " 
            << duration_cast< milliseconds >( endTime ).count() % 1000 << "ms." << endl << endl;
    
    ///////////////////////////////////////////////// Results
    cout << "finishing ... " << std::flush;

    auto x = totalSpace_.embed(res);
    
    ///////////////////////////////////////////////// VTK-Files
    
    PDE::writeVTK<Descriptions>(x,"test");
    
    cout << "done" << endl;
    
    /**/
    return 0;
}
