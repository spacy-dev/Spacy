#define FUSION_MAX_VECTOR_SIZE 15 // NOLINT(cppcoreguidelines-macro-usage)

#include <chrono>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <io/vtk.hh>
#include <utilities/kaskopt.hh>

#include <valgrind/callgrind.h>
 
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
constexpr int nT = 1;

int main( int argc, char* argv[] )
{
    cout << std::setprecision( 4 );
    cout << endl;
    
    const auto silence = 0;
    std::unique_ptr< boost::property_tree::ptree > pt = getKaskadeOptions( argc, argv, silence, false );
    
    const auto N                    = getParameter( pt, "stepsInTime", 200 );
    const auto Y                    = getParameter( pt, "stepsInTime", std::max(N/4,3) );
    const auto T                    = getParameter( pt, "endTime", 1. );
    
    const auto initialRefinements   = getParameter( pt, "initialRefinements", 3 );
    const auto iterativeRefinements = getParameter( pt, "iterativeRefinements", 0 ); //aktuell unwichtig
    
    const auto maxSteps             = getParameter( pt, "maxSteps", 5000 );
    const auto FEorder              = getParameter( pt, "FEorder", 1 );
    const auto verbose              = getParameter( pt, "verbose", 1 );
    
    double alpha                    = getParameter( pt, "alpha", 1e-2 );
    double c                        = getParameter( pt, "cPara", 0 );
    double d                        = getParameter( pt, "dPara", 1 );
    double e                        = getParameter( pt, "ePara", 0 );

    const auto modifiedCGAcc        = getParameter (pt, "modifiedCGAccuracy", 1e-1);
    const auto desiredAccuracy      = getParameter (pt, "desiredAccuracy", 1e-5);
    
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
    unsigned PDPiterations = 0, PPCGiterations = 0, MGiterations = 0;
    
    Spacy::Vector res;
        
    Spacy::KaskadeParabolic::GridManager< Spaces > gm( N, ::Spacy::Real{ T }, initialRefinements, dim, FEorder );
    
    auto startTime = high_resolution_clock::now();
    auto endTime = high_resolution_clock::now() - startTime;
    int inner_counter = 0;
        
    Spacy::KaskadeParabolic::GridManager< Spaces > gmy( Y, ::Spacy::Real{ T }, initialRefinements, dim, FEorder );
    
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

    auto totalSpace_y = makeHilbertSpace< Descriptions >( gmy, {"y","u","p"}, "Xs" );
    using XXIdxMap = boost::fusion::vector< IdxPair<0,0>, IdxPair<1,1>, IdxPair<2,2> >;
    setSubSpaceRelation< Descriptions, Descriptions, XXIdxMap >( totalSpace_y, totalSpace_ );
    
    auto primalSpace_y = makeHilbertSpace< PrimalDescriptions >( gmy, {"y","u"}, "Zs" );
    using ZZIdxMap = boost::fusion::vector< IdxPair<0,0>, IdxPair<1,1> >;
    setSubSpaceRelation< PrimalDescriptions, PrimalDescriptions, ZZIdxMap >( primalSpace_y, primalSpace_ );
    setSubSpaceRelation< PrimalDescriptions, Descriptions, ZXIdxMap >( primalSpace_y, totalSpace_y );
    setSubSpaceRelation< PrimalDescriptions, Descriptions, ZXIdxMap >( primalSpace_y, totalSpace_ );
    
    auto adjointSpace_y = makeHilbertSpace< DualDescriptions >( gmy, {"p"}, "Ps" );
    using PPIdxMap = boost::fusion::vector< IdxPair<0,0> >;
    setSubSpaceRelation< DualDescriptions, DualDescriptions, PPIdxMap >( adjointSpace_y, adjointSpace_ );
    setSubSpaceRelation< DualDescriptions, Descriptions, PXIdxMap >( adjointSpace_y, totalSpace_ );
    setSubSpaceRelation< DualDescriptions, Descriptions, PXIdxMap >( adjointSpace_y, totalSpace_y );
    
    auto stateSpace_y = makeHilbertSpace< StateDescriptions >( gmy, {"y"}, "Ys" );
    using YYIdxMap = boost::fusion::vector< IdxPair<0,0> >;
    setSubSpaceRelation< StateDescriptions, StateDescriptions, YYIdxMap >( stateSpace_y, stateSpace_ );
    setSubSpaceRelation< StateDescriptions, Descriptions, YXIdxMap >( stateSpace_y, totalSpace_y );
    setSubSpaceRelation< StateDescriptions, PrimalDescriptions, YZIdxMap >( stateSpace_y, primalSpace_y );
    setSubSpaceRelation< StateDescriptions, Descriptions, YXIdxMap >( stateSpace_y, totalSpace_ );
    setSubSpaceRelation< StateDescriptions, PrimalDescriptions, YZIdxMap >( stateSpace_y, primalSpace_ );
    
    auto controlSpace_y = makeHilbertSpace< ControlDescriptions >( gmy, {"u"}, "Us" );
    using UUIdxMap = boost::fusion::vector< IdxPair<0,0> >;
    setSubSpaceRelation< ControlDescriptions, ControlDescriptions, UUIdxMap >( controlSpace_y, controlSpace_ );
    setSubSpaceRelation< ControlDescriptions, Descriptions, UXIdxMap >( controlSpace_y, totalSpace_ );
    setSubSpaceRelation< ControlDescriptions, PrimalDescriptions, UZIdxMap >( controlSpace_y, primalSpace_ );
    setSubSpaceRelation< ControlDescriptions, Descriptions, UXIdxMap >( controlSpace_y, totalSpace_y );
    setSubSpaceRelation< ControlDescriptions, PrimalDescriptions, UZIdxMap >( controlSpace_y, primalSpace_y );
    
    ///////////////////////////////////////////////// Functional 
    cout << "create Functional" << endl;
    
    using FunctionalDefinition = TangentialStepFunctional< 0, 1, 2, double, Descriptions >;
    std::function< FunctionalDefinition( const typename Descriptions::VariableSet ) > FuncGenerator =
        [ &c, &d, &e, &alpha ]( const typename Descriptions::VariableSet y_ref ) {
            return FunctionalDefinition( alpha, y_ref, c, d, e );
        };
    

    auto fL = makeC2Functional<FunctionalDefinition>( FuncGenerator, gm, totalSpace_ ); fL.setVerbosity( verbose > 1 ? true : false );
    auto fLs = makeC2Functional<FunctionalDefinition>( FuncGenerator, gmy, totalSpace_y ); fLs.setVerbosity( verbose > 1 ? true : false );
    

    cout << "Assembling" << endl;
    
    auto x0 = zero(totalSpace_);
    auto x0s = zero(totalSpace_y);
    
    fL.hessian(x0);
    fLs.hessian(x0s);
    
    auto b = fL.d1(x0);
    auto bs = fLs.d1(x0s);
    
    auto y_ref = zero(stateSpace_);
    for(int i = 0; i<N; i++)
    {
        auto& y_refVS = Spacy::cast_ref<Vector<StateDescriptions>>(y_ref).get_nonconst(i);
        ::Kaskade::interpolateGloballyWeak<::Kaskade::PlainAverage >(
            ::boost::fusion::at_c< 0 >( y_refVS.data ),
            ::Kaskade::makeWeakFunctionView( []( auto const& cell, auto const& xLocal )
                                                -> ::Dune::FieldVector< double, 1 > {
                auto x = cell.geometry().global( xLocal );
                double ref = 12.;
                for(int j = 0; j < dim; j++)
                {
                    ref *= ( 1 - x[ j ] ) * x[ j ];
                }
                return Dune::FieldVector< double, 1 >( ref );
            } ) );
    }
    
    ///////////////////////////////////////////////// Blocks, Solvers
    cout << "create Blocks and Solvers" << endl;
    
    auto As                  = makeCallableBlock<FunctionalDefinition>(fLs,stateSpace_y,adjointSpace_y,true,nT);
    auto ATs                 = makeCallableBlock<FunctionalDefinition>(fLs,adjointSpace_y,stateSpace_y,true,nT);
    auto minusBs             = makeCallableBlock<FunctionalDefinition>(fLs,controlSpace_y,adjointSpace_y,true,nT);
    auto Ms                  = makeSymmetricCallableBlock<FunctionalDefinition>(fLs,primalSpace_y,true,nT);
    auto Mys                 = makeSymmetricCallableBlock<FunctionalDefinition>(fLs,stateSpace_y,true,nT);
    auto Mus                 = makeSymmetricCallableBlock<FunctionalDefinition>(fLs,controlSpace_y,true,nT);

    auto Precs               = makeSingleBlockPreconditioner<DualDescriptions,StateDescriptions,FunctionalDefinition>
                                (adjointSpace_y,stateSpace_y,fLs );
    auto PrecTs              = makeSingleBlockPreconditioner<StateDescriptions,DualDescriptions,FunctionalDefinition>
                                (stateSpace_y,adjointSpace_y,fLs,true );
                               
    auto controlSolver_y     = makeDirectSolver<ControlDescriptions,ControlDescriptions>(Mus,fLs,controlSpace_y, controlSpace_y);
    auto surrStateSolver_y   = makeSurrogatePDESolver<StateDescriptions,DualDescriptions>(As,Precs,gmy);
    auto surrAdjointSolver_y = makeSurrogatePDESolver<DualDescriptions,StateDescriptions>(ATs,PrecTs,gmy,false);

    auto A                   = makeCallableBlock<FunctionalDefinition>(fL,stateSpace_,adjointSpace_,true,nT);
    auto AT                  = makeCallableBlock<FunctionalDefinition>(fL,adjointSpace_,stateSpace_,true,nT);
    auto minusB              = makeCallableBlock<FunctionalDefinition>(fL,controlSpace_,adjointSpace_,true,nT);
    auto M                   = makeSymmetricCallableBlock<FunctionalDefinition>(fL,primalSpace_,true,nT);
    auto My                  = makeSymmetricCallableBlock<FunctionalDefinition>(fL,stateSpace_,true,nT);
    auto Mu                  = makeSymmetricCallableBlock<FunctionalDefinition>(fL,controlSpace_,true,nT);
    
    auto D                   = makeCallableBlock<FunctionalDefinition>(fL,stateSpace_,adjointSpace_,false,nT);
    auto DT                  = makeCallableBlock<FunctionalDefinition>(fL,adjointSpace_,stateSpace_,false,nT);
    
    auto Prec                = makeSingleBlockPreconditioner<DualDescriptions,StateDescriptions,FunctionalDefinition>
                                 (adjointSpace_,stateSpace_,fL );
    auto PrecT               = makeSingleBlockPreconditioner<StateDescriptions,DualDescriptions,FunctionalDefinition>
                                 (stateSpace_,adjointSpace_,fL,true );
                                
    auto controlSolver_      = makeDirectSolver<ControlDescriptions,ControlDescriptions>(Mu,fL,controlSpace_, controlSpace_);
    auto surrStateSolver_    = makeSurrogatePDESolver<StateDescriptions,DualDescriptions>(A,Prec,gm);
    auto surrAdjointSolver_  = makeSurrogatePDESolver<DualDescriptions,StateDescriptions>(AT,PrecT,gm,false);
    auto stateSolver_        = makePDESolver<StateDescriptions,DualDescriptions>(D,A,Prec,fL,0.1);
    auto adjointSolver_      = makePDESolver<DualDescriptions,StateDescriptions>(DT,AT,PrecT,fL,0.1,false);
    
    auto Diag                = makeDiagHandler<StateDescriptions,DualDescriptions>(adjointSolver_(y_ref),y_ref);
    auto diagStateSolver_    = makeSurrogateOperator<DualDescriptions,StateDescriptions>
                                (stateSolver_,Prec,Diag,adjointSpace_,stateSpace_);
    auto diagAdjointSolver_  = makeSurrogateOperator<StateDescriptions,DualDescriptions>
                                (adjointSolver_,PrecT,Diag,stateSpace_,adjointSpace_, false);
    
    ///////////////////////////////////////////////// Execute
    cout << "Execute" << endl;
    
    // kein Ersatzoperator
//     ::Spacy::PDP_withoutA::Solver pdp(M,My,Mu,stateSolver_,adjointSolver_,minusB,minusB,stateSolver_,adjointSolver_,controlSolver_,controlSolver_,totalSpace_,totalSpace_);
    
    // normaler Ersatzoperator
//     ::Spacy::PDP_withoutA::Solver pdp(M,My,Mu,stateSolver_,adjointSolver_,minusB,minusB,surrStateSolver_,surrAdjointSolver_,controlSolver_,controlSolver_,totalSpace_,totalSpace_);
    
    // Ersatzoperator mit grobem Gitter
    ::Spacy::PDP_withoutA::Solver pdp(M,Mys,Mus,stateSolver_,adjointSolver_,minusB,minusBs,surrStateSolver_y,surrAdjointSolver_y,controlSolver_,controlSolver_y,totalSpace_,totalSpace_y);
    
    // Diagonal-Ersatzoperator
//     ::Spacy::PDP_withoutA::Solver pdp(M,My,Mu,stateSolver_,adjointSolver_,minusB,minusB,diagStateSolver_,diagAdjointSolver_,controlSolver_,controlSolver_,totalSpace_,totalSpace_);
    
    pdp.setMaxSteps(maxSteps);
    pdp.setRelativeAccuracy(desiredAccuracy);
    pdp.setInternalAccuracies(0.1,0.1,modifiedCGAcc,0.1);
    pdp.setVerbosity(true);
    
    startTime = high_resolution_clock::now();
    res = totalSpace_.embed(pdp( b ));                         // Start PDP
    endTime = high_resolution_clock::now() - startTime;
    
    PDPiterations = pdp.getIterations();
    PPCGiterations = pdp.getPPCGIterations();
    MGiterations = pdp.getMGIterations();

    cout << "computation time Jacobi" << Y << ": " << duration_cast< seconds >( endTime ).count() << "s, " 
        << duration_cast< milliseconds >( endTime ).count() % 1000 << "ms." << endl;
    cout << "davon PPCG: " << 100.*pdp.MillisecondsPPCG / duration_cast< milliseconds >( endTime ).count() << "%" << endl;
            
    ///////////////////////////////////////////////// Results
    cout << "finishing ... " << std::flush;

    auto x = totalSpace_.embed(res);
    
    ///////////////////////////////////////////////// VTK-Files
    
    PDE::writeVTK<Descriptions>(x,"pdp");
    
    cout << "done" << endl; 
    
    return 0;
}
