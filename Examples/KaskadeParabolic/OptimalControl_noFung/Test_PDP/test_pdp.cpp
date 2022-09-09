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
#include <Spacy/Notify.h>
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
    
    auto N                    = getParameter( pt, "stepsInTime", 5 ); //cout << "N: "; std::cin >> N;
    const auto T                    = getParameter( pt, "endTime", 2. );
    
    auto initialRefinements   = getParameter( pt, "initialRefinements", 1 );
    const auto iterativeRefinements = getParameter( pt, "iterativeRefinements", 0 ); //aktuell unwichtig
    
    const auto maxSteps             = getParameter( pt, "maxSteps", 5000 );
    const auto FEorder              = getParameter( pt, "FEorder", 1 );
    const auto verbose              = getParameter( pt, "verbose", 1 );
    
    double alpha                    = getParameter( pt, "alpha", 1e-4 );
    double c                        = getParameter( pt, "cPara", 0 );
    double d                        = getParameter( pt, "dPara", 1 );
    double e                        = getParameter( pt, "ePara", 0 );

    const auto stateAcc             = getParameter (pt, "stateAccuracy", 1e-1); //aktuell unwichtig
    const auto adjointAcc           = getParameter (pt, "adjointAccuracy", 1e-1); //aktuell unwichtig
    auto modifiedCGAcc        = getParameter (pt, "modifiedCGAccuracy", 1e-1); //cout << "CG-Acc.: "; std::cin >> modifiedCGAcc;
    const auto chebyshevAcc         = getParameter (pt, "chebyshevAccuracy", 1e-2); //aktuell unwichtig
    auto desiredAccuracy      = getParameter (pt, "desiredAccuracy", 1e-7);
    
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
    int counter = 0;
    bool neu = false;
    std::string dateiName;
    unsigned OrigPDPiterations = 0, OrigPPCGiterations = 0, OrigMGiterations = 0,
     JacPDPiterations = 0, JacPPCGiterations = 0, JacMGiterations = 0,
      JacYPDPiterations = 0, JacYPPCGiterations = 0, JacYMGiterations = 0,
       DiagPDPiterations = 0, DiagPPCGiterations = 0, DiagMGiterations = 0, 
       OrigYPDPiterations = 0, OrigYPPCGiterations = 0, OrigYMGiterations = 0;
//     std::ofstream out ("/home/nuss/98/bt702998/Spacy/Examples/KaskadeParabolic/OptimalControl_noFung/Test_PDP/build/Ref" + std::to_string(initialRefinements) + "alpha" + std::to_string(alpha) + "NTest" + ".txt");
    
//     out << "N" << ", Original, Jacobi, Orig30, Diagonal, Jac30, ErrOrig, ErrJac, ErrOrig30, ErrDiag, ErrJac30" << endl;
    
//     for(alpha = 1e-2; alpha >= 1e-3; alpha *= 0.1)
    {
//     unsigned mod = 0;
//     for(desiredAccuracy = 1e-4; desiredAccuracy>=1e-10; desiredAccuracy*=0.01)
    {
     
//     std::ofstream out ("/home/nuss/98/bt702998/Spacy/Examples/KaskadeParabolic/OptimalControl_noFung/Test_PDP/build/Ref" + std::to_string(initialRefinements) + "N" + std::to_string(N) + "Acc" + std::to_string((int)log10(1./desiredAccuracy)) + "AlphaTest" + ".txt");
    
//     out << "Refinements" << ", Original, MG, Diagonal, MG30, ErrOrig, ErrMG, ErrDiag, ErrMG30, PDPstepsOrig, PDPstepsMG, PDPstepsDiag, PDPstepsMG30, PPCGstepsOrig, PPCGstepsMG, PPCGstepsDiag, PPCGstepsMG30, Orig, MG, Diag, MG30" << endl;
    
    Spacy::Vector res;
    Spacy::Vector resOrig;
    Spacy::Vector resJac;
    Spacy::Vector resJacY;
    Spacy::Vector resOrigY;
    Spacy::Vector resDiag;
    double termOrig;
    double termJac;
    double termJacY;
    double termOrigY;
    double termDiag;
    
    
    
//     int r = 2;
//     for(alpha = 1e-2; alpha >= 1e-3; alpha *= 0.1)
//     for(N = 200; N<=300; N+=100)
    {
        
//     for(desiredAccuracy = 1e-4; desiredAccuracy>=1e-10; desiredAccuracy*=1e-3)
//     for(N = 200; N<=400; N+=200)
//     for(initialRefinements = 3; initialRefinements <= 5; initialRefinements++)
//     for(alpha = 1e-2; alpha >= 1e-3; alpha *= 0.1)
//     for(modifiedCGAcc = 1e-1; modifiedCGAcc>=1e-10; modifiedCGAcc*=0.1)
    {
        
//     if (counter == 0) 
//     {
// //         desiredAccuracy = 1e-7;
// //         initialRefinements = 5;
// //         alpha = 1e-2;
// //         N = 400;
//     }
    
        
    dateiName = 
// "/home/nuss/98/bt702998/Spacy/Examples/KaskadeParabolic/OptimalControl_noFung/Test_PDP/build/R" + std::to_string(initialRefinements) + "_a" + std::to_string((int)log10(1./alpha)) + "_N" + std::to_string(N) + "_Acc" + std::to_string((int)log10(1./desiredAccuracy)) + "_CG" + std::to_string((int)log10(1./modifiedCGAcc))  + "_Test";
    "/home/nuss/98/bt702998/Spacy/Examples/KaskadeParabolic/OptimalControl_noFung/Test_PDP/build/R" + std::to_string(initialRefinements) + "_a" + std::to_string((int)log10(1./alpha)) + "_N" + std::to_string(N) + "_Acc" + std::to_string((int)log10(1./desiredAccuracy)) + "_CG" + std::to_string((int)log10(1./modifiedCGAcc))  + "_STest";
    if (dim==3) dateiName = dateiName + "_3D";
    dateiName = dateiName + ".txt";
        
    std::ofstream out;
    
    if(neu==true)
    {
        out.open(dateiName, std::ofstream::out);
        out << std::setprecision( 4 );
        out << "S" << ", Orig, MG, Diag, MGS, OrigS, Orig, MG, Diag, MGS, OrigS, Orig, MG, Diag, MGS, OrigS, Orig, MG, Diag, MGS, OrigS, Orig, MG, Diag, MGS, OrigS" << endl;
        
        out.close();
    }
    else neu = true;/**/
        
        
//     std::ofstream out ("/home/nuss/98/bt702998/Spacy/Examples/KaskadeParabolic/OptimalControl_noFung/Test_PDP/build/Ref" + std::to_string(initialRefinements) + "alpha" + std::to_string((int)log10(1./alpha)) + "N" + std::to_string(N) + "Acc" + std::to_string((int)log10(1./desiredAccuracy)) + "CGAccTest" + ".txt");
    
//     out << "modifiedCGAcc" << ", Original, MG, Diagonal, MG30, ErrOrig, ErrMG, ErrDiag, ErrMG30, PDPstepsOrig, PDPstepsMG, PDPstepsDiag, PDPstepsMG30, PPCGstepsOrig, PPCGstepsMG, PPCGstepsDiag, PPCGstepsMG30, Orig, MG, Diag, MG30" << endl;
    
//     for(N = 100; N<=1000; N*=1.2311)
    {
        
    Spacy::KaskadeParabolic::GridManager< Spaces > gm( N, ::Spacy::Real{ T }, initialRefinements, dim, FEorder );
    
//     std::ofstream outcsv ("/home/nuss/98/bt702998/Spacy/Examples/KaskadeParabolic/OptimalControl_noFung/Test_PDP/build/Ref" + std::to_string(initialRefinements) + "alpha" + std::to_string(alpha) + "Diagonal" + std::to_string(mod) + "Updated" + std::to_string(N) + ".csv");
//     outcsv << "Y" << ", " << "U" << ", " << "ms" << ", " << "Fehler" << endl;
    
    auto startTime = high_resolution_clock::now();
    auto endTime = high_resolution_clock::now() - startTime;
    auto tstartTime = high_resolution_clock::now();
    auto tendTime = high_resolution_clock::now() - tstartTime;
    auto ttstartTime = high_resolution_clock::now();
    auto ttendTime = high_resolution_clock::now() - tstartTime;
    auto tttstartTime = high_resolution_clock::now();
    auto tttendTime = high_resolution_clock::now() - tstartTime;
    auto ttttstartTime = high_resolution_clock::now();
    auto ttttendTime = high_resolution_clock::now() - tstartTime;
    int inner_counter = 0;
    int Y = N/4; if(Y<3) Y = N;
//     cout << "Y: ";
//     std::cin >> Y;
//     for(Y=10; Y<=N; Y+=10)
    {
        
    Spacy::KaskadeParabolic::GridManager< Spaces > gmy( Y, ::Spacy::Real{ T }, initialRefinements, dim, FEorder/*, "sinus"*/);
    //     gmy.getTempGrid().print();
    
    {
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
//     auto fLu = makeC2Functional<FunctionalDefinition>( FuncGenerator, gmu, totalSpace_u ); fLu.setVerbosity( verbose > 1 ? true : false );
    

    cout << "Assembling" << endl;
    auto x0 = zero(totalSpace_);
    auto x0s = zero(totalSpace_y);
//     auto x0u = zero(totalSpace_u);
    fL.hessian(x0);
    fLs.hessian(x0s);
//     fLu.hessian(x0u);
    auto b = fL.d1(x0);
    auto bs = fLs.d1(x0s);
//     auto bu = fLu.d1(x0u);
    
    auto p_ref = zero(adjointSpace_);
    for(int i = 0; i<N; i++)
    {
        auto& p_refVS = Spacy::cast_ref<Vector<StateDescriptions>>(p_ref).get_nonconst(i);
        ::Kaskade::interpolateGloballyWeak<::Kaskade::PlainAverage >(
            ::boost::fusion::at_c< 0 >( p_refVS.data ),
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
    
    auto diagStateUpdate_    = Inform();
    auto diagAdjointUpdate_  = Inform();
    
    auto As                  = makeCallableBlock<FunctionalDefinition>(fLs,stateSpace_y,adjointSpace_y,true,nT);
    auto ATs                 = makeCallableBlock<FunctionalDefinition>(fLs,adjointSpace_y,stateSpace_y,true,nT);
    auto minusBs             = makeCallableBlock<FunctionalDefinition>(fLs,controlSpace_y,adjointSpace_y,true,nT);
    auto Ms                  = makeSymmetricCallableBlock<FunctionalDefinition>(fLs,primalSpace_y,true,nT);
    auto Mys                 = makeSymmetricCallableBlock<FunctionalDefinition>(fLs,stateSpace_y,true,nT);
    auto Mus                 = makeSymmetricCallableBlock<FunctionalDefinition>(fLs,controlSpace_y,true,nT);
    
    auto Ds                  = makeCallableBlock<FunctionalDefinition>(fLs,stateSpace_y,adjointSpace_y,false,nT);
    auto DTs                 = makeCallableBlock<FunctionalDefinition>(fLs,adjointSpace_y,stateSpace_y,false,nT);

    auto Precs               = makeSingleBlockPreconditioner<DualDescriptions,StateDescriptions,FunctionalDefinition>
                                (adjointSpace_y,stateSpace_y,fLs );
    auto PrecTs              = makeSingleBlockPreconditioner<StateDescriptions,DualDescriptions,FunctionalDefinition>
                                (stateSpace_y,adjointSpace_y,fLs,true );
                               
    auto controlSolver_y     = makeDirectSolver<ControlDescriptions,ControlDescriptions>(Mus,fLs,controlSpace_y, controlSpace_y);
    auto surrStateSolver_y   = makeSurrogatePDESolver<StateDescriptions,DualDescriptions>(As,Precs,gmy);
    auto surrAdjointSolver_y = makeSurrogatePDESolver<DualDescriptions,StateDescriptions>(ATs,PrecTs,gmy,false);
    auto stateSolver_y       = makePDESolver<StateDescriptions,DualDescriptions>(Ds,As,Precs,surrStateSolver_y,fLs,stateAcc);
    auto adjointSolver_y     = makePDESolver<DualDescriptions,StateDescriptions>(DTs,ATs,PrecTs,surrAdjointSolver_y,fLs,adjointAcc,false);
    
//     auto Diags                = makeDiagHandler<DualDescriptions,StateDescriptions>(stateSolver_y(p_refy),p_refy);
//     auto diagStateSolver_y   = makeSurrogateOperator<DualDescriptions,StateDescriptions>
//                                 (stateSolver_y,Precs,Diags,adjointSpace_y,stateSpace_y, diagStateUpdate_);
//     auto diagAdjointSolver_y = makeSurrogateOperator<StateDescriptions,DualDescriptions>
//                                 (adjointSolver_y,PrecTs,Diags,stateSpace_y,adjointSpace_y, diagAdjointUpdate_, false);

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
    auto stateSolver_        = makePDESolver<StateDescriptions,DualDescriptions>(D,A,Prec,surrStateSolver_,fL,stateAcc);
    auto adjointSolver_      = makePDESolver<DualDescriptions,StateDescriptions>(DT,AT,PrecT,surrAdjointSolver_,fL,adjointAcc,false);
    
    auto Diag                = makeDiagHandler<DualDescriptions,StateDescriptions>(stateSolver_(p_ref),p_ref);
//     Diag.printDiag();
    auto DiagT                = makeDiagHandler<StateDescriptions,DualDescriptions>(adjointSolver_(y_ref),y_ref);
//     DiagT.printDiag();
    auto diagStateSolver_    = makeSurrogateOperator<DualDescriptions,StateDescriptions>
                                (stateSolver_,Prec,DiagT,adjointSpace_,stateSpace_, diagStateUpdate_);
    auto diagAdjointSolver_  = makeSurrogateOperator<StateDescriptions,DualDescriptions>
                                (adjointSolver_,PrecT,DiagT,stateSpace_,adjointSpace_, diagAdjointUpdate_, false);
//     auto diagStateSolver_    = makeSurrogateDiagPDESolver<StateDescriptions,DualDescriptions>(A,Prec,gm);
//     auto diagAdjointSolver_  = makeSurrogateDiagPDESolver<DualDescriptions,StateDescriptions>(AT,PrecT,gm,false);
    
//     auto Muu                 = makeSymmetricCallableBlock<FunctionalDefinition>(fLu,controlSpace_u,true,nT);
//     auto controlSolver_u     = makeDirectSolver<ControlDescriptions,ControlDescriptions>(Muu,fLu,controlSpace_u, controlSpace_u);
//     auto mBu                 = makeCallableBlock<FunctionalDefinition>(fLu,controlSpace_u,adjointSpace_u,true,nT);
//     auto minusBu             = makeTimeOperator( mBu, controlSpace_u, adjointSpace_y );
    
    ///////////////////////////////////////////////// Execute
    cout << "Execute" << endl;
    cout << dateiName << endl;
    /*
    cout << endl;
    startTime = high_resolution_clock::now();
    stateSolver_(p_ref);                        
    endTime = high_resolution_clock::now() - startTime;
    cout << "computation time Original: " << duration_cast< milliseconds >( endTime ).count() << "ms, " 
            << duration_cast< microseconds >( endTime ).count() % 1000 << "mics." << endl;
    cout << endl;
    startTime = high_resolution_clock::now();
    surrStateSolver_(p_ref);                        
    endTime = high_resolution_clock::now() - startTime;
    cout << "computation time Jacobi: " << duration_cast< milliseconds >( endTime ).count() << "ms, " 
            << duration_cast< microseconds >( endTime ).count() % 1000 << "mics." << endl;
    cout << endl;
    startTime = high_resolution_clock::now();
    diagStateSolver_(p_ref);                        
    endTime = high_resolution_clock::now() - startTime;
    cout << "computation time Diagonal: " << duration_cast< milliseconds >( endTime ).count() << "ms, " 
            << duration_cast< microseconds >( endTime ).count() % 1000 << "mics." << endl;
    cout << endl;
    /**/
    /*
    if( inner_counter == 0 )
//     if(counter == 0)
    {
        ::Spacy::PDP_withoutA::Solver pdpOrig(M,My,Mu,stateSolver_,adjointSolver_,minusB,minusB,stateSolver_,adjointSolver_,diagStateUpdate_,diagAdjointUpdate_,controlSolver_,controlSolver_,totalSpace_,totalSpace_,N,N,N);
        
        
        pdpOrig.setMaxSteps(maxSteps);
        pdpOrig.setRelativeAccuracy(desiredAccuracy);
        pdpOrig.setInternalAccuracies(stateAcc,adjointAcc,modifiedCGAcc,chebyshevAcc);
        pdpOrig.setVerbosity(false);
        
        ttstartTime = high_resolution_clock::now();
//         resOrig = totalSpace_.embed(pdpOrig( b ));                         // Start PDP
        ttendTime = high_resolution_clock::now() - ttstartTime;
        resOrig = zero(totalSpace_);
        
        OrigPDPiterations = pdpOrig.getIterations();
        OrigPPCGiterations = pdpOrig.getPPCGIterations();
        OrigMGiterations = pdpOrig.getMGIterations();
        termOrig = sqrt(pdpOrig.sumNormSquaredUpdate_);
        
        cout << endl;
        /* 
        for(double gen = 1e-4; gen >= 1e-5; gen *= 0.1)
        {
            pdpOrig.setMaxSteps(maxSteps);
            pdpOrig.setRelativeAccuracy(1e-10);
            pdpOrig.setInternalAccuracies(stateAcc,adjointAcc,gen,chebyshevAcc);
            pdpOrig.setVerbosity(true);
            
            startTime = high_resolution_clock::now();
            resDiag = totalSpace_.embed(pdpOrig( b ));                         // Start PDP
            endTime = high_resolution_clock::now() - ttstartTime;
            
            DiagPDPiterations = pdpOrig.getIterations();
            DiagPPCGiterations = pdpOrig.getPPCGIterations();
            DiagMGiterations = pdpOrig.getMGIterations();
            termDiag = sqrt(pdpOrig.sumNormSquaredUpdate_);
            
            auto vRDDiag = resOrig - resDiag;
            auto errRDDiag = sqrt((M(primalSpace_.project(vRDDiag)))(primalSpace_.project(vRDDiag)));
            auto errMaxDiag = desiredAccuracy * termOrig + gen * termDiag ;
            
            cout << gen << ",\t" << "exact" << ", \t" << "Gen" << endl;
            cout << gen << ",\t" << duration_cast< seconds >( ttendTime ).count()  + 0.001 * (duration_cast< milliseconds >( ttendTime ).count()%1000) << ", \t" << duration_cast< seconds >( endTime ).count()  + 0.001 * (duration_cast< milliseconds >( endTime ).count()%1000) << endl
            << gen << ",\t" << errRDDiag << " < " << errMaxDiag << endl 
            << gen << ",\t" << OrigPDPiterations << ", \t" << DiagPDPiterations << endl
            << gen << ",\t" << OrigPPCGiterations << ", \t" << DiagPPCGiterations << endl
            << gen << ",\t" << OrigMGiterations << ", \t" << DiagMGiterations << endl << endl;
        }*/
/*
        cout << "computation time Original: " << duration_cast< seconds >( ttendTime ).count() << "s, " 
            << duration_cast< milliseconds >( ttendTime ).count() % 1000 << "ms." << endl;
        cout << "davon PPCG: " << 100.*pdpOrig.MillisecondsPPCG / duration_cast< milliseconds >( ttendTime ).count() << "%" << endl;
        
        ::Spacy::PDP_withoutA::Solver pdpJac(M,My,Mu,stateSolver_,adjointSolver_,minusB,minusB,surrStateSolver_,surrAdjointSolver_,diagStateUpdate_,diagAdjointUpdate_,controlSolver_,controlSolver_,totalSpace_,totalSpace_,N,N,N);
        
        pdpJac.setMaxSteps(maxSteps);
        pdpJac.setRelativeAccuracy(desiredAccuracy);
        pdpJac.setInternalAccuracies(stateAcc,adjointAcc,modifiedCGAcc,chebyshevAcc);
        pdpJac.setVerbosity(false);
        
        tttstartTime = high_resolution_clock::now();
//         resJac = totalSpace_.embed(pdpJac( b ));                         // Start PDP
        tttendTime = high_resolution_clock::now() - tttstartTime;
        resJac = zero(totalSpace_);
        
        JacPDPiterations = pdpJac.getIterations();
        JacPPCGiterations = pdpJac.getPPCGIterations();
        JacMGiterations = pdpJac.getMGIterations();
        termJac = sqrt(pdpJac.sumNormSquaredUpdate_);
        
        cout << "computation time Jacobi: " << duration_cast< seconds >( tttendTime ).count() << "s, " 
            << duration_cast< milliseconds >( tttendTime ).count() % 1000 << "ms." << endl;
        cout << "davon PPCG: " << 100.*pdpJac.MillisecondsPPCG / duration_cast< milliseconds >( tttendTime ).count() << "%" << endl;
        
    }
        ::Spacy::PDP_withoutA::Solver pdpJacY(M,Mys,Mus,stateSolver_,adjointSolver_,minusB,minusBs,surrStateSolver_y,surrAdjointSolver_y,diagStateUpdate_,diagAdjointUpdate_,controlSolver_,controlSolver_y,totalSpace_,totalSpace_y,N,Y,Y);
        
        pdpJacY.setMaxSteps(maxSteps);
        pdpJacY.setRelativeAccuracy(desiredAccuracy);
        pdpJacY.setInternalAccuracies(stateAcc,adjointAcc,modifiedCGAcc,chebyshevAcc);
        pdpJacY.setVerbosity(true);
//         pdpJacY.MillisecondsPPCG = duration_cast< milliseconds >( 10*ttendTime ).count();
        
        ttttstartTime = high_resolution_clock::now();
        resJacY = totalSpace_.embed(pdpJacY( b ));                         // Start PDP
        ttttendTime = high_resolution_clock::now() - ttttstartTime;
//         resJacY = zero(totalSpace_);
        
        JacYPDPiterations = pdpJacY.getIterations();
        JacYPPCGiterations = pdpJacY.getPPCGIterations();
        JacYMGiterations = pdpJacY.getMGIterations();
        termJacY = sqrt(pdpJacY.sumNormSquaredUpdate_);
    
        cout << "computation time Jacobi" << Y << ": " << duration_cast< seconds >( ttttendTime ).count() << "s, " 
            << duration_cast< milliseconds >( ttttendTime ).count() % 1000 << "ms." << endl;
        cout << "davon PPCG: " << 100.*pdpJacY.MillisecondsPPCG / duration_cast< milliseconds >( ttttendTime ).count() << "%" << endl;
            
            
        ::Spacy::PDP_withoutA::Solver pdpOrigY(M,Mys,Mus,stateSolver_,adjointSolver_,minusB,minusBs,stateSolver_y,adjointSolver_y,diagStateUpdate_,diagAdjointUpdate_,controlSolver_,controlSolver_y,totalSpace_,totalSpace_y,N,Y,Y);
        
        pdpOrigY.setMaxSteps(maxSteps);
        pdpOrigY.setRelativeAccuracy(desiredAccuracy);
        pdpOrigY.setInternalAccuracies(stateAcc,adjointAcc,modifiedCGAcc,chebyshevAcc);
        pdpOrigY.setVerbosity(false);
//         pdpOrigY.MillisecondsPPCG = duration_cast< milliseconds >( 10*ttendTime ).count();
        
        tstartTime = high_resolution_clock::now();
//         resOrigY = totalSpace_.embed(pdpOrigY( b ));                         // Start PDP
        tendTime = high_resolution_clock::now() - tstartTime;
        resOrigY = zero(totalSpace_);
        
        OrigYPDPiterations = pdpOrigY.getIterations();
        OrigYPPCGiterations = pdpOrigY.getPPCGIterations();
        OrigYMGiterations = pdpOrigY.getMGIterations();
        termOrigY = sqrt(pdpOrigY.sumNormSquaredUpdate_);
    
        cout << "computation time Orig" << Y << ": " << duration_cast< seconds >( tendTime ).count() << "s, " 
            << duration_cast< milliseconds >( tendTime ).count() % 1000 << "ms." << endl;
        cout << "davon PPCG: " << 100.*pdpOrigY.MillisecondsPPCG / duration_cast< milliseconds >( tendTime ).count() << "%" << endl;
    
    if( inner_counter++ == 0 )
//     if(counter == 0)
    {
    ::Spacy::PDP_withoutA::Solver pdpDiag(M,My,Mu,stateSolver_,adjointSolver_,minusB,minusB,diagStateSolver_,diagAdjointSolver_,diagStateUpdate_,diagAdjointUpdate_,controlSolver_,controlSolver_,totalSpace_,totalSpace_,N,N,N);
    
    pdpDiag.setMaxSteps(maxSteps);
    pdpDiag.setRelativeAccuracy(desiredAccuracy);
    pdpDiag.setInternalAccuracies(stateAcc,adjointAcc,modifiedCGAcc,chebyshevAcc);
    pdpDiag.setVerbosity(false);
//     pdpDiag.MillisecondsPPCG = duration_cast< milliseconds >( 10*ttendTime ).count();
    
    startTime = high_resolution_clock::now();
//     if(alpha >= 1e-4)
//         resDiag = totalSpace_.embed(pdpDiag( b ));                         // Start PDP
    endTime = high_resolution_clock::now() - startTime;
//     if (alpha < 1e-4)
        resDiag = zero(totalSpace_);
    
    DiagPDPiterations = pdpDiag.getIterations();
    DiagPPCGiterations = pdpDiag.getPPCGIterations();
    DiagMGiterations = pdpDiag.getMGIterations();
    termDiag = sqrt(pdpDiag.sumNormSquaredUpdate_);
    
    cout << "computation time Diagonal: " << duration_cast< seconds >( endTime ).count() << "s, " 
            << duration_cast< milliseconds >( endTime ).count() % 1000 << "ms." << endl;
    cout << "davon PPCG: " << (double)pdpDiag.MillisecondsPPCG / duration_cast< milliseconds >( endTime ).count() << "%" << endl;
    }
    /**/
    
    
    auto ones = zero( stateSpace_ );
    auto co = zero(controlSpace_);
    auto& ov = Spacy::cast_ref< Vector< StateDescriptions > >( ones );
    auto& cv = Spacy::cast_ref< Vector< ControlDescriptions > >( co );
    
//     for(int t = 0; t < N; t++)
//     {
//         auto& ot = ov.getCoeffVec_nonconst(t);
//         for(int i=0; i<ot.size(); ++i)
//         {
//             ot[i] = Dune::FieldVector< double, 1 >( i*t );
//         }
//     }
//     for(int t = 0; t < N; t++)
//     {
//         auto& ct = cv.getCoeffVec_nonconst(t);
//         for(int i=0; i<ct.size(); ++i)
//         {
//             ct[i] = Dune::FieldVector< double, 1 >( i*t );
//         }
//     }
    
    
    
    auto& ct = cv.getCoeffVec_nonconst(0);
    auto& ot = ov.getCoeffVec_nonconst(0);
    for(int j=0; j<ct.size(); ++j)
    {
        for(int i=0; i<ct.size(); ++i)
        {
            if(i==j)
            {
                ct[i] = Dune::FieldVector< double, 1 >( 1 );
                ot[i] = Dune::FieldVector< double, 1 >( 1 );
            }
            else
            {
                ct[i] = Dune::FieldVector< double, 1 >( 0 );
                ot[i] = Dune::FieldVector< double, 1 >( 0 );
            }
        }
        
        auto s = - alpha * minusB(co);
        auto v = Mu(co);
        
        auto& sv = Spacy::cast_ref< Vector< ControlDescriptions > >( s );
        auto& vv = Spacy::cast_ref< Vector< ControlDescriptions > >( v );
        
        std::cout << "B(c)[" << j << "]: ( ";
        auto& st = sv.getCoeffVec_nonconst(0);
        for(int i=0; i<st.size(); ++i)
        {
            cout << st[i] << " ";
        }
        cout << ")" << endl << endl;
        
        std::cout << "Mu(c)[" << j << "]: ( ";
        auto& vt = vv.getCoeffVec_nonconst(0);
        for(int i=0; i<vt.size(); ++i)
        {
            cout << vt[i] << " ";
        }
        cout << ")" << endl << endl << endl;
    }
    /**/
    
/*
    auto vRDJac = resOrig - resJac;
    auto errRDJac = sqrt((M(primalSpace_.project(vRDJac)))(primalSpace_.project(vRDJac)));
    auto errMaxJac = desiredAccuracy * ( termOrig + termJac );
    
    auto vRDDiag = resOrig - resDiag;
    auto errRDDiag = sqrt((M(primalSpace_.project(vRDDiag)))(primalSpace_.project(vRDDiag)));
    auto errMaxDiag = desiredAccuracy * ( termOrig + termDiag );
    
    auto vRDJacY = resOrig - resJacY;
    auto errRDJacY = sqrt((M(primalSpace_.project(vRDJacY)))(primalSpace_.project(vRDJacY)));
    auto errMaxJacY = desiredAccuracy * ( termOrig + termJacY );
    
    auto vRDOrigY = resOrig - resOrigY;
    auto errRDOrigY = sqrt((M(primalSpace_.project(vRDOrigY)))(primalSpace_.project(vRDOrigY)));
    auto errMaxOrigY = desiredAccuracy * ( termOrig + termOrigY );
    
    
    double durchlaeufe = 1;
    double prz = 100.*(++counter)/(durchlaeufe);
    
    
    out.open(dateiName, std::ofstream::out | std::ofstream::app);
    
    out << Y << ", " << duration_cast< seconds >( ttendTime ).count()  + 0.001 * (duration_cast< milliseconds >( ttendTime ).count()%1000) << ", " << duration_cast< seconds >( tttendTime ).count() + 0.001 * (duration_cast< milliseconds >( tttendTime ).count()%1000) << ", " << duration_cast< seconds >( endTime ).count() + 0.001 * (duration_cast< milliseconds >( endTime ).count()%1000) << ", " << duration_cast< seconds >( ttttendTime ).count() + 0.001 * (duration_cast< milliseconds >( ttttendTime ).count()%1000) << ", " << duration_cast< seconds >( tendTime ).count()  + 0.001 * (duration_cast< milliseconds >( tendTime ).count()%1000) << ", 0, " << errRDJac << ", " << errRDDiag << ", " << errRDJacY << ", " << errRDOrigY << ", " << OrigPDPiterations << ", " << JacPDPiterations << ", " << DiagPDPiterations << ", " << JacYPDPiterations << ", " << OrigYPDPiterations << ", " << OrigPPCGiterations << ", " << JacPPCGiterations << ", " << DiagPPCGiterations << ", " << JacYPPCGiterations << ", " << OrigYPPCGiterations << ", " << OrigMGiterations << ", " << JacMGiterations << ", " << DiagMGiterations << ", " << JacYMGiterations << ", " << OrigYMGiterations << endl;
    
    out.close();
    

    cout << "Y,\t\tArt,\t\tOrig,\t\tMG,\t\tDiag,\t\tMGS,\t\tOrigS" << endl;
    cout << Y << ",\t\tSek,\t\t" << duration_cast< seconds >( ttendTime ).count()  + 0.001 * (duration_cast< milliseconds >( ttendTime ).count()%1000) << ",\t\t" << duration_cast< seconds >( tttendTime ).count()  + 0.001 * (duration_cast< milliseconds >( tttendTime ).count()%1000) << ",\t\t" << duration_cast< seconds >( endTime ).count()  + 0.001 * (duration_cast< milliseconds >( endTime ).count()%1000) << ",\t\t" << duration_cast< seconds >( ttttendTime ).count()  + 0.001 * (duration_cast< milliseconds >( ttttendTime ).count()%1000) << ",\t\t" << duration_cast< seconds >( tendTime ).count() + 0.001 * (duration_cast< milliseconds >( tendTime ).count()%1000) << endl 
    << Y << ",\t\tErr,\t\t" << "0,\t\t" << errRDJac << ",\t\t" << errRDDiag << ",\t\t" << errRDJacY << ",\t\t" << errRDOrigY << endl 
    << Y << ",\t\tErrMax,\t\t" << "0,\t\t" << errMaxJac << ",\t\t" << errMaxDiag << ",\t\t" << errMaxJacY << ",\t\t" << errMaxOrigY << endl 
    << Y << ",\t\tPDP,\t\t" << OrigPDPiterations << ",\t\t" << JacPDPiterations << ",\t\t" << DiagPDPiterations << ",\t\t" << JacYPDPiterations << ",\t\t" << OrigYPDPiterations << endl
    << Y << ",\t\tPPCG,\t\t" << OrigPPCGiterations << ",\t\t" << JacPPCGiterations << ",\t\t" << DiagPPCGiterations << ",\t\t" << JacYPPCGiterations << ",\t\t" << OrigYPPCGiterations << endl
    << Y << ",\t\tMG,\t\t" << OrigMGiterations << ",\t\t" << JacMGiterations << ",\t\t" << DiagMGiterations << ",\t\t" << JacYMGiterations << ",\t\t" << OrigYMGiterations << endl;
    /**/
    
//     cout << "\t\t\t\t\t\t\t\t Time to beat: " << duration_cast< seconds >( tendTime ).count() << "s, " << duration_cast< milliseconds >( tendTime ).count() % 1000 << "ms." << endl;
//     cout << "\t\t\t\t\t\t\t\t minimal Time left: " << duration_cast< hours >( (durchlaeufe-counter)*endTime ).count() << "h," << duration_cast< minutes >( (durchlaeufe-counter)*endTime ).count() % 60 << "min, "<< duration_cast< seconds >( (durchlaeufe-counter)*endTime ).count() % 60 << "s, " << duration_cast< milliseconds >( (durchlaeufe-counter)*endTime ).count() % 1000 << "ms." << endl;
    
    ///////////////////////////////////////////////// Results
    cout << "finishing ... " << std::flush;

    auto x = totalSpace_.embed(resJacY);
    
    ///////////////////////////////////////////////// VTK-Files
    
    PDE::writeVTK<Descriptions>(x,"testVier");
    
    cout << "done" << endl; 
    /**/
    
    } } } } } } }
    /**/
    return 0;
}
