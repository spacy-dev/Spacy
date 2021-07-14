#define FUSION_MAX_VECTOR_SIZE 15 // NOLINT(cppcoreguidelines-macro-usage)

 
#include <chrono>
#include <iostream>
 
#include <Spacy/Spacy.h>
#include <Spacy/VectorSpace.h>
#include <Spacy/Adapter/kaskadeParabolic.hh>

#include <Spacy/Algorithm/PrimalDualProjection/PDP_withoutA.hh>
#include <Spacy/Algorithm/CG/RegularizeViaPreconditioner.h>
#include <Spacy/Algorithm/CG/TerminationCriterion.h>
#include <Spacy/Algorithm/CG/TerminationCriteria.h>

#include <io/vtk.hh>
#include <utilities/kaskopt.hh>

#include <dune/grid/config.h>
#include <dune/grid/uggrid.hh>

#include "../optimal_control_withoutfung.hh"

#include "utilities/memory.hh"
#include "mg/additiveMultigrid.hh"

using std::cout;
using std::endl;
using namespace Kaskade;
using namespace Spacy::KaskadeParabolic;

constexpr int dim = 2;
constexpr int cp = 1;

int main( int argc, char* argv[] )
{
    cout << endl;
    const auto silence = 0;
    std::unique_ptr< boost::property_tree::ptree > pt = getKaskadeOptions( argc, argv, silence, false );

    const auto maxSteps             = getParameter( pt, "maxSteps", 1000 );
    const auto FEorder              = getParameter( pt, "FEorder", 1 );
    const auto verbose              = getParameter( pt, "verbose", 1 );
    
    const auto N                    = getParameter( pt, "stepsInTime", 11 );
    const auto T                    = getParameter( pt, "endTime", 1. );
    
    const auto initialRefinements   = getParameter( pt, "initialRefinements", 3 );
    const auto iterativeRefinements = getParameter( pt, "iterativeRefinements", 0 );
    
    double alpha                    = getParameter( pt, "alpha", 0.1 );
    double c                        = getParameter( pt, "cPara", 0 );
    double d                        = getParameter( pt, "dPara", 1 );
    double e                        = getParameter( pt, "ePara", 0 );
    
    const auto numberOfCubesX       = getParameter (pt, "NumberOfCubesX", 1);
    const auto numberOfCubesY       = getParameter (pt, "NumberOfCubesY", 1);
    const auto numberOfCubesZ       = getParameter (pt, "NumberOfCubesZ", 1);
    
    const auto stateAcc             = getParameter (pt, "stateAccuracy", 1e-3);
    const auto adjointAcc           = getParameter (pt, "adjointAccuracy", 1e-3);
    const auto modifiedCGAcc        = getParameter (pt, "modifiedCGAccuracy", 1e-3);
    const auto chebyshevAcc         = getParameter (pt, "chebyshevAccuracy", 1e-3);
    const auto desiredAccuracy      = getParameter (pt, "desiredAccuracy", 1e-4);
    
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
    
    auto totalSpace_ = makeHilbertSpace< VariableDescriptions, Spaces >( gm, "X" );
    
    auto primalSpace_ = makeHilbertSpace< PrimalVariables, Spaces >( gm, "Z" );
    using ZXIdxMap = boost::fusion::vector< IdxPair<0,0>, IdxPair<1,1> >;
    setSubSpaceRelation< PrimalDescriptions, Descriptions, ZXIdxMap >( primalSpace_, totalSpace_ );
    
    auto adjointSpace_ = makeHilbertSpace< DualVariables, Spaces >( gm, "P" );
    using PXIdxMap = boost::fusion::vector< IdxPair<0,2> >;
    setSubSpaceRelation< DualDescriptions, Descriptions, PXIdxMap >( adjointSpace_, totalSpace_ );
    
    auto stateSpace_ = makeHilbertSpace< StateVariables, Spaces >( gm, "Y" );
    using YXIdxMap = boost::fusion::vector< IdxPair<0,0> >;
    setSubSpaceRelation< StateDescriptions, Descriptions, YXIdxMap >( stateSpace_, totalSpace_ );
    using YZIdxMap = boost::fusion::vector< IdxPair<0,0> >;
    setSubSpaceRelation< StateDescriptions, PrimalDescriptions, YZIdxMap >( stateSpace_, primalSpace_ );
    
    auto controlSpace_ = makeHilbertSpace< ControlVariables, Spaces >( gm, "U" );
    using UXIdxMap = boost::fusion::vector< IdxPair<0,1> >;
    setSubSpaceRelation< ControlDescriptions, Descriptions, UXIdxMap >( controlSpace_, totalSpace_ );
    using UZIdxMap = boost::fusion::vector< IdxPair<0,1> >;
    setSubSpaceRelation< ControlDescriptions, PrimalDescriptions, UZIdxMap >( controlSpace_, primalSpace_ );
    
    ///////////////////////////////////////////////// Functional 
    auto resultRef = zero( stateSpace_ );
    for (int t = 0; t < N; t++)
    {
        auto& resultRefK = Spacy::cast_ref< Vector< StateDescriptions > >( resultRef );

        const Dune::FieldVector< double, 1 > reference_deformationVector{1};
        auto& refImpl = boost::fusion::at_c< 0 >( resultRefK.get_nonconst(t).data ).coefficients();
        for ( int i = 0; i < refImpl.size(); ++i )
            refImpl[ i ] = reference_deformationVector;
    }
    
    auto resultRefP = zero( adjointSpace_ );
    for (int t = 0; t < N; t++)
    {
        auto& resultRefKP = Spacy::cast_ref< Vector< DualDescriptions > >( resultRefP );

        const Dune::FieldVector< double, 1 > reference_deformationVector{1};
        auto& refImpl = boost::fusion::at_c< 0 >( resultRefKP.get_nonconst(t).data ).coefficients();
        for ( int i = 0; i < refImpl.size(); ++i )
            refImpl[ i ] = reference_deformationVector;
    }
    
    
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
    
    auto D                  = makeMassMatrixBlock<FunctionalDefinition>(fL,stateSpace_,adjointSpace_);
    auto DT                 = makeMassMatrixBlock<FunctionalDefinition>(fL,adjointSpace_,stateSpace_);

    auto BPX                = makeSingleBlockPreconditioner<DualDescriptions,StateDescriptions,FunctionalDefinition>
                                (adjointSpace_,stateSpace_,fL );
    auto BPXT               = makeSingleBlockPreconditioner<StateDescriptions,DualDescriptions,FunctionalDefinition>
                                (stateSpace_,adjointSpace_,fL );
         
    auto controlSolver_     = makeDirectSolver<ControlDescriptions,ControlDescriptions>(Mu,gm,controlSpace_, controlSpace_);
    auto stateSolver_       = makePDESolver<StateDescriptions,DualDescriptions>(D,A,BPX,gm.getTempGrid().getDtVec());
    auto adjointSolver_     = makePDESolver<DualDescriptions,StateDescriptions>(DT,AT,BPXT,gm.getTempGrid().getDtVec(),false);
    auto surrStateSolver_   = makeSurrogatePDESolver<StateDescriptions,DualDescriptions>(D,BPX,A,gm.getTempGrid().getDtVec());
    auto surrAdjointSolver_ = makeSurrogatePDESolver<DualDescriptions,StateDescriptions>(DT,BPXT,AT,gm.getTempGrid().getDtVec(),false);
    
    ///////////////////////////////////////////////// Execute
    cout << "Execute" << endl;

    ::Spacy::PDP_withoutA::Solver pdp(M,stateSolver_,adjointSolver_,minusB,surrStateSolver_,surrAdjointSolver_,controlSolver_,totalSpace_);
    
    pdp.setMaxSteps(maxSteps);
    pdp.setRelativeAccuracy(desiredAccuracy);
    pdp.setInternalAccuracies(stateAcc,adjointAcc,modifiedCGAcc,chebyshevAcc);
    
    cout << "StateSolver:" << endl;
    stateSolver_(resultRefP);
    cout << "AdjointSolver:" << endl;
    adjointSolver_(resultRef);
    cout << "StateSolver(Surr):" << endl;
    surrStateSolver_(resultRefP);
    cout << "AdjointSolver(Surr):" << endl;
    surrAdjointSolver_(resultRef);
    cout << "done" << endl;
    
    auto res = totalSpace_.embed( pdp(b) );
    
    ///////////////////////////////////////////////// Results
    cout << "finishing ... " << std::flush;

    auto x = totalSpace_.embed(res);
    
    ///////////////////////////////////////////////// VTK-Files
    
    PDE::writeVTK<Descriptions>(x,"test");
    
    cout << "done" << endl;
    
    /**/
    return 0;
}
