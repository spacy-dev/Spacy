 #define FUSION_MAX_VECTOR_SIZE 15 // NOLINT(cppcoreguidelines-macro-usage)

#include "nonlinear_elasticity.hh"
// #include "Examples/Kaskade/createSimpleGeometries.hh"

// #include <functional>
// #include <fung/examples/nonlinear_heat.hh>
// #include <fung/fung.hh>
// #include <fung/examples/volumetric_penalty_functions.hh>
// #include <fung/examples/rubber/mooney_rivlin.hh>

#include <Spacy/Spacy.h>
#include <Spacy/Algorithm/PrimalDualProjection/PDP.hh>
#include "Spacy/Algorithm/CG/RegularizeViaPreconditioner.h"
#include "Spacy/Algorithm/CG/TerminationCriterion.h"
#include "Spacy/Algorithm/CG/TerminationCriteria.h"
#include <Spacy/VectorSpace.h>

// #include <Spacy/Adapter/kaskade.hh>
#include <Spacy/Adapter/kaskadeParabolic.hh>
#include <Spacy/Adapter/KaskadeParabolic/Constraints/lqOcpHeader.hh>

// #include <fem/forEach.hh>
#include <fem/gridmanager.hh>
#include <fem/lagrangespace.hh>
// #include <fem/boundaryspace.hh>
// #include "linalg/triplet.hh"
#include <io/vtk.hh>
#include <utilities/gridGeneration.hh>
#include <utilities/kaskopt.hh>

#include <dune/grid/config.h>
#include <dune/grid/uggrid.hh>

#include <chrono>
#include <iostream>

using std::cout;
using std::endl;
using namespace Kaskade;
// using namespace Spacy::Kaskade;

constexpr int dim = 2;
constexpr int cp = 1;

int main( int argc, char* argv[] )
{
    cout << endl;
    const auto silence = 0;
    std::unique_ptr< boost::property_tree::ptree > pt = getKaskadeOptions( argc, argv, silence, false );

    const auto maxSteps             = getParameter( pt, "maxSteps", 1000 );
    const auto FEorder              = getParameter( pt, "FEorder", 1 );
    const auto verbose              = getParameter( pt, "verbose", 2 );
    
    const auto N                    = getParameter( pt, "stepsInTime", 11 );
    const auto T                    = getParameter( pt, "endTime", 30. );
    
    const auto initialRefinements   = getParameter( pt, "initialRefinements", 3 );
    const auto iterativeRefinements = getParameter( pt, "iterativeRefinements", 0 );
    
//     const auto tychonovParameter    = getParameter( pt, "tychonovParameter", 0.001 );
    double alpha                    = getParameter( pt, "alpha", 0.5 );
    double d                        = getParameter( pt, "dPara", 1e-1 );
    double mu                       = getParameter( pt, "muPara", 0 );
    
    const auto numberOfCubesX       = getParameter (pt, "NumberOfCubesX", 1);
    const auto numberOfCubesY       = getParameter (pt, "NumberOfCubesY", 1);
    const auto numberOfCubesZ       = getParameter (pt, "NumberOfCubesZ", 1);
    
    const auto stateAcc             = getParameter (pt, "stateAccuracy", 1e-3);
    const auto adjointAcc           = getParameter (pt, "adjointAccuracy", 1e-3);
    const auto modifiedCGAcc        = getParameter (pt, "modifiedCGAccuracy", 1e-3);
    const auto chebyshevAcc         = getParameter (pt, "chebyshevAccuracy", 1e-3);
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
    ::Spacy::KaskadeParabolic::GridManager< Spaces > gm( N, ::Spacy::Real{ T }, initialRefinements, FEorder );
    
    ///////////////////////////////////////////////// Spacy Vector-Spaces
    cout << "create Vector-Spaces" << endl;
    
    auto totalSpace_ = Spacy::KaskadeParabolic::makeHilbertSpace< VariableDescriptions, Spaces >( gm, "X" );
    
    auto primalSpace_ = Spacy::KaskadeParabolic::makeHilbertSpace< PrimalVariables, Spaces >( gm, "Z" );
    using ZXIdxMap = boost::fusion::vector< Spacy::KaskadeParabolic::IdxPair<0,0>, Spacy::KaskadeParabolic::IdxPair<1,1> >;
    Spacy::KaskadeParabolic::setSubSpaceRelation< PrimalDescriptions, Descriptions, ZXIdxMap >( primalSpace_, totalSpace_ );
    
    auto adjointSpace_ = Spacy::KaskadeParabolic::makeHilbertSpace< DualVariables, Spaces >( gm, "P" );
    using PXIdxMap = boost::fusion::vector< Spacy::KaskadeParabolic::IdxPair<0,2> >;
    Spacy::KaskadeParabolic::setSubSpaceRelation< DualDescriptions, Descriptions, PXIdxMap >( adjointSpace_, totalSpace_ );
    
    auto stateSpace_ = Spacy::KaskadeParabolic::makeHilbertSpace< StateVariables, Spaces >( gm, "Y" );
    using YXIdxMap = boost::fusion::vector< Spacy::KaskadeParabolic::IdxPair<0,0> >;
    Spacy::KaskadeParabolic::setSubSpaceRelation< StateDescriptions, Descriptions, YXIdxMap >( stateSpace_, totalSpace_ );
    using YZIdxMap = boost::fusion::vector< Spacy::KaskadeParabolic::IdxPair<0,0> >;
    Spacy::KaskadeParabolic::setSubSpaceRelation< StateDescriptions, PrimalDescriptions, YZIdxMap >( stateSpace_, primalSpace_ );
    
    auto controlSpace_ = Spacy::KaskadeParabolic::makeHilbertSpace< ControlVariables, Spaces >( gm, "U" );
    using UXIdxMap = boost::fusion::vector< Spacy::KaskadeParabolic::IdxPair<0,1> >;
    Spacy::KaskadeParabolic::setSubSpaceRelation< ControlDescriptions, Descriptions, UXIdxMap >( controlSpace_, totalSpace_ );
    using UZIdxMap = boost::fusion::vector< Spacy::KaskadeParabolic::IdxPair<0,1> >;
    Spacy::KaskadeParabolic::setSubSpaceRelation< ControlDescriptions, PrimalDescriptions, UZIdxMap >( controlSpace_, primalSpace_ );
    
    
    ///////////////////////////////////////////////// Functional 
    cout << "create Functional" << endl;
    auto resultRef = zero(stateSpace_);
    auto& resultRefK = Spacy::cast_ref<Spacy::KaskadeParabolic::Vector<StateDescriptions>>(resultRef);
    
    for (int t = 0; t < N; t++)
    {
        const Dune::FieldVector<double,cp> reference_deformationVector {1};
        auto& refImpl = resultRefK.getCoeffVec_nonconst(t);
        for(int i=0; i<refImpl.size(); ++i)
            refImpl[i]=reference_deformationVector;
    }
    
    

    using FunctionalDefinition = NormalStepFunctional< 0, 1, 2, double, Descriptions >;
    std::function< FunctionalDefinition( const typename Descriptions::VariableSet ) > FuncGenerator =
        [ &d, &mu, &alpha ]( const typename Descriptions::VariableSet y_ref ) { return FunctionalDefinition( alpha, y_ref, d, mu ); };

    auto fL = Spacy::KaskadeParabolic::makeC2Functional<FunctionalDefinition>( FuncGenerator, gm, totalSpace_ );
//     auto fL = Spacy::Kaskade::makeC2Functional(Functional(tychonovParameter,reference_deformation, integrand), totalSpace_);
    /*
    cout << "Assembling" << endl;
    auto x0 = zero(totalSpace_);
    fL.hessian(x0);
    auto b = fL.d1(x0);
    /*
    ///////////////////////////////////////////////// Blocks, Solvers
    cout << "create Blocks and Solvers" << endl;
    
    using KMAT = ::Kaskade::MatrixAsTriplet<double>;
    
    auto A              = makeCallableBlock<Functional>(fL,stateSpace_,adjointSpace_);
    auto AT             = makeCallableBlock<Functional>(fL,adjointSpace_,stateSpace_);
    auto minusB         = makeCallableBlock<Functional>(fL,controlSpace_,adjointSpace_);
    auto M              = makeSymmetricCallableBlock<Functional>(fL,primalSpace_);
    
    auto BPX            = makeSingleBlockPreconditioner<DualDescriptions,StateDescriptions>
                            (adjointSpace_,stateSpace_,moveUnique(makeBPX(fL.assembler().get<2,0>(), gm)));
    auto BPXT           = makeSingleBlockPreconditioner<StateDescriptions,DualDescriptions>
                            (stateSpace_,adjointSpace_,moveUnique(makeBPX(fL.assembler().get<2,0>(), gm)));
    
    auto controlSolver_ = makeDirectSolver<ControlDescriptions,ControlDescriptions>
                            (makeSymmetricMatrixBlock<decltype(fL),KMAT,ControlDescriptions>(fL,controlSpace_).get(),controlSpace_, controlSpace_);
    
    ///////////////////////////////////////////////// Execute
    cout << "Execute" << endl;
    
    Spacy::CG::NoRegularization NR;
    ::Spacy::PDP::Solver pdp(M,A,AT,minusB,BPX,BPXT,controlSolver_,NR,totalSpace_);
    
    pdp.setMaxSteps(maxSteps);
    pdp.setRelativeAccuracy(desiredAccuracy);
    pdp.setInternalAccuracies(stateAcc,adjointAcc,modifiedCGAcc,chebyshevAcc);
    
    auto res = totalSpace_.embed( pdp(b) );
    /**/
    ///////////////////////////////////////////////// Results
    cout << "finishing ... " << std::flush;

    auto x = totalSpace_.embed(resultRef);
    auto y = primalSpace_.project(x);
    auto z = stateSpace_.project(y);
    
    ///////////////////////////////////////////////// VTK-Files
    
    Spacy::KaskadeParabolic::PDE::writeVTK<StateDescriptions>(z,"test");
    
    cout << "done" << endl;
    
    /**/
    return 0;
}
