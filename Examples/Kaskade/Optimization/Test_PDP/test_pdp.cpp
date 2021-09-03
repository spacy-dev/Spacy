 #define FUSION_MAX_VECTOR_SIZE 15 // NOLINT(cppcoreguidelines-macro-usage)

#include "nonlinear_elasticity.hh"
#include "../../interpolate.h"
#include  "../../createSimpleGeometries.hh"

#include <functional>
#include <fung/examples/nonlinear_heat.hh>
#include <fung/fung.hh>
#include <fung/examples/volumetric_penalty_functions.hh>
#include <fung/examples/rubber/mooney_rivlin.hh>

#include <Spacy/Adapter/kaskade.hh>
#include <Spacy/Adapter/KaskadeParabolic/directSolver.hh>
#include <Spacy/Spacy.h>
#include <Spacy/Algorithm/PrimalDualProjection/PDP.hh>
#include "Spacy/Algorithm/CG/RegularizeViaPreconditioner.h"
#include "Spacy/Algorithm/CG/TerminationCriterion.h"
#include "Spacy/Algorithm/CG/TerminationCriteria.h"
#include <Spacy/VectorSpace.h>

#include <fem/forEach.hh>
#include <fem/gridmanager.hh>
#include <fem/lagrangespace.hh>
#include <fem/boundaryspace.hh>
#include "linalg/triplet.hh"
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
using namespace Spacy::Kaskade;

constexpr int dim = 3;

int main( int argc, char* argv[] )
{
    cout << endl;
    const auto silence = 0;
    std::unique_ptr< boost::property_tree::ptree > pt = getKaskadeOptions( argc, argv, silence, false );

    const auto maxSteps             = getParameter( pt, "maxSteps", 1000 );
    const auto FEorder              = getParameter( pt, "FEorder", 1 );
    const auto verbose              = getParameter( pt, "verbose", 2 );
    
    const auto initialRefinements   = getParameter( pt, "initialRefinements", 3 );
    const auto iterativeRefinements = getParameter( pt, "iterativeRefinements", 0 );
    
    const auto tychonovParameter    = getParameter( pt, "tychonovParameter", 0.1 );
    
    const auto numberOfCubesX       = getParameter (pt, "NumberOfCubesX", 1);
    const auto numberOfCubesY       = getParameter (pt, "NumberOfCubesY", 1);
    const auto numberOfCubesZ       = getParameter (pt, "NumberOfCubesZ", 1);
    
    const auto stateAcc             = getParameter (pt, "stateAccuracy", 1e-2);
    const auto adjointAcc           = getParameter (pt, "adjointAccuracy", 1e-2);
    const auto modifiedCGAcc        = getParameter (pt, "modifiedCGAccuracy", 1e-2);
    const auto chebyshevAcc         = getParameter (pt, "chebyshevAccuracy", 1e-2);
    const auto desiredAccuracy      = getParameter (pt, "desiredAccuracy", 1e-2);

    ///////////////////////////////////////////////// Grid
    using Grid = Dune::UGGrid< dim >;
    using GridView = Grid::LeafGridView;
    
    double factor = 1.0/std::max(numberOfCubesX,std::max(numberOfCubesY,numberOfCubesZ));
    auto grid = createPlate<Grid>(factor,numberOfCubesX,numberOfCubesY,numberOfCubesZ);
    GridManager<Grid> gm(std::move(grid),true);
    
    gm.globalRefine( initialRefinements );
    cout << "vertices: " << gm.grid().size( dim ) << endl;
    cout << std::setprecision( 10 );
    
    ///////////////////////////////////////////////// Support of control
    std::vector<int> totalIndexSet(0);
    std::vector<int> partialIndexSetId = { 1 };
    {
        const auto & gv = gm.grid().leafGridView();
        
        std::size_t numberOfBoundaryFaces = 0;
        
        forEachBoundaryFace(gv, [&numberOfBoundaryFaces](const auto & face)
        {	numberOfBoundaryFaces++;});
        totalIndexSet.resize(numberOfBoundaryFaces, 0);
        
        // write more compact with lambda
        const Dune::FieldVector<double,3> up
        { 0., 0., 1.0 };
        
        forEachBoundaryFace(gv, [&up, &totalIndexSet] (const auto & face)
        {
            if(face.centerUnitOuterNormal()*up > 0.5)
            {
                // const auto boundarySegmentId = face.boundarySegmentIndex();
                totalIndexSet.at(face.boundarySegmentIndex()) = 1;
            }
        });
    }
    
    ///////////////////////////////////////////////// Kaskade FEM-Spaces
    cout << "create FEM-Spaces" << endl;
    using L2Space = FEFunctionSpace<BoundaryMapper<ContinuousLagrangeMapper, double, GridView> >;
    using H1Space = FEFunctionSpace< ContinuousLagrangeMapper< double, GridView > >;
    using Spaces = boost::fusion::vector< FEFunctionSpace< ContinuousLagrangeMapper< double, GridView > > const*, FEFunctionSpace<BoundaryMapper<ContinuousLagrangeMapper, double, GridView> > const* >;

    using components = Components<3>;
    using VariableDescriptions = boost::fusion::vector< Variable< SpaceIndex< 0 >, components, VariableId< 0 > >,
                                                        Variable< SpaceIndex< 1 >, components, VariableId< 1 > >,
                                                        Variable< SpaceIndex< 0 >, components, VariableId< 2 > > >;
    using PrimalVariables = boost::fusion::vector< Variable< SpaceIndex< 0 >, components, VariableId< 0 > >,
                                                Variable< SpaceIndex< 1 >, components, VariableId< 1 > > >;
    using StateVariables = boost::fusion::vector< Variable< SpaceIndex< 0 >, components, VariableId< 0 > > >;
    using ControlVariables = boost::fusion::vector< Variable< SpaceIndex< 1 >, components, VariableId< 0 > > >;
    using DualVariables = boost::fusion::vector< Variable< SpaceIndex< 0 >, components, VariableId< 0 > > >;

    using Descriptions = VariableSetDescription< Spaces, VariableDescriptions >;
    using PrimalDescriptions = VariableSetDescription< Spaces, PrimalVariables >;
    using DualDescriptions = VariableSetDescription< Spaces, DualVariables >;
    using StateDescriptions = VariableSetDescription< Spaces, StateVariables >;
    using ControlDescriptions = VariableSetDescription< Spaces, ControlVariables >;
    using VarSet = Descriptions::VariableSet;
    
    L2Space l2Space(gm, gm.grid().leafGridView(), FEorder, totalIndexSet, partialIndexSetId);
    H1Space h1Space( gm, gm.grid().leafGridView(), FEorder );
    Spaces spaces( &h1Space, &l2Space);

    Descriptions desc( spaces, { "y", "u", "p" } );               // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay)
    PrimalDescriptions primalDescription( spaces, { "y", "u" } ); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay)
    DualDescriptions dualDescription( spaces, { "p" } );          // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay)
    StateDescriptions stateDescription( spaces, { "y" } );          // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay)
    ControlDescriptions controlDescription( spaces, { "u" } );          // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay)
    
    ///////////////////////////////////////////////// Spacy Vector-Spaces
    cout << "create Vector-Spaces" << endl;
    
    auto totalSpace_ = makeHilbertSpace< Descriptions >( desc, "X" );
    
    auto primalSpace_ = makeHilbertSpace< PrimalDescriptions >( primalDescription, "Z" );
    using ZXIdxMap = boost::fusion::vector< IdxPair<0,0>, IdxPair<1,1> >;
    setSubSpaceRelation< PrimalDescriptions, Descriptions, ZXIdxMap >( primalSpace_, totalSpace_ );
    
    auto adjointSpace_ = makeHilbertSpace< DualDescriptions >( dualDescription, "P" );
    using PXIdxMap = boost::fusion::vector< IdxPair<0,2> >;
    setSubSpaceRelation< DualDescriptions, Descriptions, PXIdxMap >( adjointSpace_, totalSpace_ );
    
    auto stateSpace_ = makeHilbertSpace< StateDescriptions >( stateDescription, "Y" );
    using YXIdxMap = boost::fusion::vector< IdxPair<0,0> >;
    setSubSpaceRelation< StateDescriptions, Descriptions, YXIdxMap >( stateSpace_, totalSpace_ );
    using YZIdxMap = boost::fusion::vector< IdxPair<0,0> >;
    setSubSpaceRelation< StateDescriptions, PrimalDescriptions, YZIdxMap >( stateSpace_, primalSpace_ );
    
    auto controlSpace_ = makeHilbertSpace< ControlDescriptions >( controlDescription, "U" );
    using UXIdxMap = boost::fusion::vector< IdxPair<0,1> >;
    setSubSpaceRelation< ControlDescriptions, Descriptions, UXIdxMap >( controlSpace_, totalSpace_ );
    using UZIdxMap = boost::fusion::vector< IdxPair<0,1> >;
    setSubSpaceRelation< ControlDescriptions, PrimalDescriptions, UZIdxMap >( controlSpace_, primalSpace_ );
    
    
    ///////////////////////////////////////////////// Functional 
    cout << "create Functional" << endl;
    auto resultRef = zero(stateSpace_);
    auto& resultRefK = Spacy::cast_ref<Vector<StateDescriptions>>(resultRef);
    
    const Dune::FieldVector<double,3> reference_deformationVector {1,1,1};
    auto& refImpl = boost::fusion::at_c<0>(resultRefK.get().data).coefficients();
    for(int i=0; i<refImpl.size(); ++i)
        refImpl[i]=reference_deformationVector;
    
    StateDescriptions::VariableSet reference_deformation(stateDescription);
    copy(resultRef, reference_deformation);
    
    auto y0 = FunG::LinearAlgebra::unitMatrix<Dune::FieldMatrix<double,dim,dim>>();
    auto integrand  = FunG::compressibleMooneyRivlin<FunG::Pow<2>, FunG::LN, Dune::FieldMatrix<double,dim,dim>>(3.69, 0.41, 2.09, -13.20, y0);
    
    using Functional = LagrangeFunctional<decltype(integrand),
                                               decltype(reference_deformation), 0, 1,
                                               2, double, Descriptions>;
    auto fL = makeC2Functional(Functional(tychonovParameter,reference_deformation, integrand), totalSpace_);
    
    cout << "Assembling" << endl;
    auto x0 = zero(totalSpace_);
    fL.hessian(x0);
    auto b = fL.d1(x0);
    
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
    
    ///////////////////////////////////////////////// Results
    cout << "finishing ... " << std::flush;
    
    VarSet result( desc );
    copy( res,result );
    
    VarSet x_init( desc );
    copy( b,x_init );

//     auto spacy_z = primalSpace_.project( res );
//     PrimalDescriptions::VariableSet z( primalDescription );
//     copy( spacy_z, z );
// 
//     auto spacy_y = stateSpace_.project( spacy_z );
//     StateDescriptions::VariableSet y( stateDescription );
//     copy( spacy_y, y );
//     
//     auto spacy_u = controlSpace_.project( spacy_z );
//     ControlDescriptions::VariableSet u( controlDescription );
//     copy( spacy_u, u );
//     
//     auto spacy_p = adjointSpace_.project( res );
//     DualDescriptions::VariableSet p( dualDescription );
//     copy( spacy_p, p );
    
    ///////////////////////////////////////////////// VTK-Files
    IoOptions options;
    options.outputType = IoOptions::ascii;
    writeVTKFile( gm.grid().leafGridView(), x_init, "initial_values", options, FEorder );
    writeVTKFile( gm.grid().leafGridView(), result, "variables", options, FEorder );
//     writeVTKFile( gm.grid().leafGridView(), z, "primal_variables", options, FEorder );
//     writeVTKFile( gm.grid().leafGridView(), p, "dual_variables", options, FEorder );
//     writeVTKFile( gm.grid().leafGridView(), u, "control_variables", options, FEorder );
//     writeVTKFile( gm.grid().leafGridView(), y, "state_variables", options, FEorder );
    
    cout << "done" << endl;
    
    /**/
    return 0;
}
