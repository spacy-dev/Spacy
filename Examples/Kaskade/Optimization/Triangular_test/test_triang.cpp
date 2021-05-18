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
#include <Spacy/Algorithm/PrimalDualProjection/ModifiedPPCG.h>
#include <Spacy/Algorithm/PrimalDualProjection/TriangConstraintPrec.h>
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

    const auto desiredAccuracy = getParameter( pt, "desiredAccuracy", 1e-6 );
    const auto maxSteps = getParameter( pt, "maxSteps", 500 );
    const auto initialRefinements = getParameter( pt, "initialRefinements", 3 );
    const auto iterativeRefinements = getParameter( pt, "iterativeRefinements", 0 );
    const auto FEorder = getParameter( pt, "FEorder", 1 );
    const auto verbose = getParameter( pt, "verbose", 2 );
    const auto desContr = getParameter( pt, "desiredContraction", 0.5 );
    const auto relDesContr = getParameter( pt, "relaxedContraction", desContr + 0.1 );
    const auto maxContr = getParameter( pt, "maxContraction", 0.75 );
    const auto tychonovParameter = getParameter( pt, "tychonovParameter", 1 );
    const auto numberOfCubesX = getParameter (pt, "NumberOfCubesX", 1);
    const auto numberOfCubesY = getParameter (pt, "NumberOfCubesY", 1);
    const auto numberOfCubesZ = getParameter (pt, "NumberOfCubesZ", 1);
    const auto stateAcc = getParameter (pt, "stateAccuracy", 1e-3);
    const auto controlAcc = getParameter (pt, "controlAccuracy", 1e-3);

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
    //using L2Space = FEFunctionSpace<BoundaryMapper<ContinuousLagrangeMapper, double, GridView> >;
    using H1Space = FEFunctionSpace< ContinuousLagrangeMapper< double, GridView > >;
    //using Spaces = boost::fusion::vector< H1Space const*, L2Space const* >;
    using Spaces = boost::fusion::vector< FEFunctionSpace< ContinuousLagrangeMapper< double, GridView > > const*  >;//, FEFunctionSpace<BoundaryMapper<ContinuousLagrangeMapper, double, GridView> > const* >;

    using components = Components<3>;
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
    
    //L2Space l2Space(gm, gm.grid().leafGridView(), FEorder, totalIndexSet, partialIndexSetId);
    H1Space h1Space( gm, gm.grid().leafGridView(), FEorder );
    Spaces spaces( &h1Space   );//, &l2Space);

    Descriptions desc( spaces, { "y", "u", "p" } );               // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay)
    PrimalDescriptions primalDescription( spaces, { "y", "u" } ); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay)
    DualDescriptions dualDescription( spaces, { "p" } );          // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay)
    StateDescriptions stateDescription( spaces, { "y" } );          // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay)
    ControlDescriptions controlDescription( spaces, { "u" } );          // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay)
    

    ///////////////////////////////////////////////// Reference solution (eigentlich für anderes Beispielprogramm)
    cout << "interpolate" << endl;
    VarSet x_ref( desc );
    interpolateGloballyFromFunctor< PlainAverage >(
        boost::fusion::at_c< Ids::state >( x_ref.data ), []( auto const& cell, auto const& xLocal ) -> Dune::FieldVector< double, 1 > {
            const auto x = cell.geometry().global( xLocal );
            return Dune::FieldVector< double, 1 >( 12 * ( 1 - x[ 1 ] ) * x[ 1 ] * ( 1 - x[ 0 ] ) * x[ 0 ] );
        } );

    
    ///////////////////////////////////////////////// Spacy Vector-Spaces
    cout << "create Vector-Spaces" << endl;
    
    auto totalSpace = makeHilbertSpace< Descriptions >( desc, "X" );
    
    auto primalSpace = makeHilbertSpace< PrimalDescriptions >( primalDescription, "Z" );
    using ZXIdxMap = boost::fusion::vector< IdxPair<0,0>, IdxPair<1,1> >;
    setSubSpaceRelation< PrimalDescriptions, Descriptions, ZXIdxMap >( primalSpace, totalSpace );
    
    auto dualSpace = makeHilbertSpace< DualDescriptions >( dualDescription, "P" );
    using PXIdxMap = boost::fusion::vector< IdxPair<0,2> >;
    setSubSpaceRelation< DualDescriptions, Descriptions, PXIdxMap >( dualSpace, totalSpace );
    
    auto stateSpace = makeHilbertSpace< StateDescriptions >( stateDescription, "Y" );
    using YXIdxMap = boost::fusion::vector< IdxPair<0,0> >;
    setSubSpaceRelation< StateDescriptions, Descriptions, YXIdxMap >( stateSpace, totalSpace );
    using YZIdxMap = boost::fusion::vector< IdxPair<0,0> >;
    setSubSpaceRelation< StateDescriptions, PrimalDescriptions, YZIdxMap >( stateSpace, primalSpace );
    
    auto controlSpace = makeHilbertSpace< ControlDescriptions >( controlDescription, "U" );
    using UXIdxMap = boost::fusion::vector< IdxPair<0,1> >;
    setSubSpaceRelation< ControlDescriptions, Descriptions, UXIdxMap >( controlSpace, totalSpace );
    using UZIdxMap = boost::fusion::vector< IdxPair<0,1> >;
    setSubSpaceRelation< ControlDescriptions, PrimalDescriptions, UZIdxMap >( controlSpace, primalSpace );
    
    
    
    ///////////////////////////////////////////////// Functional 
    cout << "create Functional" << endl;
    auto resultRef = zero(stateSpace);
    StateDescriptions::VariableSet reference_deformation(stateDescription);
    copy(resultRef, reference_deformation);     
    
    auto y0 = FunG::LinearAlgebra::unitMatrix<Dune::FieldMatrix<double,dim,dim>>();
    auto integrand  = FunG::compressibleMooneyRivlin<FunG::Pow<2>, FunG::LN, Dune::FieldMatrix<double,dim,dim>>(3.69, 0.41, 2.09, -13.20, y0);
    
    using Functional = LagrangeFunctional<decltype(integrand),
                                               decltype(reference_deformation), 0, 1,
                                               2, double, Descriptions>;
    auto fL = makeC2Functional(Functional(tychonovParameter,reference_deformation, integrand), totalSpace);
    
    cout << "Assembling" << endl;
    auto x0 = zero(totalSpace);
    fL.hessian(x0);
    
    ///////////////////////////////////////////////// Blocks, Solvers
    cout << "create Blocks and Solvers" << endl;
    
    using KMAT = ::Kaskade::MatrixAsTriplet<double>;
    
    auto minusB         = makeCallableBlock<Functional>(fL,controlSpace,dualSpace);
    auto stateSolver    = makeDirectSolver<DualDescriptions,StateDescriptions>
                            (makeMatrixBlock<decltype(fL),KMAT,DualDescriptions,StateDescriptions>(fL,dualSpace,stateSpace).get(),dualSpace,stateSpace);
    auto controlSolver  = makeDirectSolver<ControlDescriptions,ControlDescriptions>
                            (makeSymmetricMatrixBlock<decltype(fL),KMAT,ControlDescriptions>(fL,controlSpace).get(),controlSpace, controlSpace); 
    auto adjointSolver  = makeDirectSolver<StateDescriptions,DualDescriptions>
                            (makeMatrixBlock<decltype(fL),KMAT,StateDescriptions,DualDescriptions>(fL,stateSpace,dualSpace).get(),stateSpace,dualSpace);
                            
    ///////////////////////////////////////////////// Execute
    cout << "Execute" << endl;
    
    ::Spacy::ModifiedPPCG::TriangularConstraintPreconditioner triH(stateSolver,controlSolver,adjointSolver,minusB,totalSpace,stateSpace);
    auto temp = zero(totalSpace);
    copy(x_ref,temp);
    auto res = triH(temp);
    
    ///////////////////////////////////////////////// Results
    cout << "finishing" << endl;
    
    VarSet r( desc );
    copy( res,r );

    auto spacy_z = primalSpace.project( res );
    PrimalDescriptions::VariableSet z( primalDescription );
    copy( spacy_z, z );

    auto spacy_y = stateSpace.project( spacy_z );
    StateDescriptions::VariableSet y( stateDescription );
    copy( spacy_y, y );
    
    auto spacy_u = controlSpace.project( spacy_z );
    ControlDescriptions::VariableSet u( controlDescription );
    copy( spacy_u, u );
    
    auto spacy_p = dualSpace.project( res );
    DualDescriptions::VariableSet p( dualDescription );
    copy( spacy_p, p );
    
    ///////////////////////////////////////////////// VTK-Files
    IoOptions options;
    options.outputType = IoOptions::ascii;
    writeVTKFile( gm.grid().leafGridView(), x_ref, "reference", options, FEorder );
    writeVTKFile( gm.grid().leafGridView(), r, "variables", options, FEorder );
    writeVTKFile( gm.grid().leafGridView(), z, "primal_variables", options, FEorder );
    writeVTKFile( gm.grid().leafGridView(), p, "dual_variables", options, FEorder );
    writeVTKFile( gm.grid().leafGridView(), u, "control_variables", options, FEorder );
    writeVTKFile( gm.grid().leafGridView(), y, "state_variables", options, FEorder );
    
    /**/
    return 0;
}
