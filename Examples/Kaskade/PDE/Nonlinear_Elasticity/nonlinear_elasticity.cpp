#include <iostream>

#include <dune/grid/config.h>
#include <dune/grid/uggrid.hh>

#define SPACY_ENABLE_LOGGING
#include <Spacy/Adapter/kaskade.hh>
//#include <Spacy/Spacy.h>

// Algorithms
#include <Spacy/Algorithm/ACR/ACR.h>
#include <Spacy/Algorithm/CG/CG.h>
#include <Spacy/Algorithm/CG/LinearSolver.h>
#include <Spacy/Algorithm/CG/TerminationCriteria.h>
#include <Spacy/Algorithm/CG/TerminationCriterion.h>
#include <Spacy/Algorithm/DampingFactor.h>
#include <Spacy/Algorithm/LipschitzConstant.h>
#include <Spacy/Algorithm/Newton/Newton.h>
#include <Spacy/Algorithm/Newton/TerminationCriteria.h>

#include <Spacy/Spaces/RealSpace.h>

// Util
#include <Spacy/Util/Cast.h>
#include <Spacy/Util/Copy.h>
#include <Spacy/Util/Invoke.h>
#include <Spacy/Util/Log.h>
#include <Spacy/Util/Mixins.h>
#include <Spacy/Util/Voider.h>

// Interfaces and directly related functionality
#include <Spacy/C1Functional.h>
#include <Spacy/C1Operator.h>
#include <Spacy/C2Functional.h>
#include <Spacy/Derivative.h>
#include <Spacy/DynamicOperator.h>
#include <Spacy/Functional.h>
#include <Spacy/HilbertSpaceNorm.h>
#include <Spacy/InducedScalarProduct.h>
#include <Spacy/LinearOperator.h>
#include <Spacy/LinearSolver.h>
#include <Spacy/Norm.h>
#include <Spacy/Operator.h>
#include <Spacy/ScalarProduct.h>
#include <Spacy/Vector.h>
#include <Spacy/VectorSpace.h>
#include <Spacy/ZeroVectorCreator.h>


#include <fem/gridmanager.hh>
#include <fem/lagrangespace.hh>
#include <fem/variables.hh>
#include <io/vtk.hh>
#include <utilities/gridGeneration.hh> //  createUnitSquare

#include "../fung_functional.hh"

#include <fung/examples/biomechanics/adipose_tissue_sommer_holzapfel.hh>
#include <fung/examples/rubber/neo_hooke.hh>
#include <fung/fung.hh>

using namespace Kaskade;

int main()
{
    SET_LOG_PRINTER( Spacy::Log::StreamPrinter() );

    constexpr int dim = 3;
    int refinements = 4;
    int order = 1;

    using Grid = Dune::UGGrid< dim >;
    using H1Space = FEFunctionSpace< ContinuousLagrangeMapper< double, Grid::LeafGridView > >;
    using Spaces = boost::fusion::vector< H1Space const* >;
    using VariableDescriptions =
        boost::fusion::vector< Variable< SpaceIndex< 0 >, Components< dim >, VariableId< 0 > > >;
    using VariableSetDesc = VariableSetDescription< Spaces, VariableDescriptions >;

    GridManager< Grid > gridManager( createUnitCube< Grid >( 1., false ) );
    gridManager.globalRefine( refinements );

    H1Space temperatureSpace( gridManager, gridManager.grid().leafGridView(), order );

    Spaces spaces( &temperatureSpace );
    VariableSetDesc variableSetDesc( spaces, {"u"} );

    using Matrix = Dune::FieldMatrix< double, dim, dim >;
    // fiber tensor for 'fibers' aligned with x-axes
    auto fiberTensor = Matrix( 0 );
    fiberTensor[ 0 ][ 0 ] = 1;
    auto y0 = FunG::LinearAlgebra::unitMatrix< Matrix >();
    //    auto integrand = FunG::incompressibleNeoHooke(1.0,y0);
    auto integrand = FunG::incompressibleAdiposeTissue_SommerHolzapfel( fiberTensor, y0 );
    auto F = FungFunctional< decltype( integrand ), VariableSetDesc >( integrand );

    // compute solution
    auto domain = Spacy::Kaskade::makeHilbertSpace< VariableSetDesc >( variableSetDesc );
    auto range = Spacy::Kaskade::makeHilbertSpace< VariableSetDesc >( variableSetDesc );
    connectAsPrimalDualPair( domain, range );

    Spacy::C1Operator A = Spacy::Kaskade::makeC1Operator( F, domain, range );
    domain.setScalarProduct( Spacy::InducedScalarProduct( A.linearization( zero( domain ) ) ) );

    auto p = Spacy::Newton::Parameter{};
    p.setRelativeAccuracy( 1e-6 );

    auto x = covariantNewton( A, p );
    //  auto x = Spacy::contravariantNewton(A,p);
    //  auto x = Spacy::localNewton(A,p);

    // construct Galerkin representation
 //   writeVTK( Spacy::cast_ref< Spacy::Kaskade::Vector< VariableSetDesc > >( x ), "displacement" );
}
