#include "../fung_operator.hh"
#include <fung/examples/nonlinear_heat.hh>
#include <fung/fung.hh>

#include <Spacy/Adapter/kaskade.hh>
#include <Spacy/Spacy.h>

#include <fem/gridmanager.hh>
#include <fem/lagrangespace.hh>
#include <fem/variables.hh>
#include <io/vtk.hh>
#include <utilities/gridGeneration.hh> //  createUnitSquare

#include <dune/grid/config.h>
#include <dune/grid/uggrid.hh>

#include <iostream>

using namespace Kaskade;

int main()
{
    constexpr int dim = 2;
    int refinements = 6;
    int order = 1;

    using Grid = Dune::UGGrid< dim >;
    using H1Space = FEFunctionSpace< ContinuousLagrangeMapper< double, Grid::LeafGridView > >;
    using Spaces = boost::fusion::vector< H1Space const* >;
    using VariableDescriptions = boost::fusion::vector< Variable< SpaceIndex< 0 >, Components< 1 >, VariableId< 0 > > >;
    using VariableSetDesc = VariableSetDescription< Spaces, VariableDescriptions >;

    GridManager< Grid > gridManager( createUnitSquare< Grid >( 1., false ) );
    gridManager.globalRefine( refinements );

    H1Space temperatureSpace( gridManager, gridManager.grid().leafGridView(), order );

    Spaces spaces( &temperatureSpace );
    VariableSetDesc variableSetDesc( spaces, { "u" } );

    const auto c = 1e-1;
    const auto d = 1e2;
    Dune::FieldVector< double, 1 > u0{ 0 };
    Dune::FieldMatrix< double, 1, dim > du0{ 0 };
    auto integrand = FunG::heatModel( c, d, u0, du0 );
    auto F = FungOperator< decltype( integrand ), VariableSetDesc >( integrand );

    // compute solution
    auto domain = Spacy::Kaskade::makeHilbertSpace< VariableSetDesc >( variableSetDesc, "X" );

    Spacy::C1Operator A = Spacy::Kaskade::makeC1Operator( F, domain, domain.dualSpace() );
    domain.setScalarProduct( Spacy::InducedScalarProduct( A.linearization( zero( domain ) ) ) );

    Spacy::Newton::Parameter p{};
    p.setRelativeAccuracy( 1e-12 );
    p.setEps( 1e-15 );

    auto x = covariantNewton( A, p );
    //  auto x = Spacy::contravariantNewton(A,p);
    //  auto x = Spacy::localNewton(A,p);

    // construct Galerkin representation
    writeVTK( Spacy::cast_ref< Spacy::Kaskade::Vector< VariableSetDesc > >( x ), "temperature" );
}
