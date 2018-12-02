/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2015 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "dune/grid/config.h"
#include "dune/grid/uggrid.hh"

#include "fem/assemble.hh"
#include "fem/istlinterface.hh"
#include "fem/functional_aux.hh"
#include "fem/gridmanager.hh"
#include "fem/lagrangespace.hh"
#include "fem/variables.hh"
#include "io/vtk.hh"
#include "linalg/direct.hh"
#include "utilities/gridGeneration.hh" //  createUnitSquare

#include <Spacy/Spacy.h>
#include <Spacy/Adapter/kaskade.hh>
#include <iostream>

using namespace Kaskade;
#include "laplace.hh"

constexpr int dim         = 2;
constexpr int refinements = 5;
constexpr int order       = 2;

using Grid = Dune::UGGrid<dim>;
using LeafView = Grid::LeafGridView;
using H1Space = FEFunctionSpace<ContinuousLagrangeMapper<double,LeafView> >;
using Spaces = boost::fusion::vector<H1Space const*>;
using VariableDescriptions = boost::fusion::vector<Variable<SpaceIndex<0>,Components<1>,VariableId<0> > >;
using VariableSetDesc = VariableSetDescription<Spaces,VariableDescriptions>;
using Functional = HeatFunctional<double,VariableSetDesc>;

int main()
{
  std::cout << "Start Laplacian tutorial program" << std::endl;

  GridManager<Grid> gridManager( createUnitSquare<Grid>() );
  gridManager.globalRefine(refinements);

  // construction of finite element space for the scalar solution u.
  H1Space temperatureSpace(gridManager,gridManager.grid().leafGridView(),order);
  Spaces spaces(&temperatureSpace);
  VariableSetDesc variableSetDesc(spaces,{ "u" });
  Functional F;

  // compute solution
  const auto X = Spacy::Kaskade::makeHilbertSpace< VariableSetDesc >( variableSetDesc );
  const auto A = Spacy::Kaskade::makeC1Operator( F, X, X.dualSpace() );
  const auto solver = A.linearization(zero(X)).solver();
  const auto x = solver(-A(zero(X)));

  // output
  writeVTK( Spacy::cast_ref< Spacy::Kaskade::Vector< VariableSetDesc > >( x ), "temperature" );
  std::cout << "graphical output finished, data in VTK format is written into file temperature.vtu \n";
  std::cout << "End Laplacian tutorial program" << std::endl;
}
