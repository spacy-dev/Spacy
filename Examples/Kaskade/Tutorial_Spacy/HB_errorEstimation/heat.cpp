#include "dune/grid/config.h"
#include "dune/grid/uggrid.hh"

#include "fem/assemble.hh"
#include "fem/lagrangespace.hh"
#include "io/vtk.hh"
//#include "io/amira.hh"
#include "linalg/direct.hh"
#include "utilities/enums.hh"
#include "utilities/gridGeneration.hh" //  createUnitSquare, createUnitCube
#include "utilities/kaskopt.hh"
#include "utilities/kaskopt.hh"

#include <Spacy/Spacy.h>
#include <Spacy/Adapter/kaskade.hh>

#include <boost/timer/timer.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

using namespace Kaskade;

#include "poisson.hh"

constexpr int dim=2;

using Grid = Dune::UGGrid<dim>;
using LeafView = Grid::LeafGridView;
using H1Space = FEFunctionSpace<ContinuousLagrangeMapper<double,LeafView> >;
using Spaces = boost::fusion::vector<H1Space const*>;
using VariableDescriptions = boost::fusion::vector<Variable<SpaceIndex<0>,Components<1>,VariableId<0> >>;
using VSD = VariableSetDescription<Spaces,VariableDescriptions>;
using Functional = PoissonFunctional<double,VSD>;
using Assembler = VariationalFunctionalAssembler<LinearizationAt<Functional> >;
constexpr int neq = Functional::TestVars::noOfVariables;
using CoefficientVectors = VSD::CoefficientVectorRepresentation<0,neq>::type;

using namespace boost::fusion;

std::string getSaveFileName(int refSteps)
{
    std::ostringstream fn;
    fn << "graph/peak-grid";
    fn.width(3);
    fn.fill('0');
    fn.setf(std::ios_base::right,std::ios_base::adjustfield);
    fn << refSteps;
    fn.flush();
    return fn.str();
}

int main(int argc, char *argv[])
{
    std::cout << "Start heat transfer tutorial program using HB error estimation (SpaceDimension=2)"
              << std::endl;
    boost::timer::cpu_timer totalTimer;

    int verbosityOpt = 1;
    bool dump = true;
    std::unique_ptr<boost::property_tree::ptree> pt = getKaskadeOptions(argc, argv, verbosityOpt, dump);

    int refinements = getParameter(pt, "refinements", 5),
            order =  getParameter(pt, "order", 1),
            maxAdaptSteps =  getParameter(pt, "maxAdaptSteps", 5),
            verbosity = getParameter(pt, "verbosity", 1);
    std::cout << "original mesh shall be refined : " << refinements << " times" << std::endl;
    std::cout << "discretization order           : " << order << std::endl;
    std::cout << "max. adaptive refine steps     : " << maxAdaptSteps << std::endl;
    std::cout << "output level (verbosity)       : " << verbosity << std::endl << std::endl;

    GridManager<Grid> gridManager( createUnitSquare<Grid>() );
    gridManager.globalRefine(refinements);
    gridManager.setVerbosity(verbosity);

    // construction of finite element space for the scalar solution T
    H1Space temperatureSpace(gridManager,gridManager.grid().leafGridView(),order);
    Spaces spaces(&temperatureSpace);
    std::string varNames[1] = { "T" };
    VSD variableSetDescription(spaces,varNames);
    Functional F;

    auto X = Spacy::Kaskade::makeHilbertSpace< VSD >( variableSetDescription );
    Spacy::globalSpaceManager().add(X.index(), Spacy::Kaskade::refineCells(gridManager));

    const auto A = Spacy::Kaskade::makeC1Operator( F, X, X.dualSpace() );
    const auto minRefinementRatio = 0.2;
    const auto gamma = 1;
    const auto tolerance = std::sqrt(1-1.0/(dim*gamma))*1e-3;
    auto estimator = Spacy::Kaskade::getHierarchicalErrorEstimator(F, gridManager, variableSetDescription,
                                                                   minRefinementRatio, tolerance);

    auto estimateAndRefine = [estimator = std::move(estimator)](const Spacy::Vector& x, const Spacy::Vector& dx) mutable
    {
        static auto refSteps = 0;
        const auto doRefine = estimator.estimateError(x, dx);

        const_cast<Spacy::Vector&>(x) += dx;
        Spacy::Kaskade::writeVTK(Spacy::cast_ref<Spacy::Kaskade::Vector<VSD>>(x),
                                 getSaveFileName(refSteps++));
        if(doRefine)
            Spacy::globalSpaceManager().adjustDiscretization(x.space().index(),
                                                             estimator.getErrorIndicator());
        return doRefine;
    };

    Spacy::Newton::Parameter p;
    p.setMaxSteps(maxAdaptSteps + 1);
    Spacy::localNewton(A, p, estimateAndRefine);

    std::cout << "total computing time: " << boost::timer::format(totalTimer.elapsed()) << "\n";
    std::cout << "End heat transfer (peak source) tutorial program" << std::endl;
}
