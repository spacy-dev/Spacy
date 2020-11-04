#include "LagrangeFunctional.h"
#include "NormalStepFunctional.h"
#include <dolfin.h>

#include <Spacy/Adapter/fenics.hh>
#include <Spacy/Algorithm/CompositeStep/affineCovariantSolver.hh>

#include <chrono>

using namespace dolfin;
using namespace std::chrono;

// Source term reference solution
class Source : public Expression
{
    void eval( Array< double >& values, const Array< double >& x ) const
    {
        double dx = x[ 0 ] - 0.5;
        double dy = x[ 1 ] - 0.5;
        values[ 0 ] = 1; // 10*exp(-(dx*dx + dy*dy) / 0.02);
    }
};

// Sub domain for Dirichlet boundary condition
class DirichletBoundary : public SubDomain
{
    bool inside( const Array< double >& x, bool on_boundary ) const
    {
        return x[ 0 ] < DOLFIN_EPS or x[ 0 ] > 1.0 - DOLFIN_EPS or x[ 1 ] < DOLFIN_EPS or x[ 1 ] > 1. - DOLFIN_EPS;
    }
};

int main()
{
    // Create mesh and function space
    auto n = 256;
    auto mesh = std::make_shared< UnitSquareMesh >( n, n );
    auto V = std::make_shared< LagrangeFunctional::FunctionSpace >( mesh );

    // Define variational forms
    LagrangeFunctional::Functional L( mesh );
    LagrangeFunctional::LinearForm dL( V );
    LagrangeFunctional::BilinearForm ddL( V, V );
    NormalStepFunctional::BilinearForm Norm( V, V );

    auto dummy = std::make_shared< Function >( V );
    auto y_ref = std::make_shared< Source >();
    auto alpha = std::make_shared< Constant >( 1e-2 );
    auto c = std::make_shared< Constant >( 1e-1 );
    auto d = std::make_shared< Constant >( 1e1 );
    L.y0 = y_ref;
    L.alpha = alpha;
    L.c = c;
    L.d = d;
    dL.y0 = y_ref;
    dL.alpha = alpha;
    dL.c = c;
    dL.d = d;
    ddL.alpha = alpha;
    ddL.c = c;
    ddL.d = d;
    Norm.alpha = alpha;
    Norm.c = c;
    Norm.d = d;

    // Compute solution
    using namespace Spacy;

    std::vector< unsigned > primalSpaceIds = { 0u, 1u }, dualSpaceIds = { 2u };
    auto productSpace = FEniCS::makeHilbertSpace( V, primalSpaceIds, dualSpaceIds );

    // functionals
    Spacy::C2Functional lagrangeFunctional = FEniCS::makeC2Functional( L, dL, ddL, productSpace );
    Spacy::C2Functional normalStepFunctional = FEniCS::makeC2Functional( L, dL, Norm, productSpace );

    // composite step solve
    CompositeStep::AffineCovariantSolver solver( normalStepFunctional, lagrangeFunctional, productSpace );
    solver.setRelativeAccuracy( 1e-9 );
    solver.setEps( 1e-12 );
    solver.setMaxSteps( 500 );
    solver.setVerbosityLevel( 2 );
    solver.setIterativeRefinements( 0 );

    // solve problem
    auto startTime = high_resolution_clock::now();
    auto sol = solver();
    std::cout << "computation time: " << duration_cast< seconds >( high_resolution_clock::now() - startTime ).count() << std::endl;

    // copy back and display solution
    Function u( V );
    FEniCS::copy( sol, u );

    File filecs( "cs.pvd" );
    filecs << u[ 0 ];

    // Plot solution
    plot( u );
    interactive();

    return 0;
}
