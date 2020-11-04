/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2018 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "fem/assemble.hh"
#include "fem/spaces.hh"
#include "io/matlab.hh"
#include "io/vtk.hh"
#include "linalg/apcg.hh"
#include "linalg/direct.hh"
#include "linalg/threadedMatrix.hh"
#include "linalg/triplet.hh"
#include "mg/additiveMultigrid.hh"
#include "mg/multigrid.hh"
#include "utilities/enums.hh"
#include "utilities/gridGeneration.hh"
#include "utilities/kaskopt.hh"
#include "utilities/memory.hh"
#include "utilities/timing.hh"

#include "dune/grid/config.h"
#include "dune/grid/uggrid.hh"
#include <dune/istl/operators.hh>

#include <iostream>

using namespace Kaskade;
#include "elastomechanics.hh"

namespace Kaskade
{
    std::mutex ndLock;
}

int main( int argc, char* argv[] )
{
    using namespace boost::fusion;

    Timings& timer = Timings::instance();
    timer.start( "total computing time" );
    std::cout << "Start elastomechanics tutorial program" << std::endl;

    auto coarseGridSize = 1;
    auto maxit = 100;
    auto order = 1;
    auto refinements = 3;
    auto skipEntries = 100;
    auto skipDimension = 32;
    auto solver = 0;
    auto verbose = 0;
    auto threads = -1;
    auto additive = false;
    auto vtk = true;
    auto direct = true;
    auto atol = 1e-8;
    auto epsilon = 0.0;
    std::string material;
    std::string storageScheme;
    std::string storageType;
    std::string factorizationType;

    if ( getKaskadeOptions(
             argc, argv,
             Options( "refinements", refinements, refinements, "number of uniform grid refinements" )(
                 "coarse", coarseGridSize, coarseGridSize,
                 "number of coarse grid elements along each cube edge" )( "order", order, order, "finite element ansatz order" )(
                 "material", material, "steel", "type of material" )( "direct", direct, direct, "if true, use a direct solver" )(
                 "solver", solver, solver, "0=UMFPACK, 1=PARDISO 2=MUMPS 3=SUPERLU 4=UMFPACK32/64 5=UMFPACK64" )(
                 "additive", additive, additive, "use additive multigrid" )( "verbosity", verbose, verbose, "amount of reported details" )(
                 "vtk", vtk, true, "write solution to VTK file" )( "atol", atol, atol,
                                                                   "absolute energy error tolerance for iterative solver" )(
                 "maxit", maxit, maxit, "maximum number of iterations" )( "eps", epsilon, epsilon,
                                                                          "relative accuracy for artificial inexactness" )(
                 "storageScheme", storageScheme, "ND", "A, L, or ND" )( "storageType", storageType, "float", "float or double" )(
                 "factorizationType", factorizationType, "double",
                 "float or double" )( "threads", threads, -1, "number of threads (<0: select default based on cpu cores)" )(
                 "skipE", skipEntries, skipEntries, "number of zero block entries to skip when dissecting" )(
                 "skipD", skipDimension, skipDimension, "minimal system size for dissection" ) ) )
        return 1;

    std::cout << "refinements of original mesh : " << refinements << std::endl;
    std::cout << "discretization order         : " << order << std::endl;

    if ( threads > 0 )
        NumaThreadPool::instance( threads );
    std::cout << "NumaThreadPool::isSequential() = " << NumaThreadPool::instance().isSequential() << " [threads=" << threads << "]\n";

    constexpr int DIM = 3;
    using Grid = Dune::UGGrid< DIM >;
    using Spaces = boost::fusion::vector< H1Space< Grid > const* >;

    timer.start( "computing time for grid creation & refinement" );
    Dune::FieldVector< double, DIM > c0( 0.0 );
    Dune::FieldVector< double, DIM > dc( 1.0 );
    GridManager< Grid > gridManager( createCuboid< Grid >( c0, dc, 1.0 / coarseGridSize, true ) );
    gridManager.globalRefine( refinements );
    gridManager.enforceConcurrentReads( true );
    timer.stop( "computing time for grid creation & refinement" );

    // construction of finite element space for the scalar solution T.
    H1Space< Grid > h1Space( gridManager, gridManager.grid().leafGridView(), order );

    auto varSetDesc =
        makeVariableSetDescription( makeSpaceList( &h1Space ), make_vector( Variable< SpaceIndex< 0 >, Components< DIM > >( "u" ) ) );
    using VarSetDesc = decltype( varSetDesc );

    using Functional = ElasticityFunctional< VarSetDesc >;
    using Assembler = VariationalFunctionalAssembler< LinearizationAt< Functional > >;
    using CoefficientVectors = VarSetDesc::CoefficientVectorRepresentation< 0, 1 >::type;

    // Create the variational functional.
    Functional F( ElasticModulus::material( material ) );

    // construct Galerkin representation
    Assembler assembler( varSetDesc.spaces );
    VarSetDesc::VariableSet x( varSetDesc );
    VarSetDesc::VariableSet dx( varSetDesc );

    timer.start( "computing time for assemble" );
    assembler.assemble( linearization( F, x ) );
    timer.stop( "computing time for assemble" );

    CoefficientVectors rhs( assembler.rhs() );
    CoefficientVectors solution = varSetDesc.zeroCoefficientVector();
    timer.start( "computing time for solve" );
    if ( direct )
    {
        auto directType = static_cast< DirectType >( solver );
        AssembledGalerkinOperator< Assembler > A( assembler, directType == DirectType::MUMPS || directType == DirectType::PARDISO );
        directInverseOperator( A, directType, MatrixProperties::POSITIVEDEFINITE ).applyscaleadd( -1.0, rhs, solution );
        x.data = solution.data;
    }
    else
    {
        timer.start( "preconditioner setup" );
        using X = Dune::BlockVector< Dune::FieldVector< double, DIM > >;
        DefaultDualPairing< X, X > dp;
        using Matrix = NumaBCRSMatrix< Dune::FieldMatrix< double, DIM, DIM > >;
        using LinOp = Dune::MatrixAdapter< Matrix, X, X >;
        Matrix Amat( assembler.get< 0, 0 >(), true );
        LinOp A( Amat );
        SymmetricLinearOperatorWrapper< X, X > sa( A, dp );
        PCGEnergyErrorTerminationCriterion< double > term( atol, maxit );

        Dune::InverseOperatorResult res;
        X xi( component< 0 >( rhs ).N() );

        std::unique_ptr< SymmetricPreconditioner< X, X > > mg;
        if ( additive )
        {
            if ( order == 1 )
                mg = moveUnique( makeBPX( Amat, gridManager ) );
            else
            {
                H1Space< Grid > p1Space( gridManager, gridManager.grid().leafGridView(), 1 );
                if ( storageScheme == "A" )
                {
                    if ( storageType == "float" )
                    {
                        mg = moveUnique(
                            makePBPX( Amat, h1Space, p1Space, DenseInverseStorageTag< float >(), gridManager.grid().maxLevel() ) );
                    }
                    else
                    {
                        mg = moveUnique(
                            makePBPX( Amat, h1Space, p1Space, DenseInverseStorageTag< double >(), gridManager.grid().maxLevel() ) );
                    }
                }
                else if ( storageScheme == "L" )
                {
                    if ( storageType == "float" )
                    {
                        mg = moveUnique(
                            makePBPX( Amat, h1Space, p1Space, DenseCholeskyStorageTag< float >(), gridManager.grid().maxLevel() ) );
                    }
                    else
                    {
                        mg = moveUnique(
                            makePBPX( Amat, h1Space, p1Space, DenseCholeskyStorageTag< double >(), gridManager.grid().maxLevel() ) );
                    }
                }
                else if ( storageScheme == "ND" )
                {
                    if ( storageType == "float" )
                    {
                        if ( factorizationType == "float" )
                        {
                            mg = moveUnique( makePBPX( Amat, h1Space, p1Space, NestedDissectionStorageTag< float, float >( skipEntries ),
                                                       gridManager.grid().maxLevel() ) );
                        }
                        else
                        {
                            mg = moveUnique( makePBPX( Amat, h1Space, p1Space, NestedDissectionStorageTag< float, double >( skipEntries ),
                                                       gridManager.grid().maxLevel() ) );
                        }
                    }
                    else
                    {
                        mg = moveUnique( makePBPX( Amat, h1Space, p1Space, NestedDissectionStorageTag< double, double >( skipEntries ),
                                                   gridManager.grid().maxLevel() ) );
                    }
                }
                else
                {
                    std::cerr << "unknown storage scheme provided\n";
                    return -1;
                }
            }
        }
        else
        {
            mg = makeMultigrid( Amat, h1Space );
        }
        timer.stop( "preconditioner setup" );

        timer.start( "PCG iteration" );
        Pcg< X, X > pcg( sa, *mg, term, verbose );
        pcg.apply( xi, component< 0 >( rhs ), res );
        std::cout << "PCG iterations: " << res.iterations << "\n";
        timer.stop( "PCG iteration" );
        xi *= -1;
        component< 0 >( x ) = xi;
    }
    timer.stop( "computing time for solve" );

    // output of solution in VTK format for visualization,
    // the data are written as ascii stream into file elasto.vtu,
    // possible is also binary
    if ( vtk )
    {
        ScopedTimingSection ts( "computing time for file i/o", timer );
        writeVTKFile( x, "elasto", IoOptions().setOrder( order ).setPrecision( 7 ) );
    };

    timer.stop( "total computing time" );
    if ( verbose > 0 )
        std::cout << timer;
    std::cout << "End elastomechanics tutorial program" << std::endl;

    return 0;
}
