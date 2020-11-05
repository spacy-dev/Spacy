#include <Spacy/Adapter/kaskade.hh>
#include <Spacy/Spacy.h>

#include <fem/assemble.hh>
#include <fem/spaces.hh>
#include <io/matlab.hh>
#include <io/vtk.hh>
#include <linalg/apcg.hh>
#include <linalg/direct.hh>
#include <linalg/threadedMatrix.hh>
#include <linalg/triplet.hh>
#include <mg/multigrid.hh>
#include <utilities/enums.hh>
#include <utilities/gridGeneration.hh>
#include <utilities/kaskopt.hh>
#include <utilities/memory.hh>

#include <dune/grid/config.h>
#include <dune/grid/uggrid.hh>
#include <dune/istl/operators.hh>

#include <iostream>

using namespace Kaskade;
#include "elastomechanics.hh"

constexpr int DIM = 3;
using Grid = Dune::UGGrid< DIM >;
using Spaces = boost::fusion::vector< H1Space< Grid > const* >;
using NumaMatrix = NumaBCRSMatrix< Dune::FieldMatrix< double, DIM, DIM > >;

template < class VariableSetDescription, class Matrix >
std::function< ::Spacy::LinearSolver( const ::Spacy::Kaskade::LinearOperator< VariableSetDescription, VariableSetDescription, Matrix >& ) >
getMultiGridSolver( const H1Space< Grid >& h1Space, double atol, int maxit, int verbose )
{
    return [ & ]( const Spacy::Kaskade::LinearOperator< VariableSetDescription, VariableSetDescription, Matrix >& L ) {
        return [ &h1Space, atol, maxit, verbose, A = L.get() ]( const Spacy::Vector& x ) -> Spacy::Vector {
            using CoefficientVector = typename VariableSetDescription::template CoefficientVectorRepresentation<>::type;
            const auto spaces = ::Kaskade::makeSpaceList( &h1Space );
            CoefficientVector y_( VariableSetDescription::template CoefficientVectorRepresentation<>::init( spaces ) );
            CoefficientVector x_( VariableSetDescription::template CoefficientVectorRepresentation<>::init( spaces ) );
            Spacy::Kaskade::copyToCoefficientVector< VariableSetDescription >( x, x_ );

            using X = Dune::BlockVector< Dune::FieldVector< double, DIM > >;
            Kaskade::DefaultDualPairing< X, X > dp{};
            auto mat = Dune::MatrixAdapter< std::decay_t< std::decay_t< decltype( A ) > >, X, X >( A );
            Kaskade::SymmetricLinearOperatorWrapper< X, X > sa( mat, dp );
            Kaskade::PCGEnergyErrorTerminationCriterion< double > term( atol, maxit );
            std::unique_ptr< Kaskade::SymmetricPreconditioner< X, X > > mg =
                ::Kaskade::makeMultigrid( A, h1Space, false, Kaskade::DirectType::UMFPACK );
            Kaskade::Pcg< X, X > pcg( sa, *mg, term, verbose );
            Dune::InverseOperatorResult res{};
            pcg.apply( boost::fusion::at_c< 0 >( y_.data ), boost::fusion::at_c< 0 >( x_.data ), res );

            auto y = zero( x.space() );
            Spacy::Kaskade::copyFromCoefficientVector< VariableSetDescription >( y_, y );
            return y;
        };
    };
}

int main( int argc, char* argv[] )
{
    using namespace boost::fusion;

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

    Dune::FieldVector< double, DIM > c0( 0.0 );
    Dune::FieldVector< double, DIM > dc( 1.0 );
    GridManager< Grid > gridManager( createCuboid< Grid >( c0, dc, 1.0 / coarseGridSize, true ) );
    gridManager.globalRefine( refinements );
    gridManager.enforceConcurrentReads( true );

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

    // compute solution
    const auto X = Spacy::Kaskade::makeHilbertSpace< VarSetDesc >( varSetDesc );
    auto createMultiGridSolver = getMultiGridSolver< VarSetDesc, NumaMatrix >( h1Space, atol, maxit, verbose );
    auto A = direct ? Spacy::Kaskade::makeC2Functional< Functional, NumaMatrix >( F, X )
                    : Spacy::Kaskade::makeC2Functional< Functional, NumaMatrix >( F, X, std::move( createMultiGridSolver ) );
    const auto x0 = zero( X );
    const auto solver_ = A.hessian( x0 ).solver();
    const auto x = solver_( -A.d1( x0 ) );

    // output
    writeVTK( Spacy::cast_ref< Spacy::Kaskade::Vector< VarSetDesc > >( x ), "elasto" );
    std::cout << "graphical output finished, data in VTK format is written into file elasto.vtu \n";
    std::cout << "End elastomechanics tutorial program" << std::endl;
}
