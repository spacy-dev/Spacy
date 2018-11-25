#include <Spacy/Adapter/DealII.h>
#include <Spacy/Spacy.h>

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <fung/examples/nonlinear_heat.hh>
#include <fung/fung.hh>

using namespace dealii;

template < int dim >
class RightHandSide : public Function< dim >
{
public:
    RightHandSide() : Function< dim >()
    {
    }
    double value( const Point< dim >& p, const unsigned int = 0 ) const override
    {
        double return_value = 0.0;
        for ( unsigned int i = 0; i < dim; ++i )
            return_value += 4.0 * std::pow( p( i ), 4.0 );
        return return_value;
    }
};

template < int dim >
class BoundaryValues : public Function< dim >
{
public:
    BoundaryValues() : Function< dim >()
    {
    }
    double value( const Point< dim >& p, const unsigned int = 0 ) const override
    {
        return p.square();
    }
};

const constexpr int variable_dim = 1;
using VariableDims = Spacy::dealII::VariableDim< variable_dim >;
const constexpr auto dim = 2;

Spacy::Vector getRhs( const Spacy::VectorSpace& V )
{
    const auto& creator = Spacy::creator< Spacy::dealII::VectorCreator< dim, variable_dim > >( V );
    auto v = zero( V );
    auto& v_dealii = Spacy::cast_ref< Spacy::dealII::Vector >( v ).get();

    const auto& dof_handler = creator.dofHandler();
    const auto& fe = creator.finiteElement();
    QGauss< dim > quadrature_formula( 2 );
    const RightHandSide< dim > right_hand_side;
    FEValues< dim > fe_values( fe, quadrature_formula, update_values | update_gradients |
                                                           update_quadrature_points |
                                                           update_JxW_values );
    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points = quadrature_formula.size();
    Vector< double > cell_rhs( dofs_per_cell );
    std::vector< types::global_dof_index > local_dof_indices( dofs_per_cell );
    typename DoFHandler< dim >::active_cell_iterator cell = dof_handler.begin_active(),
                                                     endc = dof_handler.end();
    for ( ; cell != endc; ++cell )
    {
        fe_values.reinit( cell );
        cell_rhs = 0;
        for ( unsigned int q_index = 0; q_index < n_q_points; ++q_index )
            for ( unsigned int i = 0; i < dofs_per_cell; ++i )
            {
                cell_rhs( i ) += ( fe_values.shape_value( i, q_index ) *
                                   right_hand_side.value( fe_values.quadrature_point( q_index ) ) *
                                   fe_values.JxW( q_index ) );
            }
        cell->get_dof_indices( local_dof_indices );
        for ( unsigned int i = 0; i < dofs_per_cell; ++i )
        {
            v_dealii( local_dof_indices[ i ] ) += cell_rhs( i );
        }
    }

    return v;
}

template < class C1Operator >
Spacy::LinearSolver createCGSolver( const C1Operator& A, const Spacy::Vector& x )
{
    const auto& X = A.domain();
    const auto linearization = A.linearization( x );
    auto cg =
        Spacy::dealII::CGSolver< dim, VariableDims >( get( linearization ), X, X.dualSpace() );
    cg.setAbsoluteAccuracy( 1e-12 );
    return cg;
}

auto createOperator( const Spacy::VectorSpace& X )
{
    const double c = 1, d = 0;
    const double u0 = 0;
    dealii::Tensor< variable_dim, dim > du0;
    auto integrand = FunG::heatModel( c, d, u0, du0 );
    return Spacy::dealII::C1FunGOperator< decltype( integrand ), dim, VariableDims >(
        std::move( integrand ), X, X.dualSpace() );
}

int main()
{
    // dealii space setup
    dealii::Triangulation< dim > triangulation;
    dealii::GridGenerator::hyper_cube( triangulation, -1, 1 );
    triangulation.refine_global( 2 );
    const auto fe_order = 1;
    const auto boundaryConditions =
        Spacy::dealII::BoundaryPart< dim >( 0, std::make_shared< BoundaryValues< dim > >() );

    // Spacy
    const auto X = Spacy::dealII::makeHilbertSpace< variable_dim >( triangulation, fe_order, 0,
                                                                    {boundaryConditions} );
    const auto A = createOperator( X );
    auto cg = createCGSolver( A, zero( X ) );
    auto w = cg( getRhs( X ) );

    const auto& wd = Spacy::cast_ref< Spacy::dealII::Vector >( w );
    Spacy::dealII::writeVTK< dim, VariableDims >( wd, "sol" );

    return 0;
}
