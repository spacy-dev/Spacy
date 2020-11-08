#pragma once

#include <Spacy/Vector.h>

#include <linalg/apcg.hh>
#include <linalg/symmetricOperators.hh>
#include <mg/multigrid.hh>

#include <dune/common/fmatrix.hh>

/** @addtogroup KaskadeGroup
 * @{
 */
namespace Spacy::Kaskade
{
    template < class VariableSetDescription, class Grid >
    class MultiGridSolver
    {
        static constexpr auto dim = Grid::dimension;
        using Matrix = NumaBCRSMatrix< Dune::FieldMatrix< double, dim, dim > >;
        using Spaces = typename VariableSetDescription::Spaces;
        using CoefficientVector = typename VariableSetDescription::template CoefficientVectorRepresentation<>::type;
        using X = Dune::BlockVector< Dune::FieldVector< double, dim > >;

    public:
        MultiGridSolver( Matrix A, Spaces spaces, int atol, int maxit, int verbose, bool onlyLowerTriangle = false,
                         ::Kaskade::DirectType directType = ::Kaskade::DirectType::UMFPACK )
            : A_{ std::move( A ) }, spaces_{ std::move( spaces ) }, atol_{ atol }, maxit_{ maxit }, verbose_{ verbose },
              onlyLowerTriangle_{ onlyLowerTriangle }, directType_{ directType }
        {
        }

        Spacy::Vector operator()( const Spacy::Vector& x ) const
        {
            CoefficientVector y_( VariableSetDescription::template CoefficientVectorRepresentation<>::init( spaces_ ) );
            CoefficientVector x_( VariableSetDescription::template CoefficientVectorRepresentation<>::init( spaces_ ) );
            Spacy::Kaskade::copy< VariableSetDescription >( x, x_ );

            ::Kaskade::DefaultDualPairing< X, X > dp{};
            auto mat = Dune::MatrixAdapter< Matrix, X, X >( A_ );
            ::Kaskade::SymmetricLinearOperatorWrapper< X, X > sa( mat, dp );
            ::Kaskade::PCGEnergyErrorTerminationCriterion< double > term( atol_, maxit_ );
            std::unique_ptr< ::Kaskade::SymmetricPreconditioner< X, X > > mg =
                ::Kaskade::makeMultigrid( A_, *boost::fusion::at_c< 0 >( spaces_ ), onlyLowerTriangle_, directType_ );
            ::Kaskade::Pcg< X, X > pcg( sa, *mg, term, verbose_ );
            Dune::InverseOperatorResult res{};
            pcg.apply( boost::fusion::at_c< 0 >( y_.data ), boost::fusion::at_c< 0 >( x_.data ), res );

            auto y = zero( x.space() );
            Spacy::Kaskade::copy< VariableSetDescription >( y_, y );
            return y;
        }

    private:
        Matrix A_;
        Spaces spaces_;
        int atol_;
        int maxit_;
        int verbose_;
        bool onlyLowerTriangle_;
        ::Kaskade::DirectType directType_;
    };
} // namespace Spacy::Kaskade
/** @} */
