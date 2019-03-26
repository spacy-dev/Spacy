#pragma once

#include <fem/assemble.hh>
#include <fem/hierarchicErrorEstimator.hh>
#include <fem/hierarchicspace.hh>
#include <fem/lagrangespace.hh>
#include <linalg/jacobiPreconditioner.hh>

#include "Copy.h"
#include "Vector.h"

#include <Spacy/Util/Cast.h>
#include <Spacy/Vector.h>

#include <cmath>

namespace Spacy
{
    namespace Kaskade
    {
        template < class Functional, class GridManager >
        class HierarchicalErrorEstimator
        {
            static constexpr int numberOfAnsatzVariables = Functional::AnsatzVars::noOfVariables;
            static constexpr int spaceIndex = 0;
            static constexpr int extensionSpaceIndex = 1;
            static constexpr int dimension = GridManager::Grid::dimension;

            using Descriptions = typename Functional::AnsatzVars;
            using Variable = std::decay_t< typename boost::fusion::result_of::at_c<
                typename Descriptions::Variables, 0 >::type >;
            using ExtensionVariable =
                ::Kaskade::Variable<::Kaskade::SpaceIndex< extensionSpaceIndex >,
                                    ::Kaskade::Components< Variable::m >,
                                    ::Kaskade::VariableId< Variable::id > >;
            using ConstSpacePtr = std::remove_reference_t< typename boost::fusion::result_of::at_c<
                typename Descriptions::Spaces, Variable::spaceIndex >::type >;
            using ExtensionSpace =
                ::Kaskade::FEFunctionSpace<::Kaskade::ContinuousHierarchicExtensionMapper<
                    double, typename GridManager::Grid::LeafGridView > >;
            using ExtensionSpaces = boost::fusion::vector< ConstSpacePtr, ExtensionSpace const* >;
            using ExtensionVariableDescriptions = boost::fusion::vector< ExtensionVariable >;
            using ExtensionDescription =
                ::Kaskade::VariableSetDescription< ExtensionSpaces, ExtensionVariableDescriptions >;
            using CoefficientVectors =
                typename ExtensionDescription::template CoefficientVectorRepresentation<
                    0, numberOfAnsatzVariables >::type;
            using ErrorEstimator =
                ::Kaskade::HierarchicErrorEstimator<::Kaskade::LinearizationAt< Functional >,
                                                    ExtensionDescription >;
            using EstimatorAssembler = ::Kaskade::VariationalFunctionalAssembler< ErrorEstimator >;
            using EstimatorOperator = ::Kaskade::AssembledGalerkinOperator< EstimatorAssembler >;

            static constexpr int numberOfExtAnsatzVars = ErrorEstimator::AnsatzVars::noOfVariables;
            static constexpr int numberOfExtTestVars = ErrorEstimator::TestVars::noOfVariables;

        public:
            HierarchicalErrorEstimator( const Functional& F, GridManager& gridManager,
                                        const Descriptions& variableSetDescription,
                                        double minRefinementRatio, double tolerance )
                : F( F ), gridManager( gridManager ),
                  variableSetDescription( variableSetDescription ),
                  extensionSpace(
                      gridManager, gridManager.grid().leafGridView(),
                      boost::fusion::at_c< Variable::spaceIndex >( variableSetDescription.spaces )
                              ->mapper()
                              .getOrder() +
                          1 ),
                  spaces(
                      boost::fusion::at_c< Variable::spaceIndex >( variableSetDescription.spaces ),
                      &extensionSpace ),
                  minRefinementRatio( minRefinementRatio ), tolerance( tolerance )
            {
            }

            bool estimateError( const Spacy::Vector& x, const Spacy::Vector& dx )
            {
                CoefficientVectors estSol(
                    ExtensionDescription::template CoefficientVectorRepresentation< 0, 1 >::init(
                        spaces ) );
                const auto estRhs = getEstimatedRhsAndSol( x, dx, estSol );
                std::vector< double > estSolVec( estRhs.size() );
                estSol.write( estSolVec.begin() );

                if ( sqrt( innerProduct( estRhs, estSolVec ) ) < tolerance )
                {
                    std::cout << "||estim. error|| is smaller than requested" << std::endl;
                    return false;
                }

                updateErrorIndicator( estSol, norm( x + dx ).get() );
                return true;
            }

            std::vector< bool > getErrorIndicator() const
            {
                return errorIndicator;
            }

        private:
            double innerProduct( const std::vector< double >& v, const std::vector< double >& w )
            {
                assert( v.size() == w.size() );
                auto errNorm = 0.0;
                using Index = std::vector< double >::size_type;
                for ( Index k = 0; k < v.size(); ++k )
                    errNorm += v[ k ] * w[ k ];
                return errNorm;
            }

            std::vector< double > getEstimatedRhsAndSol( const Spacy::Vector& x_,
                                                         const Spacy::Vector& dx_,
                                                         CoefficientVectors& estSol )
            {
                typename Descriptions::VariableSet x( variableSetDescription ),
                    dx( variableSetDescription );
                Spacy::Kaskade::copy( x_, x );
                Spacy::Kaskade::copy( dx_, dx );

                EstimatorAssembler estAssembler{gridManager, spaces};
                estAssembler.assemble(
                    ErrorEstimator(::Kaskade::LinearizationAt< Functional >( F, x ), dx ) );

                const std::string varNames[ 2 ] = {"l", "e"};
                ExtensionDescription extensionDescription{spaces, varNames};
                size_t estSize = extensionDescription.degreesOfFreedom( 0, numberOfExtAnsatzVars );

                std::vector< double > estRhs( estSize );
                estAssembler.toSequence( 0, numberOfExtTestVars, estRhs.begin() );

                EstimatorOperator agro( estAssembler );
                CoefficientVectors rhs( estAssembler.rhs() );
                estSol = 1.0;
                ::Kaskade::JacobiPreconditioner< EstimatorOperator > jprec( agro, 1.0 );
                jprec.apply( estSol, rhs ); // single Jacobi iteration

                return estRhs;
            }

            double getErrorLevel( std::vector< double > errorDistribution, double maxErr,
                                  double norm_x )
            {
                std::sort( begin( errorDistribution ), end( errorDistribution ), std::greater<>() );
                const auto minRefineIndex = minRefinementRatio * ( errorDistribution.size() - 1 );
                const auto baseErrLevel =
                    std::max( errorDistribution[ minRefineIndex ],
                              norm_x * std::numeric_limits< double >::epsilon() );
                return std::min( baseErrLevel, 0.5 * maxErr );
            }

            void updateErrorIndicator( const CoefficientVectors& estSol, double norm_x )
            {
                std::vector< double > errorDistribution(
                    gridManager.grid().leafIndexSet().size( 0 ), 0.0 );
                std::transform(
                    gridManager.grid().leafGridView().template begin< 0 >(),
                    gridManager.grid().leafGridView().template end< 0 >(),
                    begin( errorDistribution ), [this, &estSol]( const auto& cell ) {
                        const auto indices = extensionSpace.mapper().globalIndices( cell );
                        return std::accumulate(
                            begin( indices ), end( indices ), 0.0,
                            [&estSol]( double err, const auto& index ) {
                                return err +=
                                       std::abs(::Kaskade::component< 0 >( estSol )[ index ] );
                            } );
                    } );

                const auto maxErr =
                    *std::max_element( begin( errorDistribution ), end( errorDistribution ) );
                const auto errLevel = getErrorLevel( errorDistribution, maxErr, norm_x );
                errorIndicator = std::vector< bool >( errorDistribution.size(), false );
                std::transform( begin( errorDistribution ), end( errorDistribution ),
                                begin( errorIndicator ),
                                [errLevel]( double error ) { return error > errLevel; } );
                assert( std::any_of( begin( errorIndicator ), end( errorIndicator ),
                                     []( bool value ) { return value; } ) );
            }

            const Functional& F;
            GridManager& gridManager;
            const Descriptions& variableSetDescription;

            ExtensionSpace extensionSpace;
            ExtensionSpaces spaces;

            std::vector< bool > errorIndicator{};
            double minRefinementRatio;
            double tolerance;
        };

        template < class Functional, class GridManager >
        HierarchicalErrorEstimator< Functional, GridManager > getHierarchicalErrorEstimator(
            const Functional& F, GridManager& gridManager,
            const typename Functional::AnsatzVars& variableSetDescription,
            double minRefinementRatio = 0.2, double tolerance = 1e-3 )
        {
            return HierarchicalErrorEstimator< Functional, GridManager >(
                F, gridManager, variableSetDescription, minRefinementRatio, tolerance );
        }
    }
}
