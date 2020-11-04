#include <iterator>

#include <Spacy/Adapter/Eigen/ScalarProduct.h>
#include <Spacy/Adapter/Eigen/Vector.h>
#include <Spacy/Adapter/Eigen/VectorCreator.h>
#include <Spacy/Adaptivity/SpaceManager.h>
#include <Spacy/Adaptivity/SpatialAdaptivity.h>
#include <Spacy/HilbertSpaceNorm.h>
#include <Spacy/Util/Cast.h>
#include <Spacy/Vector.h>
#include <Spacy/VectorSpace.h>
#include <Spacy/ZeroVectorCreator.h>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <vector>

using testing::Eq;

namespace Eigen
{
    auto begin( VectorXd& v )
    {
        return &v[ 0 ];
    }

    auto end( VectorXd& v ) -> decltype( &v[ 0 ] )
    {
        return nullptr;
    }
} // namespace Eigen

namespace
{
    class VectorCreatorForUnitInterval
    {
    public:
        explicit VectorCreatorForUnitInterval( unsigned dim ) : creator_( dim ), pointsInInterval_( dim )
        {
            assert( dim > 1 );
            for ( auto i = 0u; i < dim; ++i )
            {
                pointsInInterval_[ i ] = i / ( dim - 1 );
            }
        }

        Spacy::Rn::Vector operator()( const Spacy::VectorSpace* space ) const
        {
            return creator_( space );
        }

        void setDimension( unsigned dim )
        {
            creator_.setDimension( dim );
        }

        unsigned dim() const
        {
            return creator_.dim();
        }

        std::function< void( Spacy::Vector& ) > refineGrid( const Spacy::ErrorIndicator& indicator )
        {
            const auto oldPointsInInterval = pointsInInterval_;
            const auto indicators = boost::any_cast< std::vector< bool > >( indicator );

            for ( auto i = 0u; i < indicators.size(); ++i )
            {
                if ( indicators[ i ] )
                {
                    pointsInInterval_.push_back( 0.5 * ( pointsInInterval_[ i ] + pointsInInterval_[ i + 1 ] ) );
                }
            }

            using std::begin;
            using std::end;
            std::sort( begin( pointsInInterval_ ), end( pointsInInterval_ ) );
            setDimension( pointsInInterval_.size() );

            const auto newPointsInInterval = pointsInInterval_;
            return [ oldPointsInInterval, newPointsInInterval ]( Spacy::Vector& v ) mutable {
                auto& v_old = Spacy::cast_ref< Spacy::Rn::Vector >( v ).get();
                ::Eigen::VectorXd v_new( newPointsInInterval.size() );

                using std::begin;
                auto oldPoint = begin( oldPointsInInterval );
                auto oldValue = begin( v_old );
                auto newPoint = begin( newPointsInInterval );
                auto newValue = begin( v_new );

                auto equals = []( double a, double b ) {
                    const auto max = std::max( std::abs( a ), std::abs( b ) );
                    static const auto eps = std::numeric_limits< double >::epsilon();
                    if ( max < eps )
                        return true;
                    return std::abs( a - b ) < max * eps;
                };

                using std::end;
                while ( newPoint != end( newPointsInInterval ) )
                {
                    if ( equals( *newPoint, *oldPoint ) )
                    {
                        *newValue = *oldValue;
                        ++newPoint;
                        ++newValue;
                        continue;
                    }

                    if ( equals( *newPoint, *std::next( oldPoint ) ) )
                    {
                        ++oldPoint;
                        ++oldValue;
                        continue;
                    }

                    *newValue = 0.5 * ( *oldValue + *std::next( oldValue ) );
                    ++newPoint;
                    ++newValue;
                }

                get( Spacy::cast_ref< Spacy::Rn::Vector >( v ) ) = v_new;
            };
        }

        const std::vector< double >& getPoints() const
        {
            return pointsInInterval_;
        }

    private:
        Spacy::Rn::VectorCreator creator_;
        std::vector< double > pointsInInterval_;
    };
} // namespace

TEST( TestSpatialAdaptivity_1D, IF_NoAdaptivityIsDefinedForSpace_THEN_NoRefinementHappens )
{
    const auto initialDimension = 2;
    auto norm = Spacy::HilbertSpaceNorm( Spacy::Rn::EuclideanScalarProduct() );
    auto V = Spacy::VectorSpace( VectorCreatorForUnitInterval( initialDimension ), norm );
    ::Eigen::VectorXd v_e( initialDimension );
    Spacy::Vector v = Spacy::Rn::Vector( v_e, V );

    const auto errorIndicator = std::vector< bool >( { true } );
    EXPECT_THROW( Spacy::globalSpaceManager().adjustDiscretization( V.index(), errorIndicator ), std::runtime_error );

    auto& vAsEigen = get( Spacy::cast_ref< Spacy::Rn::Vector >( v ) );
    EXPECT_THAT( vAsEigen.size(), Eq( 2u ) );
}

TEST( TestSpatialAdaptivity_1D, TestRefinementFrom_2_To_3_Unknowns )
{
    const auto initialDimension = 2;
    auto norm = Spacy::HilbertSpaceNorm( Spacy::Rn::EuclideanScalarProduct() );
    auto V = Spacy::VectorSpace( VectorCreatorForUnitInterval( initialDimension ), norm );
    Spacy::globalSpaceManager().add( V.index(), [ &V ]( const Spacy::ErrorIndicator& indicator ) mutable {
        return V.creator().target< VectorCreatorForUnitInterval >()->refineGrid( indicator );
    } );
    ::Eigen::VectorXd v_e( initialDimension );
    v_e[ 0 ] = 1;
    v_e[ 1 ] = 2;
    Spacy::Vector v = Spacy::Rn::Vector( v_e, V );

    Spacy::globalSpaceManager().subscribe( &v );

    const auto errorIndicator = std::vector< bool >( { true } );
    Spacy::globalSpaceManager().adjustDiscretization( V.index(), errorIndicator );

    auto& vAsEigen = get( Spacy::cast_ref< Spacy::Rn::Vector >( v ) );
    EXPECT_THAT( vAsEigen.size(), Eq( 3u ) );
    EXPECT_THAT( vAsEigen[ 0 ], Eq( 1 ) );
    EXPECT_THAT( vAsEigen[ 1 ], Eq( 1.5 ) );
    EXPECT_THAT( vAsEigen[ 2 ], Eq( 2 ) );
}

TEST( TestSpatialAdaptivity_1D, TestRefinementFrom_2_To_3_To_4_Unknowns )
{
    const auto initialDimension = 2;
    auto norm = Spacy::HilbertSpaceNorm( Spacy::Rn::EuclideanScalarProduct() );
    auto V = Spacy::VectorSpace( VectorCreatorForUnitInterval( initialDimension ), norm );
    Spacy::globalSpaceManager().add( V.index(), [ &V ]( const Spacy::ErrorIndicator& indicator ) mutable {
        return V.creator().target< VectorCreatorForUnitInterval >()->refineGrid( indicator );
    } );
    ::Eigen::VectorXd v_e( initialDimension );
    v_e[ 0 ] = 1;
    v_e[ 1 ] = 2;
    Spacy::Vector v = Spacy::Rn::Vector( v_e, V );

    Spacy::globalSpaceManager().subscribe( &v );

    auto errorIndicator = std::vector< bool >( { true } );
    Spacy::globalSpaceManager().adjustDiscretization( V.index(), errorIndicator );

    auto& vAsEigen = get( Spacy::cast_ref< Spacy::Rn::Vector >( v ) );
    EXPECT_THAT( vAsEigen.size(), Eq( 3u ) );
    EXPECT_THAT( vAsEigen[ 0 ], Eq( 1 ) );
    EXPECT_THAT( vAsEigen[ 1 ], Eq( 1.5 ) );
    EXPECT_THAT( vAsEigen[ 2 ], Eq( 2 ) );

    errorIndicator = std::vector< bool >( { false, true } );
    Spacy::globalSpaceManager().adjustDiscretization( V.index(), errorIndicator );

    EXPECT_THAT( vAsEigen.size(), Eq( 4u ) );
    EXPECT_THAT( vAsEigen[ 0 ], Eq( 1 ) );
    EXPECT_THAT( vAsEigen[ 1 ], Eq( 1.5 ) );
    EXPECT_THAT( vAsEigen[ 2 ], Eq( 1.75 ) );
    EXPECT_THAT( vAsEigen[ 3 ], Eq( 2 ) );
}
