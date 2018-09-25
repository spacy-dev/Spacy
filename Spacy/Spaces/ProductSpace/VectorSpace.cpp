#include "VectorSpace.h"

#include <Spacy/Spaces/ScalarSpace/Real.h>
#include <Spacy/Util/cast.hh>
#include <Spacy/vectorSpace.hh>
#include <Spacy/zeroVectorCreator.hh>

#include "ScalarProduct.h"
#include "Vector.h"

#include <algorithm>
#include <cassert>

namespace Spacy
{
    namespace ProductSpace
    {
        namespace
        {
            std::vector< std::shared_ptr< VectorSpace > >
            extractSubSpaces( const std::vector< std::shared_ptr< VectorSpace > >& spaces,
                              const std::vector< unsigned >& subSpaceIds )
            {
                std::vector< std::shared_ptr< VectorSpace > > subSpaces;
                for ( unsigned i : subSpaceIds )
                    subSpaces.push_back( spaces[ i ] );
                return subSpaces;
            }

            template < class T >
            std::vector< unsigned > extractSpaceIds( const std::vector< T >& v )
            {
                std::vector< unsigned > result( v.size() );

                for ( auto i = 0u; i < result.size(); ++i )
                    result[ i ] = i;
                return result;
            }

            std::vector< unsigned > defaultGlobalIds( unsigned n )
            {
                std::vector< unsigned > globalIds( n, 0 );
                for ( auto i = 1u; i < n; ++i )
                    globalIds[ i ] = i;

                return globalIds;
            }

            auto createMap( const std::vector< unsigned >& globalIds )
            {
                std::map< unsigned, unsigned > map;
                for ( auto i = 0u; i < globalIds.size(); ++i )
                    map[ globalIds[ i ] ] = i;

                return map;
            }
        }

        VectorCreator::VectorCreator( const std::vector< std::shared_ptr< VectorSpace > >& spaces )
            : spaces_( spaces ), idMap_( createMap( defaultGlobalIds( spaces.size() ) ) )
        {
        }

        VectorCreator::VectorCreator( const std::vector< std::shared_ptr< VectorSpace > >& spaces,
                                      const std::vector< unsigned >& subSpaceIds )
            : spaces_( spaces ), idMap_( createMap( subSpaceIds ) )
        {
        }

        const std::vector< std::shared_ptr< VectorSpace > >& VectorCreator::subSpaces() const
        {
            return spaces_;
        }

        VectorSpace& VectorCreator::subSpace( unsigned k )
        {
            return *spaces_[ k ];
        }

        const VectorSpace& VectorCreator::subSpace( unsigned k ) const
        {
            return *spaces_[ k ];
        }

        unsigned VectorCreator::idMap( unsigned k ) const
        {
            return idMap_.find( k )->second;
        }

        unsigned VectorCreator::inverseIdMap( unsigned k ) const
        {
            auto found = std::find_if( begin( idMap_ ), end( idMap_ ),
                                       [k]( const auto& entry ) { return entry.second == k; } );
            assert( found != end( idMap_ ) );
            return found->first;
        }

        bool VectorCreator::hasId( unsigned k ) const
        {
            return idMap_.find( k ) != end( idMap_ );
        }

        Vector VectorCreator::operator()( const VectorSpace* space ) const
        {
            return Vector{*space};
        }

        VectorSpace makeHilbertSpace( const std::vector< std::shared_ptr< VectorSpace > >& spaces )
        {
            return ::Spacy::makeHilbertSpace( VectorCreator{spaces}, ScalarProduct{} );
        }

        VectorSpace makeHilbertSpace( const std::vector< std::shared_ptr< VectorSpace > >& spaces,
                                      const std::vector< unsigned >& globalIds )
        {
            return ::Spacy::makeHilbertSpace(
                VectorCreator{extractSubSpaces( spaces, globalIds ), globalIds}, ScalarProduct{} );
        }

        VectorSpace makeHilbertSpace( const std::vector< std::shared_ptr< VectorSpace > >& spaces,
                                      const std::vector< unsigned >& primalSubSpaceIds,
                                      const std::vector< unsigned >& dualSubSpaceIds )
        {
            return ::Spacy::makeHilbertSpace(
                VectorCreator( {std::make_shared< VectorSpace >(
                                    ProductSpace::makeHilbertSpace( spaces, primalSubSpaceIds ) ),
                                std::make_shared< VectorSpace >(
                                    ProductSpace::makeHilbertSpace( spaces, dualSubSpaceIds ) )} ),
                ScalarProduct{} );
        }
    }

    std::shared_ptr< VectorSpace > extractSubSpaceImpl( const VectorSpace& space,
                                                        unsigned global_id )
    {
        if ( is< ProductSpace::VectorCreator >( space.creator() ) )
        {
            const auto& product_space_creator =
                cast_ref< ProductSpace::VectorCreator >( space.creator() );
            const auto& subSpaces = product_space_creator.subSpaces();
            std::shared_ptr< VectorSpace > result;
            for ( auto i = 0u; i < subSpaces.size(); ++i )
            {
                if ( product_space_creator.inverseIdMap( i ) == global_id &&
                     !is< ProductSpace::VectorCreator >( subSpaces[ i ]->creator() ) )
                    return subSpaces[ i ];

                if ( is< ProductSpace::VectorCreator >( subSpaces[ i ]->creator() ) )
                {
                    auto result = extractSubSpaceImpl( *subSpaces[ i ], global_id );
                    if ( result )
                        return result;
                }
            }
        }

        return nullptr;
    }

    const VectorSpace& extractSubSpace( const VectorSpace& space, unsigned global_id )
    {
        auto sub_space = extractSubSpaceImpl( space, global_id );
        if ( sub_space )
            return *sub_space;

        return space;
    }
}
