#include "vectorCreator.hh"

#include <Spacy/Spaces/ProductSpace/VectorSpace.h>
#include <Spacy/zeroVectorCreator.hh>

#include "scalarProduct.hh"

#include <memory>
#include <vector>

namespace Spacy
{
    namespace Rn
    {
        VectorSpace makeHilbertSpace( unsigned dim )
        {
            return ::Spacy::makeHilbertSpace( VectorCreator{dim}, EuclideanScalarProduct{} );
        }

        VectorSpace makeHilbertSpace( const std::vector< unsigned >& dims )
        {
            std::vector< std::shared_ptr< VectorSpace > > spaces;
            for ( unsigned dim : dims )
                spaces.push_back( std::make_shared< VectorSpace >( makeHilbertSpace( dim ) ) );
            return ProductSpace::makeHilbertSpace( spaces );
        }

        VectorSpace makeHilbertSpace( const std::vector< unsigned >& dims,
                                      const std::vector< unsigned >& primalIds,
                                      const std::vector< unsigned >& dualIds )
        {
            std::vector< std::shared_ptr< VectorSpace > > spaces;
            for ( unsigned dim : dims )
                spaces.push_back( std::make_shared< VectorSpace >( makeHilbertSpace( dim ) ) );
            return ProductSpace::makeHilbertSpace( spaces, primalIds, dualIds );
        }
    }
}
