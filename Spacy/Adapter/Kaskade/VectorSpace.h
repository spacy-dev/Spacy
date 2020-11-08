#pragma once

#include "Copy.h"
#include "L2Product.h"
#include "Vector.h"

#include <Spacy/Spaces/ProductSpace.h>
#include <Spacy/Util/Mixins/Get.h>
#include <Spacy/Vector.h>
#include <Spacy/VectorSpace.h>

#include <type_traits>

/** @addtogroup KaskadeGroup @{ */
namespace Spacy::Kaskade
{
    /// Creator for vector space elements for %Kaskade 7
    template < class Description >
    class VectorCreator : public Mixin::Get< Description >
    {
    public:
        /**
         * @ingroup VectorSpaceGroup
         * @brief Create from %Kaskade 7 function space.
         * @param space single %Kaskade 7 function space (no product space)
         */
        template < class... Args, class = std::enable_if_t< std::is_constructible< Description, Args... >::value > >
        VectorCreator( Args&&... args ) : Mixin::Get< Description >( std::forward< Args >( args )... )
        {
        }

        /// Generate vector for %Kaskade 7.
        Vector< Description > operator()( const VectorSpace* space ) const
        {
            return Vector< Description >{ *space };
        }
    };

    /**
     * @ingroup VectorSpaceGroup
     * @brief Create single space with hilbert space structure for %Kaskade 7.
     * @param space single %Kaskade 7 function space (no product space)
     */
    template < class Description >
    auto makeHilbertSpace( const Description& description, std::string name )
    {
        return ::Spacy::makeHilbertSpace( Kaskade::VectorCreator< Description >{ description }, l2Product< Description >{},
                                          std::move( name ) );
    }

    /**
     * @ingroup VectorSpaceGroup
     * @brief Create product space with hilbert space structure for %Kaskade 7.
     * @param spaces boost fusion forward sequence of const pointers to %Kaskade 7 function
     * spaces
     * @param primalIds ids of primal variables
     * @param dualIds ids of dual variables
     */
    template < class Description >
    auto makeHilbertSpace( const Description& description, const std::vector< unsigned >& primalIds,
                           const std::vector< unsigned >& dualIds = {} )
    {
        constexpr int n = boost::fusion::result_of::size< typename Description::Variables >::value;
        std::vector< std::shared_ptr< VectorSpace > > newSpaces( n );

        std::cout << "create space with " << n << " variables." << std::endl;
        Detail::MakeSpaces< Description, 0, n >::apply( description, newSpaces );

        if ( !dualIds.empty() )
            return ::Spacy::ProductSpace::makeHilbertSpace( newSpaces, primalIds, dualIds );

        return ::Spacy::ProductSpace::makeHilbertSpace( newSpaces, primalIds );
    }
} // namespace Spacy::Kaskade
/** @} */
