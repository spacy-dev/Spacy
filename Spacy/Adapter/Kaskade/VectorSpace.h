#pragma once

#include <type_traits>

#include <Spacy/Spaces/ProductSpace.h>
#include <Spacy/Util/Mixins/Get.h>
#include <Spacy/Vector.h>
#include <Spacy/VectorSpace.h>

#include "Copy.h"
#include "L2Product.h"
#include "Vector.h"

namespace Spacy
{
    /** @addtogroup KaskadeGroup @{ */
    namespace Kaskade
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
            template < class... Args, class = std::enable_if_t<
            std::is_constructible< Description, Args... >::value > >
            VectorCreator( Args&&... args )
            : Mixin::Get< Description >( std::forward< Args >( args )... )
            {
            }
            
            /// Generate vector for %Kaskade 7.
            Vector< Description > operator()( const VectorSpace* space ) const
            {
                return Vector< Description >{*space};
            }
        };
        
        /**
         * @ingroup VectorSpaceGroup
         * @brief Create single space with hilbert space structure for %Kaskade 7.
         * @param space single %Kaskade 7 function space (no product space)
         */
        template < class Description >
        auto makeHilbertSpace( const Description& description, std::string name = "unnamed" )
        {
            return ::Spacy::makeHilbertSpace( Kaskade::VectorCreator< Description >{description},
                                              l2Product< Description >{}, name );
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
        auto makeHilbertSpace( const Description& description,
                               const std::vector< unsigned >& primalIds,
                               const std::vector< unsigned >& dualIds = {} )
        {
            constexpr int n =
            boost::fusion::result_of::size< typename Description::Variables >::value;
            std::vector< std::shared_ptr< VectorSpace > > newSpaces( n );
            
            std::cout << "create space with " << n << " variables." << std::endl;
            Detail::MakeSpaces< Description, 0, n >::apply( description, newSpaces );
            
            if ( dualIds.size() > 0 )
                return ::Spacy::ProductSpace::makeHilbertSpace( newSpaces, primalIds, dualIds );
            
            return ::Spacy::ProductSpace::makeHilbertSpace( newSpaces, primalIds );
        }
        
        
        template< class DomainDescription, class RangeDescription, int mStart, int ... m>
        class MultEmbedding : public  SubSpaceRelation 
        {
            using DomainVector = Vector<DomainDescription>;
            using RangeVector = Vector<RangeDescription>;
        public:
            
            MultEmbedding( const ::Spacy::VectorSpace& domain, const ::Spacy::VectorSpace& range ) :
            SubSpaceRelation(domain,range){ mapsOntoComponents.push_back(mStart); pushOthers<1,m...>(); };
            
            template <int count,int m0,int  ... M>
            void pushOthers()
            { 
                mapsOntoComponents.push_back(m0);
                pushOthers<count+1,M...>();
            }
            
            template <int count>
            void pushOthers()
            {
            }
            
            virtual void apply( const ::Spacy::Vector& xd, ::Spacy::Vector& xr) const
            {
                //std::cout << "Emb:" << 0 << "->" << mStart << std::endl;    
                ::Kaskade::component<mStart>(cast_ref<RangeVector>(xr).get())=::Kaskade::component<0>(cast_ref<DomainVector>(xd).get());
                applyV<1,m...>(xd,xr);
            }
            
            template<int subm,int m0, int ... M>
            void applyV( const ::Spacy::Vector& xd, ::Spacy::Vector& xr) const
            {
                //std::cout << "Emb:" << subm << "->" << m0 << std::endl;    
                ::Kaskade::component<m0>(cast_ref<RangeVector>(xr).get())=::Kaskade::component<subm>(cast_ref<DomainVector>(xd).get());
                applyV<subm+1,M...>(xd,xr);
            }
            
            template<int subm>
            void applyV( const ::Spacy::Vector& xd, ::Spacy::Vector& xr) const
            {
            }
            
            int getCBegin() const 
            {
                return  *std::min_element(mapsOntoComponents.begin(),mapsOntoComponents.end());
            }
            
            int getCEnd() const
            {
                return  *std::max_element(mapsOntoComponents.begin(),mapsOntoComponents.end())+1;
            }

            bool isInRange(int i) const
            {
               for(int k=0; k<mapsOntoComponents.size();++k) 
                   if(mapsOntoComponents[k]==i) return true;
               return false;
            }
        private:
            mutable std::vector<int> mapsOntoComponents;
        };
        
        
        template< class DomainDescription, class RangeDescription, int mStart, int ... m>
        class MultProjection : public  SubSpaceRelation 
        {
            using DomainVector = Vector<DomainDescription>;
            using RangeVector = Vector<RangeDescription>;
        public:
            
            MultProjection( const ::Spacy::VectorSpace& domain, const ::Spacy::VectorSpace& range ) :
            SubSpaceRelation(domain,range){  };
            
            
            virtual void apply( const ::Spacy::Vector& xd, ::Spacy::Vector& xr) const
            {
                //std::cout << "Proj:" << mStart << "->" << 0 << std::endl;    
                ::Kaskade::component<0>(cast_ref<RangeVector>(xr).get())=::Kaskade::component<mStart>(cast_ref<DomainVector>(xd).get());
                applyV<1,m...>(xd,xr);
            }
            
            template<int subm,int m0, int ... M>
            void applyV( const ::Spacy::Vector& xd, ::Spacy::Vector& xr) const
            {
                //std::cout << "Proj:" << m0 << "->" << subm << std::endl;    
                ::Kaskade::component<subm>(cast_ref<RangeVector>(xr).get())=::Kaskade::component<m0>(cast_ref<DomainVector>(xd).get());
                applyV<subm+1,M...>(xd,xr);
            }
            
            template<int subm>
            void applyV( const ::Spacy::Vector& xd, ::Spacy::Vector& xr) const
            {
            }
            
        };
        
        template < class DomainDescription, class RangeDescription,int ... m>
        auto makeProjection(const ::Spacy::VectorSpace& domain, const ::Spacy::VectorSpace& range )
        {
            return std::make_unique<MultProjection< DomainDescription, RangeDescription, m... > >(MultProjection< DomainDescription, RangeDescription, m... >(domain, range));
        }
        
        template < class DomainDescription, class RangeDescription,int ... m>
        auto makeEmbedding(const ::Spacy::VectorSpace& domain, const ::Spacy::VectorSpace& range )
        {
            return std::make_unique<MultEmbedding< DomainDescription, RangeDescription, m... > >(MultEmbedding< DomainDescription, RangeDescription, m... >(domain, range));
        }
        
        template<class Description, int m>
        struct SubSpace
        {
            using type=Detail::ExtractDescriptionAndSpace_t<Description,m>;
        };
        
        template<class Description, int m>
        ::Spacy::VectorSpace makeSingleEmbeddedSubspace(const ::Spacy::VectorSpace& total, const Description& totalDescription, std::string name="unamed")
        {
            using SubSpaceDescription=Detail::ExtractDescriptionAndSpace_t<Description,m>;
            auto spc=::boost::fusion::at_c<Description::template SpaceIndex<m>::spaceIndex>(totalDescription.spaces);
            typename SubSpaceDescription::Spaces subSpaceK(spc);
            SubSpaceDescription subSpaceDescription(subSpaceK,{name});
            auto subSpace=::Spacy::makeHilbertSpace(Kaskade::VectorCreator< SubSpaceDescription >{subSpaceDescription},
                                                    l2Product< SubSpaceDescription >{}, name );
            subSpace.setEmbedding(Spacy::Kaskade::makeEmbedding<SubSpaceDescription,Description,m>(subSpace,total));
            subSpace.setProjection(Spacy::Kaskade::makeProjection<Description,SubSpaceDescription,m>(total,subSpace));
            return subSpace;
        }
    }
    
    
    /** @} */
}
