#pragma once

#include "Vector.h"

#include <Spacy/VectorSpace.h>

namespace Spacy
{
    namespace Generic
    {
        namespace Detail
        {
            template < class T >
            using TryStaticMemFn_Zero = decltype( T::Zero( 0 ) );

            template < class T, class = void >
            struct HasStaticMemFn_Zero : std::false_type
            {
            };

            template < class T >
            struct HasStaticMemFn_Zero< T, std::void_t< TryStaticMemFn_Zero< T > > > : std::is_constructible< T, TryStaticMemFn_Zero< T > >
            {
            };

            template < class T, std::enable_if_t< HasStaticMemFn_Zero< T >::value >* = nullptr >
            Vector< T > create_zero_vector( unsigned dim, const VectorSpace& space )
            {
                return Vector< T >( T::Zero( dim ), space );
            }

            template < class T, std::enable_if_t< !HasStaticMemFn_Zero< T >::value >* = nullptr >
            Vector< T > create_zero_vector( unsigned dim, const VectorSpace& space )
            {
                return Vector< T >( T( dim ), space );
            }
        } // namespace Detail

        /**
         * @ingroup GenericGroup
         * @brief Generic vector creator implementation
         */
        template < class VectorImpl >
        class VectorCreator
        {
        public:
            explicit VectorCreator( unsigned dim ) : dim_( dim )
            {
            }

            Vector< VectorImpl > operator()( const VectorSpace* space ) const
            {
                return Detail::create_zero_vector< VectorImpl >( dim_, *space );
            }

            void setDimension( unsigned dim )
            {
                dim_ = dim;
            }

            unsigned dim() const
            {
                return dim_;
            }

        private:
            unsigned dim_;
        };
    } // namespace Generic
} // namespace Spacy
