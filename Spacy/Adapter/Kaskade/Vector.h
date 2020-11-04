#pragma once

#include "Copy.h"

#include <Spacy/Spaces/ScalarSpace/Real.h>
#include <Spacy/Util/Base/AddArithmeticOperators.h>
#include <Spacy/Util/Base/VectorBase.h>
#include <Spacy/Vector.h>
#include <Spacy/VectorSpace.h>

#include <boost/fusion/include/size.hpp>
#include <boost/fusion/mpl/at.hpp>

#include <string>

namespace Spacy
{
    namespace Kaskade
    {
        /// \cond
        template < class >
        class VectorCreator;
        /// \endcond

        /// Range [begin, end)
        template < class T >
        struct Range
        {
            T begin;
            T end;
        };

        template < class T >
        bool operator<( const Range< T >& lhs, const Range< T >& rhs ) noexcept
        {
            return lhs.end <= rhs.begin;
        }

        template < class BlockVector >
        class AccessBlockVectorEntry
        {
        public:
            explicit AccessBlockVectorEntry( BlockVector& v, int i = 0, int k = 0 ) : v( v ), i( i ), k( k )
            {
            }

            double& get()
            {
                return v[ i ][ k ];
            }

            const double& get() const
            {
                return v[ i ][ k ];
            }

            void increment()
            {
                ++k;
                using std::size;
                if ( k == size( v[ i ] ) )
                {
                    k = 0;
                    ++i;
                }
            }

        private:
            BlockVector& v;
            int i;
            int k;
        };

        namespace Detail
        {
            template < class ReturnT >
            struct AccessInfo
            {
                int blockSize;
                int spatialDim;
                std::function< ReturnT( int, int ) > accessor;
            };

            template < class ReturnT >
            using AccessInfos = std::vector< AccessInfo< ReturnT > >;

            template < class ReturnT, int N, int I = 0 >
            struct AddAccessInfo
            {
                template < class AccessInfos, class Container >
                AddAccessInfo( AccessInfos& accessInfos, Container& container )
                {
                    const auto blockSize = int( boost::fusion::at_c< I >( container ).size() );
                    const auto spatialDim = int( boost::fusion::at_c< I >( container )[ 0 ].size() );
                    accessInfos.push_back( AccessInfo< ReturnT >{ blockSize, spatialDim, [ &container ]( int i, int k ) -> ReturnT {
                                                                     return boost::fusion::at_c< I >( container )[ i ][ k ];
                                                                 } } );
                    AddAccessInfo< ReturnT, N, I + 1 >( accessInfos, container );
                }
            };

            template < class ReturnT, int N >
            struct AddAccessInfo< ReturnT, N, N >
            {
                template < class AccessInfos, class Container >
                AddAccessInfo( AccessInfos& /*unused*/, Container& /*unused*/ )
                {
                }
            };

            template < template < class... > class Container, class... Args >
            auto getAccessInfos( Container< Args... >& container )
            {
                using ReturnT = double&;
                AccessInfos< ReturnT > accessInfos;
                AddAccessInfo< ReturnT, sizeof...( Args ) >( accessInfos, container );
                return accessInfos;
            }

            template < template < class... > class Container, class... Args >
            auto getAccessInfos( const Container< Args... >& container )
            {
                using ReturnT = const double&;
                AccessInfos< ReturnT > accessInfos;
                AddAccessInfo< ReturnT, sizeof...( Args ) >( accessInfos, container );
                return accessInfos;
            }
        } // namespace Detail

        template < class Vector >
        struct ContiguousIterator
        {
            static constexpr auto N = boost::fusion::result_of::size< Vector >::type::value;
            using ValueRef = typename std::conditional< std::is_const< Vector >::value, const double&, double& >::type;
            using AccessInfo = Detail::AccessInfos< ValueRef >;

            explicit ContiguousIterator( Vector* v = nullptr, int i = 0, int k = 0 )
                : v( v ), i( i ), k( k ), infos( v ? Detail::getAccessInfos( *v ) : AccessInfo{} )
            {
            }

            ContiguousIterator operator++()
            {
                increment();
                return *this;
            }

            ContiguousIterator operator++( int )
            {
                ContiguousIterator tmp( v, i, k );
                ++( *this );
                return tmp;
            }

            ValueRef operator*() const
            {
                return infos[ j ].accessor( i, k );
            }

            bool operator==( const ContiguousIterator& other ) const noexcept
            {
                return isAtEnd( other ) || ( v == other.v && i == other.i && j == other.j && k == other.k );
            }

        private:
            [[nodiscard]] bool isAtEnd( const ContiguousIterator& other ) const noexcept
            {
                return ( other.v == nullptr && j == N ) || ( v == nullptr && other.j == N );
            }

            void increment()
            {
                ++k;
                if ( k == infos[ j ].spatialDim )
                {
                    ++i;
                    k = 0;
                }

                if ( i == infos[ j ].blockSize )
                {
                    ++j;
                    i = 0;
                }
            }

            Vector* v;
            int i{ 0 };
            int k{ 0 };
            AccessInfo infos;
            std::size_t j{ 0 };
        };

        /**
         * @ingroup KaskadeGroup VectorSpaceGroup
         * @brief Coefficient vector implementation for %Kaskade 7 (single space).
         * @tparam Description %Kaskade::VariableSetDescription
         */
        template < class Description >
        class Vector : public VectorBase, public AddArithmeticOperators< Vector< Description > >
        {
            using VectorImpl = typename Description::template CoefficientVectorRepresentation<>::type;
            using Variable = std::decay_t<
                std::remove_pointer_t< typename boost::fusion::result_of::value_at_c< typename Description::Variables, 0 >::type > >;
            using Space = std::decay_t< std::remove_pointer_t<
                typename boost::fusion::result_of::value_at_c< typename Description::Spaces, Variable::spaceIndex >::type > >;
            using VariableSet = typename Description::VariableSet;

        public:
            /**
             * @brief Construct vector \f$x\f$ from underlying vector space.
             * @param space underlying vector space
             */
            Vector( const VectorSpace& space )
                : VectorBase( space ), variableSet_( creator< VectorCreator< Description > >( space ).get() ),
                  description_( std::make_shared< Description >( creator< VectorCreator< Description > >( space ).get() ) ),
                  v_( Description::template CoefficientVectorRepresentation<>::init( variableSet_.descriptions.spaces ) )
            {
            }

            //      /// Assign from coefficient vector of %Kaskade 7.
            //      Vector& operator=(const typename Description::template
            //      CoefficientVectorRepresentation<>::type& v)
            //      {
            //        variableSet_
            ////        v_ = v;
            //        return *this;
            //      }

            /// Access coefficient vector.
            VariableSet& get()
            {
                return variableSet_;
            }

            /// Access coefficient vector.
            const VariableSet& get() const
            {
                return variableSet_;
            }

            /**
             * @brief Apply as dual element.
             * @param y primal vector
             * @return \f$x(y)\f$
             */
            Real operator()( const Vector& y ) const
            {
                VectorImpl v( Description::template CoefficientVectorRepresentation<>::init( variableSet_.descriptions.spaces ) );
                VectorImpl w( Description::template CoefficientVectorRepresentation<>::init( variableSet_.descriptions.spaces ) );
                variableSetToCoefficients( get(), v );
                variableSetToCoefficients( y.get(), w );
                return v * w;
                //        return get() * y.get();
                //        auto v = VectorImpl(variableSet_), w = VectorImpl(variableSet_);
                //        return v*w;
            }
            //      Vector& axpy(double a, const AbstractVector& y)
            //      {
            //        v_.axpy(a,castTo< Vector<Description> >(y).v_);
            //        return *this;
            //      }

            ContiguousIterator< typename VectorImpl::Sequence > begin()
            {
                return ContiguousIterator< typename VectorImpl::Sequence >{ &v_.data };
            }

            ContiguousIterator< const typename VectorImpl::Sequence > begin() const
            {
                return ContiguousIterator< const typename VectorImpl::Sequence >{ &v_.data };
            }

            ContiguousIterator< typename VectorImpl::Sequence > end()
            {
                return ContiguousIterator< typename VectorImpl::Sequence >{};
            }

            ContiguousIterator< const typename VectorImpl::Sequence > end() const
            {
                return ContiguousIterator< const typename VectorImpl::Sequence >{};
            }

            std::shared_ptr< Description > description() const
            {
                return description_;
            }

        private:
            VariableSet variableSet_;
            std::shared_ptr< Description > description_;
            mutable VectorImpl v_;
        };
    } // namespace Kaskade
} // namespace Spacy
