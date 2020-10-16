#pragma once

#include <Spacy/Spaces/ScalarSpace/Real.h>
#include <Spacy/Util/Base/AddArithmeticOperators.h>
#include <Spacy/Util/Base/VectorBase.h>
#include <Spacy/VectorSpace.h>

#include "Copy.h"

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

        template < class ValueRef >
        using Accessors = std::vector< std::pair< std::size_t, std::function< ValueRef( int ) > > >;

        template < int N, int I = 0 >
        struct AddAccessors
        {
            template < class Accessors, class Container >
            AddAccessors( Accessors& accessors, Container& container )
            {
                const auto size = boost::fusion::at_c< I >( container ).size();

// This is a bad hack and has to be revised. However, it compiles now                

                accessors.emplace_back( size, [&container]( int i )->double& {
                    return const_cast<double&>((boost::fusion::at_c< I >( container )[ i ])[0]);
                } );
/* Previous version:
                accessors.emplace_back( size, [&container]( int i ) {
                    return boost::fusion::at_c< I >( container )[ i ];
                } );*/


                AddAccessors< N, I + 1 >( accessors, container );
            }
        };

        template < int N >
        struct AddAccessors< N, N >
        {
            template < class Accessors, class Container >
            AddAccessors( Accessors&, Container& )
            {
            }
        };

        template < template < class... > class Container, class... Args >
        Accessors< double& > getAccessors( Container< Args... >& container )
        {
            Accessors< double& > accessors;
            AddAccessors< sizeof...( Args ) >( accessors, container );
            return accessors;
        }

        template < template < class... > class Container, class... Args >
        Accessors< const double& > getAccessors( const Container< Args... >& container )
        {
            Accessors< const double& > accessors;
            AddAccessors< sizeof...( Args ) >( accessors, container );
            return accessors;
        }

        template < class Vector >
        struct ContiguousIterator
        {
            static constexpr auto N = boost::fusion::result_of::size< Vector >::type::value;
            using ValueRef = typename std::conditional< std::is_const< Vector >::value,
                                                        const double&, double& >::type;
            using Accessor = Accessors< ValueRef >;

            explicit ContiguousIterator( Vector* v = nullptr, int i = 0 )
                : v( v ), i( i ), accessors( v ? getAccessors( *v ) : Accessor{} )
            {
            }

            ContiguousIterator operator++()
            {
                ++i;
                adjustIndices();
                return *this;
            }

            ContiguousIterator operator++( int )
            {
                ContiguousIterator tmp( v, i );
                ++( *this );
                return tmp;
            }

            ValueRef operator*() const
            {
                return accessors[ j ].second( i );
            }

            bool operator==( const ContiguousIterator& other ) const noexcept
            {
                return isAtEnd( other ) || ( v == other.v && i == other.i && j == other.j );
            }

        private:
            bool isAtEnd( const ContiguousIterator& other ) const noexcept
            {
                return ( other.v == nullptr && j == N ) || ( v == nullptr && other.j == N );
            }

            void adjustIndices()
            {
                if ( i == accessors[ j ].first )
                {
                    ++j;
                    i = 0;
                }
            }

            Vector* v;
            int i;
            Accessor accessors;
            std::size_t j{0};
        };

        /**
         * @ingroup KaskadeGroup VectorSpaceGroup
         * @brief Coefficient vector implementation for %Kaskade 7 (single space).
         * @tparam Description %Kaskade::VariableSetDescription
         */
        template < class Description >
        class Vector : public VectorBase, public AddArithmeticOperators< Vector< Description > >
        {
            using VectorImpl =
                typename Description::template CoefficientVectorRepresentation<>::type;
            using VariableSet = typename Description::VariableSet;
                
// The following is only useful, if the vector represents a single component Kaskade-VariableSet                
            using Variable =
                std::decay_t< std::remove_pointer_t< typename boost::fusion::result_of::value_at_c<
                    typename Description::Variables, 0 >::type > >;
            using Space =
                std::decay_t< std::remove_pointer_t< typename boost::fusion::result_of::value_at_c<
                    typename Description::Spaces, Variable::spaceIndex >::type > >;

        public:
            /**
             * @brief Construct vector \f$x\f$ from underlying vector space.
             * @param space underlying vector space
             */
            Vector( const VectorSpace& space )
                : VectorBase( space ),
                  description_( std::make_shared< Description >(
                      creator< VectorCreator< Description > >( space ).get() ) ),
                  v_( Description::template CoefficientVectorRepresentation<>::init(
                      description_->spaces ) )
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
            VectorImpl& get()
            {
                return v_;
            }

            /// Access coefficient vector.
            const VectorImpl& get() const
            {
                return v_;
            }

            const Description& getDescription()
            {
            return *description_;
            }
            
            /**
             * @brief Apply as dual element.
             * @param y primal vector
             * @return \f$x(y)\f$
             */
            Real operator()( const Vector& y ) const
            {
                return get() * y.get();
            }
            //      Vector& axpy(double a, const AbstractVector& y)
            //      {
            //        v_.axpy(a,castTo< Vector<Description> >(y).v_);
            //        return *this;
            //      }

            ContiguousIterator< typename VectorImpl::Sequence > begin()
            {
                return ContiguousIterator< typename VectorImpl::Sequence >{&v_.data};
            }

            ContiguousIterator< const typename VectorImpl::Sequence > begin() const
            {
                return ContiguousIterator< const typename VectorImpl::Sequence >{&v_.data};
            }

            ContiguousIterator< typename VectorImpl::Sequence > end()
            {
                return ContiguousIterator< typename VectorImpl::Sequence >{};
            }

            ContiguousIterator< const typename VectorImpl::Sequence > end() const
            {
                return ContiguousIterator< const typename VectorImpl::Sequence >{};
            }

        private:
            template < class Desc >
            friend void writeVTK( const Vector< Desc >& x, const std::string& fileName );
            std::shared_ptr< Description > description_;
            mutable VectorImpl v_;
        };
    }
}
