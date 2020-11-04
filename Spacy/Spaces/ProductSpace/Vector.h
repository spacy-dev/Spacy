#pragma once

#include <Spacy/ForwardIterator.h>
#include <Spacy/Util/Base/VectorBase.h>
#include <Spacy/Util/Cast.h>
#include <Spacy/Util/Exceptions.h>
#include <Spacy/Vector.h>

#include <algorithm>
#include <memory>
#include <numeric>
#include <utility>
#include <vector>

namespace Spacy
{
    /** @addtogroup ProductSpaceGroup
     * @{
     */
    enum
    {
        PRIMAL = 0,
        DUAL = 1
    };

    namespace ProductSpace
    {
        /// @cond
        class VectorCreator;
        /// @endcond

        template < class T >
        struct Range
        {
            T begin;
            T end;
        };

        template < class T >
        bool operator==( const Range< T >& lhs, const Range< T >& rhs )
        {
            return lhs.begin == rhs.begin && lhs.end == rhs.end;
        }

        template < class T >
        bool operator!=( const Range< T >& lhs, const Range< T >& rhs )
        {
            return !( lhs == rhs );
        }

        template < class Iterator >
        struct ProductIterator
        {
            using Iterators = std::vector< Range< Iterator > >;
            using CurrentIterator = typename Iterators::iterator;
            using Reference = decltype( *std::declval< Iterator >() );

            explicit ProductIterator( Iterators iterators ) : iterators( std::move( iterators ) ), current( std::begin( this->iterators ) )
            {
            }

            ProductIterator( const ProductIterator& other )
                : iterators( other.iterators ), current( begin( iterators ) + distance( begin( other.iterators ), other.current ) )
            {
                if ( current == end( iterators ) )
                    return;
                const auto distToEnd = distance( other.current->begin, other.current->end );
                while ( distance( current->begin, current->end ) > distToEnd )
                    ++( *this );
            }

            ProductIterator& operator=( const ProductIterator& other )
            {
                iterators = other.iterators;
                current = begin( iterators ) + distance( begin( other.iterators ), other.current );
                if ( current == end( iterators ) )
                    return *this;

                const auto distToEnd = distance( other.current->begin, other.current->end );
                while ( distance( current->begin, current->end ) > distToEnd )
                    ++( *this );
                return *this;
            }

            ProductIterator& operator++()
            {
                if ( ++( current->begin ) == current->end )
                {
                    ++current;
                }
                return *this;
            }

            ProductIterator operator++( int )
            {
                ProductIterator tmp( *this );
                ++( *this );
                return tmp;
            }

            Reference operator*() const
            {
                return *( current->begin );
            }

            bool operator==( const ProductIterator& other ) const
            {
                if ( end( iterators ) == current && end( other.iterators ) == other.current )
                    return true;
                if ( distance( begin( iterators ), current ) != distance( begin( other.iterators ), other.current ) )
                    return false;
                if ( distance( current->begin, current->end ) != distance( other.current->begin, other.current->end ) )
                    return false;
                return true;
            }

            ProductIterator& makeEnd()
            {
                current = end( iterators );
                return *this;
            }

        private:
            template < class Iter1, class Iter2 >
            [[nodiscard]] std::ptrdiff_t distance( Iter1 begin, Iter2 end ) const
            {
                std::ptrdiff_t dist{ 0 };
                while ( begin != end )
                    ++dist, ++begin;
                return dist;
            }

            Iterators iterators;
            CurrentIterator current;
        };

        /**
         * @brief Product space vector.
         *
         * Represents a vector \f$x=(x_0,x_1,\ldots,x_n)\f$ of a product space \f$X =
         * \{X_0,X_1,\ldots,X_n\}\f$.
         */
        class Vector : public VectorBase
        {
        public:
            using iterator = typename std::vector< ::Spacy::Vector >::iterator;

            using const_iterator = typename std::vector< ::Spacy::Vector >::const_iterator;

            /**
             * @brief Construct product space vector.
             * @param space associated vector space
             */
            explicit Vector( const VectorSpace& space );

            /// @return \f$x+y\f$
            Vector& operator+=( const Vector& y );

            //    /// Axpy-operation \f$x = x + ay\f$.
            //    AbstractVector& axpy(double a, const AbstractVector& y);

            /// @return \f$x-y\f$
            Vector& operator-=( const Vector& y );

            /// @return \f$ax\f$
            Vector& operator*=( double a );

            /// @return \f$-x\f$
            Vector operator-() const;

            /// Equality comparison (possibly up to the maximal attainable accuracy).
            bool operator==( const Vector& y ) const;

            [[nodiscard]] unsigned numberOfVariables() const;

            /**
             * @brief Access k-th component.
             * @return associated vector \f$x_k\f$
             */
            ::Spacy::Vector& component( unsigned k );

            /**
             * @brief Access k-th component.
             * @return associated vector \f$x_k\f$
             */
            [[nodiscard]] const ::Spacy::Vector& component( unsigned k ) const;

            /**
             * @brief Access VectorCreator object.
             * @see ProductSpace::VectorCreator
             */
            [[nodiscard]] const VectorCreator& creator() const;

            /// Apply as dual element.
            Real operator()( const Vector& y ) const;

            iterator componentBegin();

            iterator componentEnd();

            [[nodiscard]] const_iterator componentBegin() const;

            [[nodiscard]] const_iterator componentEnd() const;

            ProductIterator< ForwardIterator > begin()
            {
                return ProductIterator< ForwardIterator >( getBegin() );
            }

            ProductIterator< ForwardIterator > end()
            {
                return ProductIterator< ForwardIterator >( getEnd() ).makeEnd();
            }

            [[nodiscard]] ProductIterator< ConstForwardIterator > begin() const
            {
                return ProductIterator< ConstForwardIterator >( getBegin() );
            }

            [[nodiscard]] ProductIterator< ConstForwardIterator > end() const
            {
                return ProductIterator< ConstForwardIterator >( getEnd() ).makeEnd();
            }

        private:
            std::vector< Range< ForwardIterator > > getBegin()
            {
                std::vector< Range< ForwardIterator > > iterators( components_.size() );
                std::transform( componentBegin(), componentEnd(), std::begin( iterators ), []( Spacy::Vector& component ) {
                    return Range< ForwardIterator >{ component.begin(), component.end() };
                } );
                return iterators;
            }

            [[nodiscard]] std::vector< Range< ConstForwardIterator > > getBegin() const
            {
                std::vector< Range< ConstForwardIterator > > iterators( components_.size() );
                std::transform( componentBegin(), componentEnd(), std::begin( iterators ), []( const Spacy::Vector& component ) {
                    return Range< ConstForwardIterator >{ component.begin(), component.end() };
                } );
                return iterators;
            }

            std::vector< Range< ForwardIterator > > getEnd()
            {
                std::vector< Range< ForwardIterator > > iterators( components_.size() );
                std::transform( componentBegin(), componentEnd(), std::begin( iterators ), []( Spacy::Vector& component ) {
                    return Range< ForwardIterator >{ component.end(), component.end() };
                } );
                return iterators;
            }

            [[nodiscard]] std::vector< Range< ConstForwardIterator > > getEnd() const
            {
                std::vector< Range< ConstForwardIterator > > iterators( components_.size() );
                std::transform( componentBegin(), componentEnd(), std::begin( iterators ), []( const Spacy::Vector& component ) {
                    return Range< ConstForwardIterator >{ component.end(), component.end() };
                } );
                return iterators;
            }

            std::vector< ::Spacy::Vector > components_ = {};
        };

        template < class T >
        struct BasicComponentView
        {
            explicit BasicComponentView( std::vector< T* > components ) : components_( std::move( components ) )
            {
            }

            T& operator[]( unsigned k )
            {
                return *components_[ k ];
            }

            const T& operator[]( unsigned k ) const
            {
                return *components_[ k ];
            }

        private:
            std::vector< T* > components_;
        };

        template < class T >
        struct BasicComponentView< const T >
        {
            explicit BasicComponentView( std::vector< const T* > components ) : components_( std::move( components ) )
            {
            }

            const T& operator[]( unsigned k ) const
            {
                return *components_[ k ];
            }

        private:
            std::vector< const T* > components_;
        };

        template < class SingleSpaceVector >
        void extractSingleSpaceVectors( ::Spacy::Vector& x, std::vector< SingleSpaceVector* >& components )
        {
            if ( is< SingleSpaceVector >( x ) )
            {
                components.push_back( &cast_ref< SingleSpaceVector >( x ) );
                return;
            }

            if ( is< ProductSpace::Vector >( x ) )
            {
                auto& x_ = cast_ref< ProductSpace::Vector >( x );
                std::for_each( x_.componentBegin(), x_.componentEnd(),
                               [ &components ]( auto& x ) { extractSingleSpaceVectors( x, components ); } );
                return;
            }

            throw Exception::InvalidArgument( "ProductSpace::Vector::extractSingleSpaceVectors" );
        }

        template < class SingleSpaceVector >
        void extractConstSingleSpaceVectors( const ::Spacy::Vector& x, std::vector< const SingleSpaceVector* >& components )
        {
            if ( is< SingleSpaceVector >( x ) )
            {
                components.push_back( &cast_ref< SingleSpaceVector >( x ) );
                return;
            }

            if ( is< ProductSpace::Vector >( x ) )
            {
                const auto& x_ = cast_ref< ProductSpace::Vector >( x );
                std::for_each( x_.componentBegin(), x_.componentEnd(),
                               [ &components ]( const auto& x ) { extractConstSingleSpaceVectors( x, components ); } );
                return;
            }

            throw Exception::InvalidArgument( "ProductSpace::Vector::extractConstSingleSpaceVectors" );
        }

        template < class SingleSpaceVector >
        BasicComponentView< SingleSpaceVector > extractSingleComponentView( ::Spacy::Vector& x )
        {
            std::vector< SingleSpaceVector* > components;
            extractSingleSpaceVectors( x, components );
            return BasicComponentView< SingleSpaceVector >( components );
        }

        template < class SingleSpaceVector >
        BasicComponentView< const SingleSpaceVector > extractSingleComponentView( const ::Spacy::Vector& x )
        {
            std::vector< const SingleSpaceVector* > components;
            extractConstSingleSpaceVectors( x, components );
            return BasicComponentView< const SingleSpaceVector >( components );
        }
    } // namespace ProductSpace

    /// @return cast_ref<ProductSpace::Vector>(v).component(PRIMAL);
    ::Spacy::Vector& primalComponent( ::Spacy::Vector& v );

    /// @return cast_ref<ProductSpace::Vector>(v).component(PRIMAL);
    const ::Spacy::Vector& primalComponent( const ::Spacy::Vector& v );

    /// @return cast_ref<ProductSpace::Vector>(v).component(DUAL);
    ::Spacy::Vector& dualComponent( ::Spacy::Vector& v );

    /// @return cast_ref<ProductSpace::Vector>(v).component(DUAL);
    const ::Spacy::Vector& dualComponent( const ::Spacy::Vector& v );

    /** @} */
} // namespace Spacy
