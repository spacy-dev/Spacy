#pragma once

#include <utility>

namespace Spacy
{
    /// @cond
    class VectorSpace;
    /// @endcond

    /**
     * @brief Base class for vector implementations.
     *
     * Provides access to the underlying vector space and related operations.
     * @see Kaskade::Vector, Fenics::Vector, Rn::Vector
     */
    class VectorBase
    {
    public:
        /**
         * @brief Constructor.
         * @param space underlying vector space
         */
        VectorBase( const VectorSpace& space );

        /// Copy constructor.
        VectorBase( const VectorBase& y );

        /// Move constructor.
        VectorBase( VectorBase&& y ) noexcept;

        /// Copy assignment.
        VectorBase& operator=( const VectorBase& y );

        /// Move assignment.
        VectorBase& operator=( VectorBase&& y ) noexcept;

        /// Access underlying vector space.
        const VectorSpace& space() const;

    private:
        const VectorSpace& space_;
    };

    /// Forward iterator over contiguous memory.
    template < class T >
    struct ContiguousIterator
    {
        /// Constructor.
        explicit ContiguousIterator( T* value = nullptr ) : value( value )
        {
        }

        /// Pre-increment.
        ContiguousIterator operator++()
        {
            ++value;
            return *this;
        }

        /// Post-increment.
        ContiguousIterator operator++( int )
        {
            ContiguousIterator tmp( value );
            ++( *this );
            return tmp;
        }

        /// Access to value.
        T& operator*() const
        {
            return *value;
        }

        /// Comparison.
        bool operator==( const ContiguousIterator& other ) const
        {
            return value == other.value;
        }

    private:
        T* value;
    };

    /// Forward iterator for vectors. Only use this if the vector implementation does not provide a
    /// forward iterator itself.
    template < class Vector, class Index = std::size_t >
    struct VectorIterator
    {
        using Value = decltype( std::declval< Vector >()[ std::declval< Index >() ] );

        /// Constructor.
        explicit VectorIterator( Vector* value = nullptr, Index i = Index{0} )
            : value( value ), i( i )
        {
        }

        /// Pre-increment.
        VectorIterator operator++()
        {
            ++i;
            return *this;
        }

        /// Post-increment.
        VectorIterator operator++( int )
        {
            VectorIterator tmp( value, i );
            ++( *this );
            return tmp;
        }

        /// Access to value.
        Value& operator*() const
        {
            return ( *value )[ i ];
        }

        /// Comparison.
        bool operator==( const VectorIterator& other ) const
        {
            return value == other.value && i == other.i;
        }

    private:
        Vector* value;
        Index i;
    };
}
