#pragma once

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

    template < class T >
    struct ContiguousIterator
    {
        explicit ContiguousIterator( T* value = nullptr ) : value( value )
        {
        }

        ContiguousIterator operator++()
        {
            ++value;
            return *this;
        }

        ContiguousIterator operator++( int )
        {
            ContiguousIterator tmp( value );
            ++( *this );
            return tmp;
        }

        T& operator*() const
        {
            return *value;
        }

        bool operator!=( const ContiguousIterator& other ) const
        {
            return value != other.value;
        }

    private:
        T* value;
    };
}
