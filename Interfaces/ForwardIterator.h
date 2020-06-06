#pragma once

#include <iterator>

namespace Spacy
{
    class ForwardIterator
    {
    public:
        using difference_type = std::size_t;
        using value_type = double;
        using reference = double&;
        using pointer = double*;
        using iterator_category = std::forward_iterator_tag;

        ForwardIterator& operator++();
        ForwardIterator operator++( int offset );
        double& operator*() const;
        bool operator==( const ForwardIterator& other ) const;
    };

    class ConstForwardIterator
    {
    public:
        using difference_type = std::size_t;
        using value_type = double;
        using reference = const double&;
        using pointer = const double*;
        using iterator_category = std::forward_iterator_tag;

        ConstForwardIterator& operator++();
        ConstForwardIterator operator++( int offset );
        const double& operator*() const;
        bool operator==( const ConstForwardIterator& other ) const;
    };

    inline bool operator!=( const ForwardIterator& lhs, const ForwardIterator& rhs )
    {
        return !( lhs == rhs );
    }

    inline bool operator!=( const ConstForwardIterator& lhs, const ConstForwardIterator& rhs )
    {
        return !( lhs == rhs );
    }
}
