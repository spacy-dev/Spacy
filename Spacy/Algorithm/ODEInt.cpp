#include "ODEInt.h"

#include <Spacy/ZeroVectorCreator.h>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <iterator>

namespace boost
{
    namespace numeric
    {
        namespace odeint
        {
            Spacy::Vector operator/( const Spacy::Vector& x, const Spacy::Vector& y )
            {
                auto z( x );
                std::transform( x.begin(), x.end(), y.begin(), z.begin(),
                                []( double x, double y ) { return x / y; } );
                return z;
            }

            Spacy::Vector abs( Spacy::Vector x )
            {
                std::transform( x.begin(), x.end(), x.begin(),
                                []( double xi ) { return std::abs( xi ); } );
                return x;
            }

            vector_space_norm_inf< Spacy::Vector >::result_type
            vector_space_norm_inf< Spacy::Vector >::operator()( const Spacy::Vector& x ) const
            {
                return std::accumulate( x.begin(), x.end(), result_type{0},
                                        []( result_type max, result_type xi ) {
                                            return std::max( max, std::abs( xi ) );
                                        } );
            }

            bool same_size_impl< Spacy::Vector, Spacy::Vector >::same_size( const Spacy::Vector& x,
                                                                            const Spacy::Vector& y )
            {
                if ( !bool( x ) && !bool( y ) )
                    return true;

                if ( bool( x ) ^ bool( y ) )
                    return false;

                return std::distance( x.begin(), x.end() ) == std::distance( y.begin(), y.end() );
            }

            void resize_impl< Spacy::Vector, Spacy::Vector >::resize( Spacy::Vector& x,
                                                                      const Spacy::Vector& y )
            {
                x = zero( y.space() );
            }
        }
    }
}

namespace Spacy
{
    SpacyIterator::SpacyIterator( Vector* p ) : iterator( p->begin() )
    {
    }

    void SpacyIterator::increment()
    {
        ++iterator;
    }

    void SpacyIterator::advance( std::ptrdiff_t n )
    {
        for ( std::ptrdiff_t i = 0; i < n; ++i )
            ++iterator;
    }

    bool SpacyIterator::equal( const SpacyIterator& other ) const
    {
        return iterator == other.iterator;
    }

    double& SpacyIterator::dereference() const
    {
        return *iterator;
    }

    ConstSpacyIterator::ConstSpacyIterator( const Vector* p ) : iterator( p->begin() )
    {
    }

    ConstSpacyIterator::ConstSpacyIterator( const SpacyIterator& p ) : iterator( p.iterator )
    {
    }

    void ConstSpacyIterator::increment()
    {
        ++iterator;
    }

    void ConstSpacyIterator::advance( std::ptrdiff_t n )
    {
        for ( std::ptrdiff_t i = 0; i < n; ++i )
            ++iterator;
    }

    bool ConstSpacyIterator::equal( const ConstSpacyIterator& other ) const
    {
        return iterator == other.iterator;
    }

    bool ConstSpacyIterator::equal( const SpacyIterator& other ) const
    {
        return iterator == other.iterator;
    }

    const double& ConstSpacyIterator::dereference() const
    {
        return *iterator;
    }

    SpacyIterator end_iterator( Vector* x )
    {
        SpacyIterator iter;
        iter.iterator = x->end();
        return iter;
    }

    ConstSpacyIterator end_iterator( const Vector* x )
    {
        ConstSpacyIterator iter;
        iter.iterator = x->end();
        return iter;
    }

    SpacyIterator range_begin( Vector& x )
    {
        return SpacyIterator( &x );
    }

    SpacyIterator range_end( Vector& x )
    {
        return end_iterator( &x );
    }

    ConstSpacyIterator range_begin( const Vector& x )
    {
        return ConstSpacyIterator( &x );
    }

    ConstSpacyIterator range_end( const Vector& x )
    {
        return end_iterator( &x );
    }
}
