#include "Operator.h"

#include <Spacy/Spaces/RealSpace.h>
#include <Spacy/Util/Cast.h>
#include <Spacy/vector.hh>

#include <cassert>

namespace Spacy
{
    namespace Scalar
    {
        Operator::Operator( std::function< double( double ) > value )
            : OperatorBase( Space::R, Space::R ), value_( std::move( value ) )
        {
        }

        Operator::Operator( std::function< double( double ) > value, const VectorSpace& domain,
                            const VectorSpace& range )
            : OperatorBase( domain, range ), value_( std::move( value ) )
        {
        }

        Spacy::Vector Operator::operator()( const ::Spacy::Vector& x ) const
        {
            assert( value_ );
            return Real( value_( get( cast_ref< Real >( x ) ) ), range() );
        }
    }
}
