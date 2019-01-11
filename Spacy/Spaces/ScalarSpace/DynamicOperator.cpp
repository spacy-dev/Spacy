#include "DynamicOperator.h"

#include <Spacy/Spaces/RealSpace.h>
#include <Spacy/Util/Cast.h>
#include <Spacy/Vector.h>

#include <cassert>

namespace Spacy
{
    namespace Scalar
    {
        DynamicOperator::DynamicOperator( std::function< double( double, double ) > value )
            : OperatorBase( Space::R, Space::R ), value_( std::move( value ) )
        {
        }

        DynamicOperator::DynamicOperator( std::function< double( double, double ) > value, const VectorSpace& domain,
                            const VectorSpace& range )
            : OperatorBase( domain, range ), value_( std::move( value ) )
        {
        }

        Spacy::Vector DynamicOperator::operator()( double t,  const ::Spacy::Vector& x ) const
        {
            assert( value_ );
            return Real( value_( t, get( cast_ref< Real >( x ) ) ), range() );
        }
    }
}
