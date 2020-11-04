#include "C1Operator.h"

#include "LinearOperator.h"

#include <Spacy/Spaces/RealSpace.h>
#include <Spacy/Util/Cast.h>
#include <Spacy/Util/Exceptions.h>
#include <Spacy/Vector.h>
#include <Spacy/VectorSpace.h>
#include <Spacy/ZeroVectorCreator.h>

#include <cassert>
#include <utility>

namespace Spacy
{
    namespace Scalar
    {
        C1Operator::C1Operator( std::function< double( double ) > value, std::function< double( double ) > derivative )
            : C1Operator( std::move( value ), std::move( derivative ), Space::R, Space::R )
        {
        }

        C1Operator::C1Operator( std::function< double( double ) > value, std::function< double( double ) > derivative,
                                const VectorSpace& domain, const VectorSpace& range )
            : OperatorBase( domain, range ), value_( std::move( value ) ), derivative_( std::move( derivative ) ),
              operatorSpace_( std::make_shared< VectorSpace >(
                  []( const ::Spacy::VectorSpace* /*unused*/ ) -> Spacy::Vector {
                      throw Exception::CallOfUndefinedFunction( "OperatorCreator::operator()(const VectorSpace*)" );
                  },
                  []( const ::Spacy::Vector& /*unused*/ ) -> Spacy::Real {
                      throw Exception::CallOfUndefinedFunction( "LinearOperatorNorm" );
                  },
                  true ) )
        {
        }

        Spacy::Vector C1Operator::operator()( const ::Spacy::Vector& x ) const
        {
            assert( value_ );
            return Real( value_( get( cast_ref< Real >( x ) ) ), range() );
        }

        Spacy::Vector C1Operator::d1( const ::Spacy::Vector& x, const ::Spacy::Vector& dx ) const
        {
            assert( derivative_ );
            return Real( derivative_( get( cast_ref< Real >( x ) ) ), range() ) * get( cast_ref< Real >( dx ) );
        }

        LinearOperator C1Operator::linearization( const ::Spacy::Vector& x ) const
        {
            assert( derivative_ );
            assert( operatorSpace_ != nullptr );
            return { *operatorSpace_, derivative_( get( cast_ref< Real >( x ) ) ) };
        }
    } // namespace Scalar
} // namespace Spacy
