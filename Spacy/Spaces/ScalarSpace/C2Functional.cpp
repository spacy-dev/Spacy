#include "C2Functional.h"

#include "LinearOperator.h"

#include <Spacy/LinearOperator.h>
#include <Spacy/Spaces/RealSpace.h>
#include <Spacy/Util/Cast.h>
#include <Spacy/Util/Exceptions.h>
#include <Spacy/Vector.h>
#include <Spacy/VectorSpace.h>
#include <Spacy/ZeroVectorCreator.h>

namespace Spacy
{
    namespace Scalar
    {
        C2Functional::C2Functional( std::function< double( double ) > value, std::function< double( double ) > derivative,
                                    std::function< double( double ) > secDerivative )
            : Spacy::FunctionalBase( Space::R ), value_( std::move( value ) ), derivative_( std::move( derivative ) ),
              secDerivative_( std::move( secDerivative ) ),
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

        Spacy::Real C2Functional::operator()( const ::Spacy::Vector& x ) const
        {
            return value_( get( cast_ref< Spacy::Real >( x ) ) );
        }

        Spacy::Vector C2Functional::d1( const ::Spacy::Vector& x ) const
        {
            return Spacy::Real( derivative_( get( cast_ref< Spacy::Real >( x ) ) ) );
        }

        Spacy::Vector C2Functional::d2( const ::Spacy::Vector& x, const ::Spacy::Vector& dx ) const
        {
            return Spacy::Real( secDerivative_( get( cast_ref< Spacy::Real >( x ) ) ) * get( cast_ref< Spacy::Real >( dx ) ) );
        }
        Spacy::Scalar::LinearOperator C2Functional::hessian( const ::Spacy::Vector& x ) const
        {
            return { *operatorSpace_, secDerivative_( get( cast_ref< Spacy::Real >( x ) ) ) };
        }
    } // namespace Scalar
} // namespace Spacy
