/**
  Author: Alexander Siegl - modified from Operator.cpp
  */
#include "Operator_V2.h"

#include <Spacy/Spaces/ProductSpace/Vector.h>
#include <Spacy/Spaces/ProductSpace/VectorSpace.h>
#include <Spacy/Util/Cast.h>
#include <Spacy/ZeroVectorCreator.h>

#include <numeric>

#include <iostream>

namespace Spacy
{
    namespace ProductSpace
    {
        Operator_V2::Operator_V2( const std::vector< std::vector< CallableOperator > >& blockOperators,
                            const Spacy::VectorSpace& domain, const Spacy::VectorSpace& range )
            : OperatorBase( domain, range ), blockOperators_( blockOperators )
        {
        }


        //Todo does not work properly
        Spacy::Vector Operator_V2::operator()( const ::Spacy::Vector& x ) const
        {
            auto y = zero( range() );
            const auto domain_is_product_space = is< Vector >( x );
            const auto range_is_product_space = is< Vector >( y );

          //  std::cout << "Test Operator: " << std::endl;
            if ( domain_is_product_space && !range_is_product_space )
            {
            //    std::cout << "Test Operator: 1" << std::endl;
                const auto& xp = cast_ref< Vector >( x );
                for ( auto j = 0u; j < xp.numberOfVariables(); ++j )
                    y += blockOperators_[ 0 ][ j ]( xp.component( j ) );

              //   std::cout << "Test Operator: 1" << std::endl;

                return y;
            }

            if ( !domain_is_product_space && range_is_product_space )
            {
               // std::cout << "Test Operator: 2" << std::endl;
                auto& yp = cast_ref< Vector >( y );
                for ( auto i = 0u; i < yp.numberOfVariables(); ++i )
                    yp.component( i ) = blockOperators_[ i ][ 0 ]( x );

                // std::cout << "Test Operator: 2" << std::endl;

                return y;
            }


            // ToDo does not work for general cases
            if ( domain_is_product_space && range_is_product_space )
            {
              //  std::cout << "Test Operator: 3" << std::endl;
                const auto& xp = cast_ref< Vector >( x );
                //std::cout << "Test Operator: 3" << std::endl;
                auto& yp = cast_ref< Vector >( y );
                //std::cout << "Test Operator: 3" << std::endl;

              //  std::cout << yp.numberOfVariables() << " " << xp.numberOfVariables() << std::endl;

               // std::cout << blockOperators_.size() << std::endl;

                for ( auto i = 0u; i < yp.numberOfVariables()+1; ++i )
                    for ( auto j = 0u; j < xp.numberOfVariables()+1; ++j )
                    {
                      //std::cout << i << " " << j << std::endl;
                      if( i < 2 && j < 2 )
                       (cast_ref< Vector >(yp.component( 0 ))).component(i) += blockOperators_[ i ][ j ]( cast_ref< Vector >(xp.component( 0 )).component(j) );

                      else if( i == 2 && j  == 2 )
                         (yp.component( 1 )) += blockOperators_[ i ][ j ]( xp.component( 1 ));


                      else if( i  == 2 )
                      {
                         // blockOperators_[ i ][ j ]( cast_ref< Vector >(xp.component( 0 )).component(j));
                       // std::cout << "Test: " << std::endl;
                        (cast_ref< Vector >(yp.component( 1 ))).component(0) += blockOperators_[ i ][ j ]( cast_ref< Vector >(xp.component( 0 )).component(j));
                       }
                      else if( j == 2 )
                        ((cast_ref< Vector >(yp.component( 0 ))).component(i)) += blockOperators_[ i ][ j ]( cast_ref< Vector >(xp.component( 1 )).component(0));

                      else
                          throw;


                    }
                //  std::cout << "Test Operator: 3" << std::endl;

                return y;
            }


          //  std::cout << "Test Operator: 4" << std::endl;

            return blockOperators_[ 0 ][ 0 ]( x );
        }

        CallableOperator& Operator_V2::operator()( size_type row, size_type col )
        {
            return blockOperators_[ row ][ col ];
        }

        const CallableOperator& Operator_V2::operator()( size_type row, size_type col ) const
        {
            return blockOperators_[ row ][ col ];
        }
    }
}
