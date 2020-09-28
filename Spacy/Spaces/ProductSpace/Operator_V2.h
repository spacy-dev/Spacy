/**
  Author: Alexander Siegl - modified from Operator.hh
  */
#pragma once

#include <Spacy/Util/Base/OperatorBase.h>
#include <Spacy/Operator.h>
#include <Spacy/Util/Mixins.h>

#include <vector>

namespace Spacy
{
    namespace ProductSpace
    {
        /**
         * @brief Product space operator.
         *
         * The structure of this operator must be consistently defined through the (block-)operators
         * and
         * the domain and range space.
         */
        class Operator_V2 : public OperatorBase,
                public Mixin::AdjointIndex,
                public Mixin::ControlIndex,
                public Mixin::StateIndex
        {
        public:
            using size_type = std::vector< std::vector< ::Spacy::CallableOperator > >::size_type;

            /**
             * @brief Constructor.
             *
             * blockOperators and domain/range must be consistent, i.e. if range is a product space,
             * then the number of components in this space (ignoring subspaces) must coincide with
             * the
             * size of blockOperators. The same holds for domain wrt. each of blockOperators[i],
             * i...
             */
            Operator_V2( const std::vector< std::vector< CallableOperator > >& blockOperators,
                      const Spacy::VectorSpace& domain, const Spacy::VectorSpace& range );

            /// Compute \f$A(x)\in\mathbb{R}\f$.
            Spacy::Vector operator()( const Spacy::Vector& x ) const;

            /// Access \f$A_{row,col}\f$.
            CallableOperator& operator()( size_type row, size_type col );

            /// Access \f$A_{row,col}\f$ (read-only).
            const CallableOperator& operator()( size_type row, size_type col ) const;

        private:
            std::vector< std::vector< CallableOperator > > blockOperators_;
        };
    }
}
