#pragma once

#include <deal.II/base/tensor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>
// For boundary values
#include <deal.II/numerics/matrix_tools.h>

#include <Spacy/Util/Base/OperatorBase.h>
#include <Spacy/Util/Cast.h>
#include <Spacy/Util/Exceptions.h>
#include <Spacy/Util/Mixins/Get.h>
#include <Spacy/VectorSpace.h>
#include <Spacy/ZeroVectorCreator.h>

#include "C1FungOperatorAssembly.h"
#include "Copy.h"
#include "Init.h"
#include "LinearOperator.h"
#include "Util.h"
#include "Vector.h"
#include "VectorSpace.h"

#include <algorithm>
#include <memory>
#include <vector>

namespace Spacy
{
    /** @addtogroup dealIIGroup @{ */
    namespace dealII
    {
        template < class FunGOperator, int dim, class VariableDims = VariableDim< 1 > >
        class C1FunGOperator;

        template < class FunGOperator, int dim, int... variable_dims >
        class C1FunGOperator< FunGOperator, dim, VariableDim< variable_dims... > >
            : public OperatorBase
        {
            using VariableDims = VariableDim< variable_dims... >;
            using Value = typename Traits< VariableDims >::Vector;
            using Gradient = typename Traits< VariableDims >::Matrix;
            using SparsityPattern = typename Traits< VariableDims >::SparsityPattern;

        public:
            C1FunGOperator( FunGOperator&& operator_impl, const VectorSpace& domain,
                            const VectorSpace& range )
                : OperatorBase( domain, range ), operator_( std::move( operator_impl ) ),
                  operatorSpace_( std::make_shared< VectorSpace >(
                      LinearOperatorCreator(), []( const Spacy::Vector& ) { return Real( 0 ); },
                      true ) ),
                  value_( VariableDims::size == 1 ? GetDim< 0, VariableDims >::value
                                                  : VariableDims::size )
            {
                InitBlockVector< dim, VariableDims >::apply( domain, value_ );
                InitSparsityPattern< dim, VariableDims >::apply(
                    domain, Spacy::creator< VectorCreator< dim, VariableDims::value > >( domain )
                                .dofHandler(),
                    sparsityPattern_ );
                gradient_.reinit( sparsityPattern_ );
            }

            C1FunGOperator( const C1FunGOperator& other )
                : OperatorBase( other.domain(), other.range() ), operator_( other.operator_ ),
                  operatorSpace_( other.operatorSpace_ ),
                  value_( VariableDims::size == 1 ? GetDim< 0, VariableDims >::value
                                                  : VariableDims::size )
            {
                InitBlockVector< dim, VariableDims >::apply( domain(), value_ );
                InitSparsityPattern< dim, VariableDims >::apply(
                    domain(),
                    Spacy::creator< VectorCreator< dim, VariableDims::value > >( domain() )
                        .dofHandler(),
                    sparsityPattern_ );
                gradient_.reinit( sparsityPattern_ );
            }

            /// Compute \f$A(x)\f$.
            ::Spacy::Vector operator()( const ::Spacy::Vector& x ) const
            {
                assemble( x );

                auto y = zero( range() );
                copy( value_, y );
                return y;
            }

            /// Compute \f$A'(x)dx\f$.
            ::Spacy::Vector d1( const ::Spacy::Vector& x, const ::Spacy::Vector& dx ) const
            {
                assemble( x );

                auto dx_ = value_;
                copy( dx, dx_ );

                auto y_ = value_;
                gradient_.vmult( y_, dx_ );

                auto y = zero( range() );
                copy( y_, y );
                return y;
            }

            /**
             * @brief Access \f$A'(x)\f$ as linear operator \f$X\rightarrow Y\f$
             * @see LinearOperator, ::Spacy::LinearOperator
             */
            auto linearization( const ::Spacy::Vector& x ) const
            {
                assemble( x );
                return LinearOperator< dim, VariableDims >{gradient_, *operatorSpace_, domain(),
                                                           range()};
            }

        private:
            /// Assemble discrete representation of \f$A(x)\f$ and \f$A'(x)\f$.
            void assemble( const ::Spacy::Vector& x ) const
            {
                if ( oldX_ && ( oldX_ == x ) )
                    return;
                value_ = 0;
                gradient_ = 0;

                auto x_ = value_;
                copy( x, x_ );

                const auto& space =
                    Spacy::creator< VectorCreator< dim, VariableDims::value > >( domain() );
                Detail::OperatorUpdate< dim, VariableDims > operator_update( domain(),
                                                                             space.feSystem() );
                const auto cend = space.dofHandler().end();
                for ( auto cell = space.dofHandler().begin_active(); cell != cend; ++cell )
                {
                    operator_update.reinit( cell );
                    operator_update.init_old_solution( x_ );

                    for ( auto q_index = 0u; q_index < operator_update.impl_.n_q_points_;
                          ++q_index )
                    {
                        operator_update.update_functional( operator_, q_index );
                        operator_update.update_gradient_and_hessian( operator_, q_index );
                    }

                    operator_update.transfer_local_data( value_, gradient_ );
                }

                oldX_ = x;
            }

            mutable FunGOperator operator_;
            std::shared_ptr< VectorSpace > operatorSpace_;
            SparsityPattern sparsityPattern_;

        private:
            mutable Value value_;
            mutable Gradient gradient_;

            mutable ::Spacy::Vector oldX_;
        };
    }
    /** @} */
}
