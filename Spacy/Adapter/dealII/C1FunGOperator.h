#pragma once

#include <deal.II/base/tensor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <Spacy/Util/Base/OperatorBase.h>
#include <Spacy/Util/Cast.h>
#include <Spacy/Util/Exceptions.h>
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
#include <iostream>
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
                  operator_space_( std::make_shared< VectorSpace >(
                      LinearOperatorCreator(), []( const Spacy::Vector& ) { return Real( 0 ); },
                      true ) ),
                  fe_system_( std::make_shared< dealii::FESystem< dim > >(
                      get_finite_element_system< dim, VariableDims >( domain ) ) ),
                  dof_handler_( std::make_shared< dealii::DoFHandler< dim > >(
                      creator< VectorCreator< dim, GetDim< 0, VariableDims >::value > >(
                          extractSubSpace( domain, 0 ) )
                          .triangulation() ) ),
                  value_( VariableDims::size == 1 ? GetDim< 0, VariableDims >::value
                                                  : VariableDims::size )
            {
                dof_handler_->distribute_dofs( *fe_system_ );
                InitBlockVector< dim, VariableDims >::apply( domain, value_ );
                InitSparsityPattern< dim, VariableDims >::apply( domain, *dof_handler_,
                                                                 sparsity_pattern_ );
                gradient_.reinit( sparsity_pattern_ );
            }

            C1FunGOperator( const C1FunGOperator& other )
                : OperatorBase( other.domain(), other.range() ), operator_( other.operator_ ),
                  operator_space_( other.operator_space_ ), fe_system_( other.fe_system_ ),
                  dof_handler_( other.dof_handler_ ),
                  value_( VariableDims::size == 1 ? GetDim< 0, VariableDims >::value
                                                  : VariableDims::size )
            {
                InitBlockVector< dim, VariableDims >::apply( domain(), value_ );
                InitSparsityPattern< dim, VariableDims >::apply( domain(), *dof_handler_,
                                                                 sparsity_pattern_ );
                gradient_.reinit( sparsity_pattern_ );
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
                return LinearOperator< dim, VariableDims >{gradient_, *operator_space_, domain(),
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

                Detail::OperatorUpdate< dim, VariableDims > operator_update( domain(),
                                                                             *fe_system_ );
                std::for_each( dof_handler_->begin_active(), dof_handler_->end(),
                               [this, &x_, &operator_update]( const auto& cell ) {
                                   operator_update.reinit( cell );
                                   operator_update.init_old_solution( x_ );

                                   for ( auto q_index = 0u;
                                         q_index < operator_update.impl_.n_q_points_; ++q_index )
                                   {
                                       operator_update.update_functional( operator_, q_index );
                                       operator_update.update_gradient_and_hessian( operator_,
                                                                                    q_index );
                                   }

                                   operator_update.transfer_local_data( value_, gradient_ );
                               } );

                oldX_ = x;
            }

            mutable FunGOperator operator_;
            std::shared_ptr< VectorSpace > operator_space_;
            std::shared_ptr< dealii::FESystem< dim > > fe_system_;
            std::shared_ptr< dealii::DoFHandler< dim > > dof_handler_;
            SparsityPattern sparsity_pattern_;

            mutable Value value_;
            mutable Gradient gradient_;

            mutable ::Spacy::Vector oldX_;
        };
    }
    /** @} */
}
