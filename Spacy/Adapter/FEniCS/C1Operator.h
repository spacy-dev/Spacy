#pragma once

#include "AssignXIfPresent.h"
#include "Copy.h"
#include "LUSolver.h"
#include "LinearOperator.h"
#include "OperatorSpace.h"
#include "VectorSpace.h"
#include <dolfin.h>

#include <Spacy/Spaces/ScalarSpace/RealSpace.h>
#include <Spacy/Util/Base/OperatorBase.h>
#include <Spacy/Util/Cast.h>
#include <Spacy/VectorSpace.h>

#include <memory>

namespace Spacy
{
    /** @addtogroup FenicsGroup
     * @{
     */
    namespace FEniCS
    {
        /**
         * @brief %Operator interface for %FEniCS. Models a differentiable operator
         * \f$A:X\rightarrow Y\f$.
         * @warning In the .ufl file you have to name the argument of \f$f\f$ by "x"!
         */
        template < class ResidualForm, class JacobianForm >
        class C1Operator : public OperatorBase
        {
        public:
            /**
             * @brief Construct operator for %FEniCS.
             * @param F Residual form for the evaluation of \f$A\f$
             * @param J Jacobian form for the evaluation of \f$A'\f$
             * @param bcs Dirichlet boundary conditions
             * @param domain domain space \f$X\f$
             * @param range range space \f$Y\f$
             */
            C1Operator( const ResidualForm& F, const JacobianForm& J, const std::vector< const dolfin::DirichletBC* >& bcs,
                        const VectorSpace& domain, const VectorSpace& range )
                : OperatorBase( domain, range ), F_( F.function_space( 0 ) ), J_( J.function_space( 0 ), J.function_space( 1 ) ),
                  bcs_( bcs ),
                  operatorSpace_( std::make_shared< VectorSpace >(
                      LinearOperatorCreator( domain, range, J.function_space( 0 ) ),
                      []( const Spacy::Vector& v ) { return dolfin::Matrix( cast_ref< LinearOperator >( v ).get() ).norm( "frobenius" ); },
                      "operator space (fenics)", true ) )
            {
                copyCoefficients( F, F_ );
                copyCoefficients( J, J_ );
            }

            /**
             * @brief Construct operator without boundary conditions for %FEniCS.
             * @param F Residual form for the evaluation of \f$A\f$
             * @param J Jacobian form for the evaluation of \f$A'\f$
             * @param domain domain space \f$X\f$
             * @param range range space \f$Y\f$
             */
            C1Operator( const ResidualForm& F, const JacobianForm& J, const VectorSpace& domain, const VectorSpace& range )
                : C1Operator( F, J, {}, domain, range )
            {
            }

            C1Operator( C1Operator&& other )
                : OperatorBase( other ), F_( other.F_.function_space( 0 ) ),
                  J_( other.J_.function_space( 0 ), other.J_.function_space( 1 ) ), bcs_( std::move( other.bcs_ ) ),
                  A_( std::move( other.A_ ) ), b_( std::move( other.b_ ) ), operatorSpace_( std::move( other.operatorSpace_ ) )
            {
                copyCoefficients( other.F_, F_ );
                copyCoefficients( other.J_, J_ );
            }

            C1Operator( const C1Operator& other )
                : OperatorBase( other ), F_( other.F_.function_space( 0 ) ),
                  J_( other.J_.function_space( 0 ), other.J_.function_space( 1 ) ), bcs_( other.bcs_ ),
                  A_( ( other.A_ != nullptr ) ? /*other.A_->copy()*/ std::make_shared< dolfin::Matrix >( *other.A_ ) : nullptr ),
                  b_( ( other.b_ != nullptr ) ? other.b_->copy() : nullptr ), operatorSpace_( other.operatorSpace_ )
            {
                copyCoefficients( other.F_, F_ );
                copyCoefficients( other.J_, J_ );
            }

            C1Operator& operator=( const C1Operator& other )
            {
                OperatorBase::operator=( other );
                F_ = other.F_;
                J_ = other.J_;
                bcs_ = other.bcs_;
                A_ = ( other.A_ != nullptr ) ? /*other.A_->copy()*/ std::make_shared< dolfin::Matrix >( *other.A_ ) : nullptr;
                b_ = ( other.b_ != nullptr ) ? other.b_->copy() : nullptr;
                operatorSpace_ = other.operatorSpace_;

                copyCoefficients( other.F_, F_ );
                copyCoefficients( other.J_, J_ );
            }

            C1Operator& operator=( C1Operator&& other )
            {
                OperatorBase::operator=( other );
                F_ = other.F_;
                J_ = other.J_;
                bcs_ = std::move( other.bcs_ );
                A_ = std::move( other.A_ );
                b_ = std::move( other.b_ );
                operatorSpace_ = std::move( other.operatorSpace_ );

                copyCoefficients( other.F_, F_ );
                copyCoefficients( other.J_, J_ );
            }

            /// Compute \f$A(x)\f$.
            ::Spacy::Vector operator()( const ::Spacy::Vector& x ) const
            {
                assembleOperator( x );

                auto y = zero( range() );
                copy( *b_, y );
                return y;
            }

            /// Compute \f$A'(x)dx\f$.
            ::Spacy::Vector d1( const ::Spacy::Vector& x, const ::Spacy::Vector& dx ) const
            {
                assembleGradient( x );

                auto dx_ = dolfin::Function( cast_ref< VectorCreator >( domain().creator() ).get() );
                copy( dx, dx_ );
                auto y_ = dx_.vector()->copy();
                A_->mult( *dx_.vector(), *y_ );

                auto y = zero( range() );
                copy( *y_, y );

                return y;
            }

            /**
             * @brief Access \f$A'(x)\f$ as linear operator \f$X\rightarrow Y\f$
             * @see LinearOperator, ::Spacy::LinearOperator
             */
            auto linearization( const ::Spacy::Vector& x ) const
            {
                assembleGradient( x );
                assert( A_ != nullptr );
                assert( operatorSpace_ != nullptr );

                return LinearOperator{ A_->copy(), *operatorSpace_, J_.function_space( 0 ) };
            }

        private:
            /// Assemble discrete representation of \f$A(x)\f$.
            void assembleOperator( const ::Spacy::Vector& x ) const
            {
                if ( oldX_F && ( oldX_F == x ) )
                    return;

                auto x_ = std::make_shared< dolfin::Function >( J_.function_space( 0 ) );
                copy( x, *x_ );
                assign_x_if_present( F_, x_ );
                b_ = x_->vector()->factory().create_vector( x_->vector()->mpi_comm() );

                dolfin::Assembler{}.assemble( *b_, F_ );

                for ( const auto& bc : bcs_ )
                    bc->apply( *b_, *x_->vector() );

                oldX_F = x;
            }

            /// Assemble discrete representation of \f$A'(x)\f$.
            void assembleGradient( const ::Spacy::Vector& x ) const
            {
                if ( oldX_J && ( oldX_J == x ) )
                    return;

                auto x_ = std::make_shared< dolfin::Function >( J_.function_space( 0 ) );
                copy( x, *x_ );
                assign_x_if_present( J_, x_ );
                auto A = x_->vector()->factory().create_matrix( x_->vector()->mpi_comm() );

                dolfin::Assembler{}.assemble( *A, J_ );

                for ( const auto& bc : bcs_ )
                    bc->apply( *A, *x_->vector(), *x_->vector() );
                A_ = std::make_shared< dolfin::Matrix >( std::move( *A ) );

                oldX_J = x;
            }

            mutable ResidualForm F_;
            mutable JacobianForm J_;
            std::vector< const dolfin::DirichletBC* > bcs_ = {};
            mutable std::shared_ptr< dolfin::Matrix > A_ = nullptr;
            mutable std::shared_ptr< dolfin::GenericVector > b_ = nullptr;
            mutable ::Spacy::Vector oldX_F{}, oldX_J{};
            std::shared_ptr< VectorSpace > operatorSpace_ = nullptr;
        };

        /**
         * @brief Convenient generation of a differentiable operator \f$A: X\rightarrow Y\f$ as used
         * in %FEniCS.
         * @return @ref C1Operator "::Spacy::Fenics::C1Operator<ResidualForm,JacobianForm>( F , J ,
         * bcs , domain , range )"
         */
        template < class ResidualForm, class JacobianForm >
        auto makeC1Operator( ResidualForm& F, JacobianForm& J, const std::vector< const dolfin::DirichletBC* >& bcs,
                             const VectorSpace& domain, const VectorSpace& range )
        {
            return C1Operator< ResidualForm, JacobianForm >{ F, J, bcs, domain, range };
        }

        /**
         * @brief Convenient generation of a differentiable operator \f$A: X\rightarrow Y\f$ from
         * %FEniCS(dolfin).
         * @return @ref C1Operator "::Spacy::Fenics::C1Operator<ResidualForm,JacobianForm>( F , J ,
         * domain , range )"
         */
        template < class ResidualForm, class JacobianForm >
        auto makeC1Operator( ResidualForm& F, JacobianForm& J, const VectorSpace& domain, const VectorSpace& range )
        {
            return C1Operator< ResidualForm, JacobianForm >{ F, J, domain, range };
        }

        /**
         * @brief Convenient generation of a differentiable operator \f$A: X\rightarrow X\f$ from
         * %FEniCS(dolfin).
         * @return @ref C1Operator "::Spacy::Fenics::C1Operator<ResidualForm,JacobianForm>( F , J ,
         * bcs , space , space )"
         */
        template < class ResidualForm, class JacobianForm >
        auto makeC1Operator( ResidualForm& F, JacobianForm& J, const std::vector< const dolfin::DirichletBC* >& bcs,
                             const VectorSpace& space )
        {
            return C1Operator< ResidualForm, JacobianForm >{ F, J, bcs, space, space };
        }

        /**
         * @brief Convenient generation of a differentiable operator \f$A: X\rightarrow X\f$ from
         * %FEniCS(dolfin).
         * @return @ref C1Operator "Fenics::C1Operator<ResidualForm,JacobianForm>( F , J , space ,
         * space )"
         */
        template < class ResidualForm, class JacobianForm >
        auto makeC1Operator( ResidualForm& F, JacobianForm& J, const VectorSpace& space )
        {
            return C1Operator< ResidualForm, JacobianForm >{ F, J, space, space };
        }
    } // namespace FEniCS
    /** @} */
} // namespace Spacy
