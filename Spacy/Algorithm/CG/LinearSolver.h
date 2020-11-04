#pragma once

#include "CG.h"
#include "Regularization.h"

#include <Spacy/Util/Base/OperatorBase.h>
#include <Spacy/Util/Mixins.h>

namespace Spacy
{
    /** @addtogroup CGGroup @{ */
    namespace CG
    {
        /// Conjugate gradient solver.
        class LinearSolver : public OperatorBase,
                             public Mixin::AbsoluteAccuracy,
                             public Mixin::RelativeAccuracy,
                             public Mixin::Eps,
                             public Mixin::IterativeRefinements,
                             public Mixin::MaxSteps,
                             public Mixin::Verbosity
        {
        public:
            /**
             * @brief Set up conjugate gradient solver
             * @param A_ operator
             * @param P_ preconditioner
             * @param regularization regularization for the treatment of nonconvexities
             * @param truncated set to true if the algorithm should terminate on nonconvexities
             */
            LinearSolver( Operator A_, CallableOperator P_, bool truncated, Regularization regularization );

            LinearSolver( const LinearSolver& other );

            /// Apply conjugate gradient solver.
            Vector operator()( const Vector& y ) const;

            /// Access conjugate gradient implementation.
            Solver& impl();

            /// Checks positive definiteness of \f$A\f$.
            bool isPositiveDefinite() const noexcept;

            /// Access preconditioner \f$P\f$.
            const CallableOperator& P() const noexcept;

            /// Access operator \f$A\f$.
            const CallableOperator& A() const noexcept;

            /// Access regularization \f$R\f$.
            const Regularization& R() const noexcept;

            /// Set Termination Criterion
            template < class TermCrit >
            void setTerminationCriterion( TermCrit term_ ) const
            {
                return cg.setTerminationCriterion( term_ );
            }

            /// getter for Termination Criterion
            CG::TerminationCriterion& terminationCriterion() noexcept
            {
                return cg.terminationCriterion();
            }

            /// getter for number of iteration
            unsigned getIterations() const
            {
                return cg.getIterations();
            }

        private:
            mutable Solver cg;
        };
    } // namespace CG

    /**
     * @brief Construct conjugate gradient method.
     * @param A operator
     * @param P preconditioner
     * @param relativeAccuracy relative accuracy
     * @param eps maximal attainable accuracy
     * @param verbose verbosity
     * @return CGSolver(A,P)
     */
    CG::LinearSolver makeCGSolver( Operator A, CallableOperator P, Real relativeAccuracy = 1e-15, Real eps = 1e-15, bool verbose = false );

    /**
     * @brief Construct regularized conjugate gradient method.
     * @param A operator
     * @param P preconditioner
     * @param relativeAccuracy relative accuracy
     * @param eps maximal attainable accuracy
     * @param verbose verbosity
     * @return CGSolver(A,P,false,RegularizeViaPreconditioner())
     */
    CG::LinearSolver makeRCGSolver( Operator A, CallableOperator P, Real relativeAccuracy = 1e-15, Real eps = 1e-15, bool verbose = false );

    /**
     * @brief Construct truncated conjugate gradient method.
     * @param A operator
     * @param P preconditioner
     * @param relativeAccuracy relative accuracy
     * @param eps maximal attainable accuracy
     * @param verbose verbosity
     * @return CGSolver(A,P,true)
     */
    CG::LinearSolver makeTCGSolver( Operator A, CallableOperator P, Real relativeAccuracy = 1e-15, Real eps = 1e-15, bool verbose = false );

    /**
     * @brief Construct truncated regularized conjugate gradient method (regularized via
     * preconditioner).
     * @param A operator
     * @param P preconditioner
     * @param relativeAccuracy relative accuracy
     * @param eps maximal attainable accuracy
     * @param verbose verbosity
     * @return CGSolver(A,P,true,RegularizeViaPreconditioner())
     */
    CG::LinearSolver makeTRCGSolver( Operator A, CallableOperator P, Real relativeAccuracy = 1e-15, Real eps = 1e-15,
                                     bool verbose = false );

    /**
     * @brief Construct truncated regularized conjugate gradient method (regularized via arbitrary
     * operator).
     * @param A operator
     * @param P preconditioner
     * @param R regularization
     * @param theta_sugg suggestion for starting value of regularization paramter
     * @param relativeAccuracy relative accuracy
     * @param eps maximal attainable accuracy
     * @param verbose verbosity
     * @return CGSolver(A,P,true,RegularizeViaPreconditioner())
     */
    CG::LinearSolver makeTRCGSolver( Operator A, CallableOperator P, CallableOperator R, Real theta_sugg, Real relativeAccuracy = 1e-15,
                                     Real eps = 1e-15, bool verbose = false );
    /** @} */
} // namespace Spacy
