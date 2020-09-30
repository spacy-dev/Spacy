#pragma once

#include "../PPCG/icg.hh"
#include "../PPCG/triangularStateConstraintPreconditioner.hh"

#include <Spacy/Spaces/ProductSpace/Operator_V2.h>
#include <Spacy/Spaces/ScalarSpace/Real.h>

#include "../CG/CG.h"
#include "../CG/LinearSolver.h"
#include "../CG/RegularizeViaPreconditioner.h"
#include "../CG/TerminationCriterion.h"
#include "../CG/TerminationCriteria.h"


namespace Spacy
{
    class Vector;

    namespace ICG
    {
        class Solver
                : public OperatorBase,
                public Mixin::AdjointIndex,
                public Mixin::ControlIndex,
                public Mixin::StateIndex,
                public Mixin::Eps,
                public Mixin::MaxSteps,
                public Mixin::RelativeAccuracy,
                public Mixin::Verbosity
        {
            enum class Result
            {
                NotConverged,
                Converged,
                Failed,
                EncounteredNonConvexity,
                TruncatedAtNonConvexity,
                PreconditionerNonConvexity
            };

            enum class DefiniteNess
            {
                PositiveDefinite,
                Indefinite
            };

            public:
              /**
                 * \brief Set up projected preconditioned conjugate gradient solver.
                 *
                 * \param H linear product space operator
                 * \param BPX
                 * \param BPXT
                 * \param controlSolver
                 * \param LAPACKE_dstevr
                 */
                Solver(const ::Spacy::ProductSpace::Operator_V2& H, const ::Spacy::CallableOperator& BPX, const ::Spacy::CallableOperator& BPXT,
                       const ::Spacy::CallableOperator& controlSolver, const ::Spacy::CallableOperator& diagMuu);
               /**
                * \brief Set up projected preconditioned conjugate gradient solver.
                *
                * \param H linear product space operator
                * \param BPX
                * \param BPXT
                * \param controlSolver
                * \param LAPACKE_dstevr
                * \param Ry,Ru Regularization operators
                */
               Solver(const ::Spacy::ProductSpace::Operator_V2& H, const ::Spacy::CallableOperator& BPX, const ::Spacy::CallableOperator& BPXT,
                      const ::Spacy::CallableOperator& controlSolver,
                      const ::Spacy::CallableOperator& Ry, const ::Spacy::CallableOperator& Ru, const ::Spacy::CallableOperator& diagMuu);

                /**
                 * @brief Solving the Problem
                 * @param b right hand side
                 */
                Vector operator()(const Vector& b) const;

                /**
                 * @brief Access Operator.
                 * @return operator \f$H\f$
                 */
                const ::Spacy::ProductSpace::Operator_V2& H() const;

                /**
                 * @brief Access preconditioner.
                 * @return preconditioner \f$P\f$
                 */
                const CallableOperator& BPX() const;

                /**
                 * @brief Access preconditioner.
                 * @return preconditioner \f$PT\f$
                 */
                const CallableOperator& BPXT() const;

                unsigned getIcgIterations() const
                {
                    return icg_iterations_;
                }

                void setNorm(CallableOperator Myy , CallableOperator Muu);

                const ::Spacy::CallableOperator& controlSolver() const;

                bool isPositiveDefinite() const;
                bool isAOperatorPositiveDefinite() const noexcept;
                std::function<bool(bool setSignal, bool isConvex)> signalConvex_ = [](bool setSignal, bool isConvex){ return true; };

                std::function<void(::Spacy::Vector & x_u, ::Spacy::Vector & x_p,::Spacy::Vector & b_u, ::Spacy::Vector Mu_x_u, ::Spacy::Vector BT_x_p )> output;
                std::function<void(const ::Spacy::Vector& x)> output2;
                std::function<::Spacy::Vector(::Spacy::Vector,::Spacy::Vector)> transfer = [](::Spacy::Vector src, ::Spacy::Vector res){ return src;};
            private:
                mutable bool convex_flag_ = true;
                mutable bool energyIsConvex_ = true;

                ::Spacy::ProductSpace::Operator_V2 H_;
                ::Spacy::CallableOperator BPX_;
                ::Spacy::CallableOperator BPXT_;

                ::Spacy::CallableOperator controlSolver_;

                ::Spacy::CallableOperator Myy_ = [](const Vector & x) { return 0.0*x;};
                ::Spacy::CallableOperator Muu_ = [](const Vector & x) { return 0.0*x;};

                ::Spacy::CallableOperator Ry_ = [](const Vector & x) { return 0.0*x;};
                ::Spacy::CallableOperator Ru_ = [](const Vector & x) { return 0.0*x;};

                ::Spacy::CallableOperator diagMuu_ = [](const Vector & x) { return 0.0*x;};

                Vector outer_pcgLoop(Vector x, const Vector& b) const;
                bool convergenceTest(bool isConvex, bool energyIsConvex, CallableOperator Hyy, CallableOperator Huu, const Vector & x_y, const Vector & x_u,  const Vector & dx_y, const Vector & dx_u ) const;

                mutable Result result_ = Result::Failed;
                mutable DefiniteNess definiteness_ = DefiniteNess::PositiveDefinite;
                mutable Real theta_ = 0.0;

                mutable unsigned icg_iterations_;

        };
    }
}
