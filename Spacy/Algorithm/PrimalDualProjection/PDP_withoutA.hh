#pragma once

#include "ModifiedPPCG.h"
#include "Spacy/Algorithm/Preconditioner/Chebyshev.h"
#include "TriangConstraintPrec.h"

#include <Spacy/Spaces/ProductSpace/Operator_V2.h>
#include <Spacy/Spaces/ScalarSpace/Real.h>

#include "../CG/CG.h"
#include "../CG/RegularizeViaPreconditioner.h"
#include "../CG/TerminationCriterion.h"
#include "../CG/TerminationCriteria.h"


namespace Spacy
{
    class Vector;

    namespace PDP_withoutA
    {
        
/// Implements PDP method for optimal control problems
        class Solver :
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
               Solver(Operator M, Operator stateSolver, Operator adjointSolver, OperatorWithTranspose minusB, CallableOperator surrogateStateSolver, CallableOperator surrogateAdjointSolver, CallableOperator controlSolver, const VectorSpace& totalSpace);

                /**
                 * @brief Solving the Problem
                 * @param b right hand side
                 */
                Vector operator()(const Vector& b) const;


                unsigned getIterations() const
                {
                    return iterations_;
                }

                bool isPositiveDefinite() const;
                bool isAOperatorPositiveDefinite() const noexcept;
                std::function<bool(bool setSignal, bool isConvex)> signalConvex_ = [](bool setSignal, bool isConvex){ return true; };

                void setInternalAccuracies(double stateAcc, double adjointAcc, double modifiedCGAcc, double chebyshevAcc);

            private:
                mutable bool convex_flag_ = true;
                mutable bool energyIsConvex_ = true;

                Operator M_;
                Operator stateSolver_;
                Operator adjointSolver_;
                OperatorWithTranspose minusB_;
                CallableOperator surrogateStateSolver_;
                CallableOperator surrogateAdjointSolver_;
                CallableOperator controlSolver_;

                mutable ::Spacy::Real normDxOld_ = 0; 
                mutable ::Spacy::Real normDx_ = 0; 
                Vector PDPLoop(Vector x, const Vector& b) const;
                bool convergenceTest(bool isConvex, bool energyIsConvex, CallableOperator H , const Vector &x, const Vector &dx) const;
                
                mutable Result result_ = Result::Failed;
                mutable DefiniteNess definiteness_ = DefiniteNess::PositiveDefinite;
                mutable Real theta_ = 0.0;
                mutable bool projected_ = false;

                const VectorSpace &totalSpace_, &primalSpace_, &adjointSpace_, &stateSpace_, &controlSpace_;
                
                mutable unsigned iterations_;

                mutable unsigned dual_iterations_ = 0;
                mutable unsigned proj_iterations_ = 0;

                double stateAcc_ = 1e-2;
                double adjointAcc_ = 1e-2;
                double modifiedCGAcc_ = 1e-2;
                double chebyshevAcc_ = 1e-2;


                mutable double sumNormSquaredUpdate_ = 0.0;
        };
    }
}
