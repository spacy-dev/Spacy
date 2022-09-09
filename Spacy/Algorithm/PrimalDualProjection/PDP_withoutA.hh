#pragma once

#include "ModifiedPPCG.h"
#include "Spacy/Algorithm/Preconditioner/Chebyshev.h"
#include "TriangConstraintPrec.h"
#include "PPCG_Test.h"

#include <Spacy/Spaces/ProductSpace/Operator_V2.h>
#include <Spacy/Spaces/ScalarSpace/Real.h>

#include <Spacy/Adapter/KaskadeParabolic/PDESolverBase.h>

#include <limits>

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
               Solver(Operator M, Operator Mys, Operator Mus, Spacy::SolverBase stateSolver, Spacy::SolverBase adjointSolver, OperatorWithTranspose minusB, OperatorWithTranspose minusBs, Spacy::SolverBase surrogateStateSolver, Spacy::SolverBase surrogateAdjointSolver, Operator controlSolver, Operator controlSolverS, const VectorSpace& totalSpace, const VectorSpace& totalSpaceS);

                /**
                 * @brief Solving the Problem
                 * @param b right hand side
                 */
                Vector operator()(const Vector& b) const;


                unsigned getIterations() const
                {
                    return iterations_;
                }
                
                unsigned getPPCGIterations() const
                {
                    return PPCGiterations_;
                }
                
                unsigned getMGIterations() const
                {
                    return MGiterations_;
                }
                
                bool isPositiveDefinite() const;
                bool isAOperatorPositiveDefinite() const noexcept;
                std::function<bool(bool setSignal, bool isConvex)> signalConvex_ = [](bool setSignal, bool isConvex){ return true; };

                void setInternalAccuracies(double stateAcc, double adjointAcc, double modifiedCGAcc, double chebyshevAcc);

                mutable int MillisecondsPPCG = 0;
                
            private:
                mutable bool convex_flag_ = true;
                mutable bool energyIsConvex_ = true;

                Operator M_;
                Operator Mu_s;
                Operator My_s;
                Spacy::SolverBase stateSolver_;
                Spacy::SolverBase adjointSolver_;
                OperatorWithTranspose minusB_;
                OperatorWithTranspose minusB_s;
                Spacy::SolverBase surrogateStateSolver_;
                Spacy::SolverBase surrogateAdjointSolver_;
                Operator controlSolver_;
                Operator controlSolver_s;

                mutable ::Spacy::Real normDxOld_ = 0; 
                mutable ::Spacy::Real normDx_ = 0; 
                mutable ::Spacy::Real Theta = 0; 
                mutable ::Spacy::Real ThetaMin = 1; 
                mutable ::Spacy::Real ThetaMax = 0; 
                Vector PDPLoop(Vector x, const Vector& b) const;
                bool convergenceTest(bool isConvex, bool energyIsConvex) const;
                
                mutable Result result_ = Result::Failed;
                mutable DefiniteNess definiteness_ = DefiniteNess::PositiveDefinite;
                mutable Real theta_ = 0.0;
                mutable bool projected_ = false;

                const VectorSpace &totalSpace_, &adjointSpace_, &stateSpace_, &controlSpace_, &totalSpace_s;
                
                mutable unsigned iterations_, PPCGiterations_, MGiterations_;

                mutable unsigned dual_iterations_ = 0;
                mutable unsigned proj_iterations_ = 0;

                double stateAcc_ = 1e-2;
                double adjointAcc_ = 1e-2;
                double modifiedCGAcc_ = 1e-2;
                double chebyshevAcc_ = 1e-2;

        public:
                mutable double sumNormSquaredUpdate_ = 0.0;
        };
    }
}
