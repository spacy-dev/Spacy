/**
    ToDO: Zum Verwenden von ppcg in Spacy nötig:
            - Operator_V2 in Spacy/Spaces/ProductSpace/
            - Anpassung im cg -> callback (Für EW-Berechnung)
          Zum Testen mit einem Example:
            - C2BlockFunctional in Adapter/Kaskade/ und in Spacy/Adapter/kaskade.hh entsprechende Zeile einfügen
            - Im AffineCovariantSolver.h/.cpp die entsprechenden Anpassungen für PPCG tätigen
*/

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

    namespace PDP
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
               Solver(Operator M, OperatorWithTranspose A, OperatorWithTranspose AT, OperatorWithTranspose minusB, CallableOperator BPX, CallableOperator BPXT,
                      CallableOperator controlSolver, ::Spacy::CG::Regularization R,const VectorSpace& totalSpace);

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
                
                Vector C_(const Vector& ) const;
                Vector C_transpose(const Vector& ) const;
                
                mutable bool cheb_created_ = false;
                mutable ::Spacy::Preconditioner::Chebyshev ChebPrec_;

                mutable bool convex_flag_ = true;
                mutable bool energyIsConvex_ = true;

                Operator M_;
                OperatorWithTranspose minusB_;
                OperatorWithTranspose A_;
                OperatorWithTranspose AT_;
                CallableOperator BPX_;
                CallableOperator BPXT_;

                ::Spacy::CG::Regularization R_;
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
