#pragma once

#include "PPCG_Test.h"

#include <Spacy/Operator.h>
#include <Spacy/Vector.h>
#include <Spacy/VectorSpace.h>

#include <Spacy/Util/Mixins/Index.h>
#include <Spacy/Util/Mixins/Eps.h>
#include <Spacy/Util/Mixins/MaxSteps.h>
#include <Spacy/Util/Mixins.h>
#include <Spacy/Util/Base/OperatorBase.h>
#include <Spacy/Algorithm/CG/TerminationCriterion.h>
#include <Spacy/Algorithm/CG/TerminationCriteria.h>

namespace Spacy
{
    
    namespace ModifiedPPCG
    {
        
        class TriangularConstraintPreconditioner;
        using namespace Spacy;
        
        /// Implements projected preconditioned cg for the case of optimal control systems
        /// Requires: BlockTriangularConstraint preconditioner  P
        ///           Upper left 2x2-Block M
        ///           Control-Block minusB
        /// Does not require differential operator block A
        class Solver
        : public OperatorBase, public Mixin::Eps, public Mixin::MaxSteps, public Mixin::Verbosity, public Mixin::AbsoluteAccuracy, public Mixin::RelativeAccuracy
        
        {
            enum class Result
            {
                NotConverged,
                Converged,
                Failed,
                EncounteredNonConvexity,
                TruncatedAtNonConvexity,
                PreconditionerNonConvexity,
                TooSingular
            };
            enum class DefiniteNess
            {
                PositiveDefinite,
                Indefinite
            };
            
        public:
            Solver(const TriangularConstraintPreconditioner& P, const Operator& My, const Operator& Mu, const OperatorWithTranspose& minusB, const VectorSpace& domainIn, const PPCG_Test& test, Real theta = 0.0 );
            
            bool isPositiveDefinite() const noexcept;
            
            template < class Criterion >
            void setTerminationCriterion( Criterion newTerminate ) { terminate_ = std::move( newTerminate ); }
            
            CG::TerminationCriterion& terminationCriterion() noexcept { return terminate_; }
            
            /// Solve problem
            Vector solve( const Vector& b) const;
            
            /// Solve problem
            virtual Vector operator()( const Vector& b) const;
            
            void setRegularizationOperators(const CallableOperator& R);
            
            bool vanishingStep( unsigned step, const std::vector<Real> & alpha_vec, const std::vector<Real> & beta_vec ) const;
            
            /// Required iterations in last run
            unsigned int getIterations() const;
            
            /// Accumulates required iterations over several runs
            unsigned int getIterationsInLifeTime() const;

            bool isSingular() const;
            
            std::function<void(const ::Spacy::Vector& x)> output;
            
            double getConditionNumber();
        private:
            /// Preconditioner
            const TriangularConstraintPreconditioner& P_;

            /// Hessian of full quadratic functional            
            const Operator& My_;
            const Operator& Mu_;

            /// Negative of control operator
            const OperatorWithTranspose& minusB_;
            
            /// Regularization part
            CallableOperator R_ = [](const Vector & x) { return 0.0*x;};

             /// Required vector spaces
            const VectorSpace &stateSpace_, &controlSpace_, &adjointSpace_;
            
            const PPCG_Test& test_;
            
            
            mutable CG::TerminationCriterion terminate_={ CG::Termination::StrakosTichyEnergyError{} };
            mutable CG::Termination::StrakosTichyEnergyError stee_;
            
            
            mutable Result result_ = Result::Failed;
            mutable DefiniteNess definiteness_ = DefiniteNess::PositiveDefinite;
            
            Real theta_ = 0.0;
            
            mutable unsigned int iterations_ = 0;
            mutable unsigned int iterationsInLifeTime_ = 0;
            
            mutable double condPrecondMatrix =1;
            
            
            
        };
        
    }
}
