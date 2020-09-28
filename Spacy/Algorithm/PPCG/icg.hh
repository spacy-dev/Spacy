#pragma once

#include <Spacy/Operator.h>
#include <Spacy/Util/Mixins/Index.h>
#include <Spacy/Util/Mixins/Eps.h>
#include <Spacy/Util/Mixins/MaxSteps.h>
#include <Spacy/Util/Mixins.h>
#include <Spacy/Util/Base/OperatorBase.h>
#include <Spacy/Algorithm/CG/TerminationCriterion.h>
#include <Spacy/Spaces/ProductSpace/Operator_V2.h>

namespace Spacy
{

    class VectorSpace;
    class Vector;

        class InexactSolver
                : public OperatorBase,
                public Mixin::AdjointIndex,
                public Mixin::ControlIndex,
                public Mixin::StateIndex,
                public Mixin::Eps,
                public Mixin::MaxSteps,
                public Mixin::Verbosity,
                public Mixin::AbsoluteAccuracy,
                public Mixin::RelativeAccuracy

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
                InexactSolver(CallableOperator P, const ::Spacy::ProductSpace::Operator_V2& H,
                              CG::TerminationCriterion termination, const VectorSpace& domain, const VectorSpace& range, CallableOperator diagMuu, Real theta_ = 0.0);

                //InexactSolver(CallableOperator P,
                //              CallableOperator BT, CallableOperator Hyy, CallableOperator Hyu, CallableOperator Huu, CallableOperator Huy,
                //              CG::TerminationCriterion termination, const VectorSpace& domain, const VectorSpace& range, CallableOperator Muu,CallableOperator diagMuu, Real theta_ = 0.0);

                bool isPositiveDefinite() const noexcept;



                Vector solve( const Vector& b) const;

                Vector operator()( const Vector& b) const;

                void setRegularizationOperators(CallableOperator HyyTheta, CallableOperator HuuTheta);

                bool vanishingStep( unsigned step, const std::vector<Real> & alpha_vec, const std::vector<Real> & beta_vec ) const;

                unsigned int getIterations() const;
                unsigned int getchebIterations() const;
                unsigned int getchebTIterations() const;

                bool isSingular() const;

                std::function<void(const ::Spacy::Vector& x)> output;

                double getConditionNumber();
            private:
                CallableOperator P_;

                ::Spacy::ProductSpace::Operator_V2 H_;

                //CallableOperator BT_;
                //CallableOperator Hyy_;
                //CallableOperator Hyu_; //Für TangentialSolver (bei NormalSolver = 0)
                //CallableOperator Huu_;
                //CallableOperator Huy_; //Für TangentialSolver (bei NormalSolver = 0)


                CallableOperator HyyTheta_ = [](const Vector & x) { return 0.0*x;};
                CallableOperator HuuTheta_ = [](const Vector & x) { return 0.0*x;};


                CallableOperator diagMuu_ = [](const Vector & x) { return 0.0*x;};


                mutable CG::TerminationCriterion terminate_;

                //CallableOperator Muu_;

                mutable Result result_ = Result::Failed;
                mutable DefiniteNess definiteness_ = DefiniteNess::PositiveDefinite;

                Real theta_ = 0.0;

                mutable unsigned int iterations_ = 1;
                mutable unsigned chebIterations_;
                mutable unsigned chebTIterations_;

                mutable double condPrecondMatrix =1;



        };

}
