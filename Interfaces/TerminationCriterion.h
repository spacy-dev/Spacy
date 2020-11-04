#pragma once

#include <Spacy/Vector.h>

namespace Spacy
{
    namespace CG
    {
        /// Type-erased termination criterion for conjugate gradient methods.
        class TerminationCriterion
        {
        public:
            bool operator()() const;

            void clear();

            void update( double alpha, double qAq, double qPq, double rPINVr, const Vector& x );

            bool vanishingStep() const;

            bool minimalDecreaseAchieved() const;

            void setEps( double eps );

            void setAbsoluteAccuracy( double accuracy );

            void setMinimalAccuracy( double accuracy );

            void setRelativeAccuracy( double accuracy );
        };
    } // namespace CG
} // namespace Spacy
