#ifndef SPACY_MODIFIED_PPCG_TEST_HH
#define SPACY_MODIFIED_PPCG_TEST_HH

#include "Spacy/Operator.h"
#include "Spacy/Vector.h"
#include "Spacy/Util/Mixins/Index.h"
#include "Spacy/Util/Base/OperatorBase.h"
#include <Spacy/Adapter/KaskadeParabolic/PDESolverBase.h>

namespace Spacy
{
    namespace ModifiedPPCG
    {
        class PPCG_Test
        {
        public:
            PPCG_Test(  const ::Spacy::SolverBase& stateSolver,
                        const Operator& controlSolver,
                        const ::Spacy::SolverBase& adjointSolver,
                        const Operator& My, const Operator& Mu,
                        const OperatorWithTranspose& minusB 
            );
            
            void operator()(const Vector& dy, const Vector& du, const Vector& dp, const Vector& ru) const;
            void diffVec(const Vector& dy, const Vector& du, const Vector& dp) const;
            
        private:
            const ::Spacy::SolverBase &stateSolver_, &adjointSolver_;
            const Operator& controlSolver_, &My_, &Mu_;
            const ::Spacy::VectorSpace& stateSpace_, &controlSpace_, &adjointSpace_;
            const OperatorWithTranspose& minusB_;
            mutable Spacy::Vector yl, ul, pl;
        };
    }
}

#endif // SPACY_MODIFIED_PPCG_TEST_HH
