// Copyright (C) 2015 by Lars Lubkoll. All rights reserved.
// Released under the terms of the GNU General Public License version 3 or later.

#ifndef SPACY_CONJUGATE_GRADIENTS_TRIANGULARCONSTRAINTPRECONDITIONER_HH
#define SPACY_CONJUGATE_GRADIENTS_TRIANGULARCONSTRAINTPRECONDITIONER_HH


///---------- Originalversion in Algorithm/CG, hier nur kleinere Anpassungen für PPCG (output,...) -------------


#include "Spacy/Operator.h"
#include "Spacy/Vector.h"
#include "Spacy/Util/Mixins/Index.h"
#include "Spacy/Util/Base/OperatorBase.h"
#include <Spacy/Adapter/KaskadeParabolic/PDESolverBase.h>

namespace Spacy
{
    /// \cond
    class VectorSpace;
    /// \endcond
    
    namespace ModifiedPPCG
    {
        /**
         * @brief A triangular state constraint preconditioner used in @cite Lubkoll2015a.
         *
         * A preconditioner for optimal control problems of the form
         * \f[ P=\left( \begin{array}{ccc}  &  & A^{-T} \\  & M^{-1} & -B^T \\ A^{-1} & -B &  \end{array} \right)\f],
         * where \f$A\f$ is the state operator, \f$B\f$ the control operator and \f$M\f$ induces a scalar product on the control space.
         */
        class TriangularConstraintPreconditioner
        : public OperatorBase
        {
        public:
            /**
             * @brief Constructor.
             * @param stateSolver solver for the evaluation of \f$A^{-1}r\f$
             * @param controlSolver solver for the evaluation of \f$M^{-1}r\f$
             * @param adjointSolver solver for the evaluation of \f$A^{-T}r\f$
             * @param B control operator
             * @param BT adjoint control operator
             * @param domain domain space
             * @param range range space
             */
            TriangularConstraintPreconditioner(  const ::Spacy::SolverBase& stateSolver,
                                                 const Operator& controlSolver,
                                                 const ::Spacy::SolverBase& adjointSolver,
                                                 const OperatorWithTranspose& minusB, 
                                                 const VectorSpace& domain
            );
            
            /**
             * @brief Apply preconditioner \f$P\f$.
             * @param x argument
             * @return \f$P(x)\f$
             */
            virtual void operator()(Vector y, Vector u , Vector& outy, Vector& outu, Vector& outp) const;
                        
            unsigned getBPXCallsInLifeTime() const
            {
                return iterations_;
            }
            
            Real getSigma() const 
            { return sigma_; }
            
            void setSpectralBounds(double eigMin, double eigMax);
            
            std::function<bool(bool setSignal, bool isConvex)> signalConvex_ = [](bool setSignal, bool isConvex){ return true; };
      
            void apply( Vector& outy, Vector& outu, Vector& outp, Vector&& iny, Vector&& inu ) const;  // in will be overwritten
            
        private:
            
            
            const ::Spacy::SolverBase &stateSolver_, &adjointSolver_;
            const Operator& controlSolver_;
            const ::Spacy::VectorSpace& stateSpace_, &controlSpace_, &adjointSpace_;
            const OperatorWithTranspose& minusB_;
            
            mutable int iterations_ = 0;
            
            mutable Real sigma_ = 0.0;
        };
    }
}

#endif // SPACY_CONJUGATE_GRADIENTS_TRIANGULARCONSTRAINTPRECONDITIONER_HH
