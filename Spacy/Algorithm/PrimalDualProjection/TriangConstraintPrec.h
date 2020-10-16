// Copyright (C) 2015 by Lars Lubkoll. All rights reserved.
// Released under the terms of the GNU General Public License version 3 or later.

#ifndef SPACY_CONJUGATE_GRADIENTS_TRIANGULARCONSTRAINTPRECONDITIONER_HH
#define SPACY_CONJUGATE_GRADIENTS_TRIANGULARCONSTRAINTPRECONDITIONER_HH


///---------- Originalversion in Algorithm/CG, hier nur kleinere Anpassungen für PPCG (output,...) -------------


#include "Spacy/LinearSolver.h"
#include "Spacy/Operator.h"
#include "Spacy/Vector.h"
#include "Spacy/Util/Mixins/Index.h"
#include "Spacy/Util/Base/OperatorBase.h"
#include <Spacy/Algorithm/CG/LinearSolver.h>

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
            TriangularConstraintPreconditioner( ::Spacy::LinearSolver stateSolver,
                                                ::Spacy::LinearSolver controlSolver,
                                                ::Spacy::LinearSolver adjointSolver,
                                                OperatorWithTranspose minusB, 
                                                const VectorSpace& domain,
                                                const VectorSpace& stateSpace
            );
            
            /**
             * @brief Apply preconditioner \f$P\f$.
             * @param x argument
             * @return \f$P(x)\f$
             */
            virtual Vector operator()(Vector x) const;
            
            
            /**
             * @brief Apply preconditioner \f$P\f$ 
             * @param x argument, q auxilliary vector, needed for inexact differential block
             * @return first argument: \f$P(x)\f$
             *
             */
            Vector operator()(const Vector& x, Vector& q) const;
            
            
            std::function<void(const ::Spacy::Vector& x)> output;
            
            unsigned getBPXCallsInLifeTime() const
            {
                return iterationsInLifeTime_;
            }
            
            Real getSigma() const 
            { return sigma_; }
            
            void setSpectralBounds(double eigMin, double eigMax);
            
            std::function<bool(bool setSignal, bool isConvex)> signalConvex_ = [](bool setSignal, bool isConvex){ return true; };
            
        private:
            
            void apply(const Vector& x, Vector& y, Vector& q) const;
            
            ::Spacy::LinearSolver stateSolver_, controlSolver_, adjointSolver_;
            const ::Spacy::VectorSpace& stateSpace_, &controlSpace_, &adjointSpace_;
            OperatorWithTranspose minusB_;
            
            mutable Real sigma_ = 0.0;
            
            mutable int iterationsInLifeTime_ = 0;
        };
    }
}

#endif // SPACY_CONJUGATE_GRADIENTS_TRIANGULARCONSTRAINTPRECONDITIONER_HH
