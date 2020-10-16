// Copyright (C) 2015 by Lars Lubkoll. All rights reserved.
// Released under the terms of the GNU General Public License version 3 or later.

#pragma once


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

  namespace MINRES
  {
    /**
     * @brief A triangular state constraint preconditioner used in @cite Lubkoll2015a.
     *
     * A preconditioner for optimal control problems of the form
     * \f[ P=\left( \begin{array}{ccc}  &  & A^{-T} \\  & M^{-1} & -B^T \\ A^{-1} & -B &  \end{array} \right)\f],
     * where \f$A\f$ is the state operator, \f$B\f$ the control operator and \f$M\f$ induces a scalar product on the control space.
     */
    class BlockDiagPreconditioner 
    {
    public:
      /**
       * @brief Constructor.
       * @param stateSolver solver for the evaluation of \f$A^{-1}r\f$
       * @param controlSolver solver for the evaluation of \f$M^{-1}r\f$
       * @param domain domain space
       * @param range range space
       */
      BlockDiagPreconditioner(const ::Spacy::CallableOperator& stateAdjSolver, 
                              const ::Spacy::CallableOperator& controlSolver, 
                              const ::Spacy::VectorSpace& stateSpace, 
                              const ::Spacy::VectorSpace& controlSpace, 
                              const ::Spacy::VectorSpace& adjointSpace);
      /**
       * @brief Apply preconditioner \f$P\f$.
       * @param x argument
       * @return \f$P(x)\f$
       */
      Vector operator()(const Vector& x) const;

      
      /**
       * @brief Apply preconditioner \f$P\f$ 
       * @param x argument, q auxilliary vector, needed for inexact differential block
       * @return first argument: \f$P(x)\f$
       *
       */
      Vector operator()(const Vector& x, Vector& q) const;


    private:
    
      void apply(const Vector& x, Vector& y, Vector& q) const;
    
      ::Spacy::CallableOperator stateAdjSolver_, controlSolver_;
      
      const ::Spacy::VectorSpace& stateSpace_;
      const ::Spacy::VectorSpace& controlSpace_;
      const ::Spacy::VectorSpace&  adjointSpace_;


    };
    
    class BlockDiagSchurPreconditioner 
    {
    public:
      /**
       * @brief Constructor.
       * @param stateSolver solver for the evaluation of \f$A^{-1}r\f$
       * @param controlSolver solver for the evaluation of \f$M^{-1}r\f$
       * @param domain domain space
       * @param range range space
       */
      BlockDiagSchurPreconditioner(const ::Spacy::CallableOperator& stateAdjSolver, 
                              const ::Spacy::CallableOperator& stateL2Solver, 
                              const ::Spacy::CallableOperator& controlSolver, 
                              const ::Spacy::CallableOperator& stateOperator, 
                              const ::Spacy::VectorSpace& stateSpace, 
                              const ::Spacy::VectorSpace& controlSpace, 
                              const ::Spacy::VectorSpace& adjointSpace);
      /**
       * @brief Apply preconditioner \f$P\f$.
       * @param x argument
       * @return \f$P(x)\f$
       */
      Vector operator()(const Vector& x) const;

      
      /**
       * @brief Apply preconditioner \f$P\f$ 
       * @param x argument, q auxilliary vector, needed for inexact differential block
       * @return first argument: \f$P(x)\f$
       *
       */
      Vector operator()(const Vector& x, Vector& q) const;


    private:
    
      void apply(const Vector& x, Vector& y, Vector& q) const;
    
      ::Spacy::CallableOperator stateAdjSolver_, controlSolver_, stateL2Solver_;
      ::Spacy::CallableOperator stateOperator_;
      
      const ::Spacy::VectorSpace& stateSpace_;
      const ::Spacy::VectorSpace& controlSpace_;
      const ::Spacy::VectorSpace&  adjointSpace_;


    };

  }
}


