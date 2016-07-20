// Copyright (C) 2015 by Lars Lubkoll. All rights reserved.
// Released under the terms of the GNU General Public License version 3 or later.

#ifndef SPACY_ADAPTER_KASKADE_LINEAR_OPERATOR_HH
#define SPACY_ADAPTER_KASKADE_LINEAR_OPERATOR_HH

#include <utility>
#include <iostream>

#include "Spacy/linearSolver.hh"
#include "Spacy/vector.hh"
#include "Spacy/Util/Base/operatorBase.hh"
#include "Spacy/Util/Base/vectorBase.hh"

#include "operatorSpace.hh"
#include "directSolver.hh"
#include "util.hh" // copy

namespace Spacy
{
  /** @addtogroup KaskadeGroup
   * @{
   */

  namespace Kaskade
  {
    /**
     * @brief Linear operator interface for operators in %Kaskade 7.
     * @tparam OperatorImpl %Kaskade 7 operator, i.e. %Kaskade::AssembledGalerkinOperator or %Kaskade::MatrixRepresentedOperator.
     * @tparam AnsatzVariableSetDescription %Kaskade::VariableSetDescription for ansatz variables
     * @tparam TestVariableSetDescription %Kaskade::VariableSetDescription for test variables
     * @see ::Spacy::LinearOperator
     */
    template <class AnsatzVariableSetDescription, class TestVariableSetDescription>
    class LinearOperator :
        public OperatorBase ,
        public VectorBase ,
        public AddArithmeticOperators< LinearOperator<AnsatzVariableSetDescription,TestVariableSetDescription> >
    {
      using Spaces = typename AnsatzVariableSetDescription::Spaces;
      using Variables = typename AnsatzVariableSetDescription::Variables;
      using Domain = typename AnsatzVariableSetDescription::template CoefficientVectorRepresentation<>::type;
      using Range = typename TestVariableSetDescription::template CoefficientVectorRepresentation<>::type;
      using Matrix = ::Kaskade::MatrixAsTriplet<double>;
      using OperatorImpl = ::Kaskade::MatrixRepresentedOperator<Matrix,Domain,Range>;
      using OperatorCreator = LinearOperatorCreator<AnsatzVariableSetDescription,TestVariableSetDescription>;
      

    public:
      /**
       * @brief Construct linear operator for %Kaskade 7.
       * @param A operator from %Kaskade 7
       * @param domain domain space
       * @param range range space
       */
      LinearOperator(std::shared_ptr<OperatorImpl> A, const VectorSpace& space)
        : OperatorBase(cast_ref<OperatorCreator>(space.creator()).domain(),
                       cast_ref<OperatorCreator>(space.creator()).range()),
          VectorBase(space) ,
          A_(A),
          spaces_( extractSpaces<AnsatzVariableSetDescription>(cast_ref<OperatorCreator>(space.creator()).domain()) ),
          x_(cast_ref<OperatorCreator>(space.creator()).domain().zeroVector())
      {}

      /**
       * @brief Construct linear operator for %Kaskade 7.
       * @param A operator from %Kaskade 7
       * @param domain domain space
       * @param range range space
       * @param solverCreator function/functor implementing the creation of a linear solver
       */
      LinearOperator(std::shared_ptr<OperatorImpl> A, 
                     const VectorSpace& space, 
                     std::function<LinearSolver(const LinearOperator&, const ::Spacy::Vector&)> solverCreator,
                  const Spacy::Vector x)
        : OperatorBase(cast_ref<OperatorCreator>(space.creator()).domain(),
                       cast_ref<OperatorCreator>(space.creator()).range()),
          VectorBase(space) ,
          A_(A),
          spaces_( extractSpaces<AnsatzVariableSetDescription>(cast_ref<OperatorCreator>(space.creator()).domain()) ) ,
          solverCreator_(std::move(solverCreator)),
          x_(std::move(x))
      {  }

      /// Compute \f$A(x)\f$.
      ::Spacy::Vector operator()(const ::Spacy::Vector& x) const
      {
        Domain x_( AnsatzVariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_) );
        copyToCoefficientVector<AnsatzVariableSetDescription>(x,x_);
        Range y_( TestVariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_) );

        A_->apply(x_,y_);

        auto y = range().zeroVector();
        copyFromCoefficientVector<TestVariableSetDescription>(y_,y);

        return y;
      }

      /// Access solver representing \f$A^{-1}\f$.
      auto solver() const
      {
        return solverCreator_(*this, x_);
        //        return DirectSolver<OperatorImpl,AnsatzVariableSetDescription,TestVariableSetDescription>( A_ , spaces_, range() , domain() );
      }

      auto& get()
      {
        return A_->get_non_const();
      }

      const auto& get() const
      {
        return A_->template get<Matrix>();
      }

//      double operator()(const ::Spacy::Vector& x) const
//      {
//        const auto& x_ = cast_ref< LinearOperator<AnsatzVariableSetDescription,TestVariableSetDescription> >(x);
//        assert( impl().N() == x_.impl().N() && impl().M() == x_.impl().M() );
//        auto iend = end(impl());
//        auto iter =begin(impl()), x_iter(x_.impl());
//        auto result = 0.;
//        for( ; iter < iend ; ++iter , ++x_iter )
//          result += (*iter) * (*x_iter);
//        return result;
//      }

      const auto& A() const
      {
        return *A_;
      }

    private:
      std::shared_ptr<OperatorImpl> A_;
      Spaces spaces_;
      std::function<LinearSolver(const LinearOperator&, const ::Spacy::Vector& )> solverCreator_ = [](const LinearOperator& M, const ::Spacy::Vector& x)
        {
            return makeDirectSolver<TestVariableSetDescription,AnsatzVariableSetDescription>( M.A() ,
                                                                                              M.range() ,
                                                                                              M.domain() /*,
                                                                                              DirectType::MUMPS ,
                                                                                              MatrixProperties::GENERAL*/ );
        };

      Spacy::Vector x_;

    };



    /**
     * @brief Convenient generation of a linear operator for %Kaskade 7.
     * @param A operator from %Kaskade 7
     * @param domain domain space
     * @param range range space
     * @tparam OperatorImpl %Kaskade 7 operator, i.e. %Kaskade::AssembledGalerkinOperator or %Kaskade::MatrixRepresentedOperator.
     * @tparam AnsatzVariableSetDescription %Kaskade::VariableSetDescription for ansatz variables
     * @tparam TestVariableSetDescription %Kaskade::VariableSetDescription for test variables
     * @return LinearOperator<OperatorImpl, AnsatzVariableSetDescription, TestVariableSetDescription>( A , domain , range )
     */
    template <class AnsatzVariableSetDescription, class TestVariableSetDescription, class OperatorImpl>
    auto makeLinearOperator(const OperatorImpl& A, const VectorSpace& domain, const VectorSpace& range)
    {
      return LinearOperator<AnsatzVariableSetDescription, TestVariableSetDescription>( A , domain , range );
    }
  }
  /** @} */
}

#endif // SPACY_ADAPTER_KASKADE_LINEAR_OPERATOR_HH
