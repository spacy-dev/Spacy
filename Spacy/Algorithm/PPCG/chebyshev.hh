#pragma once

#include "Spacy/Operator.h"
#include "../../../Spacy/Util/Base/OperatorBase.h"
//#include "../../../Spacy/Util/Mixins/Eps.h"
//#include "../../../Spacy/Util/Mixins/maxSteps.h"
#include <Spacy/Util/Mixins.h>
#include "../../../Spacy/Util/Mixins/Get.h"

namespace Spacy
{
  /// \cond
  class VectorSpace;
  /// \endcond

  namespace PPCG
  {
    /**
     * @brief A Preconditioner based on the chebyshev-iteration solving A*x=b
     */
    class ChebyshevPreconditioner
            : //public OperatorBase,
              public Mixin::Eps,
              public Mixin::MaxSteps,
              public Mixin::Verbosity,
              public Mixin::AbsoluteAccuracy,
              public Mixin::RelativeAccuracy
    {
    public:
      /**
       * @brief Constructor.
       * @param A linear Operator
       * @param P preconditioner
       * @param gamma_max largest eigenvalue of A
       * @param gamma_min smallest eigenvalue of A

       */
      ChebyshevPreconditioner(::Spacy::CallableOperator A, ::Spacy::CallableOperator P,
                              ::Spacy::Real gamma_max, ::Spacy::Real gamma_min, std::function<::Spacy::Vector(::Spacy::Vector,::Spacy::Vector)> transfer
                              //,const VectorSpace& domain, const VectorSpace& range
                              );

      /**
       * @brief Apply preconditioner.
       * @param x argument
       */
      Vector operator()(const Vector& x) const;

      ::Spacy::Real acc_ = 1e-2; //public?

      unsigned getIterations() const
      {
          return iterations_;
      }

    private:
        ::Spacy::CallableOperator A_;
        ::Spacy::CallableOperator P_;
        ::Spacy::Real gamma_max_;
        ::Spacy::Real gamma_min_;

        mutable unsigned int iterations_ = 1;

        std::function<::Spacy::Vector(::Spacy::Vector,::Spacy::Vector)> transfer_ = [](::Spacy::Vector src, ::Spacy::Vector res){ return src;};

    };
  }
}

