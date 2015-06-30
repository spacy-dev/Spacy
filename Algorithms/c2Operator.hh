#ifndef ALGORITHM_C2_OPERATOR_HH
#define ALGORITHM_C2_OPERATOR_HH

#include <utility>

#include "functionSpaceElement.hh"
#include "linearOperator.hh"
#include "c1Operator.hh"

namespace Algorithm
{
  class AbstractC2Operator;
  class FunctionSpace;

  class C2Operator : public C1Operator
  {
  public:
    C2Operator(std::shared_ptr<AbstractC2Operator> impl);

    FunctionSpaceElement d2(const FunctionSpaceElement& x, const FunctionSpaceElement& dx, const FunctionSpaceElement& dy) const;
  };

  template <class Implementation, class... Args>
  C2Operator makeC2Operator(Args&&... args)
  {
    return C2Operator( std::make_shared<Implementation>(std::forward<Args>(args)...) );
  }
}
#endif // ALGORITHM_C2_OPERATOR_HH
