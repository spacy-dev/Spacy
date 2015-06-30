#include "c1Operator.hh"

#include "functionSpaceElement.hh"
#include "Interface/Operator/abstractC1Operator.hh"

namespace Algorithm
{
  C1Operator::C1Operator(std::shared_ptr<AbstractC1Operator> impl)
    : Operator(std::static_pointer_cast<AbstractOperator>(impl))
  {}

  FunctionSpaceElement C1Operator::d1(const FunctionSpaceElement& x, const FunctionSpaceElement& dx) const
  {
    return FunctionSpaceElement( dynamic_cast<const AbstractC1Operator&>( impl() ).getLinearization(x.impl())( dx.impl() ) );
  }

  LinearOperator C1Operator::getLinearization(const FunctionSpaceElement& x) const
  {
    return LinearOperator( dynamic_cast<const AbstractC1Operator&>( impl() ).getLinearization(x.impl()) );
  }
}
