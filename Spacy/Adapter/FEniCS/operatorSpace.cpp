#include "operatorSpace.hh"

#include "linearOperator.hh"
#include "Spacy/vectorSpace.hh"

namespace Spacy
{
  namespace FEniCS
  {
    LinearOperatorCreator::LinearOperatorCreator(const VectorSpace& domain, const VectorSpace& range, std::shared_ptr<const dolfin::FunctionSpace> dolfinSpace)
      : domain_(domain),
        range_(range),
        dolfinSpace_(dolfinSpace)
    {}

    LinearOperator LinearOperatorCreator::operator()(const VectorSpace* space) const
    {
      return LinearOperator{ dolfin::Matrix{}.copy() , *space , dolfinSpace_ };
    }

    const VectorSpace& LinearOperatorCreator::domain() const
    {
      return domain_;
    }

    const VectorSpace& LinearOperatorCreator::range() const
    {
      return range_;
    }
  }
}
