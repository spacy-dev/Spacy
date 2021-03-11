#include "OperatorSpace.h"

#include <Spacy/Operator.h>
#include <Spacy/VectorSpace.h>
#include <Spacy/Util/Exceptions.h>

#include "LinearOperator.h"

#include <dolfin.h>

namespace Spacy
{
    namespace FEniCS
    {
        LinearOperatorCreator::LinearOperatorCreator(
            const VectorSpace& domain, const VectorSpace& range,
            std::shared_ptr< const dolfin::FunctionSpace > dolfinSpace )
            : domain_( domain ), range_( range ), dolfinSpace_( dolfinSpace )
        {
        }

        LinearOperator LinearOperatorCreator::operator()( const VectorSpace* space ) const
        {
            throw Exception::CallOfUndefinedFunction( __func__ );
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
