#include "Derivative.h"

#include <Spacy/LinearOperator.h>
#include <Spacy/Spaces/ScalarSpace/Real.h>
#include <Spacy/Util/Base/OperatorBase.h>
#include <Spacy/VectorSpace.h>

#include <utility>

namespace Spacy
{
    namespace
    {
        class C2FunctionalDerivative : public OperatorBase
        {
        public:
            C2FunctionalDerivative( C2Functional&& f ) : OperatorBase( f.domain(), f.domain().dualSpace() ), f_( std::move( f ) )
            {
            }

            [[nodiscard]] Vector operator()( const Vector& x ) const
            {
                return f_.d1( x );
            }

            [[nodiscard]] Vector d1( const Vector& x, const Vector& dx ) const
            {
                return f_.d2( x, dx );
            }

            [[nodiscard]] LinearOperator linearization( const Vector& x ) const
            {
                return f_.hessian( x );
            }

        private:
            C2Functional f_;
        };
    } // namespace

    C1Operator derivative( C2Functional f )
    {
        return C2FunctionalDerivative( std::move( f ) );
    }
} // namespace Spacy
