#include "LinearOperatorCreator.h"

#include "VectorCreator.h"
#include <Eigen/Dense>

#include <Spacy/Util/Cast.h>
#include <Spacy/ZeroVectorCreator.h>

namespace Spacy::Rn
{

    LinearOperatorCreator::LinearOperatorCreator( const VectorSpace& X, const VectorSpace& Y )
        : Generic::LinearOperatorCreator< LinearOperator >(
              [ &X, &Y ]( const VectorSpace* space ) {
                  return LinearOperator( ::Eigen::MatrixXd::Zero( cast_ref< VectorCreator >( X.creator() ).dim(),
                                                                  cast_ref< VectorCreator >( Y.creator() ).dim() ),
                                         *space, X, Y );
              },
              X, Y )
    {
    }
} // namespace Spacy::Rn
