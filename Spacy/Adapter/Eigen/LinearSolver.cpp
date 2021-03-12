#include "LinearSolver.h"

#include "Copy.h"
#include "Vector.h"

#include <Spacy/Operator.h>
#include <Spacy/Util/Cast.h>
#include <Spacy/Vector.h>
#include <Spacy/VectorSpace.h>
#include <Spacy/ZeroVectorCreator.h>

#include <utility>

namespace Spacy::Rn
{
    LinearSolver::LinearSolver( ::Eigen::MatrixXd A, const VectorSpace& domain ) : A_( std::move( A ) ), domain_( domain )
    {
    }

    Spacy::Vector LinearSolver::operator()( const ::Spacy::Vector& y ) const
    {
        ::Eigen::VectorXd y_;
        copy( y, y_ );
        ::Eigen::VectorXd x_ = A_.lu().solve( y_ );
        ::Spacy::Vector x = zero( domain_ );
        copy( x_, x );
        return x;
    }
} // namespace Spacy::Rn
