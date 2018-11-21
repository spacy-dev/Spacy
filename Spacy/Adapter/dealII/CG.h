#pragma once

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
// For boundary values
#include <deal.II/numerics/matrix_tools.h>

#include <Spacy/Util/Base/OperatorBase.h>
#include <Spacy/Util/Cast.h>
#include <Spacy/Util/Mixins/Accuracy.h>
#include <Spacy/Util/Mixins/MaxSteps.h>
#include <Spacy/Vector.h>
#include <Spacy/VectorSpace.h>
#include <Spacy/ZeroVectorCreator.h>

#include "Util.h"
#include "Vector.h"
#include "VectorSpace.h"

namespace Spacy
{
    /** @addtogroup dealIIGroup @{ */
    namespace dealII
    {
        template < int dim, class VariableDims >
        class CGSolver : public OperatorBase, public Mixin::AbsoluteAccuracy, public Mixin::MaxSteps
        {
        public:
            CGSolver( const dealii::SparseMatrix< double >& A, const VectorSpace& domain,
                      const VectorSpace& range )
                : OperatorBase( range, domain ), A_( A.get_sparsity_pattern() )
            {
                A_.copy_from( A );
            }

            CGSolver( const CGSolver& other )
                : OperatorBase( other ), Mixin::AbsoluteAccuracy(), Mixin::MaxSteps(),
                  A_( other.A_.get_sparsity_pattern() )
            {
                checkSpaceCompatibility( domain(), other.domain() );
                checkSpaceCompatibility( range(), other.range() );

                A_.copy_from( other.A_ );
            }

            CGSolver& operator=( const CGSolver& other )
            {
                checkSpaceCompatibility( domain(), other.domain() );
                checkSpaceCompatibility( range(), other.range() );

                A_.copy_from( other.A_ );
                return *this;
            }

            ::Spacy::Vector operator()( Spacy::Vector x ) const
            {
                auto& x_ = cast_ref< Vector >( x );

                auto y = zero( range() );
                auto& y_ = cast_ref< Vector >( y );

                const auto& creator =
                    Spacy::creator< VectorCreator< dim, VariableDims::value > >( domain() );
                dealii::MatrixTools::apply_boundary_values(
                    get_boundary_map< dim, VariableDims >( domain(), creator.dofHandler() ), A_,
                    get( y_ ), get( x_ ) );

                dealii::SolverControl params( getMaxSteps(), get( getAbsoluteAccuracy() ) );
                dealii::SolverCG<> solver( params );
                solver.solve( A_, get( y_ ), get( x_ ), dealii::PreconditionIdentity() );
                return y;
            }

        private:
            mutable dealii::SparseMatrix< double > A_;
        };
    }
    /** @} */
}
