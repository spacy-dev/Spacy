#pragma once

#include <deal.II/lac/sparse_direct.h>
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

#include <map>
#include <memory>

namespace Spacy
{
    /** @addtogroup dealIIGroup @{ */
    namespace dealII
    {
        template < int dim, class VariableDims >
        class UMFPackSolver : public OperatorBase,
                              public Mixin::AbsoluteAccuracy,
                              public Mixin::MaxSteps
        {
            using Matrix = typename Traits< VariableDims >::Matrix;

        public:
            UMFPackSolver( const Matrix& A, const VectorSpace& domain, const VectorSpace& range )
                : OperatorBase( range, domain ), A_( A.get_sparsity_pattern() ),
                  solver_( std::make_shared< dealii::SparseDirectUMFPACK >() )
            {
                A_.copy_from( A );
                solver_->initialize( A_ );
            }

            UMFPackSolver( const UMFPackSolver& other )
                : OperatorBase( other ), Mixin::AbsoluteAccuracy(), Mixin::MaxSteps(),
                  A_( other.A_.get_sparsity_pattern() )
            {
                checkSpaceCompatibility( domain(), other.domain() );
                checkSpaceCompatibility( range(), other.range() );

                A_.copy_from( other.A_ );
            }

            UMFPackSolver& operator=( const UMFPackSolver& other )
            {
                checkSpaceCompatibility( domain(), other.domain() );
                checkSpaceCompatibility( range(), other.range() );

                A_.copy_from( other.A_ );
                solver_ = other.solver_;
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

                solver_->vmult( get( y_ ), get( x_ ) );

                return y + ;
            }

        private:
            mutable Matrix A_;
            std::shared_ptr< dealii::SparseDirectUMFPACK > solver_;
        };
    }
    /** @} */
}
