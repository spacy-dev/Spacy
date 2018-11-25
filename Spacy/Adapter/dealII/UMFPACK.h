#pragma once

#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_matrix.h>
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

#include <iostream>
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
        public:
            using NativeVector = typename Traits< VariableDims >::Vector;
            using Matrix = typename Traits< VariableDims >::Matrix;

            UMFPackSolver( const Matrix& A, const VectorSpace& domain, const VectorSpace& range )
                : OperatorBase( range, domain ), A_( A.get_sparsity_pattern() )
            {
                A_.copy_from( A );
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
                return *this;
            }

            ::Spacy::Vector operator()( Spacy::Vector x ) const
            {
                auto& x_ = cast_ref< Vector >( x );

                auto y = zero( range() );
                auto& y_ = cast_ref< Vector >( y );
                setBoundaryConditions< dim, VariableDims >( domain(), y_ );

                auto boundaryMap = get_boundary_map< dim, VariableDims >(
                    domain(),
                    creator< VectorCreator< dim, VariableDims::value > >( domain() ).dofHandler() );
                dealii::MatrixTools::apply_boundary_values( boundaryMap, A_, y_.get(), x_.get() );

                dealii::SparseDirectUMFPACK solver;
                solver.initialize( A_ );
                solver.vmult( get( y_ ), get( x_ ) );

                return y;
            }

        private:
            mutable Matrix A_;
        };
    }
    /** @} */
}
