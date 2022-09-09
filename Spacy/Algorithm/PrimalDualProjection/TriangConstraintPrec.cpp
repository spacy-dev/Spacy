// Copyright (C) 2015 by Lars Lubkoll. All rights reserved.
// Released under the terms of the GNU General Public License version 3 or later.

///---------- Originalversion in Algorithm/CG, hier nur kleinere Anpassungen für PPCG (output, ...) -------------

#include "TriangConstraintPrec.h"

#include <Spacy/VectorSpace.h>
#include <Spacy/ZeroVectorCreator.h>
#include <Spacy/Util/Logger.h>
#include <Spacy/Util/Cast.h>

#include "Spacy/Algorithm/Preconditioner/Chebyshev.h"

#include <utility>
#include <iostream>
#include <cmath>
#include <limits>

namespace Spacy
{
namespace
{
    Logger< unsigned > logIterationsCheb( "chebIt.log" );
    Logger< unsigned > logIterationsChebT( "chebTIt.log" );
}
    namespace ModifiedPPCG
    {
        TriangularConstraintPreconditioner::TriangularConstraintPreconditioner(
            const ::Spacy::SolverBase& stateSolver, const Operator& controlSolver,
            const ::Spacy::SolverBase& adjointSolver, const OperatorWithTranspose& minusB,
            const VectorSpace& domain)
            : OperatorBase( domain, domain ), stateSolver_( stateSolver ),
              controlSolver_( controlSolver ),
              adjointSolver_( adjointSolver ), minusB_( (minusB) ),
              stateSpace_(stateSolver.range()),controlSpace_(controlSolver.domain()),adjointSpace_(minusB.range())
        {
        }

        void TriangularConstraintPreconditioner::operator()(Vector y, Vector u , Vector& outy, Vector& outu, Vector& outp) const
        {
            apply(outy, outu, outp, std::move(y), std::move(u) );
        }

        void TriangularConstraintPreconditioner::setSpectralBounds(double eigMin, double eigMax) 
        {
//             cast_ref<::Spacy::Preconditioner::Chebyshev>(adjointSolver_).setSpectralBounds(eigMin,eigMax);
//             cast_ref<::Spacy::Preconditioner::Chebyshev>(stateSolver_).setSpectralBounds(eigMin,eigMax);
        }

        
        void TriangularConstraintPreconditioner::apply( Vector& outy, Vector& outu, Vector& outp, Vector&& iny, Vector&& inu ) const
        {
            outp = adjointSolver_( iny );
            iterations_ += adjointSolver_.getIterations();
            
            inu -= minusB_.transposed( outp );
            outu = controlSolver_( inu );
            
            sigma_ = outu( inu );
            
            auto inp = -minusB_(outu);
            outy = stateSolver_( inp );
            iterations_ += stateSolver_.getIterations();
        }


    } // namespace ModifiedPPCG
} // namespace Spacy
