// Copyright (C) 2015 by Lars Lubkoll. All rights reserved.
// Released under the terms of the GNU General Public License version 3 or later.

///---------- Originalversion in Algorithm/CG, hier nur kleinere Anpassungen für PPCG (output, ...) -------------

#include "TriangConstraintPrec.h"

#include <Spacy/VectorSpace.h>
#include <Spacy/ZeroVectorCreator.h>
#include <Spacy/Util/Logger.h>
#include <Spacy/Util/Cast.h>

#include <Spacy/Algorithm/Chebyshev/Chebyshev.h>

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
            ::Spacy::LinearSolver stateSolver, ::Spacy::LinearSolver controlSolver,
            ::Spacy::LinearSolver adjointSolver, OperatorWithTranspose minusB,
            const VectorSpace& domain,
            const VectorSpace& stateSpace)
            : OperatorBase( domain, domain ), stateSolver_( std::move( stateSolver ) ),
              controlSolver_( std::move( controlSolver ) ),
              adjointSolver_( std::move( adjointSolver ) ), minusB_( (minusB) ),
              stateSpace_(stateSpace),controlSpace_(minusB.domain()),adjointSpace_(minusB.range())
        {
        }

        Vector TriangularConstraintPreconditioner::operator()(Vector x ) const
        {
            auto out = zero( domain() );

            apply(out, std::move(x) );

            return out;
        }

        void TriangularConstraintPreconditioner::setSpectralBounds(double eigMin, double eigMax) 
        {
            cast_ref<::Spacy::Chebyshev::ChebyshevPreconditioner>(adjointSolver_).setSpectralBounds(eigMin,eigMax);
            cast_ref<::Spacy::Chebyshev::ChebyshevPreconditioner>(stateSolver_).setSpectralBounds(eigMin,eigMax);
        }

        
        void TriangularConstraintPreconditioner::apply( Vector& out, Vector&& in ) const
        {
                     
            out = adjointSpace_.embed(adjointSolver_( stateSpace_.project(in)));
            
            if(is< ::Spacy::CG::Solver >( adjointSolver_ ) && ::Spacy::cast_ref<::Spacy::CG::Solver >(adjointSolver_).indefiniteOperator())
            {
                signalConvex_(true,false);
                std::cout << "Fail: AT in Preconditioner indefinite!!!!!!!!!!!!!!!" << std::endl;
            }
            
            in -= minusB_.transposed( adjointSpace_.project(out));
            
            auto u = controlSolver_( controlSpace_.project(in));
            sigma_ = u( controlSpace_.project(in));
                        
            out+=u;
            
            in -= minusB_( u );
            out += stateSpace_.embed(stateSolver_( adjointSpace_.project(in) ));
            
            if(is< ::Spacy::CG::Solver >( stateSolver_ ) && ::Spacy::cast_ref<::Spacy::CG::Solver >(stateSolver_).indefiniteOperator())
            {
                signalConvex_(true,false);
                std::cout << "Fail: A in Preconditioner indefinite!!!!!!!!!!!!!!!" << std::endl;
            }
            
           iterationsInLifeTime_ += cast_ref<::Spacy::Chebyshev::ChebyshevPreconditioner>(adjointSolver_).getIterations()
                                  +cast_ref<::Spacy::Chebyshev::ChebyshevPreconditioner>(stateSolver_).getIterations();
        }


    }
}
