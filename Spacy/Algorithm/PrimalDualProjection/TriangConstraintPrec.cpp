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
            auto y = zero( domain() );

            apply( x, y, x );

            return y;
        }

        Vector TriangularConstraintPreconditioner::operator()( const Vector& x,
                                                                    Vector& q ) const
        {
            std::cout << "Fail TriangularStateConstraintPreconditioner operator() should not be called: " << std::endl;
            auto y = zero(domain());
            q = x;

            return y;
        }

        void TriangularConstraintPreconditioner::setSpectralBounds(double eigMin, double eigMax) 
        {
            cast_ref<::Spacy::Chebyshev::ChebyshevPreconditioner>(adjointSolver_).setSpectralBounds(eigMin,eigMax);
            cast_ref<::Spacy::Chebyshev::ChebyshevPreconditioner>(stateSolver_).setSpectralBounds(eigMin,eigMax);
        }

        
        void TriangularConstraintPreconditioner::apply( const Vector& x, Vector& y,
                                                             Vector& q ) const
        {
                     
            y += adjointSpace_.embed(adjointSolver_( stateSpace_.project(q)));
            
            if(is< ::Spacy::CG::Solver >( adjointSolver_ ) && ::Spacy::cast_ref<::Spacy::CG::Solver >(adjointSolver_).indefiniteOperator())
            {
                signalConvex_(true,false);
                std::cout << "Fail: AT in Preconditioner indefinite!!!!!!!!!!!!!!!" << std::endl;
            }
            
            q -= minusB_.transposed( adjointSpace_.project(y));
            
            auto u = controlSolver_( controlSpace_.project(q));
            sigma_ = u( controlSpace_.project(q));
                        
            y+=u;
            
            q -= minusB_( controlSpace_.project(y) );
            y += stateSpace_.embed(stateSolver_( adjointSpace_.project(q) ));
            
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
