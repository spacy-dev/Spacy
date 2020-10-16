// Copyright (C) 2015 by Lars Lubkoll. All rights reserved.
// Released under the terms of the GNU General Public License version 3 or later.

///---------- Originalversion in Algorithm/CG, hier nur kleinere Anpassungen für PPCG (output, ...) -------------

#include "blockDiagPrecond.hh"

#include <Spacy/Spaces/ProductSpace.h>
#include <Spacy/VectorSpace.h>
#include <Spacy/ZeroVectorCreator.h>
#include <Spacy/Util/Logger.h>

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
    namespace MINRES
    {
        BlockDiagPreconditioner::BlockDiagPreconditioner(
            const ::Spacy::CallableOperator& stateAdjSolver, const ::Spacy::CallableOperator& controlSolver, const ::Spacy::VectorSpace& stateSpace, const ::Spacy::VectorSpace& controlSpace, const ::Spacy::VectorSpace& adjointSpace)
        : stateAdjSolver_( std::move( stateAdjSolver ) ),
        controlSolver_( std::move( controlSolver ) ),
        stateSpace_(stateSpace), controlSpace_(controlSpace), adjointSpace_(adjointSpace)
        {
        }
        
        Vector BlockDiagPreconditioner::operator()( const Vector& x ) const
        {
            //return x;
            auto y = zero( x.space() );
            auto q = x;
            
            apply( x, y, q );
            
            return y;
        }
        
        Vector BlockDiagPreconditioner::operator()( const Vector& x, Vector& q ) const
        {
            std::cout << "Fail BlockDiagPreconditioner operator() should not be called: " << std::endl;
            auto y = zero( x.space() );
            q = x;
            
            return y;
        }
        
        void BlockDiagPreconditioner::apply( const Vector& x, Vector& y,   Vector& q ) const
        {
            
            y+= stateSpace_.embed(stateAdjSolver_( stateSpace_.project(x) ));
            
            y+= controlSolver_( controlSpace_.project(x) );
            
            y+= adjointSpace_.embed(stateAdjSolver_( adjointSpace_.project(x) ));
            
        }
        
        BlockDiagSchurPreconditioner::BlockDiagSchurPreconditioner(const ::Spacy::CallableOperator& stateAdjSolver, 
                              const ::Spacy::CallableOperator& stateL2Solver, 
                              const ::Spacy::CallableOperator& controlSolver, 
                              const ::Spacy::CallableOperator& stateOperator, 
                              const ::Spacy::VectorSpace& stateSpace, 
                              const ::Spacy::VectorSpace& controlSpace, 
                              const ::Spacy::VectorSpace& adjointSpace)
        : stateAdjSolver_( std::move( stateAdjSolver ) ), stateL2Solver_( std::move( stateL2Solver ) ),
        stateOperator_(stateOperator),
        controlSolver_( std::move( controlSolver ) ),
        stateSpace_(stateSpace), controlSpace_(controlSpace), adjointSpace_(adjointSpace)
        {
        }
        
        Vector BlockDiagSchurPreconditioner::operator()( const Vector& x ) const
        {
            //return x;
            auto y = zero( x.space() );
            auto q = x;
            
            apply( x, y, q );
            
            return y;
        }
        
        Vector BlockDiagSchurPreconditioner::operator()( const Vector& x, Vector& q ) const
        {
            std::cout << "Fail BlockDiagPreconditioner operator() should not be called: " << std::endl;
            auto y = zero( x.space() );
            q = x;
            
            return y;
        }
        
        void BlockDiagSchurPreconditioner::apply( const Vector& x, Vector& y,   Vector& q ) const
        {
            
            y+= stateL2Solver_( stateSpace_.project(x) );
            
            y+= controlSolver_( controlSpace_.project(x) );
            
            y+= adjointSpace_.embed(stateAdjSolver_(stateOperator_(stateAdjSolver_( adjointSpace_.project(x) ))));
            
        }
        
       
        
        
    }
}
