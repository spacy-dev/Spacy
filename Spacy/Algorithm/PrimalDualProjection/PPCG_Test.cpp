#include "PPCG_Test.h"

#include <Spacy/ZeroVectorCreator.h>
#include <iostream>

namespace Spacy
{
    namespace ModifiedPPCG
    {
        PPCG_Test::PPCG_Test(  const ::Spacy::SolverBase& stateSolver,
                    const Operator& controlSolver,
                    const ::Spacy::SolverBase& adjointSolver,
                    const Operator& My, const Operator& Mu,
                    const OperatorWithTranspose& minusB ) :
                    stateSolver_( stateSolver ),
                    controlSolver_( controlSolver ),
                    adjointSolver_( adjointSolver ), My_(My), Mu_(Mu), minusB_( (minusB) ),
                    stateSpace_(stateSolver.range()),controlSpace_(controlSolver.domain()),adjointSpace_(adjointSolver.range()) 
                    {
                        yl = zero(stateSpace_);
                        ul = zero(controlSpace_);
                        pl = zero(adjointSpace_);
                    }
        
        void PPCG_Test::operator()(const Vector& dy, const Vector& du, const Vector& dp, const Vector& ru) const
        {
            yl = dy; ul = du; pl = dp;
            
            auto vSAE   = adjointSolver_(My_(dy)) + dp; // M_y(dy) + A^*(dp) = 0
            auto errSAE = sqrt(vSAE(vSAE));
            std::cout << "Error Surrogate Adjoint Equation: " << errSAE << std::endl;
            auto vSC   = Mu_(du) + minusB_.transposed(dp) + ru; // M_u(du) - B^*(dp) = -ru
            auto errSC = sqrt(vSC(vSC));
            std::cout << "Error Stationary Condition: " << errSC << std::endl;
            auto vSSE   = dy + stateSolver_(minusB_(du)); // A(dy) - B(du) = 0
            auto errSSE = sqrt(vSSE(vSSE));
            std::cout << "Error Surrogate State Equation: " << errSSE << std::endl;
            std::cout << std::endl;
        }
        
        void PPCG_Test::diffVec(const Vector& dy, const Vector& du, const Vector& dp) const
        {
            auto yd = dy - yl;
            auto ud = du - ul;
            auto pd = dp - pl;
            
            std::cout << "Difference State: " << sqrt(yd(yd)) << std::endl;
            std::cout << "Difference Control: " << sqrt(ud(ud)) << std::endl;
            std::cout << "Difference Adjoint: " << sqrt(pd(pd)) << std::endl;
        }
    }
}
    
    
    
