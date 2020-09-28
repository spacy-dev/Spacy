// Copyright (C) 2015 by Lars Lubkoll. All rights reserved.
// Released under the terms of the GNU General Public License version 3 or later.

///---------- Originalversion in Algorithm/CG, hier nur kleinere Anpassungen für PPCG (output, ...) -------------

#include "triangularStateConstraintPreconditioner.hh"

#include <Spacy/Spaces/ProductSpace.h>
#include <Spacy/VectorSpace.h>
#include <Spacy/ZeroVectorCreator.h>
#include <Spacy/Util/Logger.h>

#include <Spacy/Algorithm/PPCG/chebyshev.hh>

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
    namespace PPCG
    {
        TriangularStateConstraintPreconditioner::TriangularStateConstraintPreconditioner(
            ::Spacy::LinearSolver stateSolver, ::Spacy::LinearSolver controlSolver,
            ::Spacy::LinearSolver adjointSolver, CallableOperator B, CallableOperator BT,
            const VectorSpace& domain, const VectorSpace& range )
            : OperatorBase( domain, range ), stateSolver_( std::move( stateSolver ) ),
              controlSolver_( std::move( controlSolver ) ),
              adjointSolver_( std::move( adjointSolver ) ), B_( std::move( B ) ),
              BT_( std::move( BT ) )
        {
        }

        Vector TriangularStateConstraintPreconditioner::operator()( const Vector& x ) const
        {
            auto y = zero( domain() );
            auto q = x;

            apply( x, y, q );

            return y;
        }

        Vector TriangularStateConstraintPreconditioner::operator()( const Vector& x,
                                                                    Vector& q ) const
        {
            std::cout << "Fail TriangularStateConstraintPreconditioner operator() should not be called: " << std::endl;
            auto y = zero( domain() );
            q = x;

            return y;
        }

        void TriangularStateConstraintPreconditioner::apply( const Vector& x, Vector& y,
                                                             Vector& q ) const
        {
            typedef std::numeric_limits< double > dbl;


            std::cout.precision(dbl::max_digits10);

            auto x_ = cast_ref< ProductSpace::Vector >( x );
            auto& x_prim = cast_ref< ProductSpace::Vector >( x_.component( PRIMAL ) );
            auto& x_dual = cast_ref< ProductSpace::Vector >( x_.component( DUAL ) );

            auto& y_ = cast_ref< ProductSpace::Vector >( y );
            auto& y_prim = cast_ref< ProductSpace::Vector >( y_.component( PRIMAL ) );
            auto& y_dual = cast_ref< ProductSpace::Vector >( y_.component( DUAL ) );

            const auto localAdjointIndex = creator< ProductSpace::VectorCreator >( y_dual.space() ).idMap( getAdjointIndex() );
            const auto localControlIndex = creator< ProductSpace::VectorCreator >( y_prim.space() ).idMap( getControlIndex() );
            const auto localStateIndex = creator< ProductSpace::VectorCreator >( y_prim.space() ).idMap( getStateIndex() );

            y_dual.component( localAdjointIndex ) = adjointSolver_( x_prim.component( localStateIndex ) );
            if ( is< ::Spacy::PPCG::ChebyshevPreconditioner >( adjointSolver_ ) ){
                logIterationsChebT(::Spacy::cast_ref< ::Spacy::PPCG::ChebyshevPreconditioner >( adjointSolver_ ).getIterations());
                chebTIteraions_ = ::Spacy::cast_ref< ::Spacy::PPCG::ChebyshevPreconditioner >( adjointSolver_ ).getIterations();
            }
            if(is< ::Spacy::CG::Solver >( adjointSolver_ ) && ::Spacy::cast_ref<::Spacy::CG::Solver >(adjointSolver_).indefiniteOperator())
            {
                signalConvex_(true,false);
                std::cout << "Fail: AT in Preconditioner indefinite!!!!!!!!!!!!!!!" << std::endl;
            }

            x_prim.component( localControlIndex ) -= BT_( y_dual.component( localAdjointIndex ) );

            y_prim.component( localControlIndex ) = controlSolver_( x_prim.component( localControlIndex ) );

            x_dual.component( localAdjointIndex ) -= B_( y_prim.component( localControlIndex ) );

            y_prim.component( localStateIndex ) = stateSolver_( x_dual.component( localAdjointIndex ) );
            if ( is< ::Spacy::PPCG::ChebyshevPreconditioner >( stateSolver_ ) ){
                logIterationsCheb(::Spacy::cast_ref< ::Spacy::PPCG::ChebyshevPreconditioner >( stateSolver_ ).getIterations());
                chebIterations_ = ::Spacy::cast_ref< ::Spacy::PPCG::ChebyshevPreconditioner >( stateSolver_ ).getIterations();
            }
            if(is< ::Spacy::CG::Solver >( stateSolver_ ) && ::Spacy::cast_ref<::Spacy::CG::Solver >(stateSolver_).indefiniteOperator())
            {
                signalConvex_(true,false);
                std::cout << "Fail: A in Preconditioner indefinite!!!!!!!!!!!!!!!" << std::endl;
            }

            //Warum?
            //x_dual.component( localAdjointIndex ) = cast_ref< ProductSpace::Vector >( cast_ref< ProductSpace::Vector >( x ).component( DUAL ) ).component( localAdjointIndex );
        }

//        void TriangularStateConstraintPreconditioner::setCG(::Spacy::CG::LinearSolver CgA, ::Spacy::CG::LinearSolver CgAT)
//        {
//            CgA_ = std::move(CgA);
//            CgAT_ = std::move(CgAT);
//        }

        Vector TriangularStateConstraintPreconditioner::kernelOffset( const Vector& rhs ) const
        {
            auto& rhs_ = cast_ref< ProductSpace::Vector >( rhs );
            auto& rhs_dual = cast_ref< ProductSpace::Vector >( rhs_.component( DUAL ) );

            auto y = zero( range() );
            auto& y_ = cast_ref< ProductSpace::Vector >( y );
            auto& y_prim = cast_ref< ProductSpace::Vector >( y_.component( PRIMAL ) );
            const auto localStateIndex = creator< ProductSpace::VectorCreator >( y_prim.space() ).idMap( getStateIndex() );
            const auto localAdjointIndex = creator< ProductSpace::VectorCreator >( rhs_dual.space() ).idMap( getAdjointIndex() );

            y_prim.component( localStateIndex ) = stateSolver_( rhs_dual.component( localAdjointIndex ) );
            return y;
        }
    }
}
