// Copyright (C) 2015 by Lars Lubkoll. All rights reserved.
// Released under the terms of the GNU General Public License version 3 or later.

///---------- Originalversion in Algorithm/CG, hier nur kleinere Anpassungen für PPCG (output, ...) -------------

#include "blockDiagPrecond2.hh"

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
    namespace MINRES
    {
        BlockDiagPreconditioner2::BlockDiagPreconditioner2(::Spacy::LinearSolver stateSolver, ::Spacy::LinearSolver controlSolver,
            ::Spacy::LinearSolver adjointSolver1, LinearSolver adjointSolver2, LinearSolver adjointSolver3,
            const VectorSpace& domain, const VectorSpace& range )
            : OperatorBase( domain, range ), stateSolver_( std::move( stateSolver ) ),
              controlSolver_( std::move( controlSolver ) ),
              adjointSolver1_( std::move( adjointSolver1 ) ),  adjointSolver2_( std::move( adjointSolver2 ) ),  adjointSolver3_( std::move( adjointSolver3 ) )

        {
           // std::cout << "TEST: "  << std::endl;
        }



        Vector BlockDiagPreconditioner2::operator()( const Vector& x) const
        {

            auto y = zero( domain() );
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
            const auto localStateIndex   = creator< ProductSpace::VectorCreator >( y_prim.space() ).idMap( getStateIndex() );

           // std::cout << "TestP 1: " << std::endl;
            y_prim.component( localStateIndex ) = stateSolver_( x_prim.component( localStateIndex ) );
           // std::cout << "TestP 2: " << std::endl;
//            if ( is< ::Spacy::PPCG::ChebyshevPreconditioner >( adjointSolver_ ) ){
//                logIterationsChebT(::Spacy::cast_ref< ::Spacy::PPCG::ChebyshevPreconditioner >( adjointSolver_ ).getIterations());
//                chebTIteraions_ = ::Spacy::cast_ref< ::Spacy::PPCG::ChebyshevPreconditioner >( adjointSolver_ ).getIterations();
//            }




            y_prim.component( localControlIndex ) = controlSolver_( x_prim.component( localControlIndex ) );


            y_dual.component( localAdjointIndex ) = adjointSolver1_( x_dual.component( localAdjointIndex ) );

           std::cout << "Check Component1: " << y_dual.component( localAdjointIndex )(y_dual.component( localAdjointIndex )) << std::endl;


            y_dual.component( localAdjointIndex ) = adjointSolver2_(  y_dual.component( localAdjointIndex ));

             std::cout << "Check Component2: " << y_dual.component( localAdjointIndex )(y_dual.component( localAdjointIndex )) << std::endl;

            y_dual.component( localAdjointIndex ) = adjointSolver3_(  y_dual.component( localAdjointIndex ));

             std::cout << "Check Component3: " << y_dual.component( localAdjointIndex )(y_dual.component( localAdjointIndex )) << std::endl;

         //  y_dual.component( localAdjointIndex ) = adjointSolver3_(adjointSolver2_(adjointSolver1_( x_dual.component( localAdjointIndex ))));


           std::cout << " TestComponents: " << y_prim.component( localStateIndex )(y_prim.component( localStateIndex )) << " " << y_prim.component( localControlIndex )(y_prim.component( localControlIndex )) << " " << y_dual.component( localAdjointIndex )(y_dual.component( localAdjointIndex )) << std::endl;
           std::cout << " TestComponentsRes: " << x_prim.component( localStateIndex )(x_prim.component( localStateIndex )) << " " << x_prim.component( localControlIndex )(x_prim.component( localControlIndex )) << " " << x_dual.component( localAdjointIndex )(x_dual.component( localAdjointIndex )) << std::endl;
              return y;
            //Warum?
            //x_dual.component( localAdjointIndex ) = cast_ref< ProductSpace::Vector >( cast_ref< ProductSpace::Vector >( x ).component( DUAL ) ).component( localAdjointIndex );
        }

//        void BlockDiagPreconditioner::setCG(::Spacy::CG::LinearSolver CgA, ::Spacy::CG::LinearSolver CgAT)
//        {
//            CgA_ = std::move(CgA);
//            CgAT_ = std::move(CgAT);
//        }


    }
}
