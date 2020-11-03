#include "TriangularStateConstraintPreconditioner.h"

#include <Spacy/Spaces/ProductSpace.h>
#include <Spacy/VectorSpace.h>
#include <Spacy/ZeroVectorCreator.h>

#include <utility>

namespace Spacy
{
    namespace CG
    {
        TriangularStateConstraintPreconditioner::TriangularStateConstraintPreconditioner(
            ::Spacy::LinearSolver stateSolver, ::Spacy::LinearSolver controlSolver, ::Spacy::LinearSolver adjointSolver, CallableOperator B,
            CallableOperator BT, const VectorSpace& domain, const VectorSpace& range )
            : OperatorBase( domain, range ), stateSolver_( std::move( stateSolver ) ), controlSolver_( std::move( controlSolver ) ),
              adjointSolver_( std::move( adjointSolver ) ), B_( std::move( B ) ), BT_( std::move( BT ) )
        {
        }

        Vector TriangularStateConstraintPreconditioner::operator()( const Vector& x ) const
        {
            auto y = zero( domain() );
            auto q = x;

            apply( x, y, q );

            return y;
        }

        Vector TriangularStateConstraintPreconditioner::operator()( const Vector& x, Vector& q ) const
        {
            auto y = zero( domain() );
            q = x;

            return y;
        }

        void TriangularStateConstraintPreconditioner::apply( const Vector& x, Vector& y, Vector& q ) const
        {
            auto x_ = cast_ref< ProductSpace::Vector >( x );
            auto& xp_ = cast_ref< ProductSpace::Vector >( x_.component( PRIMAL ) );
            auto& xd_ = cast_ref< ProductSpace::Vector >( x_.component( DUAL ) );

            auto& y_ = cast_ref< ProductSpace::Vector >( y );
            auto& yp_ = cast_ref< ProductSpace::Vector >( y_.component( PRIMAL ) );
            auto& yd_ = cast_ref< ProductSpace::Vector >( y_.component( DUAL ) );

            const auto localAdjointIndex = creator< ProductSpace::VectorCreator >( yd_.space() ).idMap( getAdjointIndex() );
            const auto localControlIndex = creator< ProductSpace::VectorCreator >( yp_.space() ).idMap( getControlIndex() );
            const auto localStateIndex = creator< ProductSpace::VectorCreator >( yp_.space() ).idMap( getStateIndex() );

            yd_.component( localAdjointIndex ) = adjointSolver_( xp_.component( localStateIndex ) );
            xp_.component( localControlIndex ) -= BT_( yd_.component( localAdjointIndex ) );

            yp_.component( localControlIndex ) = controlSolver_( xp_.component( localControlIndex ) );

            xd_.component( localAdjointIndex ) -= B_( yp_.component( localControlIndex ) );

            yp_.component( localStateIndex ) = stateSolver_( xd_.component( localAdjointIndex ) );

            xd_.component( localAdjointIndex ) =
                cast_ref< ProductSpace::Vector >( cast_ref< ProductSpace::Vector >( x ).component( DUAL ) ).component( localAdjointIndex );
        }

        Vector TriangularStateConstraintPreconditioner::kernelOffset( const Vector& rhs ) const
        {
            const auto& rhs_ = cast_ref< ProductSpace::Vector >( rhs );
            const auto& rhsd_ = cast_ref< ProductSpace::Vector >( rhs_.component( DUAL ) );

            auto y = zero( range() );
            auto& y_ = cast_ref< ProductSpace::Vector >( y );
            auto& yp_ = cast_ref< ProductSpace::Vector >( y_.component( PRIMAL ) );
            const auto localStateIndex = creator< ProductSpace::VectorCreator >( yp_.space() ).idMap( getStateIndex() );
            const auto localAdjointIndex = creator< ProductSpace::VectorCreator >( rhsd_.space() ).idMap( getAdjointIndex() );

            yp_.component( localStateIndex ) = stateSolver_( rhsd_.component( localAdjointIndex ) );
            return y;
        }
    } // namespace CG
} // namespace Spacy
