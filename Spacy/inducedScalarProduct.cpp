#include "inducedScalarProduct.hh"

#include <Spacy/Spaces/ProductSpace/vector.hh>
#include <Spacy/Util/cast.hh>

#include <utility>

namespace Spacy
{
    InducedScalarProduct::InducedScalarProduct( LinearOperator M ) : M_( std::move( M ) )
    {
    }

    Real InducedScalarProduct::operator()( const Vector& x, const Vector& y ) const
    {
        return M_( y )( x );
    }

    PrimalInducedScalarProduct::PrimalInducedScalarProduct( LinearOperator M )
        : M_( std::move( M ) )
    {
    }

    Real PrimalInducedScalarProduct::operator()( Vector x, Vector y ) const
    {
        auto& x_ = cast_ref< ProductSpace::Vector >( x );
        auto& y_ = cast_ref< ProductSpace::Vector >( y );
        x_.component( DUAL ) *= 0;
        y_.component( DUAL ) *= 0;

        return M_( y )( x );
    }
}
