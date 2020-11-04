#include "LipschitzConstant.h"

namespace Spacy
{
    LipschitzConstant::LipschitzConstant( Real initialOmega ) : Mixin::Get< Real >( initialOmega )
    {
    }

    LipschitzConstant& LipschitzConstant::operator=( Real newOmega )
    {
        previousOmega_ = get();
        if ( newOmega < 0 )
        {
            newOmega = minFactor_ * previousOmega_;
        }

        get() = min( max( newOmega, eps() ), previousOmega_ * maxFactor_ );
        return *this;
    }

    Real LipschitzConstant::previous() const
    {
        return previousOmega_;
    }

    void LipschitzConstant::setMaxFactor( Real factor )
    {
        maxFactor_ = factor;
    }

    void LipschitzConstant::setMinFactor( Real factor )
    {
        minFactor_ = factor;
    }

    Real operator*( const LipschitzConstant& x, Real y )
    {
        y *= get( x );
        return y;
    }

    Real operator*( Real x, const LipschitzConstant& y )
    {
        x *= get( y );
        return x;
    }
} // namespace Spacy
