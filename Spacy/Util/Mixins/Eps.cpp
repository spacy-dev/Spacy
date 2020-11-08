#include "Eps.h"

namespace Spacy::Mixin
{

    Eps::Eps( Real eps ) noexcept : eps_( eps )
    {
        setEps( eps );
    }

    void Eps::setEps( Real eps )
    {
        eps_ = eps;
        sqrtEps_ = sqrt( eps_ );
        cbrtEps_ = cbrt( eps_ );
        notify();
    }

    Real Eps::eps() const noexcept
    {
        return eps_;
    }

    Real Eps::sqrtEps() const noexcept
    {
        return sqrtEps_;
    }

    Real Eps::cbrtEps() const noexcept
    {
        return cbrtEps_;
    }

    void Eps::update( Eps* changedSubject )
    {
        setEps( changedSubject->eps() );
    }
} // namespace Spacy::Mixin
