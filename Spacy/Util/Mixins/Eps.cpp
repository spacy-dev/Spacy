#include "Eps.h"

namespace Spacy
{
    namespace Mixin
    {

        Eps::Eps( Real eps ) noexcept : eps_( eps )
        {
            set_eps( eps );
        }

        void Eps::set_eps( Real eps )
        {
            eps_ = eps;
            sqrt_eps_ = sqrt( eps_ );
            cbrt_eps_ = cbrt( eps_ );
            notify();
        }

        Real Eps::eps() const noexcept
        {
            return eps_;
        }

        Real Eps::sqrt_eps() const noexcept
        {
            return sqrt_eps_;
        }

        Real Eps::cbrt_eps() const noexcept
        {
            return cbrt_eps_;
        }

        void Eps::update( Eps* changedSubject )
        {
            set_eps( changedSubject->eps() );
        }
    }
}
