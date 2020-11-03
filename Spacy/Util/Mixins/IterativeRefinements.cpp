#include "IterativeRefinements.h"

namespace Spacy
{
    namespace Mixin
    {
        IterativeRefinements::IterativeRefinements( unsigned value ) noexcept : value_( value )
        {
        }

        void IterativeRefinements::setIterativeRefinements( unsigned value )
        {
            value_ = value;
            notify();
        }

        unsigned IterativeRefinements::getIterativeRefinements() const noexcept
        {
            return value_;
        }

        void IterativeRefinements::update( IterativeRefinements* changedSubject )
        {
            setIterativeRefinements( changedSubject->getIterativeRefinements() );
        }
    } // namespace Mixin
} // namespace Spacy
