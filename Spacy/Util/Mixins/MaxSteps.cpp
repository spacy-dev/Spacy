#include "MaxSteps.h"

namespace Spacy
{
    namespace Mixin
    {
        MaxSteps::MaxSteps( unsigned value ) noexcept : value_( value )
        {
        }

        void MaxSteps::setMaxSteps( unsigned value )
        {
            value_ = value;
            notify();
        }

        unsigned MaxSteps::getMaxSteps() const noexcept
        {
            return value_;
        }

        void MaxSteps::update( MaxSteps* changedSubject )
        {
            setMaxSteps( changedSubject->getMaxSteps() );
        }
    } // namespace Mixin
} // namespace Spacy