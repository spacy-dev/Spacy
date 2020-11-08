#include "NumberOfThreads.h"

namespace Spacy::Mixin
{
    NumberOfThreads::NumberOfThreads( unsigned value ) noexcept : value_( value )
    {
    }

    void NumberOfThreads::setNumberOfThreads( unsigned value )
    {
        value_ = value;
        notify();
    }

    unsigned NumberOfThreads::getNumberOfThreads() const noexcept
    {
        return value_;
    }

    void NumberOfThreads::update( NumberOfThreads* changedSubject )
    {
        setNumberOfThreads( changedSubject->getNumberOfThreads() );
    }
} // namespace Spacy::Mixin