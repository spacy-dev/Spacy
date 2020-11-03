#include "Index.h"

namespace Spacy
{
    namespace Mixin
    {
        StateIndex::StateIndex( unsigned value ) noexcept : value_( value )
        {
        }

        void StateIndex::setStateIndex( unsigned value )
        {
            value_ = value;
            notify();
        }

        unsigned StateIndex::getStateIndex() const noexcept
        {
            return value_;
        }

        void StateIndex::update( StateIndex* changedSubject )
        {
            setStateIndex( changedSubject->getStateIndex() );
        }

        ControlIndex::ControlIndex( unsigned value ) noexcept : value_( value )
        {
        }

        void ControlIndex::setControlIndex( unsigned value )
        {
            value_ = value;
            notify();
        }

        unsigned ControlIndex::getControlIndex() const noexcept
        {
            return value_;
        }

        void ControlIndex::update( ControlIndex* changedSubject )
        {
            setControlIndex( changedSubject->getControlIndex() );
        }

        AdjointIndex::AdjointIndex( unsigned value ) noexcept : value_( value )
        {
        }

        void AdjointIndex::setAdjointIndex( unsigned value )
        {
            value_ = value;
            notify();
        }

        unsigned AdjointIndex::getAdjointIndex() const noexcept
        {
            return value_;
        }

        void AdjointIndex::update( AdjointIndex* changedSubject )
        {
            setAdjointIndex( changedSubject->getAdjointIndex() );
        }
    } // namespace Mixin
} // namespace Spacy
