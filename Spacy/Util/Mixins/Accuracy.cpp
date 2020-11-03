#include "Accuracy.h"

namespace Spacy
{
    namespace Mixin
    {
        AbsoluteAccuracy::AbsoluteAccuracy( Real value ) noexcept : value_( value )
        {
        }

        void AbsoluteAccuracy::setAbsoluteAccuracy( Real value )
        {
            value_ = value;
            notify();
        }

        Real AbsoluteAccuracy::getAbsoluteAccuracy() const noexcept
        {
            return value_;
        }

        void AbsoluteAccuracy::update( AbsoluteAccuracy* changedSubject )
        {
            setAbsoluteAccuracy( changedSubject->getAbsoluteAccuracy() );
        }

        RelativeAccuracy::RelativeAccuracy( Real value ) noexcept : value_( value )
        {
        }

        void RelativeAccuracy::setRelativeAccuracy( Real value )
        {
            value_ = value;
            notify();
        }

        Real RelativeAccuracy::getRelativeAccuracy() const noexcept
        {
            return value_;
        }

        void RelativeAccuracy::update( RelativeAccuracy* changedSubject )
        {
            setRelativeAccuracy( changedSubject->getRelativeAccuracy() );
        }

        MinimalAccuracy::MinimalAccuracy( Real value ) noexcept : value_( value )
        {
        }

        void MinimalAccuracy::setMinimalAccuracy( Real value )
        {
            value_ = value;
            notify();
        }

        Real MinimalAccuracy::getMinimalAccuracy() const noexcept
        {
            return value_;
        }

        void MinimalAccuracy::update( MinimalAccuracy* changedSubject )
        {
            setMinimalAccuracy( changedSubject->getMinimalAccuracy() );
        }

        DampingAccuracy::DampingAccuracy( Real value ) noexcept : value_( value )
        {
        }

        void DampingAccuracy::setDampingAccuracy( Real value )
        {
            value_ = value;
            notify();
        }

        Real DampingAccuracy::getDampingAccuracy() const noexcept
        {
            return value_;
        }

        void DampingAccuracy::update( DampingAccuracy* changedSubject )
        {
            setDampingAccuracy( changedSubject->getDampingAccuracy() );
        }

        NormalAccuracy::NormalAccuracy( Real value ) noexcept : value_( value )
        {
        }

        void NormalAccuracy::setNormalAccuracy( Real value )
        {
            value_ = value;
            notify();
        }

        Real NormalAccuracy::getNormalAccuracy() const noexcept
        {
            return value_;
        }

        void NormalAccuracy::update( NormalAccuracy* changedSubject )
        {
            setNormalAccuracy( changedSubject->getNormalAccuracy() );
        }

        TangentialAccuracy::TangentialAccuracy( Real value ) noexcept : value_( value )
        {
        }

        void TangentialAccuracy::setTangentialAccuracy( Real value )
        {
            value_ = value;
            notify();
        }

        Real TangentialAccuracy::getTangentialAccuracy() const noexcept
        {
            return value_;
        }

        void TangentialAccuracy::update( TangentialAccuracy* changedSubject )
        {
            setTangentialAccuracy( changedSubject->getTangentialAccuracy() );
        }

        FallBackTangentialAccuracy::FallBackTangentialAccuracy( Real value ) noexcept : value_( value )
        {
        }

        void FallBackTangentialAccuracy::setFallBackTangentialAccuracy( Real value )
        {
            value_ = value;
            notify();
        }

        Real FallBackTangentialAccuracy::getFallBackTangentialAccuracy() const noexcept
        {
            return value_;
        }

        void FallBackTangentialAccuracy::update( FallBackTangentialAccuracy* changedSubject )
        {
            setFallBackTangentialAccuracy( changedSubject->getFallBackTangentialAccuracy() );
        }
    } // namespace Mixin
} // namespace Spacy