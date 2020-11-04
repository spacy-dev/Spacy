#include "VectorBase.h"

#include <Spacy/VectorSpace.h>

namespace Spacy
{
    VectorBase::VectorBase( const VectorSpace& space ) : space_( space )
    {
    }

    VectorBase::VectorBase( const VectorBase& /*y*/ ) = default;

    VectorBase::VectorBase( VectorBase&& y ) noexcept : space_( y.space_ )
    {
    }

    VectorBase& VectorBase::operator=( const VectorBase& y )
    {
        checkSpaceCompatibility( this->space(), y.space() );
        return *this;
    }

    VectorBase& VectorBase::operator=( VectorBase&& y ) noexcept
    {
        checkSpaceCompatibility( this->space(), y.space() );
        return *this;
    }

    const VectorSpace& VectorBase::space() const
    {
        return space_;
    }
} // namespace Spacy
