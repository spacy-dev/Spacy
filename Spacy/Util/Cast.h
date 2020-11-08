#pragma once

namespace Spacy
{
    /// Check if x can be cast to a reference of type ToType.
    template < class ToType, class FromType >
    bool is( const FromType& x )
    {
        return x.template target< ToType >() != nullptr;
    }

    /// cast x of type 'FromType&' to 'ToType&'
    template < class ToType, class FromType >
    ToType& cast_ref( FromType& x )
    {
        return *x.template target< ToType >();
    }

    /// cast x of type 'const FromType&' to 'const ToType&'
    template < class ToType, class FromType >
    const ToType& cast_ref( const FromType& x )
    {
        return *x.template target< ToType >();
    }
} // namespace Spacy
