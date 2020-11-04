// This file was automatically generated using clang-type-erase.
// Please do not modify.

#pragma once

#include <Spacy/Util/SmartPointerStorage.h>
#include <Spacy/Vector.h>
#include <Spacy/VectorSpace.h>

#include <memory>
#include <type_traits>

namespace Spacy
{
    /// Each VectorSpace needs a zero-vector creator to support generation of vector space elements.
    class ZeroVectorCreator
    {
        struct Interface
        {
            virtual ~Interface() = default;
            [[nodiscard]] virtual std::unique_ptr< Interface > clone() const = 0;
            virtual Vector call_const_VectorSpace_ptr( const VectorSpace* V ) const = 0;
        };

        template < class Impl >
        struct Wrapper : Interface
        {
            template < class T >
            Wrapper( T&& t ) : impl( std::forward< T >( t ) )
            {
            }

            [[nodiscard]] std::unique_ptr< Interface > clone() const override
            {
                return std::make_unique< Wrapper< Impl > >( impl );
            }

            Vector call_const_VectorSpace_ptr( const VectorSpace* V ) const override
            {
                return impl.operator()( V );
            }

            Impl impl;
        };

        template < class Impl >
        struct Wrapper< std::reference_wrapper< Impl > > : Wrapper< Impl& >
        {
            template < class T >
            Wrapper( T&& t ) : Wrapper< Impl& >( std::forward< T >( t ) )
            {
            }
        };

    public:
        ZeroVectorCreator() noexcept = default;

        template < class T,
                   typename std::enable_if< !std::is_same< typename std::decay< T >::type, ZeroVectorCreator >::value &&
                                            !std::is_base_of< Interface, typename std::decay< T >::type >::value >::type* = nullptr >
        ZeroVectorCreator( T&& value ) : impl_( std::forward< T >( value ) )
        {
        }

        /// Creates \f$
        Vector operator()( const VectorSpace* V ) const
        {
            assert( impl_ );
            return impl_->call_const_VectorSpace_ptr( V );
        }

        template < class T,
                   typename std::enable_if< !std::is_same< typename std::decay< T >::type, ZeroVectorCreator >::value &&
                                            !std::is_base_of< Interface, typename std::decay< T >::type >::value >::type* = nullptr >
        ZeroVectorCreator& operator=( T&& value )
        {
            return *this = ZeroVectorCreator( std::forward< T >( value ) );
        }

        explicit operator bool() const noexcept
        {
            return bool( impl_ );
        }

        template < class T >
        T* target() noexcept
        {
            return impl_.template target< T >();
        }

        template < class T >
        [[nodiscard]] const T* target() const noexcept
        {
            return impl_.template target< T >();
        }

    private:
        clang::type_erasure::polymorphic::Storage< Interface, Wrapper > impl_;
    };

    template < class T >
    T& creator( VectorSpace& space )
    {
        return *space.creator().template target< T >();
    }

    template < class T >
    const T& creator( const VectorSpace& space )
    {
        return *space.creator().template target< T >();
    }

    /// Create new vector \f$v=0\f$.
    inline Vector zero( const VectorSpace& space )
    {
        return space.creator()( &space );
    }

} // namespace Spacy
