// This file was automatically generated using clang-type-erase.
// Please do not modify.

#pragma once

#include <Spacy/Adaptivity/SpatialAdaptivity.h>
#include <Spacy/Operator.h>
#include <Spacy/Util/SmartPointerStorage.h>
#include <Spacy/Vector.h>
#include <memory>
#include <type_traits>

namespace Spacy
{
    class ErrorEstimator
    {
        struct Interface
        {
            virtual ~Interface() = default;
            virtual std::unique_ptr< Interface > clone() const = 0;
            virtual bool estimateError( const Vector& x, const Vector& dx ) = 0;
            virtual ErrorIndicator getErrorIndicator( const Vector& errorEstimate ) = 0;
        };

        template < class Impl >
        struct Wrapper : Interface
        {
            template < class T >
            Wrapper( T&& t ) : impl( std::forward< T >( t ) )
            {
            }

            std::unique_ptr< Interface > clone() const override
            {
                return std::make_unique< Wrapper< Impl > >( impl );
            }

            bool estimateError( const Vector& x, const Vector& dx ) override
            {
                return impl.estimateError( x, dx );
            }

            ErrorIndicator getErrorIndicator( const Vector& errorEstimate ) override
            {
                return impl.getErrorIndicator( errorEstimate );
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
        ErrorEstimator() noexcept = default;

        template <
            class T,
            typename std::enable_if<
                !std::is_same< typename std::decay< T >::type, ErrorEstimator >::value &&
                !std::is_base_of< Interface, typename std::decay< T >::type >::value >::type* =
                nullptr >
        ErrorEstimator( T&& value ) : impl_( std::forward< T >( value ) )
        {
        }

        bool estimateError( const Vector& x, const Vector& dx )
        {
            assert( impl_ );
            return impl_->estimateError( x, dx );
        }

        ErrorIndicator getErrorIndicator( const Vector& errorEstimate )
        {
            assert( impl_ );
            return impl_->getErrorIndicator( errorEstimate );
        }

        template <
            class T,
            typename std::enable_if<
                !std::is_same< typename std::decay< T >::type, ErrorEstimator >::value &&
                !std::is_base_of< Interface, typename std::decay< T >::type >::value >::type* =
                nullptr >
        ErrorEstimator& operator=( T&& value )
        {
            return *this = ErrorEstimator( std::forward< T >( value ) );
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
        const T* target() const noexcept
        {
            return impl_.template target< T >();
        }

    private:
        clang::type_erasure::polymorphic::Storage< Interface, Wrapper > impl_;
    };
} // namespace Spacy
