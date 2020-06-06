// This file was automatically generated using clang-type-erase.
// Please do not modify.

#pragma once

#include <Spacy/Util/SmartPointerStorage.h>
#include <memory>
#include <type_traits>
#include <iterator>

namespace Spacy
{
    class ForwardIterator
    {
        struct Interface
        {
            virtual ~Interface() = default;
            virtual std::shared_ptr< Interface > clone() const = 0;
            virtual void increment() = 0;
            virtual ForwardIterator increment_int() = 0;
            virtual double& dereference() const = 0;
            virtual bool
            compare_const_ForwardIterator_ref( const ForwardIterator& other ) const = 0;
        };

        template < class Impl >
        struct Wrapper : Interface
        {
            template < class T >
            Wrapper( T&& t ) : impl( std::forward< T >( t ) )
            {
            }

            std::shared_ptr< Interface > clone() const override
            {
                return std::make_shared< Wrapper< Impl > >( impl );
            }

            void increment() override
            {
                ++impl;
            }

            ForwardIterator increment_int() override
            {
                return impl++;
            }

            double& dereference() const override
            {
                return impl.operator*();
            }

            bool compare_const_ForwardIterator_ref( const ForwardIterator& other ) const override
            {
                return impl.operator==(
                    *other.template target< typename std::decay< Impl >::type >() );
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
        using difference_type = std::size_t;
        using value_type = double;
        using reference = double&;
        using pointer = double*;
        using iterator_category = std::forward_iterator_tag;

        ForwardIterator() noexcept = default;

        template <
            class T,
            typename std::enable_if<
                !std::is_same< typename std::decay< T >::type, ForwardIterator >::value &&
                !std::is_base_of< Interface, typename std::decay< T >::type >::value >::type* =
                nullptr >
        ForwardIterator( T&& value ) : impl_( std::forward< T >( value ) )
        {
        }

        ForwardIterator& operator++()
        {
            assert( impl_ );
            impl_->increment();
            return *this;
        }

        ForwardIterator operator++( int )
        {
            assert( impl_ );
            return impl_->increment_int();
        }

        double& operator*() const
        {
            assert( impl_ );
            return impl_->dereference();
        }

        bool operator==( const ForwardIterator& other ) const
        {
            assert( impl_ );
            return impl_->compare_const_ForwardIterator_ref( other );
        }

        template <
            class T,
            typename std::enable_if<
                !std::is_same< typename std::decay< T >::type, ForwardIterator >::value &&
                !std::is_base_of< Interface, typename std::decay< T >::type >::value >::type* =
                nullptr >
        ForwardIterator& operator=( T&& value )
        {
            return *this = ForwardIterator( std::forward< T >( value ) );
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
        clang::type_erasure::polymorphic::SBOStorage< Interface, Wrapper, 16 > impl_;
    };
    class ConstForwardIterator
    {
        struct Interface
        {
            virtual ~Interface() = default;
            virtual std::shared_ptr< Interface > clone() const = 0;
            virtual void increment() = 0;
            virtual ConstForwardIterator increment_int() = 0;
            virtual const double& dereference() const = 0;
            virtual bool
            compare_const_ConstForwardIterator_ref( const ConstForwardIterator& other ) const = 0;
        };

        template < class Impl >
        struct Wrapper : Interface
        {
            template < class T >
            Wrapper( T&& t ) : impl( std::forward< T >( t ) )
            {
            }

            std::shared_ptr< Interface > clone() const override
            {
                return std::make_shared< Wrapper< Impl > >( impl );
            }

            void increment() override
            {
                ++impl;
            }

            ConstForwardIterator increment_int() override
            {
                return impl++;
            }

            const double& dereference() const override
            {
                return impl.operator*();
            }

            bool compare_const_ConstForwardIterator_ref(
                const ConstForwardIterator& other ) const override
            {
                return impl.operator==(
                    *other.template target< typename std::decay< Impl >::type >() );
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
        using difference_type = std::size_t;
        using value_type = double;
        using reference = const double&;
        using pointer = const double*;
        using iterator_category = std::forward_iterator_tag;

        ConstForwardIterator() noexcept = default;

        template <
            class T,
            typename std::enable_if<
                !std::is_same< typename std::decay< T >::type, ConstForwardIterator >::value &&
                !std::is_base_of< Interface, typename std::decay< T >::type >::value >::type* =
                nullptr >
        ConstForwardIterator( T&& value ) : impl_( std::forward< T >( value ) )
        {
        }

        ConstForwardIterator& operator++()
        {
            assert( impl_ );
            impl_->increment();
            return *this;
        }

        ConstForwardIterator operator++( int )
        {
            assert( impl_ );
            return impl_->increment_int();
        }

        const double& operator*() const
        {
            assert( impl_ );
            return impl_->dereference();
        }

        bool operator==( const ConstForwardIterator& other ) const
        {
            assert( impl_ );
            return impl_->compare_const_ConstForwardIterator_ref( other );
        }

        template <
            class T,
            typename std::enable_if<
                !std::is_same< typename std::decay< T >::type, ConstForwardIterator >::value &&
                !std::is_base_of< Interface, typename std::decay< T >::type >::value >::type* =
                nullptr >
        ConstForwardIterator& operator=( T&& value )
        {
            return *this = ConstForwardIterator( std::forward< T >( value ) );
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
        clang::type_erasure::polymorphic::SBOStorage< Interface, Wrapper, 16 > impl_;
    };
    inline bool operator!=( const ForwardIterator& lhs, const ForwardIterator& rhs )
    {
        return !( lhs == rhs );
    }

    inline bool operator!=( const ConstForwardIterator& lhs, const ConstForwardIterator& rhs )
    {
        return !( lhs == rhs );
    }

} // namespace Spacy
