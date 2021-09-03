// This file was automatically generated using clang-type-erase.
// Please do not modify.

#pragma once

#include <Spacy/Util/SmartPointerStorage.h>
#include <Spacy/Vector.h>

#include <memory>
#include <type_traits>

namespace Spacy
{
    /// Type-erased termination criterion for conjugate gradient methods.
    class SolverBase
    {
        struct Interface
        {
            virtual ~Interface() = default;
            [[nodiscard]] virtual std::unique_ptr< Interface > clone() const = 0;
            [[nodiscard]] virtual Vector call_const_Vector_ref( const Vector& x ) const = 0;
            [[nodiscard]] virtual unsigned getIterations() const = 0;
            [[nodiscard]] virtual const VectorSpace& domain() const = 0;
            [[nodiscard]] virtual const VectorSpace& range() const = 0;
            virtual void setEps( double eps ) = 0;
            virtual void setRelativeAccuracy( double accuracy ) = 0;
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

            Vector call_const_Vector_ref( const Vector& x ) const override
            {
                return impl.operator()( x );
            }
            
            unsigned getIterations() const override
            {
                return impl.getIterations();
            }

            const VectorSpace& domain() const override
            {
                return impl.domain();
            }

            const VectorSpace& range() const override
            {
                return impl.range();
            }

            void setEps( double eps ) override
            {
                impl.setEps( std::move( eps ) );
            }

            void setRelativeAccuracy( double accuracy ) override
            {
                impl.setRelativeAccuracy( std::move( accuracy ) );
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
        SolverBase() noexcept = default;

        template < class T,
                    typename std::enable_if< !std::is_same< typename std::decay< T >::type, SolverBase >::value &&
                                            !std::is_base_of< Interface, typename std::decay< T >::type >::value >::type* = nullptr >
        SolverBase( T&& value ) : impl_( std::forward< T >( value ) )
        {
        }

        /// Apply operator.
        [[nodiscard]] Vector operator()( const Vector& x ) const
        {
            assert( impl_ );
            return impl_->call_const_Vector_ref( x );
        }
        
        /// get Iterations.
        [[nodiscard]] unsigned getIterations() const
        {
            assert( impl_ );
            return impl_->getIterations();
        }
        
        /// Access domain space \f$X\f$.
        [[nodiscard]] const VectorSpace& domain() const
        {
            assert( impl_ );
            return impl_->domain();
        }

        /// Access range space \f$Y\f$.
        [[nodiscard]] const VectorSpace& range() const
        {
            assert( impl_ );
            return impl_->range();
        }

        void setEps( double eps )
        {
            assert( impl_ );
            impl_->setEps( std::move( eps ) );
        }

        void setRelativeAccuracy( double accuracy )
        {
            assert( impl_ );
            impl_->setRelativeAccuracy( std::move( accuracy ) );
        }

        template < class T,
                    typename std::enable_if< !std::is_same< typename std::decay< T >::type, SolverBase >::value &&
                                            !std::is_base_of< Interface, typename std::decay< T >::type >::value >::type* = nullptr >
        SolverBase& operator=( T&& value )
        {
            return *this = SolverBase( std::forward< T >( value ) );
        }

        [[nodiscard]] explicit operator bool() const noexcept
        {
            return bool( impl_ );
        }

        template < class T >
        [[nodiscard]] T* target() noexcept
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
} // namespace Spacy
