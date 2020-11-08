// This file was automatically generated using clang-type-erase.
// Please do not modify.

#pragma once

#include <Spacy/LinearOperator.h>
#include <Spacy/Spaces/ScalarSpace/Real.h>
#include <Spacy/Util/SmartPointerStorage.h>
#include <Spacy/Vector.h>
#include <Spacy/VectorSpace.h>

#include <memory>
#include <type_traits>

namespace Spacy
{
    /// Type-erased twice differentiable functional \f$f:\ X \to \mathbb{R} \f$.
    class C2Functional
    {
        struct Interface
        {
            virtual ~Interface() = default;
            [[nodiscard]] virtual std::shared_ptr< Interface > clone() const = 0;
            [[nodiscard]] virtual Real call_const_Vector_ref( const Vector& x ) const = 0;
            [[nodiscard]] virtual Vector d1( const Vector& x ) const = 0;
            [[nodiscard]] virtual Vector d2( const Vector& x, const Vector& dx ) const = 0;
            [[nodiscard]] virtual LinearOperator hessian( const Vector& x ) const = 0;
            [[nodiscard]] virtual const VectorSpace& domain() const = 0;
        };

        template < class Impl >
        struct Wrapper : Interface
        {
            template < class T >
            Wrapper( T&& t ) : impl( std::forward< T >( t ) )
            {
            }

            [[nodiscard]] std::shared_ptr< Interface > clone() const override
            {
                return std::make_shared< Wrapper< Impl > >( impl );
            }

            Real call_const_Vector_ref( const Vector& x ) const override
            {
                return impl.operator()( x );
            }

            Vector d1( const Vector& x ) const override
            {
                return impl.d1( x );
            }

            Vector d2( const Vector& x, const Vector& dx ) const override
            {
                return impl.d2( x, dx );
            }

            LinearOperator hessian( const Vector& x ) const override
            {
                return impl.hessian( x );
            }

            const VectorSpace& domain() const override
            {
                return impl.domain();
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
        C2Functional() noexcept = default;

        template < class T,
                   typename std::enable_if< !std::is_same< typename std::decay< T >::type, C2Functional >::value &&
                                            !std::is_base_of< Interface, typename std::decay< T >::type >::value >::type* = nullptr >
        C2Functional( T&& value ) : impl_( std::forward< T >( value ) )
        {
        }

        /// Apply functional.
        [[nodiscard]] Real operator()( const Vector& x ) const
        {
            assert( impl_ );
            return impl_->call_const_Vector_ref( x );
        }

        /// Compute derivative as function space element in \f$X^*\f$, where \f$x\in X\f$.
        [[nodiscard]] Vector d1( const Vector& x ) const
        {
            assert( impl_ );
            return impl_->d1( x );
        }

        /// Compute second derivative as function space element in \f$X^*\f$, where \f$x,dx\in X\f$.
        [[nodiscard]] Vector d2( const Vector& x, const Vector& dx ) const
        {
            assert( impl_ );
            return impl_->d2( x, dx );
        }

        /// Access hessian as linear operator \f$ X \rightarrow X^\f$.
        [[nodiscard]] LinearOperator hessian( const Vector& x ) const
        {
            assert( impl_ );
            return impl_->hessian( x );
        }

        /// Access domain space \f$X\f$.
        [[nodiscard]] const VectorSpace& domain() const
        {
            assert( impl_ );
            return impl_->domain();
        }

        template < class T,
                   typename std::enable_if< !std::is_same< typename std::decay< T >::type, C2Functional >::value &&
                                            !std::is_base_of< Interface, typename std::decay< T >::type >::value >::type* = nullptr >
        C2Functional& operator=( T&& value )
        {
            return *this = C2Functional( std::forward< T >( value ) );
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
        clang::type_erasure::polymorphic::SBOCOWStorage< Interface, Wrapper, 16 > impl_;
    };
    /// @brief For a functional \f$ f: X\to \mathbb{R} \f$, compute \f$f'\f$ at \f$x\in X\f$ as dual
    /// element \f$ f'(x) \in X^\f$.
    ///     *
    /// Equivalent to calling f.d1(x).
    ///     *
    ///@param f twice differentiable functional
    ///@param x point of linearization
    ///@return \f$f'(x)\f$, i.e. f.d1(x).
    ///@see C2Functional
    inline Vector d1( const C2Functional& f, const Vector& x )
    {
        return f.d1( x );
    }

    /// @brief For a functional \f$ f: X\to \mathbb{R} \f$, compute \f$f''\f$ at \f$x\in X\f$ as
    /// linear operator \f$ f''(x): X \to X^\f$.
    ///     *
    /// Equivalent to calling f.hessian(x).
    ///     *
    ///@param f twice differentiable functional
    ///@param x point of linearization
    ///@return \f$f''(x)\f$, i.e. f.hessian(x).
    ///@see C2Functional
    inline LinearOperator d2( const C2Functional& f, const Vector& x )
    {
        return f.hessian( x );
    }

} // namespace Spacy
