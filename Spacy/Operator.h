// This file was automatically generated using clang-type-erase.
// Please do not modify.

#pragma once

#include <Spacy/Util/SmartPointerStorage.h>
#include <Spacy/Vector.h>
#include <Spacy/VectorSpace.h>
#include <memory>
#include <type_traits>
#include <functional>

namespace Spacy
{
    /// An operator that does not know about domain and range spaces.
    /// Good enough for some algorithms (i.e. for Krylov-methods).
    using CallableOperator = std::function< Vector( const Vector& ) >;
    /// Type-erased operator \f$A:\ X \to Y \f$.
    class Operator 
    {
        struct Interface
        {
            virtual ~Interface() = default;
            virtual std::shared_ptr< Interface > clone() const = 0;
            virtual Vector call_const_Vector_ref( const Vector& x ) const = 0;
            virtual const VectorSpace& domain() const = 0;
            virtual const VectorSpace& range() const = 0;
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

            Vector call_const_Vector_ref( const Vector& x ) const override
            {
                return impl.operator()( x );
            }

            const VectorSpace& domain() const override
            {
                return impl.domain();
            }

            const VectorSpace& range() const override
            {
                return impl.range();
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
        Operator() noexcept = default;

        template < class T, typename std::enable_if<
                                !std::is_same< typename std::decay< T >::type, Operator >::value &&
                                !std::is_base_of< Interface, typename std::decay< T >::type >::
                                    value >::type* = nullptr >
        Operator( T&& value ) : impl_( std::forward< T >( value ) )
        {
        }

        /// Apply operator.
        Vector operator()( const Vector& x ) const
        {
            assert( impl_ );
            return impl_->call_const_Vector_ref( x );
        }

        /// Access domain space \f$X\f$.
        const VectorSpace& domain() const
        {
            assert( impl_ );
            return impl_->domain();
        }

        /// Access range space \f$Y\f$.
        const VectorSpace& range() const
        {
            assert( impl_ );
            return impl_->range();
        }

        template < class T, typename std::enable_if<
                                !std::is_same< typename std::decay< T >::type, Operator >::value &&
                                !std::is_base_of< Interface, typename std::decay< T >::type >::
                                    value >::type* = nullptr >
        Operator& operator=( T&& value )
        {
            return *this = Operator( std::forward< T >( value ) );
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
        clang::type_erasure::polymorphic::SBOCOWStorage< Interface, Wrapper, 16 > impl_;
    };
    
    /// Added a new type of operator: quite often one also needs the transpose of an operator. Storing the transpose in memory would be redundant
    /// This is a straightforward way to provide this functionality. However, code-doubling occurs
    /// Alternative: if a transpose is needed, create a new, lightweight "transpose" operator whose application provied the transpose of the original one
    /// TODO: Evaluate alternatives
       class OperatorWithTranspose 
    {
        struct Interface
        {
            virtual ~Interface() = default;
            virtual std::shared_ptr< Interface > clone() const = 0;
            virtual Vector call_const_Vector_ref( const Vector& x ) const = 0;
            /// New functionality: apply transpose of the operator
            virtual Vector transposed(const Vector& x) const = 0;
            virtual const VectorSpace& domain() const = 0;
            virtual const VectorSpace& range() const = 0;
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

            Vector call_const_Vector_ref( const Vector& x ) const override
            {
                return impl.operator()( x );
            }

            Vector transposed( const Vector& x ) const override
            {
                return impl.transposed( x );
            }

            
            const VectorSpace& domain() const override
            {
                return impl.domain();
            }

            const VectorSpace& range() const override
            {
                return impl.range();
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
        OperatorWithTranspose() noexcept = default;

        template < class T, typename std::enable_if<
                                !std::is_same< typename std::decay< T >::type, OperatorWithTranspose >::value &&
                                !std::is_base_of< Interface, typename std::decay< T >::type >::
                                    value >::type* = nullptr >
        OperatorWithTranspose( T&& value ) : impl_( std::forward< T >( value ) )
        {
        }

        /// Apply operator.
        Vector operator()( const Vector& x ) const
        {
            assert( impl_ );
            return impl_->call_const_Vector_ref( x );
        }

        Vector transposed( const Vector& x ) const
        {
            assert( impl_ );
            return impl_->transposed( x );
        }
        
        /// Access domain space \f$X\f$.
        const VectorSpace& domain() const
        {
            assert( impl_ );
            return impl_->domain();
        }

        /// Access range space \f$Y\f$.
        const VectorSpace& range() const
        {
            assert( impl_ );
            return impl_->range();
        }

        template < class T, typename std::enable_if<
                                !std::is_same< typename std::decay< T >::type, OperatorWithTranspose >::value &&
                                !std::is_base_of< Interface, typename std::decay< T >::type >::
                                    value >::type* = nullptr >
        OperatorWithTranspose& operator=( T&& value )
        {
            return *this = OperatorWithTranspose( std::forward< T >( value ) );
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
        clang::type_erasure::polymorphic::SBOCOWStorage< Interface, Wrapper, 16 > impl_;
    };

} // namespace Spacy
