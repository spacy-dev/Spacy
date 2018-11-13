// This file was automatically generated using clang-type-erase.
// Please do not modify.

#pragma once

#include <Spacy/Spaces/ScalarSpace/Real.h>
#include <Spacy/Util/SmartPointerStorage.h>
#include <Spacy/vector.hh>
#include <memory>
#include <type_traits>

namespace Spacy
{
    namespace CG
    {
        /// Regularization for conjugate gradient methods.
        class Regularization
        {
            struct Interface
            {
                virtual ~Interface() = default;
                virtual std::unique_ptr< Interface > clone() const = 0;
                virtual void init() = 0;
                virtual void apply( Real& qAq, Real qPq ) const = 0;
                virtual void update( Real qAq, Real qPq ) = 0;
                virtual void adjustResidual( Real alpha, const Vector& Pq, Vector& r ) const = 0;
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

                void init() override
                {
                    impl.init();
                }

                void apply( Real& qAq, Real qPq ) const override
                {
                    impl.apply( qAq, std::move( qPq ) );
                }

                void update( Real qAq, Real qPq ) override
                {
                    impl.update( std::move( qAq ), std::move( qPq ) );
                }

                void adjustResidual( Real alpha, const Vector& Pq, Vector& r ) const override
                {
                    impl.adjustResidual( std::move( alpha ), Pq, r );
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
            Regularization() noexcept = default;

            template <
                class T,
                typename std::enable_if<
                    !std::is_same< typename std::decay< T >::type, Regularization >::value &&
                    !std::is_base_of< Interface, typename std::decay< T >::type >::value >::type* =
                    nullptr >
            Regularization( T&& value ) : impl_( std::forward< T >( value ) )
            {
            }

            /// Initialize regularization.
            void init()
            {
                impl_->init();
            }

            /// @brief Apply regularization to qAq, where q is the conjugate search direction.
            ///@param qAq \f$ qAq \f$, where \f$q\f$ is the conjugate search direction and \f$A\f$
            /// the system operator
            ///@param qPq \f$ qPq \f$, where \f$q\f$ is the conjugate search direction and \f$P\f$
            /// the preconditioner
            void apply( Real& qAq, Real qPq ) const
            {
                impl_->apply( qAq, std::move( qPq ) );
            }

            /// @brief Update regularization (parameter).
            ///@param qAq \f$ qAq \f$, where \f$q\f$ is the conjugate search direction and \f$A\f$
            /// the system operator
            ///@param qPq \f$ qPq \f$, where \f$q\f$ is the conjugate search direction and \f$P\f$
            /// the preconditioner
            void update( Real qAq, Real qPq )
            {
                impl_->update( std::move( qAq ), std::move( qPq ) );
            }

            /// @brief Adjust residual for consistency with the regularized left hand side.
            ///@param alpha step length parameter of the conjugate gradient method
            ///@param Pq \f$ Pq \f$, where \f$q\f$ is the conjugate search direction and \f$P\f$ the
            /// preconditioner
            ///@param r residual
            void adjustResidual( Real alpha, const Vector& Pq, Vector& r ) const
            {
                impl_->adjustResidual( std::move( alpha ), Pq, r );
            }

            template <
                class T,
                typename std::enable_if<
                    !std::is_same< typename std::decay< T >::type, Regularization >::value &&
                    !std::is_base_of< Interface, typename std::decay< T >::type >::value >::type* =
                    nullptr >
            Regularization& operator=( T&& value )
            {
                return *this = Regularization( std::forward< T >( value ) );
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
    }
}
