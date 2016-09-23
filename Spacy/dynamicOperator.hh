// This file was automatically generated by friendly type erasure.
// Please do not modify.

#pragma once

#include <array>
#include <functional>
#include "Spacy/Util/table_util.hh"
#include <Spacy/linearOperator.hh>
#include <Spacy/vector.hh>
#include <Spacy/vectorSpace.hh>
#include "Detail/details_for_dynamicOperator.hh"

namespace Spacy
{
    /**
    * A time-dependent operator that does not know about domain and range spaces.
    */
    using DynamicCallableOperator = std::function< Vector( double, const Vector& ) >;

    /// Type-erased time-dependent operator \f$A:\ [0,T] \times X \to Y \f$.
    class DynamicOperator
    {
    public:
        DynamicOperator() noexcept;

        template < typename T, typename std::enable_if< DynamicOperatorDetail::DynamicOperatorConcept<
                                   DynamicOperator, typename std::decay< T >::type >::value >::type* = nullptr >
        DynamicOperator( T&& value )
            : functions_(
                  {&type_erasure_table_detail::clone_into_shared_ptr< typename std::decay< T >::type >,
                   &type_erasure_table_detail::clone_into_buffer< typename std::decay< T >::type, Buffer >,
                   &DynamicOperatorDetail::execution_wrapper< DynamicOperator,
                                                              typename std::decay< T >::type >::call_const_Vector_ref,
                   &DynamicOperatorDetail::execution_wrapper< DynamicOperator, typename std::decay< T >::type >::M,
                   &DynamicOperatorDetail::execution_wrapper< DynamicOperator, typename std::decay< T >::type >::domain,
                   &DynamicOperatorDetail::execution_wrapper< DynamicOperator,
                                                              typename std::decay< T >::type >::range} ),
              type_id_( typeid( typename std::decay< T >::type ).hash_code() ), impl_( nullptr )
        {
            if ( sizeof( typename std::decay< T >::type ) <= sizeof( Buffer ) )
            {
                new ( &buffer_ ) typename std::decay< T >::type( std::forward< T >( value ) );
                impl_ = std::shared_ptr< typename std::decay< T >::type >(
                    std::shared_ptr< typename std::decay< T >::type >(),
                    static_cast< typename std::decay< T >::type* >( static_cast< void* >( &buffer_ ) ) );
            }
            else
                impl_ = std::make_shared< typename std::decay< T >::type >( std::forward< T >( value ) );
        }

        DynamicOperator( const DynamicOperator& other );

        DynamicOperator( DynamicOperator&& other ) noexcept;

        template < typename T, typename std::enable_if< DynamicOperatorDetail::DynamicOperatorConcept<
                                   DynamicOperator, typename std::decay< T >::type >::value >::type* = nullptr >
        DynamicOperator& operator=( T&& value )
        {
            return *this = DynamicOperator( std::forward< T >( value ) );
        }

        DynamicOperator& operator=( const DynamicOperator& other );

        DynamicOperator& operator=( DynamicOperator&& other ) noexcept;

        /**
         * @brief Checks if the type-erased interface holds an implementation.
         * @return true if an implementation is stored, else false
         */
        explicit operator bool() const noexcept;

        /// Apply operator.
        Vector operator()( const Vector& x ) const;

        LinearOperator M() const;

        /// Access domain space \f$X\f$.
        const VectorSpace& domain() const;

        /// Access range space \f$Y\f$.
        const VectorSpace& range() const;

        /**
        * @brief Conversion of the stored implementation to @code  T* @endcode.
        * @return pointer to the stored object if conversion was successful, else nullptr
        */
        template < class T >
        T* target() noexcept
        {
            return type_erasure_table_detail::dynamic_cast_impl< T >( type_id_, read() );
        }

        /**
        * @brief Conversion of the stored implementation to @code const T* @endcode.
        * @return pointer to the stored object if conversion was successful, else nullptr
        */
        template < class T >
        const T* target() const noexcept
        {
            return type_erasure_table_detail::dynamic_cast_impl< T >( type_id_, read() );
        }

    private:
        using Buffer = std::array< char, 64 >;

        void* read() const noexcept;

        void* write();

        DynamicOperatorDetail::Functions< DynamicOperator, Buffer > functions_;
        std::size_t type_id_;
        std::shared_ptr< void > impl_ = nullptr;
        Buffer buffer_;
    };

    /// Type-erased time-dependent linear operator \f$A:\ [0,T] \times X \to Y \f$.
    class DynamicLinearOperator
    {
    public:
        DynamicLinearOperator() noexcept;

        template < typename T, typename std::enable_if< DynamicLinearOperatorDetail::DynamicLinearOperatorConcept<
                                   DynamicLinearOperator, typename std::decay< T >::type >::value >::type* = nullptr >
        DynamicLinearOperator( T&& value )
            : functions_(
                  {&type_erasure_table_detail::clone_into_shared_ptr< typename std::decay< T >::type >,
                   &type_erasure_table_detail::clone_into_buffer< typename std::decay< T >::type, Buffer >,
                   &DynamicLinearOperatorDetail::execution_wrapper<
                       DynamicLinearOperator, typename std::decay< T >::type >::call_double_const_Vector_ref,
                   &DynamicLinearOperatorDetail::execution_wrapper<
                       DynamicLinearOperator, typename std::decay< T >::type >::add_const_DynamicLinearOperator_ref,
                   &DynamicLinearOperatorDetail::execution_wrapper<
                       DynamicLinearOperator,
                       typename std::decay< T >::type >::subtract_const_DynamicLinearOperator_ref,
                   &DynamicLinearOperatorDetail::execution_wrapper< DynamicLinearOperator,
                                                                    typename std::decay< T >::type >::multiply_double,
                   &DynamicLinearOperatorDetail::execution_wrapper< DynamicLinearOperator,
                                                                    typename std::decay< T >::type >::negate,
                   &DynamicLinearOperatorDetail::execution_wrapper<
                       DynamicLinearOperator, typename std::decay< T >::type >::compare_const_DynamicLinearOperator_ref,
                   &DynamicLinearOperatorDetail::execution_wrapper< DynamicLinearOperator,
                                                                    typename std::decay< T >::type >::solver,
                   &DynamicLinearOperatorDetail::execution_wrapper< DynamicLinearOperator,
                                                                    typename std::decay< T >::type >::domain,
                   &DynamicLinearOperatorDetail::execution_wrapper< DynamicLinearOperator,
                                                                    typename std::decay< T >::type >::range,
                   &DynamicLinearOperatorDetail::execution_wrapper< DynamicLinearOperator,
                                                                    typename std::decay< T >::type >::space} ),
              type_id_( typeid( typename std::decay< T >::type ).hash_code() ), impl_( nullptr )
        {
            if ( sizeof( typename std::decay< T >::type ) <= sizeof( Buffer ) )
            {
                new ( &buffer_ ) typename std::decay< T >::type( std::forward< T >( value ) );
                impl_ = std::shared_ptr< typename std::decay< T >::type >(
                    std::shared_ptr< typename std::decay< T >::type >(),
                    static_cast< typename std::decay< T >::type* >( static_cast< void* >( &buffer_ ) ) );
            }
            else
                impl_ = std::make_shared< typename std::decay< T >::type >( std::forward< T >( value ) );
        }

        DynamicLinearOperator( const DynamicLinearOperator& other );

        DynamicLinearOperator( DynamicLinearOperator&& other ) noexcept;

        template < typename T, typename std::enable_if< DynamicLinearOperatorDetail::DynamicLinearOperatorConcept<
                                   DynamicLinearOperator, typename std::decay< T >::type >::value >::type* = nullptr >
        DynamicLinearOperator& operator=( T&& value )
        {
            return *this = DynamicLinearOperator( std::forward< T >( value ) );
        }

        DynamicLinearOperator& operator=( const DynamicLinearOperator& other );

        DynamicLinearOperator& operator=( DynamicLinearOperator&& other ) noexcept;

        /**
         * @brief Checks if the type-erased interface holds an implementation.
         * @return true if an implementation is stored, else false
         */
        explicit operator bool() const noexcept;

        /// Apply operator.
        Vector operator()( double t, const Vector& x ) const;

        DynamicLinearOperator& operator+=( const DynamicLinearOperator& y );

        DynamicLinearOperator& operator-=( const DynamicLinearOperator& y );

        DynamicLinearOperator& operator*=( double a );

        DynamicLinearOperator operator-() const;

        bool operator==( const DynamicLinearOperator& y ) const;

        std::function< Vector( const Vector& ) > solver() const;

        /// Access domain space \f$X\f$.
        const VectorSpace& domain() const;

        /// Access range space \f$Y\f$.
        const VectorSpace& range() const;

        /// Access underlying space of linear operators.
        const VectorSpace& space() const;

        /**
        * @brief Conversion of the stored implementation to @code  T* @endcode.
        * @return pointer to the stored object if conversion was successful, else nullptr
        */
        template < class T >
        T* target() noexcept
        {
            return type_erasure_table_detail::dynamic_cast_impl< T >( type_id_, read() );
        }

        /**
        * @brief Conversion of the stored implementation to @code const T* @endcode.
        * @return pointer to the stored object if conversion was successful, else nullptr
        */
        template < class T >
        const T* target() const noexcept
        {
            return type_erasure_table_detail::dynamic_cast_impl< T >( type_id_, read() );
        }

    private:
        using Buffer = std::array< char, 64 >;

        void* read() const noexcept;

        void* write();

        DynamicLinearOperatorDetail::Functions< DynamicLinearOperator, Buffer > functions_;
        std::size_t type_id_;
        std::shared_ptr< void > impl_ = nullptr;
        Buffer buffer_;
    };

    /// Type-erased time-dependent differentiable operator \f$A:\ [0,T] \times X \to Y \f$.
    class DynamicC1Operator
    {
    public:
        DynamicC1Operator() noexcept;

        template < typename T, typename std::enable_if< DynamicC1OperatorDetail::DynamicC1OperatorConcept<
                                   DynamicC1Operator, typename std::decay< T >::type >::value >::type* = nullptr >
        DynamicC1Operator( T&& value )
            : functions_(
                  {&type_erasure_table_detail::clone_into_shared_ptr< typename std::decay< T >::type >,
                   &type_erasure_table_detail::clone_into_buffer< typename std::decay< T >::type, Buffer >,
                   &DynamicC1OperatorDetail::execution_wrapper<
                       DynamicC1Operator, typename std::decay< T >::type >::call_double_const_Vector_ref,
                   &DynamicC1OperatorDetail::execution_wrapper< DynamicC1Operator, typename std::decay< T >::type >::d1,
                   &DynamicC1OperatorDetail::execution_wrapper< DynamicC1Operator,
                                                                typename std::decay< T >::type >::linearization,
                   &DynamicC1OperatorDetail::execution_wrapper< DynamicC1Operator, typename std::decay< T >::type >::M,
                   &DynamicC1OperatorDetail::execution_wrapper< DynamicC1Operator,
                                                                typename std::decay< T >::type >::domain,
                   &DynamicC1OperatorDetail::execution_wrapper< DynamicC1Operator,
                                                                typename std::decay< T >::type >::range} ),
              type_id_( typeid( typename std::decay< T >::type ).hash_code() ), impl_( nullptr )
        {
            if ( sizeof( typename std::decay< T >::type ) <= sizeof( Buffer ) )
            {
                new ( &buffer_ ) typename std::decay< T >::type( std::forward< T >( value ) );
                impl_ = std::shared_ptr< typename std::decay< T >::type >(
                    std::shared_ptr< typename std::decay< T >::type >(),
                    static_cast< typename std::decay< T >::type* >( static_cast< void* >( &buffer_ ) ) );
            }
            else
                impl_ = std::make_shared< typename std::decay< T >::type >( std::forward< T >( value ) );
        }

        DynamicC1Operator( const DynamicC1Operator& other );

        DynamicC1Operator( DynamicC1Operator&& other ) noexcept;

        template < typename T, typename std::enable_if< DynamicC1OperatorDetail::DynamicC1OperatorConcept<
                                   DynamicC1Operator, typename std::decay< T >::type >::value >::type* = nullptr >
        DynamicC1Operator& operator=( T&& value )
        {
            return *this = DynamicC1Operator( std::forward< T >( value ) );
        }

        DynamicC1Operator& operator=( const DynamicC1Operator& other );

        DynamicC1Operator& operator=( DynamicC1Operator&& other ) noexcept;

        /**
         * @brief Checks if the type-erased interface holds an implementation.
         * @return true if an implementation is stored, else false
         */
        explicit operator bool() const noexcept;

        /// Apply operator.
        Vector operator()( double t, const Vector& x ) const;

        /// Compute directional derivative \f$A'(x)\delta x\f$.
        Vector d1( double t, const Vector& x, const Vector& dx ) const;

        /// Get linearization \f$A'(x):\ X\to Y \f$
        LinearOperator linearization( double t, const Vector& x ) const;

        LinearOperator M() const;

        /// Access domain space \f$X\f$.
        const VectorSpace& domain() const;

        /// Access range space \f$Y\f$.
        const VectorSpace& range() const;

        /**
        * @brief Conversion of the stored implementation to @code  T* @endcode.
        * @return pointer to the stored object if conversion was successful, else nullptr
        */
        template < class T >
        T* target() noexcept
        {
            return type_erasure_table_detail::dynamic_cast_impl< T >( type_id_, read() );
        }

        /**
        * @brief Conversion of the stored implementation to @code const T* @endcode.
        * @return pointer to the stored object if conversion was successful, else nullptr
        */
        template < class T >
        const T* target() const noexcept
        {
            return type_erasure_table_detail::dynamic_cast_impl< T >( type_id_, read() );
        }

    private:
        using Buffer = std::array< char, 64 >;

        void* read() const noexcept;

        void* write();

        DynamicC1OperatorDetail::Functions< DynamicC1Operator, Buffer > functions_;
        std::size_t type_id_;
        std::shared_ptr< void > impl_ = nullptr;
        Buffer buffer_;
    };
}
