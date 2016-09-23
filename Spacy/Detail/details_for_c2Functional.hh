// This file was automatically generated by friendly type erasure.
// Please do not modify.

#pragma once

#include <memory>
#include <functional>
#include <Spacy/Spaces/RealSpace/real.hh>
#include <Spacy/Util/table_util.hh>
#include <Spacy/linearOperator.hh>
#include <Spacy/vector.hh>
#include <Spacy/vectorSpace.hh>

namespace Spacy
{
    namespace C2FunctionalDetail
    {
        template < class Interface, class Buffer >
        struct Functions
        {
            using clone_function = void ( * )( void*, std::shared_ptr< void >& );
            using clone_into_function = void ( * )( void*, Buffer&, std::shared_ptr< void >& );
            using call_const_Vector_ref_function = Real ( * )( const Interface&, void*, const Vector& x );
            using d1_function = Vector ( * )( const Interface&, void*, const Vector& x );
            using d2_function = Vector ( * )( const Interface&, void*, const Vector& x, const Vector& dx );
            using hessian_function = LinearOperator ( * )( const Interface&, void*, const Vector& x );
            using domain_function = const VectorSpace& (*)( const Interface&, void* );

            clone_function clone;
            clone_into_function clone_into;
            call_const_Vector_ref_function call_const_Vector_ref;
            d1_function d1;
            d2_function d2;
            hessian_function hessian;
            domain_function domain;
        };

        template < class Interface, class Impl >
        struct execution_wrapper
        {
            static Real call_const_Vector_ref( const Interface& interface, void* impl, const Vector& x )
            {
                return static_cast< const Impl* >( impl )->operator()( x );
            }

            static Vector d1( const Interface& interface, void* impl, const Vector& x )
            {
                return static_cast< const Impl* >( impl )->d1( x );
            }

            static Vector d2( const Interface& interface, void* impl, const Vector& x, const Vector& dx )
            {
                return static_cast< const Impl* >( impl )->d2( x, dx );
            }

            static LinearOperator hessian( const Interface& interface, void* impl, const Vector& x )
            {
                return static_cast< const Impl* >( impl )->hessian( x );
            }

            static const VectorSpace& domain( const Interface& interface, void* impl )
            {
                return static_cast< const Impl* >( impl )->domain();
            }
        };

        template < class Interface, class Impl >
        struct execution_wrapper< Interface, std::reference_wrapper< Impl > >
        {
            static Real call_const_Vector_ref( const Interface& interface, void* impl, const Vector& x )
            {
                return static_cast< std::reference_wrapper< Impl >* >( impl )->get().operator()( x );
            }

            static Vector d1( const Interface& interface, void* impl, const Vector& x )
            {
                return static_cast< std::reference_wrapper< Impl >* >( impl )->get().d1( x );
            }

            static Vector d2( const Interface& interface, void* impl, const Vector& x, const Vector& dx )
            {
                return static_cast< std::reference_wrapper< Impl >* >( impl )->get().d2( x, dx );
            }

            static LinearOperator hessian( const Interface& interface, void* impl, const Vector& x )
            {
                return static_cast< std::reference_wrapper< Impl >* >( impl )->get().hessian( x );
            }

            static const VectorSpace& domain( const Interface& interface, void* impl )
            {
                return static_cast< std::reference_wrapper< Impl >* >( impl )->get().domain();
            }
        };

        template < class T >
        using TryMemFn_call_const_Vector_ref = decltype( std::declval< T >().operator()( std::declval< Vector >() ) );
        template < class T, class = void >
        struct HasMemFn_call_const_Vector_ref : std::false_type
        {
        };
        template < class T >
        struct HasMemFn_call_const_Vector_ref<
            T, type_erasure_table_detail::voider< TryMemFn_call_const_Vector_ref< T > > > : std::true_type
        {
        };
        template < class T >
        using TryMemFn_d1 = decltype( std::declval< T >().d1( std::declval< Vector >() ) );
        template < class T, class = void >
        struct HasMemFn_d1 : std::false_type
        {
        };
        template < class T >
        struct HasMemFn_d1< T, type_erasure_table_detail::voider< TryMemFn_d1< T > > > : std::true_type
        {
        };
        template < class T >
        using TryMemFn_d2 = decltype( std::declval< T >().d2( std::declval< Vector >(), std::declval< Vector >() ) );
        template < class T, class = void >
        struct HasMemFn_d2 : std::false_type
        {
        };
        template < class T >
        struct HasMemFn_d2< T, type_erasure_table_detail::voider< TryMemFn_d2< T > > > : std::true_type
        {
        };
        template < class T >
        using TryMemFn_hessian = decltype( std::declval< T >().hessian( std::declval< Vector >() ) );
        template < class T, class = void >
        struct HasMemFn_hessian : std::false_type
        {
        };
        template < class T >
        struct HasMemFn_hessian< T, type_erasure_table_detail::voider< TryMemFn_hessian< T > > > : std::true_type
        {
        };
        template < class T >
        using TryMemFn_domain = decltype( std::declval< T >().domain() );
        template < class T, class = void >
        struct HasMemFn_domain : std::false_type
        {
        };
        template < class T >
        struct HasMemFn_domain< T, type_erasure_table_detail::voider< TryMemFn_domain< T > > > : std::true_type
        {
        };

        template < class T >
        using C2FunctionalConceptImpl = std::integral_constant<
            bool, HasMemFn_call_const_Vector_ref< type_erasure_table_detail::remove_reference_wrapper_t< T > >::value &&
                      HasMemFn_d1< type_erasure_table_detail::remove_reference_wrapper_t< T > >::value &&
                      HasMemFn_d2< type_erasure_table_detail::remove_reference_wrapper_t< T > >::value &&
                      HasMemFn_hessian< type_erasure_table_detail::remove_reference_wrapper_t< T > >::value &&
                      HasMemFn_domain< type_erasure_table_detail::remove_reference_wrapper_t< T > >::value >;

        template < class Impl, class T, bool = std::is_same< Impl, T >::value >
        struct C2FunctionalConcept : std::false_type
        {
        };
        template < class Impl, class T >
        struct C2FunctionalConcept< Impl, T, false > : C2FunctionalConceptImpl< T >
        {
        };
    }
}
