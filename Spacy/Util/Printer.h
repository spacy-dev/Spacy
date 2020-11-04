// This file was automatically generated using clang-type-erase.
// Please do not modify.

#pragma once

#include <functional>
#include <ostream>

#include <Spacy/Util/SmartPointerStorage.h>

#include <memory>
#include <type_traits>

namespace Spacy
{
    namespace Log
    {
        /// Interface for logger's printers.
        class Printer
        {
            struct Interface
            {
                virtual ~Interface() = default;
                [[nodiscard]] virtual std::unique_ptr< Interface > clone() const = 0;
                virtual void call_const_char_ptr_const_char_ptr_const_std_function_void_std_ostream_ref_ref(
                    const char* tag, const char* name, const std::function< void( std::ostream& ) >& printable ) = 0;
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

                void call_const_char_ptr_const_char_ptr_const_std_function_void_std_ostream_ref_ref(
                    const char* tag, const char* name, const std::function< void( std::ostream& ) >& printable ) override
                {
                    impl.operator()( tag, name, printable );
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
            Printer() noexcept = default;

            template < class T,
                       typename std::enable_if< !std::is_same< typename std::decay< T >::type, Printer >::value &&
                                                !std::is_base_of< Interface, typename std::decay< T >::type >::value >::type* = nullptr >
            Printer( T&& value ) : impl_( std::forward< T >( value ) )
            {
            }

            /// @brief Log printable according to its level and name.
            ///@param tag specifies the algorithm/function this log originates from
            ///@param name name of the logged value
            ///@param printable wraps writing of the logged value to some stream
            void operator()( const char* tag, const char* name, const std::function< void( std::ostream& ) >& printable )
            {
                assert( impl_ );
                impl_->call_const_char_ptr_const_char_ptr_const_std_function_void_std_ostream_ref_ref( tag, name, printable );
            }

            template < class T,
                       typename std::enable_if< !std::is_same< typename std::decay< T >::type, Printer >::value &&
                                                !std::is_base_of< Interface, typename std::decay< T >::type >::value >::type* = nullptr >
            Printer& operator=( T&& value )
            {
                return *this = Printer( std::forward< T >( value ) );
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
    } // namespace Log
} // namespace Spacy
