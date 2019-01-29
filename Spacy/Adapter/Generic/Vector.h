#pragma once

#include <Spacy/Spaces/ScalarSpace/Real.h>
#include <Spacy/Util/Base/AddArithmeticOperators.h>
#include <Spacy/Util/Base/VectorBase.h>
#include <Spacy/Util/Mixins/Get.h>
#include <Spacy/Util/Voider.h>

namespace Spacy
{
    namespace Generic
    {
        /// @cond
        template < class >
        class Vector;

        namespace Detail
        {
            template < class T >
            using TryMemFn_dot = decltype( std::declval< T >().dot( std::declval< T >() ) );

            template < class T >
            using TryMemFn_inner = decltype( std::declval< T >().inner( std::declval< T >() ) );

            template < class T, class = void >
            struct HasMemFn_dot : std::false_type
            {
            };

            template < class T >
            struct HasMemFn_dot< T, voider< TryMemFn_dot< T > > > : std::true_type
            {
            };

            template < class T, class = void >
            struct HasMemFn_inner : std::false_type
            {
            };

            template < class T >
            struct HasMemFn_inner< T, voider< TryMemFn_inner< T > > > : std::true_type
            {
            };

            template < class VectorImpl, bool = HasMemFn_dot< VectorImpl >::value,
                       bool = HasMemFn_inner< VectorImpl >::value >
            struct DualPairingImpl
            {
                template < template < class > class Vector >
                static auto apply( const Vector< VectorImpl >& x, const Vector< VectorImpl >& y )
                {
                    return x.get() * y.get();
                }
            };

            template < class VectorImpl, bool hasMemFn_inner >
            struct DualPairingImpl< VectorImpl, true, hasMemFn_inner >
            {
                template < template < class > class Vector >
                static Spacy::Real apply( const Vector< VectorImpl >& x,
                                          const Vector< VectorImpl >& y )
                {
                    return x.get().dot( y.get() );
                }
            };

            template < class VectorImpl >
            struct DualPairingImpl< VectorImpl, false, true >
            {
                template < template < class > class Vector >
                static Spacy::Real apply( const Vector< VectorImpl >& x,
                                          const Vector< VectorImpl >& y )
                {
                    return x.get().inner( y.get() );
                }
            };

            template < class T >
            using TryIncrement = decltype( ++std::declval< T >() );

            template < class T >
            using TryDereference = decltype( *std::declval< T >() );

            template < class T, class = void >
            struct HasIncrement : std::false_type
            {
            };

            template < class T >
            struct HasIncrement< T, voider< TryIncrement< T > > > : std::true_type
            {
            };

            template < class T, class = void >
            struct HasDereference : std::false_type
            {
            };

            template < class T >
            struct HasDereference< T, voider< TryDereference< T > > > : std::true_type
            {
            };

            template < class T >
            constexpr bool isForwardIterator()
            {
                return HasIncrement< T >::value && HasDereference< T >::value;
            }

            template < class T >
            using TryBegin = decltype( std::declval< T >().begin() );

            template < class T >
            using TryEnd = decltype( std::declval< T >().end() );

            template < class T, class = void >
            struct HasBegin : std::false_type
            {
            };

            template < class T >
            struct HasBegin< T, voider< TryBegin< T > > >
                : std::integral_constant< bool, isForwardIterator< TryBegin< T > >() >
            {
            };

            template < class T, class = void >
            struct HasEnd : std::false_type
            {
            };

            template < class T >
            struct HasEnd< T, voider< TryEnd< T > > >
                : std::integral_constant< bool, isForwardIterator< TryEnd< T > >() >
            {
            };

            template < class T >
            constexpr bool hasForwardIterator()
            {
                return HasBegin< T >::value && HasEnd< T >::value;
            }

            template < class Derived, class VectorImpl, bool = hasForwardIterator< VectorImpl >() >
            struct IteratorBase
            {
                auto begin()
                {
                    return static_cast< Derived* >( this )->get().begin();
                }

                auto end()
                {
                    return static_cast< Derived* >( this )->get().end();
                }

                auto begin() const
                {
                    return static_cast< const Derived* >( this )->get().begin();
                }

                auto end() const
                {
                    return static_cast< const Derived* >( this )->get().end();
                }
            };

            template < class Derived, class VectorImpl >
            struct IteratorBase< Derived, VectorImpl, false >
            {
                ContiguousIterator< double > begin()
                {
                    return ContiguousIterator< double >(
                        static_cast< Derived* >( this )->get().data() );
                }

                ContiguousIterator< double > end()
                {
                    auto& v = static_cast< Derived* >( this )->get();
                    return ContiguousIterator< double >( v.data() + v.size() );
                }

                ContiguousIterator< const double > begin() const
                {
                    return ContiguousIterator< const double >(
                        static_cast< const Derived* >( this )->get().data() );
                }

                ContiguousIterator< const double > end() const
                {
                    auto& v = static_cast< const Derived* >( this )->get();
                    return ContiguousIterator< const double >( v.data() + v.size() );
                }
            };
        }
        /// @endcond

        /**
         * @ingroup GenericGroup
         * @brief Generic vector implementation for %Rn.
         */
        template < class VectorImpl >
        class Vector : public Mixin::Get< VectorImpl >,
                       public VectorBase,
                       public AddArithmeticOperators< Vector< VectorImpl > >,
                       public Detail::IteratorBase< Vector< VectorImpl >, VectorImpl >
        {
        public:
            /**
             * @brief Construct zero vector \f$x\f$ from underlying vector space.
             * @param v inital value
             * @param space underlying vector space
             */
            Vector( const VectorImpl& v, const VectorSpace& space )
                : Mixin::Get< VectorImpl >( v ), VectorBase( space )
            {
            }

            /**
             * @brief Assign from coefficient vector.
             * @param v coefficient vector
             */
            Vector& operator=( const VectorImpl& v )
            {
                this->get() = v;
                return *this;
            }

            /**
             * @brief Apply as dual element.
             * @param y primal vector
             * @return \f$x(y)\f$
             */
            Spacy::Real operator()( const Vector& y ) const
            {
                return Detail::DualPairingImpl< VectorImpl >::apply( *this, y );
            }
        };
    }
}
