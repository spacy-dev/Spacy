#pragma once

#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include <Spacy/LinearSolver.h>
#include <Spacy/Util/Base/OperatorBase.h>
#include <Spacy/Util/Base/VectorBase.h>
#include <Spacy/Util/Cast.h>
#include <Spacy/Util/Exceptions.h>
#include <Spacy/Util/Mixins/Get.h>
#include <Spacy/Vector.h>
#include <Spacy/VectorSpace.h>
#include <Spacy/ZeroVectorCreator.h>

#include "CG.h"
#include "UMFPACK.h"
#include "Vector.h"

#include <algorithm>

namespace Spacy
{
    namespace dealII
    {
        template < int dim, class VariableDims >
        class LinearOperator : public OperatorBase,
                               public VectorBase,
                               public Mixin::Get< typename Traits< VariableDims >::Matrix >
        {
        public:
            using NativeMatrix = typename Traits< VariableDims >::Matrix;

            LinearOperator( const NativeMatrix& A, const VectorSpace& operatorSpace,
                            const VectorSpace& domain, const VectorSpace& range )
                : OperatorBase( domain, range ), VectorBase( operatorSpace ),
                  Mixin::Get< NativeMatrix >( A.get_sparsity_pattern() )
            {
                this->get().copy_from( A );
            }

            LinearOperator( const LinearOperator& other )
                : OperatorBase( other.domain(), other.range() ), VectorBase( other.space() ),
                  Mixin::Get< NativeMatrix >( other.get().get_sparsity_pattern() )
            {
                checkSpaceCompatibility( domain(), other.domain() );
                checkSpaceCompatibility( range(), other.range() );

                this->get().copy_from( other.get() );
            }

            LinearOperator& operator=( const LinearOperator& other )
            {
                checkSpaceCompatibility( domain(), other.domain() );
                checkSpaceCompatibility( range(), other.range() );

                this->get().copy_from( other.get() );
                return *this;
            }

            /**
             * @brief Compute \f$A(x)\f$.
             * @param x operator argument
             */
            ::Spacy::Vector operator()( const ::Spacy::Vector& x ) const
            {
                const auto& x_ = cast_ref< Vector >( x );
                auto y = zero( range() );
                auto& y_ = cast_ref< Vector >( y );
                this->get().vmult( y_.get(), x_.get() );
                return y;
            }

            Real operator()( const LinearOperator& ) const
            {
                throw Exception::CallOfUndefinedFunction( __func__ );
            }

            LinearSolver solver() const
            {
                return UMFPackSolver< dim, VariableDims >( this->get(), domain(), range() );
                //                return CGSolver< dim, VariableDims >(this->get(), domain(),
                //                range());
            }

            /**
             * @brief In-place summation \f$ x+=y\f$.
             * @param y vector to add to this vector
             * @return \f$ x+=y\f$.
             */
            LinearOperator& operator+=( const LinearOperator& y )
            {
                checkSpaceCompatibility( static_cast< const LinearOperator* >( this )->space(),
                                         y.space() );
                this->get().add( 1, y.get() );
                return *this;
            }

            /**
             * @brief In-place subtraction \f$ x-=y\f$.
             * @param y vector to subtract from this vector
             * @return \f$ x-=y\f$.
             */
            LinearOperator& operator-=( const LinearOperator& y )
            {
                checkSpaceCompatibility( static_cast< const LinearOperator* >( this )->space(),
                                         y.space() );
                this->get().add( -1, y.get() );
                return *this;
            }

            /**
             * @brief In-place multiplication \f$ x*=a\f$.
             * @param a scaling factor
             * @return \f$ x*=a\f$.
             */
            LinearOperator& operator*=( double a )
            {
                this->get() *= a;
                return *this;
            }

            /**
             * @brief Negation \f$ -x\f$.
             * @return \f$ -x \f$.
             */
            LinearOperator operator-() const
            {
                LinearOperator y( *this );
                y.get() *= -1;
                return y;
            }

            /**
             * @brief Comparison operator \f$ x==y\f$.
             * @param y vector to compare with this vector
             * @return \f$ x==y\f$.
             */
            bool operator==( const LinearOperator& y ) const
            {
                checkSpaceCompatibility( space(), y.space() );
                auto max = std::max( ( *this )( *this ).get(), y( y ).get() );
                if ( max == 0 )
                    max = 1;
                auto dx = y;
                dx -= *this;
                return dx( dx ).get() < max * y.space().eps() * y.space().eps();
            }

            ContiguousIterator< double > begin()
            {
                return ContiguousIterator< double >( this->get().data() );
            }

            ContiguousIterator< double > end()
            {
                return ContiguousIterator< double >( this->get().data() + this->get().size() );
            }

            ContiguousIterator< const double > begin() const
            {
                return ContiguousIterator< const double >( this->get().data() );
            }

            ContiguousIterator< const double > end() const
            {
                return ContiguousIterator< const double >( this->get().data() +
                                                           this->get().size() );
            }
        };

        /// Dummy linear operator creator.
        /// @warning Throws on invocation.
        struct LinearOperatorCreator
        {
            Vector operator()( const VectorSpace* ) const
            {
                throw Exception::CallOfUndefinedFunction( __func__ );
            }
        };
    }
}
