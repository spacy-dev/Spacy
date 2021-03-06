#pragma once

#include "Copy.h"
#include "DirectSolver.h"
#include "OperatorSpace.h"
#include "Util/NumaMatrixAsOperator.h"

#include <Spacy/LinearSolver.h>
#include <Spacy/Util/Base/OperatorBase.h>
#include <Spacy/Util/Base/VectorBase.h>
#include <Spacy/Util/Exceptions.h>
#include <Spacy/Vector.h>

#include <dune/istl/operators.hh>

#include <utility>

/** @addtogroup KaskadeGroup
 * @{
 */
namespace Spacy::Kaskade
{
    /**
     * @brief Linear operator interface for operators in %Kaskade7.
     * @tparam OperatorImpl %Kaskade7 operator, i.e. %Kaskade::AssembledGalerkinOperator,
     * %Kaskade::MatrixRepresentedOperator, Dune::MatrixAdapter, ...
     * @tparam AnsatzVariableSetDescription %Kaskade::VariableSetDescription for ansatz
     * variables
     * @tparam TestVariableSetDescription %Kaskade::VariableSetDescription for test variables
     * @see ::Spacy::LinearOperator
     */
    template < class AnsatzVariableSetDescription, class TestVariableSetDescription, class Matrix = ::Kaskade::MatrixAsTriplet< double >,
               bool transposed = false, bool shared = false >
    class LinearOperator : public OperatorBase,
                           public VectorBase,
                           public AddArithmeticOperators<
                               LinearOperator< AnsatzVariableSetDescription, TestVariableSetDescription, Matrix, transposed, shared > >
    {
        using Spaces = typename AnsatzVariableSetDescription::Spaces;
        using Variables = typename AnsatzVariableSetDescription::Variables;
        using Domain = typename AnsatzVariableSetDescription::template CoefficientVectorRepresentation<>::type;
        using Range = typename TestVariableSetDescription::template CoefficientVectorRepresentation<>::type;
        using OperatorImpl = std::conditional_t< shared, std::shared_ptr< ::Kaskade::MatrixRepresentedOperator< Matrix, Domain, Range > >,
                                                 ::Kaskade::MatrixRepresentedOperator< Matrix, Domain, Range > >;
        using OperatorCreator = LinearOperatorCreator< AnsatzVariableSetDescription, TestVariableSetDescription >;

    public:
        /**
         * @brief Construct linear operator for %Kaskade 7.
         * @param A operator from %Kaskade 7
         * @param domain domain space
         * @param range range space
         */
        LinearOperator( OperatorImpl A, const VectorSpace& space )
            : OperatorBase( cast_ref< OperatorCreator >( space.creator() ).domain(),
                            cast_ref< OperatorCreator >( space.creator() ).range() ),
              VectorBase( space ), A_( std::move( A ) ),
              spaces_( extractSpaces< AnsatzVariableSetDescription >( cast_ref< OperatorCreator >( space.creator() ).domain() ) )
        {
        }

        /**
         * @brief Construct linear operator for %Kaskade 7.
         * @param A operator from %Kaskade 7
         * @param domain domain space
         * @param range range space
         * @param solverCreator function/functor implementing the creation of a linear solver
         */
        LinearOperator( OperatorImpl A, const VectorSpace& space, std::function< LinearSolver( const LinearOperator& ) > solverCreator )
            : OperatorBase( cast_ref< OperatorCreator >( space.creator() ).domain(),
                            cast_ref< OperatorCreator >( space.creator() ).range() ),
              VectorBase( space ), A_( std::move( A ) ),
              spaces_( extractSpaces< AnsatzVariableSetDescription >( cast_ref< OperatorCreator >( space.creator() ).domain() ) ),
              solverCreator_( std::move( solverCreator ) )
        {
        }

        /// Compute \f$A(x)\f$.
        [[nodiscard]] ::Spacy::Vector operator()( const ::Spacy::Vector& x ) const
        {
            Domain x_( AnsatzVariableSetDescription::template CoefficientVectorRepresentation<>::init( spaces_ ) );
            copy< AnsatzVariableSetDescription >( x, x_ );
            Range y_( TestVariableSetDescription::template CoefficientVectorRepresentation<>::init( spaces_ ) );

            if constexpr ( transposed )
            {
                if constexpr ( shared )
                {
                    A_->applyTransposed( x_, y_ );
                }
                else
                {
                    A_.applyTransposed( x_, y_ );
                }
            }
            else
            {
                if constexpr ( shared )
                {
                    A_->apply( x_, y_ );
                }
                else
                {
                    A_.apply( x_, y_ );
                }
            }

            auto y = zero( range() );
            copy< TestVariableSetDescription >( y_, y );

            return y;
        }

        ::Spacy::Real operator()( const LinearOperator& /*unused*/ ) const
        {
            return Real( 0 );
        }

        /// Access solver representing \f$A^{-1}\f$.
        [[nodiscard]] auto solver() const
        {
            return solverCreator_( *this );
        }

        //      double operator()(const ::Spacy::Vector& x) const
        //      {
        //        const auto& x_ = cast_ref<
        //        LinearOperator<AnsatzVariableSetDescription,TestVariableSetDescription> >(x);
        //        assert( impl().N() == x_.impl().N() && impl().M() == x_.impl().M() );
        //        auto iend = end(impl());
        //        auto iter =begin(impl()), x_iter(x_.impl());
        //        auto result = 0.;
        //        for( ; iter < iend ; ++iter , ++x_iter )
        //          result += (*iter) * (*x_iter);
        //        return result;
        //      }

        decltype( auto ) get()
        {
            if constexpr ( shared )
            {
                return A_->get_non_const();
            }
            else
            {
                return A_.get_non_const();
            }
        }

        [[nodiscard]] decltype( auto ) get() const
        {
            if constexpr ( shared )
            {
                return A_->template get< Matrix >();
            }
            else
            {
                return A_.template get< Matrix >();
            }
        }

        [[nodiscard]] auto& A()
        {
            if constexpr ( shared )
            {
                return *A_;
            }
            else
            {
                return A_;
            }
        }

        [[nodiscard]] const auto& A() const
        {
            if constexpr ( shared )
            {
                return *A_;
            }
            else
            {
                return A_;
            }
        }

        [[nodiscard]] Spacy::ContiguousIterator< double > begin()
        {
            return Spacy::ContiguousIterator< double >{};
        }

        [[nodiscard]] Spacy::ContiguousIterator< const double > begin() const
        {
            return Spacy::ContiguousIterator< const double >{};
        }

        [[nodiscard]] Spacy::ContiguousIterator< double > end()
        {
            return Spacy::ContiguousIterator< double >{};
        }

        [[nodiscard]] Spacy::ContiguousIterator< const double > end() const
        {
            return Spacy::ContiguousIterator< const double >{};
        }

    private:
        OperatorImpl A_;
        Spaces spaces_;
        std::function< LinearSolver( const LinearOperator& ) > solverCreator_ = []( const LinearOperator& M ) {
            if constexpr ( transposed )
            {
                throw Exception::NotImplemented( "Default solver creator", "transposed matrices" );
            }
            else
            {
                return makeDirectSolver< TestVariableSetDescription, AnsatzVariableSetDescription >( M.A(), M.range(), M.domain() );
            }
        };
    };

    /**
     * @brief Convenient generation of a linear operator for %Kaskade 7.
     * @param A operator from %Kaskade 7
     * @param domain domain space
     * @param range range space
     * @tparam OperatorImpl %Kaskade 7 operator, i.e. %Kaskade::AssembledGalerkinOperator or
     * %Kaskade::MatrixRepresentedOperator.
     * @tparam AnsatzVariableSetDescription %Kaskade::VariableSetDescription for ansatz
     * variables
     * @tparam TestVariableSetDescription %Kaskade::VariableSetDescription for test variables
     * @return LinearOperator<OperatorImpl, AnsatzVariableSetDescription,
     * TestVariableSetDescription>( A , domain , range )
     */
    template < class AnsatzVariableSetDescription, class TestVariableSetDescription, class OperatorImpl >
    auto makeLinearOperator( const OperatorImpl& A, const VectorSpace& domain, const VectorSpace& range )
    {
        return LinearOperator< AnsatzVariableSetDescription, TestVariableSetDescription, OperatorImpl >( A, domain, range );
    }
} // namespace Spacy::Kaskade
/** @} */
