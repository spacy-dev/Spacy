#pragma once

#include "Copy.h"
#include "VectorSpace.h"

#include <Spacy/Util/Base/OperatorBase.h>
#include <Spacy/Util/Cast.h>
#include <Spacy/Vector.h>

#include "linalg/direct.hh"

namespace Spacy::Kaskade
{
    /**
     * @ingroup KaskadeGroup
     * @brief Direct solver interface for %Kaskade 7.
     */
    template < class KaskadeOperator, class AnsatzVariableDescription, class TestVariableDescription >
    class DirectSolver : public OperatorBase
    {
        using Spaces = typename AnsatzVariableDescription::Spaces;
        using Domain = typename AnsatzVariableDescription::template CoefficientVectorRepresentation<>::type;
        using Range = typename TestVariableDescription::template CoefficientVectorRepresentation<>::type;

    public:
        DirectSolver() = delete;
        /**
         * @brief Constructor.
         * @param A %Kaskade operator (i.e., %AssembledGalerkinOperator or
         * %MatrixRepresentedOperator)
         * @param spaces boost fusion forward sequence of space pointers required to initialize
         * temporary storage
         * @param domain domain space of the solver
         * @param range range space of the solver
         * @param directSolver solver type (DirectType::MUMPS (default), DirectType::UMFPACK,
         * DirectType::UMFPACK3264 or DirectType::SUPERLU)
         * @param property matrix property (MatrixProperties::GENERAL (default) or
         * MatrixProperties::SYMMETRIC)
         */
        DirectSolver( KaskadeOperator A, const VectorSpace& domain, const VectorSpace& range,
                      ::Kaskade::DirectType directSolver = ::Kaskade::DirectType::UMFPACK3264,
                      ::Kaskade::MatrixProperties property = ::Kaskade::MatrixProperties::GENERAL )
            : OperatorBase( domain, range ), A_( std::move( A ) ), spaces_( extractSpaces< AnsatzVariableDescription >( domain ) ),
              directSolver_( directSolver ), property_( property )
        {
        }

        /// Compute \f$A^{-1}x\f$.
        ::Spacy::Vector operator()( const ::Spacy::Vector& x ) const
        {
            if ( solver_ == nullptr )
                solver_ = std::make_shared< ::Kaskade::InverseLinearOperator< ::Kaskade::DirectSolver< Domain, Range > > >(
                    ::Kaskade::directInverseOperator( A_, directSolver_, property_ ) );

            Range y_( TestVariableDescription::template CoefficientVectorRepresentation<>::init( spaces_ ) );
            Domain x_( AnsatzVariableDescription::template CoefficientVectorRepresentation<>::init( spaces_ ) );
            copy< AnsatzVariableDescription >( x, x_ );

            solver_->apply( x_, y_ );

            auto y = zero( range() );
            copy< TestVariableDescription >( y_, y );

            return y;
        }

    private:
        KaskadeOperator A_;
        Spaces spaces_;
        ::Kaskade::DirectType directSolver_ = ::Kaskade::DirectType::UMFPACK3264;
        ::Kaskade::MatrixProperties property_ = ::Kaskade::MatrixProperties::GENERAL;
        mutable std::shared_ptr< ::Kaskade::InverseLinearOperator< ::Kaskade::DirectSolver< Domain, Range > > > solver_ = nullptr;
    };

    /**
     * @brief Convenient generation of direct solver for %Kaskade 7.
     * @param A %Kaskade operator (i.e., AssembledGalerkinOperator or MatrixRepresentedOperator)
     * @param spaces boost fusion forward sequence of space pointers required to initialize
     * temporary storage
     * @param domain domain space of the solver
     * @param range range space of the solver
     * @param directSolver solver type (DirectType::MUMPS (default), DirectType::UMFPACK,
     * DirectType::UMFPACK3264 or DirectType::SUPERLU)
     * @param property matrix property (MatrixProperties::GENERAL (default) or
     * MatrixProperties::SYMMETRIC)
     * @return
     * DirectSolver<KaskadeOperator,AnsatzVariableSetDescription,TestVariableSetDescription>( A
     * , spaces , domain , range , directSolver , property )
     */
    template < class AnsatzVariableSetDescription, class TestVariableSetDescription, class KaskadeOperator >
    auto makeDirectSolver( KaskadeOperator A, const VectorSpace& domain, const VectorSpace& range,
                           ::Kaskade::DirectType directSolver = ::Kaskade::DirectType::UMFPACK3264,
                           ::Kaskade::MatrixProperties property = ::Kaskade::MatrixProperties::GENERAL )
    {
        return DirectSolver< KaskadeOperator, AnsatzVariableSetDescription, TestVariableSetDescription >( std::move( A ), domain, range,
                                                                                                          directSolver, property );
    }
} // namespace Spacy::Kaskade