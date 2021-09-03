#pragma once

#include <Spacy/Adapter/Kaskade/Copy.h>
#include "vectorSpace.hh"
#include "vector.hh"
#include "gridManager.hh"

#include <Spacy/Util/Base/OperatorBase.h>
#include <Spacy/Util/Cast.h>
#include <Spacy/Vector.h>

#include "linalg/direct.hh"
#include "utilities/enums.hh"

namespace Spacy::KaskadeParabolic
{
    /**
     * @ingroup KaskadeGroup
     * @brief Direct solver interface for %Kaskade 7.
     */
    template < class CallableLinearOperator, class AnsatzVariableDescription, class TestVariableDescription >
    class DirectSolver : public OperatorBase
    {
        using Spaces = typename AnsatzVariableDescription::Spaces;
        using Domain = typename AnsatzVariableDescription::template CoefficientVectorRepresentation<>::type;
        using Range = typename TestVariableDescription::template CoefficientVectorRepresentation<>::type;

    public:
        DirectSolver() = delete;
        /**
         * @brief Constructor.
         * @param A CallableLinearOperator with method get(t), to get a KaskadeOperator for each timestep
         * @param spaces boost fusion forward sequence of space pointers required to initialize
         * temporary storage
         * @param domain domain space of the solver
         * @param range range space of the solver
         * @param directSolver solver type (DirectType::MUMPS (default), DirectType::UMFPACK,
         * DirectType::UMFPACK3264 or DirectType::SUPERLU)
         * @param property matrix property (MatrixProperties::GENERAL (default) or
         * MatrixProperties::SYMMETRIC)
         */
        DirectSolver( const CallableLinearOperator& A,
                      const GridManager< Spaces >& gm,
                      const VectorSpace& domain, const VectorSpace& range,
                      ::DirectType directSolver = ::DirectType::UMFPACK3264,
                      ::MatrixProperties property = ::MatrixProperties::GENERAL )
            : OperatorBase( domain, range ), spaces_( *gm.getSpacesVec().at(0) ),
              directSolver_( directSolver ), property_( property ), tEnd_(gm.getTempGrid().getDtVec().size())
        {
            using KMAT = ::Kaskade::MatrixAsTriplet<double>;
            for (int t = 0; t < tEnd_; t++)
            {
                solver_.push_back(std::make_shared< ::Kaskade::InverseLinearOperator< ::Kaskade::DirectSolver< Domain, Range > > >(
                    ::Kaskade::directInverseOperator( ::Kaskade::MatrixRepresentedOperator<KMAT,Domain,Range>(A.template get<KMAT>(t)),
                                                      directSolver_, property_ ) ) );
            }
        }

        /// Compute \f$A^{-1}x\f$.
        ::Spacy::Vector operator()( const ::Spacy::Vector& x ) const
        {
            auto y = zero( range() );
            
            for (int t = 0; t < tEnd_; t++)
            {
                Range y_( TestVariableDescription::template CoefficientVectorRepresentation<>::init( spaces_ ) );
                Domain x_( AnsatzVariableDescription::template CoefficientVectorRepresentation<>::init( spaces_ ) );

                auto & xv = cast_ref< Vector< AnsatzVariableDescription > >( x ).get(t);
                Spacy::Kaskade::copy< AnsatzVariableDescription >( xv, x_ );
                
                solver_.at(t)->apply( x_, y_ );

                auto& yv = cast_ref< Vector< TestVariableDescription > >( y ).get_nonconst(t);
                Spacy::Kaskade::copy< TestVariableDescription >( y_, yv );
            }
            
            return y;
        }
        
        ::Spacy::Vector applyFor( const ::Spacy::Vector& x, int tBegin, int tEnd ) const
        {
            auto y = zero( range() );
            
            for (int t = tBegin; t < tEnd; t++)
            {
                Range y_( TestVariableDescription::template CoefficientVectorRepresentation<>::init( spaces_ ) );
                Domain x_( AnsatzVariableDescription::template CoefficientVectorRepresentation<>::init( spaces_ ) );

                auto & xv = cast_ref< Vector< AnsatzVariableDescription > >( x ).get(t);
                Spacy::Kaskade::copy< AnsatzVariableDescription >( xv, x_ );
                
                solver_.at(t)->apply( x_, y_ );

                auto& yv = cast_ref< Vector< TestVariableDescription > >( y ).get_nonconst(t);
                Spacy::Kaskade::copy< TestVariableDescription >( y_, yv );
            }
            
            return y;
        }

    private:
        Spaces spaces_;
        ::DirectType directSolver_ = ::DirectType::UMFPACK3264;
        ::MatrixProperties property_ = ::MatrixProperties::GENERAL;
        std::vector< std::shared_ptr< ::Kaskade::InverseLinearOperator< ::Kaskade::DirectSolver< Domain, Range > > > > solver_;
        int tEnd_;
    };

    /**
     * @brief Convenient generation of direct solver for %Kaskade 7.
     * @param A CallableLinearOperator with method get(t), to get a KaskadeOperator for each timestep
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
    template < class AnsatzVariableSetDescription, class TestVariableSetDescription, class CallableLinearOperator >
    auto makeDirectSolver( const CallableLinearOperator& A, const GridManager< typename AnsatzVariableSetDescription::Spaces >& gm,
                           const VectorSpace& domain, const VectorSpace& range,
                           ::DirectType directSolver = ::DirectType::UMFPACK3264,
                           ::MatrixProperties property = ::MatrixProperties::GENERAL )
    {
        return DirectSolver< CallableLinearOperator, AnsatzVariableSetDescription, TestVariableSetDescription >( A, gm, domain, range, directSolver, property );
    }
} // namespace Spacy::KaskadeParabolic
