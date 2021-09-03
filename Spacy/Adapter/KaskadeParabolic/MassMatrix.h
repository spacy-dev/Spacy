#pragma once

// #include "Copy.h"
// #include "DirectSolver.h"
// #include "OperatorSpace.h"
// #include "Util/NumaMatrixAsOperator.h"
#include "linalg/triplet.hh"
#include "fem/assemble.hh"
#include "fem/istlinterface.hh"

#include "c2Functional.hh"

#include <Spacy/LinearSolver.h>
#include <Spacy/Util/Base/OperatorBase.h>
#include <Spacy/Util/Base/VectorBase.h>
#include <Spacy/Util/Exceptions.h>
#include <Spacy/Vector.h>

#include <dune/istl/operators.hh>

#include <utility>

/** @addtogroup KaskadeParabolicGroup
 * @{
 */
namespace Spacy::KaskadeParabolic
{
    /** @} */
    template < class DomainDescription, class RangeDescription, class Functional>
    class MassMatrixLinearOperator
    : public OperatorBase
    {
        using DomainSpaces = typename DomainDescription::Spaces;
        using RangeSpaces = typename RangeDescription::Spaces;
        using Variables = typename DomainDescription::Variables;
        using Domain = typename DomainDescription::
        template CoefficientVectorRepresentation<>::type;
        using Range = typename RangeDescription::
        template CoefficientVectorRepresentation<>::type;
        using Matrix = ::Kaskade::MatrixAsTriplet< double >;
        using OperatorMatrixRepresentation = ::Kaskade::MatrixRepresentedOperator< Matrix, Domain, Range >;
        
        
        using DomainVector = Vector<DomainDescription>;
        using RangeVector = Vector<RangeDescription>;
        
    public:
        /**
            * @brief Construct linear operator for %Kaskade 7.
            * @param A operator from %Kaskade 7
            * @param domain domain space
            * @param range range space
            */
        MassMatrixLinearOperator(Functional& f, const VectorSpace& domain, const VectorSpace& range, int rB, int rE, int dB, int dE, int tEnd )
        : OperatorBase( domain,range),
        domainspaces_( creator< VectorCreator<DomainDescription> >(domain).getSpace(0) ),
        rangespaces_( creator< VectorCreator<RangeDescription> >( range ).getSpace(0) ), tEnd_(tEnd)
        { 
            for (int t = 0; t < tEnd_; t++)
            {
                A_.push_back( OperatorMatrixRepresentation(f.assemblerSP(t).template get<Matrix>(false,rB,rE,dB,dE)) );
            }
            
        }
        
        ::Spacy::Vector operator()( const ::Spacy::Vector& x ) const
        {
            checkSpaceCompatibility(x.space(),domain());
            
            auto y = zero( range() );
            
            for (int t=0; t<tEnd_; t++)
            {
                Domain x_( DomainDescription::template CoefficientVectorRepresentation<>::init(domainspaces_) );
                Range y_( RangeDescription::template CoefficientVectorRepresentation<>::init(rangespaces_) );
                
                boost::fusion::at_c< 0 >( x_.data ) = cast_ref< DomainVector >(x).getCoeffVec(t);
                                
                A_.at(t).apply( x_, y_ );
                
                cast_ref< RangeVector >(y).getCoeffVec_nonconst(t) = boost::fusion::at_c< 0 >( y_.data );
            }

            return y;
        }
        
        ::Spacy::Vector applyFor( const ::Spacy::Vector& x, int tBegin, int tEnd ) const
        {
            checkSpaceCompatibility(x.space(),domain());
            
            auto y = zero( range() );
            
            for (int t=tBegin; t<tEnd; t++)
            {
                Domain x_( DomainDescription::template CoefficientVectorRepresentation<>::init(domainspaces_) );
                Range y_( RangeDescription::template CoefficientVectorRepresentation<>::init(rangespaces_) );
                
                boost::fusion::at_c< 0 >( x_.data ) = cast_ref< DomainVector >(x).getCoeffVec(t);
                                
                A_.at(t).apply( x_, y_ );
                
                cast_ref< RangeVector >(y).getCoeffVec_nonconst(t) = boost::fusion::at_c< 0 >( y_.data );
            }

            return y;
        }
        
        auto& get(int t)
        {
            return A_.at(t);
        }
        
        const auto& get(int t) const
        {
            return A_.at(t);
        }
        
    private:
        
        std::vector<OperatorMatrixRepresentation> A_;
        DomainSpaces domainspaces_;
        RangeSpaces rangespaces_;
        int tEnd_;
    };
    
    template< class DomainDescription,class RangeDescription, class Functional >
    auto makeMassMatrixBlock(Functional& f,const VectorSpace& domain, const VectorSpace& range)
    {
        auto& d=domain.getEmbeddingIn(f.domain());
        auto& r=range.getEmbeddingIn(f.domain());
        return MassMatrixLinearOperator<DomainDescription,RangeDescription,Functional>
                (f,domain,range,r.getCBegin(),r.getCEnd(),d.getCBegin(),d.getCEnd(),f.getGridMan().getTempGrid().getDtVec().size());
    }
    
} // namespace Spacy::Kaskade
/** @} */
