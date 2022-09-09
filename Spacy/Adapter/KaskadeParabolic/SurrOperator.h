#pragma once

#include <Spacy/Util/Base/OperatorBase.h>

#include <Spacy/Vector.h>
#include <Spacy/VectorSpace.h>
#include <Spacy/Operator.h>

#include "gridManager.hh"

/** @addtogroup KaskadeParabolicGroup
 * @{
 */
namespace Spacy::KaskadeParabolic
{
    template <class DomainDescription, class RangeDescription, class CallableLinearOperator, class Prec,class Diag>
    class SurrOperator : 
            public OperatorBase,
            public Mixin::Eps,
            public Mixin::MaxSteps,
            public Mixin::RelativeAccuracy,
            public Mixin::Verbosity
    {
   
    public:
        typedef DomainDescription Domain;
        typedef RangeDescription Range;
        
        SurrOperator(const CallableLinearOperator& A, Prec& P, Diag& D, const Spacy::VectorSpace& domain, const Spacy::VectorSpace& range, bool forward = true)
        : OperatorBase(domain,range), A_(A), P_(P), D_(D), iterations_(0), forward_(forward)
        {
        }
        
        ::Spacy::Vector operator()( const ::Spacy::Vector& x ) const
        {
            checkSpaceCompatibility(x.space(),domain());
            
            iterations_++;
            
            if(!forward_)
            {
                return D_(x);
            }
            else
            {
                return D_.transposed(x);
            }
        }
        
        unsigned getIterations() const
        {
            return 0;//iterations_;
        }
        
        void updateDiag( const ::Spacy::Vector& y, const ::Spacy::Vector& x )
        {
            checkSpaceCompatibility(x.space(),domain());
            checkSpaceCompatibility(y.space(),range());
            
            D_.updateDiag(y,x);
        }
        
    private:
        const CallableLinearOperator& A_;
        Prec& P_;
        Diag& D_;
        mutable unsigned iterations_;
        const bool forward_;
    };
    
    template<class DomainDescription, class RangeDescription, class CallableLinearOperator, class Prec, class Diag>
    auto makeSurrogateOperator(const CallableLinearOperator& A, Prec& P, Diag& D, const Spacy::VectorSpace& domain, const Spacy::VectorSpace& range, bool forward = true)
    {
        return SurrOperator<DomainDescription,RangeDescription,CallableLinearOperator,Prec,Diag>(A,P,D,domain,range,forward);
    }
    
} // namespace Spacy::KaskadeParabolic
/** @} */
