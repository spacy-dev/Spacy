#pragma once

#include <Spacy/Util/Mixins.h>
#include <Spacy/Operator.h>
#include <Spacy/Vector.h>
#include <Spacy/ZeroVectorCreator.h>

#include "vector.hh"


namespace Spacy::KaskadeParabolic
{
    template <class DomainDescription, class RangeDescription, class CallableLinearOperator, class Prec>
    class SurrogateSolver :
            public OperatorBase,
            public Mixin::Eps,
            public Mixin::MaxSteps,
            public Mixin::RelativeAccuracy,
            public Mixin::Verbosity
    {
        
    public:
        SurrogateSolver(const CallableLinearOperator& A, Prec& P, const GridManager< typename DomainDescription::Spaces >& gm, bool forward = true)
        : OperatorBase(A.range(),A.domain()), A_(A), P_(P), gm_(gm), forward_(forward), iterations_(0)
        {
        }
        
        Spacy::Vector operator()(const Spacy::Vector& b) const
        {
            checkSpaceCompatibility(domain(),b.space());
            
            auto x = zero(range());
            auto rhs_ = zero(domain());
            
            auto dtVec_ = gm_.getTempGrid().getDtVec();
            auto N = dtVec_.size();
            
            int local_iterations_ = 0;
    
            if (forward_)
            {
                for (int t = 1; t < N; t++)
                {   
                    Real dt = dtVec_.at(t);

                    auto& xt = cast_ref< Spacy::KaskadeParabolic::Vector< DomainDescription > >( x ).getCoeffVec_nonconst(t);
                    auto xlt = cast_ref< Spacy::KaskadeParabolic::Vector< DomainDescription > >( x ).getCoeffVec(t-1);
                    
                    auto& rhs_t = cast_ref< Spacy::KaskadeParabolic::Vector< RangeDescription > >( rhs_ ).getCoeffVec_nonconst(t);
                    
                    auto Ax = A_.applyFor(x,t-1,t); 
                    const auto& Axt = cast_ref< Spacy::KaskadeParabolic::Vector< RangeDescription > >( Ax ).getCoeffVec(t-1);
                    const auto& bt = cast_ref< Spacy::KaskadeParabolic::Vector< RangeDescription > >( b ).getCoeffVec(t);
                    
                    rhs_t = bt;
                    rhs_t -= Axt;
                    
                    auto rhs = dt.get() * P_.applyFor(rhs_,t,t+1);
                    const auto& rhst = cast_ref< Spacy::KaskadeParabolic::Vector< RangeDescription > >( rhs ).getCoeffVec(t);
                    
                    //x_{t} = x_{t-1} + dt * P( b - A x_{t-1} )
                    xt = xlt;
                    xt += rhst;
                    
                    
                    local_iterations_++;
                }
            }
            else
            {
                const auto& bt = cast_ref< Spacy::KaskadeParabolic::Vector< RangeDescription > >( b ).getCoeffVec(N-1);
                auto& rhs_t = cast_ref< Spacy::KaskadeParabolic::Vector< RangeDescription > >( rhs_ ).getCoeffVec_nonconst(N-1);
                rhs_t = bt;
                
                //x_{N-1} = dt * P(b)
                x += dtVec_.at(N-1).get() * P_.applyFor(rhs_,N-1,N);
                
                for (int t = N-2; t >= 1; t--)
                {
                    
                    Real dt = dtVec_.at(t);

                    auto& xt = cast_ref< Spacy::KaskadeParabolic::Vector< DomainDescription > >( x ).getCoeffVec_nonconst(t);
                    auto xlt = cast_ref< Spacy::KaskadeParabolic::Vector< DomainDescription > >( x ).getCoeffVec(t+1);
                    
                    auto& rhs_t = cast_ref< Spacy::KaskadeParabolic::Vector< RangeDescription > >( rhs_ ).getCoeffVec_nonconst(t);
                    
                    auto Ax = A_.applyFor(x,t+1,t+2);
                    const auto& Axt = cast_ref< Spacy::KaskadeParabolic::Vector< RangeDescription > >( Ax ).getCoeffVec(t+1);
                    
                    const auto& bt = cast_ref< Spacy::KaskadeParabolic::Vector< RangeDescription > >( b ).getCoeffVec(t);
                    
                    rhs_t = bt;
                    rhs_t -= Axt;
                    
                    auto rhs = dt.get() * P_.applyFor(rhs_,t,t+1);
                    const auto& rhst = cast_ref< Spacy::KaskadeParabolic::Vector< RangeDescription > >( rhs ).getCoeffVec(t);
                    
                    //x_{t} = x_{t+1} + dt * P( b - A x_{t+1} )
                    xt = xlt;
                    xt += rhst;
                    
                    
                    local_iterations_++;
                }
            }
            
            iterations_ = local_iterations_;
            
            return x;
        }
        
        unsigned getIterations() const
        {
            return iterations_;
        }
        
    private:
        const CallableLinearOperator& A_;
        Prec& P_;
        const GridManager< typename DomainDescription::Spaces >& gm_;
        bool forward_;
        mutable int iterations_;
    };
    
    template<class DomainDescription, class RangeDescription, class CallableLinearOperator, class Prec>
    auto makeSurrogatePDESolver( const CallableLinearOperator& A, Prec& P, const GridManager< typename RangeDescription::Spaces >& gm, bool forward = true )
    {
        return SurrogateSolver<DomainDescription,RangeDescription,CallableLinearOperator,Prec>(A,P,gm,forward);
    }

} //Spacy::KaskadeParabolic
