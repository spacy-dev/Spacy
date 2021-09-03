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
            public Mixin::Eps,
            public Mixin::MaxSteps,
            public Mixin::RelativeAccuracy,
            public Mixin::Verbosity
    {
    public:
        SurrogateSolver(const CallableLinearOperator& A, Prec P, const std::vector< Real > dtVec, bool forward = true)
        : A_(A), P_(P), dtVec_(dtVec), forward_(forward)
        {
            
        }
        
        Spacy::Vector operator()(const Spacy::Vector& b) const
        {
            checkSpaceCompatibility(P_.domain(),b.space());
            
            auto x = zero(A_.domain());
    
            if (forward_)
            {
                for (int t = 1; t < dtVec_.size(); t++)
                {
                    Real dt = dtVec_.at(t);

                    auto& xt = cast_ref< Spacy::KaskadeParabolic::Vector< DomainDescription > >( x ).get_nonconst(t);
                    auto xlt = cast_ref< Spacy::KaskadeParabolic::Vector< DomainDescription > >( x ).get(t-1);
                    
                    auto rhs_ = zero(b.space());
                    auto& rhs_t = cast_ref< Spacy::KaskadeParabolic::Vector< RangeDescription > >( rhs_ ).get_nonconst(t);
                    
                    auto Ax = A_.applyFor(x,t-1,t);
                    auto Axt = cast_ref< Spacy::KaskadeParabolic::Vector< RangeDescription > >( Ax ).get(t-1);
                    auto bt = cast_ref< Spacy::KaskadeParabolic::Vector< RangeDescription > >( b ).get(t);
                    
                    
                    rhs_t = bt;
                    rhs_t -= Axt;
                    rhs_ *= dt.get();
                    
                    auto rhs = P_.applyFor(rhs_,t,t+1);
                    
                    //x_{t} = x_{t-1} + P( dt*b - dt*A x_{t-1} )
                    xt = xlt;
                    x += rhs;
                }
            }
            else
            {
                for (int t = dtVec_.size()-2; t >= 0; t--)
                {
                    Real dt = dtVec_.at(t+1);

                    auto& xt = cast_ref< Spacy::KaskadeParabolic::Vector< DomainDescription > >( x ).get_nonconst(t);
                    auto xlt = cast_ref< Spacy::KaskadeParabolic::Vector< DomainDescription > >( x ).get(t+1);
                    
                    auto rhs_ = zero(b.space());
                    auto& rhs_t = cast_ref< Spacy::KaskadeParabolic::Vector< RangeDescription > >( rhs_ ).get_nonconst(t);
                    
                    auto Ax = A_.applyFor(x,t+1,t+2);
                    auto Axt = cast_ref< Spacy::KaskadeParabolic::Vector< RangeDescription > >( Ax ).get(t+1);
                    auto bt = cast_ref< Spacy::KaskadeParabolic::Vector< RangeDescription > >( b ).get(t);
                    
                    
                    rhs_t = bt;
                    rhs_t -= Axt;
                    rhs_ *= dt.get();
                    
                    auto rhs = P_.applyFor(rhs_,t,t+1);
                    
                    //x_{t} = x_{t+1} + P( dt*b - dt*A x_{t+1} )
                    xt = xlt;
                    x += rhs;
                }
            }
            return x;
        }
        
    private:
        CallableLinearOperator A_;
        Prec P_;
        std::vector< Real > dtVec_{};
        bool forward_;
    };
    
    template<class DomainDescription, class RangeDescription, class CallableLinearOperator, class Prec>
    auto makeSurrogatePDESolver( const CallableLinearOperator& A, const Prec& P, const std::vector< Real > dtVec, bool forward = true )
    {
        return SurrogateSolver<DomainDescription,RangeDescription,CallableLinearOperator,Prec>(A,P,dtVec,forward);
    }

} //Spacy::KaskadeParabolic
