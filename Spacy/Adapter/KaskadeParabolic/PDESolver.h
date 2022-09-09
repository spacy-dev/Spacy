#pragma once

#include <Spacy/Util/Mixins.h>
#include <Spacy/Util/Base/OperatorBase.h>

#include <Spacy/Operator.h>
#include <Spacy/Vector.h>
#include <Spacy/ZeroVectorCreator.h>

#include <Spacy/Algorithm/CG/CG.h>
#include <Spacy/Algorithm/CG/RegularizeViaPreconditioner.h>

#include "vector.hh"
#include "DirectSolver.h"

namespace Spacy::KaskadeParabolic
{
    template <class DomainDescription, class RangeDescription, class MassMatrix, class CallableLinearOperator, class Prec, class Functional>
    class Solver :
            public OperatorBase,
            public Mixin::Eps,
            public Mixin::MaxSteps,
            public Mixin::RelativeAccuracy,
            public Mixin::Verbosity
    {
        using Spaces = typename DomainDescription::Spaces;
        
    public:
        Solver(const MassMatrix& D, const CallableLinearOperator& A, Prec& P, const Functional& f, Spacy::Real acc, bool forward = true)
        : OperatorBase(A.range(),A.domain()), RelativeAccuracy(acc), D_(D), A_(A), P_(P), f_(f), forward_(forward), iterations_(0)
        { 
        }
        
        Spacy::Vector operator()(const Spacy::Vector& b) const
        {
            checkSpaceCompatibility(domain(),b.space());
            
            using namespace std::chrono;
            auto startTime = high_resolution_clock::now();
            auto endTime = high_resolution_clock::now() - startTime;
            
            auto x = zero(range());
            auto rhs_ = zero(domain());

            auto dtVec_ = f_.getGridMan().getTempGrid().getDtVec();
            auto N = dtVec_.size();
            
            /// Initialize Direct Solver
            auto Oph = Op_handler(D_,A_,dtVec_);
            auto DirectSolver = makeDirectSolver<RangeDescription,DomainDescription> (Oph, f_, domain(), range());

            int local_iterations_ = 0;
            
            if (forward_)
            {
                const auto& bt = cast_ref< Spacy::KaskadeParabolic::Vector< RangeDescription > >( b ).getCoeffVec(1);
                auto& rhs_t = cast_ref< Spacy::KaskadeParabolic::Vector< RangeDescription > >( rhs_ ).getCoeffVec_nonconst(1);
                rhs_t = bt;
                
                //(D+dt*A) x_{1} = dt*b
                x += dtVec_.at(1).get() * DirectSolver.applyFor(rhs_,1,2);
                
                for (int t = 2; t < N; t++)
                {
                    
                    double dt = dtVec_.at(t).get();
                    
                    /// Initialize CG-Solver
//                     Spacy::CG::NoRegularization NR;
//                     auto CGSolver = Spacy::CG::Solver(Op<DomainDescription,RangeDescription>(D_,A_,dt,t,dtVec_.size()), Prec_Op<RangeDescription,DomainDescription>(P_,t,dtVec_.size()), NR, true);
//                     CGSolver.setEps(eps());
//                     CGSolver.setAbsoluteAccuracy(eps());
//                     CGSolver.setRelativeAccuracy(getRelativeAccuracy());
//                     CGSolver.setMaxSteps(getMaxSteps());
                    
                    
                    auto& rhs_t = cast_ref< Spacy::KaskadeParabolic::Vector< RangeDescription > >( rhs_ ).getCoeffVec_nonconst(t);
                    
                    const auto& bt = cast_ref< Spacy::KaskadeParabolic::Vector< RangeDescription > >( b ).getCoeffVec(t);
                    auto Dx = D_.applyFor(x,t-1,t);
                    const auto& Dxt = cast_ref< Spacy::KaskadeParabolic::Vector< RangeDescription > >( Dx ).getCoeffVec(t-1);
                    
                    rhs_t = bt;
                    rhs_ *= dt;
                    rhs_t += Dxt; 
                    
                    //(D+dt*A) x_{t} = D x_{t-1} + dt*b
//                     x += CGSolver(rhs_);
                    x += DirectSolver.applyFor(rhs_,t,t+1);
                    
                    
                    
//                     local_iterations_ += CGSolver.getIterations();
                    local_iterations_++;
                }
            }
            else
            {
                const auto& bt = cast_ref< Spacy::KaskadeParabolic::Vector< RangeDescription > >( b ).getCoeffVec(N-1);
                auto& rhs_t = cast_ref< Spacy::KaskadeParabolic::Vector< RangeDescription > >( rhs_ ).getCoeffVec_nonconst(N-1);
                rhs_t = bt;
                
                // (D+dt*A) x_{N-1} = dt*b
                x += dtVec_.at(N-1).get() * DirectSolver.applyFor(rhs_,N-1,N);
                
                for (int t = N-2; t >= 1; t--)
                {
                    double dt = dtVec_.at(t).get();

                    /// Initialize CG-Solver
//                     Spacy::CG::NoRegularization NR;
//                     auto CGSolver = Spacy::CG::Solver(Op<DomainDescription,RangeDescription>(D_,A_,dt,t,dtVec_.size()), Prec_Op<RangeDescription,DomainDescription>(P_,t,dtVec_.size()), NR, true);
//                     CGSolver.setEps(eps());
//                     CGSolver.setAbsoluteAccuracy(eps());
//                     CGSolver.setRelativeAccuracy(getRelativeAccuracy());
//                     CGSolver.setMaxSteps(getMaxSteps());
                    
                    
                    auto rhs_ = zero(domain());
                    auto& rhs_t = cast_ref< Spacy::KaskadeParabolic::Vector< RangeDescription > >( rhs_ ).getCoeffVec_nonconst(t);
                    auto Dx = D_.applyFor(x,t+1,t+2);
                    const auto& Dxt = cast_ref< Spacy::KaskadeParabolic::Vector< RangeDescription > >( Dx ).getCoeffVec(t+1);
                    const auto& bt = cast_ref< Spacy::KaskadeParabolic::Vector< RangeDescription > >( b ).getCoeffVec(t);
                    
                    rhs_t = bt;
                    rhs_ *= dt;
                    rhs_t += Dxt;
                    
                    // (D+dt*A) x_{t} = D x_{t+1} + dt*b
//                     x += CGSolver(rhs_);
                    x += DirectSolver.applyFor(rhs_,t,t+1);
                    
                    
//                     local_iterations_ += CGSolver.getIterations();
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
        const MassMatrix& D_;
        const CallableLinearOperator& A_;
        Prec& P_;
        const Functional& f_;
        bool forward_;
        mutable int iterations_;
        
        struct Op_handler 
        {
            
            Op_handler (const MassMatrix& D, const CallableLinearOperator& A, const std::vector< Real >& dtVec)
            : D_(D), A_(A), dtVec_(dtVec)
            { }
            
            template <class Matrix>
            Matrix get(int t) const
            {
                Matrix D = D_.template get<Matrix>(t);
                Matrix A = A_.template get<Matrix>(t);
                A *= dtVec_.at(t).get();
                D += A;
                return D;
            }
            
            const MassMatrix& D_;
            const CallableLinearOperator& A_;
            const std::vector< Real >& dtVec_;
        };
        
        template <class DomainDesc, class RangeDesc>
        struct Op : OperatorBase
        {
            Op (const MassMatrix& D, const CallableLinearOperator& A, double dt, int t, int tEnd)
            : D_(D), A_(A), dt_(dt), t_(t), tEnd_(tEnd), OperatorBase(A.domain(),A.range())
            { }
            
            Spacy::Vector operator()(const Spacy::Vector& x) const
            {
                checkSpaceCompatibility(domain(),x.space());
                
                auto y   = D_.applyFor(x,t_,t_+1);
                y       += dt_ * A_.applyFor(x,t_,t_+1);
                
                for (int t = 0; t < t_; t++)
                {
                    auto& yt = cast_ref< Spacy::KaskadeParabolic::Vector< RangeDesc > >( y ).getCoeffVec_nonconst(t);
                    auto xt = cast_ref< Spacy::KaskadeParabolic::Vector< DomainDesc > >( x ).getCoeffVec(t);
                    
                    yt = xt;
                }
                
                for (int t = t_+1; t < tEnd_; t++)
                {
                    auto& yt = cast_ref< Spacy::KaskadeParabolic::Vector< RangeDesc > >( y ).getCoeffVec_nonconst(t);
                    auto xt = cast_ref< Spacy::KaskadeParabolic::Vector< DomainDesc > >( x ).getCoeffVec(t);
                    
                    yt = xt;
                }
                
                
                return y;
            }
            
            const MassMatrix& D_;
            const CallableLinearOperator& A_;
            double dt_;
            int t_;
            int tEnd_;
        };
        
        template <class DomainDesc, class RangeDesc>
        struct Prec_Op : OperatorBase
        {
            Prec_Op (const Prec& P, int t, int tEnd)
            : P_(P), t_(t), tEnd_(tEnd), OperatorBase(P.domain(),P.range())
            { }
            
            Spacy::Vector operator()(const Spacy::Vector& x) const
            {
                checkSpaceCompatibility(domain(),x.space());
                
                auto y = P_.applyFor(x,t_,t_+1);
                
                for (int t = 0; t < t_; t++)
                {
                    auto& yt = cast_ref< Spacy::KaskadeParabolic::Vector< RangeDesc > >( y ).getCoeffVec_nonconst(t);
                    auto xt = cast_ref< Spacy::KaskadeParabolic::Vector< DomainDesc > >( x ).getCoeffVec(t);
                    
                    yt = xt;
                }
                
                for (int t = t_+1; t < tEnd_; t++)
                {
                    auto& yt = cast_ref< Spacy::KaskadeParabolic::Vector< RangeDesc > >( y ).getCoeffVec_nonconst(t);
                    auto xt = cast_ref< Spacy::KaskadeParabolic::Vector< DomainDesc > >( x ).getCoeffVec(t);
                    
                    yt = xt;
                }
                
                return y;
            }
            
            const Prec& P_;
            int t_;
            int tEnd_;
        };
    };
    
    template<class DomainDescription, class RangeDescription, class MassMatrix, class CallableLinearOperator, class Prec, class Functional>
    auto makePDESolver(const MassMatrix& D, const CallableLinearOperator& A, Prec& P, const Functional& f, Spacy::Real acc, bool forward = true)
    {
        return Solver< DomainDescription, RangeDescription, MassMatrix, CallableLinearOperator, Prec, Functional > ( D, A, P, f, acc, forward );
    }
} //Spacy::KaskadeParabolic
