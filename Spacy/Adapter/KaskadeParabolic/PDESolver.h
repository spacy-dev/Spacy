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
    template <class DomainDescription, class RangeDescription, class MassMatrix, class CallableLinearOperator, class Prec>
    class Solver :
            public OperatorBase,
            public Mixin::Eps,
            public Mixin::MaxSteps,
            public Mixin::RelativeAccuracy,
            public Mixin::Verbosity
    {
        using Spaces = typename DomainDescription::Spaces;
        
    public:
        Solver(const MassMatrix& D, const CallableLinearOperator& A, const Prec& P, const std::vector< Real > dtVec, const GridManager< Spaces >& gm, Spacy::Real acc, bool forward = true)
        : OperatorBase(A.domain(),A.range()), RelativeAccuracy(acc), D_(D), A_(A), P_(P), gm_(gm), dtVec_(dtVec), forward_(forward), iterations_(0)
        { }
        
        Spacy::Vector operator()(const Spacy::Vector& b) const
        {
            checkSpaceCompatibility(range(),b.space());
            
            auto x = zero(domain());
            
            /// Initialize Direct Solver
            auto Oph = Op_handler(D_,A_,dtVec_);
            auto DirectSolver = makeDirectSolver<RangeDescription,DomainDescription> (Oph, gm_, range(), domain());
    
            int local_iterations_ = 0;
            
            if (forward_)
            {
                for (int t = 1; t < dtVec_.size(); t++) //assume x(0) = 0
                {
                    double dt = dtVec_.at(t).get();
                    
                    /// Initialize CG-Solver
//                     Spacy::CG::NoRegularization NR;
//                     auto CGSolver = Spacy::CG::Solver(Op<DomainDescription,RangeDescription>(D_,A_,dt,t,dtVec_.size()), Prec_Op<RangeDescription,DomainDescription>(P_,t,dtVec_.size()), NR, true);
//                     CGSolver.setEps(eps());
//                     CGSolver.setAbsoluteAccuracy(eps());
//                     CGSolver.setRelativeAccuracy(getRelativeAccuracy());
//                     CGSolver.setMaxSteps(getMaxSteps());
                    
                    
                    auto rhs_ = zero(range());
                    auto& rhs_t = cast_ref< Spacy::KaskadeParabolic::Vector< RangeDescription > >( rhs_ ).get_nonconst(t);
                    
                    auto bt = cast_ref< Spacy::KaskadeParabolic::Vector< RangeDescription > >( b ).get(t);
                    
                    auto Dx = D_.applyFor(x,t-1,t);
                    auto Dxt = cast_ref< Spacy::KaskadeParabolic::Vector< RangeDescription > >( Dx ).get(t-1);
                    
                    
                    rhs_t = bt;
                    rhs_ *= dt;
                    rhs_t += Dxt;
                    
                    
                    // (D+dt*A) x_{t} = D x_{t-1} + dt*b
//                     x += CGSolver(rhs_);
                    x += DirectSolver.applyFor(rhs_,t,t+1);
                    
                    
//                     local_iterations_ += CGSolver.getIterations();
                }
//                 std::cout << std::endl << "StateSolver finished with " << local_iterations_ << " Iterations." << std::endl << std::endl;
            }
            else
            {
                for (int t = dtVec_.size()-2; t >= 0; t--) //assume x(dtVec_.size()-1) = 0
                {
                    double dt = dtVec_.at(t+1).get();

                    /// Initialize CG-Solver
//                     Spacy::CG::NoRegularization NR;
//                     auto CGSolver = Spacy::CG::Solver(Op<DomainDescription,RangeDescription>(D_,A_,dt,t,dtVec_.size()), Prec_Op<RangeDescription,DomainDescription>(P_,t,dtVec_.size()), NR, true);
//                     CGSolver.setEps(eps());
//                     CGSolver.setAbsoluteAccuracy(eps());
//                     CGSolver.setRelativeAccuracy(getRelativeAccuracy());
//                     CGSolver.setMaxSteps(getMaxSteps());
                    
                    
                    auto rhs_ = zero(range());
                    auto& rhs_t = cast_ref< Spacy::KaskadeParabolic::Vector< RangeDescription > >( rhs_ ).get_nonconst(t);
                    
                    auto Dx = D_.applyFor(x,t+1,t+2);
                    auto Dxt = cast_ref< Spacy::KaskadeParabolic::Vector< RangeDescription > >( Dx ).get(t+1);
                    auto bt = cast_ref< Spacy::KaskadeParabolic::Vector< RangeDescription > >( b ).get(t);
                    
                    
                    rhs_t = bt;
                    rhs_ *= dt;
                    rhs_t += Dxt;
                    
                    
                    // (D+dt*A) x_{t} = D x_{t+1} + dt*b
//                     x += CGSolver(rhs_);
                    x += DirectSolver.applyFor(rhs_,t,t+1);
                    
                    
//                     local_iterations_ += CGSolver.getIterations();
                }
//                 std::cout << std::endl << "AdjointSolver finished with " << local_iterations_ << " Iterations." << std::endl << std::endl;
            }
            
            iterations_ = local_iterations_;
            
            return x;
        }
        
        unsigned getIterations() const
        {
            return iterations_;
        }
        
    private:
        MassMatrix D_;
        CallableLinearOperator A_;
        Prec P_;
        std::vector< Real > dtVec_{};
        const GridManager< typename DomainDescription::Spaces > gm_;
        bool forward_;
        mutable int iterations_;
        
        struct Op_handler 
        {
            Op_handler (const MassMatrix& D, const CallableLinearOperator& A, const std::vector< Real > dtVec)
            : D_(D), A_(A), dtVec_(dtVec)
            { }
            
            template <class Matrix>
            Matrix get(int t) const
            {
                Matrix D = D_.get(t).getTriplet();
                Matrix A = A_.template get<Matrix>(t);
                A *= dtVec_.at(t).get();
                D += A;
                return D;
            }
            
            const MassMatrix D_;
            const CallableLinearOperator A_;
            const std::vector< Real > dtVec_;
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
                    auto& yt = cast_ref< Spacy::KaskadeParabolic::Vector< RangeDesc > >( y ).get_nonconst(t);
                    auto xt = cast_ref< Spacy::KaskadeParabolic::Vector< DomainDesc > >( x ).get(t);
                    
                    yt = xt;
                }
                
                for (int t = t_+1; t < tEnd_; t++)
                {
                    auto& yt = cast_ref< Spacy::KaskadeParabolic::Vector< RangeDesc > >( y ).get_nonconst(t);
                    auto xt = cast_ref< Spacy::KaskadeParabolic::Vector< DomainDesc > >( x ).get(t);
                    
                    yt = xt;
                }
                
                
                return y;
            }
            
            MassMatrix D_;
            CallableLinearOperator A_;
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
                    auto& yt = cast_ref< Spacy::KaskadeParabolic::Vector< RangeDesc > >( y ).get_nonconst(t);
                    auto xt = cast_ref< Spacy::KaskadeParabolic::Vector< DomainDesc > >( x ).get(t);
                    
                    yt = xt;
                }
                
                for (int t = t_+1; t < tEnd_; t++)
                {
                    auto& yt = cast_ref< Spacy::KaskadeParabolic::Vector< RangeDesc > >( y ).get_nonconst(t);
                    auto xt = cast_ref< Spacy::KaskadeParabolic::Vector< DomainDesc > >( x ).get(t);
                    
                    yt = xt;
                }
                
                return y;
            }
            
            Prec P_;
            int t_;
            int tEnd_;
        };
    };
    
    template<class MassMatrix, class CallableLinearOperator, class Prec>
    auto makePDESolver(const MassMatrix& D, const CallableLinearOperator& A, const Prec& P, const std::vector< Real > dtVec, const GridManager< typename Prec::RangeDesc::Spaces > gm, Spacy::Real acc, bool forward = true)
    {
        return Solver< typename Prec::RangeDesc, typename Prec::DomainDesc, MassMatrix, CallableLinearOperator, Prec > ( D, A, P, dtVec, gm, acc, forward );
    }
} //Spacy::KaskadeParabolic
