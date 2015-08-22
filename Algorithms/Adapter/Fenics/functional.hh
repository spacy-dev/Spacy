#ifndef ALGORITHMS_ADAPTER_FENICS_FUNCTIONAL_HH
#define ALGORITHMS_ADAPTER_FENICS_FUNCTIONAL_HH

#include <chrono>
#include <memory>
#include <vector>

#include "FunctionSpaces/ProductSpace/productSpaceElement.hh"
#include "Interface/abstractFunctional.hh"
#include "Interface/hessian.hh"
#include "Util/Mixins/disableAssembly.hh"
#include "Util/Mixins/primalDualSwitch.hh"
#include "Util/castTo.hh"

#include "../../vector.hh"
#include "../../vectorSpace.hh"
#include "../../functional.hh"

#include "util.hh"
#include "vector.hh"
#include "assignXIfPresent.hh"

namespace Algorithm
{
  namespace Fenics
  {
    template <class F, class DF, class DDF>
    class Functional : public Interface::AbstractFunctional , public Mixin::DisableAssembly
    {
    public:
      Functional(const F& f, const DF& J, const DDF& H,
                   const std::vector<const dolfin::DirichletBC*>& bcs, ::Algorithm::VectorSpace* space)
        : Interface::AbstractFunctional( space ),
          f_( J.function_space(0)->mesh() ),
          J_( J.function_space(0) ),
          H_( H.function_space(0) , H.function_space(1) ),
          bcs_( bcs ),
          dummy_(J_.function_space(0))
      {
        copyCoefficients(f,f_);
        copyCoefficients(J,J_);
        copyCoefficients(H,H_);
      }

      Functional(const F& f, const DF& J, const DDF& H,
                   const std::vector<const dolfin::DirichletBC*>& bcs, ::Algorithm::VectorSpace& space)
        : Functional(f,J,H,bcs,&space)
      {}

      Functional(const F& f, const DF& J, const DDF& H,
                   const std::vector<const dolfin::DirichletBC*>& bcs, ::Algorithm::VectorSpace* space,
                   const dolfin::GenericMatrix& A,
                   const ::Algorithm::Vector& oldX_H)
        : Interface::AbstractFunctional( space ),
          Mixin::DisableAssembly(true),
          f_( J.function_space(0)->mesh() ),
          J_( J.function_space(0) ),
          H_( H.function_space(0) , H.function_space(1) ),
          bcs_( bcs ),
          A_(A.copy()),
          oldX_H_(oldX_H),
          dummy_(J_.function_space(0))
      {
        copyCoefficients(f,f_);
        copyCoefficients(J,J_);
        copyCoefficients(H,H_);
      }


      double d0(const ::Algorithm::Vector& x) const final override
      {
        primalDualIgnoreReset(std::bind(&Functional::assembleFunctional,std::ref(*this), std::placeholders::_1),x);

        return value_;
      }

      ::Algorithm::Vector d1(const ::Algorithm::Vector &x) const final override
      {
        primalDualIgnoreReset(std::bind(&Functional::assembleJacobian,std::ref(*this), std::placeholders::_1),x);

        auto y = domain().dualSpace_ptr()->element();
        copy(*b_,y);
        return y;
      }

      ::Algorithm::Vector d2(const ::Algorithm::Vector &x, const ::Algorithm::Vector &dx) const final override
      {
        primalDualIgnoreReset(std::bind(&Functional::assembleHessian,std::ref(*this), std::placeholders::_1),x);

        auto x_ = std::make_shared<dolfin::Vector>(dummy_.vector()->mpi_comm(), dummy_.vector()->size());
        copy(dx,*x_);
        auto Ax = x_->copy();
        A_->mult(*x_, *Ax);

        auto result = domain().dualSpace_ptr()->element();
        copy(*Ax,result);

        return result;
    }

      void setOrigin(const ::Algorithm::Vector& x) const
      {
        auto x_ = std::make_shared<dolfin::Vector>(dummy_.vector()->mpi_comm(), dummy_.vector()->size());
        copy(x,*x_);
        *dummy_.vector() = *x_;
        assign_x0_if_present(f_,dummy_);
        assign_x0_if_present(J_,dummy_);
        assign_x0_if_present(H_,dummy_);
      }

    private:
      std::unique_ptr<Interface::Hessian> makeHessian(const ::Algorithm::Vector& x) const override
      {
        primalDualIgnoreReset(std::bind(&Functional::assembleHessian,std::ref(*this), std::placeholders::_1),x);

        assert( A_ != nullptr );
        return std::make_unique<Interface::Hessian>( std::make_unique<Functional>(f_,J_,H_,bcs_,domain_ptr(),*A_,oldX_H_), oldX_H_);
      }

      std::unique_ptr<Interface::AbstractLinearSolver> makeSolver() const
      {
        assert( A_ != nullptr );
        return std::make_unique<LUSolver>(A_,*J_.function_space(0),domain_ptr(),domain_ptr());
      }

      void assembleFunctional(const ::Algorithm::Vector& x) const
      {
        if( valueAssembled_ && (oldX_f_ == x) ) return;
        valueAssembled_ = true;
        auto x_ = std::make_shared<dolfin::Vector>(dummy_.vector()->mpi_comm(), dummy_.vector()->size());
        copy(x,*x_);

        *dummy_.vector() = *x_;

        f_.x = dummy_;
        value_ = dolfin::assemble(f_);

        oldX_f_ = x;
      }

      void assembleJacobian(const ::Algorithm::Vector& x) const
      {
        if( assemblyIsDisabled() ) return;

        if( b_ != nullptr && (oldX_J_==x) ) return;

        auto x_ = std::make_shared<dolfin::Vector>(dummy_.vector()->mpi_comm(), dummy_.vector()->size());
        copy(x,*x_);
        *dummy_.vector() = *x_;
        Assign_X_If_Present<decltype(J_)>::apply(J_,dummy_);
        b_ = dummy_.vector()->factory().create_vector();

        dolfin::assemble(*b_,J_);
        for( auto& bc : bcs_) bc->apply(*b_,*dummy_.vector());

        oldX_J_ = x;
      }

      void assembleHessian(const ::Algorithm::Vector& x) const
      {
        if( assemblyIsDisabled() ) return;

        if( A_ != nullptr && (oldX_H_==x) ) return;

        auto x_ = std::make_shared<dolfin::Vector>(dummy_.vector()->mpi_comm(), dummy_.vector()->size());
        copy(x,*x_);
        *dummy_.vector() = *x_;
        Assign_X_If_Present<decltype(H_)>::apply(H_,dummy_);
        A_ = dummy_.vector()->factory().create_matrix();
        dolfin::assemble(*A_,H_);

//        std::cout << "A: " << std::endl;
//        std::cout << A_->str(true) << std::endl;

        for( auto& bc : bcs_) bc->apply(*A_);

        oldX_H_ = x;
      }

      Functional* cloneImpl() const
      {
        if( assemblyIsDisabled() ) return new Functional(f_,J_,H_,bcs_,domain_ptr(),*A_,oldX_H_);
        return new Functional(f_,J_,H_,bcs_,domain_ptr());
      }

      mutable F f_;
      mutable DF J_;
      mutable DDF H_;
      const std::vector<const dolfin::DirichletBC*>& bcs_;
      mutable std::shared_ptr<dolfin::GenericMatrix> A_;
      mutable std::shared_ptr<dolfin::GenericVector> b_;
      mutable double value_ = 0;
      mutable bool valueAssembled_ = false;
      mutable ::Algorithm::Vector oldX_f_, oldX_J_, oldX_H_;
      mutable dolfin::Function dummy_;
    };


    template <class F, class DF, class DDF, class... Args>
    ::Algorithm::Functional makeFunctional( const F& f , const DF& J , const DDF& H , Args&&... args )
    {
      return createFromUniqueImpl< ::Algorithm::Functional , Fenics::Functional<F,DF,DDF> >( f , J , H , std::forward<Args>(args)...);
    }
  }
}

#endif // ALGORITHMS_ADAPTER_FENICS_FUNCTIONAL_HH

