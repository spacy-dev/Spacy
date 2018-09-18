/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/Numerik/numsoft/kaskade7/                        */
/*                                                                           */
/*  Copyright (C) 2002-2014 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma once

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include "fem/fixdune.hh"
#include "utilities/linalg/scalarproducts.hh"

#include "fem/variables.hh"
#include "fem/functional_aux.hh"

namespace Kaskade
{
  /**
       * @brief nonlinear constraint for optimal control problems
       */
  template <class AnsatzVars_, int stateId = 0>
  class NonlinearModelConstraint
  {
  public:
    typedef AnsatzVars_ AnsatzVars;
    typedef typename AnsatzVars::Scalar Scalar;
    static constexpr int dim = AnsatzVars::Grid::dimension;

    /**
     * @brief Construct nonlinear constraint
     * @param c (c*y*y + d) Heat conduction tensor (smaller d, larger c) -> harder problem
     * @param d see above
     * @param e nonlinear instability factor
     * @param mu linear instability factor  mu zero -> stable as laplace is stable, mu > 0 -> more stable , mu < 0 -> unstable
     */
    explicit NonlinearModelConstraint(Scalar c_=0, Scalar d_=1,Scalar e_=0, Scalar mu_=0):c(c_),d(d_),e(e_),mu(mu_) { }

    void evaluateAt(Dune::FieldVector<Scalar,AnsatzVars::template Components<stateId>::m> const& y_,
                    Dune::FieldMatrix<Scalar,AnsatzVars::template Components<stateId>::m,dim> const& dy_)
    {
      y = y_;
      dy = dy_;

      kappa   = d+c*y*y;
      dkappa  = c*2*y;
      ddkappa = c*2*unitMatrix<Scalar,AnsatzVars::template Components<stateId>::m>();

      F=0.25*e*y*y*y*y;
      dF=e*y*y*y;
      ddF=3*e*y*y;
      dddF=6*e*y;
    }


    template<int row>
    Scalar d1 (VariationalArg<Scalar,dim, AnsatzVars::template Components<row>::m> const& arg) const
    {
      if( row != stateId ) return 0;
      return kappa*sp(dy,arg.derivative) + mu*y*arg.value + dF*arg.value;

    }

    template<int row, int col>
    Scalar d2 (VariationalArg<Scalar,dim,AnsatzVars::template Components<row>::m> const &arg1, VariationalArg<Scalar,dim,AnsatzVars::template Components<col>::m> const &arg2) const
    {
      if ( row != stateId || col != stateId ) return 0;

      return ( dkappa * arg2.value ) * sp(dy,arg1.derivative)
          + kappa * sp(arg2.derivative[0],arg1.derivative[0]) + mu*arg2.value*arg1.value + ddF*arg1.value*arg2.value;
//      return  d*sp(arg2.derivative,arg1.derivative) + mu*arg2.value*arg1.value;
    }

    template <int id1, int id2, int id3>
    Scalar d3 (VariationalArg<Scalar,dim,AnsatzVars::template Components<id1>::m> const& arg1,
               VariationalArg<Scalar,dim,AnsatzVars::template Components<id2>::m> const& arg2,
               VariationalArg<Scalar,dim,AnsatzVars::template Components<id3>::m> const& arg3) const
    {
      if( id1 != stateId || id2 != stateId || id3 != stateId ) return 0;
      return ( arg3.value * ( ddkappa * arg2.value ) ) * sp(dy,arg1.derivative) + ( dkappa * arg2.value ) * sp(arg3.derivative, arg1.derivative)
          + ( dkappa * arg3.value ) * sp(arg2.derivative,arg1.derivative)+dddF*arg1.value*arg2.value*arg3.value;
    }

  private:
    Scalar c,d,e,mu,kappa,F,dF,ddF,dddF;
    Dune::FieldVector<Scalar,AnsatzVars::template Components<stateId>::m> dkappa, y;
    Dune::FieldMatrix<Scalar,AnsatzVars::template Components<stateId>::m,dim> dy;
    Dune::FieldMatrix<Scalar,AnsatzVars::template Components<stateId>::m,AnsatzVars::template Components<stateId>::m> ddkappa;
    LinAlg::EuclideanScalarProduct sp;
  };

  /**
       * @brief nonlinear pde with nonconstant source for simulation
       */
  template <class RType, class VarSet, class DescriptionControl>
  class NonlinearModelPDE : public FunctionalBase<WeakFormulation>
  {
  public:
    using Scalar = RType;
    using OriginVars = VarSet;
    using AnsatzVars = VarSet;
    using TestVars = VarSet;

    using SourceVars = DescriptionControl;

    static constexpr int dim = AnsatzVars::Grid::dimension;
    static constexpr int uIdx = 0;
    static constexpr int uSpaceIdx = boost::fusion::result_of::value_at_c<typename AnsatzVars::Variables,
                                               uIdx>::type::spaceIndex;

    class DomainCache : public CacheBase<NonlinearModelPDE,DomainCache>
    {
    public:
      DomainCache(NonlinearModelPDE const& h_,
                  typename AnsatzVars::VariableSet const& vars_,
                  int flags=7):
        h(h_), data(vars_)
      {}

      template <class Entity>
      void moveTo(Entity const& entity)
      {
        e = entity;
      }

      template <class Position, class Evaluators>
      void evaluateAt(Position const& x, Evaluators const& evaluators)
      {
        u  = boost::fusion::at_c<uIdx>(data.data).value(boost::fusion::at_c<uSpaceIdx>(evaluators));
        du = boost::fusion::at_c<uIdx>(data.data).derivative(boost::fusion::at_c<uSpaceIdx>(evaluators));
        kappa   = h.d+h.c*u*u;
        dkappa  = h.c*2*u;

        F=0.25*h.e*u*u*u*u;
        dF=h.e*u*u*u;
        ddF=3*h.e*u*u;

        f = boost::fusion::at_c<0>(h.control_.data).value(e ,x);
      }

      Scalar
      d0() const
      {
        std::cout<<"THIS SHOULD NOT BE CALLED"<<std::endl;
        return 0;
      }

      template<int row>
      Scalar d1_impl (VariationalArg<Scalar,dim,TestVars::template Components<row>::m> const& arg) const
      {
        return kappa*sp(du,arg.derivative) + h.mu*u*arg.value + dF*arg.value - f*arg.value;
      }

      template<int row, int col>
      Scalar d2_impl (VariationalArg<Scalar,dim,TestVars::template Components<row>::m> const &arg1,
                    VariationalArg<Scalar,dim,AnsatzVars::template Components<row>::m> const &arg2) const
      {
        return ( dkappa * arg2.value ) * sp(du,arg1.derivative)
            + kappa * sp(arg2.derivative[0],arg1.derivative[0]) + h.mu*arg2.value*arg1.value + ddF*arg1.value*arg2.value;
      }

    private:
      NonlinearModelPDE const& h;
      Scalar kappa,F,dF,ddF;
      typename AnsatzVars::VariableSet const& data;
      Dune::FieldVector<Scalar,AnsatzVars::template Components<uIdx>::m> u, f, dkappa;
      Dune::FieldMatrix<Scalar,AnsatzVars::template Components<uIdx>::m,dim> du;
      LinAlg::EuclideanScalarProduct sp;

      typename AnsatzVars::Grid::template Codim<0>::Entity e;
    };

    class BoundaryCache : public CacheBase<NonlinearModelPDE,BoundaryCache>
    {
    public:
      BoundaryCache(NonlinearModelPDE<RType,AnsatzVars,DescriptionControl> const&,
                    typename AnsatzVars::VariableSet const& vars_,
                    int flags=7):
        data(vars_), penalty(1e9), u(0.), uDirichletBoundaryValue(0.)
      {}

      template <class Position, class Evaluators>
      void evaluateAt(Position const&, Evaluators const& evaluators)
      {
        u = boost::fusion::at_c<uIdx>(data.data).value(boost::fusion::at_c<uSpaceIdx>(evaluators));
      }

      Scalar d0() const
      {
        return penalty*(u-uDirichletBoundaryValue)*(u-uDirichletBoundaryValue)/2;
      }

      template<int row>
      Scalar d1_impl (VariationalArg<Scalar,dim> const& arg) const
      {
        return penalty*(u-uDirichletBoundaryValue)*arg.value;
      }

      template<int row, int col>
      Scalar d2_impl (VariationalArg<Scalar,dim> const &arg1,
                      VariationalArg<Scalar,dim> const &arg2) const
      {
        return penalty*arg1.value*arg2.value;
      }

    private:
      typename AnsatzVars::VariableSet const& data;
      Scalar penalty;
      Dune::FieldVector<Scalar,AnsatzVars::template Components<uIdx>::m> u, uDirichletBoundaryValue;
    };

    /**
     * @brief Construct nonlinear pde
     * @param c (c*y*y + d) Heat conduction tensor (smaller d, larger c) -> harder problem
     * @param d see above
     * @param e nonlinear instability factor
     * @param mu linear instability factor  mu zero -> stable as laplace is stable, mu > 0 -> more stable , mu < 0 -> unstable
     * @param control_ source term of the pde
     */
    explicit NonlinearModelPDE(Scalar c_, Scalar d_, Scalar e_, Scalar mu_, typename DescriptionControl::VariableSet control): c(c_),d(d_),e(e_),mu(mu_),control_(control) {}
    template <int row>
    struct D1: public FunctionalBase<WeakFormulation>::D1<row>
    {
      static bool const present   = true;
      static bool const constant  = false;

    };

    template <int row, int col>
    struct D2: public FunctionalBase<WeakFormulation>::D2<row,col>
    {
      static bool const present = true;
      static bool const symmetric = false;
      static bool const lumped = false;
    };

    template <class Cell>
    int integrationOrder(Cell const& /* cell */, int shapeFunctionOrder, bool boundary) const
    {
      if (boundary)
        return 2*shapeFunctionOrder;
      else
      {
        int stiffnessMatrixIntegrationOrder = 2*(shapeFunctionOrder-1) + 2*shapeFunctionOrder;
        int sourceTermIntegrationOrder = shapeFunctionOrder;        // as rhs f is constant, i.e. of order 0

        return std::max(stiffnessMatrixIntegrationOrder,sourceTermIntegrationOrder);
      }
    }
    Scalar c, d, e, mu;
    typename DescriptionControl::VariableSet control_;

  };

  /**
       * @brief nonlinear boundary controlled pde with nonconstant source for simulation
       */
  template <class RType, class VarSet, class DescriptionControl>
  class BoundaryControlledNonlinearModelPDE : public FunctionalBase<WeakFormulation>
  {
  public:
    using Scalar = RType;
    using OriginVars = VarSet;
    using AnsatzVars = VarSet;
    using TestVars = VarSet;

    using SourceVars = DescriptionControl;

    static constexpr int dim = AnsatzVars::Grid::dimension;
    static constexpr int uIdx = 0;
    static constexpr int uSpaceIdx = boost::fusion::result_of::value_at_c<typename AnsatzVars::Variables,
                                               uIdx>::type::spaceIndex;

    class DomainCache : public CacheBase<BoundaryControlledNonlinearModelPDE,DomainCache>
    {
    public:
      DomainCache(BoundaryControlledNonlinearModelPDE const& h_,
                  typename AnsatzVars::VariableSet const& vars_,
                  int flags=7):
        h(h_), data(vars_)
      {}


      template <class Position, class Evaluators>
      void evaluateAt(Position const& x, Evaluators const& evaluators)
      {
        u  = boost::fusion::at_c<uIdx>(data.data).value(boost::fusion::at_c<uSpaceIdx>(evaluators));
        du = boost::fusion::at_c<uIdx>(data.data).derivative(boost::fusion::at_c<uSpaceIdx>(evaluators));
        kappa   = h.d+h.c*u*u;
        dkappa  = h.c*2*u;

        F=0.25*h.e*u*u*u*u;
        dF=h.e*u*u*u;
        ddF=3*h.e*u*u;
      }

      Scalar
      d0() const
      {
        std::cout<<"THIS SHOULD NOT BE CALLED"<<std::endl;
        return 0;
      }

      template<int row>
      Scalar d1_impl (VariationalArg<Scalar,dim,TestVars::template Components<row>::m> const& arg) const
      {
        return kappa*sp(du,arg.derivative) + h.mu*u*arg.value + dF*arg.value;
      }

      template<int row, int col>
      Scalar d2_impl (VariationalArg<Scalar,dim,TestVars::template Components<row>::m> const &arg1,
                    VariationalArg<Scalar,dim,AnsatzVars::template Components<row>::m> const &arg2) const
      {
        return ( dkappa * arg2.value ) * sp(du,arg1.derivative)
            + kappa * sp(arg2.derivative[0],arg1.derivative[0]) + h.mu*arg2.value*arg1.value + ddF*arg1.value*arg2.value;
      }

    private:
      BoundaryControlledNonlinearModelPDE const& h;
      Scalar kappa,F,dF,ddF;
      typename AnsatzVars::VariableSet const& data;
      Dune::FieldVector<Scalar,AnsatzVars::template Components<uIdx>::m> u, dkappa;
      Dune::FieldMatrix<Scalar,AnsatzVars::template Components<uIdx>::m,dim> du;
      LinAlg::EuclideanScalarProduct sp;
    };

    class BoundaryCache : public CacheBase<BoundaryControlledNonlinearModelPDE,BoundaryCache>
    {
    public:
      BoundaryCache(BoundaryControlledNonlinearModelPDE<RType,AnsatzVars,DescriptionControl> const& h_,
                    typename AnsatzVars::VariableSet const& vars_,
                    int flags=7):
        data(vars_), penalty(1e9), u(0.), uDirichletBoundaryValue(0.),h(h_)
      {}

      using FaceIterator = typename AnsatzVars::Grid::LeafIntersectionIterator;

      void moveTo(FaceIterator const& face)
      {
        e = &face;
      }

      template <class Position, class Evaluators>
      void evaluateAt(Position const& x, Evaluators const& evaluators)
      {
        u = boost::fusion::at_c<uIdx>(data.data).value(boost::fusion::at_c<uSpaceIdx>(evaluators));
        auto glob_ = (*e)->geometry().global(x);
        const auto & cell = (*e)->inside();
        const auto loc = cell.geometry().local(glob_);

        f = ::boost::fusion::at_c<0>(h.control_.data).value(cell, loc);
      }

      Scalar d0() const
      {
        std::cout<<"THIS SHOULD NOT BE CALLED"<<std::endl;
        return 0;
      }

      template<int row>
      Scalar d1_impl (VariationalArg<Scalar,dim> const& arg) const
      {
        return -f*arg.value;
      }

      template<int row, int col>
      Scalar d2_impl (VariationalArg<Scalar,dim> const &arg1,
                      VariationalArg<Scalar,dim> const &arg2) const
      {
        return 0;
      }

    private:

      BoundaryControlledNonlinearModelPDE h;
      FaceIterator const* e;
      typename AnsatzVars::VariableSet const& data;
      Scalar penalty;
      Dune::FieldVector<Scalar,AnsatzVars::template Components<uIdx>::m> u, uDirichletBoundaryValue,f;
    };

    /**
     * @brief Construct linear pde
     * @param c (c*y*y + d) Heat conduction tensor (smaller d, larger c) -> harder problem
     * @param d see above
     * @param e nonlinear instability factor
     * @param mu linear instability factor  mu > 0 -> more stable , mu < 0 -> unstable
     * @param control_ source term of the pde
     */
    explicit BoundaryControlledNonlinearModelPDE(Scalar c_, Scalar d_, Scalar e_, Scalar mu_, typename DescriptionControl::VariableSet control): c(c_),d(d_),e(e_),mu(mu_),control_(control) {}
    template <int row>
    struct D1: public FunctionalBase<WeakFormulation>::D1<row>
    {
      static bool const present   = true;
      static bool const constant  = false;

    };

    template <int row, int col>
    struct D2: public FunctionalBase<WeakFormulation>::D2<row,col>
    {
      static bool const present = true;
      static bool const symmetric = false;
      static bool const lumped = false;
    };

    template <class Cell>
    int integrationOrder(Cell const& /* cell */, int shapeFunctionOrder, bool boundary) const
    {
      if (boundary)
        return 2*shapeFunctionOrder;
      else
      {
        int stiffnessMatrixIntegrationOrder = 2*(shapeFunctionOrder-1) + 2*shapeFunctionOrder;
        int sourceTermIntegrationOrder = shapeFunctionOrder;        // as rhs f is constant, i.e. of order 0

        return std::max(stiffnessMatrixIntegrationOrder,sourceTermIntegrationOrder);
      }
    }
    Scalar c, d, e, mu;
    typename DescriptionControl::VariableSet control_;

  };

}

