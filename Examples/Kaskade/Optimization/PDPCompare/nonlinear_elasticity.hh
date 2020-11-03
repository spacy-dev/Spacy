/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/Numerik/numsoft/kaskade7/                        */
/*                                                                           */
/*  Copyright (C) 2002-2009 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#include <memory>
#include <type_traits>

#include "fem/variables.hh"
#include "utilities/linalg/scalarproducts.hh"
#include "trackingFunctional.hh"
#include "functional.hh"
// #include "fem/diffops/trackingTypeCostFunctional.hh"
// #include "fem/diffops/antonsNonlinearTestProblems.hh"


/// \cond
using namespace Kaskade;


/****************************************************************************************/
/* Boundary */
/****************************************************************************************/

// Reference Deformation has to be template ??
template <class Integrand, class Reference, int stateId, int controlId, int adjointId, class RType, class AnsatzVars_, class TestVars_=AnsatzVars_, class OriginVars_=AnsatzVars_, bool lump=false >
class LagrangeFunctional
{
public:
    typedef RType  Scalar;
    typedef AnsatzVars_ AnsatzVars;
    typedef TestVars_ TestVars;
    typedef OriginVars_ OriginVars;
    static constexpr int dim = AnsatzVars::Grid::dimension;
    static ProblemType const type = std::is_same<AnsatzVars,TestVars>::value ? VariationalFunctional : WeakFormulation;

    typedef typename AnsatzVars::Grid Grid;
    typedef typename Grid::template Codim<0>::Entity Cell;

    static int const yIdx = stateId;
    static int const uIdx = controlId;
    static int const pIdx = adjointId;

    static int const ySIdx = boost::fusion::result_of::value_at_c<typename AnsatzVars::Variables,yIdx>::type::spaceIndex;
    static int const uSIdx = boost::fusion::result_of::value_at_c<typename AnsatzVars::Variables,uIdx>::type::spaceIndex;
    static int const pSIdx = boost::fusion::result_of::value_at_c<typename AnsatzVars::Variables,pIdx>::type::spaceIndex;

    class DomainCache : public CacheBase<LagrangeFunctional,DomainCache>
    {
    public:
        DomainCache(LagrangeFunctional const& f_,
                    typename AnsatzVars::VariableSet const& vars_,int flags=7):
            f(f_), vars(vars_), c(f.integrand_), J(f.J), id_(0.0), ref_(f_.ref_)
        {
           // std::cout << std::endl;
          //  std::cout << " Domain Cache Lambda: " << lambda << std::endl;
          //  std::cout << std::endl;
            id_[0][0] = id_[1][1] = id_[2][2] = 1.0;
        }

        template <class Position, class Evaluators>
        void evaluateAt(Position const& x, Evaluators const& evaluators)
        {

            // check what is necessary Do no evaluate twice!!
            using namespace boost::fusion;
            y = at_c<yIdx>(vars.data).value(at_c<ySIdx>(evaluators));
            //             if ( y*y > 10 ) std::cout << "Y: " << y*y << std::endl;
            u = at_c<uIdx>(vars.data).value(at_c<uSIdx>(evaluators));
            // check what is really necessary
            p.value = at_c<pIdx>(vars.data).value(at_c<pSIdx>(evaluators));

            p.derivative = at_c<pIdx>(vars.data).derivative(at_c<pSIdx>(evaluators));

            c.evaluateAt(vars, evaluators);

            J.evaluateAt(y,u,evaluators);

             y_z = boost::fusion::at_c<yIdx>(ref_.data).value(boost::fusion::at_c<ySIdx>(evaluators));
        }


        Scalar d0() const
        {
            return J.d0(stateId);
        }

        template<int row>
        Scalar d1_impl (VariationalArg<Scalar,dim, TestVars::template Components<row>::m> const& arg) const
        {
            // merge control into c   //Change postion of p and arg


            if(row == yIdx) return sp(y_z,arg.value);



            return 0.0;
        }

        //Todo change position arg1arg2p
        template<int row, int col>
        Scalar d2_impl (VariationalArg<Scalar,dim,TestVars::template Components<row>::m> const &arg1, VariationalArg<Scalar,dim,TestVars::template Components<col>::m> const &arg2) const
        {
            if(row == yIdx && col == yIdx)
            {
               if( c.template d3<yIdx,yIdx,yIdx>(arg1,arg2,p) != 0.0)
                   throw;

               return J.template d2<yIdx,yIdx>(arg1,arg2);
            }

            // check
            if(row == yIdx && col == pIdx) return c.template d2<yIdx,yIdx>(arg2,arg1);
            if(row == pIdx && col == yIdx) return c.template d2<yIdx,yIdx>(arg1,arg2);
            return 0.0;
        }

    private:
        LagrangeFunctional const& f;
        typename AnsatzVars::VariableSet const& vars;
        Dune::FieldVector<Scalar,AnsatzVars::template Components<yIdx>::m> y, y_z;
        Dune::FieldVector<Scalar,AnsatzVars::template Components<uIdx>::m> u;
        NonlinearElasticityFunctionalControl<Integrand,  AnsatzVars,yIdx> c;
        TrackingFunctional<Reference,  typename AnsatzVars::VariableSet,yIdx,uIdx> J;
        // Why not a FieldVector ??
        VariationalArg<Scalar,dim,AnsatzVars::template Components<pIdx>::m> p;
        LinAlg::EuclideanScalarProduct sp;
        double lambda = 0.0;
        Dune::FieldMatrix<Scalar,AnsatzVars::template Components<yIdx>::m,dim> id_;
        Reference const& ref_;
      //  Dune::FieldMatrix<Scalar,AnsatzVars::template Components<yIdx>::m,dim> R;

    };

    class BoundaryCache : public CacheBase<LagrangeFunctional,BoundaryCache>
    {
    public:
        using FaceIterator = typename AnsatzVars::Grid::LeafIntersectionIterator;
        using Vector = Dune::FieldVector<Scalar, dim>;
        BoundaryCache(LagrangeFunctional const& f,
                      typename AnsatzVars::VariableSet const& vars_, int flags=7):
            vars(vars_), dirichlet_penalty(f.gamma),  J(f.J),   up {0, 0, 1},  side {1.0,  0, 0}
        {

        }

        void moveTo(FaceIterator const& face)
        {

            const auto scalarProdUp = face->centerUnitOuterNormal() *up;
            const auto scalarProdSide = face->centerUnitOuterNormal() *side;
            // change order for more efficiency
            //
            if (scalarProdUp > 0.5)
            {
                gamma = 0.0;
                localPenalty_ = 0.0;
                beta  = 1.0;
            }
            else if (scalarProdUp < -0.5)
            {
                gamma = 0.0;
                beta  = 0.0;
            }

            else if (scalarProdSide < -0.5)
            {
                gamma = dirichlet_penalty;
                beta  = 0.0;
                localPenalty_ = 0.0;
            }


            else if(scalarProdSide > 0.5)
            {
                gamma = dirichlet_penalty;
                beta  = 0.0;
                localPenalty_ = 0.0;
            }


            else
            {
                gamma = 0.0;
                beta = 0.0;
                localPenalty_ = 0.0;
            }
        }


        template <class Evaluators>
        void evaluateAt(Dune::FieldVector<typename AnsatzVars::Grid::ctype,AnsatzVars::Grid::dimension-1> const& x, Evaluators const& evaluators)
        {
            using namespace boost::fusion;
            y = at_c<yIdx>(vars.data).value(at_c<ySIdx>(evaluators));
            u = at_c<uIdx>(vars.data).value(at_c<uSIdx>(evaluators));
            p = at_c<pIdx>(vars.data).value(at_c<pSIdx>(evaluators));

            J.evaluateAt(y,u,evaluators);
        }

        Scalar d0() const
        {
            return beta*J.d0(controlId);
        }

        template<int row>
        Scalar d1_impl (VariationalArg<Scalar,dim,TestVars::template Components<row>::m> const& arg) const
        {
                      // c' = I''
            if(row == yIdx) return (gamma*arg.value) * p ;
            if(row == uIdx) return  beta *( J.template d1<uIdx>(arg) -p*arg.value);
            if(row == pIdx) return gamma*y*arg.value - (beta* u*arg.value);

            return 0.0;
        }

        template<int row, int col>
        Scalar d2_impl (VariationalArg<Scalar,dim,TestVars::template Components<row>::m> const &arg1, VariationalArg<Scalar,dim,AnsatzVars::template Components<col>::m> const &arg2) const
        {


            if(row == yIdx && col == yIdx)
            {

               return 0.0;
            }


            if (row == pIdx && col == pIdx) return 0.0;
            if (row==yIdx && col==pIdx)     return gamma*(arg1.value * arg2.value);

            if(row==pIdx && col==yIdx)     return gamma*(arg1.value * arg2.value);
            if(row == uIdx && col == uIdx) return beta* J.template d2<uIdx,uIdx>(arg1,arg2);

            if(row==pIdx && col==uIdx)     return - beta * arg1.value * arg2.value;
            if(row==uIdx && col==pIdx)     return - beta * arg2.value * arg1.value;

            return 0.0;
        }

    private:
        typename AnsatzVars::VariableSet const& vars;
        //         Dune::FieldVector<Scalar,AnsatzVars::template Components<yIdx>::m> u;
        Dune::FieldVector<Scalar,AnsatzVars::template Components<uIdx>::m> u;
        Dune::FieldVector<Scalar,AnsatzVars::template Components<yIdx>::m> y;
        Dune::FieldVector<Scalar,AnsatzVars::template Components<pIdx>::m> p;

        // use Reference
        TrackingFunctional<Reference,  typename AnsatzVars::VariableSet,  yIdx,uIdx> J;

        const Vector up, side;

        FaceIterator const* e;
        mutable double gamma, beta;
        const double dirichlet_penalty;
        double localPenalty_;


    };

    // use std::move by copy übergabe !!!!!  // alpha J and beta ??
    explicit LagrangeFunctional(Scalar regularization, Reference const& ref,  Integrand integrand) :
        gamma(1e9), alpha(0.), J(regularization,ref), integrand_(integrand), ref_(ref)
    {

    }


    template <class T>
    bool inDomain(T const&) const
    {
        return true;
    }

    template <int row>
    struct D1 : public FunctionalBase<WeakFormulation>::D1<row>
    {
        static bool const present = true;
        static bool const constant = false;
    };

    template <int row, int col>
    struct D2 : public FunctionalBase<WeakFormulation>::D2<row,col>
    {
        static bool const present = !( ( row == pIdx && col == pIdx ) ||
                                       ( row == yIdx && col == uIdx ) ||
                                       ( row == uIdx && col == yIdx ) );
       static bool const symmetric = row==col;
//       static bool const symmetric = (row==col || ( row == pIdx && col == yIdx ) || ( row == yIdx && col == pIdx ));
   ///    static bool const lumped =  (row==uIdx && col==uIdx && role==RoleOfFunctional::PRECONDITIONING );
       static bool const lumped = false;
    };

    template <class Cell>
    int integrationOrder(Cell const& /* cell */, int shapeFunctionOrder, bool  boundary ) const
    {
        if( boundary ) return 2*shapeFunctionOrder;
        return 4*shapeFunctionOrder - 2;
    }


    // Dirichlet and regularization control
    Scalar gamma, alpha;
    TrackingFunctional<Reference, typename AnsatzVars::VariableSet, yIdx,uIdx> J;


    // elasticity
    Integrand integrand_;
    Reference const& ref_;
};




template<typename Matrix>
void writeMatlab(const Matrix & matrix)
{
    std::ofstream myfile;
    myfile.open ("matrix.dat");
    for (int i = 0; i < matrix.ridx.size(); i++)
    {
        myfile << matrix.ridx[i]+1;
        myfile << " ";
        myfile << matrix.cidx[i]+1;
        myfile << " ";
        myfile << matrix.data[i];
        myfile << std::endl;
    }
    myfile.close();
}


/// \endcond

