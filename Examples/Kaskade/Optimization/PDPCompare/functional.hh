#pragma once

#include "../../traits.hh"
#include "../../boundary_caches.hh"
#include "domain_caches.hh"

#include <fem/functional_aux.hh>
#include <algorithm>
// Domain Cache !!!!!!!
namespace Kaskade
{
template<class Integrand, class VarSet,  class ControlVarSet,  int state_index = 0, int control_index = 0>
class NonlinearElasticityFunctionalDistributed: public FunctionalBase<VariationalFunctional>
{
public:
    using Scalar = double;
    using OriginVars = VarSet;
    using AnsatzVars = VarSet;
    using TestVars = VarSet;
    using Grid = typename AnsatzVars::Grid;
    static const int dim = AnsatzVars::Grid::dimension;
    using Vector = Dune::FieldVector<Scalar, dim>;

    double normalComplianceParameter_;
    const ControlVarSet & controlVar_;
    const double obstacle_;

    static int constexpr controlIndex_ = 0;
    static int constexpr controlSpaceIndex_ = boost::fusion::result_of::value_at_c<ControlVarSet, controlIndex_>::type::spaceIndex;

    static int constexpr u_Idx = 0;
    static int constexpr u_Space_Idx = boost::fusion::result_of::value_at_c<
            typename AnsatzVars::Variables, u_Idx>::type::spaceIndex;


    void paramUpdate(double update)
    {
        normalComplianceParameter_ = update;
    }

    template<int row>
    using D1 = NonConstD1<row>;

    template<int row, int col>
    using D2 = D2PresentAndNonSymmetric<row,col>;

    using DomainCache = FunctionalWithConstantSource<NonlinearElasticityFunctionalDistributed, Integrand, state_index>;

    class BoundaryCache: public CacheBase<NonlinearElasticityFunctionalDistributed,
            BoundaryCache>
    {
    public:
        using FaceIterator = typename AnsatzVars::Grid::LeafIntersectionIterator;

        BoundaryCache(NonlinearElasticityFunctionalDistributed const& f_,
                      typename AnsatzVars::VariableSet const& vars_, int flags = 7) :
            vars(vars_), normalComplianceParameter_(
                             f_.normalComplianceParameter_), beta(0.), up
        { 0, 0, 1 }, side
        { 1, 0, 0 }, localPenalty_(0.0), alpha(0.0), controlVarset_(f_.controlVar_), obstacle_(f_.obstacle_)
        {
        }
        void moveTo(FaceIterator const& face)
        {
            e = &face;

            const auto scalarProdUp = face->centerUnitOuterNormal() *up;
            const auto scalarProdSide = face->centerUnitOuterNormal() *side;

            if (scalarProdUp > 0.5)
            {
                alpha = 0.0;
                localPenalty_ = 0.0;

            }
            else if (scalarProdUp < -0.5)
            {
                alpha = 0.0;
                localPenalty_ = normalComplianceParameter_;
            }

            else if(scalarProdSide < -0.5)
            {
                alpha = dirichlet_penalty;
                localPenalty_ = 0.0;
            }

            else if(scalarProdSide > 0.5)
            {
                alpha = dirichlet_penalty;
                localPenalty_ = 0.0;
            }

            else
            {
                alpha = 0.0;
                localPenalty_ = normalComplianceParameter_;
            }
        }

        template<class Evaluators>
        void evaluateAt(
                Dune::FieldVector<typename AnsatzVars::Grid::ctype, dim - 1> const& x,
                Evaluators const& evaluators)
        {
            using namespace boost::fusion;
            loc_ = x;
            glob_ = (*e)->geometry().global(x);
            const auto & cell = (*e)->inside();

            // Can use evaluators ??
            const auto loc = cell.geometry().local(glob_);

            control_= at_c<0>(controlVarset_.data).value(cell, loc);
            u = at_c<u_Idx>(vars.data).value(at_c<u_Space_Idx>(evaluators));

            violation_ = std::max(obstacle_ - u[2], 0.0);

            if(violation_ == 0.0)
                contact_ = 0.0;

            else
                contact_ = 1.0;

        }

        //TODO
        Scalar d0() const
        {
            return 0.5 * alpha * (u * u) - control_ * (u  +glob_)
                    + localPenalty_ * (violation_) * violation_ * (violation_) *(1.0 /3.0);
        }

        template<int row>
        Scalar d1_impl(VariationalArg<Scalar, dim, dim> const& arg) const
        {
            return alpha * (u * arg.value) - control_ * arg.value
                    - localPenalty_ * violation_ * violation_ * arg.value[2];
        }

        template<int row, int col>
        Scalar d2_impl(VariationalArg<Scalar, dim, dim> const &arg1,
                       VariationalArg<Scalar, dim, dim> const &arg2) const
        {

            return alpha * (arg1.value * arg2.value)
                    + localPenalty_ * 2.0* violation_ * arg2.value[2] * arg1.value[2];

        }

        static constexpr Scalar dirichlet_penalty = 1e9;

        void updatePenaltyParameter(double update)
        {
            normalComplianceParameter_ = update;
        }


    private:
        typename AnsatzVars::VariableSet const& vars;
        double normalComplianceParameter_;


        Vector beta;
        const Vector up, side;
        Vector control_;
        double localPenalty_;
        Scalar alpha;
        const ControlVarSet & controlVarset_;
        const double obstacle_;
        Vector u = {0.0,0.0,0.0};
        Dune::FieldVector<typename AnsatzVars::Grid::ctype, dim - 1> loc_ = {0.0,0.0};
        Dune::FieldVector<typename AnsatzVars::Grid::ctype, dim > glob_ = {0.0,0.0,0.0};

        double violation_ = 0.0;
        double contact_ = 0.0;

        FaceIterator const* e = nullptr;
    };

    Integrand f_;

    explicit NonlinearElasticityFunctionalDistributed(Integrand f,
                                                      const double normalComplianceParameter, const ControlVarSet & controlVar, double obstacle) :
        f_(std::move(f)), normalComplianceParameter_(
                              normalComplianceParameter), controlVar_(controlVar), obstacle_(obstacle)
    {
    }

    template<class Cell>
    int integrationOrder(Cell const& /*cell*/, int shapeFunctionOrder,
                         bool boundary) const
    {
        if (boundary)
            return 2 * shapeFunctionOrder;
        else
        {
            int stiffnessMatrixIntegrationOrder = 2 * (shapeFunctionOrder - 1);
            int sourceTermIntegrationOrder = shapeFunctionOrder;
            return std::max(stiffnessMatrixIntegrationOrder,
                            sourceTermIntegrationOrder);
        }
    }
};


template<class Integrand, class VarSet,  class ControlVarSet,  int state_index = 0, int control_index = 0>
class NonlinearElasticityRegularization: public FunctionalBase<VariationalFunctional>
{
public:
    using Scalar = double;
    using OriginVars = VarSet;
    using AnsatzVars = VarSet;
    using TestVars = VarSet;
    using Grid = typename AnsatzVars::Grid;
    static const int dim = AnsatzVars::Grid::dimension;
    using Vector = Dune::FieldVector<Scalar, dim>;

    double normalComplianceParameter_;
    const ControlVarSet & controlVar_;
    const double obstacle_;

    static int constexpr controlIndex_ = 0;
    static int constexpr controlSpaceIndex_ = boost::fusion::result_of::value_at_c<ControlVarSet, controlIndex_>::type::spaceIndex;

    static int constexpr u_Idx = 0;
    static int constexpr u_Space_Idx = boost::fusion::result_of::value_at_c<
            typename AnsatzVars::Variables, u_Idx>::type::spaceIndex;


    void paramUpdate(double update)
    {
        normalComplianceParameter_ = update;
    }

    template<int row>
    using D1 = NonConstD1<row>;

    template<int row, int col>
    using D2 = D2PresentAndNonSymmetric<row,col>;

    class DomainCache : public CacheBase<NonlinearElasticityRegularization,DomainCache>
    {
    public:
        DomainCache(NonlinearElasticityRegularization const& f,
                    typename AnsatzVars::VariableSet const& vars_,int flags=7):
             c(f.f_), cLin(f.f_), id_(0.0), energyRegularization_(f.energyRegularization_), vars(vars_)
        {id_[0][0] = id_[1][1] = id_[2][2] = 1.0;}


        template <class Position, class Evaluators>
        void evaluateAt(const Position& x, Evaluators const& evaluators)
        {
            using boost::fusion::at_c;
               c.update(id_ +at_c<u_Idx>(vars.data).derivative(at_c<u_Space_Idx>(evaluators)));
               cLin.update(id_);
        //    y_ = at_c<state_index>(x_.data).value(at_c<state_space_index>(evaluators));
         //   f_.update(id_+ at_c<state_index>(x_.data).derivative(at_c<state_space_index>(evaluators)) );
        }

        Scalar d0() const
        {
            return c.d0();
        }

        template<int row>
        Scalar d1_impl (VariationalArg<Scalar,dim, TestVars::template Components<row>::m> const& arg) const
        {
            // merge control into c   //Change postion of p and arg
            return c.d1( arg.derivative );
        }

        template<int row, int col>
        Scalar d2_impl (VariationalArg<Scalar,dim,TestVars::template Components<row>::m> const &arg1, VariationalArg<Scalar,dim,TestVars::template Components<col>::m> const &arg2) const
        {
           return c.d2( arg1.derivative, arg2.derivative )+ energyRegularization_ * cLin.d2( arg1.derivative, arg2.derivative ) ;
        }

    private:

        mutable Integrand c;
        mutable Integrand cLin;
        Dune::FieldMatrix<Scalar,3,3> id_;

        double energyRegularization_ = 0.0;
        typename AnsatzVars::VariableSet const& vars;
        LinAlg::EuclideanScalarProduct sp;
    };

    class BoundaryCache: public CacheBase<NonlinearElasticityRegularization,
            BoundaryCache>
    {
    public:
        using FaceIterator = typename AnsatzVars::Grid::LeafIntersectionIterator;

        BoundaryCache(NonlinearElasticityRegularization const& f_,
                      typename AnsatzVars::VariableSet const& vars_, int flags = 7) :
            vars(vars_), normalComplianceParameter_(
                             f_.normalComplianceParameter_), beta(0.), up
        { 0, 0, 1 }, side
        { 1, 0, 0 }, localPenalty_(0.0), alpha(0.0), controlVarset_(f_.controlVar_), obstacle_(f_.obstacle_), energyRegularization_(f_.energyRegularization_)
        {
        }
        void moveTo(FaceIterator const& face)
        {
            e = &face;

            const auto scalarProdUp = face->centerUnitOuterNormal() *up;
            const auto scalarProdSide = face->centerUnitOuterNormal() *side;

            if (scalarProdUp > 0.5)
            {
                alpha = 0.0;
                localPenalty_ = 0.0;
            }

            else if (scalarProdUp < -0.5)
            {
                alpha = 0.0;
                localPenalty_ = normalComplianceParameter_;
            }


            else if(scalarProdSide < -0.5)
            {
                alpha = dirichlet_penalty;
                localPenalty_ = 0.0;
            }

            else if(scalarProdSide > 0.5)
            {
                alpha = dirichlet_penalty;
                localPenalty_ = 0.0;
            }


            else
            {
                alpha = 0.0;
                localPenalty_ = normalComplianceParameter_;
            }
        }

        template<class Evaluators>
        void evaluateAt(
                Dune::FieldVector<typename AnsatzVars::Grid::ctype, dim - 1> const& x,
                Evaluators const& evaluators)
        {
            using namespace boost::fusion;
            loc_ = x;
            glob_ = (*e)->geometry().global(x);
            const auto & cell = (*e)->inside();

            // Can use evaluators ??
            const auto loc = cell.geometry().local(glob_);

            control_= at_c<0>(controlVarset_.data).value(cell, loc);
            u = at_c<u_Idx>(vars.data).value(at_c<u_Space_Idx>(evaluators));

            violation_ = std::max(obstacle_ - u[2], 0.0);

            if(violation_ == 0.0)
                contact_ = 0.0;

            else
                contact_ = 1.0;


        }

        Scalar d0() const
        {
                return 0.5 * alpha * (u * u) - control_ * (u  +glob_)
            + localPenalty_ * (violation_)
                    * (violation_) * violation_ *(1.0 / 3.0);
        }

        template<int row>
        Scalar d1_impl(VariationalArg<Scalar, dim, dim> const& arg) const
        {
                return alpha * (u * arg.value) - control_ * arg.value
            - localPenalty_  * (violation_) *violation_ * arg.value[2];
        }

        template<int row, int col>
        Scalar d2_impl(VariationalArg<Scalar, dim, dim> const &arg1,
                        VariationalArg<Scalar, dim, dim> const &arg2) const
        {

                return alpha * (arg1.value * arg2.value)
            + 2.0 * violation_  * localPenalty_ * arg2.value[2] * arg1.value[2];
        }

        static constexpr Scalar dirichlet_penalty = 1e9;
    private:
        typename AnsatzVars::VariableSet const& vars;
        const double normalComplianceParameter_;


        Vector beta;
        const Vector up, side;
        Vector control_;
        double localPenalty_;
        Scalar alpha;
        const ControlVarSet & controlVarset_;
        const double obstacle_;
        Vector u = {0.0,0.0,0.0};
        Dune::FieldVector<typename AnsatzVars::Grid::ctype, dim - 1> loc_ = {0.0,0.0};
        Dune::FieldVector<typename AnsatzVars::Grid::ctype, dim > glob_ = {0.0,0.0,0.0};

        double violation_ = 0.0;

        FaceIterator const* e = nullptr;
        double energyRegularization_ = 0.0;
        LinAlg::EuclideanScalarProduct sp;

        double contact_ = 0.0;

    };

    Integrand f_;
    double energyRegularization_ = 0.0;

    explicit NonlinearElasticityRegularization(Integrand f,
                                                      const double normalComplianceParameter, const ControlVarSet & controlVar, double obstacle, double energyRegularization) :
        f_(std::move(f)), normalComplianceParameter_(
                              normalComplianceParameter), controlVar_(controlVar), obstacle_(obstacle), energyRegularization_(energyRegularization)
    {
    }

    template<class Cell>
    int integrationOrder(Cell const& /*cell*/, int shapeFunctionOrder,
                         bool boundary) const
    {
        if (boundary)
            return 2 * shapeFunctionOrder;
        else
        {
            int stiffnessMatrixIntegrationOrder = 2 * (shapeFunctionOrder - 1);
            int sourceTermIntegrationOrder = shapeFunctionOrder;
            return std::max(stiffnessMatrixIntegrationOrder,
                            sourceTermIntegrationOrder);
        }
    }
};








//
template<class Integrand, class VarSet, int state_index = 0>
class NonlinearElasticityFunctionalControl: public FunctionalBase<
        VariationalFunctional>
{
public:
    using Scalar = double;
    using OriginVars = VarSet;
    using AnsatzVars = VarSet;
    using TestVars = VarSet;
    static const int dim = AnsatzVars::Grid::dimension;
    using Vector = Dune::FieldVector<Scalar, dim>;
    static int constexpr y_Idx = 0;
    static int constexpr y_Space_Idx = boost::fusion::result_of::value_at_c<
            typename AnsatzVars::Variables, y_Idx>::type::spaceIndex;

    template<int row>
    using D1 = NonConstD1<row>;

    template<int row, int col>
    using D2 = D2PresentAndNonSymmetric<row,col>;

    template<class Evaluators>
    void evaluateAt(const typename AnsatzVars::VariableSet& y,
                    Evaluators const& evaluators)
    {
        using namespace boost::fusion;
        //       y_ = at_c<state_index>(y.data).value(at_c<y_Idx>(evaluators));
        f_.update(id_ );
    }

    // Problem: Boundary Conditions cannot be enforced in the functional

    Scalar d0() const
    {
        return f_();
    }

    template<int row>
    Scalar d1(
            VariationalArg<Scalar, dim, TestVars::template Components<row>::m> const &arg) const
    {
        return f_.d1(arg.derivative);
    }

    template<int row, int col>
    Scalar d2(
            VariationalArg<Scalar, dim, TestVars::template Components<row>::m> const &arg1,
            VariationalArg<Scalar, dim, TestVars::template Components<col>::m> const &arg2) const
    {
        return f_.d2(arg1.derivative, arg2.derivative);
    }
    // Check state_index and Ids in general
    template<int id1, int id2, int id3>
    Scalar d3(
            VariationalArg<Scalar, dim, TestVars::template Components<id1>::m> const &arg1,
            VariationalArg<Scalar, dim, TestVars::template Components<id2>::m> const &arg2,
            VariationalArg<Scalar, dim, TestVars::template Components<id3>::m> const &arg3) const
    {
        if (id1 != state_index || id2 != state_index || id3 != state_index)
            return 0.0;
        // Not sure about this
        return f_.d3(arg1.derivative, arg2.derivative, arg3.derivative);
    }
    // use std::move ??
    explicit NonlinearElasticityFunctionalControl(Integrand f) :id_(0.0),
        f_(std::move(f))
    {
        id_[0][0] = id_[1][1] = id_[2][2] = 1.0;
    }
private:
    Dune::FieldMatrix<Scalar,AnsatzVars::template Components<state_index>::m,dim> id_;
    Integrand f_;

};





template <class Integrand, int stateId, int controlId, int adjointId, class RType, class AnsatzVars_, class TestVars_=AnsatzVars_, class OriginVars_=AnsatzVars_,  bool lump=false >
class VolumePenaltyFunctional
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

    Integrand integrand_;

    class DomainCache : public CacheBase<VolumePenaltyFunctional,DomainCache>
    {
    public:
        DomainCache(VolumePenaltyFunctional const& f_,
                    typename AnsatzVars::VariableSet const& vars_,int flags=7):
            f(f_), vars(vars_), c(f.integrand_)
        {}

        template <class Position, class Evaluators>
        void evaluateAt(Position const& x, Evaluators const& evaluators)
        {

            // check what is necessary Do no evaluate twice!!
            using namespace boost::fusion;
            c.evaluateAt(vars, evaluators);




        }

        Scalar d0() const
        {
            return c.d0();
        }

        template<int row>
        Scalar d1_impl (VariationalArg<Scalar,dim, TestVars::template Components<row>::m> const& arg) const
        {
            // merge control into c   //Change postion of p and arg


            return 0;
        }

        template<int row, int col>
        Scalar d2_impl (VariationalArg<Scalar,dim,TestVars::template Components<row>::m> const &arg1, VariationalArg<Scalar,dim,TestVars::template Components<col>::m> const &arg2) const
        {

            return 0;
        }

    private:
        VolumePenaltyFunctional const& f;
        typename AnsatzVars::VariableSet const& vars;
        NonlinearElasticityFunctionalControl<Integrand,  AnsatzVars,yIdx> c;

    };

    class BoundaryCache : public CacheBase<VolumePenaltyFunctional,BoundaryCache>
    {
    public:

        BoundaryCache(VolumePenaltyFunctional const& f,
                      typename AnsatzVars::VariableSet const& vars_, int flags=7)

        {}




        template <class Evaluators>
        void evaluateAt(Dune::FieldVector<typename AnsatzVars::Grid::ctype,AnsatzVars::Grid::dimension-1> const& x, Evaluators const& evaluators)
        {
            using namespace boost::fusion;

        }

        Scalar d0() const
        {
            return 0.0;
        }

        template<int row>
        Scalar d1_impl (VariationalArg<Scalar,dim,TestVars::template Components<row>::m> const& arg) const
        {

            return 0.0;
        }

        template<int row, int col>
        Scalar d2_impl (VariationalArg<Scalar,dim,TestVars::template Components<row>::m> const &arg1, VariationalArg<Scalar,dim,AnsatzVars::template Components<col>::m> const &arg2) const
        {
            return 0.0;
        }

    private:


    };

    // use std::move by copy übergabe !!!!!  // alpha J and beta ??
    explicit VolumePenaltyFunctional (  Integrand integrand) :
        integrand_(integrand)
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
        static bool const present = !(
                    ( row == pIdx && col == pIdx ) ||
                    ( row == yIdx && col == uIdx ) ||
                    ( row == uIdx && col == yIdx ) );
        static bool const symmetric = row==col;
    };

    template <class Cell>
    int integrationOrder(Cell const& /* cell */, int shapeFunctionOrder, bool  boundary ) const
    {
        if( boundary ) return 2*shapeFunctionOrder;
        return 4*shapeFunctionOrder - 2;
    }



};
}
