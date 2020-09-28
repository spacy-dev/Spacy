#pragma once

#include <memory>
#include <utility>
#include <chrono>

#include "fem/assemble.hh"
#include "fem/istlinterface.hh"
#include "linalg/triplet.hh"
#include "linalg/direct.hh"

#include <linalg/lapack/lapacke.h>
#include "mg/multigrid.hh"
#include "mg/additiveMultigrid.hh"

#include <Spacy/C1Operator.h>
#include <Spacy/Vector.h>
#include <Spacy/VectorSpace.h>
#include <Spacy/ZeroVectorCreator.h>
#include <Spacy/Util/Mixins/Eps.h>
#include <Spacy/Util/Mixins/NumberOfThreads.h>
#include <Spacy/Util/Base/FunctionalBase.h>
#include <Spacy/Spaces/ProductSpace/Operator_V2.h>
#include <Spacy/Util/Mixins/Get.h>

#include <Spacy/InducedScalarProduct.h>

#include "DirectSolver.h"
#include "VectorSpace.h"
#include "Vector.h"
#include "OperatorSpace.h"
#include "LinearOperator.h"


namespace Spacy
{
/** @addtogroup KaskadeGroup
     * @{
     */
namespace Kaskade
{
/**
         * @brief %Functional interface for %Kaskade 7. Models a twice differentiable functional \f$f:X\rightarrow \mathbb{R}\f$.
         * @tparam FunctionalDefinition functional definition from %Kaskade 7
         * @see ::Spacy::C2BlockFunctional
         */
template <class FunctionalDefinition, class GridManager>
class C2BlockFunctional :
        public FunctionalBase,
        public Mixin::Eps,
        public Mixin::NumberOfThreads
{
public:
    /// %Kaskade::VariableSetDescription
    using VariableSetDescription = typename FunctionalDefinition::AnsatzVars;
    /// Coefficient vector type.
    using CoefficientVector = typename VariableSetDescription::template CoefficientVectorRepresentation<>::type;
    /// boost::fusion::vector<const Space0*,const Space1*,...>
    using Spaces = typename VariableSetDescription::Spaces;
    /// %Kaskade::VariationalFunctionalAssembler
    using Assembler = ::Kaskade::VariationalFunctionalAssembler< ::Kaskade::LinearizationAt<FunctionalDefinition> >;
    /// Matrix type
    using Matrix = ::Kaskade::MatrixAsTriplet<double>;
    /// operator for the description of the second derivative
    using KaskadeOperator = ::Kaskade::MatrixRepresentedOperator<Matrix,CoefficientVector,CoefficientVector>;



    /// VariableSetDescriptions of the variables
    using VYSetDescription = Detail::ExtractDescription_t<VariableSetDescription,0>;
    using VUSetDescription = Detail::ExtractDescription_t<VariableSetDescription,1>;
    using VPSetDescription = Detail::ExtractDescription_t<VariableSetDescription,2>;

    /// Coefficient Vector type
    using CoefficientVectorY = typename VYSetDescription::template CoefficientVectorRepresentation<>::type;
    using CoefficientVectorU = typename VUSetDescription::template CoefficientVectorRepresentation<>::type;
    using CoefficientVectorP = typename VPSetDescription::template CoefficientVectorRepresentation<>::type;

    /// Kaskade Types of the Operators held for every timestep
    template<class X, class Y>
    using KaskadeOperatorXY = ::Kaskade::MatrixRepresentedOperator<Matrix,X,Y>;
    using Mytype = KaskadeOperatorXY<CoefficientVectorY,CoefficientVectorY>;
    using Myutype = KaskadeOperatorXY<CoefficientVectorU,CoefficientVectorY>;
    using Mutype = KaskadeOperatorXY<CoefficientVectorY,CoefficientVectorY>;
    using Muytype = KaskadeOperatorXY<CoefficientVectorY,CoefficientVectorU>;
    using Atype = KaskadeOperatorXY<CoefficientVectorY,CoefficientVectorP>;
    using ATtype = KaskadeOperatorXY<CoefficientVectorP,CoefficientVectorY>;
    using Btype = KaskadeOperatorXY<CoefficientVectorU,CoefficientVectorP>;
    using BTtype = KaskadeOperatorXY<CoefficientVectorP,CoefficientVectorU>;

    using PSVC = ::Spacy::ProductSpace::VectorCreator;

    using Linearization = LinearOperator<VariableSetDescription,VariableSetDescription>;

    static constexpr unsigned bpxIndex = 3;
     bool useDirectMu = true;

    /**
             * @brief Construct a twice differentiable functional \f$f: X\rightarrow \mathbb{R}\f$ from %Kaskade 7.
             * @param f operator definition from %Kaskade 7
             * @param domain domain space
             * @param rbegin first row to be considered in the definition of f
             * @param rend one after the last row to be considered in the definition of f
             * @param cbegin first column to be considered in the definition of f
             * @param cend one after the last column to be considered in the definition of f
             *
             * The optional parameters rbegin, rend, cbegin and cend can be used to define operators that correspond to parts of
             * a system of equation.
             */
    C2BlockFunctional(const FunctionalDefinition& f, const VectorSpace& domain, const GridManager & gridManager,
                      int rbegin = 0, int rend = FunctionalDefinition::AnsatzVars::noOfVariables,
                      int cbegin = 0, int cend = FunctionalDefinition::TestVars::noOfVariables)
        : FunctionalBase(domain),
          f_(f),
          spaces_( extractSpaces<VariableSetDescription>(domain) ),
          gridManager_(gridManager),
          rhs_( zero(domain.dualSpace()) ),
          rbegin_(rbegin), rend_(rend), cbegin_(cbegin), cend_(cend),
          operatorSpace_( std::make_shared<VectorSpace>(
                              LinearOperatorCreator<VariableSetDescription,VariableSetDescription>(domain,domain.dualSpace()) ,
                              [](const ::Spacy::Vector& v)
    {

        using std::begin;
        using std::end;
        const auto& m = cast_ref<Linearization>(v).get();
        auto iend = end(m);
        auto result = 0.;
        for(auto iter = begin(m); iter!=iend; ++iter)
            result += (*iter) * (*iter);
        return Real{sqrt(result)};
    } , true ) )
    {    this->domain();}


    /**
             * @brief Copy constructor.
             * @param g functional to copy from
             */
    C2BlockFunctional(const C2BlockFunctional& g)
        : FunctionalBase(g.domain()),
          NumberOfThreads(g),
          f_(g.f_), spaces_(g.spaces_),
          gridManager_(g.gridManager_),
          A_(g.A_),
          value_(g.value_),
          old_X_f_(g.old_X_f_),
          old_X_df_(g.old_X_df_),
          old_X_ddf_(g.old_X_ddf_),
          rhs_(g.rhs_),
          onlyLowerTriangle_(g.onlyLowerTriangle_),
          rbegin_(g.rbegin_), rend_(g.rend_), cbegin_(g.cbegin_), cend_(g.cend_),
          solverCreator_(g.solverCreator_),
          operatorSpace_(g.operatorSpace_),
          My_(g.My_),
          Myu_(g.Myu_),
          Mu_(g.Mu_),
          Muy_(g.Muy_),
          AOp_(g.AOp_),
          B_(g.B_),
          AOp_t_(g.AOp_t_),
          B_t_(g.B_t_),
          MuDiag_(g.MuDiag_),
          solMu_(g.solMu_),
          solA_(g.solA_),
          solAT_(g.solAT_),
          BPXATPtr_(g.BPXATPtr_),
          BPXAPtr_(g.BPXAPtr_)
    {

        MyBlock_ = [this](const ::Spacy::Vector& y) mutable
        {
            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
            CoefficientVectorY y_coeff(boost::fusion::at_c<0>(temp.data));
            CoefficientVectorY My_coeff(boost::fusion::at_c<0>(temp.data));

            ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(y,y_coeff);
            //std::cout << "  My: y_coeff.two_norm() = " << y_coeff.two_norm() << std::endl;

            My_.apply(y_coeff,My_coeff);

            const auto& domain_ps = cast_ref<PSVC>(this->domain().creator());
            const auto& Y = (cast_ref<PSVC>(domain_ps.subSpace(0).creator())).subSpace(0);

            ::Spacy::Vector result = zero(Y);
            ::Spacy::Kaskade::copyFromCoefficientVector<VYSetDescription>(My_coeff,result);
            //std::cout << "  My: get(norm(result)) = " << get(norm(result)) << std::endl;

            return result;
        };

        MyuBlock_ = [this](const ::Spacy::Vector& u) mutable
        {
            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
            CoefficientVectorU u_coeff(boost::fusion::at_c<1>(temp.data));
            CoefficientVectorY Myu_coeff(boost::fusion::at_c<0>(temp.data));

            ::Spacy::Kaskade::copyToCoefficientVector<VUSetDescription>(u,u_coeff);
            //std::cout << "  Myu: u_coeff.two_norm() = " << u_coeff.two_norm() << std::endl;

            //Im Normalenschritt ist Myu_ leer -> gib einfach zero(...) zurück?
            if(Myu_.getTriplet().nrows() == 0)
            {
                const auto& domain_ps = cast_ref<PSVC>(this->domain().creator());
                const auto& Y = (cast_ref<PSVC>(domain_ps.subSpace(0).creator())).subSpace(0);
                return zero(Y);
            }
            Myu_.apply(u_coeff,Myu_coeff);

            const auto& domain_ps = cast_ref<PSVC>(this->domain().creator());
            const auto& Y = (cast_ref<PSVC>(domain_ps.subSpace(0).creator())).subSpace(0);

            ::Spacy::Vector result = zero(Y);
            ::Spacy::Kaskade::copyFromCoefficientVector<VYSetDescription>(Myu_coeff,result);
            //std::cout << "  Myu: get(norm(result)) = " << get(norm(result)) << std::endl;

            return result;
        };

        MuBlock_ = [this](const ::Spacy::Vector& u) mutable
        {
            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
            CoefficientVectorU u_coeff(boost::fusion::at_c<1>(temp.data));
            CoefficientVectorU Mu_coeff(boost::fusion::at_c<1>(temp.data));

            ::Spacy::Kaskade::copyToCoefficientVector<VUSetDescription>(u,u_coeff);
            //std::cout << "  Mu: u_coeff.two_norm() = " << u_coeff.two_norm() << std::endl;

            Mu_.apply(u_coeff,Mu_coeff);

            const auto& domain_ps = cast_ref<PSVC>(this->domain().creator());
            const auto& U = (cast_ref<PSVC>(domain_ps.subSpace(0).creator())).subSpace(1);

            ::Spacy::Vector result = zero(U);
            ::Spacy::Kaskade::copyFromCoefficientVector<VUSetDescription>(Mu_coeff,result);
            //std::cout << "  Mu: get(norm(result)) = " << get(norm(result)) << std::endl;

            return result;
        };

        MuyBlock_ = [this](const ::Spacy::Vector& y) mutable
        {
            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
            CoefficientVectorY y_coeff(boost::fusion::at_c<0>(temp.data));
            CoefficientVectorU Muy_coeff(boost::fusion::at_c<1>(temp.data));

            ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(y,y_coeff);
            //std::cout << "  Muy: y_coeff.two_norm() = " << y_coeff.two_norm() << std::endl;

            //Im Normalenschritt ist Muy_ leer -> gib einfach zero(...) zurück?
            if(Muy_.getTriplet().nrows() == 0)
            {
                const auto& domain_ps = cast_ref<PSVC>(this->domain().creator());
                const auto& U = (cast_ref<PSVC>(domain_ps.subSpace(0).creator())).subSpace(1);
                return zero(U);
            }
            Muy_.apply(y_coeff,Muy_coeff);

            const auto& domain_ps = cast_ref<PSVC>(this->domain().creator());
            const auto& U = (cast_ref<PSVC>(domain_ps.subSpace(0).creator())).subSpace(1);

            ::Spacy::Vector result = zero(U);
            ::Spacy::Kaskade::copyFromCoefficientVector<VUSetDescription>(Muy_coeff,result);
            //std::cout << "  Muy: get(norm(result)) = " << get(norm(result)) << std::endl;

            return result;
        };

        ATBlock_ = [this](const ::Spacy::Vector& p) mutable
        {


           // std::cout << "Test0: AT " << std::endl;

           // std::cout << "Test Input: " << p(p) << std::endl;
            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
            CoefficientVectorP p_coeff(boost::fusion::at_c<2>(temp.data));
            CoefficientVectorY ATp_coeff(boost::fusion::at_c<0>(temp.data));

            ::Spacy::Kaskade::copyToCoefficientVector<VPSetDescription>(p,p_coeff);
            //std::cout << "  AT: p_coeff.two_norm() = " << p_coeff.two_norm() << std::endl;

            unsigned int pVectorSize = boost::fusion::at_c<0>(p_coeff.data).size();
            unsigned int pBlockSize =  boost::fusion::at_c<0>(p_coeff.data).operator[](0).size();

            unsigned int ATpVectorSize = boost::fusion::at_c<0>(ATp_coeff.data).size();
            unsigned int ATpBlockSize =  boost::fusion::at_c<0>(ATp_coeff.data).operator[](0).size();


            unsigned int pSize =   pVectorSize * pBlockSize;
            unsigned int ATpSize =  ATpVectorSize * ATpBlockSize;

           // std::cout << "pSize: " << pSize << std::endl;
           // std::cout << "BTpSize: " << BTpSize << std::endl;

            std::vector<double> p_(pSize,0.0);
            std::vector<double> ATp_(ATpSize,0.0);

           // auto & matrix = AOp_t_.get_non_const();

            auto & matrix = AOp_.get_non_const();
            auto counter = 0u;

            for(auto & it: boost::fusion::at_c<0>(p_coeff.data))
                for(auto & iter: it)
                {
                    p_[counter] = iter;
                    counter++;
                }

            for(size_t i=0; i< matrix.ridx.size(); ++i)
                 ATp_[matrix.ridx[i]] += p_[matrix.cidx[i]]*matrix.data[i];

            counter = 0;

            for(auto & it: boost::fusion::at_c<0>(ATp_coeff.data))
                for(auto & iter: it)
                {
                    iter = ATp_[counter] ;
                    counter++;
                }


         //   AOp_t_.apply(p_coeff,ATp_coeff);

            const auto& domain_ps = cast_ref<PSVC>(this->domain().creator());
            const auto& Y = (cast_ref<PSVC>(domain_ps.subSpace(0).creator())).subSpace(0);

            ::Spacy::Vector result = zero(Y);
            ::Spacy::Kaskade::copyFromCoefficientVector<VYSetDescription>(ATp_coeff,result);
            //std::cout << "  AT: get(norm(result)) = " << get(norm(result)) << std::endl;

            return result;
        };

        ABlock_ = [this](const ::Spacy::Vector& y) mutable
        {

           // std::cout << "Test0: A " << std::endl;
           //  std::cout << "Test Input: " << y(y) << std::endl;
            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
            CoefficientVectorY y_coeff(boost::fusion::at_c<0>(temp.data));
            CoefficientVectorP Ay_coeff(boost::fusion::at_c<2>(temp.data));

            ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(y,y_coeff);
            //std::cout << "  A: y_coeff.two_norm() = " << y_coeff.two_norm() << std::endl;
            //std::cout << "  A: get(norm(y)) = " << get(norm(y)) << std::endl;

            // y_coeff * Ay_coeff

           // AOp_.apply(y_coeff,Ay_coeff);

            unsigned int yVectorSize = boost::fusion::at_c<0>(y_coeff.data).size();
            unsigned int yBlockSize =  boost::fusion::at_c<0>(y_coeff.data).operator[](0).size();

            unsigned int AyVectorSize = boost::fusion::at_c<0>(Ay_coeff.data).size();
            unsigned int AyBlockSize =  boost::fusion::at_c<0>(Ay_coeff.data).operator[](0).size();

            unsigned int ySize =  yVectorSize * yBlockSize;
            unsigned int AySize =  AyVectorSize * AyBlockSize;

           // std::cout << "uSize: " << uSize << std::endl;
           // std::cout << "BuSize: " << BuSize << std::endl;

            std::vector<double> y_(ySize,0.0);
            std::vector<double> Ay_(AySize,0.0);

            auto & matrix = AOp_.get_non_const();

            auto counter = 0u;

            for(auto & it: boost::fusion::at_c<0>(y_coeff.data))
                for(auto & iter: it)
                {
                    y_[counter] = iter;
                    counter++;
                }

            for(size_t i=0; i< matrix.ridx.size(); ++i)
                 Ay_[matrix.ridx[i]] += y_[matrix.cidx[i]]*matrix.data[i];

            counter = 0;

            for(auto & it: boost::fusion::at_c<0>(Ay_coeff.data))
                for(auto & iter: it)
                {
                    iter = Ay_[counter] ;
                    counter++;
                }


            const auto& domain_ps = cast_ref<PSVC>(this->domain().creator());
            const auto& P = (cast_ref<PSVC>(domain_ps.subSpace(1).creator())).subSpace(0);

            ::Spacy::Vector result = zero(P);
            ::Spacy::Kaskade::copyFromCoefficientVector<VPSetDescription>(Ay_coeff,result);
            //std::cout << "  A: get(norm(result)) = " << get(norm(result)) << std::endl;

            return result;
        };

         APrimePrime = [this](const ::Spacy::Vector& y) mutable
        {

            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
            CoefficientVectorY y_coeff(boost::fusion::at_c<0>(temp.data));
            CoefficientVectorY Ay_coeff(boost::fusion::at_c<0>(temp.data));

            ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(y,y_coeff);
            //std::cout << "  A: y_coeff.two_norm() = " << y_coeff.two_norm() << std::endl;

            AOp_.apply(y_coeff,Ay_coeff);

            const auto& domain_ps = cast_ref<PSVC>(this->domain().creator());
            const auto& P = (cast_ref<PSVC>(domain_ps.subSpace(0).creator())).subSpace(0);

            ::Spacy::Vector result = zero(P);
            ::Spacy::Kaskade::copyFromCoefficientVector<VYSetDescription>(Ay_coeff,result);
            //std::cout << "  A: get(norm(result)) = " << get(norm(result)) << std::endl;

            return result;
        };

         ADualDual = [this](const ::Spacy::Vector& y) mutable
        {

            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
            CoefficientVectorP y_coeff(boost::fusion::at_c<2>(temp.data));
            CoefficientVectorP Ay_coeff(boost::fusion::at_c<2>(temp.data));

            ::Spacy::Kaskade::copyToCoefficientVector<VPSetDescription>(y,y_coeff);
            //std::cout << "  A: y_coeff.two_norm() = " << y_coeff.two_norm() << std::endl;

            AOp_.apply(y_coeff,Ay_coeff);

            const auto& domain_ps = cast_ref<PSVC>(this->domain().creator());
            const auto& P = (cast_ref<PSVC>(domain_ps.subSpace(1).creator())).subSpace(0);

            ::Spacy::Vector result = zero(P);
            ::Spacy::Kaskade::copyFromCoefficientVector<VPSetDescription>(Ay_coeff,result);
            //std::cout << "  A: get(norm(result)) = " << get(norm(result)) << std::endl;

            return result;
        };

        BBlock_ = [this](const ::Spacy::Vector& u) mutable
        {
            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
            CoefficientVectorU u_coeff(boost::fusion::at_c<1>(temp.data));
            CoefficientVectorP Bu_coeff(boost::fusion::at_c<2>(temp.data));

            ::Spacy::Kaskade::copyToCoefficientVector<VUSetDescription>(u,u_coeff);
            // std::cout << "  B: u_coeff.two_norm() = " << Bu_coeff.two_norm() << std::endl;
          //  auto startTime = std::chrono::high_resolution_clock::now();


            unsigned int uVectorSize = boost::fusion::at_c<0>(u_coeff.data).size();
            unsigned int uBlockSize =  boost::fusion::at_c<0>(u_coeff.data).operator[](0).size();

            unsigned int BuVectorSize = boost::fusion::at_c<0>(Bu_coeff.data).size();
            unsigned int BuBlockSize =  boost::fusion::at_c<0>(Bu_coeff.data).operator[](0).size();


            unsigned int uSize =  uVectorSize * uBlockSize;
            unsigned int BuSize =  BuVectorSize * BuBlockSize;

           // std::cout << "uSize: " << uSize << std::endl;
           // std::cout << "BuSize: " << BuSize << std::endl;

            std::vector<double> u_(uSize,0.0);
            std::vector<double> Bu_(BuSize,0.0);

            auto & matrix = B_.get_non_const();

            auto counter = 0u;

            for(auto & it: boost::fusion::at_c<0>(u_coeff.data))
                for(auto & iter: it)
                {
                    u_[counter] = iter;
                    counter++;
                }

            for(size_t i=0; i< matrix.ridx.size(); ++i)
                 Bu_[matrix.ridx[i]] += u_[matrix.cidx[i]]*matrix.data[i];

            counter = 0;

            for(auto & it: boost::fusion::at_c<0>(Bu_coeff.data))
                for(auto & iter: it)
                {
                    iter = Bu_[counter] ;
                    counter++;
                }


            //            int size = B_.getTriplet().ridx.size();
            //            int blockSize = boost::fusion::at_c<0>(u_coeff.data).operator[](0).size();

            //            for(int i = 0; i < size; i++)
            //            {
            //                int row = B_.getTriplet().ridx[i];
            //                int col = B_.getTriplet().cidx[i];

            //                int index1 = col / blockSize;
            //                int index2 = col % blockSize;

            //                int index3 = row /blockSize;
            //                int index4 = row % blockSize;

            //                (boost::fusion::at_c<0>(Bu_coeff.data).operator[](index3)).operator[](index4) +=  (boost::fusion::at_c<0>(u_coeff.data).operator[](index1)).operator[](index2)*B_.getTriplet().data[i];
            //            }

            //            auto max = 0;

            //            for(int i = 0; i< B_.getTriplet().cidx.size(); ++i)
            //            {
            //                max = std::max(max, B_.getTriplet().cidx[i]);
            //            }

            //            std::cout << "B columns: " << max << std::endl;

            //            max = 0;

            //            for(int i= 0; i<B_.getTriplet().ridx.size(); ++i)
            //            {
            //                max = std::max(max, B_.getTriplet().ridx[i]);
            //            }

            //            std::cout << "B rows: " << max << std::endl;

          //  B_.apply(u_coeff,Bu_coeff);

            //  std::cout << "Test: " << boost::fusion::at_c<0>(u_coeff.data).size() << std::endl;
            //  std::cout << "Test: " << boost::fusion::at_c<0>(Bu_coeff.data).size() << std::endl;

           // auto time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - startTime).count();

            // std::cout << "Elapsed Time B: " << time <<  "ms"  << std::endl;
            // std::cout << "copy " << std::endl;


            const auto& domain_ps = cast_ref<PSVC>(this->domain().creator());
            const auto& P = (cast_ref<PSVC>(domain_ps.subSpace(1).creator())).subSpace(0);

            ::Spacy::Vector result = zero(P);
            ::Spacy::Kaskade::copyFromCoefficientVector<VPSetDescription>(Bu_coeff,result);
            //std::cout << "  B: get(norm(result)) = " << get(norm(result)) << std::endl;

            return result;
        };

        BTBlock_ = [this](const ::Spacy::Vector& p) mutable
        {
            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
            CoefficientVectorP p_coeff(boost::fusion::at_c<2>(temp.data));
            CoefficientVectorU BTp_coeff(boost::fusion::at_c<1>(temp.data));

            ::Spacy::Kaskade::copyToCoefficientVector<VPSetDescription>(p,p_coeff);
            //std::cout << "  BT: p_coeff.two_norm() = " << p_coeff.two_norm() << std::endl;

            //  B_.applyscaleaddTransposed(1.0, p_coeff, BTp_coeff);

         //   auto startTime = std::chrono::high_resolution_clock::now();


            unsigned int pVectorSize = boost::fusion::at_c<0>(p_coeff.data).size();
            unsigned int pBlockSize =  boost::fusion::at_c<0>(p_coeff.data).operator[](0).size();

            unsigned int BTpVectorSize = boost::fusion::at_c<0>(BTp_coeff.data).size();
            unsigned int BTpBlockSize =  boost::fusion::at_c<0>(BTp_coeff.data).operator[](0).size();


            unsigned int pSize =   pVectorSize * pBlockSize;
            unsigned int BTpSize =  BTpVectorSize * BTpBlockSize;

           // std::cout << "pSize: " << pSize << std::endl;
           // std::cout << "BTpSize: " << BTpSize << std::endl;

            std::vector<double> p_(pSize,0.0);
            std::vector<double> BTp_(BTpSize,0.0);

            auto & matrix = B_t_.get_non_const();

            auto counter = 0u;

            for(auto & it: boost::fusion::at_c<0>(p_coeff.data))
                for(auto & iter: it)
                {
                    p_[counter] = iter;
                    counter++;
                }

            for(size_t i=0; i< matrix.ridx.size(); ++i)
                 BTp_[matrix.ridx[i]] += p_[matrix.cidx[i]]*matrix.data[i];

            counter = 0;

            for(auto & it: boost::fusion::at_c<0>(BTp_coeff.data))
                for(auto & iter: it)
                {
                    iter = BTp_[counter] ;
                    counter++;
                }


            //            int size = B_.getTriplet().ridx.size();
            //            int blockSize = boost::fusion::at_c<0>(u_coeff.data).operator[](0).size();

            //            for(int i = 0; i < size; i++)
            //            {
            //                int row = B_.getTriplet().ridx[i];
            //                int col = B_.getTriplet().cidx[i];

            //                int index1 = col / blockSize;
            //                int index2 = col % blockSize;

            //                int index3 = row /blockSize;
            //                int index4 = row % blockSize;

            //                (boost::fusion::at_c<0>(Bu_coeff.data).operator[](index3)).operator[](index4) +=  (boost::fusion::at_c<0>(u_coeff.data).operator[](index1)).operator[](index2)*B_.getTriplet().data[i];
            //            }

            //            auto max = 0;

            //            for(int i = 0; i< B_.getTriplet().cidx.size(); ++i)
            //            {
            //                max = std::max(max, B_.getTriplet().cidx[i]);
            //            }

            //            std::cout << "B columns: " << max << std::endl;

            //            max = 0;

            //            for(int i= 0; i<B_.getTriplet().ridx.size(); ++i)
            //            {
            //                max = std::max(max, B_.getTriplet().ridx[i]);
            //            }

            //            std::cout << "B rows: " << max << std::endl;

          //  B_t_.apply(u_coeff,Bu_coeff);

            //  std::cout << "Test: " << boost::fusion::at_c<0>(u_coeff.data).size() << std::endl;
            //  std::cout << "Test: " << boost::fusion::at_c<0>(Bu_coeff.data).size() << std::endl;

        //    auto time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - startTime).count();

           // std::cout << "Elapsed Time B: " << time <<  "ms"  << std::endl;
           // std::cout << "copy " << std::endl;

            const auto& domain_ps = cast_ref<PSVC>(this->domain().creator());
            const auto& U = (cast_ref<PSVC>(domain_ps.subSpace(0).creator())).subSpace(1);

            ::Spacy::Vector result = zero(U);
            ::Spacy::Kaskade::copyFromCoefficientVector<VUSetDescription>(BTp_coeff,result);
            //std::cout << "  BT: get(norm(result)) = " << get(norm(result)) << std::endl;

            return result;
        };

        BPXAPrim_ = [this](const ::Spacy::Vector& p) mutable
        {
            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
            CoefficientVectorY p_coeff(boost::fusion::at_c<0>(temp.data));
            CoefficientVectorY y(boost::fusion::at_c<0>(temp.data));

            ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(p,p_coeff);

          //  std::cout << "  BPX: p_coeff.two_norm() = " << p_coeff.two_norm() << std::endl;
            BPXAPtr_->apply(boost::fusion::at_c<0>(y.data),boost::fusion::at_c<0>(p_coeff.data));
          //  std::cout << "  BPX: y_coeff.two_norm() = " << y.two_norm() << std::endl;


            const auto& domain_ps = ::Spacy::cast_ref<PSVC>(this->domain().creator());
            const auto& Y = (::Spacy::cast_ref<PSVC>(domain_ps.subSpace(0).creator())).subSpace(0);

            ::Spacy::Vector result = zero(Y);
            ::Spacy::Kaskade::copyFromCoefficientVector<VYSetDescription>(y,result);
            //std::cout << "  BPX: get(norm(result)) = " << get(norm(result)) << std::endl;

            return result;
        };

        BPXADual_ = [this](const ::Spacy::Vector& y) mutable
        {
            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
            CoefficientVectorP y_coeff(boost::fusion::at_c<2>(temp.data));
            CoefficientVectorP p(boost::fusion::at_c<2>(temp.data));
//std::cout << "  BPXT: y_coeff.two_norm() = " << y_coeff.two_norm() << std::endl;
            ::Spacy::Kaskade::copyToCoefficientVector<VPSetDescription>(y,y_coeff);
  //          std::cout << "  BPXT: y_coeff.two_norm() = " << y_coeff.two_norm() << std::endl;

            BPXATPtr_->apply(boost::fusion::at_c<0>(p.data),boost::fusion::at_c<0>(y_coeff.data));

    //        std::cout << "  BPXT: p.two_norm() = " << p.two_norm() << std::endl;



            const auto& domain_ps = ::Spacy::cast_ref<PSVC>(this->domain().creator());
            const auto& P = (::Spacy::cast_ref<PSVC>(domain_ps.subSpace(1).creator())).subSpace(0);

            ::Spacy::Vector result = zero(P);
            ::Spacy::Kaskade::copyFromCoefficientVector<VPSetDescription>(p,result);
            //std::cout << "  BPXT: get(norm(result)) = " << get(norm(result)) << std::endl;

            return result;
        };

         BPXAPrimDual_ = [this](const ::Spacy::Vector& y) mutable
        {
          //  std::cout << "TestNormBeforeDualBPX: " <<  sqrt(y(y)) << std::endl;
            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
            CoefficientVectorY y_coeff(boost::fusion::at_c<0>(temp.data));
            CoefficientVectorP p(boost::fusion::at_c<2>(temp.data));

            ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(y,y_coeff);
      //      std::cout << "  BPXT: y_coeff.two_norm() = " << y_coeff.two_norm() << std::endl;

         //   std::cout << " TestNormBeforeDualBPX After copy " << y_coeff.two_norm() << std::endl;
            BPXATPtr_->apply(boost::fusion::at_c<0>(p.data),boost::fusion::at_c<0>(y_coeff.data));

        //   std::cout << "  BPXT: p.two_norm() = " << p.two_norm() << std::endl;



            const auto& domain_ps = ::Spacy::cast_ref<PSVC>(this->domain().creator());
            const auto& P = (::Spacy::cast_ref<PSVC>(domain_ps.subSpace(1).creator())).subSpace(0);

            ::Spacy::Vector result = zero(P);
            ::Spacy::Kaskade::copyFromCoefficientVector<VPSetDescription>(p,result);
            //std::cout << "  BPXT: get(norm(result)) = " << get(norm(result)) << std::endl;

            return result;
        };


        BPXA_ = [this](const ::Spacy::Vector& p) mutable
        {
           // std::cout << "Test0: " << std::endl;
            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
            CoefficientVectorP p_coeff(boost::fusion::at_c<2>(temp.data));
            CoefficientVectorY y(boost::fusion::at_c<0>(temp.data));


            //std::cout << "BPX: p_coeff.two_norm() = " << p_coeff.two_norm() << std::endl;
           // std::cout << "BPX: y.two_norm() = " << y.two_norm() << std::endl;
           // std::cout << "BPX: get(norm(p)) = " << get(norm(p)) << std::endl;

            ::Spacy::Kaskade::copyToCoefficientVector<VPSetDescription>(p,p_coeff);

           // std::cout << "BPX: p_coeff.two_norm() = " << p_coeff.two_norm() << std::endl;
            BPXAPtr_->apply(boost::fusion::at_c<0>(y.data),boost::fusion::at_c<0>(p_coeff.data));
           // std::cout << "y.two_norm() = " << y.two_norm() << std::endl;


            const auto& domain_ps = ::Spacy::cast_ref<PSVC>(this->domain().creator());
            const auto& Y = (::Spacy::cast_ref<PSVC>(domain_ps.subSpace(0).creator())).subSpace(0);

            ::Spacy::Vector result = zero(Y);
            ::Spacy::Kaskade::copyFromCoefficientVector<VYSetDescription>(y,result);
            //std::cout << "  BPX: get(norm(result)) = " << get(norm(result)) << std::endl;

            return result;
        };

        BPXAT_ = [this](const ::Spacy::Vector& y) mutable
        {
            //std::cout << "Test0: " << std::endl;
            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
            CoefficientVectorY y_coeff(boost::fusion::at_c<0>(temp.data));
            CoefficientVectorP p(boost::fusion::at_c<2>(temp.data));

           // std::cout << "BPXT: y_coeff.two_norm() = " << y_coeff.two_norm() << std::endl;
           // std::cout << "BPXT: p.two_norm() = " << p.two_norm() << std::endl;
           // std::cout << "BPXT: get(norm(y)) = " << get(norm(y)) << std::endl;

            ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(y,y_coeff);
           // std::cout << "BPXT: y_coeff.two_norm() = " << y_coeff.two_norm() << std::endl;

            for(auto i=0; i<boost::fusion::at_c<0>(p.data).size(); ++i)
            {
                for(auto j=0; j<boost::fusion::at_c<0>(p.data).operator[](i).size(); ++j)
                {
                    boost::fusion::at_c<0>(p.data).operator[](i).operator[](j) = boost::fusion::at_c<0>(y_coeff.data).operator[](i).operator[](j);
                }
            }


            y_coeff = 0.0;

          //  BPXATPtr_->apply(boost::fusion::at_c<0>(p.data),boost::fusion::at_c<0>(y_coeff.data));

            BPXAPtr_->apply(boost::fusion::at_c<0>(y_coeff.data),boost::fusion::at_c<0>(p.data));


          //  std::cout << "BPXT: p.two_norm() = " << y_coeff.two_norm() << std::endl << std::endl;


            for(auto i=0; i<boost::fusion::at_c<0>(p.data).size(); ++i)
            {
                for(auto j=0; j<boost::fusion::at_c<0>(p.data).operator[](i).size(); ++j)
                {
                    boost::fusion::at_c<0>(p.data).operator[](i).operator[](j) = boost::fusion::at_c<0>(y_coeff.data).operator[](i).operator[](j);
                }
            }



            const auto& domain_ps = ::Spacy::cast_ref<PSVC>(this->domain().creator());
            const auto& P = (::Spacy::cast_ref<PSVC>(domain_ps.subSpace(1).creator())).subSpace(0);

            ::Spacy::Vector result = zero(P);
            ::Spacy::Kaskade::copyFromCoefficientVector<VPSetDescription>(p,result);
            //std::cout << "  BPXT: get(norm(result)) = " << get(norm(result)) << std::endl;

            return result;
        };

        MuDiagMatrixSolver_ = [this](const CoefficientVectorU & u)
        {
            CoefficientVectorU result(u);

            unsigned int counter = 0;

            auto resultIter =  boost::fusion::at_c<0>(result.data).begin();

            for(const auto & it : boost::fusion::at_c<0>(u.data))
            {
                for(auto i = 0u; i < it.size(); i++)
                {
                    resultIter->operator[](i) = it.operator[](i)/MuDiag_[counter];
                    counter++;
                }

                resultIter++;
            }

            return result;
        };

        diagMuOperator_ = [this](const ::Spacy::Vector& u)
        {
            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(this->spaces()));
            CoefficientVectorU u_coeff(boost::fusion::at_c<1>(temp.data));
            ::Spacy::Kaskade::copyToCoefficientVector<VUSetDescription>(u,u_coeff);


            CoefficientVectorU resultCoeff(boost::fusion::at_c<1>(temp.data));

            unsigned int counter = 0;

            auto resultIter =  boost::fusion::at_c<0>(resultCoeff.data).begin();

            for(const auto & it : boost::fusion::at_c<0>(u_coeff.data))
            {
                for(auto i = 0u; i < it.size(); i++)
                {
                    resultIter->operator[](i) = it.operator[](i) *MuDiag_[counter];
                    counter++;
                }

                resultIter++;
            }

            ::Spacy::Vector result = u;
            ::Spacy::Kaskade::copyFromCoefficientVector<VUSetDescription>(resultCoeff,result);

            return result;
        };

        MySolver_ = [this](const ::Spacy::Vector& ry)
        {


            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(this->spaces()));
            CoefficientVectorY ry_coeff(boost::fusion::at_c<0>(temp.data));
            ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(ry,ry_coeff);

            CoefficientVectorY result_coeff(boost::fusion::at_c<0>(temp.data));

            unsigned int counter = 0;

            auto resultIter =  boost::fusion::at_c<0>(result_coeff.data).begin();



            for(const auto & it : boost::fusion::at_c<0>(ry_coeff.data))
            {
                for(auto i = 0u; i < it.size(); i++)
                {

                   // std::cout << "resultIter->operator[](i) : " << resultIter->operator[](i) << std::endl;
                    resultIter->operator[](i) = it.operator[](i)/MyDiag_[counter];
                   // std::cout << "Counter  MyDiag_[counter]: " << it.operator[](i)/MyDiag_[counter] << std::endl;
                    counter++;

                }

                resultIter++;
            }

            ::Spacy::Vector result = ry;

            ::Spacy::Kaskade::copyFromCoefficientVector<VYSetDescription>(result_coeff,result);

           // std::cout << "Return: " << std::endl;
            return ry;
        };




        controlSolver_ = [this](const ::Spacy::Vector& u)
        {
            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(this->spaces()));
            CoefficientVectorU u_coeff(boost::fusion::at_c<1>(temp.data));
            ::Spacy::Kaskade::copyToCoefficientVector<VUSetDescription>(u,u_coeff);

            CoefficientVectorU result_coeff(boost::fusion::at_c<1>(temp.data));
            if(useDirectMu)
                solMu_->apply(u_coeff, result_coeff);
            else
                result_coeff = MuDiagMatrixSolver_(u_coeff);

            ::Spacy::Vector result = u;
            ::Spacy::Kaskade::copyFromCoefficientVector<VUSetDescription>(result_coeff,result);

            return result;
        };

        ASolver_ = [this](const ::Spacy::Vector& p)
        {
            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(this->spaces()));
            CoefficientVectorY p_coeff(boost::fusion::at_c<0>(temp.data));
            ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(p,p_coeff);

            CoefficientVectorY result_coeff(boost::fusion::at_c<0>(temp.data));

                solA_->apply(p_coeff, result_coeff);

            ::Spacy::Vector result = p;
            ::Spacy::Kaskade::copyFromCoefficientVector<VYSetDescription>(result_coeff,result);

            return result;
        };

        ATSolver_ = [this](const ::Spacy::Vector& y)
        {

            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(this->spaces()));
            CoefficientVectorP y_coeff(boost::fusion::at_c<2>(temp.data));
            ::Spacy::Kaskade::copyToCoefficientVector<VPSetDescription>(y,y_coeff);

            CoefficientVectorP result_coeff(boost::fusion::at_c<2>(temp.data));

                solA_->apply(y_coeff, result_coeff);

            ::Spacy::Vector result = y;
            ::Spacy::Kaskade::copyFromCoefficientVector<VPSetDescription>(result_coeff,result);

            return result;
        };

        MyDirectSolver_ = [this](const ::Spacy::Vector& y)
            {

                CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(this->spaces()));
                CoefficientVectorY y_coeff(boost::fusion::at_c<0>(temp.data));
                ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(y,y_coeff);

                CoefficientVectorY result_coeff(boost::fusion::at_c<0>(temp.data));

                    solMy_->apply(y_coeff, result_coeff);

                ::Spacy::Vector result = y;
                ::Spacy::Kaskade::copyFromCoefficientVector<VYSetDescription>(result_coeff,result);

                return result;
            };



        LAPACKE_dstevr_func = [](int matrix_layout,char jobz,char range,int n,double* d,double* e,double vl,double vu,int il,int iu,double abstol,int* m,double* eig,double* z,int ldz,int* isuppz)
        {
            LAPACKE_dstevr(matrix_layout,jobz,range,n,d,e,vl,vu,il,iu,abstol,m,eig,z,ldz,isuppz);
        };

        normY_ = [this](const Spacy::Vector& y)
        {
            Spacy::Real max = -1;
            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
            CoefficientVectorY y_coeff(boost::fusion::at_c<0>(temp.data));

            ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(y,y_coeff);

            for (auto j =0u; j < boost::fusion::at_c<0>(y_coeff.data).size(); j++)
            {
                //1 Schleife reicht falls stateDim = 1
                if(fabs(boost::fusion::at_c<0>(y_coeff.data)[j][0]) > max)
                    max = fabs(boost::fusion::at_c<0>(y_coeff.data)[j][0]);
                ;
            }
            return max;
        };

        normP_ = [this](const Spacy::Vector& p)
        {
            Spacy::Real max = -1;
            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
            CoefficientVectorP p_coeff(boost::fusion::at_c<2>(temp.data));

            ::Spacy::Kaskade::copyToCoefficientVector<VPSetDescription>(p,p_coeff);

            for (auto j =0u; j < boost::fusion::at_c<0>(p_coeff.data).size(); j++)
            {
                //1 Schleife reicht falls stateDim = 1
                if(fabs(boost::fusion::at_c<0>(p_coeff.data)[j][0]) > max)
                    max = fabs(boost::fusion::at_c<0>(p_coeff.data)[j][0]);
                ;
            }
            return max;
        };

        blockHessianPtr_ = std::make_shared<::Spacy::ProductSpace::Operator_V2 >(std::vector<std::vector<CallableOperator>>{{MyBlock_,zeroBlock_,ATBlock_},{zeroBlock_,MuBlock_,BTBlock_},{ABlock_,BBlock_,zeroBlock_}}, this->domain(), this->domain());

    }

    /**
             * @brief Copy assignment.
             * @param g functional to copy from
             */

    C2BlockFunctional& operator=(const C2BlockFunctional& g)
    {
        setNumberOfThreads(g.getNumberOfThreads());
        f_ = g.f_;
        spaces_ = g.spaces_;
        gridManager_ = g.gridManager_;
        A_ = g.A_;
        value_ = g.value_;
        old_X_f_ = g.old_X_f_;
        old_X_df_ = g.old_X_df_;
        old_X_ddf_ = g.old_X_ddf_;
        rhs_ = g.rhs_;
        onlyLowerTriangle_ = g.onlyLowerTriangle_;
        rbegin_ = g.rbegin_;
        rend_ = g.rend_;
        cbegin_ = g.cbegin_;
        cend_ = g.cend_;
        solverCreator_ = g.solverCreator_;
        operatorSpace_ = g.operatorSpace_;

        My_ = g.My_;
        Myu_ = g.Myu_;
        Mu_ = g.Mu_;
        Muy_ = g.Muy_;
        AOp_ = g.AOp_;
        B_  = g.B_;
        AOp_t_ = g.AOp_t_;
        B_t_ = g.B_t_;
        BPXATPtr_ = g.BPXATPtr_;
        BPXAPtr_ = g.BPXAPtr_;

        MuDiag_ = g.MuDiag_;
        solMu_ =  g.solMu_;
        solA_  =  g.solA_,
        solAT_ =  g.solAT_,

        MyBlock_ = [this](const ::Spacy::Vector& y) mutable
        {
            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
            CoefficientVectorY y_coeff(boost::fusion::at_c<0>(temp.data));
            CoefficientVectorY My_coeff(boost::fusion::at_c<0>(temp.data));

            ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(y,y_coeff);
            //std::cout << "  My: y_coeff.two_norm() = " << y_coeff.two_norm() << std::endl;

            My_.apply(y_coeff,My_coeff);

            const auto& domain_ps = cast_ref<PSVC>(this->domain().creator());
            const auto& Y = (cast_ref<PSVC>(domain_ps.subSpace(0).creator())).subSpace(0);

            ::Spacy::Vector result = zero(Y);
            ::Spacy::Kaskade::copyFromCoefficientVector<VYSetDescription>(My_coeff,result);
            //std::cout << "  My: get(norm(result)) = " << get(norm(result)) << std::endl;

            return result;
        };

        MyuBlock_ = [this](const ::Spacy::Vector& u) mutable
        {
            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
            CoefficientVectorU u_coeff(boost::fusion::at_c<1>(temp.data));
            CoefficientVectorY Myu_coeff(boost::fusion::at_c<0>(temp.data));

            ::Spacy::Kaskade::copyToCoefficientVector<VUSetDescription>(u,u_coeff);
            //std::cout << "  Myu: u_coeff.two_norm() = " << u_coeff.two_norm() << std::endl;

            //Im Normalenschritt ist Myu_ leer -> gib einfach zero(...) zurück?
            if(Myu_.getTriplet().nrows() == 0)
            {
                const auto& domain_ps = cast_ref<PSVC>(this->domain().creator());
                const auto& Y = (cast_ref<PSVC>(domain_ps.subSpace(0).creator())).subSpace(0);
                return zero(Y);
            }
            Myu_.apply(u_coeff,Myu_coeff);

            const auto& domain_ps = cast_ref<PSVC>(this->domain().creator());
            const auto& Y = (cast_ref<PSVC>(domain_ps.subSpace(0).creator())).subSpace(0);

            ::Spacy::Vector result = zero(Y);
            ::Spacy::Kaskade::copyFromCoefficientVector<VYSetDescription>(Myu_coeff,result);
            //std::cout << "  Myu: get(norm(result)) = " << get(norm(result)) << std::endl;

            return result;
        };

        MuBlock_ = [this](const ::Spacy::Vector& u) mutable
        {
            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
            CoefficientVectorU u_coeff(boost::fusion::at_c<1>(temp.data));
            CoefficientVectorU Mu_coeff(boost::fusion::at_c<1>(temp.data));

            ::Spacy::Kaskade::copyToCoefficientVector<VUSetDescription>(u,u_coeff);
            //std::cout << "  Mu: u_coeff.two_norm() = " << u_coeff.two_norm() << std::endl;

            Mu_.apply(u_coeff,Mu_coeff);

            const auto& domain_ps = cast_ref<PSVC>(this->domain().creator());
            const auto& U = (cast_ref<PSVC>(domain_ps.subSpace(0).creator())).subSpace(1);

            ::Spacy::Vector result = zero(U);
            ::Spacy::Kaskade::copyFromCoefficientVector<VUSetDescription>(Mu_coeff,result);
            //std::cout << "  Mu: get(norm(result)) = " << get(norm(result)) << std::endl;

            return result;
        };

        MuyBlock_ = [this](const ::Spacy::Vector& y) mutable
        {
            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
            CoefficientVectorY y_coeff(boost::fusion::at_c<0>(temp.data));
            CoefficientVectorU Muy_coeff(boost::fusion::at_c<1>(temp.data));

            ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(y,y_coeff);
            //std::cout << "  Muy: y_coeff.two_norm() = " << y_coeff.two_norm() << std::endl;

            //Im Normalenschritt ist Muy_ leer -> gib einfach zero(...) zurück?
            if(Muy_.getTriplet().nrows() == 0)
            {
                const auto& domain_ps = cast_ref<PSVC>(this->domain().creator());
                const auto& U = (cast_ref<PSVC>(domain_ps.subSpace(0).creator())).subSpace(1);
                return zero(U);
            }
            Muy_.apply(y_coeff,Muy_coeff);

            const auto& domain_ps = cast_ref<PSVC>(this->domain().creator());
            const auto& U = (cast_ref<PSVC>(domain_ps.subSpace(0).creator())).subSpace(1);

            ::Spacy::Vector result = zero(U);
            ::Spacy::Kaskade::copyFromCoefficientVector<VUSetDescription>(Muy_coeff,result);
            //std::cout << "  Muy: get(norm(result)) = " << get(norm(result)) << std::endl;

            return result;
        };

        ATBlock_ = [this](const ::Spacy::Vector& p) mutable
        {
            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
            CoefficientVectorP p_coeff(boost::fusion::at_c<2>(temp.data));
            CoefficientVectorY ATp_coeff(boost::fusion::at_c<0>(temp.data));

            ::Spacy::Kaskade::copyToCoefficientVector<VPSetDescription>(p,p_coeff);
            //std::cout << "  AT: p_coeff.two_norm() = " << p_coeff.two_norm() << std::endl;

            AOp_t_.apply(p_coeff,ATp_coeff);

            const auto& domain_ps = cast_ref<PSVC>(this->domain().creator());
            const auto& Y = (cast_ref<PSVC>(domain_ps.subSpace(0).creator())).subSpace(0);

            ::Spacy::Vector result = zero(Y);
            ::Spacy::Kaskade::copyFromCoefficientVector<VYSetDescription>(ATp_coeff,result);
            //std::cout << "  AT: get(norm(result)) = " << get(norm(result)) << std::endl;

            return result;
        };

        ABlock_ = [this](const ::Spacy::Vector& y) mutable
        {
            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
            CoefficientVectorY y_coeff(boost::fusion::at_c<0>(temp.data));
            CoefficientVectorP Ay_coeff(boost::fusion::at_c<2>(temp.data));

            ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(y,y_coeff);
            //std::cout << "  A: y_coeff.two_norm() = " << y_coeff.two_norm() << std::endl;

            AOp_.apply(y_coeff,Ay_coeff);

            const auto& domain_ps = cast_ref<PSVC>(this->domain().creator());
            const auto& P = (cast_ref<PSVC>(domain_ps.subSpace(1).creator())).subSpace(0);

            ::Spacy::Vector result = zero(P);
            ::Spacy::Kaskade::copyFromCoefficientVector<VPSetDescription>(Ay_coeff,result);
            //std::cout << "  A: get(norm(result)) = " << get(norm(result)) << std::endl;

            return result;
        };

        APrimePrime = [this](const ::Spacy::Vector& y) mutable
       {
           CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
           CoefficientVectorY y_coeff(boost::fusion::at_c<0>(temp.data));
           CoefficientVectorY Ay_coeff(boost::fusion::at_c<0>(temp.data));

           ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(y,y_coeff);
           //std::cout << "  A: y_coeff.two_norm() = " << y_coeff.two_norm() << std::endl;

           AOp_.apply(y_coeff,Ay_coeff);

           const auto& domain_ps = cast_ref<PSVC>(this->domain().creator());
           const auto& P = (cast_ref<PSVC>(domain_ps.subSpace(0).creator())).subSpace(0);

           ::Spacy::Vector result = zero(P);
           ::Spacy::Kaskade::copyFromCoefficientVector<VYSetDescription>(Ay_coeff,result);
           //std::cout << "  A: get(norm(result)) = " << get(norm(result)) << std::endl;

           return result;
       };

        ADualDual = [this](const ::Spacy::Vector& y) mutable
       {
           CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
           CoefficientVectorP y_coeff(boost::fusion::at_c<2>(temp.data));
           CoefficientVectorP Ay_coeff(boost::fusion::at_c<2>(temp.data));

           ::Spacy::Kaskade::copyToCoefficientVector<VPSetDescription>(y,y_coeff);
           //std::cout << "  A: y_coeff.two_norm() = " << y_coeff.two_norm() << std::endl;

           AOp_.apply(y_coeff,Ay_coeff);

           const auto& domain_ps = cast_ref<PSVC>(this->domain().creator());
           const auto& P = (cast_ref<PSVC>(domain_ps.subSpace(1).creator())).subSpace(0);

           ::Spacy::Vector result = zero(P);
           ::Spacy::Kaskade::copyFromCoefficientVector<VPSetDescription>(Ay_coeff,result);
           //std::cout << "  A: get(norm(result)) = " << get(norm(result)) << std::endl;

           return result;
       };


        BBlock_ = [this](const ::Spacy::Vector& u) mutable
        {
            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
            CoefficientVectorU u_coeff(boost::fusion::at_c<1>(temp.data));
            CoefficientVectorP Bu_coeff(boost::fusion::at_c<2>(temp.data));

            ::Spacy::Kaskade::copyToCoefficientVector<VUSetDescription>(u,u_coeff);
            // std::cout << "Test: " << boost::fusion::at_c<0>(u_coeff.data).size() << std::endl;
            // std::cout << "Test: " << boost::fusion::at_c<0>(Bu_coeff.data).size() << std::endl;


            // std::cout << "  B: u_coeff.two_norm() = " << Bu_coeff.two_norm() << std::endl;
         //   auto startTime = std::chrono::high_resolution_clock::now();


            unsigned int uVectorSize = boost::fusion::at_c<0>(u_coeff.data).size();
            unsigned int uBlockSize =  boost::fusion::at_c<0>(u_coeff.data).operator[](0).size();

            unsigned int BuVectorSize = boost::fusion::at_c<0>(Bu_coeff.data).size();
            unsigned int BuBlockSize =  boost::fusion::at_c<0>(Bu_coeff.data).operator[](0).size();


            unsigned int uSize =  uVectorSize * uBlockSize;
            unsigned int BuSize =  BuVectorSize * BuBlockSize;

          //  std::cout << "uSize: " << uSize << std::endl;
          //  std::cout << "BuSize: " << BuSize << std::endl;

            std::vector<double> u_(uSize,0.0);
            std::vector<double> Bu_(BuSize,0.0);

            auto & matrix = B_.get_non_const();

            auto counter = 0u;

            for(auto & it: boost::fusion::at_c<0>(u_coeff.data))
                for(auto & iter: it)
                {
                    u_[counter] = iter;
                    counter++;
                }

            for(size_t i=0; i< matrix.ridx.size(); ++i)
                 Bu_[matrix.ridx[i]] += u_[matrix.cidx[i]]*matrix.data[i];

            counter = 0;

            for(auto & it: boost::fusion::at_c<0>(Bu_coeff.data))
                for(auto & iter: it)
                {
                    iter = Bu_[counter];
                    counter++;
                }


            //            int size = B_.getTriplet().ridx.size();
            //            int blockSize = boost::fusion::at_c<0>(u_coeff.data).operator[](0).size();

            //            for(int i = 0; i < size; i++)
            //            {
            //                int row = B_.getTriplet().ridx[i];
            //                int col = B_.getTriplet().cidx[i];

            //                int index1 = col / blockSize;
            //                int index2 = col % blockSize;

            //                int index3 = row /blockSize;
            //                int index4 = row % blockSize;

            //                (boost::fusion::at_c<0>(Bu_coeff.data).operator[](index3)).operator[](index4) +=  (boost::fusion::at_c<0>(u_coeff.data).operator[](index1)).operator[](index2)*B_.getTriplet().data[i];
            //            }

            //            auto max = 0;

            //            for(int i = 0; i< B_.getTriplet().cidx.size(); ++i)
            //            {
            //                max = std::max(max, B_.getTriplet().cidx[i]);
            //            }

            //            std::cout << "B columns: " << max << std::endl;

            //            max = 0;

            //            for(int i= 0; i<B_.getTriplet().ridx.size(); ++i)
            //            {
            //                max = std::max(max, B_.getTriplet().ridx[i]);
            //            }

            //            std::cout << "B rows: " << max << std::endl;

          //  B_.apply(u_coeff,Bu_coeff);

            //  std::cout << "Test: " << boost::fusion::at_c<0>(u_coeff.data).size() << std::endl;
            //  std::cout << "Test: " << boost::fusion::at_c<0>(Bu_coeff.data).size() << std::endl;

           // auto time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - startTime).count();

          //  std::cout << "Elapsed Time B: " << time <<  "ms"  << std::endl;
          //  std::cout << "copy " << std::endl;

            const auto& domain_ps = cast_ref<PSVC>(this->domain().creator());
            const auto& P = (cast_ref<PSVC>(domain_ps.subSpace(1).creator())).subSpace(0);

            ::Spacy::Vector result = zero(P);
            ::Spacy::Kaskade::copyFromCoefficientVector<VPSetDescription>(Bu_coeff,result);
            //std::cout << "  B: get(norm(result)) = " << get(norm(result)) << std::endl;

            return result;
        };

        BTBlock_ = [this](const ::Spacy::Vector& p) mutable
        {
            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
            CoefficientVectorP p_coeff(boost::fusion::at_c<2>(temp.data));
            CoefficientVectorU BTp_coeff(boost::fusion::at_c<1>(temp.data));

            ::Spacy::Kaskade::copyToCoefficientVector<VPSetDescription>(p,p_coeff);
            //std::cout << "  BT: p_coeff.two_norm() = " << p_coeff.two_norm() << std::endl;

            // B_.applyscaleaddTransposed(1.0, p_coeff, BTp_coeff);

          //  auto startTime = std::chrono::high_resolution_clock::now();


            unsigned int pVectorSize = boost::fusion::at_c<0>(p_coeff.data).size();
            unsigned int pBlockSize =  boost::fusion::at_c<0>(p_coeff.data).operator[](0).size();

            unsigned int BTpVectorSize = boost::fusion::at_c<0>(BTp_coeff.data).size();
            unsigned int BTpBlockSize =  boost::fusion::at_c<0>(BTp_coeff.data).operator[](0).size();


            unsigned int pSize =   pVectorSize * pBlockSize;
            unsigned int BTpSize =  BTpVectorSize * BTpBlockSize;

           // std::cout << "pSize: " << pSize << std::endl;
           // std::cout << "BTpSize: " << BTpSize << std::endl;

            std::vector<double> p_(pSize,0.0);
            std::vector<double> BTp_(BTpSize,0.0);

            auto & matrix = B_t_.get_non_const();

            auto counter = 0u;

            for(auto & it: boost::fusion::at_c<0>(p_coeff.data))
                for(auto & iter: it)
                {
                    p_[counter] = iter;
                    counter++;
                }

            for(size_t i=0; i< matrix.ridx.size(); ++i)
                 BTp_[matrix.ridx[i]] += p_[matrix.cidx[i]]*matrix.data[i];

            counter = 0;

            for(auto & it: boost::fusion::at_c<0>(BTp_coeff.data))
                for(auto & iter: it)
                {
                    iter = BTp_[counter] ;
                    counter++;
                }


            //            int size = B_.getTriplet().ridx.size();
            //            int blockSize = boost::fusion::at_c<0>(u_coeff.data).operator[](0).size();

            //            for(int i = 0; i < size; i++)
            //            {
            //                int row = B_.getTriplet().ridx[i];
            //                int col = B_.getTriplet().cidx[i];

            //                int index1 = col / blockSize;
            //                int index2 = col % blockSize;

            //                int index3 = row /blockSize;
            //                int index4 = row % blockSize;

            //                (boost::fusion::at_c<0>(Bu_coeff.data).operator[](index3)).operator[](index4) +=  (boost::fusion::at_c<0>(u_coeff.data).operator[](index1)).operator[](index2)*B_.getTriplet().data[i];
            //            }

            //            auto max = 0;

            //            for(int i = 0; i< B_.getTriplet().cidx.size(); ++i)
            //            {
            //                max = std::max(max, B_.getTriplet().cidx[i]);
            //            }

            //            std::cout << "B columns: " << max << std::endl;

            //            max = 0;

            //            for(int i= 0; i<B_.getTriplet().ridx.size(); ++i)
            //            {
            //                max = std::max(max, B_.getTriplet().ridx[i]);
            //            }

            //            std::cout << "B rows: " << max << std::endl;

          //  B_t_.apply(u_coeff,Bu_coeff);

            //  std::cout << "Test: " << boost::fusion::at_c<0>(u_coeff.data).size() << std::endl;
            //  std::cout << "Test: " << boost::fusion::at_c<0>(Bu_coeff.data).size() << std::endl;

         //   auto time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - startTime).count();

           // std::cout << "Elapsed Time B: " << time <<  "ms"  << std::endl;
           // std::cout << "copy " << std::endl;



            const auto& domain_ps = cast_ref<PSVC>(this->domain().creator());
            const auto& U = (cast_ref<PSVC>(domain_ps.subSpace(0).creator())).subSpace(1);

            ::Spacy::Vector result = zero(U);
            ::Spacy::Kaskade::copyFromCoefficientVector<VUSetDescription>(BTp_coeff,result);
            //std::cout << "  BT: get(norm(result)) = " << get(norm(result)) << std::endl;

            return result;
        };

        BPXAPrim_ = [this](const ::Spacy::Vector& p) mutable
        {
            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
            CoefficientVectorY p_coeff(boost::fusion::at_c<0>(temp.data));
            CoefficientVectorY y(boost::fusion::at_c<0>(temp.data));

            ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(p,p_coeff);

          //  std::cout << "  BPX: p_coeff.two_norm() = " << p_coeff.two_norm() << std::endl;
            BPXAPtr_->apply(boost::fusion::at_c<0>(y.data),boost::fusion::at_c<0>(p_coeff.data));
          //  std::cout << "  BPX: y_coeff.two_norm() = " << y.two_norm() << std::endl;


            const auto& domain_ps = ::Spacy::cast_ref<PSVC>(this->domain().creator());
            const auto& Y = (::Spacy::cast_ref<PSVC>(domain_ps.subSpace(0).creator())).subSpace(0);

            ::Spacy::Vector result = zero(Y);
            ::Spacy::Kaskade::copyFromCoefficientVector<VYSetDescription>(y,result);
            //std::cout << "  BPX: get(norm(result)) = " << get(norm(result)) << std::endl;

            return result;
        };

        BPXADual_ = [this](const ::Spacy::Vector& y) mutable
        {
            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
            CoefficientVectorP y_coeff(boost::fusion::at_c<2>(temp.data));
            CoefficientVectorP p(boost::fusion::at_c<2>(temp.data));
//std::cout << "  BPXT: y_coeff.two_norm() = " << y_coeff.two_norm() << std::endl;
            ::Spacy::Kaskade::copyToCoefficientVector<VPSetDescription>(y,y_coeff);
  //          std::cout << "  BPXT: y_coeff.two_norm() = " << y_coeff.two_norm() << std::endl;

            BPXATPtr_->apply(boost::fusion::at_c<0>(p.data),boost::fusion::at_c<0>(y_coeff.data));

    //        std::cout << "  BPXT: p.two_norm() = " << p.two_norm() << std::endl;



            const auto& domain_ps = ::Spacy::cast_ref<PSVC>(this->domain().creator());
            const auto& P = (::Spacy::cast_ref<PSVC>(domain_ps.subSpace(1).creator())).subSpace(0);

            ::Spacy::Vector result = zero(P);
            ::Spacy::Kaskade::copyFromCoefficientVector<VPSetDescription>(p,result);
            //std::cout << "  BPXT: get(norm(result)) = " << get(norm(result)) << std::endl;

            return result;
        };


        BPXAPrimDual_ = [this](const ::Spacy::Vector& y) mutable
       {
         //  std::cout << "TestNormBeforeDualBPX: " <<  sqrt(y(y)) << std::endl;
           CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
           CoefficientVectorY y_coeff(boost::fusion::at_c<0>(temp.data));
           CoefficientVectorP p(boost::fusion::at_c<2>(temp.data));

           ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(y,y_coeff);
        //   std::cout << "  BPXT: y_coeff.two_norm() = " << y_coeff.two_norm() << std::endl;

        //   std::cout << " TestNormBeforeDualBPX After copy " << y_coeff.two_norm() << std::endl;
           BPXATPtr_->apply(boost::fusion::at_c<0>(p.data),boost::fusion::at_c<0>(y_coeff.data));

       //   std::cout << "  BPXT: p.two_norm() = " << p.two_norm() << std::endl;



           const auto& domain_ps = ::Spacy::cast_ref<PSVC>(this->domain().creator());
           const auto& P = (::Spacy::cast_ref<PSVC>(domain_ps.subSpace(1).creator())).subSpace(0);

           ::Spacy::Vector result = zero(P);
           ::Spacy::Kaskade::copyFromCoefficientVector<VPSetDescription>(p,result);
           //std::cout << "  BPXT: get(norm(result)) = " << get(norm(result)) << std::endl;

           return result;
       };

        BPXA_ = [this](const ::Spacy::Vector& p) mutable
        {
            //std::cout << "Test2: " << std::endl;

            //  (this->blockHessianPtr_) = nullptr;
            //  this->domain();

            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
            CoefficientVectorP p_coeff(boost::fusion::at_c<2>(temp.data));
            CoefficientVectorY y(boost::fusion::at_c<0>(temp.data));

            //std::cout << "BPXT: p_coeff.two_norm() = " << p_coeff.two_norm() << std::endl;
            //std::cout << "BPXT: y.two_norm() = " << y.two_norm() << std::endl;
            //std::cout << "BPXT: get(norm(p)) = " << get(norm(p)) << std::endl;

            ::Spacy::Kaskade::copyToCoefficientVector<VPSetDescription>(p,p_coeff);
           // std::cout << "  BPX: p_coeff.two_norm() = " << p_coeff.two_norm() << std::endl;

            BPXAPtr_->apply(boost::fusion::at_c<0>(y.data),boost::fusion::at_c<0>(p_coeff.data));

            // std::cout << "  y.two_norm() = " << y.two_norm() << std::endl;


            const auto& domain_ps = ::Spacy::cast_ref<PSVC>(this->domain().creator());
            const auto& Y = (::Spacy::cast_ref<PSVC>(domain_ps.subSpace(0).creator())).subSpace(0);

            ::Spacy::Vector result = zero(Y);
            ::Spacy::Kaskade::copyFromCoefficientVector<VYSetDescription>(y,result);
            //std::cout << "  BPX: get(norm(result)) = " << get(norm(result)) << std::endl;

            return result;
        };

        BPXAT_ = [this](const ::Spacy::Vector& y) mutable
        {

          //  std::cout << "Test2: " << std::endl;
            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
            CoefficientVectorY y_coeff(boost::fusion::at_c<0>(temp.data));
            CoefficientVectorP p(boost::fusion::at_c<2>(temp.data));


            //std::cout << "BPXT: y_coeff.two_norm() = " << y_coeff.two_norm() << std::endl;
            //std::cout << "BPXT: p.two_norm() = " << p.two_norm() << std::endl;
            //std::cout << "BPXT: get(norm(y)) = " << get(norm(y)) << std::endl;


            ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(y,y_coeff);
          //  std::cout << "  BPXT: y_coeff.two_norm() = " << y_coeff.two_norm() << std::endl;

            BPXATPtr_->apply(boost::fusion::at_c<0>(p.data),boost::fusion::at_c<0>(y_coeff.data));

           // std::cout << "  BPXT: p.two_norm() = " << p.two_norm() << std::endl;

            const auto& domain_ps = ::Spacy::cast_ref<PSVC>(this->domain().creator());
            const auto& P = (::Spacy::cast_ref<PSVC>(domain_ps.subSpace(1).creator())).subSpace(0);

            ::Spacy::Vector result = zero(P);
            ::Spacy::Kaskade::copyFromCoefficientVector<VPSetDescription>(p,result);
            //std::cout << "  BPXT: get(norm(result)) = " << get(norm(result)) << std::endl;

            return result;
        };

        MuDiagMatrixSolver_ = [this](const CoefficientVectorU & u)
        {
            CoefficientVectorU result(u);

            unsigned int counter = 0;

            auto resultIter =  boost::fusion::at_c<0>(result.data).begin();

            for(const auto & it : boost::fusion::at_c<0>(u.data))
            {
                for(auto i = 0u; i < it.size(); i++)
                {
                    resultIter->operator[](i) = it.operator[](i)/MuDiag_[counter];
                    counter++;
                }

                resultIter++;
            }

            return result;
        };

        diagMuOperator_ = [this](const ::Spacy::Vector& u)
        {
            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(this->spaces()));
            CoefficientVectorU u_coeff(boost::fusion::at_c<1>(temp.data));
            ::Spacy::Kaskade::copyToCoefficientVector<VUSetDescription>(u,u_coeff);


            CoefficientVectorU resultCoeff(boost::fusion::at_c<1>(temp.data));

            unsigned int counter = 0;

            auto resultIter =  boost::fusion::at_c<0>(resultCoeff.data).begin();

            for(const auto & it : boost::fusion::at_c<0>(u_coeff.data))
            {
                for(auto i = 0u; i < it.size(); i++)
                {
                    resultIter->operator[](i) = it.operator[](i) *MuDiag_[counter];
                    counter++;
                }

                resultIter++;
            }

            ::Spacy::Vector result = u;
            ::Spacy::Kaskade::copyFromCoefficientVector<VUSetDescription>(resultCoeff,result);

            return result;
        };

        MySolver_ = [this](const ::Spacy::Vector& ry)
            {


                CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(this->spaces()));
                CoefficientVectorY ry_coeff(boost::fusion::at_c<0>(temp.data));
                ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(ry,ry_coeff);

                CoefficientVectorY result_coeff(boost::fusion::at_c<0>(temp.data));

                unsigned int counter = 0;

                auto resultIter =  boost::fusion::at_c<0>(result_coeff.data).begin();



                for(const auto & it : boost::fusion::at_c<0>(ry_coeff.data))
                {
                    for(auto i = 0u; i < it.size(); i++)
                    {

                       // std::cout << "resultIter->operator[](i) : " << resultIter->operator[](i) << std::endl;
                        resultIter->operator[](i) = it.operator[](i)/MyDiag_[counter];
                       // std::cout << "Counter  MyDiag_[counter]: " << it.operator[](i)/MyDiag_[counter] << std::endl;
                        counter++;

                    }

                    resultIter++;
                }

                ::Spacy::Vector result = ry;

                ::Spacy::Kaskade::copyFromCoefficientVector<VYSetDescription>(result_coeff,result);

               // std::cout << "Return: " << std::endl;
                return ry;
            };



        controlSolver_ = [this](const ::Spacy::Vector& u)
        {
            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(this->spaces()));
            CoefficientVectorU u_coeff(boost::fusion::at_c<1>(temp.data));
            ::Spacy::Kaskade::copyToCoefficientVector<VUSetDescription>(u,u_coeff);

            CoefficientVectorU result_coeff(boost::fusion::at_c<1>(temp.data));
            if(useDirectMu)
                solMu_->apply(u_coeff, result_coeff);
            else
                result_coeff = MuDiagMatrixSolver_(u_coeff);

            ::Spacy::Vector result = u;
            ::Spacy::Kaskade::copyFromCoefficientVector<VUSetDescription>(result_coeff,result);

            return result;
        };

        ASolver_ = [this](const ::Spacy::Vector& p)
        {
            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(this->spaces()));
            CoefficientVectorY p_coeff(boost::fusion::at_c<0>(temp.data));
            ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(p,p_coeff);

            CoefficientVectorY result_coeff(boost::fusion::at_c<0>(temp.data));

                solA_->apply(p_coeff, result_coeff);

            ::Spacy::Vector result = p;
            ::Spacy::Kaskade::copyFromCoefficientVector<VYSetDescription>(result_coeff,result);

            return result;
        };

        ATSolver_ = [this](const ::Spacy::Vector& y)
        {

            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(this->spaces()));
            CoefficientVectorP y_coeff(boost::fusion::at_c<2>(temp.data));
            ::Spacy::Kaskade::copyToCoefficientVector<VPSetDescription>(y,y_coeff);

            CoefficientVectorP result_coeff(boost::fusion::at_c<2>(temp.data));

                solA_->apply(y_coeff, result_coeff);

            ::Spacy::Vector result = y;
            ::Spacy::Kaskade::copyFromCoefficientVector<VPSetDescription>(result_coeff,result);

            return result;
        };

        MyDirectSolver_ = [this](const ::Spacy::Vector& y)
            {

                CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(this->spaces()));
                CoefficientVectorY y_coeff(boost::fusion::at_c<0>(temp.data));
                ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(y,y_coeff);

                CoefficientVectorY result_coeff(boost::fusion::at_c<0>(temp.data));

                    solMy_->apply(y_coeff, result_coeff);

                ::Spacy::Vector result = y;
                ::Spacy::Kaskade::copyFromCoefficientVector<VYSetDescription>(result_coeff,result);

                return result;
            };

        LAPACKE_dstevr_func = [](int matrix_layout,char jobz,char range,int n,double* d,double* e,double vl,double vu,int il,int iu,double abstol,int* m,double* eig,double* z,int ldz,int* isuppz)
        {
            LAPACKE_dstevr(matrix_layout,jobz,range,n,d,e,vl,vu,il,iu,abstol,m,eig,z,ldz,isuppz);
        };

        normY_ = [this](const Spacy::Vector& y)
        {
            Spacy::Real max = -1;
            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
            CoefficientVectorY y_coeff(boost::fusion::at_c<0>(temp.data));

            ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(y,y_coeff);

            for (auto j =0u; j < boost::fusion::at_c<0>(y_coeff.data).size(); j++)
            {
                //1 Schleife reicht falls stateDim = 1
                if(fabs(boost::fusion::at_c<0>(y_coeff.data)[j][0]) > max)
                    max = fabs(boost::fusion::at_c<0>(y_coeff.data)[j][0]);
                ;
            }
            return max;
        };

        normP_ = [this](const Spacy::Vector& p)
        {
            Spacy::Real max = -1;
            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
            CoefficientVectorP p_coeff(boost::fusion::at_c<2>(temp.data));

            ::Spacy::Kaskade::copyToCoefficientVector<VPSetDescription>(p,p_coeff);

            for (auto j =0u; j < boost::fusion::at_c<0>(p_coeff.data).size(); j++)
            {
                //1 Schleife reicht falls stateDim = 1
                if(fabs(boost::fusion::at_c<0>(p_coeff.data)[j][0]) > max)
                    max = fabs(boost::fusion::at_c<0>(p_coeff.data)[j][0]);
                ;
            }
            return max;
        };

        blockHessianPtr_ = std::make_shared<::Spacy::ProductSpace::Operator_V2 >(std::vector<std::vector<CallableOperator>>{{MyBlock_,zeroBlock_,ATBlock_},{zeroBlock_,MuBlock_,BTBlock_},{ABlock_,BBlock_,zeroBlock_}}, this->domain(), this->domain());

    }


    void reset()
    {
        reAssembleFValue = true;
        reAssembleFd1 = true;
        reAssembleFd2 = true;
    }


    //Todo use std::move in move constructor
    /**
             * @brief Move constructor.
             * @param g functional to move from
             */
    C2BlockFunctional(C2BlockFunctional&& g):FunctionalBase(g.domain()),
        NumberOfThreads(g),
        f_(g.f_), spaces_(g.spaces_),
        gridManager_(g.gridManager_),
        A_(g.A_),
        value_(g.value_),
        old_X_f_(g.old_X_f_),
        old_X_df_(g.old_X_df_),
        old_X_ddf_(g.old_X_ddf_),
        rhs_(g.rhs_),
        onlyLowerTriangle_(g.onlyLowerTriangle_),
        rbegin_(g.rbegin_), rend_(g.rend_), cbegin_(g.cbegin_), cend_(g.cend_),
        solverCreator_(g.solverCreator_),
        operatorSpace_(g.operatorSpace_),
        My_(g.My_),
        Myu_(g.Myu_),
        Mu_(g.Mu_),
        Muy_(g.Muy_),
        AOp_(g.AOp_),
        B_(g.B_),
        AOp_t_(g.AOp_t_),
        B_t_(g.B_t_),
        MuDiag_(g.MuDiag_),
        solMu_(g.solMu_),
        solA_(g.solA_),
        solAT_(g.solAT_),
        BPXATPtr_(g.BPXATPtr_),
        BPXAPtr_(g.BPXAPtr_)
    {


        MyBlock_ = [this](const ::Spacy::Vector& y) mutable
        {
            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
            CoefficientVectorY y_coeff(boost::fusion::at_c<0>(temp.data));
            CoefficientVectorY My_coeff(boost::fusion::at_c<0>(temp.data));

            ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(y,y_coeff);
            //std::cout << "  My: y_coeff.two_norm() = " << y_coeff.two_norm() << std::endl;

            My_.apply(y_coeff,My_coeff);

            const auto& domain_ps = cast_ref<PSVC>(this->domain().creator());
            const auto& Y = (cast_ref<PSVC>(domain_ps.subSpace(0).creator())).subSpace(0);

            ::Spacy::Vector result = zero(Y);
            ::Spacy::Kaskade::copyFromCoefficientVector<VYSetDescription>(My_coeff,result);
            //std::cout << "  My: get(norm(result)) = " << get(norm(result)) << std::endl;

            return result;
        };

        MyuBlock_ = [this](const ::Spacy::Vector& u) mutable
        {
            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
            CoefficientVectorU u_coeff(boost::fusion::at_c<1>(temp.data));
            CoefficientVectorY Myu_coeff(boost::fusion::at_c<0>(temp.data));

            ::Spacy::Kaskade::copyToCoefficientVector<VUSetDescription>(u,u_coeff);
            //std::cout << "  Myu: u_coeff.two_norm() = " << u_coeff.two_norm() << std::endl;

            //Im Normalenschritt ist Myu_ leer -> gib einfach zero(...) zurück?
            if(Myu_.getTriplet().nrows() == 0)
            {
                const auto& domain_ps = cast_ref<PSVC>(this->domain().creator());
                const auto& Y = (cast_ref<PSVC>(domain_ps.subSpace(0).creator())).subSpace(0);
                return zero(Y);
            }
            Myu_.apply(u_coeff,Myu_coeff);

            const auto& domain_ps = cast_ref<PSVC>(this->domain().creator());
            const auto& Y = (cast_ref<PSVC>(domain_ps.subSpace(0).creator())).subSpace(0);

            ::Spacy::Vector result = zero(Y);
            ::Spacy::Kaskade::copyFromCoefficientVector<VYSetDescription>(Myu_coeff,result);
            //std::cout << "  Myu: get(norm(result)) = " << get(norm(result)) << std::endl;

            return result;
        };

        MuBlock_ = [this](const ::Spacy::Vector& u) mutable
        {
            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
            CoefficientVectorU u_coeff(boost::fusion::at_c<1>(temp.data));
            CoefficientVectorU Mu_coeff(boost::fusion::at_c<1>(temp.data));

            ::Spacy::Kaskade::copyToCoefficientVector<VUSetDescription>(u,u_coeff);
            //std::cout << "  Mu: u_coeff.two_norm() = " << u_coeff.two_norm() << std::endl;

            Mu_.apply(u_coeff,Mu_coeff);

            const auto& domain_ps = cast_ref<PSVC>(this->domain().creator());
            const auto& U = (cast_ref<PSVC>(domain_ps.subSpace(0).creator())).subSpace(1);

            ::Spacy::Vector result = zero(U);
            ::Spacy::Kaskade::copyFromCoefficientVector<VUSetDescription>(Mu_coeff,result);
            //std::cout << "  Mu: get(norm(result)) = " << get(norm(result)) << std::endl;

            return result;
        };

        MuyBlock_ = [this](const ::Spacy::Vector& y) mutable
        {
            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
            CoefficientVectorY y_coeff(boost::fusion::at_c<0>(temp.data));
            CoefficientVectorU Muy_coeff(boost::fusion::at_c<1>(temp.data));

            ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(y,y_coeff);
            //std::cout << "  Muy: y_coeff.two_norm() = " << y_coeff.two_norm() << std::endl;

            //Im Normalenschritt ist Muy_ leer -> gib einfach zero(...) zurück?
            if(Muy_.getTriplet().nrows() == 0)
            {
                const auto& domain_ps = cast_ref<PSVC>(this->domain().creator());
                const auto& U = (cast_ref<PSVC>(domain_ps.subSpace(0).creator())).subSpace(1);
                return zero(U);
            }
            Muy_.apply(y_coeff,Muy_coeff);

            const auto& domain_ps = cast_ref<PSVC>(this->domain().creator());
            const auto& U = (cast_ref<PSVC>(domain_ps.subSpace(0).creator())).subSpace(1);

            ::Spacy::Vector result = zero(U);
            ::Spacy::Kaskade::copyFromCoefficientVector<VUSetDescription>(Muy_coeff,result);
            //std::cout << "  Muy: get(norm(result)) = " << get(norm(result)) << std::endl;

            return result;
        };

        ATBlock_ = [this](const ::Spacy::Vector& p) mutable
        {
            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
            CoefficientVectorP p_coeff(boost::fusion::at_c<2>(temp.data));
            CoefficientVectorY ATp_coeff(boost::fusion::at_c<0>(temp.data));

            ::Spacy::Kaskade::copyToCoefficientVector<VPSetDescription>(p,p_coeff);
            //std::cout << "  AT: p_coeff.two_norm() = " << p_coeff.two_norm() << std::endl;

            AOp_t_.apply(p_coeff,ATp_coeff);

            const auto& domain_ps = cast_ref<PSVC>(this->domain().creator());
            const auto& Y = (cast_ref<PSVC>(domain_ps.subSpace(0).creator())).subSpace(0);

            ::Spacy::Vector result = zero(Y);
            ::Spacy::Kaskade::copyFromCoefficientVector<VYSetDescription>(ATp_coeff,result);
            //std::cout << "  AT: get(norm(result)) = " << get(norm(result)) << std::endl;

            return result;
        };

        ABlock_ = [this](const ::Spacy::Vector& y) mutable
        {
            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
            CoefficientVectorY y_coeff(boost::fusion::at_c<0>(temp.data));
            CoefficientVectorP Ay_coeff(boost::fusion::at_c<2>(temp.data));

            ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(y,y_coeff);
            //std::cout << "  A: y_coeff.two_norm() = " << y_coeff.two_norm() << std::endl;

            AOp_.apply(y_coeff,Ay_coeff);

            const auto& domain_ps = cast_ref<PSVC>(this->domain().creator());
            const auto& P = (cast_ref<PSVC>(domain_ps.subSpace(1).creator())).subSpace(0);

            ::Spacy::Vector result = zero(P);
            ::Spacy::Kaskade::copyFromCoefficientVector<VPSetDescription>(Ay_coeff,result);
            //std::cout << "  A: get(norm(result)) = " << get(norm(result)) << std::endl;

            return result;
        };

        BBlock_ = [this](const ::Spacy::Vector& u) mutable
        {
            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
            CoefficientVectorU u_coeff(boost::fusion::at_c<1>(temp.data));
            CoefficientVectorP Bu_coeff(boost::fusion::at_c<2>(temp.data));

            ::Spacy::Kaskade::copyToCoefficientVector<VUSetDescription>(u,u_coeff);
            //   std::cout << "  B: u_coeff.two_norm() = " << Bu_coeff.two_norm() << std::endl;

            auto startTime = std::chrono::high_resolution_clock::now();


            unsigned int uVectorSize = boost::fusion::at_c<0>(u_coeff.data).size();
            unsigned int uBlockSize =  boost::fusion::at_c<0>(u_coeff.data).operator[](0).size();

            unsigned int BuVectorSize = boost::fusion::at_c<0>(Bu_coeff.data).size();
            unsigned int BuBlockSize =  boost::fusion::at_c<0>(Bu_coeff.data).operator[](0).size();


            unsigned int uSize =  uVectorSize * uBlockSize;
            unsigned int BuSize =  BuVectorSize * BuBlockSize;

           // std::cout << "uSize: " << uSize << std::endl;
           // std::cout << "BuSize: " << BuSize << std::endl;

            std::vector<double> u_(uSize,0.0);
            std::vector<double> Bu_(BuSize,0.0);

            auto & matrix = B_.get_non_const();

            auto counter = 0u;

            for(auto & it: boost::fusion::at_c<0>(u_coeff.data))
                for(auto & iter: it)
                {
                    u_[counter] = iter;
                    counter++;
                }

            for(size_t i=0; i<matrix.ridx.size(); ++i)
                 Bu_[matrix.ridx[i]] += u_[matrix.cidx[i]]*matrix.data[i];

            counter = 0;

            for(auto & it: boost::fusion::at_c<0>(Bu_coeff.data))
                for(auto & iter: it)
                {
                    iter = Bu_[counter];
                    counter++;
                }


            //            int size = B_.getTriplet().ridx.size();
            //            int blockSize = boost::fusion::at_c<0>(u_coeff.data).operator[](0).size();

            //            for(int i = 0; i < size; i++)
            //            {
            //                int row = B_.getTriplet().ridx[i];
            //                int col = B_.getTriplet().cidx[i];

            //                int index1 = col / blockSize;
            //                int index2 = col % blockSize;

            //                int index3 = row /blockSize;
            //                int index4 = row % blockSize;

            //                (boost::fusion::at_c<0>(Bu_coeff.data).operator[](index3)).operator[](index4) +=  (boost::fusion::at_c<0>(u_coeff.data).operator[](index1)).operator[](index2)*B_.getTriplet().data[i];
            //            }

            //            auto max = 0;

            //            for(int i = 0; i< B_.getTriplet().cidx.size(); ++i)
            //            {
            //                max = std::max(max, B_.getTriplet().cidx[i]);
            //            }

            //            std::cout << "B columns: " << max << std::endl;

            //            max = 0;

            //            for(int i= 0; i<B_.getTriplet().ridx.size(); ++i)
            //            {
            //                max = std::max(max, B_.getTriplet().ridx[i]);
            //            }

            //            std::cout << "B rows: " << max << std::endl;

          //  B_.apply(u_coeff,Bu_coeff);

            //  std::cout << "Test: " << boost::fusion::at_c<0>(u_coeff.data).size() << std::endl;
            //  std::cout << "Test: " << boost::fusion::at_c<0>(Bu_coeff.data).size() << std::endl;

           // auto time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - startTime).count();

          //  std::cout << "Elapsed Time B: " << time <<  "ms"  << std::endl;
          //  std::cout << "copy " << std::endl;

            const auto& domain_ps = cast_ref<PSVC>(this->domain().creator());
            const auto& P = (cast_ref<PSVC>(domain_ps.subSpace(1).creator())).subSpace(0);

            ::Spacy::Vector result = zero(P);
            ::Spacy::Kaskade::copyFromCoefficientVector<VPSetDescription>(Bu_coeff,result);
            //std::cout << "  B: get(norm(result)) = " << get(norm(result)) << std::endl;

            return result;
        };
        APrimePrime = [this](const ::Spacy::Vector& y) mutable
       {
           CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
           CoefficientVectorY y_coeff(boost::fusion::at_c<0>(temp.data));
           CoefficientVectorY Ay_coeff(boost::fusion::at_c<0>(temp.data));

           ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(y,y_coeff);
           //std::cout << "  A: y_coeff.two_norm() = " << y_coeff.two_norm() << std::endl;

           AOp_.apply(y_coeff,Ay_coeff);

           const auto& domain_ps = cast_ref<PSVC>(this->domain().creator());
           const auto& P = (cast_ref<PSVC>(domain_ps.subSpace(0).creator())).subSpace(0);

           ::Spacy::Vector result = zero(P);
           ::Spacy::Kaskade::copyFromCoefficientVector<VYSetDescription>(Ay_coeff,result);
           //std::cout << "  A: get(norm(result)) = " << get(norm(result)) << std::endl;

           return result;
       };

        ADualDual = [this](const ::Spacy::Vector& y) mutable
       {
           CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
           CoefficientVectorP y_coeff(boost::fusion::at_c<2>(temp.data));
           CoefficientVectorP Ay_coeff(boost::fusion::at_c<2>(temp.data));

           ::Spacy::Kaskade::copyToCoefficientVector<VPSetDescription>(y,y_coeff);
           //std::cout << "  A: y_coeff.two_norm() = " << y_coeff.two_norm() << std::endl;

           AOp_.apply(y_coeff,Ay_coeff);

           const auto& domain_ps = cast_ref<PSVC>(this->domain().creator());
           const auto& P = (cast_ref<PSVC>(domain_ps.subSpace(1).creator())).subSpace(0);

           ::Spacy::Vector result = zero(P);
           ::Spacy::Kaskade::copyFromCoefficientVector<VPSetDescription>(Ay_coeff,result);
           //std::cout << "  A: get(norm(result)) = " << get(norm(result)) << std::endl;

           return result;
       };


        BTBlock_ = [this](const ::Spacy::Vector& p) mutable
        {
            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
            CoefficientVectorP p_coeff(boost::fusion::at_c<2>(temp.data));
            CoefficientVectorU BTp_coeff(boost::fusion::at_c<1>(temp.data));

            ::Spacy::Kaskade::copyToCoefficientVector<VPSetDescription>(p,p_coeff);
            //std::cout << "  BT: p_coeff.two_norm() = " << p_coeff.two_norm() << std::endl;

            auto startTime = std::chrono::high_resolution_clock::now();


            unsigned int pVectorSize = boost::fusion::at_c<0>(p_coeff.data).size();
            unsigned int pBlockSize =  boost::fusion::at_c<0>(p_coeff.data).operator[](0).size();

            unsigned int BTpVectorSize = boost::fusion::at_c<0>(BTp_coeff.data).size();
            unsigned int BTpBlockSize =  boost::fusion::at_c<0>(BTp_coeff.data).operator[](0).size();


            unsigned int pSize =   pVectorSize * pBlockSize;
            unsigned int BTpSize =  BTpVectorSize * BTpBlockSize;

          //  std::cout << "pSize: " << pSize << std::endl;
          //  std::cout << "BTpSize: " << BTpSize << std::endl;

            std::vector<double> p_(pSize,0.0);
            std::vector<double> BTp_(BTpSize,0.0);

            auto & matrix = B_t_.get_non_const();

            auto counter = 0u;

            for(auto & it: boost::fusion::at_c<0>(p_coeff.data))
                for(auto & iter: it)
                {
                    p_[counter] = iter;
                    counter++;
                }

            for(size_t i=0; i< matrix.ridx.size(); ++i)
                 BTp_[matrix.ridx[i]] += p_[matrix.cidx[i]]*matrix.data[i];

            counter = 0;

            for(auto & it: boost::fusion::at_c<0>(BTp_coeff.data))
                for(auto & iter: it)
                {
                    iter = BTp_[counter] ;
                    counter++;
                }


            //            int size = B_.getTriplet().ridx.size();
            //            int blockSize = boost::fusion::at_c<0>(u_coeff.data).operator[](0).size();

            //            for(int i = 0; i < size; i++)
            //            {
            //                int row = B_.getTriplet().ridx[i];
            //                int col = B_.getTriplet().cidx[i];

            //                int index1 = col / blockSize;
            //                int index2 = col % blockSize;

            //                int index3 = row /blockSize;
            //                int index4 = row % blockSize;

            //                (boost::fusion::at_c<0>(Bu_coeff.data).operator[](index3)).operator[](index4) +=  (boost::fusion::at_c<0>(u_coeff.data).operator[](index1)).operator[](index2)*B_.getTriplet().data[i];
            //            }

            //            auto max = 0;

            //            for(int i = 0; i< B_.getTriplet().cidx.size(); ++i)
            //            {
            //                max = std::max(max, B_.getTriplet().cidx[i]);
            //            }

            //            std::cout << "B columns: " << max << std::endl;

            //            max = 0;

            //            for(int i= 0; i<B_.getTriplet().ridx.size(); ++i)
            //            {
            //                max = std::max(max, B_.getTriplet().ridx[i]);
            //            }

            //            std::cout << "B rows: " << max << std::endl;

          //  B_t_.apply(u_coeff,Bu_coeff);

            //  std::cout << "Test: " << boost::fusion::at_c<0>(u_coeff.data).size() << std::endl;
            //  std::cout << "Test: " << boost::fusion::at_c<0>(Bu_coeff.data).size() << std::endl;

         //  auto time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - startTime).count();

           // std::cout << "Elapsed Time B: " << time <<  "ms"  << std::endl;
           // std::cout << "copy " << std::endl;

            const auto& domain_ps = cast_ref<PSVC>(this->domain().creator());
            const auto& U = (cast_ref<PSVC>(domain_ps.subSpace(0).creator())).subSpace(1);

            ::Spacy::Vector result = zero(U);
            ::Spacy::Kaskade::copyFromCoefficientVector<VUSetDescription>(BTp_coeff,result);
            //std::cout << "  BT: get(norm(result)) = " << get(norm(result)) << std::endl;

            return result;
        };

         BPXAPrim_ = [this](const ::Spacy::Vector& p) mutable
        {
            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
            CoefficientVectorY p_coeff(boost::fusion::at_c<0>(temp.data));
            CoefficientVectorY y(boost::fusion::at_c<0>(temp.data));

            ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(p,p_coeff);

          //  std::cout << "  BPX: p_coeff.two_norm() = " << p_coeff.two_norm() << std::endl;
            BPXAPtr_->apply(boost::fusion::at_c<0>(y.data),boost::fusion::at_c<0>(p_coeff.data));
          //  std::cout << "  BPX: y_coeff.two_norm() = " << y.two_norm() << std::endl;


            const auto& domain_ps = ::Spacy::cast_ref<PSVC>(this->domain().creator());
            const auto& Y = (::Spacy::cast_ref<PSVC>(domain_ps.subSpace(0).creator())).subSpace(0);

            ::Spacy::Vector result = zero(Y);
            ::Spacy::Kaskade::copyFromCoefficientVector<VYSetDescription>(y,result);
            //std::cout << "  BPX: get(norm(result)) = " << get(norm(result)) << std::endl;

            return result;
        };

         BPXADual_ = [this](const ::Spacy::Vector& y) mutable
        {
            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
            CoefficientVectorP y_coeff(boost::fusion::at_c<2>(temp.data));
            CoefficientVectorP p(boost::fusion::at_c<2>(temp.data));
//std::cout << "  BPXT: y_coeff.two_norm() = " << y_coeff.two_norm() << std::endl;
            ::Spacy::Kaskade::copyToCoefficientVector<VPSetDescription>(y,y_coeff);
  //          std::cout << "  BPXT: y_coeff.two_norm() = " << y_coeff.two_norm() << std::endl;

            BPXATPtr_->apply(boost::fusion::at_c<0>(p.data),boost::fusion::at_c<0>(y_coeff.data));

    //        std::cout << "  BPXT: p.two_norm() = " << p.two_norm() << std::endl;



            const auto& domain_ps = ::Spacy::cast_ref<PSVC>(this->domain().creator());
            const auto& P = (::Spacy::cast_ref<PSVC>(domain_ps.subSpace(1).creator())).subSpace(0);

            ::Spacy::Vector result = zero(P);
            ::Spacy::Kaskade::copyFromCoefficientVector<VPSetDescription>(p,result);
            //std::cout << "  BPXT: get(norm(result)) = " << get(norm(result)) << std::endl;

            return result;
        };

        BPXAPrimDual_ = [this](const ::Spacy::Vector& y) mutable
        {
          //  std::cout << "TestNormBeforeDualBPX: " <<  sqrt(y(y)) << std::endl;
            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
            CoefficientVectorY y_coeff(boost::fusion::at_c<0>(temp.data));
            CoefficientVectorP p(boost::fusion::at_c<2>(temp.data));

            ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(y,y_coeff);
      //      std::cout << "  BPXT: y_coeff.two_norm() = " << y_coeff.two_norm() << std::endl;

         //   std::cout << " TestNormBeforeDualBPX After copy " << y_coeff.two_norm() << std::endl;
            BPXATPtr_->apply(boost::fusion::at_c<0>(p.data),boost::fusion::at_c<0>(y_coeff.data));

        //   std::cout << "  BPXT: p.two_norm() = " << p.two_norm() << std::endl;



            const auto& domain_ps = ::Spacy::cast_ref<PSVC>(this->domain().creator());
            const auto& P = (::Spacy::cast_ref<PSVC>(domain_ps.subSpace(1).creator())).subSpace(0);

            ::Spacy::Vector result = zero(P);
            ::Spacy::Kaskade::copyFromCoefficientVector<VPSetDescription>(p,result);
            //std::cout << "  BPXT: get(norm(result)) = " << get(norm(result)) << std::endl;

            return result;
        };




        BPXA_ = [this](const ::Spacy::Vector& p) mutable
        {


           // std::cout << "Test3: " << std::endl;

            //  (this->blockHessianPtr_) = nullptr;
            //  this->domain();

            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
            CoefficientVectorP p_coeff(boost::fusion::at_c<2>(temp.data));
            CoefficientVectorY y(boost::fusion::at_c<0>(temp.data));


           // std::cout << "BPX: p_coeff.two_norm() = " << p_coeff.two_norm() << std::endl;
          //  std::cout << "BPX: y.two_norm() = " << y.two_norm() << std::endl;
          //  std::cout << "BPX: get(norm(p)) = " << get(norm(p)) << std::endl;

            ::Spacy::Kaskade::copyToCoefficientVector<VPSetDescription>(p,p_coeff);

            //std::cout << "  BPX: p_coeff.two_norm() = " << p_coeff.two_norm() << std::endl;
            BPXAPtr_->apply(boost::fusion::at_c<0>(y.data),boost::fusion::at_c<0>(p_coeff.data));
            // std::cout << "  y.two_norm() = " << y.two_norm() << std::endl;

            const auto& domain_ps = ::Spacy::cast_ref<PSVC>(this->domain().creator());
            const auto& Y = (::Spacy::cast_ref<PSVC>(domain_ps.subSpace(0).creator())).subSpace(0);

            ::Spacy::Vector result = zero(Y);
            ::Spacy::Kaskade::copyFromCoefficientVector<VYSetDescription>(y,result);
            //std::cout << "  BPX: get(norm(result)) = " << get(norm(result)) << std::endl;

            return result;
        };

        BPXAT_ = [this](const ::Spacy::Vector& y) mutable
        {

            //std::cout << "Test3: " << std::endl;
            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
            CoefficientVectorY y_coeff(boost::fusion::at_c<0>(temp.data));
            CoefficientVectorP p(boost::fusion::at_c<2>(temp.data));


           // std::cout << "BPXT: y_coeff.two_norm() = " << y_coeff.two_norm() << std::endl;
           // std::cout << "BPXT: p.two_norm() = " << p.two_norm() << std::endl;
           // std::cout << "BPXT: get(norm(y)) = " << get(norm(y)) << std::endl;


            ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(y,y_coeff);

          //  std::cout << "  BPXT: y_coeff.two_norm() = " << y_coeff.two_norm() << std::endl;

            BPXATPtr_->apply(boost::fusion::at_c<0>(p.data),boost::fusion::at_c<0>(y_coeff.data));

           // std::cout << "  BPXT: p.two_norm() = " << p.two_norm() << std::endl;

            const auto& domain_ps = ::Spacy::cast_ref<PSVC>(this->domain().creator());
            const auto& P = (::Spacy::cast_ref<PSVC>(domain_ps.subSpace(1).creator())).subSpace(0);

            ::Spacy::Vector result = zero(P);
            ::Spacy::Kaskade::copyFromCoefficientVector<VPSetDescription>(p,result);
            //std::cout << "  BPXT: get(norm(result)) = " << get(norm(result)) << std::endl;

            return result;
        };

        MuDiagMatrixSolver_ = [this](const CoefficientVectorU & u)
        {
            CoefficientVectorU result(u);

            unsigned int counter = 0;

            auto resultIter =  boost::fusion::at_c<0>(result.data).begin();

            for(const auto & it : boost::fusion::at_c<0>(u.data))
            {
                for(auto i = 0u; i < it.size(); i++)
                {
                    resultIter->operator[](i) = it.operator[](i)/MuDiag_[counter];
                    counter++;
                }

                resultIter++;
            }

            return result;
        };

        diagMuOperator_ = [this](const ::Spacy::Vector& u)
        {
            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(this->spaces()));
            CoefficientVectorU u_coeff(boost::fusion::at_c<1>(temp.data));
            ::Spacy::Kaskade::copyToCoefficientVector<VUSetDescription>(u,u_coeff);


            CoefficientVectorU resultCoeff(boost::fusion::at_c<1>(temp.data));

            unsigned int counter = 0;

            auto resultIter =  boost::fusion::at_c<0>(resultCoeff.data).begin();

            for(const auto & it : boost::fusion::at_c<0>(u_coeff.data))
            {
                for(auto i = 0u; i < it.size(); i++)
                {
                    resultIter->operator[](i) = it.operator[](i) *MuDiag_[counter];
                    counter++;
                }

                resultIter++;
            }

            ::Spacy::Vector result = u;
            ::Spacy::Kaskade::copyFromCoefficientVector<VUSetDescription>(resultCoeff,result);

            return result;
        };

        MySolver_ = [this](const ::Spacy::Vector& ry)
            {


                CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(this->spaces()));
                CoefficientVectorY ry_coeff(boost::fusion::at_c<0>(temp.data));
                ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(ry,ry_coeff);

                CoefficientVectorY result_coeff(boost::fusion::at_c<0>(temp.data));

                unsigned int counter = 0;

                auto resultIter =  boost::fusion::at_c<0>(result_coeff.data).begin();



                for(const auto & it : boost::fusion::at_c<0>(ry_coeff.data))
                {
                    for(auto i = 0u; i < it.size(); i++)
                    {

                       // std::cout << "resultIter->operator[](i) : " << resultIter->operator[](i) << std::endl;
                        resultIter->operator[](i) = it.operator[](i)/MyDiag_[counter];
                       // std::cout << "Counter  MyDiag_[counter]: " << it.operator[](i)/MyDiag_[counter] << std::endl;
                        counter++;

                    }

                    resultIter++;
                }

                ::Spacy::Vector result = ry;

                ::Spacy::Kaskade::copyFromCoefficientVector<VYSetDescription>(result_coeff,result);

               // std::cout << "Return: " << std::endl;
                return ry;
            };



        controlSolver_ = [this](const ::Spacy::Vector& u)
        {
            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(this->spaces()));
            CoefficientVectorU u_coeff(boost::fusion::at_c<1>(temp.data));
            ::Spacy::Kaskade::copyToCoefficientVector<VUSetDescription>(u,u_coeff);

            CoefficientVectorU result_coeff(boost::fusion::at_c<1>(temp.data));
            if(useDirectMu)
                solMu_->apply(u_coeff, result_coeff);
            else
                result_coeff = MuDiagMatrixSolver_(u_coeff);

            ::Spacy::Vector result = u;
            ::Spacy::Kaskade::copyFromCoefficientVector<VUSetDescription>(result_coeff,result);

            return result;
        };

        ASolver_ = [this](const ::Spacy::Vector& p)
        {
            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(this->spaces()));
            CoefficientVectorY p_coeff(boost::fusion::at_c<0>(temp.data));
            ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(p,p_coeff);

            CoefficientVectorY result_coeff(boost::fusion::at_c<0>(temp.data));

                solA_->apply(p_coeff, result_coeff);

            ::Spacy::Vector result = p;
            ::Spacy::Kaskade::copyFromCoefficientVector<VYSetDescription>(result_coeff,result);

            return result;
        };

        ATSolver_ = [this](const ::Spacy::Vector& y)
        {

            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(this->spaces()));
            CoefficientVectorP y_coeff(boost::fusion::at_c<2>(temp.data));
            ::Spacy::Kaskade::copyToCoefficientVector<VPSetDescription>(y,y_coeff);

            CoefficientVectorP result_coeff(boost::fusion::at_c<2>(temp.data));

                solA_->apply(y_coeff, result_coeff);

            ::Spacy::Vector result = y;
            ::Spacy::Kaskade::copyFromCoefficientVector<VPSetDescription>(result_coeff,result);

            return result;
        };

        MyDirectSolver_ = [this](const ::Spacy::Vector& y)
            {

                CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(this->spaces()));
                CoefficientVectorY y_coeff(boost::fusion::at_c<0>(temp.data));
                ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(y,y_coeff);

                CoefficientVectorY result_coeff(boost::fusion::at_c<0>(temp.data));

                    solMy_->apply(y_coeff, result_coeff);

                ::Spacy::Vector result = y;
                ::Spacy::Kaskade::copyFromCoefficientVector<VYSetDescription>(result_coeff,result);

                return result;
            };

        LAPACKE_dstevr_func = [](int matrix_layout,char jobz,char range,int n,double* d,double* e,double vl,double vu,int il,int iu,double abstol,int* m,double* eig,double* z,int ldz,int* isuppz)
        {
            LAPACKE_dstevr(matrix_layout,jobz,range,n,d,e,vl,vu,il,iu,abstol,m,eig,z,ldz,isuppz);
        };

        normY_ = [this](const Spacy::Vector& y)
        {
            Spacy::Real max = -1;
            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
            CoefficientVectorY y_coeff(boost::fusion::at_c<0>(temp.data));

            ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(y,y_coeff);

            for (auto j =0u; j < boost::fusion::at_c<0>(y_coeff.data).size(); j++)
            {
                //1 Schleife reicht falls stateDim = 1
                if(fabs(boost::fusion::at_c<0>(y_coeff.data)[j][0]) > max)
                    max = fabs(boost::fusion::at_c<0>(y_coeff.data)[j][0]);
                ;
            }
            return max;
        };

        normP_ = [this](const Spacy::Vector& p)
        {
            Spacy::Real max = -1;
            CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
            CoefficientVectorP p_coeff(boost::fusion::at_c<2>(temp.data));

            ::Spacy::Kaskade::copyToCoefficientVector<VPSetDescription>(p,p_coeff);

            for (auto j =0u; j < boost::fusion::at_c<0>(p_coeff.data).size(); j++)
            {
                //1 Schleife reicht falls stateDim = 1
                if(fabs(boost::fusion::at_c<0>(p_coeff.data)[j][0]) > max)
                    max = fabs(boost::fusion::at_c<0>(p_coeff.data)[j][0]);
                ;
            }
            return max;
        };

        blockHessianPtr_ = std::make_shared<::Spacy::ProductSpace::Operator_V2 >(std::vector<std::vector<CallableOperator>>{{MyBlock_,zeroBlock_,ATBlock_},{zeroBlock_,MuBlock_,BTBlock_},{ABlock_,BBlock_,zeroBlock_}}, this->domain(), this->domain());


    }

    /**
             * @brief Move assignment.
             * @param g functional to move from
             */
    C2BlockFunctional& operator=(C2BlockFunctional&& g) = delete;

    /**
             * @brief Apply functional.
             * @param x argument
             * @return \f$f(x)\f$
             */
    Real operator()(const ::Spacy::Vector& x) const
    {
        assembleFunctional(x);

        return value_;
    }

    /**
             * @brief Compute first directional derivative \f$f'(x) \in X^* \f$.
             *
             * @param x current iterate
             * @return \f$f'(x)\f$
             */
    ::Spacy::Vector d1(const ::Spacy::Vector& x) const
    {
        assembleGradient(x);

        return rhs_;
    }

    /**
             * @brief Compute second directional derivative \f$f''(x)dx\in X^* \f$.
             *
             * @param x current iterate
             * @param dx perturbation
             * @return \f$f''(x)dx\f$
             */
    ::Spacy::Vector d2(const ::Spacy::Vector& x, const ::Spacy::Vector& dx) const
    {
        assembleHessian(x);

        CoefficientVector dx_( VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_) );
        copyToCoefficientVector<VariableSetDescription>(dx,dx_);
        CoefficientVector y_( VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_) );

        A_.apply( dx_ , y_ );

        auto y = zero( domain().dualSpace() );
        copyFromCoefficientVector<VariableSetDescription>(y_,y);

        return y;
    }

    /**
             * @brief Access \f$f''(x)\f$ as linear operator \f$X\rightarrow X^*\f$.
             * @param x point of linearization
             * @see LinearOperator, ::Spacy::LinearOperator
             */
    auto hessian(const ::Spacy::Vector& x) const
    {
        assembleHessian(x);
        return Linearization{ A_ , *operatorSpace_ , solverCreator_ };
    }

    /// Access operator representing \f$f''\f$.
    const KaskadeOperator& A() const noexcept
    {
        return A_;
    }

    /// Access boost::fusion::vector of pointers to spaces.
    const Spaces& spaces() const noexcept
    {
        return spaces_;
    }

    /**
             * @brief Access onlyLowerTriangle flag.
             * @return true if only the lower triangle of a symmetric matrix is stored in the operator definition, else false
             */
    bool onlyLowerTriangle() const noexcept
    {
        return onlyLowerTriangle_;
    }

    /**
             * @brief Change solver creator.
             * @param f function/functor for the creation of a linear solver
             */
    void setSolverCreator(std::function<LinearSolver(const Linearization&)> f)
    {
        solverCreator_ = std::move(f);
    }

    auto blockHessian(const ::Spacy::Vector& x) const
    {
        assembleHessian(x);
        return *blockHessianPtr_;
    }

    const Atype& AOp() const
    {
        return AOp_;
    }

    const ATtype& AOp_t() const
    {
        return AOp_t_;
    }


    std::function<::Spacy::Vector(const ::Spacy::Vector&)> MyDiagSolver() const
    {

        return MySolver_;
    }

    std::function<::Spacy::Vector(const ::Spacy::Vector&)> MyDirectSolver() const
    {

        return MyDirectSolver_;
    }

    std::function<::Spacy::Vector(const ::Spacy::Vector&)> controlSolver() const
    {
        //return MuDiagMatrixSolver_;
        return controlSolver_;
    }

    std::function<::Spacy::Vector(const ::Spacy::Vector&)> BPXA() const
    {
        return BPXA_;
    }


    std::function<::Spacy::Vector(const ::Spacy::Vector&)> getAPrimePrime() const
    {
        return APrimePrime;
    }

    std::function<::Spacy::Vector(const ::Spacy::Vector&)> getADualDual() const
    {
        return ADualDual;
    }




    std::function<::Spacy::Vector(const ::Spacy::Vector&)> BPXAT() const
    {
        return BPXAT_;
    }

    std::function<void(int,char,char,int,double*,double*,double,double,int,int,double,int*,double*,double*,int,int*)> LAPACKE() const
    {
        return LAPACKE_dstevr_func;
    }

    std::function<Spacy::Real(const Spacy::Vector&)> NormY() const
    {
        return normY_;
    }

    std::function<Spacy::Real(const Spacy::Vector&)> NormP() const
    {
        return normP_;
    }

    std::function<::Spacy::Vector(const ::Spacy::Vector&)> MyBlock_ = [this](const ::Spacy::Vector& y) mutable
    {
        CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
        CoefficientVectorY y_coeff(boost::fusion::at_c<0>(temp.data));
        CoefficientVectorY My_coeff(boost::fusion::at_c<0>(temp.data));

        ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(y,y_coeff);
        //std::cout << "  My: y_coeff.two_norm() = " << y_coeff.two_norm() << std::endl;

        My_.apply(y_coeff,My_coeff);

        const auto& domain_ps = cast_ref<PSVC>(this->domain().creator());
        const auto& Y = (cast_ref<PSVC>(domain_ps.subSpace(0).creator())).subSpace(0);

        ::Spacy::Vector result = zero(Y);
        ::Spacy::Kaskade::copyFromCoefficientVector<VYSetDescription>(My_coeff,result);
        //std::cout << "  My: get(norm(result)) = " << get(norm(result)) << std::endl;

        return result;
    };

    std::function<::Spacy::Vector(const ::Spacy::Vector&)> MyuBlock_ = [this](const ::Spacy::Vector& u) mutable
    {
        CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
        CoefficientVectorU u_coeff(boost::fusion::at_c<1>(temp.data));
        CoefficientVectorY Myu_coeff(boost::fusion::at_c<0>(temp.data));

        ::Spacy::Kaskade::copyToCoefficientVector<VUSetDescription>(u,u_coeff);
        //std::cout << "  Myu: u_coeff.two_norm() = " << u_coeff.two_norm() << std::endl;

        //Im Normalenschritt ist Myu_ leer -> gib einfach zero(...) zurück?
        if(Myu_.getTriplet().nrows() == 0)
        {
            const auto& domain_ps = cast_ref<PSVC>(this->domain().creator());
            const auto& Y = (cast_ref<PSVC>(domain_ps.subSpace(0).creator())).subSpace(0);
            return zero(Y);
        }
        Myu_.apply(u_coeff,Myu_coeff);

        const auto& domain_ps = cast_ref<PSVC>(this->domain().creator());
        const auto& Y = (cast_ref<PSVC>(domain_ps.subSpace(0).creator())).subSpace(0);

        ::Spacy::Vector result = zero(Y);
        ::Spacy::Kaskade::copyFromCoefficientVector<VYSetDescription>(Myu_coeff,result);
        //std::cout << "  Myu: get(norm(result)) = " << get(norm(result)) << std::endl;

        return result;
    };

    std::function<::Spacy::Vector(const ::Spacy::Vector&)> MuBlock_ = [this](const ::Spacy::Vector& u) mutable
    {
        CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
        CoefficientVectorU u_coeff(boost::fusion::at_c<1>(temp.data));
        CoefficientVectorU Mu_coeff(boost::fusion::at_c<1>(temp.data));

        ::Spacy::Kaskade::copyToCoefficientVector<VUSetDescription>(u,u_coeff);
        //std::cout << "  Mu: u_coeff.two_norm() = " << u_coeff.two_norm() << std::endl;

        Mu_.apply(u_coeff,Mu_coeff);

        const auto& domain_ps = cast_ref<PSVC>(this->domain().creator());
        const auto& U = (cast_ref<PSVC>(domain_ps.subSpace(0).creator())).subSpace(1);

        ::Spacy::Vector result = zero(U);
        ::Spacy::Kaskade::copyFromCoefficientVector<VUSetDescription>(Mu_coeff,result);
        //std::cout << "  Mu: get(norm(result)) = " << get(norm(result)) << std::endl;

        return result;
    };

    std::function<::Spacy::Vector(const ::Spacy::Vector&)> MuyBlock_ = [this](const ::Spacy::Vector& y) mutable
    {
        CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
        CoefficientVectorY y_coeff(boost::fusion::at_c<0>(temp.data));
        CoefficientVectorU Muy_coeff(boost::fusion::at_c<1>(temp.data));

        ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(y,y_coeff);
        //std::cout << "  Muy: y_coeff.two_norm() = " << y_coeff.two_norm() << std::endl;

        //Im Normalenschritt ist Muy_ leer -> gib einfach zero(...) zurück?
        if(Muy_.getTriplet().nrows() == 0)
        {
            const auto& domain_ps = cast_ref<PSVC>(this->domain().creator());
            const auto& U = (cast_ref<PSVC>(domain_ps.subSpace(0).creator())).subSpace(1);
            return zero(U);
        }
        Muy_.apply(y_coeff,Muy_coeff);

        const auto& domain_ps = cast_ref<PSVC>(this->domain().creator());
        const auto& U = (cast_ref<PSVC>(domain_ps.subSpace(0).creator())).subSpace(1);

        ::Spacy::Vector result = zero(U);
        ::Spacy::Kaskade::copyFromCoefficientVector<VUSetDescription>(Muy_coeff,result);
        //std::cout << "  Muy: get(norm(result)) = " << get(norm(result)) << std::endl;

        return result;
    };

    std::function<::Spacy::Vector(const ::Spacy::Vector&)> ATBlock_ = [this](const ::Spacy::Vector& p) mutable
    {
        CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
        CoefficientVectorP p_coeff(boost::fusion::at_c<2>(temp.data));
        CoefficientVectorY ATp_coeff(boost::fusion::at_c<0>(temp.data));

        ::Spacy::Kaskade::copyToCoefficientVector<VPSetDescription>(p,p_coeff);
        //std::cout << "  AT: p_coeff.two_norm() = " << p_coeff.two_norm() << std::endl;

        AOp_t_.apply(p_coeff,ATp_coeff);

        const auto& domain_ps = cast_ref<PSVC>(this->domain().creator());
        const auto& Y = (cast_ref<PSVC>(domain_ps.subSpace(0).creator())).subSpace(0);

        ::Spacy::Vector result = zero(Y);
        ::Spacy::Kaskade::copyFromCoefficientVector<VYSetDescription>(ATp_coeff,result);
        //std::cout << "  AT: get(norm(result)) = " << get(norm(result)) << std::endl;

        return result;
    };

    std::function<::Spacy::Vector(const ::Spacy::Vector&)> ABlock_ = [this](const ::Spacy::Vector& y) mutable
    {
        CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
        CoefficientVectorY y_coeff(boost::fusion::at_c<0>(temp.data));
        CoefficientVectorP Ay_coeff(boost::fusion::at_c<2>(temp.data));

        ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(y,y_coeff);
        //std::cout << "  A: y_coeff.two_norm() = " << y_coeff.two_norm() << std::endl;

        AOp_.apply(y_coeff,Ay_coeff);

        const auto& domain_ps = cast_ref<PSVC>(this->domain().creator());
        const auto& P = (cast_ref<PSVC>(domain_ps.subSpace(1).creator())).subSpace(0);

        ::Spacy::Vector result = zero(P);
        ::Spacy::Kaskade::copyFromCoefficientVector<VPSetDescription>(Ay_coeff,result);
        //std::cout << "  A: get(norm(result)) = " << get(norm(result)) << std::endl;

        return result;
    };

    std::function<::Spacy::Vector(const ::Spacy::Vector&)> APrimePrime = [this](const ::Spacy::Vector& y) mutable
    {
       // std::cout << "TEST: " << std::endl;
        CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
        CoefficientVectorY y_coeff(boost::fusion::at_c<0>(temp.data));
        CoefficientVectorY Ay_coeff(boost::fusion::at_c<0>(temp.data));

        ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(y,y_coeff);
        //std::cout << "  A: y_coeff.two_norm() = " << y_coeff.two_norm() << std::endl;

        AOp_.apply(y_coeff,Ay_coeff);

        const auto& domain_ps = cast_ref<PSVC>(this->domain().creator());
        const auto& P = (cast_ref<PSVC>(domain_ps.subSpace(0).creator())).subSpace(0);

        ::Spacy::Vector result = zero(P);
        ::Spacy::Kaskade::copyFromCoefficientVector<VYSetDescription>(Ay_coeff,result);
        //std::cout << "  A: get(norm(result)) = " << get(norm(result)) << std::endl;

        return result;
    };

    std::function<::Spacy::Vector(const ::Spacy::Vector&)> ADualDual = [this](const ::Spacy::Vector& y) mutable
    {
        //std::cout << "ADualDual: " << y(y)<< std::endl;
        CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
        CoefficientVectorP y_coeff(boost::fusion::at_c<2>(temp.data));
        CoefficientVectorP Ay_coeff(boost::fusion::at_c<2>(temp.data));

        ::Spacy::Kaskade::copyToCoefficientVector<VPSetDescription>(y,y_coeff);
        //std::cout << "ANormAfterCopy: " << y_coeff.two_norm() << std::endl;

        AOp_.apply(y_coeff,Ay_coeff);

        const auto& domain_ps = cast_ref<PSVC>(this->domain().creator());
        const auto& P = (cast_ref<PSVC>(domain_ps.subSpace(1).creator())).subSpace(0);

        ::Spacy::Vector result = zero(P);
        ::Spacy::Kaskade::copyFromCoefficientVector<VPSetDescription>(Ay_coeff,result);
        //std::cout << "  A: get(norm(result)) = " << get(norm(result)) << std::endl;

        return result;
    };

    std::function<::Spacy::Vector(const ::Spacy::Vector&)> BBlock_ = [this](const ::Spacy::Vector& u) mutable
    {
        CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
        CoefficientVectorU u_coeff(boost::fusion::at_c<1>(temp.data));
        CoefficientVectorP Bu_coeff(boost::fusion::at_c<2>(temp.data));

        ::Spacy::Kaskade::copyToCoefficientVector<VUSetDescription>(u,u_coeff);
        // std::cout << "  B: u_coeff.two_norm() = " << Bu_coeff.two_norm() << std::endl;


        auto startTime = std::chrono::high_resolution_clock::now();


        unsigned int uVectorSize = boost::fusion::at_c<0>(u_coeff.data).size();
        unsigned int uBlockSize =  boost::fusion::at_c<0>(u_coeff.data).operator[](0).size();

        unsigned int BuVectorSize = boost::fusion::at_c<0>(Bu_coeff.data).size();
        unsigned int BuBlockSize =  boost::fusion::at_c<0>(Bu_coeff.data).operator[](0).size();


        unsigned int uSize =  uVectorSize * uBlockSize;
        unsigned int BuSize =  BuVectorSize * BuBlockSize;

       // std::cout << "uSize: " << uSize << std::endl;
       // std::cout << "BuSize: " << BuSize << std::endl;

        std::vector<double> u_(uSize,0.0);
        std::vector<double> Bu_(BuSize,0.0);

        auto & matrix = B_.get_non_const();

        auto counter = 0u;

        for(auto & it: boost::fusion::at_c<0>(u_coeff.data))
            for(auto & iter: it)
            {
                u_[counter] = iter;
                counter++;
            }

        for(size_t i=0; i< matrix.ridx.size(); ++i)
             Bu_[matrix.ridx[i]] += u_[matrix.cidx[i]]*matrix.data[i];

        counter = 0;

        for(auto & it: boost::fusion::at_c<0>(Bu_coeff.data))
            for(auto & iter: it)
            {
                iter = Bu_[counter];
                counter++;
            }


        //            int size = B_.getTriplet().ridx.size();
        //            int blockSize = boost::fusion::at_c<0>(u_coeff.data).operator[](0).size();

        //            for(int i = 0; i < size; i++)
        //            {
        //                int row = B_.getTriplet().ridx[i];
        //                int col = B_.getTriplet().cidx[i];

        //                int index1 = col / blockSize;
        //                int index2 = col % blockSize;

        //                int index3 = row /blockSize;
        //                int index4 = row % blockSize;

        //                (boost::fusion::at_c<0>(Bu_coeff.data).operator[](index3)).operator[](index4) +=  (boost::fusion::at_c<0>(u_coeff.data).operator[](index1)).operator[](index2)*B_.getTriplet().data[i];
        //            }

        //            auto max = 0;

        //            for(int i = 0; i< B_.getTriplet().cidx.size(); ++i)
        //            {
        //                max = std::max(max, B_.getTriplet().cidx[i]);
        //            }

        //            std::cout << "B columns: " << max << std::endl;

        //            max = 0;

        //            for(int i= 0; i<B_.getTriplet().ridx.size(); ++i)
        //            {
        //                max = std::max(max, B_.getTriplet().ridx[i]);
        //            }

        //            std::cout << "B rows: " << max << std::endl;

      //  B_.apply(u_coeff,Bu_coeff);

        //  std::cout << "Test: " << boost::fusion::at_c<0>(u_coeff.data).size() << std::endl;
        //  std::cout << "Test: " << boost::fusion::at_c<0>(Bu_coeff.data).size() << std::endl;

      //  auto time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - startTime).count();

      //  std::cout << "Elapsed Time B: " << time <<  "ms"  << std::endl;
      //  std::cout << "Fail: " << std::endl;

        const auto& domain_ps = cast_ref<PSVC>(this->domain().creator());
        const auto& P = (cast_ref<PSVC>(domain_ps.subSpace(1).creator())).subSpace(0);

        ::Spacy::Vector result = zero(P);
        ::Spacy::Kaskade::copyFromCoefficientVector<VPSetDescription>(Bu_coeff,result);
        //std::cout << "  B: get(norm(result)) = " << get(norm(result)) << std::endl;

        return result;
    };

    std::function<::Spacy::Vector(const ::Spacy::Vector&)> BTBlock_ = [this](const ::Spacy::Vector& p) mutable
    {
        CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
        CoefficientVectorP p_coeff(boost::fusion::at_c<2>(temp.data));
        CoefficientVectorU BTp_coeff(boost::fusion::at_c<1>(temp.data));

        ::Spacy::Kaskade::copyToCoefficientVector<VPSetDescription>(p,p_coeff);
        //std::cout << "  BT: p_coeff.two_norm() = " << p_coeff.two_norm() << std::endl;

        // B_.applyscaleaddTransposed(1.0, p_coeff, BTp_coeff); // b1_ -> b1_ +(-1.0)*(-B^T xp_)


        auto startTime = std::chrono::high_resolution_clock::now();


        unsigned int pVectorSize = boost::fusion::at_c<0>(p_coeff.data).size();
        unsigned int pBlockSize =  boost::fusion::at_c<0>(p_coeff.data).operator[](0).size();

        unsigned int BTpVectorSize = boost::fusion::at_c<0>(BTp_coeff.data).size();
        unsigned int BTpBlockSize =  boost::fusion::at_c<0>(BTp_coeff.data).operator[](0).size();


        unsigned int pSize =   pVectorSize * pBlockSize;
        unsigned int BTpSize =  BTpVectorSize * BTpBlockSize;

       // std::cout << "pSize: " << pSize << std::endl;
       // std::cout << "BTpSize: " << BTpSize << std::endl;

        std::vector<double> p_(pSize,0.0);
        std::vector<double> BTp_(BTpSize,0.0);

        auto & matrix = B_t_.get_non_const();

        auto counter = 0u;

        for(auto & it: boost::fusion::at_c<0>(p_coeff.data))
            for(auto & iter: it)
            {
                p_[counter] = iter;
                counter++;
            }

        for(size_t i=0; i< matrix.ridx.size(); ++i)
             BTp_[matrix.ridx[i]] += p_[matrix.cidx[i]]*matrix.data[i];

        counter = 0;

        for(auto & it: boost::fusion::at_c<0>(BTp_coeff.data))
            for(auto & iter: it)
            {
                iter = BTp_[counter] ;
                counter++;
            }


        //            int size = B_.getTriplet().ridx.size();
        //            int blockSize = boost::fusion::at_c<0>(u_coeff.data).operator[](0).size();

        //            for(int i = 0; i < size; i++)
        //            {
        //                int row = B_.getTriplet().ridx[i];
        //                int col = B_.getTriplet().cidx[i];

        //                int index1 = col / blockSize;
        //                int index2 = col % blockSize;

        //                int index3 = row /blockSize;
        //                int index4 = row % blockSize;

        //                (boost::fusion::at_c<0>(Bu_coeff.data).operator[](index3)).operator[](index4) +=  (boost::fusion::at_c<0>(u_coeff.data).operator[](index1)).operator[](index2)*B_.getTriplet().data[i];
        //            }

        //            auto max = 0;

        //            for(int i = 0; i< B_.getTriplet().cidx.size(); ++i)
        //            {
        //                max = std::max(max, B_.getTriplet().cidx[i]);
        //            }

        //            std::cout << "B columns: " << max << std::endl;

        //            max = 0;

        //            for(int i= 0; i<B_.getTriplet().ridx.size(); ++i)
        //            {
        //                max = std::max(max, B_.getTriplet().ridx[i]);
        //            }

        //            std::cout << "B rows: " << max << std::endl;

      //  B_t_.apply(u_coeff,Bu_coeff);

        //  std::cout << "Test: " << boost::fusion::at_c<0>(u_coeff.data).size() << std::endl;
        //  std::cout << "Test: " << boost::fusion::at_c<0>(Bu_coeff.data).size() << std::endl;

        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - startTime).count();

      //  std::cout << "Elapsed Time B: " << time <<  "ms"  << std::endl;
      //  std::cout << "copy " << std::endl;

        const auto& domain_ps = cast_ref<PSVC>(this->domain().creator());
        const auto& U = (cast_ref<PSVC>(domain_ps.subSpace(0).creator())).subSpace(1);

        ::Spacy::Vector result = zero(U);
        ::Spacy::Kaskade::copyFromCoefficientVector<VUSetDescription>(BTp_coeff,result);
        //std::cout << "  BT: get(norm(result)) = " << get(norm(result)) << std::endl;

        return result;
    };

    std::function<::Spacy::Vector(const ::Spacy::Vector&)> BPXAPrim_ = [this](const ::Spacy::Vector& p) mutable
    {

        CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
        CoefficientVectorY p_coeff(boost::fusion::at_c<0>(temp.data));
        CoefficientVectorY y(boost::fusion::at_c<0>(temp.data));

        ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(p,p_coeff);

      //  std::cout << "  BPX: p_coeff.two_norm() = " << p_coeff.two_norm() << std::endl;
        BPXAPtr_->apply(boost::fusion::at_c<0>(y.data),boost::fusion::at_c<0>(p_coeff.data));
      //  std::cout << "  BPX: y_coeff.two_norm() = " << y.two_norm() << std::endl;


        const auto& domain_ps = ::Spacy::cast_ref<PSVC>(this->domain().creator());
        const auto& Y = (::Spacy::cast_ref<PSVC>(domain_ps.subSpace(0).creator())).subSpace(0);

        ::Spacy::Vector result = zero(Y);
        ::Spacy::Kaskade::copyFromCoefficientVector<VYSetDescription>(y,result);
        //std::cout << "  BPX: get(norm(result)) = " << get(norm(result)) << std::endl;

        return result;
    };

    std::function<::Spacy::Vector(const ::Spacy::Vector&)> BPXADual_ = [this](const ::Spacy::Vector& y) mutable
    {

        //std::cout << "TestNormBeforeDualBPX: " <<  sqrt(y(y)) << std::endl;
        CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
        CoefficientVectorP y_coeff(boost::fusion::at_c<2>(temp.data));
        CoefficientVectorP p(boost::fusion::at_c<2>(temp.data));

        ::Spacy::Kaskade::copyToCoefficientVector<VPSetDescription>(y,y_coeff);
  //      std::cout << "  BPXT: y_coeff.two_norm() = " << y_coeff.two_norm() << std::endl;

      //  std::cout << " TestNormBeforeDualBPX After copy " << y_coeff.two_norm() << std::endl;
        BPXATPtr_->apply(boost::fusion::at_c<0>(p.data),boost::fusion::at_c<0>(y_coeff.data));

    //   std::cout << "  BPXT: p.two_norm() = " << p.two_norm() << std::endl;



        const auto& domain_ps = ::Spacy::cast_ref<PSVC>(this->domain().creator());
        const auto& P = (::Spacy::cast_ref<PSVC>(domain_ps.subSpace(1).creator())).subSpace(0);

        ::Spacy::Vector result = zero(P);
        ::Spacy::Kaskade::copyFromCoefficientVector<VPSetDescription>(p,result);
       // std::cout << "  BPXDual: get(norm(result)) = " << get(norm(result)) << std::endl;

        return result;
    };

    std::function<::Spacy::Vector(const ::Spacy::Vector&)> BPXAPrimDual_ = [this](const ::Spacy::Vector& y) mutable
    {
      //  std::cout << "TestNormBeforeDualBPX: " <<  sqrt(y(y)) << std::endl;
        CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
        CoefficientVectorY y_coeff(boost::fusion::at_c<0>(temp.data));
        CoefficientVectorP p(boost::fusion::at_c<2>(temp.data));

        ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(y,y_coeff);
  //      std::cout << "  BPXT: y_coeff.two_norm() = " << y_coeff.two_norm() << std::endl;

     //   std::cout << " TestNormBeforeDualBPX After copy " << y_coeff.two_norm() << std::endl;
        BPXATPtr_->apply(boost::fusion::at_c<0>(p.data),boost::fusion::at_c<0>(y_coeff.data));

    //   std::cout << "  BPXT: p.two_norm() = " << p.two_norm() << std::endl;



        const auto& domain_ps = ::Spacy::cast_ref<PSVC>(this->domain().creator());
        const auto& P = (::Spacy::cast_ref<PSVC>(domain_ps.subSpace(1).creator())).subSpace(0);

        ::Spacy::Vector result = zero(P);
        ::Spacy::Kaskade::copyFromCoefficientVector<VPSetDescription>(p,result);
        //std::cout << "  BPXT: get(norm(result)) = " << get(norm(result)) << std::endl;

        return result;
    };



    std::function<::Spacy::Vector(const ::Spacy::Vector&)> BPXA_ = [this](const ::Spacy::Vector& p) mutable
    {
        CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
        CoefficientVectorP p_coeff(boost::fusion::at_c<2>(temp.data));
        CoefficientVectorY y(boost::fusion::at_c<0>(temp.data));

        ::Spacy::Kaskade::copyToCoefficientVector<VPSetDescription>(p,p_coeff);

      //  std::cout << "  BPX: p_coeff.two_norm() = " << p_coeff.two_norm() << std::endl;
        BPXAPtr_->apply(boost::fusion::at_c<0>(y.data),boost::fusion::at_c<0>(p_coeff.data));
      //  std::cout << "  BPX: y_coeff.two_norm() = " << y.two_norm() << std::endl;


        const auto& domain_ps = ::Spacy::cast_ref<PSVC>(this->domain().creator());
        const auto& Y = (::Spacy::cast_ref<PSVC>(domain_ps.subSpace(0).creator())).subSpace(0);

        ::Spacy::Vector result = zero(Y);
        ::Spacy::Kaskade::copyFromCoefficientVector<VYSetDescription>(y,result);
        //std::cout << "  BPX: get(norm(result)) = " << get(norm(result)) << std::endl;

        return result;
    };

    std::function<::Spacy::Vector(const ::Spacy::Vector&)> BPXAT_ = [this](const ::Spacy::Vector& y) mutable
    {
        CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
        CoefficientVectorY y_coeff(boost::fusion::at_c<0>(temp.data));
        CoefficientVectorP p(boost::fusion::at_c<2>(temp.data));

        ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(y,y_coeff);
       // std::cout << "  BPXT: y_coeff.two_norm() = " << y_coeff.two_norm() << std::endl;

        BPXATPtr_->apply(boost::fusion::at_c<0>(p.data),boost::fusion::at_c<0>(y_coeff.data));

       // std::cout << "  BPXT: p.two_norm() = " << p.two_norm() << std::endl;



        const auto& domain_ps = ::Spacy::cast_ref<PSVC>(this->domain().creator());
        const auto& P = (::Spacy::cast_ref<PSVC>(domain_ps.subSpace(1).creator())).subSpace(0);

        ::Spacy::Vector result = zero(P);
        ::Spacy::Kaskade::copyFromCoefficientVector<VPSetDescription>(p,result);
        //std::cout << "  BPXT: get(norm(result)) = " << get(norm(result)) << std::endl;

        return result;
    };






    std::function<CoefficientVectorU(const CoefficientVectorU&)> MuDiagMatrixSolver_ = [this](const CoefficientVectorU & u)
    {
        CoefficientVectorU result(u);

        unsigned int counter = 0;

        auto resultIter =  boost::fusion::at_c<0>(result.data).begin();

        for(const auto & it : boost::fusion::at_c<0>(u.data))
        {
            for(auto i = 0u; i < it.size(); i++)
            {
                resultIter->operator[](i) = it.operator[](i)/MuDiag_[counter];
                counter++;
            }

            resultIter++;
        }

        return result;
    };

    //std::shared_ptr<::Kaskade::InverseLinearOperator<::Kaskade::DirectSolver<CoefficientVectorU, CoefficientVectorU> > > solMu_ =
    //std::make_shared<::Kaskade::InverseLinearOperator<::Kaskade::DirectSolver<CoefficientVectorU, CoefficientVectorU> > >(::Kaskade::directInverseOperator(Mu_, DirectType::UMFPACK3264, MatrixProperties::GENERAL));

    std::function<::Spacy::Vector(const ::Spacy::Vector&)> MySolver_ = [this](const ::Spacy::Vector& ry)
    {


        CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(this->spaces()));
        CoefficientVectorY ry_coeff(boost::fusion::at_c<0>(temp.data));
        ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(ry,ry_coeff);

        CoefficientVectorY result_coeff(boost::fusion::at_c<0>(temp.data));

        unsigned int counter = 0;

        auto resultIter =  boost::fusion::at_c<0>(result_coeff.data).begin();



        for(const auto & it : boost::fusion::at_c<0>(ry_coeff.data))
        {
            for(auto i = 0u; i < it.size(); i++)
            {

               // std::cout << "resultIter->operator[](i) : " << resultIter->operator[](i) << std::endl;
                resultIter->operator[](i) = it.operator[](i)/MyDiag_[counter];
               // std::cout << "Counter  MyDiag_[counter]: " << it.operator[](i)/MyDiag_[counter] << std::endl;
                counter++;

            }

            resultIter++;
        }

        ::Spacy::Vector result = ry;

        ::Spacy::Kaskade::copyFromCoefficientVector<VYSetDescription>(result_coeff,result);

       // std::cout << "Return: " << std::endl;
        return ry;
    };



    std::function<::Spacy::Vector(const ::Spacy::Vector&)> controlSolver_ = [this](const ::Spacy::Vector& u)
    {
        CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(this->spaces()));
        CoefficientVectorU u_coeff(boost::fusion::at_c<1>(temp.data));
        ::Spacy::Kaskade::copyToCoefficientVector<VUSetDescription>(u,u_coeff);

        CoefficientVectorU result_coeff(boost::fusion::at_c<1>(temp.data));
        if(useDirectMu)
            solMu_->apply(u_coeff, result_coeff);
        else
            result_coeff = MuDiagMatrixSolver_(u_coeff);

        ::Spacy::Vector result = u;
        ::Spacy::Kaskade::copyFromCoefficientVector<VUSetDescription>(result_coeff,result);

        return result;
    };

public:

    std::function<::Spacy::Vector(const ::Spacy::Vector&)> ASolver_ = [this](const ::Spacy::Vector& p)
    {
        CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(this->spaces()));
        CoefficientVectorY p_coeff(boost::fusion::at_c<0>(temp.data));
        ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(p,p_coeff);

        CoefficientVectorY result_coeff(boost::fusion::at_c<0>(temp.data));

            solA_->apply(p_coeff, result_coeff);

        ::Spacy::Vector result = p;
        ::Spacy::Kaskade::copyFromCoefficientVector<VYSetDescription>(result_coeff,result);

        return result;
    };

    std::function<::Spacy::Vector(const ::Spacy::Vector&)> ATSolver_ = [this](const ::Spacy::Vector& y)
    {

        CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(this->spaces()));
        CoefficientVectorP y_coeff(boost::fusion::at_c<2>(temp.data));
        ::Spacy::Kaskade::copyToCoefficientVector<VPSetDescription>(y,y_coeff);

        CoefficientVectorP result_coeff(boost::fusion::at_c<2>(temp.data));

            solA_->apply(y_coeff, result_coeff);

        ::Spacy::Vector result = y;
        ::Spacy::Kaskade::copyFromCoefficientVector<VPSetDescription>(result_coeff,result);

        return result;
    };

    std::function<::Spacy::Vector(const ::Spacy::Vector&)> MyDirectSolver_ = [this](const ::Spacy::Vector& y)
    {

        CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(this->spaces()));
        CoefficientVectorY y_coeff(boost::fusion::at_c<0>(temp.data));
        ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(y,y_coeff);

        CoefficientVectorY result_coeff(boost::fusion::at_c<0>(temp.data));

            solMy_->apply(y_coeff, result_coeff);

        ::Spacy::Vector result = y;
        ::Spacy::Kaskade::copyFromCoefficientVector<VYSetDescription>(result_coeff,result);

        return result;
    };


    std::function<::Spacy::Vector(const ::Spacy::Vector&)> diagMuOperator_ = [this](const ::Spacy::Vector& u)
    {
        CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(this->spaces()));
        CoefficientVectorU u_coeff(boost::fusion::at_c<1>(temp.data));
        ::Spacy::Kaskade::copyToCoefficientVector<VUSetDescription>(u,u_coeff);


        //  std::cout << "Test DiagMuOperator: " << std::endl;
        CoefficientVectorU resultCoeff(boost::fusion::at_c<1>(temp.data));

        unsigned int counter = 0;

        auto resultIter =  boost::fusion::at_c<0>(resultCoeff.data).begin();

        for(const auto & it : boost::fusion::at_c<0>(u_coeff.data))
        {
            for(auto i = 0u; i < it.size(); i++)
            {
                resultIter->operator[](i) = it.operator[](i) *MuDiag_[counter];
                counter++;
            }

            resultIter++;
        }

        ::Spacy::Vector result = u;
        ::Spacy::Kaskade::copyFromCoefficientVector<VUSetDescription>(resultCoeff,result);

        return result;
    };




    std::function<void(int,char,char,int,double*,double*,double,double,int,int,double,int*,double*,double*,int,int*)> LAPACKE_dstevr_func = [](int matrix_layout,char jobz,char range,int n,double* d,double* e,double vl,double vu,int il,int iu,double abstol,int* m,double* eig,double* z,int ldz,int* isuppz)
    {
        LAPACKE_dstevr(matrix_layout,jobz,range,n,d,e,vl,vu,il,iu,abstol,m,eig,z,ldz,isuppz);
    };

    std::function<Spacy::Real(const Spacy::Vector&)> normY_ = [this](const Spacy::Vector& y)
    {
        Spacy::Real max = -1;
        CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
        CoefficientVectorY y_coeff(boost::fusion::at_c<0>(temp.data));

        ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(y,y_coeff);

        for (auto j =0u; j < boost::fusion::at_c<0>(y_coeff.data).size(); j++)
        {
            //1 Schleife reicht falls stateDim = 1
            if(fabs(boost::fusion::at_c<0>(y_coeff.data)[j][0]) > max)
                max = fabs(boost::fusion::at_c<0>(y_coeff.data)[j][0]);
            ;
        }
        return max;
    };

    std::function<Spacy::Real(const Spacy::Vector&)> normP_ = [this](const Spacy::Vector& p)
    {
        Spacy::Real max = -1;
        CoefficientVector temp(VariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces_));
        CoefficientVectorP p_coeff(boost::fusion::at_c<2>(temp.data));

        ::Spacy::Kaskade::copyToCoefficientVector<VPSetDescription>(p,p_coeff);

        for (auto j =0u; j < boost::fusion::at_c<0>(p_coeff.data).size(); j++)
        {
            //1 Schleife reicht falls stateDim = 1
            if(fabs(boost::fusion::at_c<0>(p_coeff.data)[j][0]) > max)
                max = fabs(boost::fusion::at_c<0>(p_coeff.data)[j][0]);
            ;
        }
        return max;
    };

private:
    /// Assemble \f$f(x)\f$.
    void assembleFunctional(const ::Spacy::Vector& x) const
    {
        if( old_X_f_ && old_X_f_ == x && !reAssembleFValue )
            return;

        VariableSetDescription variableSet(spaces_);
        typename VariableSetDescription::VariableSet u(variableSet);

        copy(x,u);

        Assembler assembler(spaces_);
        // std::cout << "GetNumberOfThreads: " << getNumberOfThreads() << std::endl;
        assembler.assemble(::Kaskade::linearization(f_,u) , Assembler::VALUE , getNumberOfThreads() );
        value_ = assembler.functional();

        old_X_f_ = x;
        reAssembleFValue = false;
    }

    /// Assemble discrete representation of \f$f'(x)\f$.
    void assembleGradient(const ::Spacy::Vector& x) const
    {
        using boost::fusion::at_c;
        if( old_X_df_ && old_X_df_ == x && !reAssembleFd1 ) return;

        VariableSetDescription variableSet(spaces_);
        typename VariableSetDescription::VariableSet u(variableSet);

        copy(x,u);

        //std::cout << "GetNumberOfThreads: " << getNumberOfThreads() << std::endl;
        Assembler assembler(spaces_);
        assembler.assemble(::Kaskade::linearization(f_,u) , Assembler::RHS , getNumberOfThreads() );
        copyFromCoefficientVector<VariableSetDescription>(CoefficientVector(assembler.rhs()),rhs_);

        old_X_df_ = x;
        reAssembleFd1 = false;
    }
public:
    /// Assemble discrete representation of \f$f''(x)\f$.
    void assembleHessian(const ::Spacy::Vector& x) const
    {
        if( old_X_ddf_ && (old_X_ddf_==x) && !reAssembleFd2) return;

        VariableSetDescription variableSet(spaces_);
        typename VariableSetDescription::VariableSet u(variableSet);

        copy(x,u);

        Assembler assembler(spaces_);
        //std::cout << "GetNumberOfThreads: " << getNumberOfThreads() << std::endl;
        assembler.assemble(::Kaskade::linearization(f_,u) , Assembler::MATRIX , getNumberOfThreads() );
        A_ = KaskadeOperator( assembler.template get<Matrix>(onlyLowerTriangle_,rbegin_,rend_,cbegin_,cend_) );

        //NEU
        My_ = Mytype( assembler.template get<Matrix>(onlyLowerTriangle_,0,1,0,1) );
        Myu_ = Myutype( assembler.template get<Matrix>(onlyLowerTriangle_,0,1,1,2) );
        AOp_t_ = ATtype( assembler.template get<Matrix>(onlyLowerTriangle_,0,1,2,3) );
        Muy_ = Muytype( assembler.template get<Matrix>(onlyLowerTriangle_,1,2,0,1) );
        Mu_ = Mutype( assembler.template get<Matrix>(onlyLowerTriangle_,1,2,1,2) );

        B_t_ = BTtype( assembler.template get<Matrix>(onlyLowerTriangle_,1,2,2,3) );
        B_ = Btype( assembler.template get<Matrix>(onlyLowerTriangle_,2,3,1,2) );
        AOp_ = Atype( assembler.template get<Matrix>(onlyLowerTriangle_,2,3,0,1) );

       // auto stateSize = AOp_.get().N();
       // auto controlSize = Mu_.get().N();

    //    std::cout << "Test My " << My_.get().nrows() << " " << My_.get().ncols() << std::endl;
    //    std::cout << "Test My " << Mu_.get().nrows() << " " << Mu_.get().ncols() << std::endl;
   //     std::cout << "Test AOp_ " << AOp_.get().nrows() << " " << AOp_.get().ncols() << std::endl;
   //     std::cout << "Test AOp_t_ " << AOp_t_.get().nrows() << " " << AOp_t_.get().ncols() << std::endl;



      //  B_ = getBBlock<Btype>(A_, stateSize, controlSize);
      //  B_t_ = getBTBlock<Btype>(A_, stateSize, controlSize);

//        std::cout << std::endl;

//        auto max = 0;

//        for(int i=0; i<B_.getTriplet().cidx.size(); ++i)
//        {
//            max = std::max(max, B_.getTriplet().cidx[i]);
//        }

//        std::cout << "B columns max Index: " << max << std::endl;
//        std::cout << "B col max Index Orig: " << max+stateSize << std::endl;
//        std::cout << "B getColumns(): " << B_.getTriplet().ncols() << std::endl;
//        std::cout << std::endl;

//        max = 0;

//        for(int i=0; i<B_.getTriplet().ridx.size(); ++i)
//        {
//            max = std::max(max, B_.getTriplet().ridx[i]);
//        }

//        std::cout << "B rows max Index: " << max << std::endl;
//        std::cout << "B rows max Index Orig: " << max+stateSize+controlSize << std::endl;
//        std::cout << "B getRows()" << B_.getTriplet().nrows() << std::endl;
//        std::cout << std::endl;

//        max = 0;

//        for(int i=0; i<B_t_.getTriplet().cidx.size(); ++i)
//        {
//            max = std::max(max, B_t_.getTriplet().cidx[i]);
//        }

//        std::cout << "BT column max Index: " << max << std::endl;
//        std::cout << "Bt col max Index Orig: " << max+stateSize+controlSize << std::endl;
//        std::cout << "Bt getColumns(): " << B_t_.getTriplet().ncols() << std::endl;
//        std::cout << std::endl;
//        max = 0;

//        for(int i=0; i<B_t_.getTriplet().ridx.size(); ++i)
//        {
//            max = std::max(max, B_t_.getTriplet().ridx[i]);
//        }

//        std::cout << "Bt rows max Index: " << max << std::endl;
//        std::cout << "Bt rows max Index Orig: " << max+stateSize << std::endl;
//        std::cout << "Bt getRows(): " << B_t_.getTriplet().nrows() << std::endl;
//        std::cout << std::endl;


//        //  std::cout << "Size: " << B_.getTriplet().ridx.size() << std::endl;


//        std::cout << "B largest index row should be: " << stateSize + stateSize + controlSize -1 << std::endl;
//        std::cout << "B largest index col should be: " << stateSize + controlSize -1 << std::endl;

//        std::cout << "BT largest index row should be: " << stateSize + controlSize -1 << std::endl;
//        std::cout << "BT largest index col should be: " << stateSize + controlSize + stateSize -1 << std::endl;

//        std::cout << std::endl;

        blockHessianPtr_ = std::make_shared<::Spacy::ProductSpace::Operator_V2 >(std::vector<std::vector<CallableOperator>>{{MyBlock_,MyuBlock_,ATBlock_},{MuyBlock_,MuBlock_,BTBlock_},{ABlock_,BBlock_,zeroBlock_}}, this->domain(), this->domain());

        const auto & AOptemp = AOp_.get();
        auto AOpbcrs = AOptemp.template toBCRS<bpxIndex>(); //TODO: statt 1 ...?
        BPXAPtr_ = std::make_shared<BPX>(::Kaskade::makeBPX( ::Kaskade::NumaBCRSMatrix<Dune::FieldMatrix<double, bpxIndex, bpxIndex>>(*AOpbcrs,false), gridManager_));

        const auto & ATOptemp = AOp_t_.get();
        auto  ATOpbcrs = ATOptemp.template toBCRS<bpxIndex>();//TODO: statt 1 ...?
        BPXATPtr_ = std::make_shared<BPX>(::Kaskade::makeBPX( ::Kaskade::NumaBCRSMatrix<Dune::FieldMatrix<double, bpxIndex, bpxIndex>>(*ATOpbcrs,false), gridManager_));


        auto degreesOfFreedomControl = 0;
        degreesOfFreedomControl = u.descriptions.degreesOfFreedom(1, 2);
        MuDiag_.resize(degreesOfFreedomControl);
        for(size_t i = 0; i < Mu_.get().nnz(); i++)
        {
            if(Mu_.get().ridx[i] == Mu_.get().cidx[i]){
                MuDiag_[Mu_.get().ridx[i]] = Mu_.get().data[i];
            }
        }


        auto degreesOfFreedomState = 0;
        degreesOfFreedomState = u.descriptions.degreesOfFreedom(0, 1);
        MyDiag_.resize(degreesOfFreedomState);

        for(size_t i = 0; i < My_.get().nnz(); i++)
        {
            if(My_.get().ridx[i] == My_.get().cidx[i]){
                MyDiag_[My_.get().ridx[i]] = My_.get().data[i];
            }
        }






        solMu_ = std::make_shared
                < ::Kaskade::InverseLinearOperator<
                ::Kaskade::DirectSolver<CoefficientVectorU, CoefficientVectorU> >
                > (::Kaskade::directInverseOperator(Mu_, DirectType::UMFPACK3264,
                                                    MatrixProperties::GENERAL));


        solA_ = std::make_shared
                < ::Kaskade::InverseLinearOperator<
                ::Kaskade::DirectSolver<CoefficientVectorY, CoefficientVectorY> >
                > (::Kaskade::directInverseOperator(AOp_, DirectType::UMFPACK3264,
                                                    MatrixProperties::GENERAL));


        solAT_ = std::make_shared
                < ::Kaskade::InverseLinearOperator<
                ::Kaskade::DirectSolver<CoefficientVectorP, CoefficientVectorP> >
                > (::Kaskade::directInverseOperator(AOp_t_, DirectType::UMFPACK3264,
                                                    MatrixProperties::GENERAL));

        solMy_ = std::make_shared
                < ::Kaskade::InverseLinearOperator<
                ::Kaskade::DirectSolver<CoefficientVectorY, CoefficientVectorY> >
                > (::Kaskade::directInverseOperator(My_, DirectType::UMFPACK3264,
                                                    MatrixProperties::GENERAL));



        old_X_ddf_ = x;
        reAssembleFd2 = false;

        auto & matrixA = AOp_.get_non_const();

        auto & matrixAT = AOp_t_.get_non_const();

        if(matrixA.ridx.size() != matrixAT.ridx.size())
            std::cout << "Fail A is not symmetric size: " << std::endl;

  //      for(int i = 0; i < matrixA.ridx.size(); i++)
    //    {
      //      if(matrixA.data[i] != matrixAT.data[i] /*|| matrixA.ridx[i] != matrixAT.ridx[i] || matrixA.cidx[i] != matrixAT.cidx[i]*/ )
        //    {
          //       std::cout << "Fail A is not symmetric: " << std::endl;
            //    break;
           // }
      //  }



//        std::cout << "test 1 " << std::endl;
//        if(!matrix.isSymmetric())
//            std::cout << " Fail HesseMatrix is not Symmetric: " << std::endl;

//        std::cout << "test 2 " << std::endl;
        /*
        std::cout << "Myu_.cidx: ";
        for(int i=0; i<Myu_.getTriplet().cidx.size(); ++i)
        {
            std::cout << Myu_.getTriplet().cidx[i] << ", ";
        }
        std::cout << std::endl;
        std::cout << "Myu_.ridx: ";
        for(int i=0; i<Myu_.getTriplet().ridx.size(); ++i)
        {
            std::cout << Myu_.getTriplet().ridx[i] << ", ";
        }
        std::cout << std::endl;

        std::cout << "Myu_: ";
        for(int i=0; i<Myu_.getTriplet().data.size(); ++i)
        {
            std::cout << std::setprecision(20) << Myu_.getTriplet().data[i] << ", ";
        }
        std::cout << std::endl;
        */

    }
    mutable std::shared_ptr<::Kaskade::InverseLinearOperator<::Kaskade::DirectSolver<CoefficientVectorY, CoefficientVectorY> > > solMy_ = nullptr;
    mutable std::shared_ptr<::Kaskade::InverseLinearOperator<::Kaskade::DirectSolver<CoefficientVectorY, CoefficientVectorY> > > solA_ = nullptr;
    mutable std::shared_ptr<::Kaskade::InverseLinearOperator<::Kaskade::DirectSolver<CoefficientVectorP, CoefficientVectorP> > > solAT_ = nullptr;


private:
    FunctionalDefinition f_;
    Spaces spaces_;
    const GridManager & gridManager_;
    mutable KaskadeOperator A_ = {};

    mutable Mytype My_{};
    mutable Myutype Myu_{};
    mutable Mutype Mu_{};
    mutable Muytype Muy_{};
    mutable Atype  AOp_{};
    mutable Btype B_{};
    mutable ATtype AOp_t_{};
    mutable BTtype B_t_{};
    mutable std::shared_ptr<::Spacy::ProductSpace::Operator_V2> blockHessianPtr_ = nullptr;
    mutable std::vector<double> MuDiag_ = {};
    mutable std::vector<double> MyDiag_ = {};
    mutable std::shared_ptr<::Kaskade::InverseLinearOperator<::Kaskade::DirectSolver<CoefficientVectorU, CoefficientVectorU> > > solMu_ = nullptr;

    using BPX = decltype(::Kaskade::makeBPX(std::declval<::Kaskade::NumaBCRSMatrix<Dune::FieldMatrix<double, bpxIndex, bpxIndex>>>(),std::declval<GridManager>() ) );
    mutable std::shared_ptr<BPX> BPXAPtr_ = nullptr;
    mutable std::shared_ptr<BPX> BPXATPtr_ = nullptr;


    mutable double value_ = 0;
    mutable ::Spacy::Vector old_X_f_{}, old_X_df_{}, old_X_ddf_{}, rhs_{};
    bool onlyLowerTriangle_ = false;
    int rbegin_=0, rend_=FunctionalDefinition::AnsatzVars::noOfVariables;
    int cbegin_=0, cend_=FunctionalDefinition::TestVars::noOfVariables;

    mutable bool reAssembleFValue = false;
    mutable bool reAssembleFd1 = false;
    mutable bool reAssembleFd2 = false;


    std::function<LinearSolver(const Linearization&)> solverCreator_ = [](const Linearization& f)
    {
        return makeDirectSolver<VariableSetDescription,VariableSetDescription>( f.A() ,
                                                                                f.range() ,
                                                                                f.domain() /*,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                DirectType::MUMPS ,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                MatrixProperties::GENERAL */);

    };
    std::shared_ptr<VectorSpace> operatorSpace_ = nullptr;

    //Überarbeitung notwendig
    std::function<::Spacy::Vector(const ::Spacy::Vector&)> zeroBlock_ = [this](const ::Spacy::Vector& x) mutable
    {
        return x*0.0;
    };

    template<class Result>
    Result getBBlock(const KaskadeOperator& K, unsigned int stateSize,
                     unsigned int controlSize) const
    {
        const Matrix& A(K.get());
        Result B;
        auto maxRows = 0;
        auto maxCols = 0;

        auto counter = 0u;
        for (int i = 0; i < A.nnz(); ++i)
        {
            if (A.cidx[i] >= stateSize
                    && A.ridx[i] >= controlSize + stateSize
                    && A.cidx[i] < controlSize + stateSize)
            {
                maxRows = std::max(maxRows, A.ridx[i]);
                maxCols = std::max(maxCols, A.cidx[i]);
                B.get_non_const().addEntry(A.ridx[i]- controlSize -stateSize, A.cidx[i]- stateSize, A.data[i]);
                counter++;
            }
        }

        //  std::cout << "A.nnz() " << counter << std::endl;

        // std::cout << "B Orig max rows: " << maxRows << std::endl;
        //    std::cout << "B Orig max cols: " << maxCols << std::endl;

        return B;
    }

    template<class Result>
    Result getBTBlock(const KaskadeOperator& K, unsigned int stateSize,
                      unsigned int controlSize) const
    {
        const Matrix& A(K.get());
        Result B;
        auto maxRows = 0;
        auto maxCols = 0;
        for (int i = 0; i < A.nnz(); ++i)
        {


            if (A.ridx[i] >= stateSize
                    && A.cidx[i] >= controlSize + stateSize
                    && A.ridx[i] < controlSize + stateSize)
            {
                maxRows = std::max(maxRows, A.ridx[i]);
                maxCols = std::max(maxCols, A.cidx[i]);

                B.get_non_const().addEntry(A.ridx[i]- stateSize, A.cidx[i]- stateSize- controlSize, A.data[i]);

            }
        }

        //  std::cout << "BT Orig max rows: " << maxRows << std::endl;
        //   std::cout << "BT Orig max cols: " << maxCols << std::endl;

        return B;
    }



};

/**
         * @brief Convenient generation of a twice differentiable functional \f$f: X\rightarrow \mathbb{R}\f$ for %Kaskade 7.
         * @param f operator definition from %Kaskade 7
         * @param domain domain space
         * @param rbegin first row to be considered in the definition of f
         * @param rend one after the last row to be considered in the definition of f
         * @param cbegin first column to be considered in the definition of f
         * @param cend one after the last column to be considered in the definition of f
         * @return @ref C2BlockFunctional "::Spacy::Kaskade::C2BlockFunctional<FunctionalDefinition>( f, domain, rbegin, rend, cbegin, cend )"
         *
         * The optional parameters rbegin, rend, cbegin and cend can be used to define operators that correspond to parts of
         * a system of equation.
         */
template <class FunctionalDefinition, class GridManager>
auto makeC2BlockFunctional(const FunctionalDefinition& f, const VectorSpace& domain, const GridManager & gridManager,
                           int rbegin = 0, int rend = FunctionalDefinition::AnsatzVars::noOfVariables,
                           int cbegin = 0, int cend = FunctionalDefinition::TestVars::noOfVariables)
{
    return C2BlockFunctional<FunctionalDefinition, GridManager>( f, domain , gridManager, rbegin, rend , cbegin , cend );
}
}
/** @} */
}
