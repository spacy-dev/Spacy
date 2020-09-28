
// Copyright (C) 2015 by Lars Lubkoll. All rights reserved.
// Released under the terms of the GNU General Public License version 3 or later.

#ifndef SPACY_KASKADE_DIRECTBLOCKPRECONDITIONER_HH
#define SPACY_KASKADE_DIRECTBLOCKPRECONDITIONER_HH


#include <utility>


#include "fem/variables.hh"

#include "Spacy/LinearSolver.h"
#include "Spacy/Operator.h"
#include "Spacy/C2Functional.h"
#include "Spacy/Vector.h"
#include "Spacy/Spaces/ProductSpace/Vector.h"
#include "Spacy/Util/Base/OperatorBase.h"
//#include "Spacy/Adapter/Kaskade/Util.hh"
#include "Spacy/Adapter/Kaskade/Operator.h"
#include "Spacy/Adapter/Kaskade/preconditioner.hh"
#include <Spacy/ZeroVectorCreator.h>
#include <Spacy/VectorSpace.h>


//Added by Stoecklein Matthias
#include <Spacy/Algorithm/CG/CG.h>
#include <Spacy/Algorithm/CG/RegularizeViaPreconditioner.h>
#include <Spacy/Algorithm/CG/TerminationCriteria.h>
#include <iostream>
#include <boost/type_index.hpp>
#include <linalg/apcg.hh>
#include "mg/multigrid.hh"
#include "mg/additiveMultigrid.hh"

#include "linalg/symmetricOperators.hh"
#include "mg/pcg.hh"

//////////////////////////////////////

namespace Spacy
{
/// \cond
class VectorSpace;
/// \endcond

namespace Kaskade
{
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

template<typename Vector>
void writeMatlabVector(const Vector & vector)
{


    std::ofstream myfile;
    myfile.open ("rhs.dat");
    for (int i = 0; i < boost::fusion::at_c<0>(vector.data).size(); i++)
    {
        for(int j = 0; j < boost::fusion::at_c<0>(vector.data).operator[](i).size(); j++)
        {
            myfile << boost::fusion::at_c<0>(vector.data).operator[](i).operator[](j) << std::endl;
        }

    }
    for (int i = 0; i < boost::fusion::at_c<1>(vector.data).size(); i++)
    {
        for(int j = 0; j < boost::fusion::at_c<1>(vector.data).operator[](i).size(); j++)
        {
            myfile << boost::fusion::at_c<1>(vector.data).operator[](i).operator[](j) << std::endl;
        }

    }
    for (int i = 0; i < boost::fusion::at_c<2>(vector.data).size(); i++)
    {
        for(int j = 0; j < boost::fusion::at_c<2>(vector.data).operator[](i).size(); j++)
        {
            myfile << boost::fusion::at_c<2>(vector.data).operator[](i).operator[](j) << std::endl;
        }

    }
    myfile.close();
}



/**
 * @brief A Preconditioner that uses a direct solver
 *
 */
template<class FunctionalDefinition/*, class Preconditioner, class TestMatrix, class GridManager*/>
class DirectBlockPreconditioner: public OperatorBase
{
public:

    using Description = typename FunctionalDefinition::AnsatzVars;
    using SpacyFunctional = typename ::Spacy::Kaskade::C2Functional<FunctionalDefinition>;
    using AVD = typename SpacyFunctional::VariableSetDescription;
    using VariableSet = typename Description::VariableSet;
    using Spaces = typename AVD::Spaces;

    using KaskadeOperator = typename SpacyFunctional::KaskadeOperator;
    using Domain = typename AVD::template CoefficientVectorRepresentation<>::type;
    using Range = typename AVD::template CoefficientVectorRepresentation<>::type;
    template<class X, class Y>
    using M = ::Kaskade::MatrixRepresentedOperator<MatrixAsTriplet<double>,X,Y>;

    template<class X, class Y>
    using MyMatrixType = ::Kaskade::MatrixRepresentedOperator<NumaBCRSMatrix<Dune::FieldMatrix<double,3,3>>,X,Y>;

    using VYSetDescription = Detail::ExtractDescription_t<AVD,0>;
    using VUSetDescription = Detail::ExtractDescription_t<AVD,1>;
    using VPSetDescription = Detail::ExtractDescription_t<AVD,2>;
    using DomainY = typename VYSetDescription::template CoefficientVectorRepresentation<>::type;
    using DomainU = typename VUSetDescription::template CoefficientVectorRepresentation<>::type;
    using DomainP = typename VPSetDescription::template CoefficientVectorRepresentation<>::type;
    using Atype = M<DomainY,DomainP>;
    using ATtype = M<DomainP,DomainY>;
    using Btype = M<DomainU,DomainP>;
    using BTtype = M<DomainP,DomainU>;
    using Mutype = M<DomainU,DomainU>;

    using PSV = ::Spacy::ProductSpace::Vector;
    using PSVC = ::Spacy::ProductSpace::VectorCreator;


    /**
         * @brief Constructor.
         * @param domain domain space
         * @param range range space
         */
    DirectBlockPreconditioner(const KaskadeOperator& K,
                              /*  const ::Spacy::Vector& origin, */
                              const Description& desc, const VectorSpace& domain,
                              const VectorSpace& range/*, Preconditioner & preconditioner,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                TestMatrix & testMatrix, GridManager & gridManager*/) :
        OperatorBase(domain, range), desc_(desc), spaces_(
                                                      ::Spacy::Kaskade::extractSpaces<AVD>(domain)),

        /*, preconditioner_(
                                                                                                                        preconditioner), testMatrix_(testMatrix),
                                                                                                                        gridManager_(gridManager),*/ K_(K)

      /*     origin_(origin)          */
    {
        blockbounds.push_back(0);
        {
            typename Description::VariableSet originK(desc);
            // 								copy(origin_,originK);
            for (int i = 0; i < originK.descriptions.noOfVariables; i++)
            {
                blockbounds.push_back(
                            originK.descriptions.degreesOfFreedom(i, i + 1)
                            + blockbounds[i]);
            }
        }

        //~ ATtype AT;
        //~ AT=getBlock<ATtype>(K,0,2);
        //~ AT.get_non_const().setStartToZero();
        //~ solAT=std::make_shared< ::Kaskade::InverseLinearOperator< ::Kaskade::DirectSolver<DomainP,DomainY> > >
        //~ ( ::Kaskade::directInverseOperator(AT, DirectType::UMFPACK3264, MatrixProperties::GENERAL) );
        auto indexShift = 0;

        //  A Quadratic ???
        {
            std::cout << "  Factorization .." << std::flush;
            Atype A;
            A = getBlock<Atype>(K, 2, 0);
            A.get_non_const().setStartToZero();
            indexShift = A.get().N();

            //  writeMatlab(A.get_non_const());
            solA =
                    std::make_shared
                    < ::Kaskade::InverseLinearOperator<
                    ::Kaskade::DirectSolver<DomainY, DomainP> >
                    > (::Kaskade::directInverseOperator(A,
                                                        DirectType::UMFPACK3264,
                                                        MatrixProperties::GENERAL));
            std::cout << "." << std::endl;
        }

        Mutype Mu;
        Mu = getBlockDiag<Mutype>(K, 1, 1);

        // Todo check is this really is not necessary
        // Run with debugger

        for(size_t i = 0; i < Mu.get().nnz(); i++)
        {
            Mu.get_non_const().ridx[i] -= indexShift;
            Mu.get_non_const().cidx[i] -= indexShift;
        }


        const auto & MuTemp = Mu.get();


        //        std::cout << "MuTemp.N(): " << MuTemp.N() << std::endl;

        //        for(int i = 0; i < MuTemp.ridx.size(); i++)
        //        {
        //            if(MuTemp.ridx[i] == MuTemp.cidx[i])
        //            {
        //               std::cout << "ridx[i]: " << MuTemp.ridx[i] << " " << "data[i]: " << MuTemp.data[i] << std::endl;
        //            }
        //        }

        std::cout << std::endl;

        // Why does this work ?? Test Again if this is correct ??
        // Test U MassMatrix !!!
        // Mu.get_non_const().setStartToZero();
        // Mu_ = getBlockDiag<Mutype>(K, 1, 1);
        // Mu_.get_non_const().setStartToZero();

        solMu = std::make_shared
                < ::Kaskade::InverseLinearOperator<
                ::Kaskade::DirectSolver<DomainU, DomainU> >
                > (::Kaskade::directInverseOperator(Mu, DirectType::UMFPACK3264,
                                                    MatrixProperties::GENERAL));

        B = getBlock<Btype>(K, 2, 1);

        // check if this works
        for(size_t i = 0; i < B.get().nnz(); i++)
        {
            B.get_non_const().ridx[i] -= (indexShift + Mu.get().N());
            B.get_non_const().cidx[i] -= indexShift;
        }


        //     const auto & BTemp = B.get();


        //        std::cout << "BTemp.N(): " << MuTemp.N() << std::endl;
        //        std::cout << "BTemp.M(): " << MuTemp.M() << std::endl;

        //        int max = 0;
        //        for(int i = 0; i < BTemp.ridx.size(); i++)
        //        {
        //            max = std::max(max, BTemp.ridx[i]);

        //            if(BTemp.ridx[i] == BTemp.cidx[i])
        //            {
        //               std::cout << "ridx[i]: " << BTemp.ridx[i] << " " << "data[i]: " << BTemp.data[i] << std::endl;
        //            }
        //        }

        //        std::cout << "Max row: " <<  max << std::endl << std::endl;


        // Does not always work
        //B.get_non_const().setStartToZero();
    }

    /**
         * @brief Apply preconditioner \f$P\f$.
         * @param x argument
         * @return \f$P(x)\f$
         */
    ::Spacy::Vector operator()(const ::Spacy::Vector& b) const
    {

        using namespace boost::fusion;

        Domain b_(
                    AVD::template CoefficientVectorRepresentation<>::init(spaces_));
        Domain x_(
                    AVD::template CoefficientVectorRepresentation<>::init(spaces_));

        /*Domain x_dummy(AVD::template CoefficientVectorRepresentation<>::init(spaces_));*/

        copyToCoefficientVector<AVD>(b, b_);

        DomainY b0_(at_c<0>(b_.data));
        DomainU b1_(at_c<1>(b_.data));
        DomainP b2_(at_c<2>(b_.data));

        DomainY xy_(at_c<0>(x_.data));
        DomainU xu_(at_c<1>(x_.data));
        DomainP xp_(at_c<2>(x_.data));

        //solAT->apply(b0_,xp_);

        std::vector<double> p(xp_.dim());
        std::vector<double> b0(b0_.dim());

        b0_.write(b0.begin());

        // Solve A^T p=b0

        solA->op.solver->solve(b0, p, true);   // true means transposed

        //solAT->apply(b0_,xp_);
        xp_.read(p.begin());

        B.applyscaleaddTransposed(-1.0, xp_, b1_); // b1_ -> b1_ +(-1.0)*(-B^T xp_)

        // Solve M_u xu_ = b1_
        solMu->apply(b1_, xu_);


        B.applyscaleadd(-1.0, xu_, b2_); // b2_ -> b2_ +(-1.0)*(-B xu_)

        // Solve Axy_ = b2_
        solA->apply(b2_, xy_);

        at_c<0>(x_.data) = at_c<0>(xy_.data);
        at_c<1>(x_.data) = at_c<0>(xu_.data);
        at_c<2>(x_.data) = at_c<0>(xp_.data);

        //  Domain x_dummy(AVD::template CoefficientVectorRepresentation<>::init(spaces_));
        //  K_.apply(x_,x_dummy);

        // Check Positivity for x_ and kaskade operator

        //~ b2_ *= 0.0;
        //~
        //~ A.applyscaleadd(1.0,xy_,b2_);
        //~ B.applyscaleadd(1.0,xu_,b2_);
        //~
        //~
        //~ DomainY b0__(at_c<0>(b_.data));
        //~ DomainU b1__(at_c<1>(b_.data));
        //~ DomainP b2__(at_c<2>(b_.data));
        //~
        //~ std::cout << "Scalar products:" << xy_*b0__ << " " << xu_*b1__ << " " << xp_*b2__ << std::endl;
        //~ std::cout << "U Norms:" << xy_.two_norm() << " " << xu_.two_norm() << " " << xp_.two_norm() << std::endl;
        //~ std::cout << "R Norms:" << b0__.two_norm() << " " << b1__.two_norm() << " " << b2__.two_norm() << std::endl;

        //         auto x  = return space.creator()( &space )


        // Todo check domain() function !!
        //sfgs
        auto x = zero(domain());

        copyFromCoefficientVector<AVD>(x_, x);

        return x;
    };

private:
    const Description& desc_;
    Spaces spaces_;
    //         const ::Spacy::Vector& origin_;
    std::vector<int> blockbounds =
    { };
    Btype B;
    //Atype A;
    //	Mutype Mu_;
    mutable std::shared_ptr<
    ::Kaskade::InverseLinearOperator<
    ::Kaskade::DirectSolver<DomainY, DomainP> > > solA = nullptr;
    //~ mutable std::shared_ptr< ::Kaskade::InverseLinearOperator< ::Kaskade::DirectSolver<DomainP,DomainY> > > solAT = nullptr;
    mutable std::shared_ptr<
    ::Kaskade::InverseLinearOperator<
    ::Kaskade::DirectSolver<DomainU, DomainU> > > solMu =
            nullptr;

    //	Preconditioner & preconditioner_;
    //	TestMatrix testMatrix_;
    //	GridManager & gridManager_;
    const KaskadeOperator& K_;

    template<class Result>
    Result getBlock(const KaskadeOperator& K, unsigned int row,
                    unsigned int col)
    {
        const MatrixAsTriplet<double>& A(K.get());
        Result B;
        for (int i = 0; i < A.nnz(); ++i)
        {
            if (A.ridx[i] >= blockbounds[row]
                    && A.ridx[i] < blockbounds[row + 1]
                    && A.cidx[i] >= blockbounds[col]
                    && A.cidx[i] < blockbounds[col + 1])
                B.get_non_const().addEntry(A.ridx[i], A.cidx[i], A.data[i]);
        }
        return B;

    }
    template<class Result>
    Result getBlockDiag(const KaskadeOperator& K, unsigned int row,
                        unsigned int col)
    {
        const MatrixAsTriplet<double>& A(K.get());
        Result B;
        for (int i = 0; i < A.nnz(); ++i)
        {
            if (A.ridx[i] >= blockbounds[row]
                    && A.ridx[i] < blockbounds[row + 1]
                    && A.cidx[i] >= blockbounds[col]
                    && A.cidx[i] < blockbounds[col + 1]
                    && A.ridx[i] == A.cidx[i])
                B.get_non_const().addEntry(A.ridx[i], A.cidx[i], A.data[i]);
        }
        return B;

    }
};


template<class FunctionalDefinition, class GridManager>
class InexactDirectBlockPreconditioner: public OperatorBase
{
public:

    // todo more generic
    static constexpr int dim = 3;

    using Description = typename FunctionalDefinition::AnsatzVars;
    using SpacyFunctional = typename ::Spacy::Kaskade::C2Functional<FunctionalDefinition>;
    using AVD = typename SpacyFunctional::VariableSetDescription;
    using VariableSet = typename Description::VariableSet;
    using Spaces = typename AVD::Spaces;

    using KaskadeOperator = typename SpacyFunctional::KaskadeOperator;
    using Domain = typename AVD::template CoefficientVectorRepresentation<>::type;
    using Range = typename AVD::template CoefficientVectorRepresentation<>::type;

    template<class X, class Y>
    using M = ::Kaskade::MatrixRepresentedOperator<MatrixAsTriplet<double>,X,Y>;

    template<class X, class Y>
    using MyMatrixType = ::Kaskade::MatrixRepresentedOperator<NumaBCRSMatrix<Dune::FieldMatrix<double,dim,dim>>,X,Y>;

    using Vector = Dune::FieldVector<double,dim>;
    using BlockVector = Dune::BlockVector<Dune::FieldVector<double,dim>>;
    using Matrix = NumaBCRSMatrix<Dune::FieldMatrix<double,dim,dim>>;

    using VYSetDescription = Detail::ExtractDescription_t<AVD,0>;
    using VUSetDescription = Detail::ExtractDescription_t<AVD,1>;
    using VPSetDescription = Detail::ExtractDescription_t<AVD,2>;

    using CoefficientVectorY = typename VYSetDescription::template CoefficientVectorRepresentation<>::type;
    using CoefficientVectorU = typename VUSetDescription::template CoefficientVectorRepresentation<>::type;
    using CoefficientVectorP = typename VPSetDescription::template CoefficientVectorRepresentation<>::type;

    using DomainY = typename VYSetDescription::template CoefficientVectorRepresentation<>::type;
    using DomainU = typename VUSetDescription::template CoefficientVectorRepresentation<>::type;
    using DomainP = typename VPSetDescription::template CoefficientVectorRepresentation<>::type;

    using Atype = M<DomainY,DomainP>;
    using ATtype = M<DomainP,DomainY>;
    using Btype = M<DomainU,DomainP>;
    using BTtype = M<DomainP,DomainU>;
    using Mutype = M<DomainU,DomainU>;

    using PSV  = ::Spacy::ProductSpace::Vector;
    using PSVC = ::Spacy::ProductSpace::VectorCreator;

    using BPX = decltype(makeBPX(std::declval<Matrix>(),std::declval<GridManager>()));

    using LinOp  = Dune::MatrixAdapter<Matrix,BlockVector,BlockVector>;
    using SymOp  = SymmetricLinearOperatorWrapper<BlockVector,BlockVector>;


    /**
         * @brief Constructor.
         * @param domain domain space
         * @param range range space
         */
    InexactDirectBlockPreconditioner(const KaskadeOperator& K,
                                     const Description& desc, const VectorSpace& domain,
                                     const VectorSpace& range, const GridManager & gridManager) :
        OperatorBase(domain, range), desc_(desc), spaces_(
                                                      ::Spacy::Kaskade::extractSpaces<AVD>(domain)), K_(K), gridManager_(gridManager), domain_(domain)


    {

        //std::cout << std::endl;
        //std::cout << "Use InexactBlockPreconditioner: " << std::endl;
        //std::cout << std::endl;
        auto degreesOfFreedomControl = 0u;

        blockbounds.push_back(0);
        {
            typename Description::VariableSet originK(desc);
            degreesOfFreedomControl = originK.descriptions.degreesOfFreedom(1, 2);
            for (int i = 0; i < originK.descriptions.noOfVariables; i++)
            {
                blockbounds.push_back(
                            originK.descriptions.degreesOfFreedom(i, i + 1)
                            + blockbounds[i]);
            }
        }



        auto indexShift = 0;

        //  A Quadratic ???
        {
            //std::cout << "  Factorization .." << std::flush;

            A = getBlock<Atype>(K, 2, 0);
            A.get_non_const().setStartToZero();
            indexShift = A.get().N();

            const auto & ATemp = A.get();

            A_ = std::make_shared<Matrix>(NumaBCRSMatrix<Dune::FieldMatrix<double, dim, dim>>(*(ATemp.template toBCRS<dim>()),true));
            bpxPtr_ = std::make_shared<BPX>(makeBPX( *A_, gridManager_));

            //            solA =
            //                    std::make_shared
            //                    < ::Kaskade::InverseLinearOperator<
            //                    ::Kaskade::DirectSolver<DomainY, DomainP> >
            //                    > (::Kaskade::directInverseOperator(A,
            //                                                        DirectType::UMFPACK3264,
            //                                                        MatrixProperties::GENERAL));
            std::cout << "." << std::endl;
        }

        Mutype Mu;
        Mu = getBlockDiag<Mutype>(K, 1, 1);

        // Todo check is this really is not necessary
        // Run with debugger

        MuDiag_.resize(degreesOfFreedomControl);


        for(size_t i = 0; i < Mu.get().nnz(); i++)
        {
            Mu.get_non_const().ridx[i] -= indexShift;
            Mu.get_non_const().cidx[i] -= indexShift;

            if(Mu.get().ridx[i] == Mu.get().cidx[i])
                MuDiag_[Mu.get().ridx[i]] = Mu.get().data[i];
        }



        std::cout << std::endl;

        solMu = std::make_shared
                < ::Kaskade::InverseLinearOperator<
                ::Kaskade::DirectSolver<DomainU, DomainU> >
                > (::Kaskade::directInverseOperator(Mu, DirectType::UMFPACK3264,
                                                    MatrixProperties::GENERAL));

        B = getBlock<Btype>(K, 2, 1);

        // check if this works
        for(size_t i = 0; i < B.get().nnz(); i++)
        {
            B.get_non_const().ridx[i] -= (indexShift + Mu.get().N());
            B.get_non_const().cidx[i] -= indexShift;
        }

    }


    std::function<bool(bool isConvex, bool setConvex)> signalConvex_ = [] (bool isConvex, bool setConvex) {return true;};

    /**
         * @brief Apply preconditioner \f$P\f$.
         * @param x argument
         * @return \f$P(x)\f$
         */
    ::Spacy::Vector operator()(const ::Spacy::Vector& b) const
    {
        using namespace boost::fusion;


        // !!!!!!!! Always use the DomainY because otherwise CopyToCoefficientVector is not going to work
        // Spacy distinguishes between the same vector spaces

        AOperator_ = [this](const ::Spacy::Vector& x) mutable
        {
            Domain temp( AVD::template CoefficientVectorRepresentation<>::init(spaces_));

            CoefficientVectorY p(at_c<0>(temp.data));
            CoefficientVectorY ATp(at_c<0>(temp.data));
            CoefficientVectorP var(at_c<0>(temp.data));

            ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(x,p);

            ATp = 0.0;
            A.apply(p,var);

            boost::fusion::at_c<0>(ATp.data) = boost::fusion::at_c<0>(var.data);
            auto result = x;
            ::Spacy::Kaskade::copyFromCoefficientVector<VYSetDescription>(ATp,result);

            return result;
        };

        BPXOperator_ = [this](const ::Spacy::Vector& x) mutable
        {

            Domain temp( AVD::template CoefficientVectorRepresentation<>::init(spaces_));

            CoefficientVectorP p(at_c<2>(temp.data));
            CoefficientVectorY y(at_c<0>(temp.data));

            ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(x,y);

            bpxPtr_->apply(::Kaskade::component<0>(p),::Kaskade::component<0>(y));

            boost::fusion::at_c<0>(y.data) = boost::fusion::at_c<0>(p.data);
            auto result = x;

            ::Spacy::Kaskade::copyFromCoefficientVector<VYSetDescription>(y,result);

            return result;

        };


        MuDiagMatrixSolver_ = [this](const DomainU & u)
        {

            DomainU result(u);
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


        MuDiagMatrixOperator_ = [this](const DomainU & u)
        {
            DomainU result(u);
            unsigned int counter = 0;

            auto resultIter =  boost::fusion::at_c<0>(result.data).begin();

            for(const auto & it : boost::fusion::at_c<0>(u.data))
            {
                for(auto i = 0u; i < it.size(); i++)
                {
                    resultIter->operator[](i) = it[i] * MuDiag_[counter];
                    counter++;
                }

                resultIter++;
            }

            return result;
        };


        Domain b_(
                    AVD::template CoefficientVectorRepresentation<>::init(spaces_));

        Domain x_(
                    AVD::template CoefficientVectorRepresentation<>::init(spaces_));

        copyToCoefficientVector<AVD>(b, b_);

        DomainY b0_(at_c<0>(b_.data));
        DomainU b1_(at_c<1>(b_.data));
        DomainP b2_(at_c<2>(b_.data));

        DomainY xy_(at_c<0>(x_.data));
        DomainU xu_(at_c<1>(x_.data));
        DomainP xp_(at_c<2>(x_.data));

        std::vector<double> p(xp_.dim());
        std::vector<double> b0(b0_.dim());

        b0_.write(b0.begin());

        CoefficientVectorY y(at_c<0>(x_.data));

        CoefficientVectorY p_(at_c<0>(x_.data));

        CoefficientVectorY Ay(at_c<0>(x_.data));
        CoefficientVectorY rhsy_(at_c<0>(x_.data));
        CoefficientVectorU u(at_c<1>(x_.data));


        // Solve A^T p=b0

        //solA->op.solver->solve(b0, p, true);   // true means transposed



        // Use Cg instead of direct solver
        const auto& domain_ps = cast_ref<PSVC>(domain_.creator());
        const auto& Y = (cast_ref<PSVC>(domain_ps.subSpace(0).creator())).subSpace(0);


        ::Spacy::Vector rhs = zero(Y);
        auto rhsTemp = b0_;

        ::Spacy::Kaskade::copyFromCoefficientVector<VYSetDescription>(b0_,rhs);
        ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(rhs, rhsTemp);

        auto cg = ::Spacy::CG::Solver(AOperator_, BPXOperator_,  ::Spacy::CG::NoRegularization(), true);
        cg.set_eps(1e-11);
        cg.setRelativeAccuracy(1e-11);
        //cg.setVerbosityLevel(1);
        cg.setTerminationCriterion(CG::Termination::AdaptiveRelativeEnergyError());

        //std::cout << std::endl;
      //  std::cout << "Start CG in DBP: " << std::endl;

        auto result = cg.solve(zero(Y),rhs);

      //  std::cout << "Finished CG in DBP: " << std::endl;
      //  std::cout << std::endl;

        if(cg.indefiniteOperator())
        {
            if(get(norm(rhs)) != 0.0)
            {
                signalConvex_(true,false);
                std::cout << "norm(rhs): " << norm(rhs) << std::endl;
                std::cout << "Fail In indirectBlockPreconditioner: A is not positive definite: " << std::endl << std::endl;
                return -b;
            }
        }

        ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(result,p_);
        boost::fusion::at_c<0>(xp_.data) = boost::fusion::at_c<0>(p_.data);


        /////////////////
        // xp_.read(p.begin());

        //  std::cout << "TestTwoNorm 1: " << std::endl;

        //        p_ -= xp_;
        //        std::cout << boost::fusion::at_c<0>(p_.data).two_norm() << std::endl;

        //        std::cout << std::endl;

        B.applyscaleaddTransposed(-1.0, xp_, b1_); // b1_ -> b1_ +(-1.0)*(-B^T xp_)

        // Mu xu_ = b1_
        xu_ = MuDiagMatrixSolver_(b1_);
        // Solve M_u xu_ = b1_
        // solMu->apply(b1_, xu_);

        ////u = MuDiagMatrixOperator_(xu_);

        ////std::cout << "u = Mu * xu_ twoNorm: " << u.two_norm() << std::endl;
        ////std::cout << "b1_: " << b1_.two_norm() << std::endl;

        ////  u -= b1_;
        //// std::cout << "u -= b1_: " << u.two_norm() << std::endl;


        // b2_ -> b2_ +(-1.0)*(-B xu_)
        B.applyscaleadd(-1.0, xu_, b2_);

        // Solve Axy_ = b2_
        // solA->apply(b2_, xy_);



        // Solve Axy_ = b2_
        cg = ::Spacy::CG::Solver(AOperator_, BPXOperator_,  ::Spacy::CG::NoRegularization(), true);
        cg.set_eps(1e-11);
        cg.setRelativeAccuracy(1e-11);
        //cg.setVerbosityLevel(1);
        cg.setTerminationCriterion(CG::Termination::AdaptiveRelativeEnergyError());

        // b2_ has different type than VYSetDescription. Thus, copyFromCoefficientVector does not work
        boost::fusion::at_c<0>(rhsy_.data) = boost::fusion::at_c<0>(b2_.data);
        ::Spacy::Kaskade::copyFromCoefficientVector<VYSetDescription>(rhsy_,rhs);


        // Solve Axy_ = b2_
        //std::cout << "Starting CG in DBP: " << std::endl;
        result = cg.solve(zero(Y),rhs);

       // std::cout << "Finished CG in DBP: " << std::endl;
       // std::cout << std::endl;


        if(cg.indefiniteOperator())
        {
            if(get(norm(rhs)) != 0.0)
            {
                signalConvex_(true,false);
                std::cout << "norm(rhs): " << norm(rhs) << std::endl;
                std::cout << "Fail In indirectBlockPreconditioner: A is not positive definite: " << std::endl << std::endl;
                return -b;
            }
        }

        ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(result,y);

        boost::fusion::at_c<0>(xy_.data) = boost::fusion::at_c<0>(y.data);


        //        std::cout << "TestTwoNorm 2: " << std::endl;
        //        y -= xy_;
        //        std::cout << boost::fusion::at_c<0>(y.data).two_norm() << std::endl;
        //        std::cout << std::endl;


        at_c<0>(x_.data) = at_c<0>(xy_.data);
        at_c<1>(x_.data) = at_c<0>(xu_.data);
        at_c<2>(x_.data) = at_c<0>(xp_.data);

        auto x = zero(domain());
        copyFromCoefficientVector<AVD>(x_, x);

        // test nur x zurück geben

        return x;
    };

private:
    const Description& desc_;
    Spaces spaces_;

    std::vector<int> blockbounds =
    { };
    Btype B;
    Atype A;
    mutable std::shared_ptr<
    ::Kaskade::InverseLinearOperator<
    ::Kaskade::DirectSolver<DomainY, DomainP> > > solA = nullptr;

    mutable std::shared_ptr<
    ::Kaskade::InverseLinearOperator<
    ::Kaskade::DirectSolver<DomainU, DomainU> > > solMu =
            nullptr;

    const KaskadeOperator& K_;
    const GridManager & gridManager_;
    const VectorSpace& domain_;

    std::shared_ptr<Matrix> A_ = nullptr;
    mutable std::shared_ptr<BPX> bpxPtr_ = nullptr;

    std::vector<double> MuDiag_ = {};

    mutable std::function<DomainU(const DomainU & u)> MuDiagMatrixOperator_ = [this](const DomainU & u) { return u; };
    mutable std::function<DomainU(const DomainU & u)> MuDiagMatrixSolver_ = [this](const DomainU & u) { return u; };

    mutable std::function<::Spacy::Vector(const ::Spacy::Vector&)> AOperator_ = [this](const ::Spacy::Vector & x) {return x;};
    mutable std::function<::Spacy::Vector(const ::Spacy::Vector&)> BPXOperator_ = [this](const ::Spacy::Vector& x) mutable
    {
        return x;
    };

    template<class Result>
    Result getBlock(const KaskadeOperator& K, unsigned int row,
                    unsigned int col)
    {
        const MatrixAsTriplet<double>& A(K.get());
        Result B;
        for (int i = 0; i < A.nnz(); ++i)
        {
            if (A.ridx[i] >= blockbounds[row]
                    && A.ridx[i] < blockbounds[row + 1]
                    && A.cidx[i] >= blockbounds[col]
                    && A.cidx[i] < blockbounds[col + 1])
                B.get_non_const().addEntry(A.ridx[i], A.cidx[i], A.data[i]);
        }
        return B;
    }


    template<class Result>
    Result getBlockDiag(const KaskadeOperator& K, unsigned int row,
                        unsigned int col)
    {
        const MatrixAsTriplet<double>& A(K.get());
        Result B;
        for (int i = 0; i < A.nnz(); ++i)
        {
            if (A.ridx[i] >= blockbounds[row]
                    && A.ridx[i] < blockbounds[row + 1]
                    && A.cidx[i] >= blockbounds[col]
                    && A.cidx[i] < blockbounds[col + 1]
                    && A.ridx[i] == A.cidx[i])
                B.get_non_const().addEntry(A.ridx[i], A.cidx[i], A.data[i]);
        }
        return B;

    }
};



//template<class FunctionalDefinition, class GridManager, class RegOperator>
//class InexactNonConvexDirectBlockPreconditioner: public OperatorBase
//{
//public:

//    // todo more generic
//    static constexpr int dim = 3;

//    using Description = typename FunctionalDefinition::AnsatzVars;
//    using SpacyFunctional = typename ::Spacy::Kaskade::C2Functional<FunctionalDefinition>;
//    using AVD = typename SpacyFunctional::VariableSetDescription;
//    using VariableSet = typename Description::VariableSet;
//    using Spaces = typename AVD::Spaces;

//    using KaskadeOperator = typename SpacyFunctional::KaskadeOperator;
//    using Domain = typename AVD::template CoefficientVectorRepresentation<>::type;
//    using Range = typename AVD::template CoefficientVectorRepresentation<>::type;

//    template<class X, class Y>
//    using M = ::Kaskade::MatrixRepresentedOperator<MatrixAsTriplet<double>,X,Y>;

//    template<class X, class Y>
//    using MyMatrixType = ::Kaskade::MatrixRepresentedOperator<NumaBCRSMatrix<Dune::FieldMatrix<double,dim,dim>>,X,Y>;

//    using Vector = Dune::FieldVector<double,dim>;
//    using BlockVector = Dune::BlockVector<Dune::FieldVector<double,dim>>;
//    using Matrix = NumaBCRSMatrix<Dune::FieldMatrix<double,dim,dim>>;

//    using VYSetDescription = Detail::ExtractDescription_t<AVD,0>;
//    using VUSetDescription = Detail::ExtractDescription_t<AVD,1>;
//    using VPSetDescription = Detail::ExtractDescription_t<AVD,2>;

//    using CoefficientVectorY = typename VYSetDescription::template CoefficientVectorRepresentation<>::type;
//    using CoefficientVectorU = typename VUSetDescription::template CoefficientVectorRepresentation<>::type;
//    using CoefficientVectorP = typename VPSetDescription::template CoefficientVectorRepresentation<>::type;

//    using DomainY = typename VYSetDescription::template CoefficientVectorRepresentation<>::type;
//    using DomainU = typename VUSetDescription::template CoefficientVectorRepresentation<>::type;
//    using DomainP = typename VPSetDescription::template CoefficientVectorRepresentation<>::type;

//    using Atype = M<DomainY,DomainP>;
//    using ATtype = M<DomainP,DomainY>;
//    using Btype = M<DomainU,DomainP>;
//    using BTtype = M<DomainP,DomainU>;
//    using Mutype = M<DomainU,DomainU>;

//    using PSV  = ::Spacy::ProductSpace::Vector;
//    using PSVC = ::Spacy::ProductSpace::VectorCreator;

//    using BPX = decltype(makeBPX(std::declval<Matrix>(),std::declval<GridManager>()));

//    using LinOp  = Dune::MatrixAdapter<Matrix,BlockVector,BlockVector>;
//    using SymOp  = SymmetricLinearOperatorWrapper<BlockVector,BlockVector>;


//    /**
//         * @brief Constructor.
//         * @param domain domain space
//         * @param range range space
//         */
//    InexactNonConvexDirectBlockPreconditioner(const KaskadeOperator& K,
//                                              const Description& desc, const VectorSpace& domain,
//                                              const VectorSpace& range, const GridManager & gridManager,  BPX & positiveBPX, const RegOperator &regOperator) :
//        OperatorBase(domain, range), desc_(desc), spaces_(
//                                                      ::Spacy::Kaskade::extractSpaces<AVD>(domain)), K_(K), gridManager_(gridManager), domain_(domain), positiveBPX_(positiveBPX), regOperator_(regOperator)


//    {
//        auto degreesOfFreedomControl = 0u;

//        blockbounds.push_back(0);
//        {
//            typename Description::VariableSet originK(desc);
//            degreesOfFreedomControl = originK.descriptions.degreesOfFreedom(1, 2);
//            for (int i = 0; i < originK.descriptions.noOfVariables; i++)
//            {
//                blockbounds.push_back(
//                            originK.descriptions.degreesOfFreedom(i, i + 1)
//                            + blockbounds[i]);
//            }
//        }



//        auto indexShift = 0;


//        {
//           // std::cout << "  Factorization .." << std::flush;

//            A = getBlock<Atype>(K, 2, 0);
//            A.get_non_const().setStartToZero();
//            indexShift = A.get().N();

//            const auto & ATemp = A.get();

//            A_ = std::make_shared<Matrix>(NumaBCRSMatrix<Dune::FieldMatrix<double, dim, dim>>(*(ATemp.template toBCRS<dim>()),true));
//            bpxPtr_ = std::make_shared<BPX>(makeBPX( *A_, gridManager_));

//            std::cout << "." << std::endl;
//        }

//        Mutype Mu;
//        Mu = getBlockDiag<Mutype>(K, 1, 1);

//        // Todo check is this really is not necessary
//        // Run with debugger

//        MuDiag_.resize(degreesOfFreedomControl);


//        for(size_t i = 0; i < Mu.get().nnz(); i++)
//        {
//            Mu.get_non_const().ridx[i] -= indexShift;
//            Mu.get_non_const().cidx[i] -= indexShift;

//            if(Mu.get().ridx[i] == Mu.get().cidx[i])
//                MuDiag_[Mu.get().ridx[i]] = Mu.get().data[i];
//        }



//        std::cout << std::endl;

//        solMu = std::make_shared
//                < ::Kaskade::InverseLinearOperator<
//                ::Kaskade::DirectSolver<DomainU, DomainU> >
//                > (::Kaskade::directInverseOperator(Mu, DirectType::UMFPACK3264,
//                                                    MatrixProperties::GENERAL));

//        B = getBlock<Btype>(K, 2, 1);

//        // check if this works
//        for(size_t i = 0; i < B.get().nnz(); i++)
//        {
//            B.get_non_const().ridx[i] -= (indexShift + Mu.get().N());
//            B.get_non_const().cidx[i] -= indexShift;
//        }

//    }


//    std::function<bool(bool isConvex, bool setConvex)> signalConvex_ = [] (bool isConvex, bool setConvex) {return true;};

//    /**
//         * @brief Apply preconditioner \f$P\f$.
//         * @param x argument
//         * @return \f$P(x)\f$
//         */
//    ::Spacy::Vector operator()(const ::Spacy::Vector& b) const
//    {
//        using namespace boost::fusion;


//        // !!!!!!!! Always use the DomainY because otherwise CopyToCoefficientVector is not going to work
//        // Spacy distinguishes between the same vector spaces

//        AOperator_ = [this](const ::Spacy::Vector& x) mutable
//        {
//            Domain temp( AVD::template CoefficientVectorRepresentation<>::init(spaces_));

//            CoefficientVectorY p(at_c<0>(temp.data));
//            CoefficientVectorY ATp(at_c<0>(temp.data));
//            CoefficientVectorP var(at_c<0>(temp.data));

//            ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(x,p);

//            ATp = 0.0;
//            A.apply(p,var);

//            boost::fusion::at_c<0>(ATp.data) = boost::fusion::at_c<0>(var.data);
//            auto result = x;
//            ::Spacy::Kaskade::copyFromCoefficientVector<VYSetDescription>(ATp,result);

//            return result;
//        };

//        BPXOperator_ = [this](const ::Spacy::Vector& x) mutable
//        {

//            Domain temp( AVD::template CoefficientVectorRepresentation<>::init(spaces_));

//            CoefficientVectorP p(at_c<2>(temp.data));
//            CoefficientVectorY y(at_c<0>(temp.data));

//            ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(x,y);

//            bpxPtr_->apply(::Kaskade::component<0>(p),::Kaskade::component<0>(y));

//            boost::fusion::at_c<0>(y.data) = boost::fusion::at_c<0>(p.data);
//            auto result = x;

//            ::Spacy::Kaskade::copyFromCoefficientVector<VYSetDescription>(y,result);

//            return result;

//        };


//        MuDiagMatrixSolver_ = [this](const DomainU & u)
//        {

//            DomainU result(u);
//            unsigned int counter = 0;

//            auto resultIter =  boost::fusion::at_c<0>(result.data).begin();

//            for(const auto & it : boost::fusion::at_c<0>(u.data))
//            {
//                for(auto i = 0u; i < it.size(); i++)
//                {
//                    resultIter->operator[](i) = it.operator[](i)/MuDiag_[counter];
//                    counter++;
//                }

//                resultIter++;
//            }

//            return result;
//        };


//        MuDiagMatrixOperator_ = [this](const DomainU & u)
//        {
//            DomainU result(u);
//            unsigned int counter = 0;

//            auto resultIter =  boost::fusion::at_c<0>(result.data).begin();

//            for(const auto & it : boost::fusion::at_c<0>(u.data))
//            {
//                for(auto i = 0u; i < it.size(); i++)
//                {
//                    resultIter->operator[](i) = it[i] * MuDiag_[counter];
//                    counter++;
//                }

//                resultIter++;
//            }

//            return result;
//        };


//        positivePrecond_ = [this](const ::Spacy::Vector& x) mutable
//        {

//            Domain temp( AVD::template CoefficientVectorRepresentation<>::init(spaces_));

//            CoefficientVectorP p(at_c<2>(temp.data));
//            CoefficientVectorY y(at_c<0>(temp.data));

//            ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(x,y);

//            positiveBPX_.apply(::Kaskade::component<0>(p),::Kaskade::component<0>(y));

//            boost::fusion::at_c<0>(y.data) = boost::fusion::at_c<0>(p.data);
//            auto result = x;

//            ::Spacy::Kaskade::copyFromCoefficientVector<VYSetDescription>(y,result);

//            return result;

//        };




//        Domain b_(
//                    AVD::template CoefficientVectorRepresentation<>::init(spaces_));

//        Domain x_(
//                    AVD::template CoefficientVectorRepresentation<>::init(spaces_));

//        copyToCoefficientVector<AVD>(b, b_);

//        DomainY b0_(at_c<0>(b_.data));
//        DomainU b1_(at_c<1>(b_.data));
//        DomainP b2_(at_c<2>(b_.data));

//        DomainY xy_(at_c<0>(x_.data));
//        DomainU xu_(at_c<1>(x_.data));
//        DomainP xp_(at_c<2>(x_.data));

//        std::vector<double> p(xp_.dim());
//        std::vector<double> b0(b0_.dim());

//        b0_.write(b0.begin());

//        CoefficientVectorY y(at_c<0>(x_.data));

//        CoefficientVectorY p_(at_c<0>(x_.data));

//        CoefficientVectorY Ay(at_c<0>(x_.data));
//        CoefficientVectorY rhsy_(at_c<0>(x_.data));
//        CoefficientVectorU u(at_c<1>(x_.data));






//        // Use Cg instead of direct solver
//        const auto& domain_ps = cast_ref<PSVC>(domain_.creator());
//        const auto& Y = (cast_ref<PSVC>(domain_ps.subSpace(0).creator())).subSpace(0);


//        ::Spacy::Vector rhs = zero(Y);
//        auto rhsTemp = b0_;

//        ::Spacy::Kaskade::copyFromCoefficientVector<VYSetDescription>(b0_,rhs);
//        ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(rhs, rhsTemp);

//        auto cg = ::Spacy::CG::Solver(AOperator_, BPXOperator_,  ::Spacy::CG::NoRegularization(), true);
//        cg.set_eps(1e-8);
//        cg.setRelativeAccuracy(1e-8);
//        //cg.setVerbosityLevel(1);
//        cg.setTerminationCriterion(CG::Termination::AdaptiveRelativeEnergyError());

//        //std::cout << std::endl;
//     //   std::cout << "Start CG in DBP: " << std::endl;

//        auto result = cg.solve(zero(Y),rhs);

//       // std::cout << "Finished CG in DBP: " << std::endl;
//      //  std::cout << std::endl;


//        // Todo A might get solved with different regularizations
//        if((cg.indefiniteOperator()))
//        {
//         //   std::cout << "CG not definit: " << std::endl;
//            if(get(norm(rhs)) != 0.0)
//            {
//           //     std::cout << "CG not definit2: " << std::endl;
//                auto newCg = ::Spacy::CG::Solver( AOperator_,positivePrecond_, ::Spacy::CG::RegularizeViaCallableOperator(regOperator_, 0.0),false );
//                newCg.setRelativeAccuracy(1e-6);
//                newCg.setVerbosity(2);
//                newCg.setTerminationCriterion(CG::Termination::AdaptiveRelativeEnergyError());
//                result = newCg.solve(zero(Y),rhs);
//            }
//        }



//        ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(result,p_);
//        boost::fusion::at_c<0>(xp_.data) = boost::fusion::at_c<0>(p_.data);


//        /////////////////


//        B.applyscaleaddTransposed(-1.0, xp_, b1_); // b1_ -> b1_ +(-1.0)*(-B^T xp_)

//        // Mu xu_ = b1_
//        xu_ = MuDiagMatrixSolver_(b1_);

//        // b2_ -> b2_ +(-1.0)*(-B xu_)
//        B.applyscaleadd(-1.0, xu_, b2_);

//        // Solve Axy_ = b2_
//        cg = ::Spacy::CG::Solver(AOperator_, BPXOperator_,  ::Spacy::CG::NoRegularization(), true);
//        cg.set_eps(1e-8);
//        cg.setRelativeAccuracy(1e-8);
//       // cg.setVerbosityLevel(2);
//        cg.setTerminationCriterion(CG::Termination::AdaptiveRelativeEnergyError());

//        // b2_ has different type than VYSetDescription. Thus, copyFromCoefficientVector does not work
//        boost::fusion::at_c<0>(rhsy_.data) = boost::fusion::at_c<0>(b2_.data);
//        ::Spacy::Kaskade::copyFromCoefficientVector<VYSetDescription>(rhsy_,rhs);

//        // Solve Axy_ = b2_
//        result = cg.solve(zero(Y),rhs);

//        //std::cout << "Finished CG in DBP: " << std::endl;
//        //std::cout << std::endl;


//        if((cg.indefiniteOperator()))
//        {
//             //  std::cout << "CG not definit: " << std::endl;
//            if(get(norm(rhs)) != 0.0)
//            {
//                    //          std::cout << "CG not definit3: " << std::endl;
//                auto newCg = ::Spacy::CG::Solver( AOperator_,positivePrecond_, ::Spacy::CG::RegularizeViaCallableOperator(regOperator_, 0.0),false );
//                newCg.setRelativeAccuracy(1e-6);
//                newCg.setVerbosity(2);
//                newCg.setTerminationCriterion(CG::Termination::AdaptiveRelativeEnergyError());
//                result = newCg.solve(zero(Y),rhs);
//            }
//        }

//        ::Spacy::Kaskade::copyToCoefficientVector<VYSetDescription>(result,y);
//        boost::fusion::at_c<0>(xy_.data) = boost::fusion::at_c<0>(y.data);

//        at_c<0>(x_.data) = at_c<0>(xy_.data);
//        at_c<1>(x_.data) = at_c<0>(xu_.data);
//        at_c<2>(x_.data) = at_c<0>(xp_.data);

//        auto x = zero(domain());
//        copyFromCoefficientVector<AVD>(x_, x);

//        // test nur x zurück geben

//        return x;
//    };

//private:
//    const Description& desc_;
//    Spaces spaces_;

//    std::vector<int> blockbounds =
//    { };
//    Btype B;
//    Atype A;
//    mutable std::shared_ptr<
//    ::Kaskade::InverseLinearOperator<
//    ::Kaskade::DirectSolver<DomainY, DomainP> > > solA = nullptr;

//    mutable std::shared_ptr<
//    ::Kaskade::InverseLinearOperator<
//    ::Kaskade::DirectSolver<DomainU, DomainU> > > solMu =
//            nullptr;

//    const KaskadeOperator& K_;
//    const GridManager & gridManager_;
//    const VectorSpace& domain_;

//    std::shared_ptr<Matrix> A_ = nullptr;
//    mutable std::shared_ptr<BPX> bpxPtr_ = nullptr;

//    std::vector<double> MuDiag_ = {};
//    mutable std::function<::Spacy::Vector(const ::Spacy::Vector&)> positivePrecond_= [this](const ::Spacy::Vector& x) mutable
//    {
//        return x;
//    };

//    BPX & positiveBPX_;
//    RegOperator     regOperator_;

//    mutable std::function<DomainU(const DomainU & u)> MuDiagMatrixOperator_ = [this](const DomainU & u) { return u; };
//    mutable std::function<DomainU(const DomainU & u)> MuDiagMatrixSolver_ = [this](const DomainU & u) { return u; };

//    mutable std::function<::Spacy::Vector(const ::Spacy::Vector&)> AOperator_ = [this](const ::Spacy::Vector & x) {return x;};
//    mutable std::function<::Spacy::Vector(const ::Spacy::Vector&)> BPXOperator_ = [this](const ::Spacy::Vector& x) mutable
//    {
//        return x;
//    };

//    template<class Result>
//    Result getBlock(const KaskadeOperator& K, unsigned int row,
//                    unsigned int col)
//    {
//        const MatrixAsTriplet<double>& A(K.get());
//        Result B;
//        for (int i = 0; i < A.nnz(); ++i)
//        {
//            if (A.ridx[i] >= blockbounds[row]
//                    && A.ridx[i] < blockbounds[row + 1]
//                    && A.cidx[i] >= blockbounds[col]
//                    && A.cidx[i] < blockbounds[col + 1])
//                B.get_non_const().addEntry(A.ridx[i], A.cidx[i], A.data[i]);
//        }
//        return B;

//    }
//    template<class Result>
//    Result getBlockDiag(const KaskadeOperator& K, unsigned int row,
//                        unsigned int col)
//    {
//        const MatrixAsTriplet<double>& A(K.get());
//        Result B;
//        for (int i = 0; i < A.nnz(); ++i)
//        {
//            if (A.ridx[i] >= blockbounds[row]
//                    && A.ridx[i] < blockbounds[row + 1]
//                    && A.cidx[i] >= blockbounds[col]
//                    && A.cidx[i] < blockbounds[col + 1]
//                    && A.ridx[i] == A.cidx[i])
//                B.get_non_const().addEntry(A.ridx[i], A.cidx[i], A.data[i]);
//        }
//        return B;

//    }
//};





}
}

#endif // SPACY_KASKADE_DIRECTBLOCKPRECONDITIONER_HH
