// Copyright (C) 2015 by Lars Lubkoll. All rights reserved.
// Released under the terms of the GNU General Public License version 3 or later.

#pragma once

#include "linalg/direct.hh"
#include <Spacy/ZeroVectorCreator.h>
#include "Spacy/Vector.h"
#include "Spacy/Util/Cast.h"
#include "Spacy/Util/Base/OperatorBase.h"
// #include "Copy.h"
#include "vectorSpace.hh"
#include "vector.hh"
// #include "linalg/symmetricOperators.hh"
#include "linalg/threadedMatrix.hh"

#include "c2Functional.hh"

#include "fem/assemble.hh"
#include "utilities/memory.hh"
#include "mg/additiveMultigrid.hh"
#include "mg/multiplicativeMultigrid.hh"
#include "linalg/jacobiPreconditioner.hh"
#include "linalg/icc0precond.hh"

namespace Spacy
{
namespace KaskadeParabolic
{
/**
 * @ingroup KaskadeGroup
 * @brief Wrapper class to call apply-routines from different Kaskade preconditioners the same way (is to be unified in Kaskade)
 * @param Prec  %Kaskade Preconditioner (i.e., %Kaskade::AdditiveMultiGrid < > )
 * @param range range space of the preconditioner
 * @param domain domain space of the preconditioner
 */

template<class Domain, class Range, class Prec>
class Apply
{
public:
	// const Refs ?
    template<class Dom, class Rge>
	static void wrap(Prec & prec, Dom & x, Rge & y)
	{
		prec.apply(y,x);
	}
};

/**
 * @ingroup KaskadeGroup
 * @brief Preconditioner interface for %Kaskade 7.
 * First template argument is domain which corresonds to the rhs Mx = b
 */
template<class AnsatzVariableDescription, class TestVariableDescription,
		class FunctionalDefinition, class Prec>
class SingleBlockPreconditioner: public OperatorBase
{
	using Spaces = typename AnsatzVariableDescription::Spaces;
	using Domain = typename AnsatzVariableDescription::template CoefficientVectorRepresentation<>::type;
	using Range = typename TestVariableDescription::template CoefficientVectorRepresentation<>::type;
    
    using DomainVector = Vector<AnsatzVariableDescription>;
    using RangeVector = Vector<TestVariableDescription>;

    using NumaMat = ::Kaskade::NumaBCRSMatrix<Dune::FieldMatrix<double,1,1>>;
  
public:
    using DomainDesc = AnsatzVariableDescription;
    using RangeDesc = TestVariableDescription;
    
	SingleBlockPreconditioner() = delete;
	/**
	 * @brief Constructor.
	 * @param spaces boost fusion forward sequence of space pointers required to initialize temporary storage
	 * @param domain domain space of the solver
	 * @param range range space of the solver
	 */
	SingleBlockPreconditioner(const VectorSpace& domain, const VectorSpace& range,
				  Spacy::KaskadeParabolic::C2Functional<FunctionalDefinition>& f, bool transposed = false) :
			OperatorBase(domain, range)
	{
        auto dtVec_ = f.getGridMan().getTempGrid().getDtVec();
        tEnd_ = dtVec_.size();
        if(!transposed)
        {
            for (int t = 0; t < tEnd_; t++)
            {
                auto StiffMat = f.assembler(t).template get<NumaMat>(false,2,3,0,1);
                StiffMat    *= dtVec_.at(t).get();
                
                auto MassMat = f.assemblerSP(t).template get<NumaMat>(false,2,3,0,1);
                
                auto Mat = MassMat + StiffMat;

                
                prec_.push_back(moveUnique(::Kaskade::makeJacobiMultiGrid(Mat, *f.getGridMan().getKaskGridMan().at(t) )));
            }
        }
        else
        {
            for (int t = 0; t < tEnd_-1; t++)    
            {
                auto StiffMat = f.assembler(t).template get<NumaMat>(false,0,1,2,3);
                StiffMat    *= dtVec_.at(t+1).get();

                auto MassMat = f.assemblerSP(t).template get<NumaMat>(false,0,1,2,3);

                auto Mat = MassMat + StiffMat;


                prec_.push_back(moveUnique(::Kaskade::makeJacobiMultiGrid(Mat, *f.getGridMan().getKaskGridMan().at(t) )));
            }
        }
	}

	/// Compute \f$A^{-1}x\f$.
	::Spacy::Vector operator()(const ::Spacy::Vector& x) const
	{
        checkSpaceCompatibility(x.space(),domain());
        
        auto y = ::Spacy::zero( range() );
    
        const auto& x__ = cast_ref<DomainVector>(x);
        auto& y__ = cast_ref<RangeVector>(y);
                
        for (int t = 0; t < tEnd_; t++)
        {
            Apply<Domain, Range, Prec>::wrap( *prec_.at(t), ::boost::fusion::at_c<0>(x__.get(t).data).coefficients(), ::boost::fusion::at_c<0>(y__.get_nonconst(t).data).coefficients() );
        }


		return y;
	}

	::Spacy::Vector applyFor(const ::Spacy::Vector& x, int tBegin, int tEnd) const
	{
        checkSpaceCompatibility(x.space(),domain());

        auto y = ::Spacy::zero( range() );
    
        const auto& x__ = cast_ref<DomainVector>(x);
        auto& y__ = cast_ref<RangeVector>(y);

        for (int t = tBegin; t < tEnd; t++)
            Apply<Domain, Range, Prec>::wrap( *prec_.at(t), ::boost::fusion::at_c<0>(x__.get(t).data).coefficients(), ::boost::fusion::at_c<0>(y__.get_nonconst(t).data).coefficients() );
        
		return y;
	}

private:
        
    // This is a pImpl, thus should be a std::unique_ptr, problem with that is that copying the Object is difficulat
    std::vector< std::shared_ptr<Prec> > prec_;
    int tEnd_;
};

/**
 * @brief Convenient generation of a Spacy preconditioner from a Kaskade preconditioner
 * @param A %Kaskade operator (i.e., AssembledGalerkinOperator or MatrixRepresentedOperator)
 * @param spaces boost fusion forward sequence of space pointers required to initialize temporary storage
 * @param domain domain space of the solver
 * @param range range space of the solver
 * @param Prec Kaskade preconditioner (e.g. Kaskade::AdditiveMultiGrid < > )
 * @return Preconditioner<AnsatzVariableSetDescription, TestVariableSetDescription, Preconditioner  >( domain , range , preconditioner )
 */
template<class AnsatzVariableSetDescription, class TestVariableSetDescription,
		class FunctionalDefinition>
auto makeSingleBlockPreconditioner(const VectorSpace& domain, const VectorSpace& range,
				   Spacy::KaskadeParabolic::C2Functional<FunctionalDefinition>& f, bool transposed = false)
{
	return SingleBlockPreconditioner< AnsatzVariableSetDescription,
			TestVariableSetDescription, FunctionalDefinition,
            decltype(::Kaskade::makeJacobiMultiGrid(f.assembler(0).template get<2,0>(), *f.getGridMan().getKaskGridMan().at(0) )) >
	  (domain, range, f, transposed);
}

}
}
