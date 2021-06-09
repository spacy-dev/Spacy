// Copyright (C) 2015 by Lars Lubkoll. All rights reserved.
// Released under the terms of the GNU General Public License version 3 or later.

#pragma once

#include "linalg/direct.hh"
#include <Spacy/ZeroVectorCreator.h>
#include "Spacy/Vector.h"
#include "Spacy/Util/Cast.h"
#include "Spacy/Util/Base/OperatorBase.h"
#include "Copy.h"
#include "VectorSpace.h"
#include "linalg/symmetricOperators.hh"



namespace Spacy
{
namespace Kaskade
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
		class Prec>
class SingleBlockPreconditioner: public OperatorBase
{
	using Spaces = typename AnsatzVariableDescription::Spaces;
	using Domain = typename AnsatzVariableDescription::template CoefficientVectorRepresentation<>::type;
	using Range = typename TestVariableDescription::template CoefficientVectorRepresentation<>::type;
    
    using DomainVector = Vector<AnsatzVariableDescription>;
    using RangeVector = Vector<TestVariableDescription>;
public:
	SingleBlockPreconditioner() = delete;
	/**
	 * @brief Constructor.
	 * @param spaces boost fusion forward sequence of space pointers required to initialize temporary storage
	 * @param domain domain space of the solver
	 * @param range range space of the solver
	 */
	SingleBlockPreconditioner(const VectorSpace& domain, const VectorSpace& range,
			std::unique_ptr<Prec> prec) :
			OperatorBase(domain, range), prec_(std::move(prec)), spaces_(
					extractSpaces<AnsatzVariableDescription>(domain))
	{
	}

	/// Compute \f$A^{-1}x\f$.
	::Spacy::Vector operator()(const ::Spacy::Vector& x) const
	{
        checkSpaceCompatibility(x.space(),domain());

		// Changed Range and Domain in Kaskade
                      
                auto y = ::Spacy::zero( range() );
           
                const auto& x__ = cast_ref<DomainVector>(x);
                auto& y__ = cast_ref<RangeVector>(y);

                
		Apply<Domain, Range, Prec>::wrap( *prec_, ::boost::fusion::at_c<0>(x__.get().data).coefficients(), ::boost::fusion::at_c<0>(y__.get().data).coefficients() );


		return y;
	}

	void apply(::Spacy::Vector& y,const ::Spacy::Vector &x) const
	{
        checkSpaceCompatibility(x.space(),domain());
        checkSpaceCompatibility(y.space(),range());

		// Changed Range and Domain in Kaskade
                      
                const auto& x__ = cast_ref<DomainVector>(x);
                auto& y__ = cast_ref<RangeVector>(y);

                
		Apply<Domain, Range, Prec>::wrap( *prec_, ::boost::fusion::at_c<0>(x__.get().data).coefficients(), ::boost::fusion::at_c<0>(y__.get().data).coefficients() );
	}

private:
        
    // This is a pImpl, thus should be a std::unique_ptr, problem with that is that copying the Object is difficulat
    std::shared_ptr<Prec> prec_;
	Spaces spaces_;

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
		class Prec>
auto makeSingleBlockPreconditioner(const VectorSpace& domain, const VectorSpace& range,
		std::unique_ptr<Prec> prec)
{
	return SingleBlockPreconditioner<AnsatzVariableSetDescription,
			TestVariableSetDescription, Prec>(domain, range,std::move(prec));
}

}
}
