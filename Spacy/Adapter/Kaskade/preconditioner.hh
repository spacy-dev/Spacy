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

template<class Range, class Domain, class Prec>
class Apply
{

public:
	// const Refs ?
    template<class Rge, class Dom>
	static void wrap(Prec & prec, Rge & y, Dom & x)
	{
		prec.apply(y, x);
	}

};
/**
 * @ingroup KaskadeGroup
 * @brief Wrapper class to call apply-routines from different Kaskade preconditioners the same way (is to be unified in Kaskade)
 * @param Prec  %Kaskade Preconditioner (i.e., %Kaskade::AdditiveMultiGrid < > )
 * @param range range space of the preconditioner
 * @param domain domain space of the preconditioner
 */

// TODO
// Domain Range Vertauscht ???
template<class Range, class Domain, class T, class U, class V>
class Apply<Range, Domain, ::Kaskade::AdditiveMultiGrid<T, U, V> >
{
public:
    template<class Rge, class Dom>
	static void wrap(::Kaskade::AdditiveMultiGrid<T, U, V> & prec, Rge & y,
			Dom & x)
	{
		prec.apply(y, x);
//		prec.apply(::Kaskade::component<0>(y), ::Kaskade::component<0>(x));
	}

};





/** NOT TESTED
 * @ingroup KaskadeGroup
 * @brief Simple Jacobi Preconditioner
 * @param range range space of the preconditioner
 * @param domain domain space of the preconditioner
 */
template<class Domain, class Range, class KaskadeOperator>
class SimpleJacobiPreconditioner: public ::Kaskade::SymmetricPreconditioner<Domain, Range>
{
public:

	using Base = ::Kaskade::SymmetricPreconditioner<Domain,Range>;

	SimpleJacobiPreconditioner() = default;

	SimpleJacobiPreconditioner(const SimpleJacobiPreconditioner &) = default;

	// Not nice
	SimpleJacobiPreconditioner(SimpleJacobiPreconditioner &&) = default;

	SimpleJacobiPreconditioner(const KaskadeOperator &A) :
			init_(true)
	{

		const auto n = A.get().N();
		diag_.resize(n, 0.0);
		const auto size = A.get().size();
		const auto & rows = A.get().ridx;
		const auto & cols = A.get().cidx;
		const auto & data = A.get().data;

		size_t i = 0;

		while (i < size)
		{
			// use iterators
			if (rows[i] == cols[i])
			{
				diag_[rows[i]] = data[i];
			}

			i++;
		}

	}

	 virtual void apply( Domain & domain, const Range & range) override
	{

//


	auto & x = domain;

	const auto & b = range;

		size_t size = x.size();
		assert(init_ && size == b.size() && x.dim() == b.dim() == diag_.size());

		auto xIter = x.begin();
		auto diagIter = diag_.cbegin();

		for (const auto & it : b)
		{
			auto vecIter = xIter->begin();

			for (const auto & entry : it)
			{
				*(vecIter) = entry / (*diagIter);
				diagIter++;
				vecIter++;
			}
			xIter++;
		}


	}

	  typename Base::field_type applyDp(Domain& x, Range const& r) override
	        {

	         apply(x,r);

	         return x*r;
	        }

	  bool requiresInitializedInput() const
	         {



	           return false;
	         }




	void init(const KaskadeOperator &A)
	{
;
		init_ = true;
		const auto n = A.get().N();
		diag_.resize(n, 0.0);
		const auto size = A.get().size();
		const auto & rows = A.get().ridx;
		const auto & cols = A.get().cidx;
		const auto & data = A.get().data;

		size_t i = 0;

		while (i < size)
		{
			// use iterators
			if (rows[i] == cols[i])
			{
				diag_[rows[i]] = data[i];
			}
			i++;
		}

	}

private:
// not well implemented

	std::vector<double> diag_;
	bool init_ = false;

};




template<class Range, class Domain, class KaskadeOperator>
class Apply<Range, Domain, ::Spacy::Kaskade::SimpleJacobiPreconditioner< Range, Domain, KaskadeOperator > >
{
public:
    template<class Rge, class Dom>
	static void wrap(::Spacy::Kaskade::SimpleJacobiPreconditioner< Range, Domain, KaskadeOperator > & prec, Rge & y,
			Dom & x)
	{
		prec.apply(::Kaskade::component<0>(y), ::Kaskade::component<0>(x));
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

		// Changed Range and Domain in Kaskade
                      
                auto y = ::Spacy::zero( range() );
           
                const auto& x__ = cast_ref<DomainVector>(x);
                auto& y__ = cast_ref<RangeVector>(y);

                
		Apply<Range, Domain, Prec>::wrap(*prec_, ::boost::fusion::at_c<0>(y__.get().data), ::boost::fusion::at_c<0>(x__.get().data));


		return y;
	}

	void apply(::Spacy::Vector& y,const ::Spacy::Vector &x) const
	{

		// Changed Range and Domain in Kaskade
                      
                const auto& x__ = cast_ref<DomainVector>(x);
                auto& y__ = cast_ref<RangeVector>(y);

                
		Apply<Range, Domain, Prec>::wrap(*prec_, ::boost::fusion::at_c<0>(y__.get().data), ::boost::fusion::at_c<0>(x__.get().data));
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

// template<class AnsatzVariableSetDescription, class TestVariableSetDescription,
// 		class Prec>
// auto makeSingleBlockPreconditioner(const VectorSpace& domain, const VectorSpace& range,
// 		const Prec& prec)
// {
// 	return SingleBlockPreconditioner<AnsatzVariableSetDescription,
// 			TestVariableSetDescription, Prec>(domain, range,prec);
// }

}
}
