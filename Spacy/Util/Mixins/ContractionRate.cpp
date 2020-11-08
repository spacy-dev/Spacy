#include "ContractionRate.h"

#include <cassert>
#include <utility>

namespace Spacy::Mixin
{
    ContractionRate::ContractionRate( Real desiredContraction, Real relaxedDesiredContraction, Real maximalContraction ) noexcept
        : desiredContraction_( std::move( desiredContraction ) ), relaxedDesiredContraction_( std::move( relaxedDesiredContraction ) ),
          maximalContraction_( std::move( maximalContraction ) )
    {
        assert( desiredContraction_ < relaxedDesiredContraction_ );
        assert( relaxedDesiredContraction_ < maximalContraction_ );
    }

    void ContractionRate::setContraction( Real contraction ) noexcept
    {
        contraction_ = contraction;
    }

    void ContractionRate::setDesiredContraction( Real desiredContraction ) noexcept
    {
        desiredContraction_ = desiredContraction;
    }

    void ContractionRate::setRelaxedDesiredContraction( Real relaxedDesiredContraction ) noexcept
    {
        relaxedDesiredContraction_ = relaxedDesiredContraction;
    }

    void ContractionRate::setMaximalContraction( Real maximalContraction ) noexcept
    {
        maximalContraction_ = maximalContraction;
    }

    Real ContractionRate::getContraction() const noexcept
    {
        return contraction_;
    }

    Real ContractionRate::getDesiredContraction() const noexcept
    {
        return desiredContraction_;
    }

    Real ContractionRate::getRelaxedDesiredContraction() const noexcept
    {
        return relaxedDesiredContraction_;
    }

    Real ContractionRate::getMaximalContraction() const noexcept
    {
        return maximalContraction_;
    }

    bool ContractionRate::contractionIsAdmissible() const noexcept
    {
        return getContraction() < getMaximalContraction();
    }
} // namespace Spacy::Mixin
