#pragma once

#include <boost/any.hpp>

#include <functional>

namespace Spacy
{
    /// @cond
    class Vector;
    /// @endcond

    /// Transfer a vector to a refined and or coarsened grid
    using TransferVector = std::function<void(Vector&)>;

    /// Error indicator for grid modifications.
    using ErrorIndicator = boost::any;

    /**
     * @brief Modify a vector space according to an error indicator.
     * @return function to transfer vectors of the modified vector space to the new discretization
     */
    using SpatialAdaptivity = std::function< TransferVector(const ErrorIndicator& indicator)>;

    /**
     * @brief A simplified interface for error estimation.
     *
     * @param x current iterate
     * @param dx current correction
     * @return true if grid was refined, else false
     */
    using EstimateAndRefine = std::function<bool(const Vector& /*x*/, const Vector& /*dx*/)>;
}
