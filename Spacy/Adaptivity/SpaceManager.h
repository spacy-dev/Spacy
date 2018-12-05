#pragma once

#include <Spacy/Adaptivity/SpatialAdaptivity.h>
#include <Spacy/VectorSpace.h>

#include <boost/any.hpp>

#include <set>
#include <functional>
#include <unordered_map>

namespace Spacy
{
    /// @cond
    class Vector;
    /// @endcond

    class SpaceManager
    {
    public:
        using SpaceIndex = VectorSpace::Index;

        /**
         * @brief Add handling of adaptive grid adjustments for space idx
         * @param idx space for which adaptive grid adjustments are registered
         * @param adaptivity implementation of adaptive grid adjustments
         */
        void add( SpaceIndex idx, SpatialAdaptivity adaptivity );

        /**
         * @brief Remove handling of adaptive grid adjustments for space idx
         * @param idx space for which adaptive grid adjustments are registered
         */
        void remove( SpaceIndex idx );

        /// Subscribe vector to automatic grid transfers of its underlying space.
        void subscribe( Vector* v );

        /// Unsubscribe vector from automatic grid transfers of its underlying space.
        void unsubscribe( Vector* v );

        /**
         * @brief Adjust the discretization of space idx according to the indicator errorIndicator
         * @param idx index of the space to be refined
         * @param errorIndicator error indicator describing how the discretization should be
         * adjusted
         */
        void adjustDiscretization( SpaceIndex idx, const ErrorIndicator& errorIndicator );

    private:
        struct AdaptivityAndVectors
        {
            SpatialAdaptivity adaptivity;
            std::set< Vector* > vectors;
        };
        using SubscribedVectors = std::unordered_map< SpaceIndex, AdaptivityAndVectors >;
        using SubscribedVectorsIterator = typename SubscribedVectors::iterator;

        SubscribedVectorsIterator find( SpaceIndex spaceIndex );

        SubscribedVectors subscribedVectors_;
    };

    /**
     * @brief Global space manager for handling grid modifications
     * @return
     */
    SpaceManager& globalSpaceManager();

    /**
     * @brief Combine error estimation and grid refinement.
     * @param estimator error estimator
     * @param index index of the space that should be refined
     */
    template <class Estimator>
    EstimateAndRefine makeEstimateAndRefine(Estimator&& estimator, VectorSpace::Index index)
    {
        return [estimator = std::forward<Estimator>(estimator), index](const Vector& x, const Vector& dx)
        {
            if(estimator.estimateError(x, dx))
            {
                globalSpaceManager().adjustDiscretization(index, estimator.getErrorIndicator());
                return true;
            }

            return false;
        };
    }
}
