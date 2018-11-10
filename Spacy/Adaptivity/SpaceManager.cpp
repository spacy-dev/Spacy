#include "SpaceManager.h"

#include <Spacy/vector.hh>

#include <algorithm>
#include <string>
#include <utility>

namespace Spacy
{
    namespace
    {
        std::string notFound(SpaceManager::SpaceIndex spaceIndex)
        {
            return "Space " + std::to_string(spaceIndex) + " not found.";
        }
    }

    void SpaceManager::add(SpaceIndex idx, SpatialAdaptivity adaptivity)
    {
        using std::end;
        if(subscribedVectors_.find(idx) == end(subscribedVectors_))
            subscribedVectors_[idx] = AdaptivityAndVectors{adaptivity, {}};
    }

    void SpaceManager::remove(SpaceIndex idx)
    {
        const auto subscribedVectorsIterator = subscribedVectors_.find(idx);
        using std::end;
        if(subscribedVectorsIterator != end(subscribedVectors_))
            subscribedVectors_.erase(subscribedVectorsIterator);
    }

    void SpaceManager::subscribe(Vector* v)
    {
        using std::end;
        const auto subscribedVectorsIterator = find(v->space().index());
        const auto containsSpace = subscribedVectorsIterator != end(subscribedVectors_);
        if(!containsSpace)
            throw std::runtime_error(notFound(v->space().index()).append(" Can't subscribe vector."));

        subscribedVectorsIterator->second.vectors.insert(v);
    }

    void SpaceManager::unsubscribe(Vector* v)
    {
        using std::end;
        const auto subscribedVectorsIterator = find(v->space().index());
        const auto containsSpace = subscribedVectorsIterator != end(subscribedVectors_);
        if(!containsSpace)
            throw std::runtime_error(notFound(v->space().index()).append(" Can't unsubscribe vector."));

        subscribedVectorsIterator->second.vectors.erase(v);
    }

    void SpaceManager::adjustDiscretization(SpaceIndex idx, const ErrorIndicator& errorIndicator)
    {
        using std::begin;
        using std::end;
        const auto subscribedVectorsIterator = find(idx);
        const auto containsSpace = subscribedVectorsIterator != end(subscribedVectors_);
        if(!containsSpace)
            throw std::runtime_error(notFound(idx).append(" Ignoring error indicator."));

        auto& adjustGrid = subscribedVectorsIterator->second.adaptivity;
        const auto transfer = adjustGrid(errorIndicator);
        const auto transferPtr = [&transfer](Vector* v) { transfer(*v); };

        auto& vectors = subscribedVectorsIterator->second.vectors;
        std::for_each(begin(vectors), end(vectors), transferPtr);
    }

    SpaceManager::SubscribedVectorsIterator SpaceManager::find(unsigned spaceIndex)
    {
        using std::begin;
        using std::end;
        return std::find_if(begin(subscribedVectors_), end(subscribedVectors_),
                            [spaceIndex](const SubscribedVectors::value_type& subscribedVectorsForVectorSpace)
        {
            return spaceIndex == subscribedVectorsForVectorSpace.first;
        });
    }

    SpaceManager& globalSpaceManager()
    {
        static SpaceManager spaceManager;
        return spaceManager;
    }
}
