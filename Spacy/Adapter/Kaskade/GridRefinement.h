#pragma once

#include <Spacy/Adaptivity/SpatialAdaptivity.h>

#include <algorithm>

namespace Spacy
{
    namespace Kaskade
    {
        template <class GridManager>
        struct RefineCells
        {
        public:
            explicit RefineCells(GridManager& gridManager) noexcept
                : gridManager(gridManager)
            {}

            TransferVector operator()(const Spacy::ErrorIndicator& indicator)
            {
                const auto cellIndicator = boost::any_cast<std::vector<bool>>(indicator);
                std::for_each(gridManager.grid().leafGridView().template begin<0>(),
                              gridManager.grid().leafGridView().template end<0>(),
                              [this, &cellIndicator](const auto& cell)
                {
                    if (cellIndicator[gridManager.grid().leafIndexSet().index(cell)])
                        gridManager.mark(1, cell);
                });
                gridManager.adaptAtOnce();

                return [](const Spacy::Vector&){};
            }

        private:
            GridManager& gridManager;
        };

        template <class GridManager>
        RefineCells<GridManager> refineCells(GridManager& gridManager) noexcept
        {
            return RefineCells<GridManager>(gridManager);
        }
    }
}
