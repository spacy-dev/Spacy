#pragma once

#include <Spacy/Spaces/ScalarSpace/Real.h>

namespace Spacy
{
    namespace Mixin
    {
        /// %Mixin class for contraction rates.
        class ContractionRate
        {
        public:
            /**
             * @brief Constructor.
             * @param desiredContraction desired contraction rate
             * @param relaxedDesiredContraction relaxed contraction rate
             * @param maximalContraction maximal allowed contraction rate
             */
            explicit ContractionRate( Real desiredContraction = 0.25, Real relaxedDesiredContraction = 0.5,
                                      Real maximalContraction = 0.75 ) noexcept;

            void setContraction( Real contraction ) noexcept;

            void setDesiredContraction( Real desiredContraction ) noexcept;

            void setRelaxedDesiredContraction( Real relaxedDesiredContraction ) noexcept;

            void setMaximalContraction( Real maximalContraction ) noexcept;

            [[nodiscard]] Real getContraction() const noexcept;

            [[nodiscard]] Real getDesiredContraction() const noexcept;

            [[nodiscard]] Real getRelaxedDesiredContraction() const noexcept;

            [[nodiscard]] Real getMaximalContraction() const noexcept;

            [[nodiscard]] bool contractionIsAdmissible() const noexcept;

        private:
            Real contraction_ = 1;
            Real desiredContraction_;
            Real relaxedDesiredContraction_;
            Real maximalContraction_;
        };
    } // namespace Mixin
} // namespace Spacy
