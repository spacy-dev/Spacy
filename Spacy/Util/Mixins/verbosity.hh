#pragma once

#include "MixinConnection.h"

namespace Spacy
{
    namespace Mixin
    {
        /// %Mixin class for verbosity.
        class Verbosity : public MixinConnection< Verbosity >
        {
        public:
            /**
             * @brief Constructor.
             * @param verbosityLevel verbosity level (0 = silent, 1 = print information, 2 = print
             * detailed information)
             */
            explicit Verbosity( unsigned verbosityLevel = 0 ) noexcept;

            /**
             * @brief Enable/disable verbosity.
             * @param verbose true: if verbosityLevel = 0, set verbosityLevel = 1; false: if set
             * verbosityLevel = 0
             */
            void setVerbosity( bool verbose );

            /// Check if verbosityLevel > 0
            bool verbose() const noexcept;

            void setVerbosityLevel( unsigned level ) noexcept;

            /// Access verbosity level.
            unsigned getVerbosityLevel() const noexcept;

            /// update function for observer pattern.
            void update( Verbosity* changedSubject );

        private:
            unsigned verbosityLevel_;
        };
    }
}
