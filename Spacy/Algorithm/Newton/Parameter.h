#pragma once

#include <Spacy/Util/Mixins.h>

namespace Spacy
{
    namespace Newton
    {
        struct Parameter : Mixin::Eps,
                           Mixin::MaxSteps,
                           Mixin::RegularityTest,
                           Mixin::RelativeAccuracy,
                           Mixin::Timer< std::chrono::milliseconds >
        {
        };
    }
}
