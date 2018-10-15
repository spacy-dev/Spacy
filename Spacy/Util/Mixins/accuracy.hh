#pragma once

#include "macros.hh"

#include <Spacy/Spaces/ScalarSpace/Real.h>

GENERATE_MIXIN_HEADER( Real, AbsoluteAccuracy, 1e-15 )

GENERATE_MIXIN_HEADER( Real, RelativeAccuracy, 1e-15 )

GENERATE_MIXIN_HEADER( Real, MinimalAccuracy, 1e-15 )

GENERATE_MIXIN_HEADER( Real, DampingAccuracy, 5e-2 )

GENERATE_MIXIN_HEADER( Real, NormalAccuracy, 0.1 )

GENERATE_MIXIN_HEADER( Real, TangentialAccuracy, 0.1 )

GENERATE_MIXIN_HEADER( Real, FallBackTangentialAccuracy, 0.25 )

#include "undefMacros.hh"
