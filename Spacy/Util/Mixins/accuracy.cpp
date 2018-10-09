#include "accuracy.hh"

#include "macros.hh"

GENERATE_MIXIN_SOURCE( Spacy::Real, AbsoluteAccuracy )

GENERATE_MIXIN_SOURCE( Spacy::Real, RelativeAccuracy )

GENERATE_MIXIN_SOURCE( Spacy::Real, MinimalAccuracy )

GENERATE_MIXIN_SOURCE( Spacy::Real, DampingAccuracy )

GENERATE_MIXIN_SOURCE( Spacy::Real, NormalAccuracy )

GENERATE_MIXIN_SOURCE( Spacy::Real, TangentialAccuracy )

GENERATE_MIXIN_SOURCE( Spacy::Real, FallBackTangentialAccuracy )

#include "undefMacros.hh"
