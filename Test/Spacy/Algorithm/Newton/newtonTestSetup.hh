#pragma once

#include <functional>

#include <Spacy/Algorithm/Newton/Newton.h>
#include <Spacy/C1Operator.h>
#include <Spacy/InducedScalarProduct.h>
#include <Spacy/Spaces/RealSpace.h>
#include <Spacy/ZeroVectorCreator.h>

using namespace Spacy;

auto createLocalNewton( C1Operator F, double relativeAccuracy, double eps, unsigned maxSteps )
{
    auto p = Newton::Parameter{};
    p.setRelativeAccuracy( relativeAccuracy );
    p.setEps( eps );
    p.setMaxSteps( maxSteps );

    return [ F, p ]( const Vector& x0 ) { return localNewton( F, x0, p ); };
}

auto createCovariantNewton( C1Operator F, double relativeAccuracy, double eps, unsigned maxSteps )
{
    auto p = Newton::Parameter{};
    p.setRelativeAccuracy( relativeAccuracy );
    p.setEps( eps );
    p.setMaxSteps( maxSteps );
    Space::R.setScalarProduct( InducedScalarProduct( F.linearization( zero( F.domain() ) ) ) );

    return [ F, p ]( const Vector& x0 ) { return covariantNewton( F, x0, p ); };
}

auto createContravariantNewton( C1Operator F, double relativeAccuracy, double eps, unsigned maxSteps )
{
    auto p = Newton::Parameter{};
    p.setRelativeAccuracy( relativeAccuracy );
    p.setEps( eps );
    p.setMaxSteps( maxSteps );

    return [ F, p ]( const Vector& x0 ) { return contravariantNewton( F, x0, p ); };
}
