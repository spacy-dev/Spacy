#pragma once

// Algorithms
#include <Spacy/Algorithm/ACR/ACR.h>
#include <Spacy/Algorithm/CG/CG.h>
#include <Spacy/Algorithm/CG/LinearSolver.h>
#include <Spacy/Algorithm/CG/TerminationCriteria.h>
#include <Spacy/Algorithm/CG/TerminationCriterion.h>
#include <Spacy/Algorithm/CompositeStep/AffineCovariantSolver.h>
#include <Spacy/Algorithm/DampingFactor.h>
#include <Spacy/Algorithm/LipschitzConstant.h>
#include <Spacy/Algorithm/Newton/Newton.h>
#include <Spacy/Algorithm/Newton/TerminationCriteria.h>
#include <Spacy/Algorithm/TrustRegion/TrustRegionSolver.h>

// Spaces
#include <Spacy/Spaces/ProductSpace.h>
#include <Spacy/Spaces/RealSpace.h>

// Util
#include <Spacy/Util/Cast.h>
#include <Spacy/Util/Copy.h>
#include <Spacy/Util/Invoke.h>
#include <Spacy/Util/Log.h>
#include <Spacy/Util/Mixins.h>
#include <Spacy/Util/Voider.h>

// Interfaces and directly related functionality
#include <Spacy/Derivative.h>
#include <Spacy/DynamicOperator.h>
#include <Spacy/HilbertSpaceNorm.h>
#include <Spacy/Norm.h>
#include <Spacy/c1Functional.hh>
#include <Spacy/c1Operator.hh>
#include <Spacy/c2Functional.hh>
#include <Spacy/functional.hh>
#include <Spacy/inducedScalarProduct.hh>
#include <Spacy/linearOperator.hh>
#include <Spacy/linearSolver.hh>
#include <Spacy/operator.hh>
#include <Spacy/scalarProduct.hh>
#include <Spacy/vector.hh>
#include <Spacy/vectorSpace.hh>
#include <Spacy/zeroVectorCreator.hh>
