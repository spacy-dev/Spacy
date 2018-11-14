#pragma once

// Algorithms
#include <Spacy/Algorithm/ACR/acr.hh>
#include <Spacy/Algorithm/CG/CG.h>
#include <Spacy/Algorithm/CG/LinearSolver.h>
#include <Spacy/Algorithm/CG/TerminationCriteria.h>
#include <Spacy/Algorithm/CG/TerminationCriterion.h>
#include <Spacy/Algorithm/CompositeStep/affineCovariantSolver.hh>
#include <Spacy/Algorithm/DampingFactor.h>
#include <Spacy/Algorithm/LipschitzConstant.h>
#include <Spacy/Algorithm/Newton/Newton.h>
#include <Spacy/Algorithm/Newton/TerminationCriteria.h>
#include <Spacy/Algorithm/TrustRegion/TrustRegionSolver.h>
#include <Spacy/Algorithm/parameter.hh>

// Spaces
#include <Spacy/Spaces/ProductSpace.h>
#include <Spacy/Spaces/RealSpace.h>

// Util
#include <Spacy/Util/Invoke.h>
#include <Spacy/Util/Log.h>
#include <Spacy/Util/Mixins.h>
#include <Spacy/Util/cast.hh>
#include <Spacy/Util/copy.hh>
#include <Spacy/Util/voider.hh>

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
