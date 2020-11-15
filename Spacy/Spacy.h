#pragma once

// Algorithms
#include <Spacy/Algorithm/ACR/ACR.h>
#include <Spacy/Algorithm/CG/CG.h>
#include <Spacy/Algorithm/CG/LinearSolver.h>
#include <Spacy/Algorithm/CG/TerminationCriteria.h>
#include <Spacy/Algorithm/CG/TerminationCriterion.h>
#include <Spacy/Algorithm/CompositeStep/AffineCovariantSolver.h>
#include <Spacy/Algorithm/DampingFactor.h>
#include <Spacy/Algorithm/Krylov/MinRes.h>
#include <Spacy/Algorithm/LipschitzConstant.h>
#include <Spacy/Algorithm/Newton/Newton.h>
#include <Spacy/Algorithm/Newton/TerminationCriteria.h>
#include <Spacy/Algorithm/Preconditioner/Chebyshev.h>
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

// Interfaces and directly related functionality
#include <Spacy/C1Functional.h>
#include <Spacy/C1Operator.h>
#include <Spacy/C2Functional.h>
#include <Spacy/Derivative.h>
#include <Spacy/DynamicOperator.h>
#include <Spacy/Functional.h>
#include <Spacy/HilbertSpaceNorm.h>
#include <Spacy/InducedScalarProduct.h>
#include <Spacy/LinearOperator.h>
#include <Spacy/LinearSolver.h>
#include <Spacy/Norm.h>
#include <Spacy/Operator.h>
#include <Spacy/ScalarProduct.h>
#include <Spacy/Vector.h>
#include <Spacy/VectorSpace.h>
#include <Spacy/ZeroVectorCreator.h>
