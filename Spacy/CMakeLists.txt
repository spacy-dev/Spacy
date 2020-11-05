if(Dolfin)
  aux_source_directory(Adapter/FEniCS SRC_LIST)
endif()

if(Eigen3_FOUND)
  aux_source_directory(Adapter/Eigen SRC_LIST)
endif()

if(Kaskade)
    aux_source_directory(Adapter/Kaskade SRC_LIST)
endif()


if(dealII)
    aux_source_directory(Adapter/dealII SRC_LIST)
endif()
list(APPEND SRC_LIST
)
aux_source_directory(Adaptivity SRC_LIST)
aux_source_directory(Adapter/Scalar SRC_LIST)
aux_source_directory(Adapter/Generic SRC_LIST)
aux_source_directory(Spaces/ScalarSpace SRC_LIST)
aux_source_directory(Spaces/ProductSpace SRC_LIST)

add_library(spacy ${SRC_LIST}
    C1Functional.h
    C1Operator.h
    C2Functional.h
    Derivative.h
    Derivative.cpp
    DynamicOperator.h
    ForwardIterator.h
    Functional.h
    HilbertSpaceNorm.h
    HilbertSpaceNorm.cpp
    InducedScalarProduct.h
    InducedScalarProduct.cpp
    LinearOperator.h
    LinearSolver.h
    Norm.h
    Operator.h
    ScalarProduct.h
    Spacy.h
    Vector.h
    VectorSpace.h
    VectorSpace.cpp
    ZeroVectorCreator.h
    ### Algorithm ###
    Algorithm/DampingFactor.h
    Algorithm/DampingFactor.cpp
    Algorithm/LipschitzConstant.h
    Algorithm/LipschitzConstant.cpp
    Algorithm/ODEInt.h
    Algorithm/ODEInt.cpp
    Algorithm/ACR/ACR.h
    Algorithm/ACR/ACR.cpp
    Algorithm/CG/CG.h
    Algorithm/CG/CG.cpp
    Algorithm/CG/LinearSolver.h
    Algorithm/CG/LinearSolver.cpp
    Algorithm/CG/Regularization.h
    Algorithm/CG/RegularizeViaPreconditioner.h
    Algorithm/CG/RegularizeViaPreconditioner.cpp
    Algorithm/CG/TerminationCriteria.h
    Algorithm/CG/TerminationCriteria.cpp
    Algorithm/CG/TerminationCriterion.h
    Algorithm/CG/TriangularStateConstraintPreconditioner.h
    Algorithm/CG/TriangularStateConstraintPreconditioner.cpp
    Algorithm/CompositeStep/AffineCovariantSolver.h
    Algorithm/CompositeStep/AffineCovariantSolver.cpp
    Algorithm/CompositeStep/QuadraticModel.h
    Algorithm/CompositeStep/QuadraticModel.cpp
    Algorithm/Newton/DampingStrategies.h
    Algorithm/Newton/DampingStrategies.cpp
    Algorithm/Newton/Newton.h
    Algorithm/Newton/Newton.cpp
    Algorithm/Newton/Parameter.h
    Algorithm/Newton/TerminationCriteria.h
    Algorithm/Newton/TerminationCriteria.cpp
    Algorithm/Preconditioner/Chebyshev.h
    Algorithm/Preconditioner/Chebyshev.cpp
    Algorithm/Rothe/GenericRothe.h
    Algorithm/Rothe/GenericRothe.cpp
    Algorithm/Scalar/FindGlobalMinimizer.h
    Algorithm/Scalar/FindGlobalMinimizer.cpp
    Algorithm/Scalar/Fmin.h
    Algorithm/Scalar/Models.h
    Algorithm/Scalar/Models.cpp
    Algorithm/TrustRegion/TrustRegionSolver.h
    Algorithm/TrustRegion/TrustRegionSolver.cpp
    ### Util ###
    Util/Cast.h
    Util/CompareVTK.h
    Util/Copy.h
    Util/Exceptions.h
    Util/Exceptions.cpp
    Util/ForEach.h
    Util/Invoke.h
    Util/Log.h
    Util/Log.cpp
    Util/Logger.h
    Util/Mixins.h
    Util/Printer.h
    Util/SmartPointerStorage.h
    Util/Timer.h
    Util/Timer.cpp
    Util/Base/AddArithmeticOperators.h
    Util/Base/FunctionalBase.h
    Util/Base/OperatorBase.h
    Util/Base/OperatorBase.cpp
    Util/Base/VectorBase.h
    Util/Base/VectorBase.cpp
    Util/Mixins/Accuracy.h
    Util/Mixins/Accuracy.cpp
    Util/Mixins/ContractionRate.h
    Util/Mixins/ContractionRate.cpp
    Util/Mixins/DecreaseCondition.h
    Util/Mixins/DecreaseCondition.cpp
    Util/Mixins/Domain.h
    Util/Mixins/Domain.cpp
    Util/Mixins/Eps.h
    Util/Mixins/Eps.cpp
    Util/Mixins/Get.h
    Util/Mixins/Index.h
    Util/Mixins/Index.cpp
    Util/Mixins/IterativeRefinements.h
    Util/Mixins/IterativeRefinements.cpp
    Util/Mixins/MaxSteps.h
    Util/Mixins/MaxSteps.cpp
    Util/Mixins/MixinConnection.h
    Util/Mixins/NumberOfThreads.h
    Util/Mixins/NumberOfThreads.cpp
    Util/Mixins/Range.h
    Util/Mixins/Range.cpp
    Util/Mixins/RegularityTest.h
    Util/Mixins/RegularityTest.cpp
    Util/Mixins/Timer.h
    Util/Mixins/Timer.cpp
    Util/Mixins/Verbosity.h
    Util/Mixins/Verbosity.cpp
)
target_include_directories(spacy PUBLIC
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
)

target_link_libraries(spacy PUBLIC lapacke)


target_compile_features(spacy PUBLIC cxx_std_14)
set_target_properties(spacy PROPERTIES CXX_EXTENSIONS OFF)
#This is required so that the exported target has the name Spacy and not spacy
set_target_properties(spacy PROPERTIES EXPORT_NAME Spacy)
add_library(Spacy::Spacy ALIAS spacy)

if(dealII)
    target_include_directories(spacy PUBLIC ${deal.II_INCLUDE_DIRS})
    DEAL_II_SETUP_TARGET(spacy)
endif()

include(GNUInstallDirs)
set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/Spacy)

install(TARGETS spacy
    EXPORT spacy-targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(DIRECTORY . DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Spacy FILES_MATCHING PATTERN "*.h")
install(DIRECTORY . DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Spacy FILES_MATCHING PATTERN "*.hh")
install(EXPORT spacy-targets
  FILE
    SpacyTargets.cmake
  NAMESPACE
    Spacy::
  DESTINATION
    ${CMAKE_INSTALL_LIBDIR}/cmake/Spacy
)

#Create a ConfigVersion.cmake file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/SpacyConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

configure_package_config_file(${PROJECT_SOURCE_DIR}/cmake/SpacyConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/SpacyConfig.cmake
    INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
)

#Install the config, configversion and custom find modules
install(FILES
#    ${PROJECT_SOURCE_DIR}/cmake/FindSpacy.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/SpacyConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/SpacyConfigVersion.cmake
    DESTINATION ${INSTALL_CONFIGDIR}
)

##############################################
## Exporting from the build tree
#configure_file(${PROJECT_SOURCE_DIR}/cmake/FindSpacy.cmake
#    ${CMAKE_CURRENT_BINARY_DIR}/FindSpacy.cmake
#    COPYONLY)

export(EXPORT spacy-targets FILE ${CMAKE_CURRENT_BINARY_DIR}/SpacyTargets.cmake NAMESPACE Spacy::)
#Register package in user's package registry
export(PACKAGE Spacy)
