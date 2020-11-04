find_package(GTest REQUIRED)
find_package(Threads REQUIRED)

set(CMAKE_CXX_CLANG_TIDY "")

aux_source_directory(Mock MOCK_SRC_LIST)
if(Eigen3_FOUND)
  list(APPEND TEST_SRC_LIST
      Spacy/Adaptivity/Refinement_1D.cpp
      Spacy/Algorithm/ODEInt.cpp
  )
endif()

list(APPEND TEST_SRC_LIST
    Spacy/C1Functional.cpp
    Spacy/C1Operator.cpp
    Spacy/C2Functional.cpp
    Spacy/Functional.cpp
    Spacy/IndefiniteLinearSolver.cpp
    Spacy/LinearOperator.cpp
    Spacy/Operator.cpp
    Spacy/VectorSpace.cpp
    Spacy/Algorithm/DampingFactor.cpp
    Spacy/Algorithm/LipschitzConstant.cpp
    Spacy/Algorithm/CompositeStep/Models.cpp
    Spacy/Algorithm/Newton/Scalar/Newton.cpp
    Spacy/Algorithm/Scalar/FindGlobalMinimizer.cpp
    Spacy/Algorithm/Scalar/Models.cpp
    Spacy/Spaces/ProductSpace.cpp
    Spacy/Spaces/RealSpace.cpp
    Spacy/Spaces/ScalarSpace/C1Operator.cpp
    Spacy/Spaces/ScalarSpace/LinearOperator.cpp
    Spacy/Spaces/ScalarSpace/LinearSolver.cpp
    Spacy/Spaces/ScalarSpace/Real.cpp
    Spacy/Spaces/ProductSpace/C1Operator.cpp
    Spacy/Spaces/ProductSpace/Creator.cpp
    Spacy/Spaces/ProductSpace/Norm.cpp
    Spacy/Spaces/ProductSpace/Operator.cpp
    Spacy/Spaces/ProductSpace/ScalarProduct.cpp
    Spacy/Spaces/ProductSpace/Vector.cpp
    Spacy/Spaces/ProductSpace/VectorSpace.cpp
    Spacy/Spaces/PrimalDualProductSpace/Creator.cpp
    Spacy/Spaces/PrimalDualProductSpace/Vector.cpp
    Spacy/Util/Create.cpp
    Spacy/Util/Base/VectorBase.cpp
    Spacy/Util/Mixins/DecreaseCondition.cpp
)

enable_testing()
set(TEST_NAME ${CMAKE_PROJECT_NAME}-unit-test)
add_executable(${TEST_NAME} ${TEST_SRC_LIST} ${MOCK_SRC_LIST})
target_link_libraries(${TEST_NAME} Spacy::Spacy GTest::GTest GTest::Main Threads::Threads)
target_compile_options(${TEST_NAME} PRIVATE -Wno-unused-result)
add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
