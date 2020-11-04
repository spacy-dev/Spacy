cmake_minimum_required(VERSION 3.8)

find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
endif()

option(BuildTest "BuildTest" OFF)
option(Kaskade "Kaskade" OFF)
option(Dolfin "Dolfin" OFF)
option(dealII "dealII" OFF)
option(Logging "Logging" ON)

set(CMAKE_CXX_FLAGS_COVERAGE
    "-O0 -g --coverage" CACHE STRING "Flags used by the CXX compiler during COVERAGE builds" FORCE
)
set(CMAKE_CXX_FLAGS_ASAN
    "-O1 -fsanitize=address -fno-omit-frame-pointer -fno-optimize-sibling-calls"
    CACHE STRING "Flags used by the CXX compiler during ASAN builds" FORCE
)
set(CMAKE_CXX_FLAGS_TSAN
    "-O1 -fsanitize=thread -fno-omit-frame-pointer -fno-optimize-sibling-calls"
    CACHE STRING "Flags used by the CXX compiler during TSAN builds" FORCE
)
set(CMAKE_CXX_FLAGS_UBSAN
    "-O1 -fsanitize=undefined -fno-omit-frame-pointer -fno-optimize-sibling-calls"
    CACHE STRING "Flags used by the CXX compiler during UBSAN builds" FORCE
)

mark_as_advanced (CMAKE_CXX_FLAGS_COVERAGE CMAKE_CXX_FLAGS_ASAN CMAKE_CXX_FLAGS_TSAN CMAKE_CXX_FLAGS_UBSAN)

set(CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}"
    CACHE STRING 
    "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel Coverage ASAN TSAN UBSAN."
    FORCE
)

list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_SOURCE_DIR}/cmake)

if(dealII)
    find_package(deal.II 8.4 QUIET
      HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
      )

  DEAL_II_INITIALIZE_CACHED_VARIABLES()
endif()

project(spacy VERSION 0.4.0 LANGUAGES CXX)

if(Logging)
    add_definitions(-DSPACY_ENABLE_LOGGING)
endif()

if(Dolfin)
  find_package(DOLFIN CONFIG REQUIRED)
  include_directories(${DOLFIN_INCLUDE_DIRS})
  include_directories(SYSTEM ${DOLFIN_3RD_PARTY_INCLUDE_DIRS})
endif()

find_package(Eigen3)
if(Eigen3_FOUND)
  include_directories(${EIGEN3_INCLUDE_DIR})
endif()

add_subdirectory(Spacy)
if(BuildTest)
  add_subdirectory(Test)
endif()

# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
    set(DOC_DIR "doc" CACHE STRING "Output directory for the doxygen documentation")
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(doc
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen" VERBATIM
    )
endif(DOXYGEN_FOUND)

if(Kaskade)
    add_subdirectory(Examples/Kaskade)
endif()