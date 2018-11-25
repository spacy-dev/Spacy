cmake_minimum_required(VERSION 3.5)

find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
endif()

option(BuildTest "BuildTest" OFF)
option(Coverage "Coverage" OFF)
option(Kaskade "Kaskade" OFF)
option(Dolfin "Dolfin" OFF)
option(dealII "dealII" OFF)
option(Logging "Logging" ON)

list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_SOURCE_DIR}/cmake)

if(dealII)
    find_package(deal.II 8.4 QUIET
      HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
      )
  include_directories(${deal.II_INCLUDE_DIRS})

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

if(Coverage)
  set(CMAKE_BUILD_TYPE Debug)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -coverage")
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
