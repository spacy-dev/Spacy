cmake_minimum_required(VERSION 3.1)

FIND_PACKAGE(deal.II 8.4 REQUIRED
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )

DEAL_II_INITIALIZE_CACHED_VARIABLES()

project(dealii-tutorial-step-4)
include_directories(${deal.II_INCLUDE_DIRS})

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
    "Choose the type of build, options are: Debug MinSizeRel Release RelWithDebInfo." FORCE)
endif()

set(CMAKE_CXX_STANDARD 14)

add_executable(run step4.cpp)
DEAL_II_SETUP_TARGET(run)
