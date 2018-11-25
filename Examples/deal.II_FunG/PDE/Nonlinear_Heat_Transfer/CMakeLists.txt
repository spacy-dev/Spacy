cmake_minimum_required(VERSION 3.5)

FIND_PACKAGE(deal.II 8.4 REQUIRED
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )

DEAL_II_INITIALIZE_CACHED_VARIABLES()

project(Spacy-Examples-dealII-Nonlinear-Heat)

list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_CURRENT_SOURCE_DIR}/../../../../cmake)

find_package(FunG CONFIG REQUIRED)
find_package(Spacy CONFIG REQUIRED)

# Executable
add_executable(run nonlinear_heat.cpp)

DEAL_II_SETUP_TARGET(run)

# Target libraries
target_link_libraries(run Spacy::Spacy FunG::FunG)
