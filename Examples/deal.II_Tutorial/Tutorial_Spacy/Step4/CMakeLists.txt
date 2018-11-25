cmake_minimum_required(VERSION 3.5)

FIND_PACKAGE(deal.II 8.4 REQUIRED
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )

DEAL_II_INITIALIZE_CACHED_VARIABLES()

project(deal.II-stacy-tutorial-step-4)

find_package(FunG CONFIG REQUIRED)
find_package(Spacy CONFIG REQUIRED)

# Executable
add_executable(run step4.cpp)
DEAL_II_SETUP_TARGET(run)
target_link_libraries(run FunG::FunG Spacy::Spacy)
