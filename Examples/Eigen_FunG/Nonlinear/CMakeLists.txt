cmake_minimum_required(VERSION 3.1)

project(spacy-eigen-fung-nonlinear-example)

find_package(FunG REQUIRED)
find_package(Eigen3 REQUIRED)

# Include path, link path, library and compiler flags for Spacy
if(NOT TARGET Spacy::Spacy)
    find_package(Spacy REQUIRED)
endif()

# Executable
add_executable(${PROJECT_NAME} main.cpp)

# Target libraries
target_link_libraries(${PROJECT_NAME} Spacy::Spacy Eigen3::Eigen FunG::FunG)
