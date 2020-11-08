cmake_minimum_required(VERSION 3.1)

project(spacy-examples-scalar-newton)

if(NOT TARGET Spacy::Spacy)
    find_package(Spacy REQUIRED)
endif()

add_executable(${PROJECT_NAME} newton.cpp)
target_link_libraries(${PROJECT_NAME} Spacy::Spacy)
