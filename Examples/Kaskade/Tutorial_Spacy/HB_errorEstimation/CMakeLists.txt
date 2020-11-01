cmake_minimum_required(VERSION 3.5)

project(spacy-kaskade-tutorial-hierarchical-error-estimation)

find_package(Kaskade7 REQUIRED)
if(NOT TARGET Spacy::Spacy)
    find_package(Spacy REQUIRED)
endif()

add_executable(${PROJECT_NAME} heat.cpp poisson.hh ${Kaskade7_UMFPACK_SRC})
target_link_libraries(${PROJECT_NAME}
    Spacy::Spacy Kaskade7::Kaskade7 ${Kaskade7_UMFPACK_LIB}
)

file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/graph)
