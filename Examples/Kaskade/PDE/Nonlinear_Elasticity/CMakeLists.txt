cmake_minimum_required(VERSION 3.5)

project(spacy-kaskade-pde-nonlinear-elasticity)

find_package(Kaskade7 REQUIRED)
if(NOT TARGET Spacy::Spacy)
    find_package(Spacy REQUIRED)
endif()

add_executable(${PROJECT_NAME} nonlinear_elasticity.cpp ${Kaskade7_UMFPACK_SRC})
target_link_libraries(${PROJECT_NAME}
    Spacy::Spacy Kaskade7::Kaskade7 ${Kaskade7_UMFPACK_LIB}
)