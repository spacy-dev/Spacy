cmake_minimum_required(VERSION 3.5)

project(spacy-kaskade-pde-nonlinear-heat-transfer)

find_package(Kaskade7 REQUIRED)
if(NOT TARGET Spacy::Spacy)
    find_package(Spacy REQUIRED)
endif()

add_executable(${PROJECT_NAME} nonlinear_heat.cpp ${Kaskade7_UMFPACK_SRC})
target_link_libraries(${PROJECT_NAME}
    Spacy::Spacy Kaskade7::Kaskade7 ${Kaskade7_UMFPACK_LIB}
)