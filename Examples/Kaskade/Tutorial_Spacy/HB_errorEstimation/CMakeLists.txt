cmake_minimum_required(VERSION 3.5)

project(spacy-kaskade-tutorial-hierarchical-error-estimation)

find_package(Kaskade7 REQUIRED)
find_package(Spacy REQUIRED)
find_package(Boost REQUIRED system timer thread QUIET)

add_executable(run heat.cpp)
target_link_libraries(run 
    Spacy::Spacy 
    Kaskade7::Kaskade7 Kaskade7::umfpack_solve
    Boost::system Boost::timer Boost::thread
)

file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/graph)
