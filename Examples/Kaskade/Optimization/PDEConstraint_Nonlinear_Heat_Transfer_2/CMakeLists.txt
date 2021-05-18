#cmake_minimum_required(VERSION 3.2)

project(test_spaces)

set(CMAKE_VERBOSE_MAKEFILE ON)

set(CMAKE_CXX_COMPILER "g++-8")
set(CMAKE_CXX_FLAGS "-std=c++17 -w -DHAVE_NUMA -DHAVE_LIBAMIRA=1 -DENABLE_UG=1 -DUG_USE_NEW_DIMENSION_DEFINES -DINDEXSET_HAS_ITERATORS -DHYPRE_SEQUENTIAL=1")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG ")
set(CMAKE_CXX_FLAGS_DEBUG " -g -O0 -DDEBUG")
#set(CMAKE_CXX_FLAGS_RELEASE " -g -O0 -DDEBUG")

if(NOT CMAKE_CONFIGURATION_TYPE AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

### FILE(GLOB SRCS "*.cpp" "$ENV{HOME}/Kaskade7.4_svn_2019-04-02_blas_lapack/libs/*.o")
FILE(GLOB SRCS "nonlinear_heat.cpp" "$ENV{HOME}/Kaskade7.4_svn_2019-04-02_blas_lapack/libs/mumps_solve.o" "$ENV{HOME}/Kaskade7.4_svn_2019-04-02_blas_lapack/libs/umfpack_solve.o")
FILE(GLOB HDRS "*.hh, *.h")

include_directories($ENV{HOME}/Kaskade7.4_svn_2019-04-02_blas_lapack)
include_directories($ENV{HOME}/Programmierung/ITSOL_2/include)
include_directories($ENV{HOME}/Programmierung/libamiramesh/include)
include_directories(/usr/include/hypre)
include_directories($ENV{HOME}/Programmierung/Dune-2.6/include)
include_directories($ENV{HOME}/Spacy)
include_directories($ENV{HOME}/Programmierung/FunG/include)
include_directories(/usr/include/eigen3)
### include_directories(/usr/include/superlu)
#include_directories($ENV{HOME}/Boost-1.69/include)

set(MYLIBS "-L$ENV{HOME}/Kaskade7.4_svn_2019-04-02_blas_lapack/libs -lkaskade \ ")
set(MYLIBS "${MYLIBS} -L$ENV{HOME}/Programmierung/Dune-2.6/lib -ldunegeometry -ldunegrid -ldunecommon ")
set(MYLIBS "${MYLIBS} -lugS2 -lugS3 -lugL -lmpiseq_seq")
set(MYLIBS "${MYLIBS} -lboost_program_options")
set(MYLIBS "${MYLIBS} -lboost_system -lboost_timer -lboost_thread -lboost_chrono -lpthread")
set(MYLIBS "${MYLIBS} -lumfpack -lamd -ldmumps_seq -lmumps_common_seq -lpord_seq -lmpiseq_seq -lmetis")
#set(MYLIBS "${MYLIBS} -L$ENV{HOME}/Programmierung/ITSOL_2/lib -litsol")
#set(MYLIBS "${MYLIBS} -lHYPRE -lHYPRE_parcsr_ls -lHYPRE_IJ_mv")
#set(MYLIBS "${MYLIBS} -L$ENV{HOME}/Programmierung/libamiramesh/lib -lamiramesh")
set(MYLIBS "${MYLIBS} -lgfortran -llapacke -lcblas -lblas -llapack -lnuma")
#set(MYLIBS "${MYLIBS} -lgfortran -llapacke -lcblas -lblas -llapack_atlas -llapack -lnuma")
set(MYLIBS "${MYLIBS} -L$ENV{HOME}/Spacy/build/Spacy -lspacy")
#message("MYLIBS = ${MYLIBS}")


#set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS ${MYLIBS})

string(STRIP ${MYLIBS} MYLIBS)

add_custom_target(headers SOURCES ${HDRS})
add_executable(nonlinear_heat ${SRCS})
target_link_libraries(nonlinear_heat ${MYLIBS})
