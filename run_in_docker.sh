#!/bin/bash

C_COMPILER=$1
CXX_COMPILER=$2
GCOV=$3
COVERALLS_TOKEN=$4
COMPUTE_COVERAGE=$5
GENERATE_DOCUMENTATION=$6

DEPS=/home/deps
SHARED=/home/shared

mkdir -p $DEPS && cd $DEPS

export CXX=$(which $CXX_COMPILER)
export CC=$(which $C_COMPILER)
export PATH="$PATH:$DEPS"

apt install wget

cd $DEPS
wget http://github.com/Kitware/CMake/releases/download/v3.13.1/cmake-3.13.1.tar.gz
tar xzf cmake-3.13.1.tar.gz && cd cmake-3.13.1 && ./bootstrap.sh --prefix=$DEPS && make && make install

cd $DEPS
git clone https://github.com/google/googletest.git
cd googletest
mkdir -p build && cd build && cmake .. -DCMAKE_INSTALL_PREFIX=$DEPS && cmake --build . && cmake --build . --target install

cd $DEPS
git clone https://github.com/eigenteam/eigen-git-mirror.git && cd eigen-git-mirror && mkdir build && cd build && cmake .. -DCMAKE_INSTALL_PREFIX=$DEPS && make install

cd $SHARED
mkdir -p build && cd build && rm -rf *
if [ "$COMPUTE_COVERAGE" == "true" ]; then
  $DEPS/bin/cmake -DBuildTest=ON -DCoverage=ON ..
else
  $DEPS/bin/cmake -DBuildTest=ON ..
fi
$DEPS/bin/cmake --build .
cd Test && ctest

if [ "$COMPUTE_COVERAGE" == "true" ]; then
  lcov --gcov-tool $GCOV --capture --no-external --directory .. --base-directory ../../Spacy --output-file coverage.info
  lcov --remove coverage.info '*/Spacy/Adapter/*' -o coverage_without_adapter.info
  coveralls-lcov --repo-token ${COVERALLS_TOKEN} coverage_without_adapter.info
fi

if [ "$GENERATE_DOCUMENTATION" == "true" ]; then
  cd $SHARED/build
  $DEPS/bin/cmake --build . --target doc
  ../deploy_doc.sh
fi
