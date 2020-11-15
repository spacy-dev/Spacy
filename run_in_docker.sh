#!/bin/bash

C_COMPILER=$1
CXX_COMPILER=$2
GCOV=$3
COVERALLS_TOKEN=$4
GENERATE_DOCUMENTATION=$5

DEPS=/home/deps
SHARED=/home/shared

mkdir -p $DEPS && cd $DEPS || exit 1

export CXX=$(which $CXX_COMPILER)
export CC=$(which $C_COMPILER)
export PATH="$PATH:$DEPS"

cd $DEPS || exit 1
git clone https://github.com/google/googletest.git
cd googletest || exit 1
mkdir -p build && cd build && cmake .. -DCMAKE_INSTALL_PREFIX=$DEPS && cmake --build . && cmake --build . --target install

cd $SHARED || exit 1
mkdir -p build && cd build && rm -rf *
if [ "$COMPUTE_COVERAGE" == "true" ]; then
  $DEPS/bin/cmake -DBuildTest=ON -DCoverage=ON ..
else
  $DEPS/bin/cmake -DBuildTest=ON ..
fi
$DEPS/bin/cmake --build .
cd Test && ctest

if [ "$GCOV" != "" ]; then
  lcov --gcov-tool "$GCOV" --capture --no-external --directory .. --base-directory ../../Spacy --output-file coverage.info
  lcov --remove coverage.info '*/Spacy/Adapter/*' -o coverage_without_adapter.info
  coveralls-lcov --repo-token "${COVERALLS_TOKEN}" coverage_without_adapter.info
fi

if [ "$GENERATE_DOCUMENTATION" == "true" ]; then
  cd $SHARED/build || exit 1
  $DEPS/bin/cmake --build . --target doc
  ../deploy_doc.sh
fi
