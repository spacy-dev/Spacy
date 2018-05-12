#!/bin/bash

C_COMPILER=$1
CXX_COMPILER=$2
GCOV=$3
COVERALLS_TOKEN=$4
COMPUTE_COVERAGE=$5
GENERATE_DOCUMENTATION=$6

DEPS=/home/deps
SHARED=/home/shared

mkdir $DEPS && cd $DEPS

export CXX=$(which $CXX_COMPILER)
export CC=$(which $C_COMPILER)

git clone https://github.com/google/googletest.git && cd googletest 
mkdir -p build && cd build && cmake .. && cmake --build .
cp -r ../googletest/include/gtest /usr/local/include/
cp -r ../googlemock/include/gmock /usr/local/include/
cp googlemock/gtest/lib*.a /usr/local/lib

cd $SHARED
git checkout master && git pull
mkdir -p build && cd build && rm -rf *
if [ "$COMPUTE_COVERAGE" == "true" ]; then
  cmake -DBuildTest=ON -DCoverage=ON ..
else
  cmake -DBuildTest=ON ..
fi
cmake --build .
cd Test && ctest

if [ "$COMPUTE_COVERAGE" == "true" ]; then
  lcov --gcov-tool $GCOV --capture --no-external --directory .. --base-directory ../../Spacy --output-file coverage.info
  lcov --remove coverage.info '*/Spacy/Adapter/*' -o coverage_without_adapter.info
  coveralls-lcov --repo-token ${COVERALLS_TOKEN} coverage_without_adapter.info
fi

if [ "$GENERATE_DOCUMENTATION" == "true" ]; then
  cmake --build . --target doc
  ../deploy_doc.sh
fi
