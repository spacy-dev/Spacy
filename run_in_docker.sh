#!/bin/bash

C_COMPILER=$1
CXX_COMPILER=$2
GCOV=$3
COVERALLS_TOKEN=$4

DEPS=/home/deps
SHARED=/home/shared

mkdir $DEPS && cd $DEPS

apt-get update && apt-get install -y $C_COMPILER $CXX_COMPILER
export CXX=/usr/bin/$CXX_COMPILER
export CC=/usr/bin/$C_COMPILER

git clone https://github.com/google/googletest.git && cd googletest 
mkdir -p build && cd build && cmake .. && cmake --build .
cp -r ../googletest/include/gtest /usr/local/include/
cp -r ../googlemock/include/gmock /usr/local/include/
cp googlemock/gtest/lib*.a /usr/local/lib

cd $SHARED
echo "<$GCOV>"
git checkout master && git pull
mkdir -p build && cd build && rm -rf *
if [ "$GCOV" == "" ]; then
  cmake -DBuildTest=ON ..
else
  cmake -DBuildTest=ON -DCoverage=ON ..
fi
cmake --build .
cd Test && ctest


if [ "$GCOV" != "" ]; then
  lcov --gcov-tool $GCOV --capture --no-external --directory .. --base-directory ../../Spacy --output-file coverage.info
  lcov --remove coverage.info '*/Spacy/Adapter/*' -o coverage_without_adapter.info
  coveralls-lcov --repo-token ${COVERALLS_TOKEN} coverage_without_adapter.info
fi
