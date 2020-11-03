#!/usr/bin/env bash

set -e

mkdir -p /tmp/deps
cd /tmp/deps || exit 1

# Kaskade
tar xzf Kaskade7.4.tar.gz
cd Kaskade7.4 || exit 1
mkdir -p build
cd build || exit 1
conan install ..
cmake .. -GNinja -DCMAKE_TOOLCHAIN_FILE=conan_paths.cmake -DCMAKE_BUILD_TYPE=Release
cmake --build . --target install
cd -

