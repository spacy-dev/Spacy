#!/usr/bin/env bash

mkdir -p "$1"/build-"$2"
cd "$1"/build-"$2" || exit 1
conan install .. -s build_type="$2"
 