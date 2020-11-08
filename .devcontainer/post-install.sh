#!/usr/bin/env bash

REPO_ROOT=$(git rev-parse --show-toplevel)
mkdir -p "$REPO_ROOT"/build-Release
cd "$REPO_ROOT"/build-Release || exit 1
conan install .. --build=missing
