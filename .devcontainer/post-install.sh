#!/usr/bin/env bash

REPO_ROOT=$(git rev-parse --show-toplevel)
mkdir -p "$REPO_ROOT"/build
cd "$REPO_ROOT"/build || exit 1
conan install .. --build=missing
