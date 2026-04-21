#!/bin/bash
mkdir -p build
cd build
cmake -D compiler=intel \
    ..
make