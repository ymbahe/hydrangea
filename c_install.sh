#!/bin/bash

# Compile and install C libraries

cd ./hydrangea/clib/src/ckat
./build.sh

cd ../../src/sumbins
./build.sh

