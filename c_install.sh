#!/bin/bash

# Compile and install C libraries

cd ./hydrangea/clib/src/HsmlAndProjectYB
./build.sh

cd ../../src/ckat
./build.sh
