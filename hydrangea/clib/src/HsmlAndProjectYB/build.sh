#!/bin/bash

gcc -fPIC -g -O3 -fopenmp -DDOUBLEPRECISION -c -Wall -lgomp main.c tree.c peano.c \
    allvars.c map.c
gcc -shared -lgomp -lm -lc -o HsmlAndProjectYB.so main.o peano.o tree.o allvars.o map.o
cp HsmlAndProjectYB.so ../../
