#!/bin/bash

clang src/dadac.c src/kdtree.c src/kdtreeV2.c src/heap.c src/vptree.c src/vptreeV2.c -Wall -Wextra -O3 -lm -fopenmp=libomp -I./include -fpic -shared -o bin/libdadac.so 


