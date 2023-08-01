LIBRARIES=-lm -fopenmp 
OPTIM=-O4 -march=native -DUSE_INT32 -DUSE_FLOAT32 
DEBUG=-g 
SRC="src"
VERBOSE=-DVERBOSE

DADAC=bin/clustering.o bin/kdtree.o 

driver: driver.o libclustering.so
	gcc driver.o -L./bin -lclustering  ${LIBRARIES} ${DEBUG} -o driver

driver.o: driver.c
	gcc -c driver.c -o driver.o ${OPTIM} -fopenmp ${DEBUG}

libclustering.so: clustering.o kdtree.o 
	gcc -shared bin/clustering.o bin/kdtree.o ${DEBUG} -o bin/libclustering.so 

clustering.o: src/clustering.c
	gcc -c src/clustering.c -o bin/clustering.o ${OPTIM} -fopenmp ${DEBUG} -fpic ${VERBOSE}

kdtree.o: src/kdtree.c
	gcc -c src/kdtree.c -o bin/kdtree.o ${OPTIM} -fopenmp ${DEBUG} -fpic


clean:

	rm driver
	rm driver.o
	rm bin/*.so
	rm bin/*.o
