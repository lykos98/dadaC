LIBRARIES=-lm 
OPTIM=-O3 -march=native
DEBUG=-g
VERBOSE=-DVERBOSE

all:libclustering.so 

libclustering.so: clustering.o kdtree.o read_fof_snapshot.o
	gcc -shared ../bin/clustering.o ../bin/kdtree.o ../bin/read_fof_snapshot.o -o ../bin/libclustering.so 
	rm ../bin/*.o

clustering.o: clustering.c
	gcc -c clustering.c -o ../bin/clustering.o ${OPTIM} -fopenmp ${DEBUG} -fpic ${VERBOSE}

read_fof_snapshot.o: read_fof_snapshot.c
	gcc -c read_fof_snapshot.c -o ../bin/read_fof_snapshot.o ${OPTIM} -fpic
kdtree.o: kdtree.c
	mkdir bin
	gcc -c kdtree.c -o ../bin/kdtree.o ${OPTIM} -fopenmp ${DEBUG} -fpic

clean:
	rm ../bin/*.so
	rm ../bin/*.o
