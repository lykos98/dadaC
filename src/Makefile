LIBRARIES=-lm 
OPTIM=-O3 -march=native
DEBUG=-g 

all: clustering.o read_fof_snapshot.o kdtree.o

clustering.o: clustering.c
	gcc -c clustering.c -o ../bin/clustering.o ${OPTIM} -fopenmp ${DEBUG}

read_fof_snapshot.o: read_fof_snapshot.c
	gcc -c read_fof_snapshot.c -o ../bin/read_fof_snapshot.o ${OPTIM} 
kdtree.o: kdtree.c
	gcc -c kdtree.c -o ../bin/kdtree.o ${OPTIM} -fopenmp ${DEBUG}

clean:
	rm ../bin/*.o
