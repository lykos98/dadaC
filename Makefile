LIBRARIES=-lm -fopenmp 
OPTIM=-O3 -march=native
DEBUG=-g 

DADAC=bin/clustering.o bin/kdtree.o bin/read_fof_snapshot.o

driver: driver.o 
	gcc driver.o ${DADAC} -o driver ${LIBRARIES}
driver.o: driver.c
	gcc -c driver.c -o driver.o ${OPTIM} -fopenmp ${DEBUG}
