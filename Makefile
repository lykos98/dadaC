LIBRARIES=-lm -fopenmp 
OPTIM=-O3 -march=native
DEBUG= 

DADAC=bin/clustering.o bin/kdtree.o bin/read_fof_snapshot.o

driver: driver.o 
	#gcc driver.o ${DADAC} -o driver ${LIBRARIES}
	gcc driver.o -L./bin -lclustering  ${LIBRARIES} -o driver
driver.o: driver.c
	gcc -c driver.c -o driver.o ${OPTIM} -fopenmp ${DEBUG}
clean:
	rm driver
	rm driver.o
