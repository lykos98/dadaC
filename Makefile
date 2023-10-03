LIBRARIES=-lm -fopenmp 
OPTIM=-O4 -march=native  -Wall -Wextra 
DEBUG=-g 
SRC="src"
VERBOSE=-DVERBOSE

DADAC=bin/clustering.o bin/kdtree.o 

driver: driver.o libdadac.so
	gcc driver.o -L./bin -ldadac ${LIBRARIES} ${DEBUG} -o driver

driver.o: driver.c
	gcc -c driver.c -o driver.o ${OPTIM} -fopenmp ${DEBUG}

libdadac.so: dadac.o kdtree.o 
	gcc -shared bin/dadac.o bin/kdtree.o ${DEBUG} ${OPTIM} ${LIBRARIES} -o bin/libdadac.so 

dadac.o: src/dadac.c
	gcc -c src/dadac.c -o bin/dadac.o ${OPTIM} -fopenmp ${DEBUG} -fpic ${VERBOSE}

kdtree.o: src/kdtree.c
	gcc -c src/kdtree.c -o bin/kdtree.o ${OPTIM} -fopenmp ${DEBUG} -fpic


clean:

	rm driver
	rm driver.o
	rm bin/*.so
	rm bin/*.o
