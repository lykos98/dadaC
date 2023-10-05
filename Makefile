LIBRARIES=-lm -fopenmp 
OPTIM=-O4 -march=native  -Wall -Wextra 
DEBUG= 
SRC="src"
VERBOSE=-DVERBOSE

CC=gcc

DADAC=bin/clustering.o bin/kdtree.o 

driver: driver.o libdadac.so
	${CC} driver.o -L./bin -ldadac ${LIBRARIES} ${DEBUG} -o driver

driver.o: driver.c
	${CC} -c driver.c -o driver.o ${OPTIM} -fopenmp ${DEBUG}

libdadac.so: dadac.o kdtree.o 
	${CC} -shared bin/dadac.o bin/kdtree.o ${DEBUG} ${OPTIM} ${LIBRARIES} -o bin/libdadac.so 

dadac.o: src/dadac.c
	${CC} -c src/dadac.c -o bin/dadac.o ${OPTIM} -fopenmp ${DEBUG} -fpic ${VERBOSE}

kdtree.o: src/kdtree.c
	${CC} -c src/kdtree.c -o bin/kdtree.o ${OPTIM} -fopenmp ${DEBUG} -fpic


clean:

	rm driver
	rm driver.o
	rm bin/*.so
	rm bin/*.o
