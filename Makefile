LIBRARIES=-lm -fopenmp 
OPTIM=-O3 -march=native  -Wall -Wextra 
DEBUG=-ggdb
SRC="src"
VERBOSE=-DVERBOSE

CC=gcc

all: driver test

test: test.o bin/libdadac.so
	${CC} test.o -L./bin -ldadac ${LIBRARIES} ${DEBUG} -o test

driver: driver.o bin/libdadac.so
	${CC} driver.o -L./bin -ldadac ${LIBRARIES} ${DEBUG} -o driver

test.o: test.c
	${CC} -c test.c -o test.o ${OPTIM} ${LIBRARIES} ${DEBUG} ${VERBOSE}

driver.o: driver.c
	${CC} -c driver.c -o driver.o ${OPTIM} -fopenmp ${DEBUG}

bin/libdadac.so: bin/dadac.o bin/kdtree.o bin/heap.o bin/vptree.o
	${CC} -shared bin/dadac.o bin/kdtree.o bin/heap.o bin/vptree.o ${DEBUG} ${OPTIM} ${LIBRARIES} -o bin/libdadac.so 

bin/dadac.o: src/dadac.c
	${CC} -c src/dadac.c -o bin/dadac.o ${OPTIM} -fopenmp ${DEBUG} -fpic ${VERBOSE}

bin/kdtree.o: src/kdtree.c
	${CC} -c src/kdtree.c -o bin/kdtree.o ${OPTIM} -fopenmp ${DEBUG} -fpic ${VERBOSE}

bin/vptree.o: src/vptree.c
	${CC} -c src/vptree.c -o bin/vptree.o ${OPTIM} -fopenmp ${DEBUG} -fpic ${VERBOSE}

bin/heap.o: src/heap.c
	${CC} -c src/heap.c -o bin/heap.o ${OPTIM} -fopenmp ${DEBUG} -fpic ${VERBOSE}


clean:
	rm bin/*.so bin/*.o *.o driver test
