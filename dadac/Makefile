LDFLAGS=-lm -fopenmp 
CFLAGS=-O3 -march=native -Wall -Wextra 
DEBUG=-ggdb
SRC="src"
VERBOSE=-DVERBOSE

# ----------------------------------------------------------------------
# /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ 
# uncomment the following lines to compile against a BLAS implementation 
# you can install libblas from system package manager or set your implementation to openblas, mkl, or anything else
# e.g.

#BLAS_INCLUDE=""
#BLAS_LIB=""
#LDFLAGS += -lyourblas  -L${BLAS_LIB}
#CFLAGS += -DUSE_BLAS -I${BLAS_INCLUDE}




CC=@gcc
COMPILING_INFO=@echo "[\e[36mCOMPILING\e[0m]: $@"
LINKING_INFO=@echo "[\e[34mLINKING\e[0m]: $@"



all: driver test
	@echo "[\e[31mFLAGS\e[0m]: ${CFLAGS}"

test: test.o bin/libdadac.so
	${CC} test.o -L./bin -ldadac ${LDFLAGS} ${DEBUG} -o test
	${LINKING_INFO}

driver: driver.o bin/libdadac.so
	${CC} driver.o -L./bin -ldadac ${LDFLAGS} ${DEBUG} -o driver
	${LINKING_INFO}

test.o: test.c
	${CC} -c test.c -o test.o ${CFLAGS} ${LDFLAGS} ${DEBUG} ${VERBOSE}
	${COMPILING_INFO}

driver.o: driver.c
	${CC} -c driver.c -o driver.o ${CFLAGS} -fopenmp ${DEBUG}
	${COMPILING_INFO}

bin/libdadac.so: bin/dadac.o bin/kdtree.o bin/kdtreeV2.o bin/heap.o bin/vptree.o bin/vptreeV2.o
	${CC} -shared bin/dadac.o bin/kdtree.o bin/heap.o bin/vptree.o bin/vptreeV2.o bin/kdtreeV2.o ${DEBUG} ${CFLAGS} ${LDFLAGS}  -o bin/libdadac.so 
	${LINKING_INFO}

bin/%.o : src/%.c
	${CC} -c $< -o $@ ${CFLAGS} -fopenmp ${DEBUG} -fpic ${VERBOSE} 
	${COMPILING_INFO}

clean:
	rm bin/*.so bin/*.o *.o driver test
