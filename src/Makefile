LIBRARIES=-lm -fopenmp 
OPTIM=-O1 -march=native
DEBUG=-g 
id: read_fof_snapshot.o kdtree.o id.o
	gcc kdtree.o read_fof_snapshot.o id.o -o id ${LIBRARIES} 
	rm *.o

id.o: id.c
	gcc -c id.c -o id.o ${OPTIM} -fopenmp ${DEBUG}

read_fof_snapshot.o: read_fof_snapshot.c
	gcc -c read_fof_snapshot.c -o read_fof_snapshot.o ${OPTIM} 
kdtree.o: kdtree.c
	gcc -c kdtree.c -o kdtree.o ${OPTIM} -fopenmp ${DEBUG}
