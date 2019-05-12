NODES := 4
MPICC := mpiCC
MPIXX := mpic++
CFLAGS := -std=c++14 -Wall -Wextra
SOURCES := task.cpp
TARGET := run.bin
NUMBERS := 40000

all: clean compile run check

check:
	g++ ${CFLAGS} check.cpp -o check.o
	./check.o

run: 
	mpirun -n ${NODES} ./run.bin ${NUMBERS}

compile: clean
	${MPIXX} ${CFLAGS} ${SOURCES} -o ${TARGET}

clean:
	rm -rf *.o ${TARGET}
