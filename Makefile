NODES := 1
MPICC := mpiCC
MPIXX := mpic++
CFLAGS := -std=c++14 -Wall -Wextra
SOURCES := task.cpp
TARGET := run.bin
NUMBERS := 10

all: clean compile run

run: 
	mpirun -n ${NODES} ./run.bin ${NUMBERS}

compile: clean
	${MPIXX} ${CFLAGS} ${SOURCES} -o ${TARGET}

clean:
	rm -rf *.o ${TARGET}
