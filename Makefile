NODES := 3
MPICC := mpiCC
MPIXX := mpic++
CFLAGS := -std=c++14 -O3 #-Wall -Wextra
SOURCES := task.cpp
TARGET := run.bin
NUMBERS := 100000

all: clean compile run

run: 
	mpirun -n ${NODES} ./run.bin ${NUMBERS}

compile: clean
	${MPIXX} ${CFLAGS} ${SOURCES} -o ${TARGET}

clean:
	rm -rf *.o ${TARGET} *.txt ppi ppg
