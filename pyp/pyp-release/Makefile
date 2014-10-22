ROOTDIR=${HOME}/work/src/pyp
SRC=${ROOTDIR}/pyp-release/src
INC=${SRC}/include
LIBS=-lpthread -lxml2 -lm #-licui18n -licuuc -licudata 
BINDIR=${ROOTDIR}/bin
DEFS=-D_FILE_OFFSET_BITS=64 -D_LARGE_FILES -O3
CCFLAGS=${DEFS} -Wall

file.o: ${INC}/file.h ${INC}/file.cpp ${INC}/fdstream.h
	g++ -c ${CCFLAGS} ${INC}/file.cpp 
params.o: ${INC}/params.h ${INC}/params.cpp file.o
	g++ -c ${CCFLAGS} ${INC}/params.cpp
vocab.o: ${INC}/vocab.h ${INC}/vocab.cpp
	g++ -c ${CCFLAGS} ${INC}/vocab.cpp

binary: src/multiNTs.h vocab.o file.o params.o 
	g++ ${CCFLAGS} -I${INC} -lz file.o vocab.o params.o -pthread \
	src/mt19937ar.c src/multiNTs.cpp -o multiNTs

block: src/multiNTs.h vocab.o file.o params.o 
	g++ ${CCFLAGS} -I${INC} -I/Users/ablev/work/cdec/include -L/Users/ablev/work/cdec/include \
 	-lz -lcdec -lutils file.o vocab.o params.o -pthread \
	src/mt19937ar.c src/multiNTs.cpp -o multiNTs
