ROOTDIR=${HOME}/work/src/workfusion/
SRC=${ROOTDIR}/ibcc-c++
INC=${SRC}/include
DEFS=-D_FILE_OFFSET_BITS=64 -D_LARGE_FILES -O3
CCFLAGS=${DEFS} -Wall
INCPATHS=-I${INC} -L/usr/local/lib

#export LD_RUN_PATH=$LD_RUN_PATH:LIBS
all: main 
#tests: all 
file.o: ${INC}/file.h ${INC}/file.cpp ${INC}/fdstream.h 
	g++ -c ${CCFLAGS} ${INC}/file.cpp 
params.o: ${INC}/params.h ${INC}/params.cpp file.o
	g++ -c ${CCFLAGS} ${INC}/params.cpp 
main: distributions.h file.o params.o 
	g++ ${SRC}/solid/file.o ${SRC}/solid/params.o ${LIBS} ${INCPATHS} \
	${IBAPIFLAGS} -I${IBAPIINC} -I${IBAPITEST} -L${IBAPITEST} -libapi \
	goldlabels.cpp ibcc.cpp main.cpp
