ROOTDIR=${HOME}/work/webspam
INC=${ROOTDIR}/src/include
DEFS=-D_FILE_OFFSET_BITS=64 -D_LARGE_FILES #-O3
CCFLAGS=${DEFS} -Wall

file.o: ${INC}/file.h ${INC}/file.cpp ${INC}/fdstream.h
	g++ -c ${CCFLAGS} ${INC}/file.cpp 
params.o: ${INC}/params.h ${INC}/params.cpp file.o
	g++ -c ${CCFLAGS} ${INC}/params.cpp
vocab.o: ${INC}/vocab.h ${INC}/vocab.cpp
	  g++ -c ${CCFLAGS} ${INC}/vocab.cpp
nb: params.o file.o features.h base.h hostGraph.h evaluate.h classify.h
	g++ ${CCFLAGS} -I${INC} -DBIN params.o file.o features.cpp base.cpp \
		hostGraph.cpp naiveBayes.cpp main.cpp -o nb.out
pa: params.o file.o features.h base.h hostGraph.h evaluate.h classify.h vocab.o
	g++ ${CCFLAGS} -I${INC} -I/usr/include/libxml2/ params.o file.o features.cpp base.cpp \
		hostGraph.cpp onlinePA.cpp main.cpp -o pa.out vocab.o -lxml2
html: file.o ${INC}/htmlParser.h vocab.o   
	g++ ${CCFLAGS} -I${INC} -I/usr/include/libxml2/ html.cpp file.o vocab.o -lxml2
