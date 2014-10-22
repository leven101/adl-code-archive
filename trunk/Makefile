ROOTDIR=${HOME}/work
SRC=${ROOTDIR}/src
INC=${SRC}/include
LIBS=-lpthread -lxml2 -lm #-licui18n -licuuc -licudata 
BINDIR=${ROOTDIR}/bin
DEFS=-D_FILE_OFFSET_BITS=64 -D_LARGE_FILES -O3
CCFLAGS=${DEFS} -Wall
SRILM=FIXME

#export LD_RUN_PATH=$LD_RUN_PATH:LIBS
file.o: ${INC}/file.h ${INC}/file.cpp ${INC}/fdstream.h
	g++ -c ${CCFLAGS} ${INC}/file.cpp 
params.o: ${INC}/params.h ${INC}/params.cpp file.o
	g++ -c ${CCFLAGS} ${INC}/params.cpp
vocab.o: ${INC}/vocab.h ${INC}/vocab.cpp
	g++ -c ${CCFLAGS} ${INC}/vocab.cpp
streamcounter.o: streamCounter.h streamCounter.cpp ${INC}/types.h ${INC}/utils.h
	g++ -c ${CCFLAGS} streamCounter.cpp
reuters: reuters.h reuters.cpp ${INC}/xmlParser.h ${INC}/dirList.h ${INC}/utils.h \
         streamCounter.o file.o
	g++ ${CCFLAGS} reuters.cpp streamCounter.o file.o \
        -I/usr/include/libxml2 -o reuters
weeks: weekCounts.cpp file.o ${INC}/utils.h ${INC}/dirList.h
	g++ ${CCFLAGS} weekCounts.cpp file.o -o weeks
countmin: countmin.cpp countmin.h hash.h ${INC}/utils.h ${INC}/file.h
	g++ ${CCFLAGS} countmin.cpp -o cmsketch
hash: ${INC}/hash.h hash.cpp ${INC}/types.h ${INC}/utils.h params.o file.o vocab.o 
	g++ ${CCFLAGS} -I${INC} hash.cpp file.o params.o vocab.o
lossless: losslessLM.cpp params.o file.o vocab.o ${INC}/types.h
	g++ losslessLM.cpp vocab.o file.o params.o \
	-L${SRILM}/lib/i686 -ldstruct -I${INC} \
	-I${SRILM}/include -I${SRILM}/dstruct/src \
	-DUSE_SARRAY_TRIE -DINSTANTIATE_TEMPLATES \
	-o main
phash: ${INC}/hash.h ${INC}/perfectHash.h ${INC}/onlineRLM.h ${INC}/types.h ${INC}/utils.h \
    params.o file.o vocab.o ${INC}/RandLMFilter.h ${INC}/quantizer.h ${INC}/RandLMCache.h
	g++ ${CCFLAGS} -I${INC} file.o params.o vocab.o -o orlm/orlm.exe orlm/onlineRLM.cpp 
bloomier: ${INC}/hash.h ${INC}/types.h ${INC}/utils.h ${INC}/RandLMCache.h \
    params.o file.o vocab.o ${INC}/RandLMFilter.h ${INC}/quantizer.h ${INC}/bloomier.h 
	g++ ${CCFLAGS} -I${INC} bloomier.cpp vocab.o file.o params.o -o bloomier.exe
entropy: ${INC}/hash.h ${INC}/types.h ${INC}/utils.h ${INC}/countmin.h \
    params.o file.o vocab.o 
	g++ ${CCFLAGS} -I${INC} vocab.o file.o params.o  entropy.cpp -o entropy.exe
model1: file.o params.o ${INC}/utils.h vocab.o model1.h 
	g++ ${CCFLAGS} -I${INC} model1.cpp file.o params.o vocab.o -o model1.exe
hmm: file.o params.o ${INC}/utils.h vocab.o model1.h model1.cpp hmm.h
	g++ ${CCFLAGS} -I${INC} hmm.cpp file.o params.o vocab.o \
	  -DUSEMODEL1FORHMM -o hmm.exe
epstuff: file.o params.o ${INC}/dirList.h ${INC}/utils.h
	g++ ${CCFLAGS} -I${INC} epstuff.cpp file.o params.o -o epstuff.exe
multi: ${INC}/hash.h ${INC}/perfectHash.h ${INC}/onlineRLM.h ${INC}/types.h ${INC}/utils.h \
    params.o file.o vocab.o ${INC}/RandLMFilter.h ${INC}/quantizer.h ${INC}/RandLMCache.h \
		${INC}/multiOnlineRLM.h
	g++ ${CCFLAGS} -I${INC} orlm/multiOnlineRLM.cpp file.o params.o vocab.o -o orlm/multiOrlm.exe
ppl: ${INC}/TextStats.h ${INC}/onlineRLM.h
	g++ ${CCFLAGS} -I${INC} orlm/multiOrlmPPL.cpp file.o params.o vocab.o ${INC}/TextStats.cpp -o orlm/multiOrlmPPL.out
itg: vocab.o file.o 
	g++ ${CCFLAGS} -I${INC}  file.o vocab.o itg.cpp -o itg.exe 
mnt: pyp/multiNTs.h vocab.o file.o params.o 
	g++ ${CCFLAGS} -I${INC} -I/Users/ablev/work/cdec/include -L/Users/ablev/work/cdec/include \
 	-lz -lcdec -lutils file.o vocab.o params.o -pthread \
	pyp/mt19937ar.c pyp/multiNTs.cpp -o pyp/multiNTs
#clean: 
#	rm -f *.o *.a ${LIBS}/* ${BINDIR}/*
#hiero: pyp/hiero.h vocab.o file.o params.o 
#	g++ ${CCFLAGS} -I${INC} file.o vocab.o params.o pyp/hiero.cpp pyp/mt19937ar.c -pthread -o pyp/hiero
