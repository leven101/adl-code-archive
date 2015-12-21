ROOTDIR=${HOME}/work
SRC=${ROOTDIR}/src
INC=${SRC}/include
DEFS=-D_FILE_OFFSET_BITS=64 -D_LARGE_FILES -O3
CCFLAGS=${DEFS} -Wall
ALGLIB=${INC}/alglib/src
LIBS=-lpthread -lxml2 -lm -lalglib -lcurl -lboost_regex
INCPATHS=-I${INC} -I${ALGLIB} -L${ALGLIB} -L/usr/local/lib

IBAPI=${SRC}/solid/IBAPI/cpp
IBAPIINC=${IBAPI}/Shared
IBAPISRC=${IBAPI}/PosixSocketClient
IBAPITEST=${IBAPI}/PosixSocketClientTest
IBAPIFLAGS = ${CCFLAGS} -DIB_USE_STD_STRING -Wno-switch
#export LD_RUN_PATH=$LD_RUN_PATH:LIBS
all: main 
#tests: all 
file.o: ${INC}/file.h ${INC}/file.cpp ${INC}/fdstream.h 
	g++ -c ${CCFLAGS} ${INC}/file.cpp 
params.o: ${INC}/params.h ${INC}/params.cpp file.o
	g++ -c ${CCFLAGS} ${INC}/params.cpp 
ibapi.o: ${IBAPITEST}/PosixTestClient.*
	g++ $(IBAPIFLAGS) -I$(IBAPIINC) -o $(IBAPITEST)/SocketBase.o -c $(IBAPISRC)/EClientSocketBase.cpp 
	g++ $(IBAPIFLAGS) -I$(IBAPIINC) -o $(IBAPITEST)/ClientSocket.o -c $(IBAPISRC)/EPosixClientSocket.cpp 
	g++ $(IBAPIFLAGS) -I$(IBAPIINC) -I$(IBAPITEST) -I$(IBAPISRC) -I$(INC) -o $(IBAPITEST)/TestClient.o -c $(IBAPITEST)/PosixTestClient.cpp
	ar rvs $(IBAPITEST)/libibapi.a $(IBAPITEST)/*.o
#ibapi.o ib.h 
main: header.h timeseries.h distributions.h file.o params.o tsreg.h txtreg.h \
		txtfeats.h tsfeats.h curl_theysay.h goldlabels.h 
	g++ ${SRC}/solid/file.o ${SRC}/solid/params.o ${LIBS} ${INCPATHS} \
	${IBAPIFLAGS} -I${IBAPIINC} -I${IBAPITEST} -L${IBAPITEST} -libapi \
	goldlabels.cpp ibcc.cpp main.cpp
