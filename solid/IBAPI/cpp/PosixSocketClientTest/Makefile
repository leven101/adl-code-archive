CXX=g++
CXXFLAGS=-DIB_USE_STD_STRING -Wall -Wno-switch
ROOT_DIR=..
BASE_SRC_DIR=${ROOT_DIR}/PosixSocketClient
INCLUDES=-I${ROOT_DIR}/Shared/ -I${BASE_SRC_DIR} -I/Users/ablev/work/src/include -I /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.9.sdk/usr/include/
TARGET=PosixSocketClientTest

$(TARGET): ../*/*.h ../*/*.cpp 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o EClientSocketBase.o -c $(BASE_SRC_DIR)/EClientSocketBase.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o EPosixClientSocket.o -c $(BASE_SRC_DIR)/EPosixClientSocket.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o PosixTestClient.o -c PosixTestClient.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o Main.o -c Main.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ EClientSocketBase.o EPosixClientSocket.o PosixTestClient.o Main.o

clean:
	rm -f $(TARGET) *.o

