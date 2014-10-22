#ifndef ibwrapper_solid_h
#define ibwrapper_solid_h

#include "PosixTestClient.h"
#include "file.h"

class IBWrapper {
public:
  void data(threeD_t* retData) { // get historical streams
    FileHandler fcodes("data/qcodes.txt", std::ios::in); 
    string line;
    vector<string> vl, vcodes;
    while(getline(fcodes, line)) {
      Utils::trim(line);
      vcodes.push_back(line);
    }
    fcodes.close();
    threeD_t tsfeats;
    for(size_t i=0; i < vcodes.size(); ++i) {
      PosixTestClient client(ST_REQ_HIST_DATA);
      client.connect();
      while(client.isConnected()) {
        client.processMessages(vcodes[i]);
      }
      twoD_t data = *client.getHistData();
      retData->push_back(data);
    }
  }
  void order(const string action) { // put an order through IB API 
    PosixTestClient client;
    client.connect();
    while(client.isConnected()) {
      client.processMessages(action);
    }
  }
};
#endif
