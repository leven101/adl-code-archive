#ifndef tsfeats_solid_h
#define tsfeats_solid_h

#include <algorithm>
#include "goldlabels.h"

class TSFeats {
public:
  TSFeats(Parameters& p, GoldLabels& gls): params_(&p), gldLbls_(gls) {
  }
  int readTSData(); // read timeseries streams
  void extraFeatures();
  threeD_t tsfeats_, xfeats_; // source -> epoch -> features 
  Parameters* params_;
  GoldLabels gldLbls_;
};
void TSFeats::extraFeatures() {
  const size_t nosrcs = tsfeats_.size();
  const size_t noeps = tsfeats_[0].size();
  const size_t nofeats = tsfeats_[0][0].size();
  cerr << "Computing extra features..." << endl;
  xfeats_.resize(nosrcs);
  for(size_t i=0; i < nosrcs; ++i) { // for each source
    xfeats_[i].resize(noeps);
    for(size_t k=0; k < nofeats; ++k) {
      vf_t a(noeps); // array for each feature set
      for(size_t j=0; j < noeps; ++j) { // for each epoch
        a[j] = tsfeats_[i][j][k];
      }
      // get all extra features
      vf_t bb = TS::binary_band(a); // compute its binary band
      vf_t absDiff = TS::diff(a); // compute the absolute vaue of change 
      vf_t smav3 = TS::sma(absDiff, 3); // compute period SMA 
      vector<vf_t> xf(noeps); // store extra features
      for(size_t j=0; j < noeps; ++j) { // for each epoch
        xf[j].push_back(bb[j]);
        xf[j].push_back(absDiff[j]);
        xf[j].push_back(smav3[j]);
      }
      for(size_t j=0; j < noeps; ++j) { // for each epoch
        iterate(xf[j], f) {
          xfeats_[i][j].push_back(*f);
        }
      }
    }
  }
  /*iterate(xfeats_, it) { // each source
    piterate(it, iit) { // each day
      cerr << "[";
      piterate(iit, iiit) { // all features
        cerr << *iiit << " | ";
      }
      cerr << "] ";
    }
    cerr << endl;
  }*/
}
int TSFeats::readTSData() {
  cerr << "Reading in timeseries...";
  const string fname = params_->getParam("ts-file");
  const int skip = atoi(params_->getParam("skip").c_str());
  const int delay = atoi(params_->getParam("delay").c_str());
  int num_ts=-1, numEpochs = StrUtils::numlines(fname) - (skip + delay);
  int curEp(-1), tmp(0);
  vector<string> vline;
  vector<float> vfeats;
  string line;
  FileHandler fin(fname, std::ios::in);
  const int tsidx = atoi(params_->getParam("ts-index").c_str());
  while(getline(fin, line)) {
    if(++tmp <= skip) continue;
    if(++curEp >= numEpochs) break; 
    // skipped subsampled epochs
    if(gldLbls_.eps2Skip_.find(curEp) != gldLbls_.eps2Skip_.end())
      continue;
    Utils::splitToStr(line, vline, " ");
    if(num_ts==-1) { // set up data structures
      assert(tsfeats_.size()==0);
      num_ts = tsidx > -1 ? 1 : vline.size(); 
      vector<vf_t> empty(numEpochs);
      for(int i=0; i < num_ts; ++i) {
        tsfeats_.push_back(empty);
      }
    }
    if(num_ts == 1) { // using a single timeseries 
      Utils::splitToFloat(vline[tsidx], vfeats, "/");
      for(size_t j=0; j < vfeats.size(); ++j) {
        tsfeats_[0][curEp].push_back(vfeats[j]);
      }
    }
    else { // using multiple timeseries
      for(int i=0; i < num_ts; ++i) {
        Utils::splitToFloat(vline[i], vfeats, "/");
        for(size_t j=0; j < vfeats.size(); ++j) {
          tsfeats_[i][curEp].push_back(vfeats[j]);
        }
      }
    }
  }
  fin.close();
  cerr << "Total timeseries streams: " << tsfeats_.size() << endl;
  cerr << "Total data points per stream: " << numEpochs << endl; 
  return tsfeats_.size();
}
#endif
