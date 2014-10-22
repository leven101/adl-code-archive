#ifndef INC_WEBSPAM_FEATURES_H
#define INC_WEBSPAM_FEATURES_H
#include "file.h"
#include "params.h"
#include "utils.h"
#include "types.h"

using std::vector;
using std::map;
using std::set;
class Base;
class HostGraph;
enum LABEL_T {SPAM = 1, NONSPAM = -1, UNKNOWN = 0};
class Features {
public:
  Features(Parameters& params): params_(&params) {
    bool doBinary = params_->getBoolValue("load-raw-feats");
    lmFeatures_.resize(114529);
    loadLMFeatures();
    if(doBinary) {
      loadRawContentFeatures();
      loadRawLinkFeatures();
    }
    else {
      loadNormalizedContentFeatures();
      loadNormalizedLinkFeatures(); 
    }
  }
  //~Features();
  friend class Base;
  friend class HostGraph;
private:
  vector<vector<float> > contentFeatures_, linkFeatures_, lmFeatures_;
  set<int> hostsMissingFeats_, smallLMSet;
  Parameters* const params_;
  void loadRawContentFeatures();
  void loadRawLinkFeatures();
  void loadNormalizedContentFeatures();
  void loadNormalizedLinkFeatures();
  void loadHostsMissingFeats();
  void normalizeContentFeatures();
  void normalizeLinkFeatures();
  void normalize(int, vector<vector<float> >&, string);
  void printStuff(); 
  void loadLMFeatures();
  void streamLMFeatures(int);
  float logisticF(float z) {
    z *= -1;
    return 1.0f / (1.0f + exp(z));
  }
};
#endif 
