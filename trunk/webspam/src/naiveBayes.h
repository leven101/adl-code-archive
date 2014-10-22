#ifndef INC_WEBSPAM_ONLINENB_H
#define INC_WEBSPAM_ONLINENB_H
#include "base.h"

class NaiveBayes {
public:
  NaiveBayes(Parameters& params, vector<vector< float> >& data,
    map<int, LABEL_T>& train_set, map<int, LABEL_T>& test_set): 
    params_(&params), binary_(params_->getBoolValue("load-raw-feats")),  
    data_(data), train_set_(train_set), test_set_(test_set), 
    noTrSpamLbls_(0), noTrNonspamLbls_(0) { 
  }
  ~NaiveBayes() {
  }
  LABEL_T predict(int, float*, float*, bool =false);
  void train();
  void normalizedHistograms(int, LABEL_T);
private:
  Parameters* const params_;
  const bool binary_;
  bool newstuff_;
  vector<vector<float> >& data_;
  map<int, LABEL_T> train_set_, &test_set_;
  float noTrSpamLbls_, noTrNonspamLbls_;
  map<int, map<float, float> > spamFeats_, nonspamFeats_;  //featureID->{value, count of value}
  map<int, int> binSpamFeats_, binNonspamFeats_;
  map<int, pair<float, float> > spamStats_, nonspamStats_;
  void binaryHistograms(int, LABEL_T);
  float normalPredict(int, float*, float*, bool);
  void gaussPredict(int, float*, float*);
  float binaryPredict(int, float*, float*);
  void computeStats();
  float gaussian(float mu, float sigma, float x) {
    float term1 = 1 / sqrt(2* 3.14159 * sigma);
    float term2 = exp(-(pow(x-mu, 2) / (2 * sigma)));
    return term1 * term2;
  }
};
#endif
