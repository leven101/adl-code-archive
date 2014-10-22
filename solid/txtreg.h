#ifndef txtreg_solid_h
#define txtreg_solid_h

#include <iomanip>
#include "log_add.h"
#include "ibcc.h"

class TxtReg { 
public:
  TxtReg(const vector<DocFeats*>& docs, Parameters& p, 
    GoldLabels& gls, threeD_t& weakClfs): 
    docs_(docs), params_(&p), goldLbls_(gls), output_(&weakClfs)  {
    if(goldLbls_.noClasses_ != 2) {
      cerr << "WARNING: Need to setup text regression to handle > 2 classes.\n";
    }
    assert(goldLbls_.noClasses_ == 2);
    if(docs_.size() % goldLbls_.totEps() != 0) {
      cerr << "WARNING: Skipping test epochs with no gold labels.\n";
    }
    counts();
  }
  ~TxtReg() {
    iterate(docs_, it) {
      delete *it;
    }
  }
private:
  const vector<DocFeats*>& docs_;
  vector<vf_t> sntcDist_;
  Parameters* params_;
  GoldLabels goldLbls_;
  threeD_t* output_;  
  void counts();
  vf_t TheySayFeatures(int);
  vf_t FSWireFeatures(int);
  void testEpoch(real_2d_array&, int);
  void printTrainingProbs(real_2d_array&, int, logitmodel&); 
  void results(real_1d_array&);
  void results(vf_t&);
  float maxLikelihood(double, double);
};
vf_t TxtReg::FSWireFeatures(int doc) {
  static float poslast(0), neglast(0);
  vf_t feats(2,0.5);
  float noSnts = docs_[doc]->noSnts();
  float pos(0), neg(0);
  if(noSnts) {
    for(int snt=0; snt < noSnts; ++snt) {
      float f = docs_[doc]->txtFeats_[snt][0];
      pos += f;
      neg += (1 - f);
      //if(f > 0.5) ++pos;
      //else ++neg;
    }
    //feats[0] = neg / (neg+pos);
    //feats[1] = pos / (neg+pos);
    feats[0] = neg - neglast; 
    feats[1] = pos - poslast;
  }
  neglast = neg;
  poslast = pos;
  return feats;
}
vf_t TxtReg::TheySayFeatures(int doc) {
  vf_t feats(2,0.5);
  static float poslast(0), neglast(0);
  float noSnts = docs_[doc]->noSnts();
  if(noSnts) {
    feats[0] = docs_[doc]->pcntNeg();
    feats[1] = docs_[doc]->pcntPos();
  }
  return feats;
  float pos = 0.001, neg = 0.001;
  for(int snt=0; snt < noSnts; ++snt) {
    float norm = docs_[doc]->txtFeats_[snt][4];
    pos += docs_[doc]->txtFeats_[snt][1] * norm;
    pos += docs_[doc]->txtFeats_[snt][3] * norm;
    neg += docs_[doc]->txtFeats_[snt][2] * norm;
  }
  /*float d = pos + neg;
  feats[0] = neg / d;
  feats[1] = pos / d;
  return feats;*/
  if(noSnts) {
    if(doc % goldLbls_.lbls().size() == 0) {
      neglast=0;
      poslast=0;
    }
    else {
      //double ndelta = fabs(neg-neglast);
      //double pdelta = fabs(pos-poslast); 
      feats[0] = neg - neglast; 
      feats[1] = pos - poslast;
    }
  }
  neglast = neg;
  poslast = pos;
  return feats;
}
void TxtReg::testEpoch(real_2d_array& trdata, int row) {
  logitmodel lm2;
  mnlreport mr;
  ae_int_t info;
  // train model
  mnltrainh(trdata, row, trdata.cols()-1, 2, info, lm2, mr);
  if(info != 1) {
    cerr << "ERROR: mnltrainh(..) returned error value: " << info << endl;
    exit(-1);
  }
  if(row == goldLbls_.noTrainEp_) {
    printTrainingProbs(trdata, row, lm2); // 'test' all data
  }
  // test epoch
  real_1d_array r1d, pred;
  r1d.setlength(trdata.cols()-1);
  for(int i=0; i < r1d.length(); ++i) {
    r1d[i] = trdata[row][i];
  }
  mnlprocess(lm2, r1d, pred);
  results(pred);
}
void TxtReg::counts() {
  real_2d_array trdata; 
  int ordTrEp(0), ordTstEp(goldLbls_.noTrainEp_);
  const bool reg = params_->getParam("run-reg") == Parameters::kTrueValue;
  const int delay = atoi(params_->getParam("delay").c_str());
  for(size_t doc=0; doc < docs_.size(); ++doc) {
    const int epoch = doc % (goldLbls_.totEps() + delay); 
    if(epoch == goldLbls_.totEps()) continue; // handles delay 
    if(goldLbls_.eps2Skip_.find(epoch) != goldLbls_.eps2Skip_.end()) {
      continue;
    }
    if(epoch == 0) { // when source changes
      if(reg && (doc != 0)) {
        assert(false); 
        assert(ordTrEp == goldLbls_.noTrainEp_);
        assert(ordTstEp == goldLbls_.totEps());
        // all features now in trData. test (and retrain) model for all test points
        for(size_t i=goldLbls_.noTrainEp_; i < goldLbls_.totEps(); ++i) {
          testEpoch(trdata, i); 
        }
      }
      ordTrEp=0;
      ordTstEp=goldLbls_.noTrainEp_;
      output_->push_back(vector<vf_t>()); // new source
    }
    vf_t feats = TheySayFeatures(doc);
    if(reg) {
      if(trdata.rows() == 0) { // set matrix dimensions
        trdata.setlength(goldLbls_.totEps(), feats.size()+1);
      }
      // add features to matrix 
      int row;
      if(goldLbls_.tstEpochs_.find(epoch) != goldLbls_.tstEpochs_.end()) {
        row = ordTstEp;
        ++ordTstEp;
      }
      else { // training epoch
        row = ordTrEp;
        ++ordTrEp;
      }
      for(size_t i=0; i < feats.size(); ++i) {
        trdata[row][i] = feats[i];
      }
      trdata[row][feats.size()] = goldLbls_.lbls()[epoch];
    }
    else {
      assert(feats.size() == 2);
      results(feats);
    }
  }
  if(reg) { // test 
    // all features now in trData. test (and retrain) model for all test points
    for(size_t i=goldLbls_.noTrainEp_; i < goldLbls_.totEps(); ++i) {
      testEpoch(trdata, i); 
    }
  }
  cerr << "Finished counting text stats." << endl;
}
void TxtReg::printTrainingProbs(real_2d_array& data, int row, logitmodel& lm2) {
  assert(row == (int)goldLbls_.trEpochs_.size());
  for(int i=0; i < row; ++i) {
    real_1d_array r1d, pred;
    r1d.setlength(data.cols()-1);
    for(int j=0; j < data.cols()-1; ++j) {
      r1d[j] = data[i][j];
    }
    mnlprocess(lm2, r1d, pred);
    results(pred);
  }
}
void TxtReg::results(real_1d_array& pred) {
  output_->back().push_back(vf_t()); // next epoch
  for(int k=0; k < pred.length(); ++k) {
    output_->back().back().push_back(pred[k]);
  }
  if(pred[0] > pred[1]) {
    output_->back().back().push_back(0);
  }
  else {
    output_->back().back().push_back(1);
  }
}
void TxtReg::results(vf_t& pred) {
  output_->back().push_back(vf_t()); // next epoch
  for(int k=0; k < pred.size(); ++k) {
    output_->back().back().push_back(pred[k]);
  }
  if(pred[0] > pred[1]) {
    output_->back().back().push_back(0);
  }
  else {
    output_->back().back().push_back(1);
  }
}
float TxtReg::maxLikelihood(double a1, double a2) {
  float maxLike(Log<float>::zero()), maxProb;
  for(float i=0.1; i < 1; i += 0.001) {
    float bdist = Dist::lbetadist(i, a1, a2); 
    if(bdist > maxLike) {
      maxLike = bdist;
      maxProb = i;
    }
  }
  return maxProb;
}
#endif 
