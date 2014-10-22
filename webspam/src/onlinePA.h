#ifndef INC_WEBSPAM_ONLINEPA_H
#define INC_WEBSPAM_ONLINEPA_H

#include <numeric>
#include "base.h"
/* Optimization needed for variables: 
 * C_spam, C_nonspam
 * gamma, lambda2, labmda1
 */
class OnlinePAAlg{
public:
  OnlinePAAlg(Parameters& params, vector<vector<float> >& data, 
    map<int, LABEL_T>& trainSet, map<int, LABEL_T>& testSet, 
    HostGraph* graph): params_(&params), data_(data), 
    train_set_(trainSet), test_set_(testSet), graph_(graph) {
    // init weight vectors to zero
    assert(data_.size() > 0);
    weights_.resize(data[0].size(), 0); 
    cumm_weights_ = weights_;
    C_ = atof(params_->getParam("slack").c_str());
    C_s = atof(params_->getParam("c-spam").c_str());
    C_ns = atof(params_->getParam("c-nonspam").c_str());
    assert(C_ > 0);
    gamma_ = atof(params_->getParam("gamma").c_str());
    lambda_ = atof(params_->getParam("lambda").c_str());
    numUpdates_ = 1;
  }
  ~OnlinePAAlg() {
  }
  void train();
  void batchTrain();
  LABEL_T predict(int, float*);
  void optimizeWeights(int, float, bool);
private:
  float C_, numUpdates_; // aggressiveness parameter
  float gamma_, lambda_;
  float C_s, C_ns;
  Parameters* const params_;
  vector<float> weights_, cumm_weights_;
  vector<vector<float> >& data_;
  map<int, LABEL_T> train_set_, &test_set_;
  HostGraph* graph_;
  float getScore(int);
  void updateWeights(const vector<float>&, float, LABEL_T, bool);
  void graphRegularizer(bool);
  void graphRegularizer(int);
  float computeV_i(int);
  float computeMSProd(int);
  float computeM_ij(int, int, int);
  bool innerConvergence();
  bool outerConvergence();
  float biDirLinkWeights(int node); 
  void batchOptimizeWeights(bool);
  void binaryUpdate(int, LABEL_T);  
  void regressionUpdate(int, float, bool);  
  float dotProduct(const vector<float>& x) {
    assert(weights_.size() == x.size());
    float dp = std::inner_product(x.begin(), x.end(), weights_.begin(), 0.0f);
    //cerr << "dp = " << dp << " vs scalar = " << scalar << endl;
    return dp;
  }
  float dotProduct2(const vector<float>& x) {
    assert(cumm_weights_.size() == x.size());
    float scalar(0);
    for(int i=0; i < (int)x.size(); ++i) {
      scalar += ((cumm_weights_[i] / numUpdates_) * x[i]);
    }
    return scalar;
  }
  LABEL_T sign(float n) {
    return n > 0 ? SPAM : NONSPAM;
  }
  float norm(const vector<float>& v) {
    float result(0);
    iterate(v, itr) result += pow(*itr, 2.0f);
    //cerr << "Norm = " <<  sqrt(result) << endl;
    return sqrt(result);
  }
  float getTau(float loss, const vector<float>& v) {
    float normSqrd = pow(norm(v), 2);
    //return loss / normSqrd; 
    //return std::min(C_, loss / normSqrd);
    return loss / (normSqrd + (1.0f / (2.0f*C_))); 
  }
  int numTrainExamples() {
    int numTrain = atoi(params_->getParam("percent-data").c_str()); 
    assert(numTrain <= 100);
    numTrain = int(((float)numTrain / 100.0f) * train_set_.size()); //* data_.size()); 
    cerr << "Training on " << numTrain << " examples.\n";
    return numTrain;
  }
};
#endif
