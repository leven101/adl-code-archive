#ifndef INC_WEBSPAM_BASECLASS_H
#define INC_WEBSPAM_BASECLASS_H
#include "features.h"
#include "evaluate.h"
#include "hostGraph.h"

using std::pair;

class Base {
public:
  Base(Parameters& params): params_(&params), numTotNodes_(-1), 
    numTotFeats_(-1) { 
    features_ = new Features(params);
    graph_ = new HostGraph(params, *features_);
    loadLabels();
    selectFeatSet();
    delete features_; // don't need anymore
    features_ = 0;
  }
  virtual ~Base() {
    if(features_)
      delete features_;
    if(graph_)
      delete graph_;
  }
  void printTrainTestFeatures();
  map<int, LABEL_T> train_labels_, test_labels_;
  vector<vector<float> > data_;
  HostGraph* graph_;
protected: 
  Features* features_;
  Parameters* const params_;
  int numTotNodes_, numTotFeats_;
  void loadTestLabels();
private:
  void selectFeatSet(); 
  void loadLabels(); 
  vector<int> fisherDiscriminant(vector<vector<float> >&);
};
#endif
