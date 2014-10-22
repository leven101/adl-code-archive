#ifndef goldlabels_solid_h
#define goldlabels_solid_h

#include "header.h"
#include "timeseries.h"

class GoldLabels {
public:
  GoldLabels(Parameters&);
  set<int> eps2Skip_, trEpochs_, tstEpochs_;
  int noTrainEp_, noTestEp_, noClasses_, delay_; 
  const vf_t lbls() { return lbls_; }
  const vf_t orgLbls() { return org_; }
  const size_t totEps() { return lbls_.size();} 
//private:
  Parameters* params_;
  vf_t org_, lbls_;
  bool ss_;
  void load();
  void subsample();
  void stats();
  void buildClasses();
};
#endif 
