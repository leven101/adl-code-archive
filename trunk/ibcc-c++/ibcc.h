#ifndef ibcc_solid_h 
#define ibcc_solid_h 

#include <algorithm>
#include "distributions.h"
#include "log_add.h"

class IBCC {
public:
  //IBCC(Parameters, threeD_t&, GoldLabels&);
  IBCC(Parameters);
  void train(int);
  float backtest();
protected:
  Parameters params_;
  bool unsprvd_;
  // data: source -> epoch -> features
  threeD_t data_; 
  // pi is confusion matrix/probs. alpha/alpha0 are sufficient stats for pi
  // pi/alpha/alpha0: source -> classes -> probs/counts/priors 
  threeD_t pi_, alpha_, alpha0_; 
  twoD_t rho_; // rho_: epoch -> p(label)
  // kappa is log p(labels). nu/nu0 are priors on kappa
  vf_t goldLbls_, trainLbls_, kappa_, nu_, nu0_;
  double noSrcs_, noTstEps_, noTrEps_, noClasses_, noOutputs_; 
  virtual void baseCnts(bool=false)=0; // equation (1.10)
  virtual void expctRho()=0; // equation (1.7) & (1.8)
private:
  void matrixFromFile(const string&);
  void initialize();
  void expctKappa(); // equation (1.12)
  void expctPi(); // equation (1.15)
  void updateNu();
  void updateAlpha(); 
  void lblCnts(bool=false); // equation (1.9)
  void initRho(); 
  void initAlpha0();
  void initNu0();
};
class SIBCC: public IBCC {
public:
  SIBCC(Parameters p, threeD_t& d, GoldLabels& g): IBCC(p, d, g) {
  }
  SIBCC(Parameters p): IBCC(p) {
  }
private:
  void baseCnts(bool); 
  void expctRho(); 
};

class CIBCC: public IBCC {
public:
  CIBCC(Parameters p, threeD_t& d, GoldLabels& g): IBCC(p, d, g) {
    cerr << "noClasses: " << noClasses_ << "\tnoOutputs: " << noOutputs_ << endl;
    assert(noClasses_ == noOutputs_);
  }
  CIBCC(Parameters p): IBCC(p) {
    assert(noClasses_ == noOutputs_);
  }
private:
  void baseCnts(bool); 
  void expctRho();
};
#endif
