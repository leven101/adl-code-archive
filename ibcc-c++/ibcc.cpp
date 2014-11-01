#include "ibcc.h"

bool debug_ = false;

IBCC::IBCC(Parameters p): params_(p) {
  cout.precision(4);
  cerr.precision(4);
  string fname = params_.getParam("input");
  assert(Utils::fileExists(fname));
  matrixFromFile(fname);
  noTstEps_ = atoi(params_.getParam("test-epochs").c_str());
  noClasses_ = noOutputs_; // fix this! 
  unsprvd_ = params_.getParam("unsupervised") == Parameters::kTrueValue;
  noSrcs_ = data_.size();
  noTrEps_ = data_[0].size() - noTstEps_;
  trainLbls_.assign(goldLbls_.begin(), goldLbls_.begin()+noTrEps_);
  cerr << "Number of sources: " << noSrcs_ << endl;
  cerr << "Number of total epochs: " << goldLbls_.size() << endl;
  cerr << "Number of training epochs: " << noTrEps_ << endl;
  cerr << "Number of test epochs: " << noTstEps_ << endl;
  cerr << "Unsuperised: " << unsprvd_ << endl;
}
void IBCC::initialize() {
  cerr << "Initialising priors..." << endl;
  debug_ = true;
  initRho();
  initNu0();
  initAlpha0();
  debug_ = false;
}
void IBCC::initRho() {
  rho_.resize(goldLbls_.size());
  for(size_t i=0; i < goldLbls_.size(); ++i) {
    rho_[i].resize(noClasses_);
  }
  for(size_t i=0; i < goldLbls_.size(); ++i) {
    for(size_t j=0; j < noClasses_; ++j) {
      if(unsprvd_ || (i >= noTrEps_)) {
        rho_[i][j] = log(0.5);
      }
      else {
        rho_[i][j] = trainLbls_[i] == j ? 0 : Log<double>::zero();
      }
    }
  }
}
void IBCC::initNu0() {
  kappa_.resize(noClasses_);
  nu_.resize(noClasses_); 
  nu0_.resize(noClasses_);
  vector<float> iniVals;
  string nu0 = params_.getParam("nu0");
  if(!nu0.empty()) { // initialize based on user parameters
    Utils::splitToFloat(nu0, iniVals, " ");
    bool bNum = iniVals.size() > 1;
    if(bNum) {
      assert(iniVals.size() == noClasses_);
    } 
    for(int i=0; i < noClasses_; ++i) {
      nu0_[i] = log(iniVals[(bNum ? i : 0)]);
    }
  }
  else { // use count of each class from training labels 
    for(int i=0; i < noClasses_; ++i) {
      nu0_[i] = log(std::count(trainLbls_.begin(), trainLbls_.end(), i));
    }
  }
  if(debug_) {
    for(int i=0; i < noClasses_; ++i) {
      cerr << "\tnu0_[" << i << "]=" << exp(nu0_[i]) << "\t";
    }
    cerr << endl;
  }
}
void IBCC::initAlpha0() {
  // setup matrix dimensions
  pi_.resize(noSrcs_);
  alpha_.resize(noSrcs_);
  alpha0_.resize(noSrcs_);
  for(int i=0; i < noSrcs_; ++i) {
    pi_[i].resize(noClasses_);
    alpha_[i].resize(noClasses_);
    alpha0_[i].resize(noClasses_);
    for(int j=0; j < noClasses_; ++j) {
      pi_[i][j].resize(noOutputs_);
      alpha0_[i][j].resize(noOutputs_);
      alpha_[i][j].resize(noOutputs_);
    }
  }
  string alp0 = params_.getParamValue("alp0");
  if(!alp0.empty()) { // initialize based on user parameters
    vector<float> iniVals;
    Utils::splitToFloat(alp0, iniVals, " ");
    bool bNum = iniVals.size() > 1;
    if(bNum) {
      assert(iniVals.size() == (noOutputs_*noClasses_));
    }
    for(size_t i=0; i < noClasses_; ++i) {
      for(size_t j=0; j < noOutputs_; ++j) {
        alpha0_[0][i][j] = log(iniVals[(bNum ? (j+(i*noClasses_)) : 0)]);
      }
      for(size_t i=1; i < noSrcs_; ++i) {
        alpha0_[i] = alpha0_[0];
      }
    }
  }
  else { // initialize based on weak classifier counts
    for(int i=0; i < noSrcs_ ; ++i) { // for each source
      for(int j=0; j < noTrEps_; ++j) { // for each training epoch
         // CIBCC - count partial predictions
        for(size_t k=0; k < noOutputs_; ++k) {
          float *p = &alpha0_[i][trainLbls_[j]][k]; 
          *p = Log<double>::add(*p, log(data_[i][j][k])); 
        }
        /*
        // SIBCC - count class predictions 
        float *p = &alpha0_[i][trainLbls_[j]][data_[i][j].back()];
        *p = Log<double>::add(*p, 0); // log(1)
        */
      }
    }
  }
  { // add true label counts
    for(int i=0; i < noSrcs_ ; ++i) { // for each source
      for(int j=0; j < noTrEps_; ++j) { // for each training epoch
        // prior is correct answer
        float *p = &alpha0_[i][trainLbls_[j]][trainLbls_[j]]; 
        *p = Log<double>::add(*p, 0); // log(1)
      }
    }
  }
  if(debug_) {
    for(int i=0; i < noSrcs_; ++i) {
      for(int j=0; j < noClasses_; ++j) {
        for(size_t k=0; k < noOutputs_; ++k) {
          cerr << "\talpha0_[" << i << "][" << j << "][" << k << "]=" 
            << exp(alpha0_[i][j][k]) << "\t";
        }
        cerr << endl;
      }
    }
  }
}
void IBCC::train(int noItr) {
  initialize();
  for(int i=0; i < noItr; ++i) {
    cerr << "VB Iteration " << (i+1) << "..." << endl;
    // E-step
    lblCnts(i==0); // get class counts  
    baseCnts(i==0); // get base classifiers counts
    // M-step
    updateNu(); // update class priors 
    expctKappa(); // compute class expectations 
    updateAlpha(); // update base classifiers priors 
    expctPi(); // compute base classifiers expectations
    // E-step
    expctRho();
  }
}
float IBCC::backtest() {
  cerr << "Class probabilities: " << endl;
  for(size_t i=0; i < kappa_.size(); ++i) {
    cerr << i << ": " << exp(kappa_[i]) << "\t";
  }
  cerr << endl;
  for(size_t i=noTrEps_; i < goldLbls_.size(); ++i) {
    cout << goldLbls_[i] << "\t";
    for(int j=0; j < noClasses_; ++j) {
      cout << exp(rho_[i][j]) << "\t"; 
    }
    cout << endl;
  }
  return exp(rho_[goldLbls_.size()-1][1]);
}
void IBCC::updateNu() {
  // add the prior counts nu0 to current counts in nu
  if(debug_) cerr << "updateNu()..." << endl;
  for(int i=0; i < noClasses_; ++i) {
    nu_[i] = Log<double>::add(nu0_[i], nu_[i]);
    if(debug_) cerr << "\tnu_[" << i << "]=" << exp(nu_[i]) << endl;
  }
}
void IBCC::expctKappa() { 
  if(debug_) cerr << "expctKappa()..." << endl;
  double lgsum(Log<double>::zero());
  for(int i=0; i < noClasses_; ++i) {
    lgsum = Log<double>::add(lgsum, nu_[i]);
  }
  double psiSum = Dist::psi(exp(lgsum));
  for(size_t i=0; i < noClasses_; ++i) {
    kappa_[i] = Dist::psi(exp(nu_[i])) - psiSum; 
    if(debug_) cerr << "\tkappa_[" << i << "]=" << exp(kappa_[i]) << endl;
  }
}
void CIBCC::baseCnts(bool init) {
  // counts the partial classifications from the base classifiers 
  if(debug_) {
    cerr << "CIBCC baseCnts()..." << endl;
  }
  // clear the counts
  for(int i=0; i < noSrcs_; ++i) {
    for(int j=0; j < noClasses_; ++j) {
      for(int k=0; k < noOutputs_; ++k) {
        alpha_[i][j][k] = Log<double>::zero();
      }
    }
  }
  int startEp(0);
  if(!unsprvd_) {
    for(int i=0; i < noSrcs_; ++i) {
      for(size_t j=0; j < noTrEps_; ++j) { // training points 
        int gold = trainLbls_[j];
        for(size_t k=0; k < noOutputs_; ++k) {
          double prb = log(data_[i][j][k]);
          alpha_[i][gold][k] = Log<double>::add(alpha_[i][gold][k], prb); 
        }
      }
    }
    startEp = noTrEps_;
  }
  for(int i=0; i < noSrcs_; ++i) {
    for(int j=0; j < noClasses_; ++j) {
      for(int k=0; k < noOutputs_; ++k) {
        double lgsum(Log<double>::zero());
        for(size_t l=startEp; l < goldLbls_.size(); ++l) {
          double lgcnt = log(data_[i][l][k]) + rho_[l][j];
          lgsum = Log<double>::add(lgsum, lgcnt); 
        }
        alpha_[i][j][k] = Log<double>::add(alpha_[i][j][k], lgsum);
      }
    }
  }
  if(debug_) {
    for(int i=0; i < noSrcs_; ++i) {
      for(int j=0; j < noClasses_; ++j) {
        for(int k=0; k < noOutputs_; ++k) {
          cerr << "\talpha_[" << i << "][" << j << "][" << k << 
            "]=" << exp(alpha_[i][j][k]) << "\t";
        }
        cerr << endl;
      }
      cerr << endl;
    }
  }
}
void SIBCC::baseCnts(bool init) {
  // only counts the hard classifications from the base classifiers 
  if(debug_) {
    cerr << "SIBCC baseCnts()..." << endl;
  }
  // clear the counts
  for(int i=0; i < noSrcs_; ++i) {
    for(int j=0; j < noClasses_; ++j) {
      for(int k=0; k < noOutputs_; ++k) {
        alpha_[i][j][k] = Log<double>::zero();
      }
    }
  }
  int startEp(0);
  if(!unsprvd_) {
    for(int i=0; i < noSrcs_; ++i) {
      for(size_t j=0; j < noTrEps_; ++j) { // training points 
        int pred = data_[i][j].back();
        int gold = trainLbls_[j];
        alpha_[i][gold][pred] = Log<double>::add(alpha_[i][gold][pred], 0);
      }
    }
    startEp = noTrEps_;
  }
  for(int i=0; i < noSrcs_; ++i) {
    for(int j=0; j < noClasses_; ++j) {
      for(int k=0; k < noOutputs_; ++k) {
        double lgsum(Log<double>::zero());
        for(size_t l=startEp; l < goldLbls_.size(); ++l) {
          if(data_[i][l].back() == k) {
            lgsum = Log<double>::add(lgsum, rho_[l][j]);
          }
        }
        alpha_[i][j][k] = Log<double>::add(alpha_[i][j][k], lgsum);
      }
    }
  }
  if(debug_) {
    for(int i=0; i < noSrcs_; ++i) {
      for(int j=0; j < noClasses_; ++j) {
        for(int k=0; k < noOutputs_; ++k) {
          cerr << "\talpha_[" << i << "][" << j << "][" << k << 
            "]=" << exp(alpha_[i][j][k]) << "\t";
        }
        cerr << endl;
      }
      cerr << endl;
    }
  }
}
void IBCC::lblCnts(bool init) {
  // get the expected counts of the labels
  if(debug_) {
    cerr << "lblCnts()..." << endl; 
  }
  for(int i=0; i < noClasses_; ++i) {
    if(unsprvd_) {
      nu_[i] = Log<double>::zero();
      for(int j=0; j < noTrEps_; ++j) {
        nu_[i] = Log<double>::add(nu_[i], rho_[j][i]);
      }
    }
    else { // supervised 
      // count the gold labels for the training points
      nu_[i] = log(std::count(trainLbls_.begin(), trainLbls_.end(), i));
      if(init) {
        cerr << exp(nu_[i]) << " train instances of class " << i << endl;
      }
      // get the expectations for the test points 
      for(size_t j=noTrEps_; j < goldLbls_.size(); ++j) {
        nu_[i] = Log<double>::add(nu_[i], rho_[j][i]);
      }
    }
    if(debug_) {
      cerr << "\tnu_[" << i << "]=" << exp(nu_[i]) << endl;
    }
  }
}
void IBCC::updateAlpha() {
  if(debug_) cerr << "updateAlpha()..." << endl;
  for(size_t i=0; i < noSrcs_; ++i) {
    for(size_t j=0; j < noClasses_; ++j) {
      for(size_t k=0; k < noOutputs_; ++k) {
        alpha_[i][j][k] = Log<double>::add(alpha0_[i][j][k], alpha_[i][j][k]);
        if(debug_) {
          cerr << "\talpha_[" << i << "][" << j << "][" << k << "]=" <<
            exp(alpha_[i][j][k]) << "\t";
        }
      }
      if(debug_) {
        cerr << endl;
      }
    }
    if(debug_) {
      cerr << endl;
    }
  }
}
void IBCC::expctPi() {  
  if(debug_) {
    cerr << "expctPi()..." << endl;
  }
  for(int i=0; i < noSrcs_; ++i) {
    for(int j=0; j < noClasses_; ++j) {
      double lgsum(Log<double>::zero());
      for(int k=0; k < noOutputs_; ++k) {
        lgsum = Log<double>::add(lgsum, alpha_[i][j][k]);
      }
      double psiSum = Dist::psi(exp(lgsum));
      for(int k=0; k < noOutputs_; ++k) {
        pi_[i][j][k] = Dist::psi(exp(alpha_[i][j][k])) - psiSum; 
        if(debug_) {
          cerr << "\tpi_[" << i << "][" << j << "][" << k << "]=" <<
            exp(pi_[i][j][k]) << "\t";
        }
      }
      if(debug_) {
        cerr << endl;
      }
    }
    if(debug_) {
      cerr << endl;
    }
  }
}
void CIBCC::expctRho() { 
  if(debug_) {
    cerr << "expctRho()..." << endl;
  }
  for(size_t i=0; i < goldLbls_.size(); ++i) {
    if(unsprvd_ || (i >= noTrEps_)) {
      double lgdenom(Log<double>::zero());
      for(size_t j=0; j < noClasses_; ++j) {
        rho_[i][j] = 0; 
        for(int k=0; k < noSrcs_; ++k) {
          for(size_t l=0; l < noOutputs_; ++l) {
            double prb = data_[k][i][l] * pi_[k][j][l];
            rho_[i][j] += prb;  
          }
        }
        rho_[i][j] += kappa_[j]; 
        lgdenom = Log<double>::add(lgdenom, rho_[i][j]);
      }
      for(size_t j=0; j < noClasses_; ++j) {
        rho_[i][j] -= lgdenom;
      }
    }
    // else do nothing 
    /*if(debug_) {
      for(size_t j=0; j < noClasses_; ++j) {
        cerr << "\trho_[" << i << "][" << j << "]=" 
          << exp(rho_[i][j]) << "\t";
      }
      cerr << endl;
    }*/
  }
}
void SIBCC::expctRho() { 
  if(debug_) cerr << "expctRho()..." << endl;
  for(size_t i=0; i < goldLbls_.size(); ++i) {
    if(unsprvd_ || (i >= noTrEps_)) {
      double lgdenom(Log<double>::zero());
      for(size_t j=0; j < noClasses_; ++j) {
        rho_[i][j] = kappa_[j]; 
        for(int k=0; k < noSrcs_; ++k) {
            int pred = data_[k][i].back(); 
            rho_[i][j] += pi_[k][j][pred];
        }
        lgdenom = Log<double>::add(lgdenom, rho_[i][j]);
      }
      for(size_t j=0; j < noClasses_; ++j) {
        rho_[i][j] -= lgdenom;
      }
    }
    // else do nothing 
    /*(if(debug_) {
      for(size_t j=0; j < noClasses_; ++j) {
        cerr << "\trho_[" << i << "][" << j << "]=" 
          << exp(rho_[i][j]) << "\t";
      }
      cerr << endl;
    }*/
  }
}

void IBCC::matrixFromFile(const string& fname) {
  // load data from file
  FileHandler fin(fname, std::ios::in);
  string line;
  vector<string> v;
  int cursrc(-1);
  while(getline(fin, line)) {
    Utils::trim(line);
    Utils::splitToStr(line, v, "\t");
    if(cursrc == -1) {
      noOutputs_ = v.size() - 4;
      cerr << "noOutputs: " << noOutputs_ << endl;
    }
    int src = atoi(v[0].c_str());
    if(cursrc != src) {
      data_.push_back(twoD_t());
      cursrc = src;
    }
    vf_t ptData(noOutputs_+1);
    for(size_t i=2; i < noOutputs_+3; ++i) {
      ptData[i-2] = atof(v[i].c_str());
    }
    data_[src].push_back(ptData);
    assert((int)data_[src].size() == atoi(v[1].c_str())+1);
    if(src == 0) {
      goldLbls_.push_back(atoi(v.back().c_str()));
    }
  }
  fin.close();
}
