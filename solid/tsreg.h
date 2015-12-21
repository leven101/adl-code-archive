#ifndef linreg_solid_h
#define linreg_solid_h

#include "ibcc.h"

class TSReg {
public:
  TSReg(const TSFeats* tsf, GoldLabels& gls, threeD_t& wcs): 
    tsf_(tsf), goldLbls_(gls), output_(&wcs) {
    xfeats_=0;
    if(tsf_->xfeats_.size()) 
      xfeats_=tsf_->xfeats_[0][0].size();
    formatData();
    cerr << "Done formatting data." << endl;
  }
  void runLogReg();
  void runRandForest();
  void formatData();
  void runBatch();
  void baselines();
private:
  real_2d_array addTrainingPoints(const int, const int);
  real_2d_array addTrainingPoint(const real_2d_array&, const real_1d_array&, f_t);
  const TSFeats* tsf_;
  // separate models per source 
  vector<real_2d_array> trdata_; // each vector is source with features|gold label
  vector<vector<real_1d_array> > tstdata_; // source -> epochs -> features
  int xfeats_;
  GoldLabels goldLbls_;
  threeD_t* output_; // source -> epochs -> predictions 
  void clear_alglib_matrix(real_2d_array&); 
  void updateOutput(real_1d_array&);
};
void TSReg::formatData() {
  // each row of data is a feature, each column of data is a data point
  const vector<vector<vf_t> >& features = tsf_->tsfeats_;
  int noStreams = features.size();  
  // setup alglib matrices with preset sizes
  int baseFeats = features[0][0].size(); 
  int totFeats = baseFeats + (xfeats_ ? xfeats_ : 0);
  cerr << "Number of features per timeseries: " << totFeats << endl;
  real_2d_array r2d; // training data 
  r2d.setlength(goldLbls_.noTrainEp_, totFeats+1); // features + gold label 
  real_1d_array r1d;
  r1d.setlength(totFeats);
  int nTest = goldLbls_.totEps(); // noTestEp_;  // change this to test all
  vector<real_1d_array> vr1d(nTest); 
  for(int i=0; i < nTest; ++i) {
    vr1d[i] = r1d;
  }
  for(int i=0; i < noStreams; ++i) {
    trdata_.push_back(r2d);
    tstdata_.push_back(vr1d);
  }
  // add all features to alblib training and test arrays
  for(int i=0; i < noStreams; ++i) {
    //cerr << "Processing stream " << i << endl;
    for(int j=0; j < goldLbls_.totEps(); ++j) {
      if(j < goldLbls_.noTrainEp_) {
        for(int k=0; k < baseFeats; ++k) {
          trdata_[i][j][k] = features[i][j][k]; // raw feat 
        }
        if(xfeats_) {
          for(int k=0; k < xfeats_; ++k) {
            trdata_[i][j][k+baseFeats] = tsf_->xfeats_[i][j][k]; 
          }
        }
        trdata_[i][j][totFeats] = goldLbls_.lbls()[j];
      }
      //else {  // comment out this else statement to test all
        real_1d_array tstpt;
        tstpt.setlength(totFeats);
        for(int k=0; k < baseFeats; ++k) {
          tstpt[k] = features[i][j][k];
        }
        if(xfeats_) {
          for(int k=0; k < xfeats_; ++k) {
            tstpt[k+baseFeats] = tsf_->xfeats_[i][j][k]; 
          }
        }
        //tstdata_[i][(j - goldLbls_.noTrainEp_)] = tstpt;
        tstdata_[i][j] = tstpt;
      //}
    }
  }
}
void TSReg::runRandForest() {
  for(size_t i=0; i < trdata_.size(); ++i) {  // process each stream separately     
    /*if((i != 9) && (i != 12) && (i != 15) && (i != 32) && (i != 33)) {
      continue;
    }*/
    output_->push_back(vector<vf_t>());
    ae_int_t info;
    decisionforest df;
    dfreport rep;
    dfbuildrandomdecisionforest(trdata_[i], trdata_[i].rows(), trdata_[i].cols()-1, 
      goldLbls_.noClasses_, 50, 0.2, info, df, rep);
    if(info != 1) {
      cerr << "ERROR: mnltrainh(..) returned error value: " << info << endl;
      exit(-1);
    }
    // predict/output test points in matrix
    for(size_t j=0; j < tstdata_[i].size(); ++j) {
      output_->back().push_back(vf_t());
      real_1d_array pred;
      dfprocess(df, tstdata_[i][j], pred);
      updateOutput(pred);
      /*if((int)j > trdata_[i].rows()) { // retrain model 
        if(j % 20 == 0) {
          trdata_[i] = addTrainingPoints(i, j);
          dfbuildrandomdecisionforest(trdata_[i], trdata_[i].rows(), 
            trdata_[i].cols()-1, goldLbls_.noClasses_, 50, 0.2, info, df, rep);
        }
      }*/
    }
  }
}
void TSReg::runLogReg() {
  // output matrix of predictions - 1 row per stream
  assert(trdata_.size() == tstdata_.size());
  for(size_t i=0; i < trdata_.size(); ++i) {  // process each stream separately 
    /*if((i != 2) && (i != 13) && (i != 17) && (i != 21)) { 
      continue;
    }*/
    output_->push_back(vector<vf_t>());
    logitmodel lm2;
    ae_int_t info;
    mnlreport mr;
    mnltrainh(trdata_[i], trdata_[i].rows(), trdata_[i].cols()-1, 
      goldLbls_.noClasses_, info, lm2, mr);
    if(info != 1) {
      cerr << "ERROR: mnltrainh(..) returned error value: " << info << endl;
      exit(-1);
    }
    // predict/output test points in matrix
    for(size_t j=0; j < tstdata_[i].size(); ++j) {
      output_->back().push_back(vf_t());
      real_1d_array pred;
      mnlprocess(lm2, tstdata_[i][j], pred);
      updateOutput(pred);
      /*if((int)j > trdata_[i].rows()) { // retrain model 
        if(j % 20 == 0) { // retrain monthly
          trdata_[i] = addTrainingPoints(i, j);
          mnltrainh(trdata_[i], trdata_[i].rows(), trdata_[i].cols()-1, 
            goldLbls_.noClasses_, info, lm2, mr);
        }
      }*/
    }
  }
}
real_2d_array TSReg::addTrainingPoints(const int idx, const int tstpt) {
  // add test points and gold labels from end of current training data up to present 
  real_2d_array tmp;
  int curNoRows = trdata_[idx].rows();
  int newNoRows = tstpt;
  int noCols = trdata_[idx].cols();
  tmp.setlength(newNoRows, noCols);
  // copy old training points
  for(int i=0; i < curNoRows; ++i) {
    for(int j=0; j < noCols; ++j) {
      tmp[i][j] = trdata_[idx][i][j];
    }
  }
  // add new training points
  for(int i=curNoRows; i < tstpt; ++i) {
    for(int j=0; j < noCols; ++j) {
      if(j == noCols-1) {
        tmp[i][j] = goldLbls_.lbls()[tstpt];
      }
      else {
        tmp[i][j] = tstdata_[idx][tstpt][j];
      }
    }
  }
  return tmp;
}
real_2d_array TSReg::addTrainingPoint(const real_2d_array& train, 
  const real_1d_array& test, f_t goldlbl) {
  // copy old training data
  real_2d_array tmp;
  int noRows = train.rows() + 1; 
  int noCols = train.cols(); 
  tmp.setlength(noRows, noCols);
  for(int i=0; i < noRows-1; ++i) {
    for(int j=0; j < noCols; ++j) {
      tmp[i][j] = train[i][j];
    }
  }
  // add test point
  for(int i=0; i < test.length(); ++i) {
    tmp[noRows-1][i] = test[i];
  }
  // add gold label
  tmp[noRows-1][noCols-1] = goldlbl;
  return tmp;
}
void TSReg::baselines() {
  // NOTE: delay=0 for this since goldlabels are training data
  cerr << "\nRunning back returns baseline..." << endl;
  //srand48(time(0));
  int last = 0; 
  output_->push_back(vector<vf_t>());
  for(size_t i=0; i < goldLbls_.totEps(); ++i) {
    output_->back().push_back(vf_t());
    float pDown = last ? 0 : 1;  // "back return"
    //float pDown = Utils::dRand(); // random
    //float pDown = 0; // constant up/down
    float pUp = 1-pDown;
    output_->back().back().push_back(pDown);
    output_->back().back().push_back(pUp);
    if(pDown > pUp) {
      output_->back().back().push_back(0);
    }
    else {
      output_->back().back().push_back(1);
    }
    last = goldLbls_.lbls()[i];
  }
}
void TSReg::runBatch() {
  cerr << "Running batch baseline... " << endl;
  // note: doesn't use tsf_->xfeats
  const vector<vector<vf_t> >& matrix = tsf_->tsfeats_;
  // read in data
  const size_t num_tr = goldLbls_.noTrainEp_; 
  //const size_t num_tst = goldLbls_.noTestEp_;
  const size_t streams = matrix.size();
  const size_t noEpochs = matrix[0].size();
  const size_t nofeats = matrix[0][0].size();
  const size_t combNoFeats = (nofeats*streams);
  real_2d_array train;
  train.setlength(num_tr, combNoFeats + 1);  // plus DV
  //vector<real_1d_array> test(num_tst);
  vector<real_1d_array> test(noEpochs);
  real_1d_array tmp;
  tmp.setlength(combNoFeats);
  for(size_t i=0; i < test.size(); ++i)
    test[i] = tmp;
  // setup training/testing matrix
  for(size_t i=0; i < noEpochs; ++i) {
    for(size_t j=0; j <= streams; ++j) {
      if(i < num_tr) {
        if(j < streams) {
          for(size_t k=0; k < nofeats; ++k) {
            train[i][(j*nofeats)+k] = matrix[j][i][k];
          }
        }
        else { // add gold label for epoch i
          train[i][combNoFeats] = goldLbls_.lbls().at(i); 
        }
      }
      if(j < streams) {
        for(size_t k=0; k < nofeats; ++k) {
          test[i][(j*nofeats)+k] = matrix[j][i][k];
        }
      }
    }
  }
  // train the model
  logitmodel lm2;
  mnlreport mr;
  ae_int_t info;
  mnltrainh(train, train.rows(), train.cols()-1, 2, info, lm2, mr);
  if(info==1) 
    cerr << "Batch training complete." << endl; 
  cerr << "Testing model. " << test.size() << " testing points with " <<  
    combNoFeats << " features..." << endl;
  // test the model
  output_->push_back(vector<vf_t>());  // new source
  for(size_t i=0; i < noEpochs; ++i) { 
    output_->back().push_back(vf_t()); // new epoch
    real_1d_array pred;
    mnlprocess(lm2, test[i], pred);
    updateOutput(pred);
    /*// retrain "online" model
    addTrainingPoint(train, test[i], goldLbls_.lbls()[(i+num_tr)]);
    mnltrainh(train, train.rows(), train.cols()-1, 2, info, lm2, mr);*/
  }
}
void TSReg::updateOutput(real_1d_array& pred) {
  float maxPrb(-1), maxIdx(-1);
  for(int k=0; k < pred.length(); ++k) {
    output_->back().back().push_back(pred[k]);
    if(pred[k] > maxPrb) {
      maxIdx = k;
      maxPrb = pred[k];
    }
  }
  output_->back().back().push_back(maxIdx);
}
void TSReg::clear_alglib_matrix(real_2d_array& a) {
  alglib_impl::ae_matrix_clear(const_cast<alglib_impl::ae_matrix*>(a.c_ptr()));
}
#endif 
