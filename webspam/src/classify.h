#ifndef INC_WEBSPAM_WEAKCLASSIFIER_H
#define INC_WEBSPAM_WEALCLASSIFIER_H
#include "base.h"
#include "onlinePA.h"
#include "naiveBayes.h"

class Classifier: public Base {
public:
  Classifier(Parameters& params, int numCfiers): Base(params),
  numCfiers_(numCfiers) {
    srand(time(NULL));
    //srand(0);
    cf_errs_ = new float[numCfiers_];
    for(int i=0; i < numCfiers_; ++i) {
      cf_errs_[i] = 1;
    }
#ifdef BIN
    classifiers_ = new NaiveBayes*[numCfiers_];
    for(int i=0; i < numCfiers_; ++i) {
      classifiers_[i] = buildNBClassifier(cf_errs_ + i); 
    }
#else
    classifiers_ = new OnlinePAAlg*[numCfiers_];
    for(int i=0; i < numCfiers_; ++i) {
      classifiers_[i] = buildOPAClassifier(cf_errs_ + i); 
    }
#endif
  }
  ~Classifier() {
    for(int i=0; i < numCfiers_; ++i)
      delete classifiers_[i];
    delete[] classifiers_;
    if(cf_errs_) 
      delete[] cf_errs_;
  }
#ifdef BIN
  void testWeightedAvg() {
    loadTestLabels(); // load correct labels here
    getErrorWeights();
    map<int, LABEL_T> results;
    //int testcnt(0);
    iterate(test_labels_, tstItr) { // for each test instance
      float tps(0), tpns(0), ps, pns;
      // for all classifiers get weighted prob of each class
      for(int i=0; i < numCfiers_; ++i) {
        classifiers_[i]->predict(tstItr->first, &ps, &pns, false);
        tps += cf_errs_[i] * ps;
        tpns += cf_errs_[i] * pns;
      }
      //cerr << "tps = " << tps << "\ttpsn = " << tpns << "\tcorrect = " << tstItr->second << endl;
      float scalar = std::max(tps, tpns); // choose highest prob
      results[tstItr->first] = scalar == tps ? SPAM : NONSPAM; 
      //if(++testcnt >= 300) break; 
    }
    Evaluate::percentCorrect(test_labels_, results);
  }
  void testBagging() {
    loadTestLabels(); // load correct labels here
    map<int, LABEL_T> results;
    //int testcnt(0);
    iterate(test_labels_, tstItr) { // for each test instance
      float tps(0), tpns(0), ps, pns;
      // for all classifiers get weighted prob of each class
      for(int i=0; i < numCfiers_; ++i) {
        classifiers_[i]->predict(tstItr->first, &ps, &pns, false);
        tps += ps;
        tpns += pns;
      }
      //cerr << "tps = " << tps << "\ttpsn = " << tpns << "\tcorrect = " << tstItr->second << endl;
      tps /= (float)numCfiers_;
      tpns /= (float)numCfiers_;
      float scalar = std::max(tps, tpns); // choose highest prob
      results[tstItr->first] = scalar == tps ? SPAM : NONSPAM; 
      //if(++testcnt >= 300) break; 
    }
    Evaluate::percentCorrect(test_labels_, results);
  }
#else
  void testWeightedAvg() {
    loadTestLabels(); // load correct labels here
    getErrorWeights();
    map<int, LABEL_T> results;
    int testcnt(0);
    iterate(test_labels_, tstItr) { // for each test instance
      float tps(0), scalar; 
      // for all classifiers get weighted prob of each class
      for(int i=0; i < numCfiers_; ++i) {
        classifiers_[i]->predict(tstItr->first, &scalar);
        tps += cf_errs_[i] * scalar;
      }
      //cerr << "tps = " << tps << "\ttpsn = " << tpns << "\tcorrect = " << tstItr->second << endl;
      results[tstItr->first] = sign(tps); 
      if(++testcnt >= 30) break; 
    }
    Evaluate::percentCorrect(test_labels_, results);
  }
  void testBagging() {
    loadTestLabels(); // load correct labels here
    map<int, LABEL_T> results;
    float scalar;
    //int testcnt(0);
    iterate(test_labels_, tstItr) { // for each test instance
      //cerr << "\nNew test point...\n";
      float all(0);
      int vsp(0), vnsp(0);
      for(int i=0; i < numCfiers_; ++i) {
        LABEL_T pred = classifiers_[i]->predict(tstItr->first, &scalar);
        //cerr << "\tPred label for classifier " << i << " = " << pred << endl; 
        if(pred == SPAM) {
          ++vsp;
        }
        else {
          ++vnsp;
        }
        all += scalar;
      }
      //cerr << "\tscalar= " << scalar << endl;
      scalar = fabs(all / (float)numCfiers_);
      if(vnsp >= vsp) scalar *= -1;
      //cerr << "\tscalar= " << scalar << endl;
      //cerr << "\tspam votes = " << vsp << "  nonspam votes = " << vnsp << endl;
      results[tstItr->first] = sign(scalar); 
      LABEL_T correctLbl = tstItr->second;
      if(scalar < -1 ) scalar = -1;
      if(scalar > 1) scalar = 1;
      cout << (correctLbl == -1 ? 0 : 1) << "  " ;//<< scalar << endl;
      cout << ((scalar * .5) + .5) << endl;
      //if(++testcnt >= 100) break; 
    }
    Evaluate::percentCorrect(test_labels_, results);
  }
#endif
private:
  int numCfiers_;
#ifdef BIN
  NaiveBayes** classifiers_;
#else
  OnlinePAAlg** classifiers_;
#endif
  float* cf_errs_;
  NaiveBayes* buildNBClassifier(float* err) {
    NaiveBayes* nb(0);
    if(numCfiers_ == 0) {
      nb = new NaiveBayes(*params_, data_, train_labels_, test_labels_); 
      nb->train();
      return nb;
    } 
    else {
      map<int, LABEL_T> randSample = selectRandSamp();
      nb = new NaiveBayes(*params_, data_, randSample, test_labels_); 
      nb->train();
      //return nb;
      // now use the rest of the training data to validate and add examples based on error rates
      float ps, pns;
      LABEL_T pred;
      *err = 0;
      std::multimap<float, int> m_errSort; // holds the errors in ascending order
      bool validate(true);
      int inLoop(0);
      while(validate && inLoop++ < 0) {
        iterate(train_labels_, tr) {
          if(randSample.find(tr->first) == randSample.end()) { // for unused training examples
            nb->predict(tr->first, &ps, &pns); // get prediction error
            //cerr << "ps = " << ps << "\tpns = " << pns << endl;
            pred = ps > pns ? SPAM : NONSPAM; 
            if(pred != tr->second) { // if label is correct don't add 
              float lerr = fabs(ps - pns); 
              //cerr << "lerr = " << lerr << endl;
              m_errSort.insert(pair<float, int>(lerr, tr->first)); // get sorted list 
            }
          }
        } // end train_labels_ for loop 
        if(m_errSort.size()) {
          int cnt(0);
          riterate(m_errSort, mmitr) {
            if(++cnt >= 100) break;
            nb->normalizedHistograms(mmitr->second, train_labels_[mmitr->second]);
            //cerr << mmitr->first << "-->" << mmitr->second << endl;
            *err += mmitr->first;
            randSample[mmitr->second] = UNKNOWN; // value doesn't matter
          }
          if(randSample.size() == train_labels_.size()) { //used all training data
            validate = false;
          }
          m_errSort.clear();
        }
        else {
          validate = false;
        }
      } // end while
    }
    return nb;
  }
  OnlinePAAlg* buildOPAClassifier(float* err) {
    OnlinePAAlg* opa;
    map<int, LABEL_T> randSample;
    if(numCfiers_ == 1) {
      opa = new OnlinePAAlg(*params_, data_, train_labels_, test_labels_, graph_);
      opa->train(); // do initial training
      return opa;
    } 
    else {
      randSample = selectRandSamp();
      opa = new OnlinePAAlg(*params_, data_, randSample, test_labels_, graph_);
      opa->train();
      return opa;
      // now use the rest of the training data to validate and add examples based on error rates
      float scalar(0);
      LABEL_T pred;
      *err = 0;
      std::multimap<float, int> m_errSort; // holds the errors in ascending order
      bool validate(true);
      while(validate) {
        iterate(train_labels_, tr) {
          if(randSample.find(tr->first) == randSample.end()) { // for unused training examples
            pred = opa->predict(tr->first, &scalar); // get prediction error
            //cerr << "ps = " << ps << "\tpns = " << pns << endl;
            if(pred != tr->second) { // if label is correct don't add 
              float lerr = fabs(tr->second - scalar); 
              //cerr << "lerr = " << lerr << endl;
              m_errSort.insert(pair<float, int>(lerr, tr->first)); // get sorted list 
            }
          }
        } // end train_labels_ for loop 
        if(m_errSort.size()) {
          int cnt(0);
          riterate(m_errSort, mmitr) {
            if(++cnt >= 500) break;
            opa->optimizeWeights(mmitr->second, train_labels_[mmitr->second], true); 
            //cerr << mmitr->first << "-->" << mmitr->second << endl;
            *err += mmitr->first;
            randSample[mmitr->second] = UNKNOWN; // value doesn't matter
          }
          if(randSample.size() == train_labels_.size()) { //used all training data
            validate = false;
          }
          m_errSort.clear();
        }
        else {
          validate = false;
        }
      } // end while
    }
    return opa;
  }
  map<int, LABEL_T> selectRandSamp() {
    /* // use to test ERUS with subset % of data
    int numTrain = numTrainExamples(); 
    int cnt(0);
    map<int, LABEL_T> trSubset;
    iterate(train_labels_, tr) {
      trSubset[tr->first] = tr->second;
      if(++cnt >= numTrain) break;
    }*/
    static bool rev(false);
    rev = !rev;
    const int targ(135);
    int nsp(0), sp(0);
    map<int, LABEL_T> trSet;
    if(rev) {
      //riterate(trSubset, tr) {
      riterate(train_labels_, tr) {
        float chance = Utils::rand<unsigned>(RAND_MAX); 
        chance /= (float)RAND_MAX;
        if(tr->second == SPAM && (sp < targ)) { 
          ++sp;
          trSet.insert(*tr);
        }
        else if(nsp < targ) {
          if(chance > .85) {
            ++nsp;
            trSet.insert(*tr);
          }
        }
      }
    }
    else {
      //iterate(trSubset, tr) {
      iterate(train_labels_, tr) {
        float chance = Utils::rand<unsigned>(RAND_MAX); 
        chance /= (float)RAND_MAX;
        if(tr->second == SPAM && (sp < targ)) { 
          ++sp;
          trSet.insert(*tr);
        }
        else if(nsp < targ) {
          if(chance > .85) {
            ++nsp;
            trSet.insert(*tr);
          }
        }
      }
    }
//    cerr << "Example set size = " << trSet.size() << endl;
    return trSet;
  }
  LABEL_T sign(float n) {
    return n > 0 ? SPAM : NONSPAM;
  }
  void getErrorWeights() {
    float esum(0);
    for(int i=0; i < numCfiers_; ++i) {
      cf_errs_[i] /= train_labels_.size();
      esum += exp(-cf_errs_[i]);
    }
    for(int i=0; i < numCfiers_; ++i) {
      cf_errs_[i] = exp(-cf_errs_[i]) / esum;
      cerr << "cf_err_[" << i << "] = " << cf_errs_[i] << endl;
    }

  }
  int numTrainExamples() {
    int numTrain = atoi(params_->getParam("percent-data").c_str()); 
    assert(numTrain <= 100);
    numTrain = int(((float)numTrain / 100.0f) * train_labels_.size()); //* data_.size()); 
    cerr << "Training on " << numTrain << " examples.\n";
    return numTrain;
  }
};
#endif
