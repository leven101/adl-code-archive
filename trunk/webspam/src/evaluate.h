#ifndef INC_WEBSPAM_EVAL_H
#define INC_WEBSPAM_EVAL_H
#include "base.h"

class Evaluate {
public:
  static void Fmeasure(float trueNonspam, float falseNonspam, 
    float falseSpam, float trueSpam) {
    float recall = trueSpam / (falseSpam + trueSpam); // also true positive rate
    float precision = trueSpam / (falseNonspam + trueSpam);
    float falsePosRate = falseNonspam / (falseNonspam + trueNonspam);
    float F = 2 * ((precision * recall) / (precision + recall));
    cerr << "RESULTS: Precision=" << precision << "\tRecall=" << recall << "\tFP%=" << falsePosRate << "\tF=" << F << endl; 
  }
  static void Fmeasure(const map<int, LABEL_T>& realLbls, 
    const map<int, LABEL_T>& predicted) {
    float tp(0), fp(0), tn(0), fn(0); 
    map<int, LABEL_T>::const_iterator ritr;
    iterate(predicted,pitr) {
      ritr = realLbls.find(pitr->first);
      assert(ritr != realLbls.end());
      if(ritr->second == SPAM) {
        if(pitr->second == SPAM) ++tp;
        else ++fn;
      }
      else if(ritr->second == NONSPAM) {
        if(pitr->second == NONSPAM) ++tn;
        else ++fp;
      }
    }
    Evaluate::Fmeasure(tn, fn, fp, tp);
  }
  static void percentCorrect(const map<int, LABEL_T>& realLbls, 
    const map<int, LABEL_T>& predicted) {
    // get sensitivity (#-of-true-pos/total-pos) and specificity (#-of-true-neg/total-neg)
    float truePos(0), trueNeg(0);
    float falsePos(0), falseNeg(0);
    map<int, LABEL_T>::const_iterator mitr;
    iterate(predicted, itr) {
      mitr = realLbls.find(itr->first);
      assert(mitr != realLbls.end());
      if(itr->second == mitr->second) { // then it's correct
        if(itr->second == SPAM) ++truePos;
        else if(itr->second == NONSPAM) ++trueNeg;
      }
      else { // incorrect
        if(itr->second == SPAM) ++falsePos;
        else if(itr->second == NONSPAM) ++falseNeg;
      }
    }
    // return values
    float sens = truePos / (truePos + falseNeg); 
    float spec = trueNeg / (trueNeg + falsePos); 
    Evaluate::Fmeasure(trueNeg, falseNeg, falsePos, truePos);
    cerr << "RESULTS: Sensitivity = " << sens << "\tSpecificity = " << spec << endl;
    cerr << "Number of TP = " << truePos << "\tNumber of TN = " << trueNeg << endl; 
    cerr << "Number of FP = " << falsePos << "\tNumber of FN = " << falseNeg << endl; 
  }
};
#endif
