#include "tsfeats.h"
#include "txtfeats.h"
#include "tsreg.h"
#include "txtreg.h"
#include "ib.h"

using namespace std;
// paramter definitions
const ParamDefs paramdefs[] = {
  //name, default value, abbrev, type, msg
  {"top-dir", "~/work/solid/data/text/nfp/sources/*/*/*", "d", Parameters::kStringValue, ""},
  //{"top-dir", "~/work/solid/data/text/fswire/data/gold/*/", "d", Parameters::kStringValue, ""},
  {"ext-name", "*.feats", "e", Parameters::kStringValue, ""},
  {"test-epochs", "1", "t", Parameters::kIntValue, ""},
  {"skip", "625", "", Parameters::kIntValue, ""}, // 13471 for corn, 9538 for gold, 5684 for crude, 625 for nfp data
  {"ts-file", "./data/ts.nsa.all", "tsf", Parameters::kStringValue, ""},
  {"ts-index", "-1", "tsi", Parameters::kIntValue, "column in ts-file to read in (zero-indexed)"},
  {"delay", "1", "", Parameters::kIntValue, "epochs to delay prediction by"},
  {"run-reg", Parameters::kFalseValue, "reg", Parameters::kBoolValue, ""},
  {"no-ss", Parameters::kTrueValue, "ns", Parameters::kBoolValue, ""},
  {"vb-itr", "7", "vb", Parameters::kIntValue, ""},
  {"unsupervised", Parameters::kFalseValue, "unsup", Parameters::kBoolValue, ""},
  {"nu0", "", "", Parameters::kFloatValue, ""},
  {"alp0", "", "", Parameters::kFloatValue, ""},
  {"interval", Parameters::kFalseValue, "int", Parameters::kBoolValue, "do interval robustness testing"},
  {"open", Parameters::kFalseValue, "", Parameters::kBoolValue, "predict open versus close prices"},
};
const string tradingDatesPath = "/Users/ablev/work/solid/data/timeseries/commodities/FUTURE_GC1.csv"; 
void rearrange(GoldLabels& gls, Parameters& params, TSFeats* tsf=NULL) {
  // assumes gls.tstEpochs_ and gls.trEpochs_ are set correctly
  cerr << "Rearranging data for current robustness test..." << endl;
  cerr << "WARNING: orgLbls have no delay\n";
  gls.lbls_ = gls.orgLbls();
  if(*gls.tstEpochs_.rbegin() != gls.lbls_.size()) { // else already in correct order
    size_t streams = tsf ? tsf->tsfeats_.size() : 1;
    for(size_t i=0; i < streams; ++i) {
      set<int>::iterator it = gls.tstEpochs_.begin();
      for(; it != gls.tstEpochs_.end() ; ++it) { // add test epochs to end of tsf.tsfeats_
        if(tsf)
          tsf->tsfeats_[i].push_back(tsf->tsfeats_[i][*it]);
        if(i==0) {
          gls.lbls_.push_back(gls.lbls_.at(*it));
        }
      }
      set<int>::reverse_iterator rit = gls.tstEpochs_.rbegin();
      for(; rit != gls.tstEpochs_.rend(); ++rit) { // remove test epochs from training phase
        if(tsf)
          tsf->tsfeats_[i].erase(tsf->tsfeats_[i].begin() + *rit);
        if(i==0) {
          gls.lbls_.erase(gls.lbls_.begin() + *rit);
        }
      }
    }
  }
  const int delay = atoi(params.getParam("delay").c_str());
  gls.lbls_.erase(gls.lbls_.begin(), gls.lbls_.begin()+delay);
  cerr << "Done rearranging data." << endl;
}
void baselines(Parameters& params, threeD_t& weakClfs, 
  GoldLabels& gls, bool rearr=false) {
  TSFeats tsf(params, gls);
  tsf.readTSData();
  if(rearr) rearrange(gls, params, &tsf);
  TSReg rr(&tsf, gls, weakClfs);
  rr.runBatch();
  //rr.baselines();
}
void getTxtFeats(Parameters& params) {
  vector<DocFeats*> vtxf;
  DirList dir(params.getParam("top-dir"), "*.tok"); 
  vector<string> vs;
  iterate(dir.cFileList_, itr) { // for each source file 
    cerr << "Processing file " << *itr << endl;
    if(Utils::fileExists(*itr + ".feats")) {
      cerr << "WARNING: feature file already exists. Skipping file." << endl;
      continue;
    }
    DocFeats txf(*itr, false, params);  // get its features
    //if(*itr != dir.cFileList_.back())
      //sleep(60);
  }
}
set<string> getTradingDates() {
  // read in FUTURE_GC1.csv and get all trading dates from 2013-01-01 onwards
  set<string> trDates;
  FileHandler fdates(tradingDatesPath, ios::in);
  string line;
  vector<string> vtmp;
  bool keep=false;
  while(getline(fdates, line)) {
    if(!keep) {
      if(line.find("2013-01") == 0) {
        keep=true;
      }
    }
    if(keep) {
      Utils::splitToStr(line, vtmp, ",");
      trDates.insert(vtmp[0]);
    }
  }
  fdates.close();
  return trDates;
}
string extractDateFromPath(const string& path) {
  vector<string> vtmp;
  // extract date from file name
  Utils::splitToStr(path, vtmp, "/"); 
  time_t date = atoi(vtmp.back().substr(0,10).c_str());
  struct tm * tinfo = localtime(&date);
  char buff[15];
  strftime(buff, 15, "%F", tinfo);
  return string(buff);
}
void runTxtReg(Parameters& params, GoldLabels& gls, 
  threeD_t& weakClfs, bool rearr=false) {
  if(rearr) rearrange(gls, params);
  //set<string> trDates = getTradingDates(); // use for fswire
  vector<DocFeats*> vtxf;
  DirList dir(params.getParam("top-dir"), params.getParam("ext-name")); 
  vector<string> vs;
  cerr << "Loading data from disk..." << endl; 
  const vector<string>& cfl = dir.cFileList_;
  int numZero(0);
  for(size_t i=0; i < cfl.size(); ++i) { // for each source input file 
    /*string date = extractDateFromPath(cfl[i]);
    if(trDates.find(date) == trDates.end()) { // skip nontrading days 
      continue; 
    }*/
    DocFeats *txf = new DocFeats(cfl[i], true, params);  // get its features
    if(txf->noSnts() == 0) {
      //cerr << "WARNING: No features from file: " << cfl[i] << endl;
      ++numZero;
    }
    else {
      //cerr << "Loaded " << txf->noSnts() << " from " << cfl[i] << endl;
    }
    vtxf.push_back(txf);
  }
  cerr << "Have " << vtxf.size() << " feature sets. (";
  cerr << numZero << " are empty.)" << endl;
  TxtReg txtreg(vtxf, params, gls, weakClfs);
}
void runTSReg(Parameters& params, GoldLabels& gls, 
  threeD_t& weakClfs, bool rearr=false) {
  TSFeats tsf(params, gls);
  tsf.readTSData();
  if(rearr) rearrange(gls, params, &tsf);
  tsf.extraFeatures();
  TSReg rr(&tsf, gls, weakClfs);
  // !!can't run both logReg and RandForest models with online updates 
  //rr.runLogReg();
  rr.runRandForest(); 
}
void saveWeakCls(threeD_t& weakClfs, GoldLabels& gls, Parameters &params) {
  FileHandler fout("./wcOutput", ios::out, false);
  for(size_t i=0; i < weakClfs.size(); ++i) { // for each source
    for(size_t j=0; j < weakClfs[i].size(); ++j) { // for each epochs
      fout << i << "\t" << j << "\t";
      for(size_t k=0; k < weakClfs[i][j].size(); ++k) { // for each pred
        fout << weakClfs[i][j][k] << "\t";
      }
      fout << (gls.lbls()[j]) << "\t";
      fout << endl;
    }
  }
  fout.close();
}
void testAcc(threeD_t& weakClfs, GoldLabels& gls, Parameters &params) {
  float mean(0);
  cerr << "Number of outputs: " << weakClfs.size() << endl;
  cerr << "source\t\%right\t\%wrong\n";
  for(size_t i=0; i < weakClfs.size(); ++i) { // for each source
    float right(0); // count the number correct
    for(size_t j=gls.noTrainEp_; j < weakClfs[i].size(); ++j) { 
      right += (weakClfs[i][j].back() == gls.lbls_.at(j));
    }
    float pctRight = right / (float)gls.noTestEp_;
    mean += pctRight;
    cerr << i << "\t" << pctRight << "\t" << (1-pctRight) << endl; 
  }
  mean /= (float)weakClfs.size();
  cerr << "Average \% correct: " << mean << endl; 
}
void trade(Parameters &params, const float prb1) {
  // execute trade with IBAPI 
  string action = prb1 > 0.5 ? "BUY" : "SELL";
  IBWrapper ib;
  ib.order(action);
}
void testIntervals(Parameters& params) {
  /* For text, testing with intervals only works with -reg set */
  GoldLabels gldLbls(params);
  int interval(gldLbls.noTestEp_), numofPts(gldLbls.totEps());
  const int itr = atoi(params.getParamValue("vb-itr").c_str());
  const string cmd = "./perf.sh " + Utils::IntToStr(interval);
  int expr(0);
  for(int i=0; i < numofPts; i+=interval) {
    gldLbls.trEpochs_.clear();
    gldLbls.tstEpochs_.clear();
    int begTest = i; 
    int endTest = i+interval-1;
    if(endTest >= numofPts) {
      begTest = numofPts - interval;
      endTest = numofPts;
    }
    cerr << "Experiment " << ++expr << ":\tTest epochs are " << begTest << " - " << endTest << endl;
    for(int j=0; j < numofPts; ++j) {
      if((j < begTest) || (j > endTest)) {
        gldLbls.trEpochs_.insert(j);
      }
      else {
        gldLbls.tstEpochs_.insert(j);
      }
    }
    assert(gldLbls.noTrainEp_ == (int)gldLbls.trEpochs_.size());
    assert(gldLbls.noTestEp_ == (int)gldLbls.tstEpochs_.size());
    /*threeD_t weakClfs;  
    baselines(params, weakClfs, gldLbls, true);
    runTxtReg(params, gldLbls, weakClfs, false);
    saveWeakCls(weakClfs, gldLbls, params); //save test probs to file 
    system(cmd.c_str());
    CIBCC ibcc(params, weakClfs, gldLbls);
    ibcc.train(itr);
    ibcc.backtest();*/
  }
}
int main(int argc, char** argv) {
  Parameters params(argc, argv, paramdefs, NumOfParams(paramdefs));
  bool interval = params.getParam("interval") == Parameters::kTrueValue;
  if(interval) {
    testIntervals(params);
    return 1;
  }
  GoldLabels gldLbls(params);
  threeD_t weakClfs;
  //baselines(params, weakClfs, gldLbls);
  runTSReg(params, gldLbls, weakClfs);
  testAcc(weakClfs, gldLbls, params);
  //saveWeakCls(weakClfs, gldLbls, params);
  //runTxtReg(params, gldLbls, weakClfs);
  cerr << endl << "Starting IBCC..." << endl;
  const int itr = atoi(params.getParam("vb-itr").c_str());
  CIBCC ibcc(params, weakClfs, gldLbls);
  ibcc.train(itr);
  ibcc.backtest();
  return 1;
}
