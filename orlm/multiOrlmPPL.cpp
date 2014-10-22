#include "onlineRLM.h"
#include "file.h"
#include "params.h"
#include "utils.h"
#include "TextStats.h"

// paramter definitions
const ParamDefs paramdefs[] = {
  //name, default value, abbrev, type, msg
  {"orlm", "_NO_ORLM_", "", Parameters::kStringValue, "Path to read in trained ORLM from file"},
  {"pplFile", "_NO_PPLFILE_", "", Parameters::kStringValue, "Path to read in ppl test file"},
};


int main(int argc, char** argv) {
  Parameters params(argc, argv, paramdefs, NumOfParams(paramdefs));
  TextStats docStats;
  const count_t order = 5;
  // load ORLM
  string fname = params.getParamValue("orlm"); 
  FileHandler fin(fname, ios::in|ios::binary);
  OnlineRLM<count_t>* ph = new OnlineRLM<count_t>(&fin, order);
  fin.close();
  // load pplFile
  string pplFile = params.getParamValue("pplFile"); 
  FileHandler fppl(pplFile, ios::in);
  vector<string> vwords;
  string line;
  // for each sentence in pplFile
  while(getline(fppl, line)) {
    TextStats sentStats;
    Utils::trim(line);
    if(line.empty()) continue;
    // print sentence
    cout << line << endl;
    // split, get vocab IDs (add BOS and EOS IDs)
    Utils::splitToStr(line, vwords, " "); 
    count_t sentLength = vwords.size() + 2; // 1-indexed
    wordID_t wids[sentLength];
    wids[0] = ph->vocab_->getWordID(Vocab::kBOS);
    for(count_t i=1; i < sentLength-1; ++i) {
      wids[i] = ph->vocab_->getWordID(vwords[i-1]);
    }
    wids[vwords.size()+1] = ph->vocab_->getWordID(Vocab::kEOS);
    // get probability for each n-gram in sentence 
    int found(0), len(1), startIndex(0);
    for(count_t i=0; i < sentLength-1; ++i) {
      // find prob for current point in sentence to either 0-index of sentence or order
      if(i > (order-2)) ++startIndex;  
      if(len < order) ++len;
      //cout << "i = " << i << "\tstartIndex = " << startIndex << "\tlen = " << len << endl;
      float prob = ph->getProb(&wids[startIndex], len, NULL, &found);
      // print it as SRILM does
      cout << "\tp( " << ph->vocab_->getWord(wids[i+1]) << " | " 
        //<< (wids[i] != Vocab::kOOVWordID ? ph->vocab_->getWord(wids[i]) : "")
        << ph->vocab_->getWord(wids[i])
        << (i != 0 ? " ..." : " ") << ")\t= ";
      if(found) {
        cout << "[" << found << "gram]";
      } else {
        cout << "[1gram]";
        //cout << "[OOV]";
      }
      cout << " " << LogPtoProb(prob) << " [ " << prob << " ]" << endl;
      sentStats.prob += prob;
    }
      // print sentence stats
      sentStats.numSentences = 1;
      sentStats.numWords = sentLength - 2;
      cout << sentStats << endl;
      docStats.increment(sentStats);
  }
  cout << "file " << pplFile << ": " << docStats;
  delete ph;
  return 0;
}

