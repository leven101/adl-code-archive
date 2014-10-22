#include "model1.h"
#ifndef USEMODEL1FORHMM 
#include "params.h"

const ParamDefs paramdefs[] = {
  //name, default value, abbrev, type, msg
  {"source", "./source", "src", Parameters::kStringValue, ""},
  {"target", "./target", "trg", Parameters::kStringValue, ""},
  {"EMIterations", "1", "em", Parameters::kIntValue, ""},
};
int main(int argc, char**argv) {
  Parameters params(argc, argv, paramdefs, NumOfParams(paramdefs));
  Model1 model1(params.getParam("source"), params.getParam("target"));
  model1.EM(atoi(params.getParam("EMIterations").c_str()));
  model1.printTTable();
  return 1;
}
#endif

int Model1::EM(const int noIters) {
  initializeUniformly();
  std::vector<wordID_t> *src, *trg;
  std::map<wordID_t, float> trgSntTot, totalSrc; 
  int noSents = src_sents.size();
  unsigned int srcVcbSize = srcVcb->size(); 
  unsigned int trgVcbSize = trgVcb->size();
  cerr << "srcVcbSize = " << srcVcbSize << "\ttrgVcbSize = " << trgVcbSize << endl;
  ttable_t::iterator pairPrb;
  wordID_t tID, sID;
  for(int iter = 1; iter <= noIters; ++iter) {
    initZeroCounts();
    totalSrc.clear();
    for(int snt = 0; snt < noSents; ++snt) { // for all sentences
      cerr << "NEW sentence number: " << (snt + 1) << endl;
      trg = &trg_sents[snt];
      src = &src_sents[snt];
      cerr << "target size = " << trg->size() << "\tsrc size = " << src->size() << endl; 
      
      for(size_t i = 0; i < trg->size(); ++i) { // for all trg words in sentence 
        tID = trg->at(i);
        trgSntTot[tID] = 0; // resets for repeat words
        for(size_t j = 0; j < src->size(); ++j) { // for all src words in sentence
          pairPrb = wrdProbs_->find(wrdPair_t(src->at(j), tID));
          assert(pairPrb != wrdProbs_->end());
          trgSntTot[tID] += pairPrb->second; //get total sum_j Pr(tID_i | sID_j) 
        } // src words
      } // trg words
      iterate(trgSntTot, itr)
        cerr << trgVcb->getWord(itr->first) << " " << itr->second << endl; 
      
      for(size_t i = 0; i < trg->size(); ++i) { // for all trg words 
        tID = trg->at(i);
        for(size_t j = 0; j < src->size(); ++j) { // for all src words
          sID = src->at(j);
          pairPrb = wrdProbs_->find(wrdPair_t(sID, tID));
          assert(pairPrb != wrdProbs_->end());
          float x = pairPrb->second / trgSntTot[tID]; // sentence-level normalized probability
          (*counts_)[pairPrb->first] += x;  // collect counts evidence per sentence
          totalSrc[sID] += x;  // store for normalization   
        } // end src words
      } // end trg words
      cerr << "______" << endl;
    } // end all sentences
    piterate(counts_, itr)
      cerr << "count of (" << srcVcb->getWord(itr->first.first) << "," << 
        trgVcb->getWord(itr->first.second) << ") " << itr->second << endl;
    iterate(totalSrc, itr)
      cerr << "totalSrc[" << srcVcb->getWord(itr->first) << "]=" << itr->second << endl;
    
    for(sID = 0; sID <= srcVcbSize; ++sID) // include NULL word
      for(tID = 1; tID <= trgVcbSize; ++tID) {
        pairPrb = wrdProbs_->find(wrdPair_t(sID, tID));
        if(pairPrb != wrdProbs_->end())
          pairPrb->second = (*counts_)[pairPrb->first] / totalSrc[sID];
      }
  } // EM Iterations
  return 1;
}
void Model1::initializeUniformly() {
  // for each sentence pair set probability each word goes to 
  // each word in the sentence + NULL word
  std::vector<wordID_t> *src, *trg;
  int noSents = src_sents.size();
  assert(noSents = trg_sents.size());
  for(int i = 0; i < noSents; ++i) {
    src = &src_sents[i];
    trg = &trg_sents[i];
    float uniform = 1.0 / (float)src->size();
    piterate(src, sID)  // for each src word
      piterate(trg, tID)  // add link to each coocurring trg word
        wrdProbs_->insert(pairStat_t(wrdPair_t(*sID, *tID), uniform));
  }
  //printTTable();
}
void Model1::initZeroCounts() {
  assert(wrdProbs_->size() > 0);
  counts_->clear();
  piterate(wrdProbs_, itr) 
    counts_->insert(pairStat_t(itr->first, 0));
}
void Model1::loadCorpora(string srcCorpus, string trgCorpus) {
// for prototype loads all parallel sentences into memory
  totSrcToks = 0; totTrgToks = 0;
  string line;
  std::vector<string> v;
  std::vector<wordID_t> vi;
  FileHandler fsrc(srcCorpus, std::ios::in);
  while(getline(fsrc, line)) {
    vi.clear();
#ifndef USEMODEL1FORHMM 
    vi.push_back(0); // add NULL word to each source sentence
#endif
    Utils::trim(line);
    Utils::splitToStr(line, v, " ");
    iterate(v, itr)
      vi.push_back(srcVcb->getWordID(*itr));
    src_sents.push_back(vi);
    totSrcToks += vi.size();
  }
  FileHandler ftrg(trgCorpus, std::ios::in);
  while(getline(ftrg, line)) {
    vi.clear();
    Utils::trim(line);
    Utils::splitToStr(line, v, " ");
    iterate(v, itr)
      vi.push_back(trgVcb->getWordID(*itr));
    trg_sents.push_back(vi);
    totTrgToks += vi.size();
  }
  fsrc.close();
  ftrg.close();
  srcVcb->makeClosed();
  trgVcb->makeClosed();
}
void Model1::printTTable(string fname) {
  if(fname.empty()) {
    piterate(wrdProbs_, pitr) {
      cerr << srcVcb->getWord(pitr->first.first) << " ";
      //cerr << "(" << pitr->first.first << ") "; 
      cerr << trgVcb->getWord(pitr->first.second) << " ";
      //cerr << "(" << pitr->first.second << ") ";
      cerr << pitr->second << endl;
    }
  }
}
