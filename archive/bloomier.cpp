#include <map>
#include <stack>
#include <list>
#include <set>
#include <vector>
#include <iomanip>
#include "types.h"
#include "hash.h"
#include "vocab.h"
#include "file.h"
#include "utils.h"
#include "params.h"
#include "RandLMFilter.h"
#include "quantizer.h"
#include "bloomier.h"
#include "dirList.h"

using randlm::Filter;

// paramter definitions
const ParamDefs paramdefs[] = {
  //name, default value, abbrev, type, msg
  {"memory", "2", "m", Parameters::kIntValue, "Total megabytes for filter"},
  {"width", "16", "w", Parameters::kIntValue, "Width in bits for each ngram"},
  {"order", "5", "o", Parameters::kIntValue, "N-gram order"},
  {"graph", Parameters::kFalseValue, "", Parameters::kBoolValue, "Signal graphing"}, 
  {"query", Parameters::kFalseValue, "", Parameters::kBoolValue, "Signal querying"},
  {"train", Parameters::kFalseValue, "", Parameters::kBoolValue, "Signal training"},
  {"input", "_NO_FILE_", "", Parameters::kStringValue, ""},
  {"get-counts", Parameters::kFalseValue, "gc", Parameters::kBoolValue, ""},
  {"data-dir", "_NO_DIR_", "dir", Parameters::kStringValue, "Directory with training data."},
  {"qbase", "8", "", Parameters::kFloatValue, "Used as 2^(1/qBase)"},
  {"target", "200000000", "xx", Parameters::kIntValue, ""}
};
Vocab* vocab_;
UnivHash_linear<count_t>* h_;
count_t cells_ = 0;
const int k_ = 3;
int order_(0);
int items_(0);
int tot_subsets_(0), subset_size_(0), cur_subset_(0);
int** left2right(0);
static int last_(0);  // state pointer to last ngram processed 
std::map<int, std::set<int> >* right2left = new std::map<int, std::set<int> >();
typedef std::map<int, std::set<int> >::iterator r2l_itr_t;
wordID_t** lastGrams(0);
string sdir_("");
void deconstruct2() {
  if(vocab_) delete vocab_;
  for(int i = 0; i < subset_size_; ++i) 
    delete[] left2right[i];
  if(left2right) delete[] left2right;
  for(int i = 0; i < tot_subsets_; ++i)
    delete[] lastGrams[i];
  if(lastGrams) delete[] lastGrams;
  delete right2left;
  if(h_) delete h_;
}
bool setupGraph2(Parameters* params, const int lines2read) {
  cerr << "lines2read = " << lines2read << endl; 
  // instantiate new hash function for each graph 
  h_ = new UnivHash_linear<count_t>(cells_, k_, PRIME);
  FileHandler fin(params->getParamValue("input"), std::ios::in);
  string line;
  int i(0); 
  std::vector<string> vwords, vtmp;
  for(int j = last_; j > 0; --j)
    getline(fin, line);  // get to last place processed in file
  while(i < lines2read) {
    assert(getline(fin, line));
    // get word ID's
    Utils::splitToStr(line, vtmp);
    Utils::splitToStr(vtmp[0], vwords, " ");
    wordID_t ids[vwords.size()];
    for(int j = 0; j < (int)vwords.size(); ++j) 
      ids[j] = vocab_->getWordID(vwords[j]);
    for(int k = 0; k < k_; ++k) { // through each hash function k
      int arow = h_->hash(ids, vwords.size(), k); // current row mapped to
      // keep track of which rows are associated with this ngram
      left2right[i][k] = arow;
      // track which left nodes are matched to this right one
      (*right2left)[arow].insert(i);
    }
    if(i % 1000000 == 0 && (i > 0)) 
      cerr << "loaded " << i << " ngrams\n";
    ++i; // increment index pointer
    if(i == lines2read) {  // track last ngram of each subset
      for(int j = 0; j < (int)vwords.size(); ++j) 
        lastGrams[cur_subset_ -1][j] = ids[j];
    }
  }
  bool done = fin.eof() || (last_ >= items_);
  fin.close(); // close file between subsets 
  return done;
}
bool verifyGraph2(int size2match) {
  //push all degree one nodes to list
  std::list<int> degree_one;
  std::stack<std::pair<int, int> > matched;
  r2l_itr_t itr1_r2l(0), itr2_r2l(0);
  piterate(right2left, it) {
    if(it->second.size() == 1) 
      degree_one.push_back(it->first);
  }
  //cerr << "degree_one size = " << degree_one.size() << endl;
  while(!degree_one.empty()) {
    int rhs = degree_one.front(); 
    degree_one.pop_front(); // remove from list
    itr1_r2l = right2left->find(rhs);
    //cerr << "rhs = " << rhs << endl;
    //assert((*right2left)[rhs].size() <= 1);
    if(itr1_r2l->second.size() == 1) {
      int lhs = *itr1_r2l->second.begin();
      // push onto stack "matched"
      matched.push(std::pair<int, int>(lhs, rhs));
      if(matched.size() % 1000000 == 0) 
        cerr << "matched " << matched.size() << endl;
      //cerr << "matched lhs " << lhs << " with rhs " << rhs << endl;
      // use left2right[] to remove this ngram index from matching right2left[] nodes
      for(int k = 0; k < k_; ++k) {
        int rnode = left2right[lhs][k];
        itr2_r2l = right2left->find(rnode);
        itr2_r2l->second.erase(lhs);  // (*right2left)[rnode].erase(lhs);
        if(itr2_r2l->second.size() == 1) {// if this node is of degree one add to list 
          //cerr << "pushing " << rnode << " onto degree_one list\n";
          degree_one.push_back(rnode);
        }
      }
    }
  }
  //cerr << "matched size = " << matched.size() << endl;
  bool success = (int)matched.size() == size2match ? true : false;
  //OLD::bool success = (int)matched.size() == items_ ? true : false;
  cerr << "subset " << cur_subset_ << (success ? " SUCCESS" : " FAILURE") << endl;
  if(success) { // save order and hash params
    string fprefix = sdir_; 
    if(cur_subset_ < 10) fprefix += "0";
    fprefix += Utils::IntToStr(cur_subset_);
    FileHandler fData(fprefix + ".matched.gz", std::ios::out);
    while(!matched.empty()) {
      std::pair<int, int> p = matched.top();
      //cerr << p.first << " goes to " << p.second << endl;
      // write out "filter index \t ngram \t count"
      fData << p.second << "\t" << p.first << endl;
      matched.pop();
    }
    fData.close();
    FileHandler fInfo(fprefix + ".hash.gz", std::ios::out);
    h_->save(&fInfo);
    fInfo.close();
  }
  delete h_;
  h_ = 0;
  return success;
}
void outputOrder(Parameters* params, const int lines2read) {
  // read in ngrams
  std::vector<string>* ngrams = new std::vector<string>();
  FileHandler fin(params->getParamValue("input"), std::ios::in);
  string line;
  int i(0); 
  std::vector<string> vwords, vtmp;
  for(int j = last_; j > 0; --j)
    getline(fin, line);  // get to last place processed in file
  while(i++ < lines2read) {
    assert(getline(fin, line));
    ngrams->push_back(line);
  }
  fin.close();
  // loop through matched 
  string fprefix = sdir_;
  if(cur_subset_ < 10) fprefix += "0";
  fprefix += Utils::IntToStr(cur_subset_);
  FileHandler fData(fprefix + ".matched.gz", std::ios::in);
  FileHandler fout(fprefix + ".ordered.gz", std::ios::out);
  std::vector<int> vec;
  // print ngrams in order of matched
  while(getline(fData, line)) {
    assert(Utils::splitToInt(line, vec, "\t") == 2);
    // vec[0] is array index, vec[1] is file index 
    fout << vec[0] << "\t" << (*ngrams)[vec[1]] << "\n";
  }
  fData.close();
  fout.close();
  delete ngrams;
}
void graph(Parameters* params) {
  vocab_ = new Vocab();
  FileHandler fin(params->getParamValue("input"), std::ios::in);
  string line;
  order_ = atoi(params->getParamValue("order").c_str());
  int target = atoi(params->getParamValue("target").c_str());
  sdir_ = params->getParamValue("data-dir");
  while(getline(fin, line) && (items_ < target)) 
    ++items_;
  cerr << "items = " << items_ << endl;
  fin.close();
  tot_subsets_ = 10;
  subset_size_ = static_cast<int>(ceil(static_cast<float>(items_) / 
    static_cast<float>(tot_subsets_)));
  cerr << "subset_size = " << subset_size_ << endl;
  // handle last subset size of data in case different than other subsets
  int lastsize = items_ % subset_size_ == 0 ? subset_size_ : 
    items_ % subset_size_;
  cerr << "lastsize = " << lastsize << endl;
  left2right = new int*[subset_size_]; 
  for(int i = 0; i < subset_size_; ++i)
    left2right[i] = new int[k_];  // initialized below
  lastGrams = new wordID_t*[tot_subsets_]; 
  for(int i = 0; i < tot_subsets_; ++i) {
    lastGrams[i] = new wordID_t[order_];
    for(int j = 0; j < order_; ++j)
      lastGrams[i][j] = 0;  // initilize with OOV ID
  }
  cells_ = (count_t)ceil(subset_size_* 1.23);
  for(cur_subset_ = 1; cur_subset_ <= tot_subsets_; ++cur_subset_) {
    cerr << "processing subset = " << cur_subset_ << endl;
    int num2process = cur_subset_ == tot_subsets_ ? lastsize : subset_size_;
    do {  // loop until matched graph found 
      for(int i = 0; i < subset_size_; ++i)
        for(int k = 0; k < k_; ++k) left2right[i][k] = -1;
      right2left->clear(); // clear data
      setupGraph2(params, num2process);
    } while(!verifyGraph2(num2process));
    right2left->clear();
    outputOrder(params, num2process);
    last_ += num2process; 
  }
  assert(last_ == items_);
  // save vocab 
  FileHandler fsubinfo(sdir_ + "subinfo.gz", std::ios::out);
  vocab_->save(&fsubinfo);
  fsubinfo.write((char*)&cells_, sizeof(cells_));
  fsubinfo.write((char*)&tot_subsets_, sizeof(tot_subsets_));
  fsubinfo.write((char*)&order_, sizeof(order_));
  for(int i = 0; i < tot_subsets_; ++i)
    for(int j = 0; j < order_; ++j)
      fsubinfo.write((char*)&lastGrams[i][j], sizeof(lastGrams[i][j]));
  fsubinfo.close();
}
void train2(Parameters* params) {
  float base = atof(params->getParamValue("qbase").c_str());
  LogQtizer qtizer(base);
  string sdir = params->getParamValue("data-dir");
  // get all subset information 
  FileHandler fsubinfo(sdir + "subinfo.gz", std::ios::in);
  vocab_ = new Vocab(&fsubinfo);  // read in vocab
  fsubinfo.read((char*)&cells_, sizeof(cells_));
  fsubinfo.read((char*)&tot_subsets_, sizeof(tot_subsets_));
  cerr << "tot_subsets_ = " << tot_subsets_ << endl;
  fsubinfo.read((char*)&order_, sizeof(order_));
  lastGrams = new wordID_t*[tot_subsets_];   // read in lastGrams
  for(int i = 0; i < tot_subsets_; ++i)
    lastGrams[i] = new wordID_t[order_];
  for(int i = 0; i < tot_subsets_; ++i)
    for(int j = 0; j < order_; ++j)
      fsubinfo.read((char*)&lastGrams[i][j], sizeof(lastGrams[i][j]));
  fsubinfo.close();
  cerr << "cells = " << cells_ << endl;
  // setup structures and LM tracking 
  int maxCode(0), corpusSize(0), sofar(0);
  int cellWidth = atoi(params->getParamValue("width").c_str());
  Filter<count_t>* filter = new Filter<count_t>(cells_, cellWidth);
  cerr << "filter memusage = " << filter->size() << endl;
  UnivHash_linear<count_t> fghash(pow(2, cellWidth), 10, PRIME); // fingerprint hash function
  string cur_file(""), line;
  std::vector<string> vline, vwords;
  DirList dir(sdir, "*.ordered.gz");  // ordered file extensions
  assert((int)dir.cFileList_.size() == tot_subsets_);
  iterate(dir.cFileList_, m_itr) {         // open each file
    // load this subset's hash params
    cur_file = *m_itr;
    Utils::rtrim(cur_file, ".ordered.gz");
    FileHandler fhash(cur_file + ".hash.gz", std::ios::in); // read in hash params
    h_ = new UnivHash_linear<count_t>(&fhash);
    fhash.close();
    // fill this subset's filter
    FileHandler fsubdata(*m_itr, std::ios::in);
    while(getline(fsubdata, line)) {
      // split apart index, ngram, value
      assert(Utils::splitToStr(line, vline) == 3);
      int idx = atoi(vline[0].c_str());
      int v = qtizer.code(atoi(vline[2].c_str())); // quantize value
      if(v > maxCode) maxCode = v;  // track max value
      Utils::splitToStr(vline[1], vwords, " ");

      // keep track of total items from training data minus "<s>"
      if(vwords.size() == 1)  // sum of unigram counts is total corpus size
        corpusSize += (vwords[0] != Vocab::kBOS ? atoi(vline[2].c_str()) : 0);

      wordID_t IDs[vwords.size()]; // get each vocab ID
      for(int i = 0; i < (int)vwords.size(); ++i)
        IDs[i] = vocab_->getWordID(vwords[i]);
      count_t f = BloomierXX::nonZeroSignature(IDs, vwords.size(), &fghash); // get fingerprint
      assert(log2(f) < cellWidth);
      count_t b(cells_ + 1);
      bool fnd(false);
      for(int k = 0; k < k_; ++k) {  // for each k hash function
        int h = h_->hash(IDs, vwords.size(), k);
        if(h != idx) {
          count_t arrVal = filter->read(h);
          b = (b == cells_ + 1 ? arrVal : (b ^ arrVal)); // xor dependent nodes values 
        }
        else fnd = true;
      }
      assert(fnd);
      filter->write(idx, (v ^ f ^ b));
      if(++sofar % 5000000 == 0) cerr << "processed " << sofar << endl;
    }
    fsubdata.close();
    //save bloomier filter for each subset file
    FileHandler ffilter(cur_file + ".bfilter.gz" , std::ios::out|std::ios::binary);
    filter->save(&ffilter);
    h_->save(&ffilter);
    ffilter.close();

    delete h_;
    h_ = 0;
    filter->reset();  // clear filter and reuse for next subset
  }
  // save other info
  FileHandler finfo(sdir + "filter.info.gz", std::ios::out);
  vocab_->save(&finfo);
  fghash.save(&finfo);
  finfo.write((char*)&cells_, sizeof(cells_));
  finfo.write((char*)&cellWidth, sizeof(cellWidth));
  finfo.write((char*)&maxCode, sizeof(maxCode));
  cerr << "maxcode = " << maxCode << endl;
  finfo.write((char*)&corpusSize, sizeof(corpusSize));
  cerr << "corpusSize = " << corpusSize << endl;
  finfo.write((char*)&tot_subsets_, sizeof(tot_subsets_));
  finfo.write((char*)&order_, sizeof(order_));
  qtizer.save(&finfo);
  for(int i = 0; i < tot_subsets_; ++i)
    for(int j = 0; j < order_; ++j)
      finfo.write((char*)&lastGrams[i][j], sizeof(lastGrams[i][j]));
  finfo.close();
  delete filter;
  filter = 0;
}
void test(Parameters* params) {
  BloomierXX* b;
  b = new BloomierXX(params->getParamValue("data-dir"), 5);
  FileHandler fin2(params->getParamValue("input"), std::ios::in);
  string line;
  std::vector<string> vwords;
  while(getline(fin2, line)) {
    Utils::splitToStr(line, vwords, " ");
    wordID_t ngram[vwords.size()];
    for(int i = 0; i < (int)vwords.size(); ++i)
      ngram[i] = b->vocab_->getWordID(vwords[i]);
    if(params->getParamValue("get-counts") == Parameters::kTrueValue)
      cout << b->query(ngram, vwords.size()) << endl;
    else
      cout << b->getProb(ngram, vwords.size(), NULL) << endl;
  }
  delete b;
}
int main(int argc, char** argv) {
  Parameters params(argc, argv, paramdefs, NumOfParams(paramdefs));
  time_t start, finish;
  time(&start);
  if(params.getParamValue("graph") == Parameters::kTrueValue) {
    srand(time(NULL));
    graph(&params);
  }
  else if(params.getParamValue("train") == Parameters::kTrueValue) {
    train2(&params);
  }
  else if(params.getParamValue("query") == Parameters::kTrueValue) {
    test(&params); 
  }
  time(&finish);
  cerr << "time taken was " << difftime(finish, start) << " seconds\n";
  deconstruct2();
  return 1;
}
/*
 * //MY ATTEMPT TO DO IT THE PERFECT HASHING WITHOUT SPLITTING THE DATA INTO SUBSETS
void deconstruct() {
  if(vocab_) delete vocab_;
  for(int i = 0; i < items_; ++i)
    delete[] left2right[i];
  if(left2right) delete[] left2right;
  delete right2left;
  if(h_) delete h_;
}
void preprocess(Parameters* params) {  
  FileHandler fin(params->getParamValue("input"), std::ios::in);
  string line;
  int target = atoi(params->getParamValue("target").c_str());
  while(getline(fin, line) && (items_ < target)) 
    ++items_; 
  fin.close();
  left2right = new int*[items_]; 
  for(int i = 0; i < items_; ++i)
    left2right[i] = new int[3];
  cells_ = (count_t)ceil(items_ * 1.23);
  cerr << cells_ << " rows for " << items_ << " ngrams\n";
}
void setupGraph(Parameters* params) {
  vocab_ = new Vocab();
  h_ = new UnivHash_linear<count_t>(cells_, k_, PRIME);
  FileHandler fin(params->getParamValue("input"), std::ios::in);
  string line;
  int i(0); 
  std::vector<string> vwords, vtmp;
  // this loop sets up the graph
  while(getline(fin, line) && (i < items_)) { 
    left2right[i][0] = 0, left2right[i][1] = 0, left2right[i][2] = 0; 
    // get word ID's
    Utils::splitToStr(line, vtmp);
    Utils::splitToStr(vtmp[0], vwords, " ");
    wordID_t ids[vwords.size()];
    for(int j = 0; j < (int)vwords.size(); ++j) 
      ids[j] = vocab_->getWordID(vwords[j]);
    for(int k = 0; k < k_; ++k) { // through each hash function k
      int arow = h_->hash(ids, vwords.size(), k); // current row mapped to
      // keep track of which rows are associated with this ngram
      left2right[i][k] = arow;
      //cerr << "left2right[" << i << "][" << k << "] = " << arow << endl;
      // track which left nodes are matched to this right one
      (*right2left)[arow].insert(i);
    }
    if(++i % 1000000 == 0) cerr << "processed " << i << " ngrams\n";
  }
  fin.close();
}
void train(Parameters* params) {
  LogQtizer qtizer(8);
  FileHandler fInfo("../ds/data/bl.info.gz", std::ios::in);
  vocab_ = new Vocab(&fInfo);  // read in vocab
  h_ = new UnivHash_linear<count_t>(&fInfo);  // read in hash params
  fInfo.close();
  int maxCode(0), corpusSize(0);
  int cellWidth = atoi(params->getParamValue("width").c_str());
  Filter<count_t>* filter = new Filter<count_t>(cells_, cellWidth); // main array
  UnivHash_linear<count_t> fghash(cellWidth, 10, PRIME); // fingerprint hash function
  // get file with ngrams and indexes in it
  FileHandler fData("../ds/data/bl.data.gz", std::ios::in); 
  string line;
  std::vector<string> vline, vwords;
  int sofar(0);
  while(getline(fData, line)) { // for each line
    if(++sofar % 100000 == 0) cerr << sofar << endl;
    // split apart index, ngram, value
    assert(Utils::splitToStr(line, vline) == 3);
    int idx = atoi(vline[0].c_str());
    int v = qtizer.code(atoi(vline[2].c_str())); // quantize value
    assert(log2(v) < cellWidth);
    if(v > maxCode) maxCode = v;  // track max value
    Utils::splitToStr(vline[1], vwords, " ");

    // keep track of total items from training data minus "<s>"
    if(vwords.size() == 1)  // sum of unigram counts is total corpus size
      corpusSize += (vwords[0] != Vocab::kBOS ? atoi(vline[2].c_str()) : 0);

    wordID_t IDs[vwords.size()]; // get each vocab ID
    for(int i = 0; i < (int)vwords.size(); ++i)
      IDs[i] = vocab_->getWordID(vwords[i]);
    count_t f = BloomierXX::nonZeroSignature(IDs, vwords.size(), &fghash); // get fingerprint
    assert(log2(f) < cellWidth);
    count_t b(cells_ + 1);
    bool fnd(false);
    for(int k = 0; k < k_; ++k) {  // for each k hash function
      int h = h_->hash(IDs, vwords.size(), k);
      if(h != idx) {
        count_t arrVal = filter->read(h);
        b = (b == cells_ + 1 ? arrVal : (b ^ arrVal)); 
      }
      else fnd = true;
    }
    assert(fnd);
    filter->write(idx, (v ^ f ^ b));
  }
  fData.close();
  // save bloomier filter
  FileHandler fsave("../ds/data/bl.lm.gz", std::ios::out|std::ios::binary);
  vocab_->save(&fsave);
  filter->save(&fsave);
  h_->save(&fsave);
  fghash.save(&fsave);
  fsave.write((char*)&maxCode, sizeof(maxCode));
  fsave.write((char*)&cells_, sizeof(cells_));
  fsave.write((char*)&cellWidth, sizeof(cellWidth));
  fsave.write((char*)&corpusSize, sizeof(corpusSize));
  delete filter;
}
 */
