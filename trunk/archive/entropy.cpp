#include "file.h"
#include "vocab.h"
#include "params.h"
#include "countmin.h"
#include <tr1/unordered_map>
#include "google/dense_hash_map"

using google::dense_hash_map;

typedef std::tr1::unordered_multimap<wordID_t, count_t> MM_t; 

// paramter definitions
const ParamDefs paramdefs[] = {
  //name, default value, abbrev, type, msg
  {"input", "../ds/data/1996-09.sents.tok.gz", "i", Parameters::kStringValue, ""},
  {"bl-input", "../ds/data/reu.96-09.unigrams.gz", "bi", Parameters::kStringValue, ""},
  {"stream-len", "200", "sl", Parameters::kIntValue, ""},
  {"epsilon", ".5", "e", Parameters::kFloatValue, ""},
  {"delta", ".1", "d", Parameters::kFloatValue, ""}
};
/*
Unigram Entropy Estimate:
Input a file stream 
Get unigram vocab id
Keep track via 3 components
Output entropy after x items

Components:
A1. Basic estimator (more later)
A2. Vector for first m^(2alpha) items 
A3. DONE! Count-Min Sketch (implement separately)
(A2 and A3 not needed - always high enough entropy) 
*/

float eps = .001;
float alp = .02;
float delta = .01;
float m = 200;
Vocab vocab;
float X(int count) {
  assert(count > 0);
  float r = static_cast<float>(count);
  float X = (r/m) * log2(m/r);
  X -= r > 1 ? ((r-1)/m) * log2(m/(r-1)) : 0;
  X *= m;
  return X;
}
float baselineH(Parameters& params) {
  //read in file 
  FileHandler fblin(params.getParamValue("bl-input"), std::ios::in);
  string line;
  int freq(0), totfreq(0);
  int max(0);
  while(getline(fblin, line)) {
    //get sum of unigram counts
    std::istringstream um(line.c_str());
    um >> line;
    um >> freq;
    if(freq > max) max = freq;
    totfreq += freq;
  }
  cerr << "total freq = " << totfreq << endl;
  cerr << "max = " << max << endl;

  fblin.reset();
  float H(0);
  float tf = (float)totfreq;
  while(getline(fblin, line)) {
    std::istringstream um(line.c_str());
    um >> line;
    um >> freq;
    float m_i = (float)freq;
    H += (m_i / tf) * log2(tf / m_i);
  }
  fblin.close();
  return H; 
}
float getH(MM_t* estmap, int s1, int s2){
  //assert(estmap->size() == s1 * s2);
  cerr << "F_0 = " << vocab.size() << endl;
  cerr << "# of estimators = (actual) " << estmap->size() << 
    " vs (defined) " << s1 * s2 << endl;
  float ests[s1*s2];
  //iterate through basic estimators
  int idx(0);
  piterate(estmap, itr) { // should be unordered!
    ests[idx++] = X(itr->second); //store X() for each sample in the map
    //cerr << ests[idx -1] << endl;
  }
  //get mean of each batch of s1 samples -- should give s2 numbers
  float means[s2];
  for(int j = 0; j < s2; ++j) {
    float sum(0);
    for(idx = j * s1; idx < (j * s1) + s1; ++idx)
      sum += ests[idx]; 
    means[j] = sum / s1;
  }
  for(int i = 0; i < s2; i++)
    cerr << means[i] << endl;
  //get median of the s2 numbers
  int mid = (int)floor(s2/2);
  std::nth_element(means, means + mid, means + s2);
  float H = means[mid];
  return H;
}
int main(int argc, char** argv) {
  srand(time(NULL));
  Parameters params(argc, argv, paramdefs, NumOfParams(paramdefs));
  // initialize stuff
  m = atof(params.getParamValue("stream-len").c_str());
  eps = atof(params.getParamValue("epsilon").c_str());
  delta = atof(params.getParamValue("delta").c_str());
  float s1  = ceil((8/pow(eps, 2)) * pow(log2(m),2));  // force integers 
  float s2 = ceil(4 * log2(1/delta)); 
  int s1s2 = (int)(s1 * s2); 
  cerr << "s1 = " << s1 << ", s2 = " << s2 << ", s1s2 = " << s1s2<< endl;
  MM_t basicEsts;
  MM_t::iterator beItr;
  std::pair<MM_t::iterator, MM_t::iterator> ret;
  int tokTot(0), basEstTot(0);
  // get stream
  FileHandler fin(params.getParamValue("input"), std::ios::in);
  string line, token;
  wordID_t id(Vocab::kOOVWordID);
  // maintain samples
  while(getline(fin, line)) {
    std::istringstream stream(line);
    while((stream >> token) && ++tokTot <= m) {
      id = vocab.getWordID(token);
      // sample first s1s2 items (should sample uar here!)
      if(basEstTot < s1s2) { 
        basicEsts.insert(MM_t::value_type(id,0));
        ++basEstTot;
        /*float chance = Utils::rand<unsigned>(RAND_MAX); 
        chance /= (float)RAND_MAX;
        if(chance > .5) {
          basicEsts.insert(MM_t::value_type(id,0));
          ++basEstTot;
        }*/
      }
      // increment all estimators that track this id
      ret = basicEsts.equal_range(id);
      for(beItr = ret.first; beItr != ret.second; ++beItr)
        ++beItr->second;
    }
    if(tokTot >= m) break; // stop processing stream 
    if(tokTot % 100000 == 0) cerr << "processed " << tokTot << endl;
  }
  fin.close();
  m = tokTot ; // length of stream is total tokens seen 
  //output H
  cerr << "total tokens processed= " << tokTot << endl;
  cerr << "entropy= " << getH(&basicEsts, (int)s1, (int)s2) << endl;
  cerr << "baseline H = " << baselineH(params) << endl;
  return 1;
}
