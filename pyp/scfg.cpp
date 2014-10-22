// g++ scfg.cpp ../file.o ../params.o -I../include
#include <algorithm>
#include <queue>
#include "file.h"
#include "params.h"
#include "utils.h"
#include "types.h"

using std::vector;
using std::queue;
const ParamDefs paramdefs[] = {
  //name, default value, abbrev, type, msg
  {"sout", "source.txt", "", Parameters::kStringValue, ""},
  {"tout", "target.txt", "", Parameters::kStringValue, ""},
  {"grammarPath", "./miniscfg.txt", "g", Parameters::kStringValue, ""},
  {"sentences", "1", "n", Parameters::kIntValue, ""},
};

struct scfgRule {
  vector<string> s,t;
  float prob;
};
vector<scfgRule> rules;
vector<float> distr;
Parameters* params;
FileHandler* sout, *tout;
const string delm = "|||";
bool isTerm(string i) {
  return atoi(i.c_str()) >= 0;
}
void readGrammarFile(string path) {
  // read in SCFG file 
  FileHandler fin(path, std::ios::in);
  string line;
  vector<string> vline, vrule;
  while(getline(fin, line)) {
    scfgRule rule;
    // split the rule into a vector
    Utils::splitToStrMD(line, vline, delm);
    // use 1=source, 2=target, 3=prob
    Utils::splitToStr(vline[1], vrule, " "); 
    int sarity(0);
    iterate(vrule, v) {
      rule.s.push_back(*v);
      if(!isTerm(rule.s.back())) ++sarity;
    }
    Utils::splitToStr(vline[2], vrule, " "); 
    int tarity(0);
    iterate(vrule, v) {
      rule.t.push_back(*v);
      if(!isTerm(rule.t.back())) ++tarity;
    }
    if(sarity != tarity) {
      cout << "ERROR!! Bad grammar rule: " << line << endl;
      assert(false);
    }
    rule.prob = atof(vline[3].c_str());
    rules.push_back(rule);
  }
  // build rule distribution 
  distr.resize(rules.size());
  distr[0] = rules[0].prob;
  for(int i=1; i < rules.size(); ++i) {
    distr[i] = rules[i].prob + distr[i-1]; 
  }
  assert(distr.back() <= 1);
  if(distr.back() != 1)
    cerr << "WARNING!! grammar ratios don't sum to 1.\n";
}
void generateSourceSide(queue<int>& cache) {
  // choose a rule randomly
  vector<float>::iterator it;
  float rs = drand48();
  it = std::lower_bound(distr.begin(), distr.end(), rs);
  int index = int(it - distr.begin());
  cache.push(index);
  scfgRule& rule = rules[index];
  // print source side
  iterate(rule.s, sit) {
    if(isTerm(*sit))
      *sout << *sit << " ";
    else
      generateSourceSide(cache);
  }
}
void generateTargetSide(queue<int>& cache) {
  assert(!cache.empty());
  scfgRule& rule = rules[cache.front()];
  iterate(rule.t, tit) {
    if(isTerm(*tit))
      *tout << *tit << " ";
    else {
      cache.pop();
      generateTargetSide(cache);
    }
  }
}
void generateCorpus() {
  int numSents = atoi(params->getParam("sentences").c_str());
  for(int i=0; i < numSents; ++i) {
    queue<int> cache;
    generateSourceSide(cache);
    *sout << endl;
    generateTargetSide(cache);
    *tout << endl;
  }
}
int main(int argc, char** argv) {
  params = new Parameters(argc, argv, paramdefs, NumOfParams(paramdefs));
  srand48(time(NULL));
  readGrammarFile(params->getParam("grammarPath"));
  sout = new FileHandler(params->getParam("sout"), std::ios::out);
  tout = new FileHandler(params->getParam("tout"), std::ios::out);
  generateCorpus(); // generate corpus recursively
  delete params;
  delete sout;
  delete tout;
  return 1;
}
