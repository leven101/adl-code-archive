// g++ aer.cpp -I ../include/ ../file.o -o aer.out
#include "file.h"
#include "types.h"
#include <algorithm>
using namespace std;

int readAlignments(string line, map<size_t, set<size_t> >& algMap) {
	vector<int> vidx;
  Utils::splitToInt(line, vidx, "- ");
  for(size_t i=0; i < vidx.size(); i+=2) {
    size_t srcIdx = vidx[i];
    size_t trgIdx = vidx[i+1];
    algMap[srcIdx].insert(trgIdx);
  }
  return vidx.size()/2; 
}
void phraseAER(int argc, char** argv) {
  string path1 = argv[1];
  string path2 = argv[2];
  FileHandler fin1(path1, ios::in);
  FileHandler fin2(path2, ios::in);
  string line;
  vector<string> v;
  set<string> gram1, gram2, intr1; 
  while(getline(fin1, line)) {
    Utils::splitToStrMD(line, v, "|||");
    string rule = v[1] + " ||| " + v[2]; 
    gram1.insert(rule);
  }
  while(getline(fin2, line)) {
    Utils::splitToStrMD(line, v, "|||");
    string rule = v[1] + " ||| " + v[2]; 
    gram2.insert(rule);
  }
  fin1.close();
  fin2.close();
  insert_iterator<set<string> > insIt (intr1, intr1.begin());
  set_intersection(gram1.begin(), gram1.end(), gram2.begin(), gram2.end(), insIt); 
  // calculate total stats
  float sizeOfIntr = intr1.size();
  float sizeOfGA = gram1.size();
  float sizeOfTA = gram2.size();
  float recall = sizeOfIntr / sizeOfGA; 
  float precision =  sizeOfIntr / sizeOfTA;
  float aer = 1 - (sizeOfIntr*2 / (sizeOfGA + sizeOfTA));
  cout << "R=" << recall << ", P=" << precision << ", AER=" << aer << endl;
}
void grammarDiffs(int argc, char** argv) {
  string path1 = argv[1];
  string path2 = argv[2];
  FileHandler fin1(path1, ios::in);
  FileHandler fin2(path2, ios::in);
  string line;
  vector<string> v;
  vector<float> vf;
  map<string, float> gram1, gram2, diff1, diff2;
  float totPrb_gram1(0), totPrb_gram2(0);
  while(getline(fin1, line)) {
    Utils::splitToStrMD(line, v, "|||");
    string rule = v[0] + " ||| " + v[1] + " ||| " + v[2]; 
    Utils::splitToFloat(v[3], vf, " ");
    gram1[rule] = pow(10, vf[2]);
    totPrb_gram1 += pow(10, vf[2]);
  }
  while(getline(fin2, line)) {
    Utils::splitToStrMD(line, v, "|||");
    string rule = v[0] + " ||| " + v[1] + " ||| " + v[2]; 
    Utils::splitToFloat(v[3], vf, " ");
    gram2[rule] = pow(10, vf[2]);
    totPrb_gram2 += pow(10, vf[2]);
  }
  fin1.close();
  fin2.close();
  set_difference(gram1.begin(), gram1.end(), gram2.begin(), gram2.end(), 
    inserter(diff1, diff1.begin()), gram1.value_comp());
  //cout << "Size of set difference = " << diff1.size() << endl;
  iterate(diff1, dit)
    cout << dit->first << "\t" << ((float)dit->second/totPrb_gram1) << endl;
  set_difference(gram2.begin(), gram2.end(), gram1.begin(), gram1.end(), 
    inserter(diff2, diff2.begin()), gram2.value_comp());
  //cout << "Size of set difference = " << diff2.size() << endl;
  iterate(diff2, dit)
    cerr << dit->first << "\t" << ((float)dit->second/totPrb_gram2) << endl;
}
void calculateAER(char** argv) {
  map<size_t, set<size_t> > goldAlgs;
  map<size_t, set<size_t> > testAlgs;
  map<size_t, set<size_t> >::iterator ait; 
  string goldFile(argv[1]), testFile(argv[2]);
  FileHandler goldIn(goldFile, ios::in);
  FileHandler testIn(testFile, ios::in);
  string line, line2;
	vector<int> vidx, vidx2;
  float sizeOfGA(0), sizeOfTA(0), sizeOfIntr(0);
  int totSnts(0);
	while(getline(goldIn, line)) {
    ++totSnts;
    goldAlgs.clear();
    testAlgs.clear();
    sizeOfGA += readAlignments(line, goldAlgs);
    getline(testIn, line2);
    sizeOfTA += readAlignments(line2, testAlgs);
    iterate(goldAlgs, ga) { // for each src word in gold alignment
      ait = testAlgs.find(ga->first);
      if(ait != testAlgs.end()) { // if this src word is aligned in sampled sentence  
        set<size_t> intersect, s1 = ga->second, s2 = ait->second;
        // get the intersection of aligned words
        insert_iterator<set<size_t> > insIt (intersect, intersect.begin());
        set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(), insIt); 
        sizeOfIntr += intersect.size();
      }
    }
  }
  // calculate total stats
  float recall = sizeOfIntr / sizeOfGA; 
  float precision =  sizeOfIntr / sizeOfTA;
  float aer = 1 - (sizeOfIntr*2 / (sizeOfGA + sizeOfTA));
  cout << "Total sentences in parallel corpus = " << totSnts << endl;
  cout << "R=" << recall << ", P=" << precision << ", AER=" << aer << endl;
  goldIn.close();
  testIn.close();
}

void cutoffSample(char** argv) {
  string alPath = argv[1];
  FileHandler fOldAls(alPath, ios::in);
  FileHandler fNewAls(alPath + ".new", ios::out, false);
  string line;
  while(getline(fOldAls, line)) {
    map<size_t, set<size_t> > algs1, algs2;
    readAlignments(line, algs1);
    set<size_t> srcIdxs;
    iterate(algs1, alit) {  
      srcIdxs.insert(alit->first);
    }
    iterate(srcIdxs, sit) {  // for each source word
      if(algs1.find(*sit) == algs1.end()) continue;  // if deleted
      set<size_t> phrase;
      set<size_t>& tidxs = algs1[*sit]; 
      //size_t sright(0), sleft(10000);
      iterate(algs1, alit) { 
        if(alit->second == tidxs) { //  if matching target side
          phrase.insert(alit->first); // rebuild source side phrase
          //if(alit->first < sleft) sleft = alit->first; // get leftmost index
          //if(alit->first > sright) sright = alit->first; // get rightmost index
        }
      }
      if(phrase.size() + tidxs.size() < 7) { // total number of terminals in phrase
        iterate(phrase,pit) {  // add full phrase
          algs2[*pit] = tidxs;
        }
      }
      else { 
        // add left and right most links only
        //size_t tright(0), tleft(10000);
        //iterate(tidxs, tit) {
          //if(*tit < tleft) tleft = *tit;
          //if(*tit > tright) tright = *tit;
        //}
        //algs2[sleft].insert(tleft);
        //algs2[sright].insert(tright);
        iterate(phrase,pit) // delete nodes
          algs1.erase(*pit);
      }
    }
    // print out new alignments
    iterate(algs2, alit) {
      iterate(alit->second, tit) {
        fNewAls << alit->first << "-" << *tit << " ";
      }
    }
    fNewAls << endl;
  }
  fOldAls.close();
  fNewAls.close();
}

int main(int argc, char** argv) {
  //grammarDiffs(argc, argv);
  if(argc == 2)
    cutoffSample(argv);
  else if(argc == 3) {
    //phraseAER(argc, argv);
    calculateAER(argv);
  }
  else { // display baseline vs. sampler alignment differences 
    string srcPath = argv[1];
    string trgPath = argv[2];
    string algPath = argv[3];
    int to(atoi(argv[4])), at(0);
    FileHandler fsrc(srcPath, ios::in);
    FileHandler ftrg(trgPath, ios::in);
    FileHandler falg(algPath, ios::in);
    string sline, tline, aline;
    vector<string> svec, tvec;
    map<size_t, set<size_t> > mAlgs;
    while(getline(fsrc, sline) && (++at <= to)) {
      getline(ftrg, tline);
      getline(falg, aline);
      if(at == to) {
        Utils::splitToStr(sline, svec, " ");
        Utils::splitToStr(tline, tvec, " ");
        string aline2 = aline.substr();
        mAlgs.clear();
        readAlignments(aline, mAlgs);
        cout << sline << endl << tline << endl << aline2 << endl;
        iterate(mAlgs, algitr) {
          cout << svec[algitr->first] << ":";
          iterate(algitr->second, trgitr) {
            cout << tvec[*trgitr] << ", ";
          }
          cout << endl;
        }
      }
    }
    fsrc.close();
    ftrg.close();
    falg.close();
  }
  return 1;
}
