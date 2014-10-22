//g++ reorder.cpp -I ~/src/include/ ~/src/file.o
#include "file.h"
#include "utils.h"
#include "types.h"

using namespace std;

int main(int argc, char** argv) {
  FileHandler fin1(argv[1], ios::in);
  FileHandler fin2(argv[2], ios::in);
  FileHandler fin3(argv[3], ios::in);
  string line1, line2, line3;
  vector<string> v, v2, v3;
  while(getline(fin1, line1)) {
    getline(fin2, line2);
    getline(fin3, line3);
    v.push_back(line1);
    v2.push_back(line2);
    v3.push_back(line3);
  }
  FileHandler fout1("./urdu.reordered", ios::out, false);
  FileHandler fout2("./english.reordered", ios::out, false);
  FileHandler fout3("./hmm.reordered", ios::out, false);
  while(v.size() > 0) {
    assert(v2.size() == v.size());
    assert(v3.size() == v.size());
    int rnd = rand() % v.size();
    fout1 << v[rnd] << endl;
    fout2 << v2[rnd] << endl;
    fout3 << v3[rnd] << endl;
    v.erase(v.begin() + rnd);
    v2.erase(v2.begin() + rnd);
    v3.erase(v3.begin() + rnd);
  }
  fin1.close();
  fin2.close();
  fin3.close();
  fout1.close();
  fout2.close();
  fout3.close();
  return 1;
}

