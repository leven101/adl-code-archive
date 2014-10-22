/* pattern1 looks like: 
0 / 0 
0 / 0 1 
0 / 0 1 2 
0 / 1 
0 / 1 2 
0 / 2 
0 1 / 0 
0 1 / 0 1 
0 1 / 0 1 2 
0 1 / 1 
0 1 / 1 2 
0 1 / 2 
0 1 2 / 0 
0 1 2 / 0 1 
0 1 2 / 0 1 2 
0 1 2 / 1 
0 1 2 / 1 2 
0 1 2 / 2 
1 / 0 
1 / 0 1 
1 / 0 1 2 
...
*/
void traverseGridPattern1() {
  const int srcsize=3, trgsize=3;
  for(int sss=0; sss < srcsize; ++sss) {
    vector<int> src;
    for(int ss=sss; ss < srcsize; ++ss) {
      src.push_back(ss); 
      for(int ttt=0; ttt < trgsize; ++ttt) {
        vector<int> trg;
        for(int tt=ttt; tt < trgsize; ++tt) {
          trg.push_back(tt);
          iterate(src, sit)
            cout << *sit << " ";
          cout << "/ ";
          iterate(trg, tit)
            cout << *tit << " ";
          cout << endl;
        }
      }
    }
  }
}
/* pattern 2 
0 / 0 
0 / 0 1 
0 / 0 1 2 
0 / 1 
0 / 1 2 
0 / 2 
1 / 0 
1 / 0 1 
1 / 0 1 2 
1 / 1 
1 / 1 2 
1 / 2 
2 / 0 
2 / 0 1 
2 / 0 1 2 
2 / 1 
2 / 1 2 
2 / 2 
3 / 0 
3 / 0 1 
3 / 0 1 2 
3 / 1 
3 / 1 2 
3 / 2 
0 1 / 0 
0 1 / 0 1 
0 1 / 0 1 2 
...
*/
void traverseGridPattern2() {
  const int srcsize=4, trgsize=3;
  for(int span=1; span <= srcsize; ++span) { // current span
    for(int i=0; i < srcsize && (i+span <= srcsize); ++i) {
      vector<int> src;
      for(int j=i; j-i < span; ++j) {
        src.push_back(j); 
      }
      for(int ttt=0; ttt < trgsize; ++ttt) {
        vector<int> trg;
        for(int tt=ttt; tt < trgsize; ++tt) {
          trg.push_back(tt);
          iterate(src, sit)
            cout << *sit << " ";
          cout << "/ ";
          iterate(trg, tit)
            cout << *tit << " ";
          cout << endl;
        }
      }
    }
  }
}
/* pattern 3 - large to small
0 1 2 / 0 1 
0 1 2 / 0 
0 1 2 / 1 
0 1 / 0 1 
0 1 / 0 
0 1 / 1 
1 2 / 0 1 
1 2 / 0 
1 2 / 1 
0 / 0 1 
0 / 0 
0 / 1 
1 / 0 1 
1 / 0 
1 / 1 
2 / 0 1 
2 / 0 
2 / 1 
 */
void traverseGridPattern3() {
  const int srcsize=4, trgsize=3;
  for(int sspan=srcsize; sspan > 0; --sspan) { // current span
    for(int i=0; i+sspan <= srcsize; ++i) {
      vector<int> src;
      for(int j=i; j-i < sspan; ++j) {
        src.push_back(j); 
      }
      for(int tspan=trgsize; tspan > 0; --tspan) {
        for(int t=0; t+tspan <= trgsize; ++t) {
          vector<int> trg;
          for(int tt=t; tt-t < tspan; ++tt) {
            trg.push_back(tt);
          }
          iterate(src, sit)
            cout << *sit << " ";
          cout << "/ ";
          iterate(trg, tit)
            cout << *tit << " ";
          cout << endl;
        }
      }
    }
  }
}

