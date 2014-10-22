#include <algorithm>
#include <string>
#include <iostream>
#include <types.h>
#include "tree.hh"
#include <pthread.h>
#include <cmath>
#define NUM_THREADS     5
using namespace std;
void printTree(tree<string>& tr) {
  tree<string>::iterator sib2=tr.begin();
  tree<string>::iterator end2=tr.end();
  while(sib2!=end2) {
    for(int i=0; i<tr.depth(sib2)-2; ++i) 
      cout << "--";
    cout << (*sib2) << endl;
    ++sib2;
  }
  cout << endl;
}
void treeStuff() {
  tree<string> tr;
  tree<string>::iterator top, one, two, loc, banana;

  top=tr.begin();
  one=tr.insert(top, "one");
  two=tr.append_child(one, "two");
  tr.append_child(two, "apple");
  banana=tr.append_child(two, "banana");
  tr.append_child(banana,"cherry");
  tr.append_child(two, "peach");
  tr.append_child(one,"three");
  printTree(tr);
  const tree<string>::iterator it1 = tr.begin();
  tree<string>::iterator it2 = tr.end();
  while(it2 != it1) {
    --it2;
    cout << *it2 << endl;
  }
  cout << ":)\n";
  cout << *it1 << endl;

  loc=find(tr.begin(), tr.end(), "two");
  cout << "Printing only siblings of " << *loc << endl;
  if(loc!=tr.end()) {
    tree<string>::sibling_iterator sib=tr.begin(loc);
    while(sib!=tr.end(loc)) {
      cout << (*sib) << endl;
      ++sib;
    }
    cout << endl;
  }
  cout << endl << "my thing\n";
  loc=find(tr.begin(), tr.end(), "apple");
  assert(loc != tr.end());
  one = tr.append_child(find(tr.begin(), tr.end(), "cherry"), loc);
  printTree(tr);
  if(one == loc) 
    cout << "they do match so not a copy\n";
  else cout << "inserts a copy\n";
  tr.erase(loc);
  printTree(tr);
  cout << endl;
  it2 = tr.end();
  while(it2 != it1) {
    --it2;
    cout << *it2 << endl;
  }
}

void* PrintHello(void *threadid)
{
  long tid;
  tid = (long)threadid;
  printf("Hello World! It's me, thread #%ld!\n", tid);
  pthread_exit(NULL);
}

int main (int argc, char *argv[])
{
  treeStuff();
  return 1;
  pthread_t threads[NUM_THREADS];
  int rc;
  long t;
  for(t=0; t<NUM_THREADS; t++){
    printf("In main: creating thread %ld\n", t);
    rc = pthread_create(&threads[t], NULL, PrintHello, (void *)t);
    if (rc){
      printf("ERROR; return code from pthread_create() is %d\n", rc);
      exit(-1);
    }
  }
  /* Last thing that main() should do */
  pthread_exit(NULL);
}
