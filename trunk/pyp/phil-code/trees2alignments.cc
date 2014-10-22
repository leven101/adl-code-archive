#include "text.hh"
#include "alignment.hh"
#include "rule.hh"

// STL
#include <iostream>
#include <fstream>

using namespace std;


int main(int argc, char* argv[])
{
  while(cin)
  {
    // read the alignment tree
    RuleTree rule_tree;
    if (!read_tree(cin, rule_tree))
      break;

    Alignment alignment = rule_tree_to_alignment(rule_tree);
    print_alignment(cout, alignment);
    cout << endl;
  }
  return 0;
}
