#ifndef itg2trees_h 
#define itg2trees_h 

#include <stack>
#include "multiNTs.h"
#include "sampler_utils.h"

string itglinee = "([-1]->[0,1|||-1,-3] ([-1]->[0,1|||-1,-2] ([-1]->[0,1|||-1,-2] ([-1]->[0,1|||-1,-3] ([-1]->[0,1|||-2,-2] ([-2]->[1,0|||-3,-2] [-3]->[121|||18] ([-2]->[1,0|||-3,-3] [-3]->[219|||1264] [-3]->[170|||184])) ([-2]->[1,0|||-3,-2] [-3]->[367|||325] ([-2]->[1,0|||-3,-1] [-3]->[40|||238] ([-1]->[0,1|||-3,-3] [-3]->[9|||27] [-3]->[2176|||2727])))) [-3]->[368|||145]) ([-2]->[1,0|||-1,-2] ([-1]->[0,1|||-1,-2] ([-1]->[0,1|||-2,-3] ([-2]->[1,0|||-3,-2] [-3]->[2725|||2728] ([-2]->[1,0|||-3,-1] [-3]->[13|||76] ([-1]->[0,1|||-1,-3] ([-1]->[0,1|||-3,-3] [-3]->[2723|||1070] [-3]->[159|||50]) [-3]->[2724|||2729]))) [-3]->[156|||50]) ([-2]->[1,0|||-3,-1] [-3]->[2726|||2730] ([-1]->[0,1|||-2,-4] ([-2]->[1,0|||-3,-3] [-3]->[134|||76] [-3]->[2697|||1770]) [-4]->[|||151]))) ([-2]->[1,0|||-3,-3] [-3]->[11|||18] [-3]->[2722|||269]))) ([-2]->[1,0|||-1,-1] ([-1]->[0,1|||-3,-4] [-3]->[16|||456] [-4]->[|||411]) ([-1]->[0,1|||-3,-3] [-3]->[931|||2731] [-3]->[932|||1009]))) [-3]->[17|||17])";


class ITG2Trees {
public:
  ITG2Trees(MultiNT* mnt): mnt_(mnt) {
    loadVocab();
    FileHandler fin(mnt_->params_->getParam("init-aligns"), std::ios::in); 
    string line;
    int tot=0;
    while(getline(fin, line) && (++tot <= mnt_->src_sents.size())) {
      tree<node_t>* tr = itg2tree(line);
      int sstr=0, send=0;
      getSrcSpans(tr->begin(), sstr, send); 
      //handleInsertionRules(tr); // for block sampler 
      mnt_->trees_.push_back(tr);
    }
    mnt_->saveTreeState();
    fin.close();
  }
private: 
  MultiNT* const mnt_; 
  map<int, string> itgvcb_;
  void getSrcSpans(const tree<node_t>::iterator node, int& sstr, int& send) {
    node->sstr = sstr;
    for(int i=0; i < node->rule.t.size(); ++i) {
      if(!isTerm(node->rule.t[i])) {
        tree<node_t>::iterator child = std::find(node.begin(), node.end(), node->rule.t[i]);
        getSrcSpans(child, sstr, send);
      }
    }
    for(int s=0; s < node->rule.s.size(); ++s) {
      if(isTerm(node->rule.s[s])) {
        ++send;
      }
    }
    // done recursion
    node->send = send;
    sstr = send;
  }
  void loadVocab() {
    // load itg vocab
    string line;
    vector<string> v1;
    FileHandler fin("/data/taipan/ablev/data/ur-en/itg-trees/ur-en.t.itg.rl10.vocab", std::ios::in); 
    while(getline(fin, line)) {
      Utils::splitToStr(line, v1, " ");
      itgvcb_[atoi(v1[0].c_str())] = v1[1]; 
    }
    fin.close();
    cerr << "itgvcb_ has " << itgvcb_.size() << " items" << endl;
  }
  tree<node_t>* itg2tree(const string& itgline) {
    tree<node_t>* tr = new tree<node_t>;
    std::stack<int> nt_stack;
    map<int, int> node_prts;
    vector<string> v_rules;
    Utils::splitToStr(itgline, v_rules, " ");
    // build scfg rules from itg rules
    int nt_lbl = -1;
    nt_stack.push(nt_lbl);
    iterate(v_rules, vit) {
      int label = nt_stack.top();
      nt_stack.pop();
      hieroRule hr;
      hr.NT = 0;
      vector<string> vs, vs2;
      vector<int> vi;
      bool b_new_node(true);
      Utils::splitToStrMD(*vit, vs, "->");
      Utils::trim(vs[0], "([]");
      int ruleType = atoi(vs[0].c_str());
      switch(ruleType) {
        case -1: { // monotonic ITG rule
          for(int i=0; i < 2; ++i) {
            hr.s.push_back(--nt_lbl);
            hr.t.push_back(nt_lbl);
            node_prts[nt_lbl] = label;
          }
          nt_stack.push(nt_lbl);
          nt_stack.push(nt_lbl+1);
          break;
        }
        case -2: { // reordered ITG rule
          for(int i=0; i < 2; ++i) {
            hr.s.push_back(--nt_lbl);
            hr.t.push_back(nt_lbl);
            node_prts[nt_lbl] = label;
          }
          int tmp = hr.s[0];
          hr.s[0] = hr.s[1];
          hr.s[1] = tmp;
          nt_stack.push(nt_lbl);
          nt_stack.push(nt_lbl+1);
          break;
        }
        case -3: { // lexical itg rule 
          Utils::trim(vs[1], "[]()");
          Utils::splitToStrMD(vs[1], vs2, "|||");
          hr.s.push_back(getWordID(itgvcb_[atoi(vs2[0].c_str())], srcVcb_));
          hr.t.push_back(getWordID(itgvcb_[atoi(vs2[1].c_str())], trgVcb_));
          break;
        }
        default: { // case 4
          int type4 = getType4Rule(vs[1]);
          switch(type4) {
            case -1: { // double null rule 
              tree<node_t>::iterator prt = std::find(tr->begin(), tr->end(), node_prts[label]);
              // erase parent label for both source and target
              for(int i=0; i < prt->rule.s.size(); ++i) {
                if(prt->rule.s[i] == label) 
                  prt->rule.s.erase(prt->rule.s.begin()+i);
              }
              for(int i=0; i < prt->rule.t.size(); ++i) {
                if(prt->rule.t[i] == label) 
                  prt->rule.t.erase(prt->rule.t.begin()+i);
              }
              // add place holders for child null aligned rules
              for(int i=0; i < 2; ++i) {
                prt->rule.s.push_back(--nt_lbl);
                prt->rule.t.push_back(nt_lbl);
                node_prts[nt_lbl] = prt->label;
              }
              nt_stack.push(nt_lbl);
              nt_stack.push(nt_lbl+1);
              b_new_node = false;
              break;
            }
            case -2: {  // trg word is null aligned
              tree<node_t>::iterator prt = std::find(tr->begin(), tr->end(), node_prts[label]);
              hieroRule& prtRule = prt->rule;
              Utils::trim(vs[1], "[|]");
              int tid = getWordID(itgvcb_[atoi(vs[1].c_str())], trgVcb_); 
              for(int i=0; i < prtRule.t.size(); ++i) {
                if(prtRule.t[i] == label) 
                  prtRule.t[i] = tid;
              }
              for(int i=0; i < prtRule.s.size(); ++i) {
                if(prtRule.s[i] == label) 
                  prtRule.s.erase(prtRule.s.begin()+i);
              }
              b_new_node = false;
              break;
            }
            case -3: { // src word is null aligned
              tree<node_t>::iterator prt = std::find(tr->begin(), tr->end(), node_prts[label]);
              hieroRule& prtRule = prt->rule;
              Utils::trim(vs[1], "[|]");
              int sid = getWordID(itgvcb_[atoi(vs[1].c_str())], srcVcb_); 
              for(int i=0; i < prtRule.s.size(); ++i) {
                if(prtRule.s[i] == label) 
                  prtRule.s[i] = sid;
              }
              for(int i=0; i < prtRule.t.size(); ++i) {
                if(prtRule.t[i] == label) 
                  prtRule.t.erase(prtRule.t.begin()+i);
              }
              b_new_node = false;
              break;
            }
            default: { // lexical rule
              Utils::trim(vs[1], "[]()");
              Utils::splitToStrMD(vs[1], vs2, "|||");
              hr.s.push_back(getWordID(itgvcb_[atoi(vs2[0].c_str())], srcVcb_));
              hr.t.push_back(getWordID(itgvcb_[atoi(vs2[1].c_str())], trgVcb_));
            }
          }
        } //  end default case
      } // end switch
      if(b_new_node) {
        node_t node;
        node.label = label;
        node.rule = hr;
        if(tr->empty()) {
          tr->set_head(node);
        }
        else {
          tree<node_t>::iterator prt = std::find(tr->begin(), tr->end(), node_prts[label]);
          tr->append_child(prt, node); // add new node to tree
        }
      }
    }
    tr->begin()->nodeIndex = nt_lbl;
    return tr;
  }
  int getType4Rule(string rule) {
    Utils::trim(rule, "()");
    int type(0);
    if(rule == "[0,1|||-4,-4]") {
      type = -1;
    }
    else if(rule.find("[|||") == 0) {
      type = -2;
    }
    else if(rule.rfind("|||]") == rule.size()-4) {
      type = -3;
    }
    else type = -4;
    return type;
  }
  void handleInsertionRules(tree<node_t>* tr) {
    tree<node_t>::iterator trit = tr->begin();
    tree<node_t>::iterator childnode;
    while(trit != tr->end()) {
      if(isInsertionRule(trit->rule)) {
        //cout << "Before: " << endl;
        //mnt_->printTree(tr);
        hieroRule& insrule = trit->rule; // reference 
        do {
          //cerr << "insertion rule: " << insrule << endl;
          childnode = std::find(tr->begin(), tr->end(), insrule.s[0]); // get its child 
          //cerr << "child rule: " ;
          //mnt_->printTreeNode(childnode);
          assert(childnode != tr->end());
          const hieroRule& childrule = childnode->rule;
          int tidx(-1);
          tidx = std::find(insrule.t.begin(), insrule.t.end(), insrule.s[0]) - insrule.t.begin(); 
          insrule = combineChildParent(insrule, childrule, 0, tidx, false);
          // attach children of childnode to new parent
          iterate(childrule.s, sit) {
            if(!isTerm(*sit)) {
              tree<node_t>::iterator newchild = std::find(tr->begin(), tr->end(), *sit); // get its children 
              assert(newchild != tr->end());
              tr->append_child(trit, newchild); 
            }
          }
          //cerr << "deleting : " << childrule << endl;
          tr->erase(childnode);
          //cerr << "new rule: " << insrule << endl;
        } while(isInsertionRule(insrule));
        //cout << "After: " << endl;
        //mnt_->printTree(tr);
        trit = tr->begin(); // reset tree iterator
      }
      else ++trit;
    }
  }
};

#endif
