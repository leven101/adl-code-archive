#include "hiero.h"
using std::make_pair;
double get_p0(const hieroRule* hr, const bool debug) {
  static const double usp = log(1.0 / (double)srcVcb_->size());
  static const double utp = log(1.0 / (double)trgVcb_->size());
  const int sl = hr->s.size();
  const int tl = hr->t.size();
  const int arity = hr->arity(); 
  const int nst = sl - arity;
  const int ntt = tl - arity;
  // probability of source length
  static const float pm = atof(params_->getParam("poisson-mean").c_str());
  double psl = log_poissDist(sl, pm);
  // probability of source side configuration
  double tpp = pow(termPenalty(), sl); // term penalty param
  double ssc = (nst * log(tpp)) + (arity * log(1-tpp));
  double ptl = log_poissDist(ntt, nst+0.01); 
  double tsc(0), n = arity;
  for(int i=0; i < tl; ++i) {
    double ntp = n / (tl-i); // non-term param
    if(isTerm(hr->t[i])) {
      tsc += log(1-ntp); 
    }
    else {
      tsc += log(ntp);   // probability of drawing NT
      tsc -= log(n);    // uniform permutation of NTs
      --n;
    }
  }
  static const wordID_t nullID = getWordID("NULL", NULL);
  // model1 p(t|s) p(s)
  double s2t = nst * usp;   // p(s)
  if(ntt > 0) { //get m1 prob for p(t|s)
    iterate(hr->t, t) {
      if(isTerm(*t)) { // for each target terminal 
        double m1prb = m1_[make_pair(nullID,*t)].first;
        double tsum = Log<double>::add(Log<double>::zero(), m1prb); // get P(t|NULL)
        if(nst > 0) {
          iterate(hr->s, s) { // for each source word
            if(isTerm(*s)) {
              m1prb =  m1_[make_pair(*s,*t)].first;
              tsum = Log<double>::add(tsum, m1prb); // get P(t|s)
            }
          }
        }
        s2t += tsum;
      }
    }
    s2t -= log(pow(nst+1, ntt));  
  }
  // model 1 p(s|t)p(t)
  double t2s = ntt * utp; // p(t)
  if(nst > 0) {
    iterate(hr->s, s) {
      if(isTerm(*s)) { // for each source terminal 
        double m1prb = m1_[make_pair(*s, nullID)].second;
        double ssum = Log<double>::add(Log<double>::zero(), m1prb); // get P(s|NULL)
        if(ntt > 0) {
          iterate(hr->t, t) { // for each target word
            if(isTerm(*t)) {
              m1prb =  m1_[make_pair(*s,*t)].second;
              ssum = Log<double>::add(ssum, m1prb); // get sum(P(s|t))
            }
          }
        }
        t2s += ssum;
      }
    }
    t2s -= log(pow(ntt+1, nst)); // normalize 
  }
  double joint_ts = Log<double>::add(s2t, t2s) - log(2);

  double p0 = psl + ptl + ssc + tsc + joint_ts;
  if(debug) {
    cout << "prior for rule " << *hr << " = psl: " << psl << "\tptl: " << ptl
      << "\tssc: " << ssc << "\ttsc: " << tsc << "\tjoint_ts: " << joint_ts
      << "\tp0: " << p0 << endl;
  }
  return p0; 
}
void loadM1Params() {
  cerr << "Loading model1 parameters (src->trg)...\n";
  string line; 
  wordID_t sid, tid;
  vector<string> v;
  FileHandler fm1(params_->getParam("model1"), std::ios::in); 
  while(getline(fm1, line)) {
    Utils::splitToStr(line, v, " ");
    sid = getWordID(v[0], srcVcb_);
    tid = getWordID(v[1], trgVcb_);
    m1_[make_pair(sid, tid)].first = log(atof(v[3].c_str()));
  }
  fm1.close();
  cerr << "Loading model1 parameters (trg->src)...\n";
  FileHandler fm2(params_->getParam("model1-inv"), std::ios::in); 
  while(getline(fm2, line)) {
    Utils::splitToStr(line, v, " ");
    tid = getWordID(v[0], trgVcb_);
    sid = getWordID(v[1], srcVcb_);
    m1_[make_pair(sid, tid)].second = log(atof(v[3].c_str()));
  }
  fm2.close();
  cerr << "Finished loading model1 parameters.\n";
}
void loadSentences() {
  cerr << "Loading source sentences...\n";
  size_t maxSnts = atoi(params_->getParam("maxSnts").c_str());
  FileHandler fsrc(params_->getParam("source"), std::ios::in);
  string line;
  vector<string> vline; 
  float srcwords(0), trgwords(0);
  while(getline(fsrc, line)) {
    srcwords += Utils::splitToStr(line,vline, " ");
    vector<wordID_t> lineIDs; 
    iterate(vline, w) { 
      lineIDs.push_back(getWordID(*w, srcVcb_));
    }
    src_sents.push_back(lineIDs);
    if(src_sents.size() == maxSnts) break;
  }
  fsrc.close();
  srcVcb_->makeClosed();
  cerr << "Loading target sentences...\n";
  FileHandler ftrg(params_->getParam("target"), std::ios::in);
  while(getline(ftrg, line) && (trg_sents.size() != src_sents.size())) {
    trgwords += Utils::splitToStr(line,vline, " ");
    vector<wordID_t> lineIDs; 
    iterate(vline, w) { 
      lineIDs.push_back(getWordID(*w, trgVcb_));
    }
    trg_sents.push_back(lineIDs);
  }
  ftrg.close();
  assert(trg_sents.size() == src_sents.size());
  trgVcb_->makeClosed();
  cout << "Source/target words: " << srcwords << " / " << trgwords << endl;
  cout << "Average source/target sentence length: ";
  cout << (srcwords/(float)src_sents.size()) << " / "; 
  cout << (trgwords/(float)trg_sents.size()) << endl;
  cerr << "Loaded " << src_sents.size() << " sentence pairs.\n";
}
void loadData() {
  srcVcb_ = new Vocab(false);
  trgVcb_ = new Vocab(false);
  loadM1Params();
  loadSentences();
  cerr << "Source vocabularly size = " << srcVcb_->size() << endl;
  cerr << "Target vocabularly size = " << trgVcb_->size() << endl;
}
void fullSentenceInit() {
  for(size_t i=0; i < src_sents.size(); ++i) {
    bool debug = (int)i == snt2debug_;
    hieroRule sntRule;
    sntRule.s.assign(src_sents[i].begin(), src_sents[i].end());
    sntRule.t.assign(trg_sents[i].begin(), trg_sents[i].end());
    tree<node_t>* tr = new tree<node_t>;
    node_t topNode;
    topNode.rule = sntRule;
    topNode.label = 1; // set label;
    topNode.NTIndex = 0; // total NTs in tree;
    topNode.sstr = 0;
    topNode.send = src_sents[i].size();
    //topNode.tstr = 0;
    //topNode.tend = trg_sents[i].size();
    tr->set_head(topNode); // set top node to sentence rule 
    //increment(sntRule, exp(get_p0(&sntRule, debug)));
    trees_.push_back(tr);
    if(debug) printTree(i);
  } 
}
void initPYP() {
  float a = atof(params_->getParam("discount").c_str()),
    b = atof(params_->getParam("strength").c_str());
  restr_ = new PYP<hieroRule, hieroRule>(a, b);
  bool restart = params_->getParam("restart") ==
      Parameters::kTrueValue ? true: false;
  if(restart) {
    readTreeState();
    assert(trees_.size() == src_sents.size());
  }
  else if(noInit_) // add each sentence as a rule 
    fullSentenceInit();
  else { // add probable word alignments
    initPhrs_ = params_->getParam("phr-init") ==
      Parameters::kTrueValue ? true: false;
    FileHandler m1alg(params_->getParam("init-aligns"), std::ios::in);
    string line;
    vector<int> v;
    for(size_t i=0; i < src_sents.size(); ++i) {
      bool debug = (int)i == snt2debug_;
      vector<int> alGrid(trg_sents[i].size(), -1);
      vector<bool> srcIdxUsed(src_sents[i].size(), false);
      getline(m1alg, line);
      Utils::splitToInt(line, v, "- ");
      for(int j=0; j < (int)v.size(); j+=2) {
        int srcIdx = v[j];
        int trgIdx = v[j+1];
        if(srcIdxUsed[srcIdx] == false || (srcIdx == alGrid[trgIdx-1])) {
          alGrid[trgIdx] = srcIdx;
          srcIdxUsed[srcIdx] = true;
        }
      }
      tree<node_t>* tr = buildXRuleGrid(alGrid, i); // expects vector trgIdx->srcIdx
      tree<node_t>::iterator trit = tr->begin();
      while(trit != tr->end()) {
        const hieroRule& hr = trit->rule; 
        double p0 = get_p0(&hr, debug);
        hieroRule transf = transform(hr);
        increment(transf, exp(p0)); // add rule to rest
        ++trit;
      }
      trees_.push_back(tr);
      if(debug) printTree(i);
    } // end corpus
    m1alg.close();
  }
  printAlignments(-1);
  restr_->debug_info(cout);
  cerr << "Initial data LL = " << dataLLAndOtherStuff(-1) << endl;
}
tree<node_t>* buildXRuleGrid(const vector<int>& grid, const int snt) {
  tree<node_t>* tr = new tree<node_t>;
  node_t dummy;
  tree<node_t>::iterator top = tr->set_head(dummy);
  hieroRule topX;
  vector<wordID_t>& src = src_sents[snt];
  vector<wordID_t>& trg = trg_sents[snt];
  assert(grid.size() == trg.size());
  map<pair<int, int> , int> scover; // helper structure since source rules are found in arbitrary order 
  int ntIdx(-1);
  if(snt==snt2debug_) { // debug
    for(int i=0; i < (int)grid.size(); ++i) {  // bottom up
      cout << "grid[" << i << "] = " << grid[i] << endl;
    }
  }
  for(int i = grid.size()-1; i >= 0; --i) {
    if(grid[i] != -1) {
      if(initPhrs_) {
        for(int j=0; j <= i; ++j) {  // find 'initial' phrase pairs
          if(grid[j] == -1) {
            continue;
          }
          bool ok(true);
          vector<wordID_t> tmp(grid.begin()+j, grid.begin()+i+1); // get this subsequence
          std::sort(tmp.begin(), tmp.end()); // sort to check for contiguous sequence
          tmp.resize(std::unique(tmp.begin(), tmp.end()) - tmp.begin()); // only keep unique indices
          for(int k=1; k < (int)tmp.size(); ++k) {
            if(tmp[k] != tmp[k-1]+1) { // check in order  
              ok = false;
              break;
            }
          }
          if(ok) { // found longer phrase
            hieroRule X;
            X.t.assign(trg.begin()+j, trg.begin()+i+1); // add target  
            for(int k=0; k < (int)tmp.size(); ++k) { // add source side
              X.s.push_back(src[tmp[k]]);
            }
            scover[make_pair(tmp.front(), tmp.back())] = ntIdx; // source indexes covered by nt
            node_t node;
            node.rule = X;
            node.label = ntIdx;
            //node.tstr = j;
            //node.tend = i+1;
            node.sstr = tmp.front();
            node.send = tmp.back() + 1;
            tr->append_child(top, node); // add rule to tree
            topX.t.insert(topX.t.begin(), ntIdx--); // add nt to target side of full sentence rule
            i=j; // only get longest phrase pair
            break; 
          }
        }
      }
      else {
        hieroRule X;
        X.s.push_back(src[grid[i]]);
        X.t.push_back(trg[i]);
        int j=i-1;
        while((j > -1) && (grid[j] == grid[i])) {
          X.t.insert(X.t.begin(), trg[j]);
          --j;
        }
        i = j+1; 
        node_t node;
        node.rule = X;
        node.label = ntIdx;
        //node.tstr = i;
        //node.tend = i+1;
        node.sstr = grid[i];
        node.send = grid[i] + 1; 
        tr->append_child(top, node); // add rule to tree
        scover[make_pair(grid[i], grid[i])] = ntIdx; // source word index covered by nt
        topX.t.insert(topX.t.begin(), ntIdx--); // add nt to target side of full sentence rule
      }
    }
    else { // nonterminal
      topX.t.insert(topX.t.begin(), trg[i]); // building backwards
    }
  }
  // have full target sentence rule -- need source 
  int curIdx(0);
  iterate(scover, scit) {
    while(curIdx < scit->first.first) {
      topX.s.push_back(src[curIdx++]);
    }
    topX.s.push_back(scit->second);
    curIdx = scit->first.second+1;
  }
  while(curIdx < (int)src.size()) { // add remaining target words
    topX.s.push_back(src[curIdx++]);
  }
  node_t topNode;
  topNode.rule = topX;
  topNode.label = 1; // set label;
  topNode.NTIndex = ntIdx+1; // total NTs in tree;
  topNode.sstr = 0;
  topNode.send = src.size();
  //topNode.tstr = 0;
  //topNode.tend = trg.size();
  top.node->data = topNode; // set top node to sentence rule 
  return tr;
}
bool sampleAll(hieroRule& rule, const vector<pair<vector<int>, int> >& trgSpans, 
  const int sIdx, const double prtSpanPrb, tree<node_t>::iterator node, tree<node_t>* fullTree, 
  const float annealT, const bool debug, const int ss, const int se) {
  bool change(false);
  vector<double> probs;
  map<int, pair<double, double> > mchildPrbs;
  double allChildPrbs(0);
  mchildPrbs = cacheChildPrbs(node, &allChildPrbs, debug) ;
  double currConfigPrb = (prtSpanPrb + allChildPrbs) * annealT; // combine and anneal 
  probs.push_back(currConfigPrb);
  if(debug)
    cout << endl << node->rule << "\t" << probs.back() << endl;
  // for each rule in trgSpan, get probability
  iterate(trgSpans, ts) {
    rule.t = ts->first;
    if(debug) cout << rule << "\t";
    if(alreadyOn(rule)) {
      tree<node_t>::iterator child = std::find(node.begin(), node.end(), rule.s[0]);
      hieroRule ocr = child->rule;
      // combine child with parent
      hieroRule npr = combineChildParent(node->rule, ocr, sIdx, ts->second, debug);
      // get probs for potential rules 
      double nprP0 = get_p0(&npr, debug);
      double nprPrb = restr_->log_prob(transform(npr), nprP0) * annealT; // anneal 
      // get prob of all children minus this one
      double remChildPrbs = allChildPrbs - mchildPrbs[rule.s[0]].first;
      double combPrb = (remChildPrbs + nprPrb) * annealT; // combine and anneal 
      probs.push_back(Log<double>::add(probs.back(), combPrb));
      if(debug) cout << "[" << combPrb << ", ";
    }
    else {
      // split out rule from parent
      int lbl = fullTree->begin()->NTIndex; // get current label from full tree
      hieroRule npr = newParentRule1(node->rule, rule, sIdx, ts->second, lbl-1); // create new parent rule 
      // get probs for potential rules 
      double nprP0 = get_p0(&npr, debug);
      double ncrP0 = get_p0(&rule,debug);
      double nprPrb = restr_->log_prob(transform(npr), nprP0);
      double ncrPrb = restr_->log_prob(transform(rule), ncrP0); // assume independence for now
      double combPrb = (nprPrb + ncrPrb + allChildPrbs) * annealT; // normalize and anneal 
      probs.push_back(Log<double>::add(probs.back(), combPrb));
      if(debug) cout << "[" << combPrb << ", ";
    }
    if(debug) cout << probs.back() << "]" << endl;
  }
  double rSample = probs.back() + log(dRand());
  vector<double>::iterator pit = std::lower_bound(probs.begin(), probs.end(), rSample);
  int index = int(pit - probs.begin());
  if(debug) cout << "Selected index " << index << endl;
  if(index > 0) { 
    change = true;
    rule.t = trgSpans[index-1].first;
    int tStr = trgSpans[index-1].second;
    if(alreadyOn(rule)) {
      tree<node_t>::iterator child = std::find(node.begin(), node.end(), rule.s[0]);
      hieroRule ocr = child->rule;
      decrement(transform(ocr));
      // combine child with parent
      hieroRule npr = combineChildParent(node->rule, ocr, sIdx, tStr, debug);
      // restructure tree
      node->rule = npr; 
      fullTree->reparent(node, child); // move any children of 'child' to parent 
      fullTree->erase(child); 
    }
    else {
      int& lbl = fullTree->begin()->NTIndex; // get current label from full tree
      hieroRule npr = newParentRule1(node->rule, rule, sIdx, tStr, lbl-1); // create new parent rule 
      updateTree(node, fullTree, npr, rule, --lbl, ss, se);
      increment(transform(rule), exp(get_p0(&rule, debug))); // add new child rule to restaurant
    }
  }
  return change;
}
void* sampleMN(void* threadArgs) {
  threadData_t* td = (threadData_t*)threadArgs; 
  set<int> srcNTCache, trgNTCache; // thread specific caching
  int lftRuleCache; 
  for(int t = td->sntBegin; t < td->sntEnd; ++t) {
    tree<node_t>::iterator trit = trees_[t]->begin();
    bool debug = t == snt2debug_;
    while(trit != trees_[t]->end()) { // iterate top-down through all legal spans in alignment grid 
      hieroRule prtRule = trit->rule, prtTrsfm = transform(prtRule); 
      if((!td->skipDecr) || (trit != trees_[t]->begin())) 
        decrement(prtTrsfm); // remove current rule from restr configuration
      double prtP0 = get_p0(&prtRule, debug);
      double prtProb = restr_->log_prob(prtTrsfm, prtP0) * td->anneal;
      vector<int> symSpans = precomputeSpans(trit, debug); 
      int maxSrcSpan = std::min(trit->send - trit->sstr, td->maxspan);
      for(int srcspan=maxSrcSpan; srcspan > 0; --srcspan) { // top level loop from biggest to smallest 
        for(int srIdxStart=0; srIdxStart < (int)prtRule.s.size(); ++srIdxStart) { // start index of rule lhs 
          bool newSpan(false);
          hieroRule rule;
          bool allNewSrc(true); 
          int sidx=srIdxStart, ss, se;
          int curspan(0);
          bool checkTarget(false);
          do { // build source rule
            rule.s.push_back(prtRule.s[sidx]); //add to the rule
            curspan += (symSpans[sidx+1] - symSpans[sidx]);
            if(curspan == srcspan) {
              ss = symSpans[srIdxStart];
              se = ss + curspan;
              checkTarget = true; // sample something
              break;
            }
          } while(curspan < srcspan && (++sidx < (int)prtRule.s.size()));
          if(checkTarget) { // build all target side rules
            vector<pair<vector<int>, int> > trgSpans;
            int maxTrgSpan = std::min((int)prtRule.t.size(), td->maxspan);
            for(int trgspan = maxTrgSpan; trgspan > 0; --trgspan) {
              bool allNewTrg(true);
              for(int tIdxStart=0; tIdxStart + trgspan <= (int)prtRule.t.size(); ++tIdxStart) {
                rule.t.clear();
                for(int tidx=tIdxStart; tidx-tIdxStart < trgspan; ++tidx) {
                  rule.t.push_back(prtRule.t[tidx]);
                }
                if(validRule(rule, prtRule, &allNewSrc, &allNewTrg, 
                  srcNTCache, trgNTCache, &lftRuleCache)) {
                  trgSpans.push_back(make_pair(rule.t, tIdxStart));
                } // end validRule
              } // end tIdxStart
            } // end trgspan
            if(trgSpans.size()) {
              newSpan = sampleAll(rule, trgSpans, srIdxStart, 
                prtProb, trit, trees_[t], td->anneal, debug, ss, se); 
            }
            if(newSpan) {
              prtRule = trit->rule;
              prtTrsfm = transform(prtRule);
              prtP0 = get_p0(&prtRule, debug);
              prtProb = restr_->log_prob(prtTrsfm, prtP0) * td->anneal;
              symSpans = precomputeSpans(trit, debug);
              //if(alreadyOn(rule))
                //--srIdxStart; //keep src start index the same
            }
          }  // end checkTarget
        }  // end srIdxStart 
      }  // end maxSrcSpan
      increment(prtTrsfm, exp(prtP0)); // add back current rule into restr
      ++trit;
    } // end sentence 
    if(debug) printTree(t);
  } // end batch of sentences 
  pthread_exit((void*) threadArgs);
}
void* sample(void* threadArgs) {
  threadData_t* td = (threadData_t*)threadArgs; 
  set<int> srcNTCache, trgNTCache; // thread specific caching
  int lftRuleCache; 
  for(int t = td->sntBegin; t < td->sntEnd; ++t) {
    tree<node_t>::iterator trit = trees_[t]->begin();
    bool debug = t == snt2debug_;
    while(trit != trees_[t]->end()) { // iterate top-down through all legal spans in alignment grid 
      hieroRule prtRule = trit->rule, prtTrsfm = transform(prtRule); 
      if((!td->skipDecr) || (trit != trees_[t]->begin())) 
        decrement(prtTrsfm); // remove current rule from restr configuration
      double prtP0 = get_p0(&prtRule, debug);
      double prtProb = restr_->log_prob(prtTrsfm, prtP0) * td->anneal;
      vector<int> symSpans = precomputeSpans(trit, debug); 
      int maxSrcSpan = std::min(trit->send - trit->sstr, td->maxspan);
      for(int srcspan=maxSrcSpan; srcspan > 0; --srcspan) { // top level loop from biggest to smallest 
        for(int srIdxStart=0; srIdxStart < (int)prtRule.s.size(); ++srIdxStart) { // start index of rule lhs 
          bool newPrtSpan(false);
          hieroRule rule;
          bool allNewSrc(true); 
          int sidx=srIdxStart, ss, se;
          int curspan(0);
          bool checkTarget(false);
          do { // build source rule
            rule.s.push_back(prtRule.s[sidx]); //add to the rule
            curspan += (symSpans[sidx+1] - symSpans[sidx]);
            if(curspan == srcspan) {
              ss = symSpans[srIdxStart];
              se = ss + curspan;
              checkTarget = true; // sample something
              break;
            }
          } while(curspan < srcspan && (++sidx < (int)prtRule.s.size()));
          if(checkTarget) {
            int maxTrgSpan = std::min((int)prtRule.t.size(), td->maxspan);
            for(int trgspan = maxTrgSpan; trgspan > 0; --trgspan) {
              bool allNewTrg(true);
              for(int tIdxStart=0; tIdxStart + trgspan <= (int)prtRule.t.size(); ++tIdxStart) {
                rule.t.clear();
                for(int tidx=tIdxStart; tidx-tIdxStart < trgspan; ++tidx) {
                  rule.t.push_back(prtRule.t[tidx]);
                }
                //if(debug) cout << "Checking validity: " << rule << endl;
                //++noSamplesBuilt_;
                if(validRule(rule, prtRule, &allNewSrc, &allNewTrg, 
                  srcNTCache, trgNTCache, &lftRuleCache)) {
                  //++noSamplesChecked_;
                  if(debug) { 
                    //cout << "rule is valid" << endl;
                    //cout << "Current parent rule : " << prtRule << endl;
                    cout << "Checking valid rule : " << rule << endl; 
                  }
                  if(alreadyOn(rule)) {
                    newPrtSpan = sampleChildOff(rule.s[0], srIdxStart, tIdxStart, prtP0,
                      prtProb, trit, trees_[t], td->anneal, debug);
                    if(newPrtSpan) {
                      --srIdxStart; //keep src start index the same
                      //++timesSampled[transform(trit->rule)];
                    }
                  }
                  else {
                    newPrtSpan = sampleChildOn(rule, srIdxStart, tIdxStart, prtP0, 
                      prtProb, trit, trees_[t], td->anneal, debug, ss, se);
                    /*if(newPrtSpan) {
                      ++timesSampled[transform(rule)];
                    }*/
                  }
                  if(newPrtSpan) { // get new parent's stats
                    prtRule = trit->rule;
                    prtTrsfm = transform(prtRule);
                    prtP0 = get_p0(&prtRule, debug);
                    prtProb = restr_->log_prob(prtTrsfm, prtP0) * td->anneal;
                    //tIdxStart =  (int)prtRule.t.size(); // kill inner target loop
                    trgspan = -1; // kill outer target loop
                    symSpans = precomputeSpans(trit, debug);
                    break;
                  }
                } // end validRule
              } // end tIdxStart
            } // end trgspan
          }  // end checkTarget
        }  // end srIdxStart 
      }  // end maxSrcSpan
      increment(prtTrsfm, exp(prtP0)); // add back current rule into restr
      ++trit;
    } // end sentence 
    if(debug) printTree(t);
  } // end batch of sentences 
  pthread_exit((void*) threadArgs);
}
bool sampleChildOff(const int nodeLbl, const int sIdx, const int tIdx, 
  const double prtSpanP0, const double  prtSpanPrb, 
  tree<node_t>::iterator node, tree<node_t>* fullTree, const float annealT,
  const bool debug) {
  bool newSample(false); 
  hieroRule prtSpan = node->rule;
  // find child span to merge
  tree<node_t>::iterator child = std::find(node.begin(), node.end(), nodeLbl);
  assert(child != node.end());
  const hieroRule& ocr = child->rule; // old child rule
  // have child span so merge with parent span 
  hieroRule npr = prtSpan;
  npr.s.erase(npr.s.begin() + sIdx);
  npr.s.insert(npr.s.begin() + sIdx, ocr.s.begin(),
    ocr.s.end());
  npr.t.erase(npr.t.begin() + tIdx);
  npr.t.insert(npr.t.begin() + tIdx, ocr.t.begin(),
    ocr.t.end());
  decrement(transform(ocr));
  // get probs for potential rules 
  double nprP0 = get_p0(&npr, debug);
  double ocrP0 = get_p0(&ocr, debug);
  double nprPrb = restr_->log_prob(transform(npr), nprP0) * annealT; // anneal 
  double ocrPrb = restr_->log_prob(transform(ocr), ocrP0);
  double currRulesProb = (prtSpanPrb + ocrPrb) * annealT; // combine and anneal 
  double log_normalizer = Log<double>::add(currRulesProb, nprPrb);
  double rSample = log_normalizer + log(dRand());
  if(debug) {
    cout << "sampling off...\n";
    cout << "currparentRule : " << prtSpan << " [" << prtSpanP0 << ", " 
      << restr_->count(transform(prtSpan)) << ", " << prtSpanPrb << "]" << endl;
    cout << "currchildRule : " << ocr << " [" << ocrP0 << ", " 
      << restr_->count(transform(ocr)) << ", " << ocrPrb << "]" << endl;
    cout << "current Rules Prob (combined and annealed) : " << currRulesProb << endl;
    cout << "newParentRule : " << npr << " [" << nprP0 << ", " << 
      restr_->count(transform(npr)) << ", " << nprPrb << "]" << endl;
    cout << "normalizer : " << log_normalizer << "\tsampled point : " << rSample << endl;
  }
  if(rSample > currRulesProb) {
    if(debug) {
      cout << " sample accepted\n\n";
    }
    node->rule = npr; // restructure tree
    fullTree->reparent(node, child); // move any children of 'child' to parent 
    fullTree->erase(child); 
    newSample = true;
  }
  else {
    if(debug) {
      cout << " sample rejected\n\n";
    }
    increment(transform(ocr), exp(ocrP0)); // put child rule back into restr
  }
  return newSample;
}
bool sampleChildOn(const hieroRule& ncr, const int sIdx, const int tIdx, 
  const double prtSpanP0, const double prtSpanPrb,
  tree<node_t>::iterator node, tree<node_t>* fullTree, const float annealT,
  const bool debug, const int ss, const int se) {
  bool newSample(false);
  hieroRule prtSpan = node->rule;
  int& lbl = fullTree->begin()->NTIndex; // get current label from full tree
  hieroRule npr = newParentRule1(prtSpan, ncr, sIdx, tIdx, lbl-1); // create new parent rule 
  // get probs for potential rules 
  double nprP0 = get_p0(&npr, debug);
  double ncrP0 = get_p0(&ncr,debug);
  double nprPrb = restr_->log_prob(transform(npr), nprP0);
  double ncrPrb = restr_->log_prob(transform(ncr), ncrP0); // assume independence for now
  double newRulesProb = (nprPrb + ncrPrb) * annealT; // normalize and anneal 
  double log_normalizer = Log<double>::add(prtSpanPrb, newRulesProb);
  double rSample = log_normalizer + log(dRand());
  if(debug) {
    cout << "sampling on...\n";
    cout << "oldparentRule : " << prtSpan << " [" << prtSpanP0 << ", " 
      << restr_->count(transform(prtSpan)) << ", " << prtSpanPrb << "]" << endl;
    cout << "newParentRule : " << npr << " [" << nprP0 << ", " 
      << restr_->count(transform(npr)) << ", " << nprPrb << "]" << endl;
    cout << "newChildRule : " << ncr << " [" << ncrP0 << ", " 
      << restr_->count(transform(ncr)) << ", " << ncrPrb << "]" << endl;
    cout << "newRulesProb (combined and annealed) : " << newRulesProb << endl;
    cout << "normalizer : " << log_normalizer << "\tsampled point : " << rSample << endl;
  }
  if(rSample > prtSpanPrb) {
    if(debug) {
      cout << " sample accepted\n\n";
    }
    updateTree(node, fullTree, npr, ncr, --lbl, ss, se);
    increment(transform(ncr), exp(ncrP0)); // add new child rule to restaurant
    newSample = true;
  }
  else {
    if(debug) {
      cout << " sample rejected\n\n";
    }
  }
  return newSample;
}
void printTreeNode(tree<node_t>::iterator node) {
  cout << "X(" << node->label << ")-> " << node->rule;
  cout << "[" << node->sstr << "-" << node->send << " | ?-? ]";
  //out << " | " << X.tstr << "-" << X.tend << "]";
  if(node.node->parent)
    cout << "  (parent = " << node.node->parent->data.label << ")";
  cout << endl;
}
std::ostream& operator<<(std::ostream& out, const hieroRule& X) {
  out << " < ";
  iterate(X.s, sit) {
    if(isTerm(*sit))
      out << srcVcb_->getWord(*sit) << " ";
    else
      out << *sit << " ";
  }
  out << "||| "; 
  iterate(X.t, tit) {
    if(isTerm(*tit))
      out << trgVcb_->getWord(*tit) << " ";
    else
      out << *tit << " ";
  }
  out << " > ";
  return out;
}
void getTermIdxs(const tree<node_t>::iterator iter, map<int, set<int> >& lblidx, 
    int* currIdx, const tree<node_t>& tr, bool target) {
  if(iter == tr.end())
    return;
  const vector<int>& rule = target ? iter->rule.t : iter->rule.s; // source or target??
  iterate(rule, r) {
    if(isTerm(*r)) { //  found a terminal symbol
      lblidx[iter->label].insert(*currIdx);
      ++(*currIdx);
    }
    else { // non-term 
      tree<node_t>::iterator child = std::find(iter.begin(), iter.end(), *r);  // get child 
      assert(child != iter.end());
      getTermIdxs(child, lblidx, currIdx, tr, target);
    }
  } 
}
void saveTreeState() {
  std::stringstream fname;
  fname << outputDir_.str() << "/trees.state.gz";
  FileHandler fout(fname.str(), std::ios::out, false); 
  for(size_t t=0; t< trees_.size(); ++t) {
    const tree<node_t>& tr = *trees_[t];  // get tree for snt
    tree<node_t>::iterator itr = tr.begin();
    fout << "sentence " << t << endl;
    for(; itr != tr.end(); ++itr) {
      if(itr != tr.begin())
        fout << itr.node->parent->data.label;
      fout << " ||| ";
      const hieroRule r = itr->rule;
      iterate(r.s, s) fout << *s << " ";
      fout << " ||| ";
      iterate(r.t, t) fout << *t << " ";
      fout << " ||| ";
      fout << itr->label << " " << itr->NTIndex << " " 
        << itr->sstr << " " << itr->send << endl;
    }
  }
  fout.close();
}
void readTreeState() {
  cerr << "Loading former sentence derivations...\n";
  FileHandler fin(params_->getParam("init-aligns"), std::ios::in);
  string line;
  vector<string> v, v2;
  tree<node_t>* tr; 
  while(getline(fin, line)) {
    Utils::trim(line);
    if(line.find("sentence ") == 0) { // start of new sentence
      if(line != "sentence 0")
        trees_.push_back(tr);
      tr = 0;
    }
    else {
      Utils::splitToStrMD(line, v, "|||");
      hieroRule r;
      Utils::splitToStr(v[1], v2, " ");
      iterate(v2, vit)
        r.s.push_back(atoi(vit->c_str()));
      Utils::splitToStr(v[2], v2, " ");
      iterate(v2, vit)
        r.t.push_back(atoi(vit->c_str()));
      // add rule to restr
      increment(transform(r), exp(get_p0(&r, false)));
      node_t node;
      node.rule = r;
      Utils::splitToStr(v[3], v2, " ");
      node.label = atoi(v2[0].c_str());
      node.NTIndex = atoi(v2[1].c_str());
      node.sstr = atoi(v2[2].c_str());
      node.send = atoi(v2[3].c_str());
      // add rule to tree
      if(tr == 0) {
        tr = new tree<node_t>;
        tr->set_head(node);
      }
      else {
        assert(tr);
        int prtLbl = atoi(v[0].c_str());
        tree<node_t>::iterator parent;
        parent = std::find(tr->begin(), tr->end(), prtLbl);  // get child 
        assert(parent != tr->end()); 
        tr->append_child(parent, node); 
      }
    }
  }
  trees_.push_back(tr); // /last sentence
  cerr << "Loaded " << trees_.size() << " sentences\n";
  //assert(trees_.size() == src_sents.size());
}
void printAlignments(int iter) {
  cerr << "Saving current alignments for iteration " << iter << "...\n";
  std::stringstream fname;
  fname << outputDir_.str() << "/alignment." << iter << ".grow-diag-final-and";
  //FileHandler fout("___stdout___", std::ios::out, false); 
  FileHandler fout(fname.str(), std::ios::out, false); 
  for(size_t snt=0; snt < trees_.size(); ++snt) {
    map<int, set<int> > trgIndexes, srcIndexes; // index->label 
    const tree<node_t>& tr = *trees_[snt];  // get tree for snt
    int currIdx=0;
    getTermIdxs(tr.begin(), trgIndexes, &currIdx, tr, true);
    currIdx=0;
    getTermIdxs(tr.begin(), srcIndexes, &currIdx, tr, false);
    iterate(srcIndexes, sit) { // for each label 
      map<int, set<int> >::iterator it = trgIndexes.find(sit->first); 
      if(it == trgIndexes.end()) continue;
      iterate(sit->second, ssit) {
        iterate(it->second, ttit) {
          fout << *ssit << "-" << *ttit << " ";
        }
      }
    }
    fout << endl;
  }
  fout.close();
  cerr << "Alignments saved to " << fname.str() << endl;
}
void runSampling() {
  const int numThreads = atoi(params_->getParam("threads").c_str());
  const int batchSize = trees_.size() / numThreads;
  time_t iterFinish, iterStart;
  double totSecs(0);
  for(int iter=1; iter <= totIter_; ++iter) {
    vector<threadData_t> threadData(numThreads);
    pthread_t threads[numThreads];
    currIter_ = iter;
    cout << "\nSampling Iteration " << iter << "\n";
    time(&iterStart);
    for(int i=0; i < numThreads; ++i) { // create threaded jobs
      threadData[i].sntBegin = batchSize * i;
      threadData[i].sntEnd = (i == numThreads - 1) ? trees_.size() 
        : (batchSize * i) + batchSize;
      threadData[i].anneal = annealT();
      threadData[i].maxspan = noInit_ ? maxSpan() : MAX_SENT_LEN;
      threadData[i].skipDecr = noInit_ && (iter==1);
      //if(iter % 12 == 1 || (iter % 12 == 2))
        //pthread_create(&threads[i], &attr, sampleMN, (void*)&threadData[i]);
      //else
      pthread_create(&threads[i], &attr, sample, (void*)&threadData[i]);
    }
    for(int i=0; i < numThreads; ++i) { // wait on the other threads 
      pthread_join(threads[i], NULL);
    }
    time(&iterFinish);
    double iterSecs = difftime(iterFinish, iterStart);
    cerr << "Total seconds for iteration " << iter << " : " << iterSecs;
    cerr << " (" << (iterSecs / float(trees_.size())) << " per sentence.)\n";
    totSecs += iterSecs;
    if(iter % 10 == 0) {
      printAlignments(iter);
      restr_->debug_info(cout);
      cerr << "Resampling hyperparameters..\n";
      restr_->resample_prior(mt_genrand_res53);
      cerr << "a=" << restr_->a() << "\tb=" << restr_->b() << endl;
      resampleTermPenalty();
      cerr << "Data LL = " << dataLLAndOtherStuff(iter) << endl;
      saveTreeState();
    }
  } // end iterations
  printAlignments(currIter_);
  cerr << "sampled " << trees_.size() << " sentences for " << totIter_ << " iterations\n";
  cerr << "total seconds taken : " << totSecs << " seconds\n";
  cerr << "time per sentence : " << (totSecs / float(trees_.size() * totIter_)) << " seconds" << endl;
  /*cerr << "total number of samples built: " << noSamplesBuilt_ << endl;
  cerr << "total number of samples checked : " << noSamplesChecked_ << endl;
  std::stringstream fname;
  fname << outputDir_.str() << "/samples";
  FileHandler f(fname.str(), std::ios::out, false);
  iterate(timesSampled, it)
    if(it->second > 1)
      f << it->first << "\t" << it->second << endl;
  f.close();*/
}
int main(int argc, char** argv) {
  time_t curtime; 
  time(&curtime);
  srand48(curtime);
  params_ = new Parameters(argc, argv, paramdefs, NumOfParams(paramdefs));
  if(params_->getParam("outputdir").empty())
    outputDir_ << "pyp/outputfiles/" << curtime;
  else
    outputDir_ << params_->getParam("outputdir");
  if(!Utils::fileExists(outputDir_.str())) {
    string cmd = "mkdir -p " + outputDir_.str();
    cerr << "Making output directory " << outputDir_.str() << endl;
    system(cmd.c_str());
  }
  snt2debug_ = atoi(params_->getParam("snt2Debug").c_str());
  /* Initialize and set thread detached attribute */
  cout.precision(10);
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  pthread_mutex_init(&mutexsum, NULL);
  // set class variables
  noInit_ = params_->getParam("no-init") ==
    Parameters::kTrueValue ? true: false;
  tp_ = atof(params_->getParam("terminal-penalty").c_str());
  currIter_=0;
  totIter_ = atoi(params_->getParam("iterations").c_str());
  noSamplesChecked_ = 0;
  noSamplesBuilt_ = 0;
  loadData();
  initPYP();
  runSampling();
  printArityHist();
  freeStuff();
  pthread_attr_destroy(&attr);
  pthread_mutex_destroy(&mutexsum);
  pthread_exit(NULL);
  return 1;
}
