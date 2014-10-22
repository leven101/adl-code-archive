#include "multiNTs.h"
#include "binary_sampler.h"
#include "block_sampler.h"
#include "multinomial_sampler.h"
#include "itg2trees.h"

MultiNT::p0_t MultiNT::get_p0(const hieroRule* hr, const bool debug) {
  static const double usntp = log(1.0 / pow(numSrcWrdCls_, 2)); // uniform src NT prob
  double arity, baseP0 = base_p0(hr, debug, arity);
  p0_t p0;
  if(hierch_) { 
    p0.first = base_restr_->log_prob(transform(*hr, false), baseP0);
  }
  else p0.first = baseP0;
  double pnts = usntp * arity; // probability of individual NTs
  p0.second = p0.first + pnts;
  if(debug) {
    cout << "\tpnts: " << pnts << "\tp0: " << p0.second << endl;
  }
  return p0;
}
double MultiNT::base_p0(const hieroRule* hr, const bool debug, double& arity) {
  static const double usp = log(1.0 / (double)srcVcb_->size()); // uniform src prob
  static const double utp = log(1.0 / (double)trgVcb_->size()); // uniform trg prob
  const double sl = hr->s.size();
  const double tl = hr->t.size();
  arity = hr->arity(); 
  const double nst = sl - arity;
  const double ntt = tl - arity;
  // probability of source length
  static const float pm = atof(params_->getParam("poisson-mean").c_str());
  double psl = log_poissDist(sl, pm);
  // probability of source side configuration
  double tpp = pow(termPenalty(), sl); // term penalty param
  double ssc = (nst * log(tpp)) + (arity * log(1-tpp));
  double ptl = log_poissDist(ntt, nst+0.01); 
  // probability of target side configuration
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

  //double p0 = psl + ptl + ssc + tsc + joint_ts + pnts;
  double p0 = psl + ptl + ssc + tsc + joint_ts;
  if(debug) {
    cout << "prior for rule " << *hr << " = psl: " << psl << "\tptl: " << ptl
      << "\tssc: " << ssc << "\ttsc: " << tsc << "\tjoint_ts: " << joint_ts;
      //<< "\tpnts: " << pnts << "\tp0: " << p0 << endl;
  }
  return p0; 
}
void MultiNT::fullSentenceInit() {
  for(size_t i=0; i < src_sents.size(); ++i) {
    bool debug = (int)i == snt2debug_;
    hieroRule sntRule;
    sntRule.NT = 0;
    sntRule.s.assign(src_sents[i].begin(), src_sents[i].end());
    sntRule.t.assign(trg_sents[i].begin(), trg_sents[i].end());
    tree<node_t>* tr = new tree<node_t>;
    node_t topNode;
    topNode.rule = sntRule;
    topNode.label = 1; // set label;
    topNode.nodeIndex = 0; // total NTs in tree;
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
void MultiNT::loadLongSnts() {
  int longSnts(0);
  const string dir = "/data/taipan/ablev/data/ur-en/";
  FileHandler fhmm(dir + "hmm.justLongSnts", std::ios::in);
  FileHandler fsrc(dir + "urdu.justLongSnts", std::ios::in);
  FileHandler ftrg(dir + "english.justLongSnts", std::ios::in);
  /* Load the sentences */
  cerr << "Loading long source sentences...\n";
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
    ++longSnts;
  }
  fsrc.close();
  srcVcb_->makeClosed();
  cerr << "Loading long target sentences...\n";
  while(getline(ftrg, line)) { 
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
  cerr << "Loaded " << longSnts << " more sentence pairs.\n";
  /* Load the alignments */
  vector<int> v;
  for(size_t i=src_sents.size() - longSnts; i < src_sents.size(); ++i) {
    bool debug = (int)i == snt2debug_;
    vector<int> alGrid(trg_sents[i].size(), -1);
    vector<bool> srcIdxUsed(src_sents[i].size(), false);
    getline(fhmm, line);
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
      const hieroRule& r = trit->rule;
      increment(transform(r, trit), get_p0(&r, debug)); // add rule to rest
      ++trit;
    }
    trees_.push_back(tr);
    if(debug) printTree(i);
  } // end corpus
  cerr << "Finished loading long sentences." << endl;
  fhmm.close();
}
void MultiNT::initPYP() {
  a_ = atof(params_->getParam("discount").c_str());
  b_ = atof(params_->getParam("strength").c_str());
  restr_ = new vector<PYP<hieroRule, hieroRule> >;
  if(singleNT_) {
    restr_->push_back(PYP<hieroRule, hieroRule>(a_,b_));
  }
  if(hierch_) {
    base_restr_ = new PYP<hieroRule, hieroRule>(a_, b_);
  }
  const bool restart = params_->getParam("restart") ==
      Parameters::kTrueValue ? true: false;
  if(restart) {
    int maxNT = readTreeState();
    assert(trees_.size() == src_sents.size());
    // create old restaurants
    while(maxNT >= (int)restr_->size()) { // ensure PYP for every NT
      restr_->push_back(PYP<hieroRule, hieroRule>(a_,b_));
    }
    int maxarity=0;
    // add all rules to restr
    for(size_t snt=0; snt < trees_.size(); ++snt) {
      const tree<node_t>& tr = *trees_[snt];  // get tree for snt
      tree<node_t>::iterator trit = tr.begin();
      while(trit != tr.end()) {
        const hieroRule& r = trit->rule;
        increment(transform(r, trit), get_p0(&r, false));
        if(r.arity() > maxarity) maxarity = r.arity();
        ++trit;
      }
      if((int)snt == snt2debug_) printTree(snt);
    }
    cerr << "Maximum rule arity: " << maxarity << endl;
    loadLongSnts();  // load long HMM sentences here
    resampleHyperParams(); // get good hyper-parameters 
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
        const hieroRule& r = trit->rule;
        increment(transform(r, trit), get_p0(&r, debug)); // add rule to rest
        ++trit;
      }
      trees_.push_back(tr);
      if(debug) printTree(i);
    } // end corpus
    m1alg.close();
  }
  printAlignments(-1);
  printRules(-1);
  dataLLAndOtherStuff(-1);
}
tree<node_t>* MultiNT::buildXRuleGrid(const vector<int>& grid, const int snt) {
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
            X.NT = NTSymbol(src[tmp[0]], src[tmp.back()]); // add NT symbol to rule
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
        X.NT = NTSymbol(src[grid[i]], src[grid[i]]);
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
    else { // terminal
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
  topX.NT = NTSymbol(src[0], src.back());
  node_t topNode;
  topNode.rule = topX;
  topNode.label = 1; // set label
  topNode.nodeIndex = ntIdx+1; // track nodes in tree
  topNode.sstr = 0;
  topNode.send = src.size();
  //topNode.tstr = 0;
  //topNode.tend = trg.size();
  top.node->data = topNode; // set top node to sentence rule 
  return tr;
}
void MultiNT::loadM1Params() {
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
void MultiNT::loadNonterminals() {
  numSrcWrdCls_ = 0;
  if(singleNT_) {
    numSrcWrdCls_ = 1;
    cerr << "Using restricted SCFG with single nonterminal X\n" << endl;
    return;
  }
  FileHandler sfin(params_->getParam("src-nts"), std::ios::in);
  string line;
  vector<string> vline; 
  while(getline(sfin, line)) {
    Utils::splitToStr(line, vline, "\t");
    wordID_t sid = srcVcb_->getWordID(vline[0]); 
    assert(sid != Vocab::kOOVWordID);
    int cls = atoi((vline[1].c_str()));
    srcWrdCls_[sid] = cls;
    if(cls > numSrcWrdCls_) numSrcWrdCls_ = cls;
  }
  --numSrcWrdCls_; // always 1 more than specified from mkcls
  sfin.close();
  /*FileHandler tfin(params_->getParam("trg-nts"), std::ios::in);
  while(getline(tfin, line)) {
    Utils::splitToStr(line, vline, "\t");
    wordID_t tid = trgVcb_->getWordID(vline[0]); 
    assert(tid != Vocab::kOOVWordID);
    int cls = atoi((vline[1].c_str()));
    trgWrdCls_[tid] = cls;
  }
  tfin.close();*/
  cerr << "Size of srcWrdCls map: " << srcWrdCls_.size() << endl;
  //cerr << "Number of target classes: " << trgWrdCls_.size() << endl;
}
void MultiNT::loadSentences() {
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
void MultiNT::loadData() {
  srcVcb_ = new Vocab(false);
  trgVcb_ = new Vocab(false);
  nts_ = new Vocab(false);
  loadM1Params();
  loadSentences();
  loadNonterminals();
  cerr << "Source vocabularly size = " << srcVcb_->size() << endl;
  cerr << "Target vocabularly size = " << trgVcb_->size() << endl;
}
void MultiNT::printAlignments(int iter) {
  // for each sentence pair get map of index->node label for source and 
  // target sides recursively and then print matching links 
  cerr << "Saving current alignments for iteration " << iter << "...\n";
  std::stringstream fname;
  fname << outputDir_.str() << "/alignment." << iter << ".grow-diag-final-and";
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
void MultiNT::getTermIdxs(const tree<node_t>::iterator iter, map<int, set<int> >& lblidx, 
    int* currIdx, const tree<node_t>& tr, bool target) {
  if(iter == tr.end())
    return;
  vector<int> rule = target ? iter->rule.t : iter->rule.s; // source or target??
  iterate(rule, r) {
    if(isTerm(*r)) { //  found a terminal symbol
      lblidx[iter->label].insert(*currIdx);
      ++(*currIdx);
    }
    else { // non-term 
      tree<node_t>::iterator child = std::find(iter.begin(), iter.end(), *r);  // get child 
      getTermIdxs(child, lblidx, currIdx, tr, target);
    }
  } 
}
void* sample(void *threadArgs) {
  BinarySampler smp(threadArgs);
  pthread_exit((void*) threadArgs);
}
void* blocksample(void *threadArgs) {
  BlockSampler bs(threadArgs);
  pthread_exit((void*) threadArgs);
}
void* sampleMN(void *threadArgs) {
  MultinomialSampler mns(threadArgs); 
  pthread_exit((void*) threadArgs);
}
void MultiNT::runSampling() {
  const int numThreads = atoi(params_->getParam("threads").c_str());
  const int maxThreads = atoi(params_->getParam("threads").c_str());
  const int bs_every = atoi(params_->getParam("bs-every").c_str());
  const int batchSize = trees_.size() / numThreads;
  const bool restart = params_->getParam("restart") ==
      Parameters::kTrueValue ? true: false;
  time_t iterFinish, iterStart;
  double totSecs(0);
  for(int iter=1; iter <= totIter_; ++iter) {
    vector<threadData_t> threadData(numThreads);
    pthread_t threads[numThreads];
    currIter_ = iter;
    cout << "\nSampling Iteration " << iter << "\n";
    time(&iterStart);
    bool doblock = ((blockSample_) && (iter % bs_every == 0)) && iter > 10;
    GrammarCache* cache = doblock ? new GrammarCache(this) : NULL;
    for(int i=0; i < numThreads; ++i) { // create threaded jobs
      threadData[i].sntBegin = batchSize * i;
      threadData[i].sntEnd = (i == numThreads - 1) ? trees_.size() 
        : (batchSize * i) + batchSize;
      threadData[i].anneal = annealT();
      threadData[i].maxspan = noInit_ ? maxSpan() : MAX_SENT_LEN();
      threadData[i].skipDecr = noInit_ && (iter==1);
      threadData[i].mnt = this;
      threadData[i].cache = cache;
      if(doblock) {
        //BlockSampler bs((void*)&threadData[0]);
        pthread_create(&threads[i], &attr, blocksample, (void*)&threadData[i]);
      }
      else {
        //BinarySampler((void*)&threadData[0]);
        pthread_create(&threads[i], &attr, sample, (void*)&threadData[i]);
      }
    }
    for(int i=0; i < numThreads; ++i) { // wait on the other threads 
      pthread_join(threads[i], NULL);
    }
    if(cache) delete cache;
    time(&iterFinish);
    double iterSecs = difftime(iterFinish, iterStart);
    cerr << "Total seconds for iteration " << iter << " : " << iterSecs;
    cerr << " (" << (iterSecs / float(trees_.size())) << " per sentence.)\n";
    totSecs += iterSecs;
    dataLLAndOtherStuff(iter);
    if(iter % 10 == 0) {
      resampleHyperParams();
      printAlignments(iter);
      printRules(iter);
      saveTreeState();
    }
  } // end iterations
  cerr << "sampled " << trees_.size() << " sentences for " << totIter_ << " iterations\n";
  cerr << "total seconds taken : " << totSecs << " seconds\n";
  cerr << "time per sentence : " << (totSecs / float(trees_.size() * totIter_)) << " seconds" << endl;
}
void MultiNT::saveTreeState() {
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
      fout << r.NT << " ||| ";
      iterate(r.s, s) fout << *s << " ";
      fout << " ||| ";
      iterate(r.t, t) fout << *t << " ";
      fout << " ||| ";
      fout << itr->label << " " << itr->nodeIndex << " " 
        << itr->sstr << " " << itr->send << endl;
    }
  }
  fout.close();
}
int MultiNT::readTreeState() {
  cerr << "Loading former sentence derivations\n";
  FileHandler fin(params_->getParam("init-aligns"), std::ios::in);
  string line;
  vector<string> v, v2;
  tree<node_t>* tr; 
  int maxNT(0);
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
      r.NT = atoi(v[1].c_str());
      maxNT = r.NT > maxNT ? r.NT : maxNT;
      Utils::splitToStr(v[2], v2, " ");
      iterate(v2, vit)
        r.s.push_back(atoi(vit->c_str()));
      Utils::splitToStr(v[3], v2, " ");
      iterate(v2, vit)
        r.t.push_back(atoi(vit->c_str()));
      node_t node;
      node.rule = r;
      Utils::splitToStr(v[4], v2, " ");
      node.label = atoi(v2[0].c_str());
      node.nodeIndex = atoi(v2[1].c_str());
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
  return maxNT;
}
void MultiNT::initialize(int argc, char** argv) {
  params_ = new Parameters(argc, argv, paramdefs, NumOfParams(paramdefs));
  time_t curtime; 
  time(&curtime);
  srand48(curtime);
  mt_init_genrand(curtime);
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
  noInit_ = params_->getParam("no-init") ==
    Parameters::kTrueValue ? true: false;
  tp_ = atof(params_->getParam("terminal-penalty").c_str());
  hierch_ = params_->getParam("hierarchical") ==
    Parameters::kTrueValue ? true: false;
  singleNT_ = params_->getParam("restrict-scfg") ==
    Parameters::kTrueValue ? true: false;
  blockSample_ = params_->getParam("block-sample") ==
    Parameters::kTrueValue ? true: false;
  assert(!(hierch_ && singleNT_));
  currIter_ = 0;
  totIter_ = atoi(params_->getParam("iterations").c_str());
}
std::ostream& operator<<(std::ostream& out, const hieroRule& X) {
  out << X.NT << " --> < ";
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
int main(int argc, char** argv) {
  /* Initialize and set thread detached attribute */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  pthread_mutex_init(&mutexsum, NULL);
  pthread_mutex_init(&mutexsum_NT, NULL);
  MultiNT mnt;
  mnt.initialize(argc, argv);
  mnt.loadData();
  //ITG2Trees itg2trees(&mnt); 
  mnt.initPYP();
  mnt.runSampling(); 
  pthread_attr_destroy(&attr);
  pthread_mutex_destroy(&mutexsum);
  pthread_mutex_destroy(&mutexsum_NT);
  pthread_exit(NULL);
}
