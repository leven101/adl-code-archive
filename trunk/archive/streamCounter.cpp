#include "streamCounter.h"
#include "utils.h"


StreamCounter::StreamCounter() : total_tokens_(0), period_(""),
 uni_idx_(0), bi_idx_(1), tri_idx_(2) {
  for(int i=0; i < NGRAM_ORDER; i++) {
    ngrm_cntrs_[i] = new date_cnt_t;
    longstring_ ="";
  }
}

StreamCounter::~StreamCounter() {
  clearCntrs();
  for(int i=0; i < NGRAM_ORDER; i++)
    delete ngrm_cntrs_[i];
}

void StreamCounter::extractRawData(std::istringstream& stream) {
  word_t sdate;
  stream >> sdate;
  if( sdate != period_ ) {
    if(! longstring_.empty() ) {
      printRaw();
      longstring_.clear();
    }
    period_ = sdate;
  }
  std::vector<string> tmp;
  Utils::tokenizeToStr(stream.str(), tmp, "\n");
  assert(tmp.size() >= 2);
  //prints everything minus date and 
  // "Reuters Limited 1997" etc
  for(int i = 1; i < tmp.size() - 1; i++)
    longstring_ += tmp[i] + "\n";
}

void StreamCounter::printRaw() {
  string fname = "../output/" + period_ + ".rawtext.gz";
  FileHandler fout(fname, std::ios::out, false);
  fout << longstring_;
}

void StreamCounter::countStreamItems(std::istringstream& stream) {
  word_t sdate, tokens[NGRAM_ORDER], token, ngram;
  //first item is date
  stream >> sdate;
  if( sdate != period_ ) {
    if( total_tokens_ > 0 ) { // not the very beginning
      // print out contents of ngrm_cntrs
      print(period_); 
      clearCntrs(); // clear ngrm_cntrs
    }
    period_ = sdate; // assign new period
  }
  // initialize token array
  for( int n = 0; n < NGRAM_ORDER - 1; n++ )
    stream >> tokens[n]; 
  size_t tokenCnt = NGRAM_ORDER - 1;
  // while stream exists build ngrams from order 
  // 1 to NGRAM_ORDER - 1.
  while( stream >> token ) {
    //remove punctuation
    Utils::trim(token, " \",+/*:-()'<>;#"); 
    if( token.empty() ) continue;
    Utils::strToLowercase(token);
    tokens[tokenCnt % NGRAM_ORDER] = token;
    ngram = "";
    //build and record all order ngrams
    for( int n = NGRAM_ORDER - 1; n >= 0; n-- ) {
      int idx = (tokenCnt - n) % NGRAM_ORDER;
      ngram += " " + tokens[idx];
      //add to appropriate ngram container by date
      if( n == 2 ) {
        // ???Utils::trim(token, " .?!"); //remove eos punctuation
        if( token.empty() ) continue;
        (*ngrm_cntrs_[uni_idx_])[sdate][token]++;
      }
      else if( n == 1 ) //add to bigram
        (*ngrm_cntrs_[bi_idx_])[sdate][ngram]++;
      else if( n == 0 ) //add to trigram
        (*ngrm_cntrs_[tri_idx_])[sdate][ngram]++;
    }
    tokenCnt++; 
  }
  total_tokens_ += tokenCnt;
}

void StreamCounter::print(string prefix)
{
  cerr << "Total items = " << total_tokens_ << endl;
  string path = "../output/" +  
                (prefix.empty() ? period_ : prefix)
                + ".";
  printCounts("unigrams", path + "unigrams.gz",
                *ngrm_cntrs_[uni_idx_]);
  printCounts("bigrams", path + "bigrams.gz",
                *ngrm_cntrs_[bi_idx_]);
  printCounts("trigrams", path + "trigrams.gz",
                *ngrm_cntrs_[tri_idx_]);
}

void StreamCounter::printCounts(const string& what, const string& filename,
                                date_cnt_t& cntr) {
  cerr << "Printing " << what << " distribution to ";
  cerr << filename << ".\n";
  {
  FileHandler out(filename, std::ios::out, false);
  iterate(cntr, i) {
    out << i->first << endl;
    iterate(i->second, p) {
      out << p->second;
      out << "\t" << p->first << "\n";
    }
  }
  out.flush();
  }
  sort(filename);
}

void StreamCounter::sort(const string& filename) {
  // issues a cmd to shell to sort file
  string tmp = filename;
  tmp.replace(filename.find(".gz"), 3, ".sorted.gz"); 
  string cmd = "zcat " + filename + " | sort -n -r | gzip -f > ";
  cmd += tmp + " ;rm -f " + filename;
  cerr << cmd << endl;
  int r = system(cmd.c_str());
  if( r != 0 ) perror("sort: ");
}

void StreamCounter::clearCntrs() {
  for( int i = NGRAM_ORDER - 1; i > -1; i-- ) {
    date_cnt_t& tmpCntr = *ngrm_cntrs_[i];
    iterate(tmpCntr, itr) {
      itr->second.clear();
    }
    tmpCntr.clear();
  }
}

