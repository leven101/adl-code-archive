#ifndef INC_VOCAB_H
#define INC_VOCAB_H

#include <map>
#include <string>
#include "types.h"
#include "file.h"
#include "utils.h"

//namespace randlm {
  // Vocab maps between strings and uint32 ids.

  class Vocab {
  public:
    static const wordID_t kOOVWordID = 0;   // out of vocabulary word id
    static const wordID_t kBOSWordID = 1;    
    static const word_t kBOS;  // beginning of sentence marker
    static const word_t kEOS;  // end of sentence marker
    static const word_t kOOVWord;  // <unk>
    Vocab(bool sntMarkers = true):closed_(false) {
      if(sntMarkers) {
        getWordID(kBOS);  // added in case not observed in corpus
        getWordID(kEOS);
      }
    }  
    // if no file then must allow new words
    // specify whether more words can be added via 'closed'
    // assume that if a vocab is loaded from file then it should be closed.
    Vocab(const std::string & vocab_path, bool closed = true) {
      assert(load(vocab_path, closed));  
    }
    Vocab(FileHandler* fin, bool closed = true) {
      assert(load(fin, closed));  
    }
    ~Vocab() {}
    wordID_t getWordID(const word_t & word);
    word_t getWord(wordID_t id);
    bool inVocab(wordID_t id);
    bool inVocab(const word_t & word);
    uint32_t size() { return words2ids_.size(); }
    void makeClosed() { closed_ = true; }
    void makeOpen() { closed_ = false; }
    bool isClosed() { return closed_; }
    bool save(const std::string & vocab_path);
    bool save(FileHandler* fout);
    bool load(const std::string & vocab_path, bool closed = true);
    bool load(FileHandler* fin, bool closed = true);
    void printVocab();
    std::map<word_t, wordID_t>::const_iterator vocabStart() {
      return words2ids_.begin();
    }
    std::map<word_t, wordID_t>::const_iterator vocabEnd() {
      return words2ids_.end();
    }
  private:
    std::map<word_t, wordID_t> words2ids_;  // map from strings to word ids
    std::map<wordID_t, word_t> ids2words_;  // map from ids to strings
    bool closed_;  // can more words be added
  };
//}

#endif // INC_VOCAB_H
