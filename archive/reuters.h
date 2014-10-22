#ifndef REUTERS_STREAM_H
#define REUTERS_STREAM_H

#include <iostream>
#include <sstream>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <cassert>
#include <cstdlib>
#include <pthread.h>
#include <cstdio>
#include "types.h"
#include "dirList.h"
#include "xmlParser.h"
#include "utils.h"

int total_tokens = 0;

class ReutersCorpusStream {
  public:
    ReutersCorpusStream(const string path, const string tmpDir, 
                        const bool threaded);
    ~ReutersCorpusStream();
    template<typename C, typename F>
    bool stream(C &obj, F func);
  private:
    template<typename C, typename F>
      bool singleStream(C &obj, F func);
    template<typename C, typename F>
      bool multiStream(C &obj, F func);
    template<typename C, typename F>
      static void* threadEntry(void* pArgs);
    template<typename C, typename F>
      void threadWorker(C &obj, F func, 
                        const string tmpDir,
                        const std::vector<string> &list);
    bool isDirectory(string path);
    int unpack(const string& dir, const string& zipFile);
    int cleanTmpDir(const string& dir);
    void testIt(std::istringstream& stream);
    const string path_, tmp_dir_;
    const bool threaded_;
    pthread_mutex_t mutex_1;
    size_t fileIndex_;
    template<typename C, typename F>
      struct threadArgs {
        ReutersCorpusStream* rcs_;
        C *obj_;
        F func_;
        string tmpDir_;
        const std::vector<string>* list_;
      };
};

void ReutersCorpusStream::testIt(std::istringstream& stream) {
  word_t word, date;
  //date is first item in each stream 
  stream >> date;
  // for entire stream
  while( stream >> word ) {
    Utils::trim(word);
    // count it
    //unigram_cnts_[word]++;
    total_tokens++;
  }
}

template<typename C, typename F>
bool ReutersCorpusStream::multiStream(C &obj, F func) {
  // get list of zip files in reuters directory
  // only find files that begin with numbers
  DirList reutersDir(path_, "[0-9]*.zip"); 
  assert(reutersDir.cFileList_.size() > 0);
  //setup threads 
  pthread_t threads[THREAD_MAX];
  // a structure to pass the args
  threadArgs<C,F> args_array[THREAD_MAX];
  // spawn some threads
  for( int i=0; i < THREAD_MAX; ++i ) {
    args_array[i].rcs_ = this;
    args_array[i].func_ = func;
    args_array[i].tmpDir_ = tmp_dir_ + Utils::IntToStr(i) + "/";
    args_array[i].list_ = &reutersDir.cFileList_;
    pthread_create(&threads[i], 0, 
                   ReutersCorpusStream::threadEntry<C,F>, 
                   (void*) &args_array[i]);
    sleep(2);  //stagger the thread starts
  }
  for( int i=0; i < THREAD_MAX; ++i ) {
    pthread_join(threads[i], NULL);
  }
  cerr << "Total items = " << total_tokens << std::endl;
  return true;
}

/*static*/
template<typename C, typename F>
void* ReutersCorpusStream::threadEntry(void* pArgs) {
  threadArgs<C,F>* args;
  args = (threadArgs<C,F>*) pArgs;
  // pass process to non-static method of class
  args->rcs_->threadWorker<C,F>(*args->obj_, args->func_,
                                args->tmpDir_, *args->list_);
  return (void*) true; // to shut up compiler
}

template<typename C, typename F>
void ReutersCorpusStream::threadWorker(C &obj, F func, 
     const string tmpDir, const std::vector<string> &list) {
  cerr << "I am threadWorker id: " << (unsigned) pthread_self() << endl;
  int idx = 0;
  XMLParser xml_handler;
  cerr << "add of my xmlHandler is " << &xml_handler << endl;
  cleanTmpDir(tmpDir);
  while( fileIndex_ < /*list.size()*/2 ) {
    // grabs a zip file from index and 
    // increment global counter
    pthread_mutex_lock(&mutex_1);
    idx = fileIndex_++;
    pthread_mutex_unlock(&mutex_1);
    unpack(tmpDir, list[idx]);
    // process each xml file in tmp dir
    DirList tmp(tmpDir, "*.xml");
    assert(tmp.cFileList_.size() > 0);
    iterate(tmp.cFileList_, tmpItr) {
      //std::istringstream packet(xml_handler.streamXMLFile(*tmpItr));
      // pass to stream handler object entry-point 
      //(obj.*func)(packet);
      //testIt(packet);
      //cerr << packet.str() << endl;
    }
    // delete xml files from tmp dir
    cleanTmpDir(tmpDir);
  }
}

ReutersCorpusStream::ReutersCorpusStream(const string path, const string tmpDir, 
                                         const bool threaded) : 
  path_(path), tmp_dir_(tmpDir), threaded_(threaded), fileIndex_(0) {
  // verify the directorys exist
  assert(isDirectory(path_));
  assert(isDirectory(tmp_dir_));
  //if( threaded_ )
  //  pthread_mutex_init(&mutex_1, NULL);
}

ReutersCorpusStream::~ReutersCorpusStream() {
  if( threaded_ )
    pthread_mutex_destroy(&mutex_1);
}

bool ReutersCorpusStream::isDirectory(string path) {
  struct stat s;
  // if the path exists with no error
  if( stat(path.c_str(), &s) == 0 )
    return S_ISDIR(s.st_mode); // return if it's a directory 
  else
    return false;
}

template<typename C, typename F>
bool ReutersCorpusStream::stream(C &obj, F func) {
 // if( threaded_ ) 
 //   return multiStream<C,F>(obj, func);
 // else
  return singleStream<C,F>(obj, func);
}

template<typename C, typename F>
bool ReutersCorpusStream::singleStream(C &obj, F func) {
  // get list of zip files in reuters directory
  // only find files that begin with numbers
  cleanTmpDir(tmp_dir_);
  XMLParser xml_handler;
  DirList reutersDir(path_, "[0-9]*.zip", NULL, false); 
  assert(reutersDir.cFileList_.size() > 0);
  iterate(reutersDir.cFileList_, reuItr) {
    std::cerr << "processing file " << *reuItr << endl;
    // unpack each zipped file into tmp dir
    unpack(tmp_dir_, *reuItr);
    // process each xml file in tmp dir
    DirList tmpDir(tmp_dir_, "*.xml", NULL, false);
    assert(tmpDir.cFileList_.size() > 0);
    iterate(tmpDir.cFileList_, tmpItr) {
      std::istringstream packet(xml_handler.streamXMLFile(*tmpItr));
      // pass to stream handler object 
      (obj.*func)(packet);
    }
    // delete xml files from tmp dir
    cleanTmpDir(tmp_dir_);
   // cerr << "finished with " << *reuItr << endl;
  }
  // process each xml node
  return true;
}

int ReutersCorpusStream::cleanTmpDir(const string& dir) {
  string cmd = "rm -f " + dir + "*.xml";
  //cerr << cmd << endl;
  int r = system(cmd.c_str());
  if( r != 0) perror("cleanTmpDir: ");
  return r;
}

int ReutersCorpusStream::unpack(const string& dir,
                                const string& zipFile) {
  //unzip files silently with overwrites and to specific dir
  string cmd = "unzip -qq -o -d " + dir + " " + zipFile;
  //cerr << cmd << endl;
  int r = system(cmd.c_str());
  if( r != 0) perror("unpack: ");
  return r;
}

#endif //REUTERS_STREAM_H
