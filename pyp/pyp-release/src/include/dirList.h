#ifndef DIRLIST_H
#define DIRLIST_H

#include <vector>
#include <string>
#include <cstdio>
#include <cstdlib>

class DirList {
  private:
    typedef void (*FUNC_PTR)(const std::string &fileName);
    const std::string fileTypes_;
    FUNC_PTR funcPtr_;
    std::vector<std::string> fileList_;
  public:
    const std::vector<std::string> &cFileList_;

    DirList(const std::string path, const std::string &types, 
            FUNC_PTR fptr = NULL, bool searchSubDirs = false)
    :fileTypes_(types), funcPtr_(fptr), cFileList_(fileList_) {
      findInDir(path, searchSubDirs);
    }
    DirList(const std::string path )
    :fileTypes_(""), funcPtr_(0), cFileList_(fileList_) {
      findDirs(path);
    }
    virtual ~DirList() {
      fileList_.clear();
      if( funcPtr_ ) funcPtr_ = NULL;
    }
    virtual void execFunc(FUNC_PTR fptr) {
      const int count = fileList_.size();
      for( int i = 0; i < count; ++i )
	fptr(fileList_[i]);
     }
        
    template<typename T, typename F>
    void execClassFunc(T &srcObj, F func) {
      const int count = fileList_.size();
      for (int i = 0;i < count;++i)
      {
        (srcObj.*func)(fileList_[i]);
      }
    }
  protected:
    virtual void findDirs(const std::string &path) {
      std::string cmd = "find " + path + " -type d -print"; 
      executeCmd(cmd);
    }
    virtual void findInDir(const std::string &path, bool searchSubDirs) {
      std::string cmd = (searchSubDirs)?(std::string("find ") + path + 
                         " -name '" + fileTypes_ + "' -print")
                         :(std::string("ls -1 ") + path + "/" + fileTypes_ 
                           /*+ " 2> /dev/null"*/);
      //std::cerr << cmd << std::endl;
      executeCmd(cmd);
    }
  private:
    virtual void executeCmd(const std::string& cmd) {
      FILE *fp = popen(cmd.c_str(), "r");
      if( !fp ) return;
      char psBuffer[1024];
      while( !feof(fp) ) {
        if( fgets(psBuffer, sizeof(psBuffer), fp ) != NULL ) {
          if( psBuffer[strlen(psBuffer)-1] == '\n' )
            psBuffer[strlen(psBuffer)-1] = 0;
          fileList_.push_back(psBuffer);
          if( funcPtr_ ) 
            funcPtr_(psBuffer);
        }
      }
      pclose(fp);
    }
};
#endif // DIRLIST_H
