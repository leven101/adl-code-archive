#ifndef STREAM_UTILS_H
#define STREAM_UTILS_H

#include <cstdlib>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <sstream>
#include <cctype>
#include <cmath>
#include <cstring>
#include <sys/stat.h>
#include <unistd.h>

using std::vector;
using std::map;
using std::set;

class Utils {
public: 
  static void trim(std::string& str, const std::string dropChars = " \t\n\r") {
    str.erase(str.find_last_not_of(dropChars)+1);
    str.erase(0, str.find_first_not_of(dropChars));
  }
  static void rtrim(std::string& str, const std::string dropChars = " \t\n\r") {
    str.erase(str.find_last_not_of(dropChars)+1);
  }
  static void ltrim(std::string& str, const std::string dropChars = " \t\n\r") {
    str.erase(0, str.find_first_not_of(dropChars));
  }
  static std::string IntToStr(int integer) {
    std::ostringstream stream;
    stream << integer;
    return stream.str();
  }
  static int splitToStr(const char * str, 
                           std::vector<std::string> & items, 
                           const char * delm = "\t") {
    char * buff = const_cast<char *>(str);
    items.clear();
    char * pch = strtok(buff, delm);
    while( pch != NULL ) {
      items.push_back(pch);
      pch = strtok(NULL, delm);
    }
    return items.size();
  }
  static int splitToStr(const std::string& buff, 
                           std::vector<std::string> & items, 
                           std::string delm = "\t") {
    std::string cp = buff.substr();
    return splitToStr(cp.c_str(), items, delm.c_str());
  }
  static int splitToStrMD(const std::string& buff,
      std::vector<std::string>& items, const std::string delm) {
    // find 'delm' in string
    items.clear();
    size_t lindex(0), rindex = buff.find(delm);
    while(rindex != std::string::npos) {
      items.push_back(buff.substr(lindex,rindex-lindex));
      lindex = rindex + delm.size();
      rindex = buff.find(delm, lindex);
    }
    if(lindex < buff.size())
      items.push_back(buff.substr(lindex));
    return items.size();
  }
  static int splitToIntMD(const std::string& buff, std::vector<int>& items, 
                           std::string delm = ",") {
    items.clear();
    std::vector<std::string> tmpVector;
    /*char copy[buff.size()];
    buff.copy(copy, buff.size());
    copy[buff.size()] = '\0';
    int i = splitToStrMD(copy, tmpVector, delm.c_str());*/
    int i = splitToStrMD(buff.c_str(), tmpVector, delm.c_str());
    if( i > 0 )
      for( int j = 0; j < i; j++ )
        items.push_back(atoi(tmpVector[j].c_str()));
    return i;
  }
  static int splitToInt(const std::string& buff, std::vector<int>& items, 
                           std::string delm = ",") {
    items.clear();
    std::vector<std::string> tmpVector;
    /*char copy[buff.size()];
    buff.copy(copy, buff.size());
    copy[buff.size()] = '\0';
    int i = splitToStr(copy, tmpVector, delm.c_str());*/
    int i = splitToStr(buff.c_str(), tmpVector, delm.c_str());
    if( i > 0 )
      for( int j = 0; j < i; ++j )
        items.push_back(atoi(tmpVector[j].c_str()));
    return i;
  }
  static int splitToFloat(const std::string& buff, std::vector<float>& items, 
                           std::string delm = ",") {
    items.clear();
    std::vector<std::string> tmpVector;
    char copy[buff.size()];
    buff.copy(copy, buff.size());
    copy[buff.size()] = '\0';
    int i = splitToStr(copy, tmpVector, delm.c_str());
    //int i = splitToStr(buff.c_str(), tmpVector, delm.c_str());
    if( i > 0 )
      for(int j=0; j < i; j++ )
        items.push_back(atof(tmpVector[j].c_str()));
    return i;
  }
  static std::string lowercase(const std::string& str) {
    std::string retstr(str);
   for(size_t i=0; i < str.length(); i++) {
      retstr[i] = tolower(str[i]);
   }
   return retstr;
  }
  static std::string uppercase(const std::string& str) {
    std::string retstr(str);
   for(size_t i=0; i < str.length(); i++) {
      retstr[i] = toupper(str[i]);
   }
   return retstr;
  }
  // TODO: interface with decent PRG
  template<typename T>
  static T rand(T mod_bnd = 0) {
    T random = 0;
    if(sizeof(T) <= 4) {
      random = static_cast<T>(std::rand());
    }
    else if(sizeof(T) == 8) {
      random = static_cast<T>(std::rand());
      random <<= 31; random <<= 1;
      random |= static_cast<T>(std::rand());
    }
    if(mod_bnd != 0) 
      return random % mod_bnd;
    else return random;
  }
  static double dRand() {
    return drand48();
    //return rand() * (1/double(RAND_MAX));
  }
  static int iRand(int modval=1) {
    return lrand48() % modval;
    //return rand() % modval;
  }
  static void inline sleep(const size_t seconds) {
    usleep(seconds);
  }
  static bool fileExists(const std::string path) {
    bool exists = false;
    struct stat f_info;
    if( stat(path.c_str(), &f_info) == 0 ) //if stat() returns no errors
      exists = true;
    return( exists );
  }
  static double choose(unsigned n, unsigned k) { // binomial coefficient
    return exp(lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1));
  }
};
#endif //STREAM_UTILS_H
