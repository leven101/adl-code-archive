#ifndef STREAM_TYPES_H
#define STREAM_TYPES_H
#include <iostream>
#include <map>
#include <set>
#include <vector>
#include <typeinfo>
#include <stdint.h>

#define iterate(c, i) for(typeof(c.begin()) i = c.begin(); i != c.end(); ++i)
#define piterate(c, i) for(typeof(c->begin()) i = c->begin(); i != c->end(); ++i)
#define riterate(c, i) for(typeof(c.rbegin()) i = c.rbegin(); i != c.rend(); ++i)
#define THREADED false
#define THREAD_MAX 2
#define MAX_NGRAM_ORDER 8
#define MAX_STR_LEN 300
#define PRIME 8589935681ULL
#define MAX_HASH_FUNCS 1000
//#define PRIME 409 

using std::string;
using std::cout;
using std::cerr;
using std::endl;

//typedefs for projects
typedef std::string word_t;     // word as string 
typedef uint32_t wordID_t;      // word mapped to integer
typedef std::string date_t;     // a date marker
typedef unsigned int count_t;   // for 64-bit to 32-bit compatibility 
#endif // STREAM_TYPES_H
