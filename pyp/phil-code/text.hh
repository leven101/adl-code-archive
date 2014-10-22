#ifndef GTL_VOCAB_H
#define GTL_VOCAB_H

// STL
#include <vector>
#include <iostream>

// Boost
#include <boost/tuple/tuple.hpp>
#include <boost/bimap.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/functional/hash.hpp>


typedef std::string             Token;
typedef std::vector<Token>      TokenString;
typedef int                     IndexedToken;
typedef std::vector<IndexedToken> IndexedString;

class Vocab
{
  public:
    Vocab() {}

    IndexedString
      index(const TokenString &in)
      {
        // lookup the index for each word. If the word has not been seen, insert
        // adds it to the dictionary
        IndexedString result;
        for(TokenString::const_iterator it=in.begin(); it != in.end(); ++it)
          result.push_back(index(*it));
        return result;
      }

    IndexedToken
      index(const Token &in)
      {
        // lookup the index for the word. If the word has not been seen, insert
        // adds it to the dictionary
        bool inserted;
        TokenDictionary::iterator result;
        boost::tie(result, inserted) = m_token_index.insert(TokenDictionary::relation(in, m_token_index.size()));
        return result->right;
      }

    TokenString
      lookup(const IndexedString &in) const
      {
        // lookup the index for each word.
        TokenString result;
        for(IndexedString::const_iterator it=in.begin(); it != in.end(); ++it)
          result.push_back(lookup(*it));
        return result;
      }

    Token
      lookup(const IndexedToken &in) const
      {
        TokenDictionary::right_const_iterator result;
        result = m_token_index.right.find(in);
        if (result == m_token_index.right.end())
        {
          std::cerr << "<UNKNOWN_TOKEN_INDEX," << boost::lexical_cast<std::string>(in) << ">" << std::endl;
          assert(false);
          return "<UNKNOWN_TOKEN_INDEX," + boost::lexical_cast<std::string>(in) + ">";
        }
        return result->second;
      }

    bool 
      read(std::istream& in)
      {
        // table format is index tokend
        std::string buf;
        while (std::getline(in, buf))
        {
          if (buf.empty()) return false;

          // split line into fields
          std::istringstream ss(buf);

          Token token;
          IndexedToken indexed_token;
          ss >> indexed_token >> token;
          if (!ss) return false;

          bool inserted;
          TokenDictionary::iterator result;
          boost::tie(result, inserted) = m_token_index.insert(TokenDictionary::relation(token, indexed_token));
          if (inserted == false)
            return false;
        }
        return true;
      }

    std::ostream&
      write(std::ostream& out) {
        for (TokenDictionary::const_iterator it=m_token_index.begin(); it != m_token_index.end(); ++it)
          out << it->right << " " << it->left << std::endl;
        return out;
      }

  protected:
    typedef boost::bimap<Token,IndexedToken> TokenDictionary;

    TokenDictionary m_token_index;
};
typedef boost::shared_ptr<Vocab> VocabPtr;


class Sentence
{
  public:
    friend std::istream& 
      operator>>(std::istream& in, Sentence& s);
    friend std::ostream& 
      operator<<(std::ostream& out, const Sentence& s);

  public:
    Sentence(VocabPtr v) : m_vocab(v) {}
    ~Sentence() {}

    IndexedToken
      word(unsigned i) const
      {
        return m_indexed_tokens.at(i);
      }

    // get the sub-sequence [i,j), i.e. span(0,1) returns the first token of the sentence.  
    IndexedString
      span(unsigned i, unsigned j) const
      {
        if (i > j) { std::cerr << i << " " << j << std::endl; }
        assert(i <= j);

        IndexedString result;
        IndexedString::const_iterator it=m_indexed_tokens.begin() + i;
        IndexedString::const_iterator end
          = (j < m_indexed_tokens.size() ? it + (j-i) : m_indexed_tokens.end());

        for(; it != end; ++it)
          result.push_back(*it);

        return result;
      }

    bool
      empty() const
      { return m_indexed_tokens.empty(); }

    unsigned
      length() const
      { return m_indexed_tokens.size(); }

    VocabPtr
      vocab() const
      { return m_vocab; }

    IndexedString
      indexed_string() const
      { return m_indexed_tokens; }

  protected:
    IndexedString   m_indexed_tokens;
    VocabPtr        m_vocab;
};

inline std::istream& operator>>(std::istream& in, Sentence& s)
{
  std::string buf, token;
  TokenString tokens;

  if (!std::getline(in, buf)) 
    return in;

  std::istringstream ss(buf);

  while(ss >> token) 
    tokens.push_back(token);

  s.m_indexed_tokens = s.m_vocab->index(tokens);

  return in;
}

inline std::ostream& operator<<(std::ostream& out, const Sentence& s)
{
  IndexedString::const_iterator it=s.m_indexed_tokens.begin();
  if (it != s.m_indexed_tokens.end())
  {
    out << s.m_vocab->lookup(*it);

    ++it;
    for (; it != s.m_indexed_tokens.end(); ++it)
      out << " " << s.m_vocab->lookup(*it);
  }

  return out;
}

inline std::ostream& operator<<(std::ostream& out, const IndexedString& s)
{
  IndexedString::const_iterator it=s.begin();
  if (it != s.end())
  {
    out << *it;

    ++it;
    for (; it != s.end(); ++it)
      out << " " << *it;
  }

  return out;
}

inline std::ostream& operator<<(std::ostream& out, const TokenString& s)
{
  TokenString::const_iterator it=s.begin();
  if (it != s.end())
  {
    out << *it;

    ++it;
    for (; it != s.end(); ++it)
      out << " " << *it;
  }

  return out;
}

inline std::string str(const TokenString& s)
{
  std::string result;
  TokenString::const_iterator it=s.begin();
  if (it != s.end())
  {
    result += *it;

    ++it;
    for (; it != s.end(); ++it)
      result += " " + *it;
  }

  return result;
}

// Container for a pair of sentences
class SentencePair
{
  public:
    SentencePair(const Sentence& s, const Sentence& t)
      : source(s), target(t) {}

    Sentence source;
    Sentence target;
};

inline std::size_t 
hash_value(IndexedString const& s) {
  static boost::hash<IndexedString> indexed_string_hasher;

  return indexed_string_hasher(s);
}

inline std::string
map_brakets(const std::string &str) {
  std::string new_str;
  for (std::string::const_iterator it=str.begin(); it != str.end(); ++it) {
    if (*it == '(') new_str.append("-LRB-");
    else if (*it == ')') new_str.append("-RRB-");
    else new_str.push_back(*it);
  }
  return new_str;
}

inline IndexedString
span(const IndexedString& s, unsigned i, unsigned j)
{
  if (i > j) { std::cerr << i << " " << j << std::endl; }
  assert(i <= j);

  IndexedString result;
  IndexedString::const_iterator it=s.begin() + i;
  IndexedString::const_iterator end = (j < s.size() ? it + (j-i) : s.end());

  for(; it != end; ++it)
    result.push_back(*it);

  return result;
}

#endif // GTL_VOCAB_H
