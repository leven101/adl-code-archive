#include <libxml/HTMLparser.h>
#include <ctype.h>
#include "file.h"
#include "utils.h"
using namespace std;

#define COMPARE(a, b) (!strcasecmp((a), (b)))

struct Context
{
  Context(): addBody(false), addTitle(false) { }
  bool addBody, addTitle;
  std::string title;
  std::string body;
};
static void StartElement(void *voidContext, const xmlChar *name,
  const xmlChar **attributes) {
  Context *context = (Context *)voidContext;
  if(COMPARE((char *)name, "BODY")) {
    context->body = "";
    context->addBody = true;
  }
  else if(COMPARE((char *)name, "TITLE")) {
    context->title = "";
    context->addTitle = true;
  }
  (void) attributes;
}
static void EndElement(void *voidContext, const xmlChar *name)
{
  Context *context = (Context *)voidContext;
  if(COMPARE((char *)name, "BODY")) {
    context->addBody = false;
  }
  if(COMPARE((char *)name, "TITLE")) {
    context->addTitle = false;
  }
}
static void handleCharacters(Context *context, const xmlChar *chars, int length)
{
  if(context->addTitle) {
    context->title.append((char *)chars, length);
  }
  if(context->addBody) {
    size_t oldsize = context->body.size();
    string str = (char *)chars;
    Utils::trim(str);
    if(str.size() && (str.size() < 400)) {
      //context->body += str + "\n";
      for(size_t i=0; i < str.size(); ++i) {
        if(isalnum(str[i]) || isspace(str[i]) || ispunct(str[i]))
          context->body += str[i];
      }
      if(oldsize < context->body.size())
        context->body += "\n";
    }
  }
}
static void Characters(void *voidContext, const xmlChar *chars, int length){
  Context *context = (Context *)voidContext;
  handleCharacters(context, chars, length);
}
static void cdata(void *voidContext, const xmlChar *chars, int length)
{
  Context *context = (Context *)voidContext;
  handleCharacters(context, chars, length);
}

static htmlSAXHandler saxHandler =
{
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  StartElement,
  EndElement,
  NULL,
  Characters,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  cdata,
  NULL
};
static string parseHtml(const std::string &html)
{
  htmlParserCtxtPtr ctxt;
  Context context;

  ctxt = htmlCreatePushParserCtxt(&saxHandler, &context, "", 0, "",
        XML_CHAR_ENCODING_NONE);
  htmlCtxtUseOptions(ctxt, HTML_PARSE_NOBLANKS | HTML_PARSE_NOERROR | HTML_PARSE_NOWARNING | HTML_PARSE_NONET);

  htmlParseChunk(ctxt, html.c_str(), html.size(), 0);
  htmlParseChunk(ctxt, "", 0, 1);

  htmlFreeParserCtxt(ctxt);

  return (context.title + "\n" + context.body + "\n");
}


