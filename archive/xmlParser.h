#ifndef INC_MYXMLPARSER_H
#define INC_MYXMLPARSER_H
/**
 * section: xmlReader
 * synopsis: Parse an XML file with an xmlReader
 * purpose: Demonstrate the use of xmlReaderForFile() to parse an XML file
 *          and dump the informations about the nodes found in the process.
 *          (Note that the XMLReader functions require libxml2 version later
 *          than 2.6.)
 * usage: reader1 <filename>
 * test: reader1 test2.xml > reader1.tmp ; diff reader1.tmp reader1.res ; rm reader1.tmp
 * author: Daniel Veillard
 * copy: see Copyright for the status of this software.
 */
#include <iostream>
#include <cstdio>
#include <libxml/xmlreader.h>

#ifdef LIBXML_READER_ENABLED

/**
 * processNode:
 * @reader: the xmlReader
 *
 * Dump information about the current node
 */
class XMLParser {
public:
  XMLParser() {
    // initializes and checks the xml library 
    LIBXML_TEST_VERSION
  }
  ~XMLParser() {
    // Cleanup function for the XML library.
    xmlCleanupParser();
    // this is to debug memory for regression tests
    // xmlMemoryDump();
  }

  const std::string
  processNode(xmlTextReaderPtr reader) {
    static bool inText = false, read = false,
      headline = false;
    const xmlChar *name;
    std::string date(""), value("");
    name = xmlTextReaderConstName(reader);
    if (name == NULL)
        name = BAD_CAST "--";
    
    if( xmlStrEqual(name, (const xmlChar *)"newsitem") &&
        xmlTextReaderHasAttributes(reader) &&
        xmlTextReaderNodeType(reader) == 1 ) {
      date = (const char *)xmlTextReaderGetAttribute(reader, 
                                   (const xmlChar *)"date");
      return date + "\n";
    }
    if( xmlStrEqual(name, (const xmlChar *)"headline"))
      headline = !headline;
    if(headline && xmlStrEqual(name, (const xmlChar *)"#text")) {
      value = (const char *) xmlTextReaderConstValue(reader);
      return value + "\n";
    }
    
    if( xmlStrEqual(name, (const xmlChar *)"text") ) 
      inText = !inText;
    if( xmlStrEqual(name, (const xmlChar *)"p") ) 
      read = !read;
    
    if( inText && read && 
        xmlStrEqual(name, (const xmlChar *)"#text") ) {
      value = (const char *) xmlTextReaderConstValue(reader);
      return value + "\n";
    }
    return "";
  }

  /**
   * streamFile:
   * @filename: the file name to parse
   *
   * Parse and print information about an XML file.
   */
  std::string 
  streamXMLFile(std::string filename) {
    xmlTextReaderPtr reader;
    int ret;
    std::string packet = "";
    reader = xmlReaderForFile(filename.c_str(), NULL, 0);
    if (reader != NULL) {
      ret = xmlTextReaderRead(reader);
      while (ret == 1) {
        packet += processNode(reader);
        ret = xmlTextReaderRead(reader);
      }
      xmlFreeTextReader(reader);
      if (ret != 0) {
        fprintf(stderr, "%s : failed to parse\n", filename.c_str());
      }
    } 
    else {
      fprintf(stderr, "Unable to open %s\n", filename.c_str());
    }
    return packet;
  } 
};
/*
int main(int argc, char **argv) {
    if (argc != 2)
        return(1);

    // * this initialize the library and check potential ABI mismatches
    // * between the version it was compiled for and the actual shared
    // * library used.
    LIBXML_TEST_VERSION

    streamFile(argv[1]);

    // * Cleanup function for the XML library.
    xmlCleanupParser();
    // * this is to debug memory for regression tests
    xmlMemoryDump();
    return(0);
}
*/
#endif
#endif
