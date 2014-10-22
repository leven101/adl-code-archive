#ifndef curl_solid_h 
#define curl_solid_h
/***************************************************************************
 *                                  _   _ ____  _
 *  Project                     ___| | | |  _ \| |
 *                             / __| | | | |_) | |
 *                            | (__| |_| |  _ <| |___
 *                             \___|\___/|_| \_\_____|
 *
 * Copyright (C) 1998 - 2011, Daniel Stenberg, <daniel@haxx.se>, et al.
 *
 * This software is licensed as described in the file COPYING, which
 * you should have received as part of this distribution. The terms
 * are also available at http://curl.haxx.se/docs/copyright.html.
 *
 * You may opt to use, copy, modify, merge, publish, distribute and/or sell
 * copies of the Software, and permit persons to whom the Software is
 * furnished to do so, under the terms of the COPYING file.
 *
 * This software is distributed on an "AS IS" basis, WITHOUT WARRANTY OF ANY
 * KIND, either express or implied.
 *
 ***************************************************************************/
#include <curl/curl.h>

static size_t sentiment(char *ptr, size_t size, size_t nmemb, void *userdata) {
  int strsz = nmemb * size;
  string* data = (string*) userdata;
  for(int i=0; i < strsz; ++i) {
    data->push_back(ptr[i]);
  }
  return strsz; 
}

static string theySayAPI(const string& sntn, const int service=0) {
  CURL *curl;
  CURLcode res;
  string retstr;
  struct curl_slist *slist=NULL;
  //string postthis = "{\"text\":\"" + sntn +  "\",\"level\":\"sentence\"}";
  string postthis = "{\"text\":\"" + sntn +  "\"";
  //postthis += ",\"bias\":{\"positive\":1, \"negative\":1, \"neutral\":0}";
  //postthis += ",\"bias\":{\"neutral\":30.5}";
  postthis += "}";
  //cerr << "Sending \"" + sntn + "\" to TheySay API." << endl;
  curl = curl_easy_init();
  if(curl) {
    slist = curl_slist_append(slist, "Content-Type: application/json");
    slist = curl_slist_append(slist, "Accept: application/json");
    curl_easy_setopt(curl, CURLOPT_HTTPHEADER, slist);
    curl_easy_setopt(curl, CURLOPT_USERPWD, "levenberga@theysayanalytics.com:Disent-Ils8835");
    string call = "http://api.theysay.io/v1/";
    size_t (*fncptr)(char*, size_t, size_t, void*);
    switch(service) {
      case 1:
        call += "risk";
        break;
      default:
        call += "sentiment";
        fncptr = &sentiment;
    }
    if(service == -1) {
      curl_easy_setopt(curl, CURLOPT_HEADER, 1); // returns header information 
    }
    else {
      curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, fncptr);
      curl_easy_setopt(curl, CURLOPT_WRITEDATA, (void*)&retstr);
    }
    curl_easy_setopt(curl, CURLOPT_URL, call.c_str());
    curl_easy_setopt(curl, CURLOPT_POSTFIELDS, postthis.c_str());

    /* if we don't provide POSTFIELDSIZE, libcurl will strlen() by itself */
    curl_easy_setopt(curl, CURLOPT_POSTFIELDSIZE, (long)strlen(postthis.c_str()));

    /* Perform the request, res will get the return code */
    res = curl_easy_perform(curl);
    /* Check for errors */
    if(res != CURLE_OK) {
      fprintf(stderr, "curl_easy_perform() failed: %s\n", 
        curl_easy_strerror(res));
    }
    /* always cleanup */
    curl_easy_cleanup(curl);
    curl_slist_free_all(slist); 
  }
  return retstr;
}
#endif
