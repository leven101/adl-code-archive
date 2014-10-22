#!/usr/bin/python

import os;
import sys;
import time;
from selenium import webdriver
from selenium.common.exceptions import TimeoutException
from selenium.webdriver.support.ui import WebDriverWait # available since 2.4.0
from selenium.webdriver.support import expected_conditions as EC # available since 2.26.0
from selenium.webdriver.support.ui import Select

searchQuery = "NONE";
months = {1:31, 2:28, 3:31, 4:30, 5:31, 6:30, 7:31, 8:31, 9:30, 10:31, 11:30, 12:31};

def createDataFolder():
  global folder; 
  folder = './data/' + str(time.time()).replace('.', '');
  try:
    os.makedirs(folder)
  except OSError:
    raise 

def save(fname, text):
  f = open(folder+"/"+fname, 'w');
  f.write(text.encode('utf8'));
  f.close();

def login(driver):
  fin = open('params.txt', 'r').readlines();
  driver.get("http://oxford1-direct.hosted.exlibrisgroup.com/V/3PV6CKJ26USICFNEVXV5KG6K1KL5SHQRFCGH32DNECR6PDT5YN-39390?func=native-link&resource=OXF00954");
  driver.find_element_by_id("username").send_keys(fin[0].strip());
  driver.find_element_by_id("password").send_keys(fin[1].strip()); 
  driver.find_element_by_name("Submit").click();
  driver.find_element_by_xpath("/html/body/div/div/div[2]/a").click();    
  global searchQuery; 
  searchQuery = fin[2].strip();

def login2(driver):
  driver.get("http://oxford1-direct.hosted.exlibrisgroup.com/V/CXYFNXG9EEY7QMRBJSFALGLYRSAS31F9TFH9CRRHL2Q4FB1CP1-10639?func=native-link&resource=OXF00954")
  time.sleep(1)

def crawl(year, frm, tom, frd, tod):
  try:
    createDataFolder();
    driver = webdriver.Firefox(); # Create a new instance of the Firefox driver
    driver.implicitly_wait(10);
    login(driver);
    # keywords
    driver.find_element_by_xpath("//*[@id='ftx']").send_keys(searchQuery);
    # date ranges
    selectDate = Select(driver.find_element_by_xpath("//*[@id='dr']"));
    selectDate.select_by_visible_text("Enter data range...");
    driver.find_element_by_id('frd').send_keys(frd);
    driver.find_element_by_id('frm').send_keys(frm);
    driver.find_element_by_id('fry').send_keys(year);
    driver.find_element_by_id('tod').send_keys(tod);
    driver.find_element_by_id('tom').send_keys(tom);
    driver.find_element_by_id('toy').send_keys(year);
    # select region
    driver.find_element_by_id("reTab").click();
    driver.find_element_by_xpath("/html/body/form[2]/div[3]/div[2]/table/tbody/tr[16]/td[2]/div[2]/ul/li[15]/span").click();
    driver.find_element_by_xpath("/html/body/form[2]/div[3]/div[2]/table/tbody/tr[16]/td[2]/div[2]/ul/li[15]/div/ul/li[6]/a[1]").click();
    time.sleep(3);
    # click search
    driver.find_element_by_xpath("/html/body/form[2]/div[3]/div[2]/table/tbody/tr[2]/td/div/div/table/tbody/tr/td[2]/div[3]/div[2]/ul/li/div/span").click();
    # get total items to crawl
    sAllTxt = driver.find_element_by_xpath("//*[@id='headlineHeader33']/table/tbody/tr/td/span[2]").text;
    totArticles = int(sAllTxt.split(" ")[-1]);
    totPages = totArticles / 100;
    if(totArticles % 100 != 0):
      totPages +=1;
    curPage=1;
    while(curPage <= totPages):
      alltxt = "";
      a_list = driver.find_elements_by_xpath("//a[@class='enHeadline']");
      cnt = 0;
      for a in a_list:
        cnt += 1;
        if(cnt > 5):
          break;
        a.click(); 
        txt = driver.find_element_by_xpath("/html/body/form[2]/div[2]/div[2]/div[3]/div[3]/div/div[2]/div/div[2]/div[3]").text;
        alltxt += txt; 
      save("page" + str(curPage), alltxt);
      if(curPage < totPages): # click "Next 100"
        driver.find_element_by_xpath("//div[@class='headlineHeader']//a[last()]").click();
        time.sleep(2);
      curPage += 1;
  finally:
    driver.quit()

def setupCrawler():
  fyr = int(sys.argv[1]);
  tyr = int(sys.argv[2]);
  timeframe = sys.argv[3];
  for y in range(fyr, tyr+1):
    if(timeframe == 'quarter'):
      for m in range(1, 13, 3):
        tom = m+2; 
        tod = months[tom];
        print 'Year:', y, 'From Month:', m, 'To Month:', tom, 'From Day:', frd, 'To Day:', tod;
        #crawl(y, m, tom, frd, tod);
    for m in range(1, 13):
      tod = months[m];
      if((m==2) and (y % 4 == 0)):
        tod+=1; 
      if(timeframe == 'day'): 
        for frd in range(1, tod):
          print 'Year:', y, 'From Month:', m, 'To Month:', m, 'From Day:', frd, 'To Day:', frd+1;
        #crawl(y, m, m, frd, frd+1);
      else: #default is monthly
        frd=1;
        print 'Year:', y, 'From Month:', m, 'To Month:', m, 'From Day:', frd, 'To Day:', tod;
        #crawl(y, m, m, frd, tod);
  
def main():
  assert len(sys.argv) == 4;
  setupCrawler();
  
main();

