#!/usr/bin/python

import os;
from selenium import webdriver
from selenium.common.exceptions import TimeoutException
from selenium.webdriver.support.ui import WebDriverWait # available since 2.4.0
from selenium.webdriver.support import expected_conditions as EC # available since 2.26.0
from selenium.webdriver.support.ui import Select

#os.system('cp data/comm.dat data/comm.backup');
urls = ['http://www.ofdp.org/continuous_contracts/data?exchange=CBT&symbol=O&depth=1', 'http://www.ofdp.org/continuous_contracts/data?exchange=CBT&symbol=S&depth=1', 'http://www.ofdp.org/continuous_contracts/data?exchange=CBT&symbol=W&depth=1', 'http://www.ofdp.org/continuous_contracts/data?exchange=CBT&symbol=C&depth=1']
cmd = 'echo ';
for url in urls:
  driver = webdriver.Firefox(); # Create a new instance of the Firefox driver
  driver.implicitly_wait(20);
  driver.get(url);
  text = driver.find_element_by_xpath("/html/body/div/div/table/tbody/tr[2]").text;
  print text;
  text = text.split(' ');
  cmd += text[1] + '/';
  cmd += text[4] + '/';
  cmd += text[5] + '\t';
  driver.quit();
cmd += ' >> data/comm.dat';
#print cmd;
os.system(cmd);
