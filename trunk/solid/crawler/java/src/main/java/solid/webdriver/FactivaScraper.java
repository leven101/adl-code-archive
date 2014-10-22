package solid.webdriver;

import java.util.Iterator;
import java.util.Set;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.TimeUnit;
import java.util.Date;
import java.io.File;
import java.io.FileWriter;
import java.io.BufferedWriter;
import java.io.IOException;

import org.openqa.selenium.By;
import org.openqa.selenium.WebDriver;
import org.openqa.selenium.WebElement;
import org.openqa.selenium.firefox.FirefoxDriver;
import org.openqa.selenium.firefox.FirefoxProfile;
import org.openqa.selenium.support.ui.ExpectedCondition;
import org.openqa.selenium.support.ui.Select;
import org.openqa.selenium.support.ui.WebDriverWait;

public class FactivaScraper {
  static String foldername_, year_;
  final static String searchQuery = "nonfarm near20 (predict$ OR forecast$)";
  // ("(unemploy$ OR employ$) near50 economy");
  final static FirefoxProfile firefoxProfile = new FirefoxProfile();
  public static void createDataFolder() {
    Date date = new Date();
    foldername_ = "data/Factiva/nonfarm-pred/" + year_
      + "/" + String.valueOf(date.getTime());
    File dir = new File(foldername_);
    if(dir.exists()) {
      System.out.println("ERROR: Output directory " + foldername_ + " already exists.");
      System.exit(1);
    }
    boolean success = dir.mkdirs();
    if(!success) {
      System.out.println("ERROR: failed to create output directory " + foldername_ + ".");
      System.exit(1);
    }
    else {
      System.out.println("Created output directory " + foldername_ + ".");
    }
  }
  public static void saveResults(String fileName, String s) {
    try {
      BufferedWriter out = new BufferedWriter(new FileWriter(foldername_ + "/" + fileName));
      out.write(s);
      out.close();
    } catch (IOException e) {
      e.printStackTrace();
    }
  }
  public static void crawlArticles(String month, String frd, String tod) throws InterruptedException {
    FirefoxDriver driver = new FirefoxDriver(firefoxProfile);
    driver.manage().deleteAllCookies();
    driver.manage().timeouts().implicitlyWait(10, TimeUnit.SECONDS); //http://seleniumhq.org/docs/04_webdriver_advanced.jsp
    try {
      // get through login page
      /*driver.get("http://oxford1-direct.hosted.exlibrisgroup.com/V/3PV6CKJ26USICFNEVXV5KG6K1KL5SHQRFCGH32DNECR6PDT5YN-39390?func=native-link&resource=OXF00954");
      driver.findElement(By.id("username")).sendKeys("coml0399"); 
      driver.findElement(By.id("password")).sendKeys("coml0399Password1"); 
      driver.findElement(By.name("Submit")).click();
      driver.findElementByXPath("/html/body/div/div/table/tbody/tr[3]/td/p/span/a").click();*/
      // setup search
      // keywords
      driver.findElementByXPath("//*[@id='ftx']").sendKeys(searchQuery);
      // date range
      Select selectDate = new Select(driver.findElementByXPath("//*[@id='dr']"));
      //selectDate.selectByVisibleText("In the last day");
      selectDate.selectByVisibleText("Enter date range...");
      driver.findElement(By.id("frd")).sendKeys(frd);
      driver.findElement(By.id("frm")).sendKeys(month);
      driver.findElement(By.id("fry")).sendKeys(year_);
      String toYear;
      if(Integer.valueOf(month) == 12) { 
        month = "1"; 
        toYear = String.valueOf(Integer.valueOf(year_)+1);
      }
      else {
        month = String.valueOf(Integer.valueOf(month)+1);
        toYear = year_;
      }
      driver.findElement(By.id("tod")).sendKeys(tod);
      driver.findElement(By.id("tom")).sendKeys(month);
      driver.findElement(By.id("toy")).sendKeys(toYear);
      // select region
      driver.findElement(By.id("reTab")).click();
      driver.findElement(By.xpath("/html/body/form[2]/div[3]/div[2]/table/tbody/tr[16]/td[2]/div[2]/ul/li[15]/span")).click();
      driver.findElement(By.xpath("/html/body/form[2]/div[3]/div[2]/table/tbody/tr[16]/td[2]/div[2]/ul/li[15]/div/ul/li[6]/a[1]")).click();
      TimeUnit.SECONDS.sleep(3);
      
      // click Search
      driver.findElementByXPath("/html/body/form[2]/div[3]/div[2]/table/tbody/tr[2]/td/div/div/table/tbody/tr/td[2]/div[3]/div[2]/ul/li/div/span").click();
      // skip this click to get Publications  
      //driver.findElement(By.xpath("//*[@id='headlineTabs']/table[1]/tbody/tr/td/span[1]/a")).click();
      // get number of total items 
      String sAllTxt = driver.findElement(By.xpath("//*[@id='headlineHeader33']/table/tbody/tr/td/span[2]")).getText(); 
      int start = sAllTxt.indexOf(" of ");
      String snum = sAllTxt.substring(start+4, sAllTxt.length()).replace(",","");
      // parse out number of documents here
      int totArticles = Integer.parseInt(snum);
      int totPages = (totArticles / 100);
      if(totArticles % 100 != 0) ++totPages;
      int curPage=1;
      System.out.println("btnAll text:" + sAllTxt + "  totArticles:" + totArticles + "  totPages:" + totPages);

      // iterate through all pages of news 
      HashSet<String> dups= new HashSet<String>();
      while(curPage <= totPages)  {
        Iterator<WebElement> a_itr = driver.findElements(By.xpath("//a[@class='enHeadline']")).iterator();
        StringBuilder builder = new StringBuilder(20000000);
        int stsize=0;
        while(a_itr.hasNext()) {
          WebElement a = a_itr.next();
          a.click();
          try {
            String text = driver.findElement(By.xpath("//*[@id='articleFrame']/div[2]/div[2]")).getText();
            text.trim();
            stsize += text.length();
            //System.out.println("total string size: " + stsize);
            builder.ensureCapacity(stsize);
            //System.out.println("Capacity: " + builder.capacity());
            builder.append(text);
            TimeUnit.SECONDS.sleep(1);
          } catch (Exception NoSuchElementException) {
            // just skip this article
          }
        }
        String fname = "data/Factiva/page" + String.valueOf(curPage) + ".txt"; 
        if(curPage==1) {
          createDataFolder();
        }
        saveResults("page" + String.valueOf(curPage) + ".txt", builder.toString());
        if(curPage < totPages) {
          // click "Next 100" 
          driver.findElement(By.xpath("//div[@class='headlineHeader']//a[last()]")).click();
          TimeUnit.SECONDS.sleep(2);
        }
        ++curPage;
      }
      driver.findElement(By.xpath("/html/body/div[2]/div/ul[2]/li/a")).click();
      driver.findElement(By.xpath("/html/body/div[2]/div/ul[2]/li/div/ul/li[4]/a")).click();
    }
    finally {
      final Set<String> wins = driver.getWindowHandles();
      for (final String string : wins)
        driver.switchTo().window(string).close();
    }
  }
  public static void main(final String[] args) throws InterruptedException {
    // instantiate a browser either mozilla or htmlunit
    //logger.info("Instantiate browser");
    firefoxProfile.setEnableNativeEvents(true);
    // disable cache
    firefoxProfile.setPreference("browser.cache.disk.enable", false);
    firefoxProfile.setPreference("browser.cache.memory.enable", false);
    firefoxProfile.setPreference("browser.cache.offline.enable", false);
    firefoxProfile.setPreference("network.http.use-cache", false);

    /*firefoxProfile.setPreference("browser.download.folderList", 2);
    firefoxProfile.setPreference("browser.download.manager.showWhenStarting", false);
    firefoxProfile.setPreference("browser.download.dir", "/home/ablev/workspace/solid.oxpathclient/data/Factiva/");
    firefoxProfile
        .setPreference("browser.helperApps.neverAsk.saveToDisk",
           "application/rtf,application/x-rtf,text/rtf,text/richtext,application/msword,application/doc,application/x-soffice");*/
    
    try {
      HashMap<Integer, Integer> months = new HashMap<Integer, Integer>();
      months.put(1, 31);
      months.put(2, 28);
      months.put(3, 31);
      months.put(4, 30);
      months.put(5, 31);
      months.put(6, 30);
      months.put(7, 31);
      months.put(8, 31);
      months.put(9, 30);
      months.put(10, 31);
      months.put(11, 30);
      months.put(12, 31);
      for(int y=2010; y <= 2012; ++y) { 
        year_ = String.valueOf(y);
        for(int m=1; m <= 12; ++m) { 
          String frd = "1";
          String tod = String.valueOf(months.get(m));
          System.out.println("Year:" + y + "\tMonth:" + m + "\tFrom: " + frd + "\tTo: " + tod);
          crawlArticles(String.valueOf(m), frd, tod);
        }
      }
    }
    finally {
    }
  }
}

/*
      //This requires CAPTCHA breaking
      // iterate through all pages of news 
      while(curPage <= totPages)  {
        // select all headlines 
        driver.findElement(By.xpath("/html/body/form[2]/div[3]/div[2]/div[3]/div[3]/div/div/div/div/table/tbody/tr/td/span/input")).click();
        // click rtfbutton
        driver.findElement(By.xpath("/html/body/form[2]/div[3]/div[2]/div[2]/table[2]/tbody/tr/td/div/div/span[3]/ul/li[5]/a")).click();
        // click 'article format'
        driver.findElement(By.xpath("/html/body/form[2]/div[3]/div[2]/div[2]/table[2]/tbody/tr/td/div/div/span[3]/ul/li[5]/ul/li[2]/a")).click();
        TimeUnit.SECONDS.sleep(10);
        // deselect all headlines 
        driver.findElement(By.xpath("/html/body/form[2]/div[3]/div[2]/div[3]/div[3]/div/div/div/div/table/tbody/tr/td/span/input")).click();
        TimeUnit.SECONDS.sleep(2);
        if(curPage < totPages) {
          // click "Next 100" 
          driver.findElement(By.xpath("//div[@class='headlineHeader']//a[last()]")).click();
        }
        ++curPage;
      }
 */
