package solid.webdriver;

import org.openqa.selenium.By;
import org.openqa.selenium.WebDriver;
import org.openqa.selenium.WebElement;
import org.openqa.selenium.firefox.FirefoxDriver;
import org.openqa.selenium.support.ui.ExpectedCondition;
import org.openqa.selenium.support.ui.WebDriverWait;
import java.util.concurrent.TimeUnit;
import java.util.Date;
import java.io.File;

public class WebdriverTest {
  static String foldername_;
  public static void createDataFolder() {
    Date date = new Date();
    foldername_ = "data/Factiva/" + String.valueOf(date.getTime());
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
  }
  public static void main(String[] args) throws InterruptedException {
    // Create a new instance of the Firefox driver
    // Notice that the remainder of the code relies on the interface, 
    // not the implementation.
    WebDriver driver = new FirefoxDriver();

    // And now use this to visit Google
    driver.get("http://www.google.com");
    // Alternatively the same thing can be done like this
    // driver.navigate().to("http://www.google.com");

    // Find the text input element by its name
    WebElement element = driver.findElement(By.name("q"));

    // Enter something to search for
    element.sendKeys("Cheese!");

    // Now submit the form. WebDriver will find the form for us from the element
    element.submit();

    // Check the title of the page
    System.out.println("Page title is: " + driver.getTitle());
    
    // Google's search is rendered dynamically with JavaScript.
    // Wait for the page to load, timeout after 10 seconds
    (new WebDriverWait(driver, 10)).until(new ExpectedCondition<Boolean>() {
        public Boolean apply(WebDriver d) {
            return d.getTitle().toLowerCase().startsWith("cheese!");
        }
    });

    // Should see: "cheese! - Google Search"
    System.out.println("Page title is: " + driver.getTitle());

    TimeUnit.SECONDS.sleep(2);
    //Close the browser
    driver.quit();
  }
}
