package cseis.test;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.logging.Level;
import java.util.logging.Logger;

public class TestProcess {
  public static void main( String[] args ) {
    Runtime runtime = Runtime.getRuntime();
    Process processTEMP = null;
    try {
      processTEMP = runtime.exec("cd");
    }
    catch( IOException ex) {
      System.out.println( ex.getMessage() );
    }
    catch (Exception err) {
      err.printStackTrace();
    }
    InputStream stderr = processTEMP.getErrorStream();
    InputStream stdout = processTEMP.getInputStream();
    String line;
    BufferedReader brCleanUp = new BufferedReader (new InputStreamReader (stdout));
    try {
      while( (line = brCleanUp.readLine ()) != null ) {
        System.out.println(line);
      }
      brCleanUp.close();
      brCleanUp = new BufferedReader (new InputStreamReader (stderr));
      while( (line = brCleanUp.readLine ()) != null ) {
          System.out.println(line);
      }
      brCleanUp.close();
    }
    catch (IOException ex) {
      Logger.getLogger(TestProcess.class.getName()).log(Level.SEVERE, null, ex);
    }
    
  }
}
