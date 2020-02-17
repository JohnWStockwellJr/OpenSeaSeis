/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.general;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

import javax.swing.JOptionPane;
import javax.swing.filechooser.FileSystemView;

/**
 * Base class for user preferences/configuration file.
 * @author Bjorn Olofsson
 */
public abstract class csAbstractPreference {
  public static final boolean IS_WINDOWS = ( java.io.File.separatorChar == '\\' );
  public static final String FILE_NAME_START = "preferences_";

  private static File myCseisDirectory;

  private PrintWriter myErrorLogWriter;
  protected File myPreferencesFile;
  protected File myErrorFile;
  private ArrayList<String> myExceptionErrorList;
  
  public static File getCseisDirectory() {
    return myCseisDirectory;
  }
  public static File getCseisDirectoryFile() {
    FileSystemView fsv = FileSystemView.getFileSystemView();
    String folder = ".cseis";
    File homeDirectory = fsv.getHomeDirectory();
    if( IS_WINDOWS ) {
      // Trick to get to C:\\Documents and Settings\\user_name :
      File dirFullPath = new File( homeDirectory, "" );
      homeDirectory = dirFullPath.getParentFile();
    }
    File cseisDirectory = new File( homeDirectory, folder );
    return cseisDirectory;
  }
  public static boolean setCseisDirectory() throws Exception {
    myCseisDirectory = getCseisDirectoryFile();

    if( !myCseisDirectory.exists() || !myCseisDirectory.isDirectory() ) {
      if( !myCseisDirectory.mkdir() ) {
        JOptionPane.showMessageDialog( null,
            "Error occurred when trying to create Cseis configuration directory\n" +
            myCseisDirectory + "\nProgram will be terminated.",
            "", JOptionPane.ERROR_MESSAGE );
        throw new Exception( "Cannot create configuration directory" );
      }
      return true;
    }
    return false;
  }

  //--------------------------------------
  public csAbstractPreference( String name ) {
    myPreferencesFile = new File( getCseisDirectory(), "pref_" + name + ".txt" );
    myExceptionErrorList = new ArrayList<String>();
    myErrorFile = new File( getCseisDirectory(), "pref_" + name + ".err" );
  }
  public boolean readPreferences() {
    if( myPreferencesFile.exists() ) {
      readPreferences( myPreferencesFile );
      return true;
    }
    else {
      return false;
    }
  }
  public void readPreferences( File file ) {
    openErrorLogWriter();
    myExceptionErrorList.clear();
    BufferedReader reader = null;
    try {
      reader = new BufferedReader( new FileReader( file ) );
      readPreferences(reader);
      reader.close();
    }
    catch( IOException e ) {
      if( myErrorLogWriter != null ) e.printStackTrace(myErrorLogWriter);
      if( reader != null ) {
        try {
          reader.close();
        }
        catch( IOException exc ) {
          if( myErrorLogWriter != null ) exc.printStackTrace(myErrorLogWriter);
        }
      }
    }
    closeErrorLogWriter();
  }
//----------------------------------------------
  public void writePreferences() {
    writePreferences( myPreferencesFile );
  }
  public void writePreferences( File file ) {
    openErrorLogWriter();
    myExceptionErrorList.clear();

    BufferedWriter writer = null;
    try {
      writer = new BufferedWriter( new FileWriter( file ) );
      writePreferences(writer);
      writer.close();
    }
    catch( IOException e ) {
      if( myErrorLogWriter != null ) e.printStackTrace(myErrorLogWriter);
      if( writer != null ) {
        try {
          writer.close();
        }
        catch( IOException exc ) {
          if( myErrorLogWriter != null ) exc.printStackTrace(myErrorLogWriter);
        }
      }
    }
    closeErrorLogWriter();
  }
  protected abstract void readPreferences(BufferedReader reader) throws IOException;
  protected abstract void writePreferences(BufferedWriter writer) throws IOException;
//----------------------------------------------
  private void closeErrorLogWriter() {
    if( myErrorLogWriter != null ) {
      myErrorLogWriter.close();
      myErrorLogWriter = null;
    }
  }
//----------------------------------------------
  private void openErrorLogWriter() {
    try {
      myErrorLogWriter = new PrintWriter( new BufferedWriter( new FileWriter( myErrorFile ) ) );
    }
    catch( Exception e ) {
      if( myErrorLogWriter != null ) myErrorLogWriter = null;
    }
    myErrorLogWriter.close();
    myErrorLogWriter = null;
  }  

}


