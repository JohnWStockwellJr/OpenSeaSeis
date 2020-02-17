/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.general;

/**
 * Defines one data file.<br>
 * Splits into three strings: The file name, the path name, and the full path name including the file name.
 * @author 2013 Felip Punto
 */
public class csFilename {
  /// File name, without path
  public String filename;
  /// File name with full path
  public String filenamePath;
  /// Path name / directory where file is located, without trailing file separator
  public String path;
  
  public csFilename( csFilename obj ) {
    filename = obj.filename;
    filenamePath = obj.filenamePath;
    path = obj.path;
  }
  public csFilename( String filenamePath_in ) {
    setFilenamePath( filenamePath_in );
  }
  public csFilename( String path_in, String filename_in ) {
    setFilenamePath( path_in, filename_in );
  }
  public csFilename() {
    filename = "";
    filenamePath = "";
    path = "";
  }
  public void setFilenamePath( String filenamePath_in ) {
    filenamePath = filenamePath_in;
    int index = filenamePath.lastIndexOf(java.io.File.separator);
    if( index < 0 || index >= filenamePath.length() ) {
      filename = filenamePath;
      path = "";
    }
    else {
      filename = filenamePath.substring( index+1, filenamePath.length() );
      path = filenamePath.substring( 0, index );
    }
  }
  public void setFilenamePath( String path_in, String filename_in ) {
    path = path_in;
    filename = filename_in;
    filenamePath = path + java.io.File.separator + filename;
  }
  public void dump() {
    System.out.println("FilenamePath:    " + filenamePath);
    System.out.println("Path + filename: '" + path + "' + '" + filename + "'");
  }
}

