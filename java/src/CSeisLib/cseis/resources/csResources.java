/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.resources;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URL;
import javax.swing.ImageIcon;

/**
 * Access to external "resources", i.e. image files
 * @author Bjorn Olofsson
 */
public class csResources {
  public static ImageIcon getIcon( String name ) {
    return( getIcon(name,"") );
  }
  public static ImageIcon getIcon( String name, String subDirectory ) {
    String fullName;
    if( subDirectory.compareTo("") == 0 ) {
      fullName = name;
    }
    else {
      fullName = subDirectory + java.io.File.separatorChar + name;
    }
    URL url = csResources.class.getResource( fullName );
    if( url == null ) return null;
    return( new ImageIcon(url) );
  }
  public static String getText( String fileName ) {
    return( getText( fileName, "" ) );
  }
  public static String getText( String fileName, String subDirectory ) {
    String fullName;
    if( subDirectory.compareTo("") == 0 ) {
      fullName = fileName;
    }
    else {
      fullName = subDirectory + java.io.File.separatorChar + fileName;
    }
    URL url = csResources.class.getResource( fullName );
    if( url == null ) return null;

    try {
      BufferedReader reader = new BufferedReader( new InputStreamReader( url.openStream() ) );
      String text = "";
      String line;
      while( (line = reader.readLine()) != null ) {
        text += line + "\n";
      }
      reader.close();
      return text;
    }
    catch( IOException e ) {
    }
    return null;
  }
}


