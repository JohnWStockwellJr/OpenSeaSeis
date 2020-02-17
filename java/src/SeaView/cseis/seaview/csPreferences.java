/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seaview;

import cseis.general.csAbstractPreference;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;

/**
 * Preferences specific to seismic viewer application SeaView.<br>
 * Preferences are saved when the application is quit, and reloaded during next start-up.
 * @author 2008 Bjorn Olofsson
 */
public class csPreferences extends csAbstractPreference {
  public static final String FILE_HEADER = "### SEAVIEW v" + SeaView.VERSION + " preferences file ###";
  public static final String NODE_END            = "<end>";
  public static final String NODE_RECENT_FILES   = "<recent_files>";

  private java.util.List<String> myRecentFileList = null;
  
  public csPreferences() {
    super( "seaview" );
  }
  public java.util.List<String> getRecentFileList() {
    return myRecentFileList;
  }
  public void setRecentFileList( java.util.List<String> list ) {
    myRecentFileList = list;
  }
//----------------------------------------------
  @Override
  protected void readPreferences(BufferedReader reader) throws IOException {
    String line;
    myRecentFileList = new ArrayList<String>();
    while( (line = reader.readLine()) != null ) {
      if( line.equalsIgnoreCase(NODE_RECENT_FILES) ) {
        myRecentFileList = new ArrayList<String>();
        while( (line = reader.readLine()) != null && !line.trim().equalsIgnoreCase(NODE_END)) {
          myRecentFileList.add( line.trim() );
        }
      }

    }
  }
  @Override
  protected void writePreferences( BufferedWriter writer ) throws IOException {
    writer.write(FILE_HEADER);
    writer.newLine();
    writer.newLine();
    writer.write(NODE_RECENT_FILES);
    writer.newLine();
    if( myRecentFileList != null ) {
      for( int i = 0; i < myRecentFileList.size(); i++ ) {
        writer.write( (String)myRecentFileList.get(i) );
        writer.newLine();
      }
    }
    writer.write(NODE_END);
    writer.newLine();
    writer.newLine();
  }
}


