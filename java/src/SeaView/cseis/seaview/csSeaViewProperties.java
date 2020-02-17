/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seaview;

import cseis.general.csAbstractPreference;
import cseis.seisdisp.csAnnotationAttributes;
import cseis.general.csColorMap;
import cseis.general.csCustomColorMap;
import cseis.segy.csSegyAttr;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Properties;

/**
 * This class is used to store user specified parameters and settings to/from Java Properties file.<br>
 * Properties are saved when the application is quit, and reloaded during next start-up.
 * @author 2009 Bjorn Olofsson
 */
public class csSeaViewProperties {
  public String dataDirectoryPath;
  public String screenDumpDirectoryPath;
  public String propertiesDirectoryPath;
  public String velFieldDirectoryPath;
  public String wellPathDirectoryPath;
  public int windowLayout;

  public csSegyAttr segyAttr;
  public csAnnotationAttributes annAttr;
  public int numTraces;
  public int maxTabWidth;
  public boolean showFilename;
  public boolean showActivePaneFeedback;
  public int colorBitType;
  
  private File myPropertiesFile;

  public static final String PROPERTY_colorBitType     = "colorBitType";
  public static final String PROPERTY_endianSEGYHdr    = "endianSEGYHdr";
  public static final String PROPERTY_endianSEGYData   = "endianSEGYData";
  public static final String PROPERTY_endianSUHdr      = "endianSUHdr";
  public static final String PROPERTY_endianSUData     = "endianSUData";
  public static final String PROPERTY_isCustomMapSEGY  = "isCustomMapSEGY";
  public static final String PROPERTY_hdrMapSEGY       = "hdrMapSEGY";
  public static final String PROPERTY_autoScaleSEGY    = "autoScaleSEGY";
  public static final String PROPERTY_numTraces        = "numTraces";
  public static final String PROPERTY_maxTabWidth      = "maxTabWidth";
  public static final String PROPERTY_showFilename     = "showFilename";
  public static final String PROPERTY_showActivePaneFeedback = "showActivePaneFeedback";

  public static final String PROPERTY_windowLayout     = "windowLayout";
  
  public static final String PROPERTY_annotationOmitRepeating  = "omitRepeating";
  public static final String PROPERTY_annotationShowSequential  = "showSequential";
  public static final String PROPERTY_annotationFixedTraceLabelStep  = "fixedTraceLabelStep";
  public static final String PROPERTY_annotationTraceLabelStep  = "traceLabelStep";
  
  public static final String PROPERTY_dataDir          = "dataDir";
  public static final String PROPERTY_screenDumpDir    = "screenDumpDir";
  public static final String PROPERTY_propertyDir      = "propertyDir";
  public static final String PROPERTY_velFieldDir      = "velFieldDir";
  public static final String PROPERTY_wellPathDir      = "wellPathDir";
  
  private ArrayList<csCustomColorMap> myCustomColorMapList;
          
  public csSeaViewProperties() {
    myPropertiesFile  = createPropertyFile();
    //new File( csAbstractPreference.getCseisDirectory(), "prop_seaview_"+SeaView.VERSION+".txt" );
    dataDirectoryPath = csAbstractPreference.getCseisDirectory().getAbsolutePath();
    screenDumpDirectoryPath = csAbstractPreference.getCseisDirectory().getAbsolutePath();
    propertiesDirectoryPath = csAbstractPreference.getCseisDirectory().getAbsolutePath();
    windowLayout   = cseis.swing.csDockPaneManager.LAYOUT_ONE_ROW;
    
    colorBitType   = csColorMap.COLOR_MAP_TYPE_32BIT;
    
    numTraces = csTraceSelectionParam.NUM_TRACES_DEFAULT;
    maxTabWidth = 9999;
    showFilename = true;
    showActivePaneFeedback = true;
    segyAttr = new csSegyAttr();
    annAttr = new csAnnotationAttributes();

    myCustomColorMapList = new ArrayList<csCustomColorMap>();
    loadColorMaps();
  }
  private void loadColorMaps() {
    File[] allFiles = csAbstractPreference.getCseisDirectory().listFiles();
    ArrayList<String> cmapFilenameList = new ArrayList<String>();
    for( File file : allFiles ) {
      if( file.getName().endsWith(".cmap") ) {
        cmapFilenameList.add( file.getAbsolutePath() );
      }
    }
    // Sort CMAP files by name
    Collections.sort( cmapFilenameList );
    for( String filename : cmapFilenameList ) {
      csCustomColorMap cmap = new csCustomColorMap();
      try {
        cmap.load( new BufferedReader( new FileReader( new File(filename) ) ) );
        myCustomColorMapList.add( cmap );
      } catch( IOException ex ) {}
    }
  }
  public void updateColorMaps( java.util.List<csCustomColorMap> list ) {
    myCustomColorMapList.clear();
    for( int i = 0; i < list.size(); i++ ) {
      myCustomColorMapList.add( list.get(i) );
    }
  }
  public java.util.List<csCustomColorMap> retrieveColorMaps() {
    return myCustomColorMapList;
  }
  public String filename() {
    return myPropertiesFile.getAbsolutePath();
  }
  public static File createPropertyFile() {
    File file = csAbstractPreference.getCseisDirectoryFile();
    return new File( file, "prop_seaview_"+SeaView.VERSION+".txt" );
  }
  public static boolean fileExists() {
    File properties = createPropertyFile();
    return properties.exists();
  }
  public void load() throws IOException {
    Properties p = new Properties();
    p.load(new FileInputStream( myPropertiesFile ));
    String value;
    value = p.getProperty( PROPERTY_colorBitType );
    if( value != null ) colorBitType  = Integer.parseInt( value );
    value = p.getProperty( PROPERTY_endianSEGYHdr );
    if( value != null ) {
      segyAttr.endianSEGYHdr = value.compareTo("LITTLE_ENDIAN") == 0 ?
              java.nio.ByteOrder.LITTLE_ENDIAN : java.nio.ByteOrder.BIG_ENDIAN;
    }
    value = p.getProperty( PROPERTY_endianSEGYData );
    if( value != null ) {
      segyAttr.endianSEGYData = value.compareTo("LITTLE_ENDIAN") == 0 ?
              java.nio.ByteOrder.LITTLE_ENDIAN : java.nio.ByteOrder.BIG_ENDIAN;
    }
    value = p.getProperty( PROPERTY_endianSUHdr );
    if( value != null ) {
      segyAttr.endianSUHdr = value.compareTo("LITTLE_ENDIAN") == 0 ?
              java.nio.ByteOrder.LITTLE_ENDIAN : java.nio.ByteOrder.BIG_ENDIAN;
    }
    value = p.getProperty( PROPERTY_endianSUData );
    if( value != null ) {
      segyAttr.endianSUData = value.compareTo("LITTLE_ENDIAN") == 0 ?
              java.nio.ByteOrder.LITTLE_ENDIAN : java.nio.ByteOrder.BIG_ENDIAN;
    }
    value = p.getProperty( PROPERTY_isCustomMapSEGY );
    if( value != null ) segyAttr.isCustomMap  = Boolean.parseBoolean( value );
    value = p.getProperty( PROPERTY_hdrMapSEGY );
    if( value != null ) segyAttr.hdrMap  = Integer.parseInt( value );
    value = p.getProperty( PROPERTY_autoScaleSEGY );
    if( value != null ) segyAttr.autoScaleCoord  = Boolean.parseBoolean( value );
    value = p.getProperty( PROPERTY_windowLayout );
    if( value != null ) windowLayout  = Integer.parseInt( value );

    value = p.getProperty( PROPERTY_annotationOmitRepeating );
    if( value != null ) annAttr.omitRepeating  = Boolean.parseBoolean( value );
    value = p.getProperty( PROPERTY_annotationShowSequential );
    if( value != null ) annAttr.showSequential  = Boolean.parseBoolean( value );
    value = p.getProperty( PROPERTY_annotationFixedTraceLabelStep );
    if( value != null ) annAttr.fixedTraceLabelStep  = Boolean.parseBoolean( value );
    value = p.getProperty( PROPERTY_annotationTraceLabelStep );
    if( value != null ) annAttr.traceLabelStep  = Integer.parseInt( value );

    value = p.getProperty( PROPERTY_numTraces );
    if( value != null ) numTraces  = Integer.parseInt( value );
    value = p.getProperty( PROPERTY_maxTabWidth );
    if( value != null ) maxTabWidth = Integer.parseInt( value );
    value = p.getProperty( PROPERTY_showFilename );
    if( value != null ) showFilename = Boolean.parseBoolean( value );
    value = p.getProperty( PROPERTY_showActivePaneFeedback );
    if( value != null ) showActivePaneFeedback = Boolean.parseBoolean( value );
    
    propertiesDirectoryPath = p.getProperty( PROPERTY_propertyDir );
    dataDirectoryPath = p.getProperty( PROPERTY_dataDir );
    screenDumpDirectoryPath = p.getProperty( PROPERTY_screenDumpDir );
    velFieldDirectoryPath = p.getProperty( PROPERTY_velFieldDir );
    wellPathDirectoryPath = p.getProperty( PROPERTY_wellPathDir );
  }
  public void write() throws IOException {
    Properties p = new Properties();
    if( dataDirectoryPath != null ) p.setProperty( PROPERTY_dataDir, dataDirectoryPath );
    if( screenDumpDirectoryPath != null ) p.setProperty( PROPERTY_screenDumpDir, screenDumpDirectoryPath );
    if( propertiesDirectoryPath != null ) p.setProperty( PROPERTY_propertyDir, propertiesDirectoryPath );
    if( velFieldDirectoryPath != null ) p.setProperty( PROPERTY_velFieldDir, velFieldDirectoryPath );
    if( wellPathDirectoryPath != null ) p.setProperty( PROPERTY_wellPathDir, wellPathDirectoryPath );
    p.setProperty( PROPERTY_colorBitType, colorBitType+"" );
    
    p.setProperty( PROPERTY_endianSEGYHdr, segyAttr.endianSEGYHdr.toString() );
    p.setProperty( PROPERTY_endianSEGYData, segyAttr.endianSEGYData.toString() );
    p.setProperty( PROPERTY_endianSUHdr, segyAttr.endianSUHdr.toString() );
    p.setProperty( PROPERTY_endianSUData, segyAttr.endianSUData.toString() );
    p.setProperty( PROPERTY_isCustomMapSEGY, segyAttr.isCustomMap+"" );
    p.setProperty( PROPERTY_hdrMapSEGY, segyAttr.hdrMap+"" );
    p.setProperty( PROPERTY_autoScaleSEGY, segyAttr.autoScaleCoord+"" );
    p.setProperty( PROPERTY_windowLayout, windowLayout+"" );
    p.setProperty( PROPERTY_annotationOmitRepeating, annAttr.omitRepeating+"" );
    p.setProperty( PROPERTY_annotationShowSequential, annAttr.showSequential+"" );
    p.setProperty( PROPERTY_annotationFixedTraceLabelStep, annAttr.fixedTraceLabelStep+"" );
    p.setProperty( PROPERTY_annotationTraceLabelStep, annAttr.traceLabelStep+"" );

    p.setProperty( PROPERTY_numTraces, numTraces+"" );
    p.setProperty( PROPERTY_maxTabWidth, maxTabWidth+"" );
    p.setProperty( PROPERTY_showFilename, showFilename+"" );
    p.setProperty( PROPERTY_showActivePaneFeedback, showActivePaneFeedback+"" );
    
    p.store(new FileOutputStream( myPropertiesFile ), null);
  }
}


