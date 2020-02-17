/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seaview;

import cseis.general.csColorBarPanel;
import cseis.general.csColorMap;
import cseis.general.csUnits;
import cseis.jni.csNativeSegyReader;
import cseis.jni.csNativeSeismicReader;
import cseis.seis.*;
import cseis.seisdisp.csSeisDispProperties;
import cseis.seisdisp.csSeisDispSettings;
import cseis.seisdisp.csWellPathAttr;
import cseis.seisdisp.csWellPathOverlay;
import cseis.seisdisp.csWellPathOverlayAttr;

import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.imageio.ImageIO;

/**
 * Plot seismic view to external image file.<br>
 * This class provides mostly the same seismic display functionality as SeaView.
 * It is run in the background to create an external image file displaying seimic data.
 * Annotations of the vertical and horizontal axes are using nicer fonts than SeaView, and there are
 * more options to decorate color bars, add axis labels and a title.
 * @author 2007 Bjorn Olofsson
 */
public class PlotImage {
  public static final int COLORBAR_NONE   = 0;
  public static final int COLORBAR_RIGHT  = 5;
  public static final int COLORBAR_BOTTOM = 10;
  public static final int WIDTH_COLORBAR  = 40;

  private String myPropertiesFilename;
  private boolean myIsProperties;
  private csSeaViewProperties mySeaViewProperties;
//  private csSeisDispSettings myDispSettings;

  private csPlotImageSeisView mySeisView;
  private csISeismicReader myReader = null;
  private csTraceBuffer myTraceBuffer;
  private csHeaderDef[] myTraceHeaders;
  private String mySeismicFilename;
  private String myImageFilename;

  private int myMarginLeftRight;
  private int myMarginTopBottom;
  private int myFontSize;
  private int myLabelMargin = 6;  // Margin/pixels between annotation label and image
  private int myHeaderIndex = -1;
  private String myHeaderName = null;
  private int myNumDecimals = 1;
  private int myNumMaxDecimals = 3;
  private DecimalFormat myDecimalFormat;
  private int mySampleIntScalar;
  private float myMinTime_ms;
  private float myMaxTime_ms;
  private int myColorBarHeight;
  private int myColorBarWidth;

  public static String IMAGE_FORMAT_PNG = "PNG";
  public static String IMAGE_FORMAT_JPG = "JPEG";
  public static String IMAGE_FORMAT_GIF = "GIF";
  public static String IMAGE_FORMAT_BMP = "BMP";

  public static int DATE_NONE  = 0;
  public static int DATE_HOURS = 1;
  public static int DATE_DAYS  = 2;
  public static int DATE_DAY_MONTH = 3;
  public static int DATE_DAY_MONTH_YEAR = 4;
  //------------------------------------------------------------------------
  public PlotImage( String filenameIn, String filenameOut, String propertiesFilename ) throws Exception {
    myPropertiesFilename = propertiesFilename;
    myIsProperties = myPropertiesFilename != null ? true : false;
    mySeaViewProperties = new csSeaViewProperties();
    try {
      mySeaViewProperties.load();
    }
    catch( IOException ex ) {
    }
    mySeismicFilename = filenameIn;
    myImageFilename   = filenameOut;
    myReader    = null;
    myFontSize  = 11;
    mySampleIntScalar   = 1000;
    myMinTime_ms = 0;
    myMaxTime_ms = 0;
    myColorBarHeight = 0;
    myColorBarWidth  = 0;

    setNumDecimals( 0 );

    File file = new File( mySeismicFilename );
    if( !file.exists() ) {
      throw( new Exception("File not found: " + mySeismicFilename) );
    }
    file = null;
    try {
      if( mySeismicFilename.toLowerCase().endsWith("cseis") || mySeismicFilename.toLowerCase().endsWith("oseis") ) {
        myReader = new csNativeSeismicReader( mySeismicFilename );
      }
      else {
        myReader = new csNativeSegyReader( mySeismicFilename, csNativeSegyReader.HDR_MAP_STANDARD, 200, false, false, false );
      }
    }
    catch( IOException e ) {
      throw( new Exception("Error opening seismic data input file: " + mySeismicFilename + ". System message: " + e.getMessage() ) );
    }

    try {
      file = new File( myImageFilename );
    }
    catch( Exception e ) {
      throw( new Exception("int Error opening image output file: " + myImageFilename ) );
    }
    file = null;

    myMarginLeftRight = 40;
    myMarginTopBottom = 20;
  }
  //------------------------------------------------------------------------
  //
  //
  private void setProperties( String propertiesFilename, csSeisDispSettings settings ) throws Exception {
    myIsProperties = true;
    myPropertiesFilename = propertiesFilename;
    csSeisDispProperties sp = new csSeisDispProperties( myPropertiesFilename, mySeaViewProperties.retrieveColorMaps() );
    try {
      sp.load( settings );
//      settings.dump();
    }
    catch( IOException e ) {
      throw( new Exception("Error while opening/reading properties file " + propertiesFilename +
        ". System message:\n" + e.getMessage()) );
    }
  }
  //------------------------------------------------------------------------
  //
  //
  private boolean readData() {
    int numHeaders = myReader.numHeaders();
    myTraceBuffer = new csTraceBuffer( myReader.numSamples(), numHeaders );
    if( myReader.isFrequencyDomain() ) mySampleIntScalar = 1;

    int numSamples = myReader.numSamples();
    csSeismicTrace trace = new csSeismicTrace( numSamples, numHeaders );
    int currentTraceIndex = 0;
    try {
      myTraceBuffer.clear();
      while( myReader.getNextTrace( trace ) ) {
        trace.setOriginalTraceNumber( currentTraceIndex+1 );
        myTraceBuffer.addTrace( trace );
        trace = new csSeismicTrace( numSamples, numHeaders );
        currentTraceIndex += 1;
      }
    }
    catch( Exception e ) {
      System.out.println("Error occurred when reading next trace from input file. System message:\n\n" + e.getMessage() );
      return false;
    }

    myTraceHeaders = new csHeaderDef[numHeaders];
    for( int i = 0; i < numHeaders; i++ ) {
      myTraceHeaders[i] = new csHeaderDef( myReader.headerName(i), myReader.headerDesc(i), myReader.headerType(i) );
    }
    myMaxTime_ms = myReader.sampleInt() * ( myReader.numSamples() -1 );

    return true;
  }
  public void closeFiles() {
    if( myReader != null ) myReader.closeFile();
    myReader    = null;
  }
  //------------------------------------------------------------------------
  //
  //
  public void convert_db() {
    float epsilon = myTraceBuffer.samples( 0 )[0] / 1e6f;
    for( int itrc = 0; itrc < myTraceBuffer.numTraces(); itrc++ ) {
      float [] samples = myTraceBuffer.samples( itrc );
      for( int isamp = 0; isamp < myTraceBuffer.numSamples(); isamp++ ) {
        samples[isamp] = (float)( 10.0 * Math.log10( Math.abs( samples[isamp] + epsilon ) ) );
      }
    }
  }
  //------------------------------------------------------------------------
  //
  //
  public void generatePlot( int width, int height, boolean evenTraceSpacing, int fontsize, String headerName,
      float minTime, float maxTime, boolean convert2db, int colorBarOption, int colorBarSize, int colorBarAnnOption, int colorBarAnnSize,
      int maxDecimals, int hdrDateType, csWellPathOverlayAttr wellPathOverlayAttr, ArrayList<csWellPathAttr> wellPathAttrList ) throws Exception {
    myNumMaxDecimals = maxDecimals;
    myFontSize = fontsize;

    if( !readData() ) return;

    if( convert2db ) convert_db();

    if( minTime > 0 ) {
      if( minTime >= myMaxTime_ms ) {
        throw( new Exception("Inconsistent minimum time specified: " + minTime + " > maximum time in input data: " + myMaxTime_ms) );
      }
      myMinTime_ms = minTime;
      // Change minimum time to full multiple of sample interval, to match current constraint of Seaview
      int numTimes = (int)( minTime / myReader.sampleInt() + 0.51 );
      myMinTime_ms = numTimes * myReader.sampleInt();
    }
    myMaxTime_ms = myReader.sampleInt()*(myReader.numSamples()-1);
    if( maxTime > 0 ) {
      if( maxTime <= myMinTime_ms ) {
        throw( new Exception("Inconsistent maximum time specified: " + maxTime + " > specified minimum time: " + myMinTime_ms) );
      }
      if( maxTime < myMaxTime_ms ) myMaxTime_ms = maxTime;
      // Change maximum time to full multiple of sample interval, to match current constraint of Seaview
      int numTimes = (int)( maxTime / myReader.sampleInt() + 0.51 );
      myMaxTime_ms = numTimes * myReader.sampleInt();
    }

    String imageFormat = IMAGE_FORMAT_JPG;
    if( myImageFilename.endsWith( ".JPG" ) || myImageFilename.endsWith( ".jpg" ) || myImageFilename.endsWith( ".JPEG" ) || myImageFilename.endsWith( ".jpeg" ) ) {
      imageFormat = IMAGE_FORMAT_JPG;
    }
    else if( myImageFilename.endsWith( ".PNG" ) || myImageFilename.endsWith( ".png" ) ) {
      imageFormat = IMAGE_FORMAT_PNG;
    }
    else if( myImageFilename.endsWith( ".GIF" ) || myImageFilename.endsWith( ".gif" ) ) {
      imageFormat = IMAGE_FORMAT_GIF;
    }
    else if( myImageFilename.endsWith( ".BMP" ) || myImageFilename.endsWith( ".bmp" ) ) {
      imageFormat = IMAGE_FORMAT_BMP;
    }
    else if( myImageFilename.contains( "." ) ) {
      throw( new Exception("Image output file does not have recognized extension: " + myImageFilename ) );
    }
    else {
      myImageFilename += ".jpg";
    }

    // Check annotation header, if specified
    if( headerName != null ) {
      myHeaderName = headerName;
      myHeaderIndex = -1;
      for( int i = 0; i < myTraceHeaders.length; i++ ) {
        if( headerName.compareTo( myTraceHeaders[i].name ) == 0 ) {
          myHeaderIndex = i;
          break;
        }
      }
      if( myHeaderIndex < 0 ) {
        throw( new Exception("Trace header not defined in input file: " + headerName ) );
      }
    }

    if( myHeaderIndex >= 0 && hdrDateType != DATE_NONE ) {
      csHeader hdrValue1 = myTraceBuffer.headerValues( 0 )[myHeaderIndex];
      int time0 = hdrValue1.intValue();
      for( int i = 0; i < myTraceBuffer.numTraces(); i++ ) {
        csHeader hdrValue2 = myTraceBuffer.headerValues( i )[myHeaderIndex];
        if( hdrDateType == DATE_HOURS || hdrDateType == DATE_DAYS ) {
          double period = hdrValue2.intValue() - time0;
          hdrValue2.setValue( 0.0 );
          period /= 3600.0;
          if( hdrDateType == DATE_DAYS ) {
            myHeaderName = "Days";
            period /= 24.0;
          }
          else {
            myHeaderName = "Hours";
          }
          hdrValue2.setValue( (long)(period*10.0)/10.0 );
        }
        else {
          SimpleDateFormat dateFormat = null;
          if( hdrDateType == DATE_DAY_MONTH ) {
            dateFormat = new SimpleDateFormat("dd/MM");
            myHeaderName = "Date(d/m)";
          }
          else {
            dateFormat = new SimpleDateFormat("dd/MM/yyyy");
            myHeaderName = "Date(d/m/y)";
          }
          java.util.Date date2 = new java.util.Date((long)hdrValue2.intValue()*1000L);
          hdrValue2.setValue( dateFormat.format( date2 ));
        }
      }
    }

    mySeisView = new csPlotImageSeisView( null, myTraceBuffer, myReader.sampleInt() );
    
    csSeisDispSettings settings = mySeisView.getDispSettings();

    if( myIsProperties ) {
      setProperties( myPropertiesFilename, settings );
      if( settings.title != null && settings.title.compareTo("") != 0 ) {
        exchangeTitleHeaders( settings );
      }
      else {
        settings.title = "";
      }
      mySeisView.updateDispSettings( settings );
    }
    else {
    }

    recomputeDecimals( mySeisView.getDispSettings().timeLineMinorInc );

    csPlotImageColorBar colorBarPanel = null;
    if( colorBarOption != PlotImage.COLORBAR_NONE ) {
      if( colorBarOption == PlotImage.COLORBAR_BOTTOM ) {
         colorBarPanel = new csPlotImageColorBar( new csColorMap(), csColorBarPanel.ORIENT_HORIZONTAL, width, colorBarSize, colorBarAnnOption );
         myColorBarHeight = colorBarSize+colorBarAnnSize;
      }
      else {
         colorBarPanel = new csPlotImageColorBar( new csColorMap(), csColorBarPanel.ORIENT_VERTICAL, colorBarSize, height, colorBarAnnOption );
         myColorBarWidth = colorBarSize+colorBarAnnSize;
      }
      if( settings.isVIDisplay ) {
        colorBarPanel.setColorMap( settings.viColorMap );
      }
      else if( settings.isVariableColor ) {
        colorBarPanel.setColorMap( settings.wiggleColorMap );
      }
      else {
        colorBarPanel.setColorMap( new csColorMap( Color.black ) );
      }
      colorBarPanel.setFontSize(myFontSize);
      colorBarPanel.setMinMax( settings.minValue, settings.maxValue );
    }

    if( wellPathAttrList.size()> 0 ) {
      csWellPathOverlay wellPathOverlay = new csWellPathOverlay( mySeisView );
      int hdrIndex_x = getTraceHeaderIndex( wellPathOverlayAttr.hdrName_binx );
      int hdrIndex_y = getTraceHeaderIndex( wellPathOverlayAttr.hdrName_biny );
      wellPathOverlay.updateHeaders( hdrIndex_x, hdrIndex_y );
      wellPathOverlay.update( wellPathOverlayAttr, wellPathAttrList );
      mySeisView.addOverlay( wellPathOverlay );
    }
    
    double minorInc = mySeisView.getDispSettings().timeLineMinorInc;
    minorInc = ( (int)( minorInc / (double)mySeisView.getSampleInt() + 0.51 ) ) * (double)mySeisView.getSampleInt();
    recomputeDecimals( minorInc );

    int newWidth = mySeisView.set( width, height, myMarginTopBottom, myMarginLeftRight, evenTraceSpacing,
      myMinTime_ms, myMaxTime_ms, myColorBarWidth, myColorBarHeight );
    if( newWidth != width ) {
//      System.out.println("+++++++++++ New width: " + width + " " + newWidth );
      width = newWidth;
//      if( colorBarOption == PlotImage.COLORBAR_RIGHT ) {
//        width += WIDTH_COLORBAR;
//      }
    }

    height += myColorBarHeight;
    width  += myColorBarWidth;

/*    if( colorBarOption == PlotImage.COLORBAR_BOTTOM ) {
    }
    else {
    } */

    BufferedImage bitmap = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
    Graphics2D graphics = bitmap.createGraphics();


    int axisLabelWidth = (int)( myMarginLeftRight / 1.5 );
    graphics.translate( axisLabelWidth, 0 );

    mySeisView.paintImage( graphics );

    graphics.setStroke( new BasicStroke(2) );
    graphics.setColor( Color.white );
    graphics.fillRect( 0, 0, myMarginLeftRight-2, height );
    graphics.fillRect( width-myMarginLeftRight+4-myColorBarWidth, 0, myMarginLeftRight-4+myColorBarWidth, height );

    paintAnnotation( graphics, width-myColorBarWidth, height-myColorBarHeight, settings.title, hdrDateType );
    if( colorBarOption != PlotImage.COLORBAR_NONE ) {
      if( colorBarOption == PlotImage.COLORBAR_BOTTOM ) {
        graphics.translate( 0, height-myColorBarHeight );
//        graphics.translate( 10, 0 );
      }
      else {
        graphics.translate( width-myColorBarWidth-(myMarginLeftRight-2), 0 );
      }
      colorBarPanel.paintComponent( graphics );
      if( colorBarOption == PlotImage.COLORBAR_BOTTOM ) {
        graphics.translate( 0, myColorBarHeight-height );
      }
      else {
        graphics.translate( myColorBarWidth-width+(myMarginLeftRight-2), 0 );
      }
    }

    // Add time/frequency axis label
    graphics.setColor( Color.white );
    graphics.fillRect( -axisLabelWidth, 0, axisLabelWidth, height );
    graphics.setColor( Color.black );
    String label = "Time (seconds)";
    if( settings.titleVertAxis.compareToIgnoreCase("DEFAULT") != 0 ) {
      label = settings.titleVertAxis;
    }
    else if( myReader.isFrequencyDomain() ) {
      label = "Frequency [Hz]";
    }
    else if( myReader.verticalDomain() == csUnits.DOMAIN_DEPTH ) {
      label = "Depth (kilometers)";
    }
    FontMetrics metrics = graphics.getFontMetrics( graphics.getFont() );
    graphics.rotate( -Math.PI/2 );
    int labelWidth  = metrics.stringWidth( label );
    graphics.drawString( label, -height/2 - labelWidth/2, -3 );

    try {
      ImageIO.write( bitmap, imageFormat, new File( myImageFilename ) );
    }
    catch( IOException ex ) {
      ex.printStackTrace();
      Logger.getLogger( PlotImage.class.getName() ).log( Level.SEVERE, null, ex );
    }
  }
  private int getTraceHeaderIndex( String headerName ) throws Exception {
    for( int i = 0; i < myTraceHeaders.length; i++ ) {
      if( headerName.compareTo( myTraceHeaders[i].name ) == 0 ) {
        return i;
      }
    }
    throw( new Exception("Trace header not defined in input file: " + headerName ) );
  }
  //---------------------------------------------------
  //
  private void paintAnnotation( Graphics2D g2, int width, int height, String title, int hdrDateType ) {
    g2.setColor( Color.black );
    g2.drawRect( myMarginLeftRight, myMarginTopBottom, width-2*myMarginLeftRight, height-2*myMarginTopBottom );
    g2.setFont( new Font("SansSerif",Font.PLAIN,myFontSize) );
    g2.setStroke( new BasicStroke(2) );
    FontMetrics metrics = g2.getFontMetrics( g2.getFont() );
    int labelHeightHalf = metrics.getHeight() / 2 - 3;

    // Title
    int titleWidth = metrics.stringWidth( title );
    g2.setFont( new Font("SansSerif",Font.BOLD,myFontSize) );
    g2.drawString( title, width/2-titleWidth/2, height-myMarginTopBottom+2*labelHeightHalf+myLabelMargin );
    g2.setFont( new Font("SansSerif",Font.PLAIN,myFontSize) );

    // Time axis annotation
    double minorInc = mySeisView.getDispSettings().timeLineMinorInc;
    minorInc = ( (int)( minorInc / (double)mySeisView.getSampleInt() + 0.51 ) ) * (double)mySeisView.getSampleInt();
    if( minorInc == 0 ) minorInc = (double)mySeisView.getSampleInt();

    recomputeDecimals( minorInc );
    int numTimes = (int)( (myMaxTime_ms-myMinTime_ms) / minorInc + 0.5 );
    float zoomVert = mySeisView.getDispSettings().zoomVert;
    double pixelsPerMS = zoomVert / (double)mySeisView.getSampleInt();
    int numPreSteps = (int)( myMinTime_ms / minorInc );
    double firstTimeToAnnotate_ms = (numPreSteps+1) * minorInc;
    for( int i = 0; i < numTimes; i++ ) {
      double time_ms = firstTimeToAnnotate_ms + minorInc * i;
      int ypixel = (int)( (time_ms-myMinTime_ms) * pixelsPerMS ) + myMarginTopBottom;
      String label = myDecimalFormat.format( (time_ms)/(double)mySampleIntScalar );
      int labelWidth  = metrics.stringWidth( label );
      g2.drawString( label, myMarginLeftRight-labelWidth-myLabelMargin, ypixel+labelHeightHalf );
      // g2.drawString( label, width-myMarginLeftRight+myLabelMargin, ypixel+labelHeightHalf );
    }

    // Trace axis annotation
    if( myHeaderIndex >= 0 ) {
      csHeader hdrValue1 = myTraceBuffer.headerValues( 0 )[myHeaderIndex];
      csHeader hdrValue2 = myTraceBuffer.headerValues( myTraceBuffer.numTraces()-1 )[myHeaderIndex];
      int labelWidthMax  = metrics.stringWidth( hdrValue1.toString() );
      labelWidthMax  = Math.max( labelWidthMax, metrics.stringWidth( hdrValue2.toString() ) );
      float numFit = (float)( width - 2*myMarginLeftRight ) / (float)labelWidthMax;
      if( numFit <= 0 ) numFit = 2;
      int traceStep = (int)( 2.0f*myTraceBuffer.numTraces() / numFit );
      if( traceStep >= myTraceBuffer.numTraces() ) traceStep = (int)( myTraceBuffer.numTraces() / 3 );
      if( traceStep <= 0 ) traceStep = 1;
      String label;
      for( int itrc = 0; itrc < myTraceBuffer.numTraces(); itrc += traceStep ) {
        csHeader[] hdrValues = myTraceBuffer.headerValues( itrc );
        float xView = mySeisView.xModel2View( itrc );
        label = hdrValues[myHeaderIndex].toString();
        int labelWidth  = metrics.stringWidth( label );
        if( itrc == 0 ) {
          g2.drawString( myHeaderName, myLabelMargin, myLabelMargin+labelHeightHalf*2 );
        }
        else {
          g2.drawString( label, xView-labelWidth/2, myLabelMargin+labelHeightHalf*2 );
        }
        g2.drawLine( (int)xView, myMarginTopBottom-3, (int)xView, myMarginTopBottom );
      }
    }
  }
  private void exchangeTitleHeaders( csSeisDispSettings s ) {
    String titleIn  = s.title;
    String titleOut = s.title;
    int counterIn  = 0;
    int counterOut = 0;
    titleOut.replace( '_', ' ' );
    while( counterIn < titleIn.length() ) {
      if( titleIn.charAt( counterIn ) == '%' ) {
        int pos1 = counterIn;
        counterIn += 1;
        while( counterIn < titleIn.length() && titleIn.charAt( counterIn ) != '%') {
          counterIn += 1;
        }
        if( counterIn == titleIn.length() ) return;
        String headerName = titleIn.substring( pos1+1, counterIn );
        for( int ihdr = 0; ihdr < myTraceHeaders.length; ihdr++ ) {
          if( myTraceHeaders[ihdr].name.compareTo( headerName ) == 0 ) {
            String newText = myTraceBuffer.headerValues( 0)[ihdr].toString();
            titleOut = titleOut.substring( 0, counterOut ) + newText;
            counterOut += newText.length() - 1;
            if( counterIn < titleIn.length()-1 ) {
              titleOut += titleIn.substring( counterIn+1 );
            }
            break;
          }
        }
      }
      counterIn  += 1;
      counterOut += 1;
    }
    s.title = titleOut;
  }
  private void setNumDecimals( int numDecimals ) {
    myNumDecimals = numDecimals;
    String s = new String( "0" );
    if( myNumDecimals > 0 ) {
      s += ".0";
    }
    for( int i = 1; i < myNumDecimals; i++ ) {
      s += "0";
    }
    myDecimalFormat = new DecimalFormat( s );
  }
  private void recomputeDecimals( double minorIncMS ) {
    long number = (long)( 1000*1000/(double)mySampleIntScalar*(double)minorIncMS );
    int minDecimal = 1;
    while( number != 0 && number % 10 == 0 ) {
      number /= 10;
      minDecimal += 1;
    }
    int numDecimals = 7 - minDecimal;
    if( numDecimals > myNumMaxDecimals ) numDecimals = myNumMaxDecimals;

    if( numDecimals != myNumDecimals ) setNumDecimals( numDecimals );

    number = (long)( (1.0/(double)mySampleIntScalar)*(double)myMaxTime_ms );
    int counter = 0;
    while( number != 0 ) {
      number /= 10;
      counter += 1;
    }
    if( counter > 0 ) counter -= 1;
    myMarginLeftRight = 15 + (myNumDecimals+counter)*(myFontSize);
    if( myMarginLeftRight < 2*myFontSize ) myMarginLeftRight = (int)( 2*myFontSize );
    myMarginTopBottom = 12 + myFontSize;
  }

  public static void help() {
    System.out.println("PlotImage -f <seismicFilename> -o <imageFilename> [options]");
    System.out.println(" -f <seismicFilename>:   Seismic input file ");
    System.out.println(" -o <imageFilename>:     Image output file");
    System.out.println("    Options:");
    System.out.println(" -p <propertiesFilename>:  Display properties file");
    System.out.println(" -height <height>:         Height of image in pixels");
    System.out.println(" -width <width>:           Width of image on pixels");
    System.out.println(" -header <headerName>:     Name of trace header to use for trace annotation");
    System.out.println(" -fontsize <fontsize>:     Font size in pixels");
    System.out.println(" -trace_spacing [auto|even]:  Automatic trace spacing, or integer number of pixels per trace (even)");
    System.out.println("                              For 'even' trace spacing, image width may be adjusted");
    System.out.println(" -min_time <minTime>:      Minimum time to display [ms] (or [Hz] for f-x plot)");
    System.out.println(" -max_time <maxTime>:      Maximum time to display [ms] (or [Hz] for f-x plot)");
    System.out.println(" -db [yes,no]:             Convert to db, yes/no");
    System.out.println(" -vert_title <title>:             Vertical axis title");
    System.out.println(" -color_bar <bottom|right><size> : Display color bar");
    System.out.println(" -color_bar <simple|fancy><size> : Color bar annotation");
    System.out.println(" -hdr_date_type [h|d|m|y]:        Convert trace header value to date. h: Hours, d: Julian day, m: Day/Month, Day/Month/Year");
    System.out.println(" -well_param <lineSize>,<maxDistance>,<transparency>,<hdrX>,<hdrY>,<plotShadow[yes|no]>,<shadowDepth> : Well path plot parameters");
    System.out.println("    lineSize: Line size [pixels], maxDistance: Max distance from seismic to plot [m], transparency: 0: Min, 100: Max");
    System.out.println("    hdrX,hdrY: Trace header names/XY coordinates, plotShadow: Plot well 'shadow'?, shadowDepth: 0: Min, 100: Max");
    System.out.println(" -well_path <filename>,<color> : Well path plot, file name and plot color. Can specify multiple well paths");
  }
//--------------------------------------------------------------------------------
//
  public static void main( String[] args ) {
    String path = System.getProperty( "java.library.path" );
    System.setProperty("java.awt.headless", "true");
    int counterArg = 0;
    try {
      // Create CSEIS configuration folder, if it doesn't exist yet
      cseis.general.csAbstractPreference.setCseisDirectory();
    }
    catch( Exception e ) {
      System.err.println( "FATAL ERROR occurred in PlotImage. System message: " + e.getMessage() );
      System.exit( 0 );
    }
  /*
    try {
      File tempErr = File.createTempFile("stderr", ".log");
      File tempOut = File.createTempFile("stdut", ".log");
//      tempErr.deleteOnExit();
//      tempOut.deleteOnExit();
      FileOutputStream err = new FileOutputStream( tempErr );
      FileOutputStream std = new FileOutputStream( tempOut );
      System.setErr(new PrintStream(err));
      System.setOut(new PrintStream(std));
    }
    catch( FileNotFoundException e ) {
      e.printStackTrace();
      return;
    }
    catch( IOException e ) {
      e.printStackTrace();
      return;
    }
*/

//    for( int i = 0; i < args.length; i++ ) {
//      System.out.println("OPTION #" + i + ":  " + args[i]);
//    }

    String filenameIn  = null;
    String filenameOut = null;;
    String propertiesFilename = null;
    String headerName = null;
    int width  = 2000;
    int height = 800;
    boolean evenTraceSpacing = true;
    int fontsize = 11;
    float minTime = 0;
    float maxTime = 0;
    boolean convert2db = false;
    int colorBarOption = PlotImage.COLORBAR_NONE;
    int colorBarSize = 0;
    int colorBarAnnOption = csColorBarPanel.ANNOTATION_SIMPLE;
    int colorBarAnnSize = 0;
    int maxDecimals = 3;
    int hdrDateType = PlotImage.DATE_NONE;
    csWellPathOverlayAttr wellPathOverlayAttr = new csWellPathOverlayAttr();
    ArrayList<csWellPathAttr> wellPathAttrList = new ArrayList<csWellPathAttr>();

    while( counterArg < args.length ) {
      char letter = args[counterArg].charAt(0);
      if( letter == '-' ) {
        if( args[counterArg].length() <= 1 ) {
          System.out.println("Error in command line input");
          help();
          System.exit(0);
        }
        String word = args[counterArg].substring( 1 );
        if( counterArg >= args.length-1 ) {
          System.out.println("Error in command line input.");
          help();
          System.exit(0);
        }
        if( word.compareTo( "fontsize" ) == 0 ) {
          counterArg += 1;
          fontsize = Integer.parseInt( args[counterArg] );
        }
        else if( word.compareTo( "width" ) == 0 ) {
          counterArg += 1;
          width = Integer.parseInt( args[counterArg] );
        }
        else if( word.compareTo( "height" ) == 0 ) {
          counterArg += 1;
          height = Integer.parseInt( args[counterArg] );
        }
        else if( word.compareTo( "trace_spacing" ) == 0 ) {
          counterArg += 1;
          if( args[counterArg].compareTo( "even" ) == 0 ) {
            evenTraceSpacing = true;
          }
          else if( args[counterArg].compareTo( "auto" ) == 0 ) {
            evenTraceSpacing = false;
          }
          else {
            System.out.println("Unknown value for option -trace_spacing: " + args[counterArg]);
            help();
            System.exit(0);
          }
        }
        else if( word.compareTo( "db" ) == 0 ) {
          counterArg += 1;
          if( args[counterArg].compareTo( "yes" ) == 0 ) {
            convert2db = true;
          }
          else if( args[counterArg].compareTo( "no" ) == 0 ) {
            convert2db = false;
          }
          else {
            System.out.println("Unknown value for option -db: " + args[counterArg]);
            help();
            System.exit(0);
          }
        }
        else if( word.compareTo( "color_bar" ) == 0 ) {
          counterArg += 1;
          if( args[counterArg].substring( 0, 6 ).compareTo( "bottom" ) == 0 ) {
            colorBarOption = PlotImage.COLORBAR_BOTTOM;
            colorBarSize   = Integer.parseInt( args[counterArg].substring( 6 ) );
          }
          else if( args[counterArg].substring( 0, 5 ).compareTo( "right" ) == 0 ) {
            colorBarOption = PlotImage.COLORBAR_RIGHT;
            colorBarSize   = Integer.parseInt( args[counterArg].substring( 5 ) );
          }
          else {
            System.out.println("Unknown value for option -color_bar: " + args[counterArg] + " " + args[counterArg].substring( 0, 5 ));
            help();
            System.exit(0);
          }
        }
        else if( word.compareTo( "color_bar_ann" ) == 0 ) {
          counterArg += 1;
          if( args[counterArg].substring( 0, 6 ).compareTo( "simple" ) == 0 ) {
            colorBarAnnOption = csColorBarPanel.ANNOTATION_SIMPLE;
            colorBarAnnSize   = Integer.parseInt( args[counterArg].substring( 6 ) );
          }
          else if( args[counterArg].substring( 0, 5 ).compareTo( "fancy" ) == 0 ) {
            colorBarAnnOption = csColorBarPanel.ANNOTATION_ADVANCED;
            colorBarAnnSize   = Integer.parseInt( args[counterArg].substring( 5 ) );
          }
          else {
            System.out.println("Unknown value for option -color_bar_ann: " + args[counterArg] + " " + args[counterArg].substring( 0, 5 ));
            help();
            System.exit(0);
          }
        }
        else if( word.compareTo( "well_path" ) == 0 ) {
          counterArg += 1;
          String[] wellPathParamList = args[counterArg].split(",");
          if( wellPathParamList.length < 2 ) {
            System.out.println("Too few parameters provided for well path" + args[counterArg] );
            help();
            System.exit(0);
          }
          csWellPathAttr attr = new csWellPathAttr();
          wellPathAttrList.add( attr );
          int yp = 0;
          attr.showWell = true;
          attr.filename = wellPathParamList[yp++];
          attr.path     = attr.filename;
          String colorStr = wellPathParamList[yp++];
          if( colorStr.compareTo("red") == 0 ) attr.color = Color.red;
          else if( colorStr.compareTo("red") == 0 ) attr.color = Color.blue;
          else if( colorStr.compareTo("yellow") == 0 ) attr.color = Color.yellow;
          else if( colorStr.compareTo("orange") == 0 ) attr.color = Color.orange;
          else if( colorStr.compareTo("black") == 0 ) attr.color = Color.black;
          else if( colorStr.compareTo("white") == 0 ) attr.color = Color.white;
          else if( colorStr.compareTo("gray") == 0 ) attr.color = Color.gray;
          else if( colorStr.compareTo("green") == 0 ) attr.color = Color.green;
          else if( colorStr.compareTo("blue") == 0 ) attr.color = Color.blue;
          if( !csWellPathDialog.readFile( null, attr ) ) { System.out.println("Cannot read well path ASCII file " + attr.filename); System.exit(0); }
        }
        else if( word.compareTo( "well_path_param" ) == 0 ) {
          counterArg += 1;
          String[] wellPathParamList = args[counterArg].split(",");
          if( wellPathParamList.length < 7 ) {
            System.out.println("Too few parameters provided for well_path_param" + args[counterArg] );
            help();
            System.exit(0);
          }
          int yp = 0;
          wellPathOverlayAttr.lineSize = Integer.parseInt( wellPathParamList[yp++] );
          wellPathOverlayAttr.maxDisplayDistance = Integer.parseInt( wellPathParamList[yp++] );
          wellPathOverlayAttr.transparency = ( Integer.parseInt( wellPathParamList[yp++] ) * (csWellPathOverlayAttr.MAX_TRANSPARENCY - csWellPathOverlayAttr.MIN_TRANSPARENCY) ) / 100 + csWellPathOverlayAttr.MIN_TRANSPARENCY;
          wellPathOverlayAttr.hdrName_binx = wellPathParamList[yp++];
          wellPathOverlayAttr.hdrName_biny = wellPathParamList[yp++];
          String showShadowStr = wellPathParamList[yp++];
          if( showShadowStr.compareTo("yes") == 0 ) wellPathOverlayAttr.showShadow = true;
          if( showShadowStr.compareTo("no") == 0 ) wellPathOverlayAttr.showShadow = false;
          wellPathOverlayAttr.shadowDepth = ( Integer.parseInt( wellPathParamList[yp++] ) * (csWellPathOverlayAttr.MAX_SHADOW_DEPTH- csWellPathOverlayAttr.MIN_SHADOW_DEPTH) ) / 100 + csWellPathOverlayAttr.MIN_SHADOW_DEPTH;
        }
        else if( word.compareTo( "min_time" ) == 0 ) {
          counterArg += 1;
          minTime = Float.parseFloat( args[counterArg] );
        }
        else if( word.compareTo( "max_time" ) == 0 ) {
          counterArg += 1;
          maxTime = Float.parseFloat( args[counterArg] );
        }
        else if( word.compareTo( "o" ) == 0 ) {
          counterArg += 1;
          filenameOut = args[counterArg];
        }
        else if( word.compareTo( "f" ) == 0 ) {
          counterArg += 1;
          filenameIn = args[counterArg];
        }
        else if( word.compareTo( "p" ) == 0 ) {
          counterArg += 1;
          propertiesFilename = args[counterArg];
        }
        else if( word.compareTo( "header" ) == 0 ) {
          counterArg += 1;
          headerName = args[counterArg];
        }
        else if( word.compareTo( "hdr_date_type" ) == 0 ) {
          counterArg += 1;
          hdrDateType = PlotImage.DATE_NONE;
          if( args[counterArg].compareTo("h") == 0 ) {
            hdrDateType = PlotImage.DATE_HOURS;
          }
          else if( args[counterArg].compareTo("d") == 0 ) {
            hdrDateType = PlotImage.DATE_DAYS;
          }
          else if( args[counterArg].compareTo("m") == 0 ) {
            hdrDateType = PlotImage.DATE_DAY_MONTH;
          }
          else if( args[counterArg].compareTo("y") == 0 ) {
            hdrDateType = PlotImage.DATE_DAY_MONTH_YEAR;
          }
          else {
            hdrDateType = PlotImage.DATE_DAY_MONTH_YEAR;
          }
        }
        else if( word.compareTo( "max_decimals" ) == 0 ) {
          counterArg += 1;
          maxDecimals = Integer.parseInt( args[counterArg] );
        }
      }
      counterArg += 1;
    } // END while
    if( filenameIn == null ) {
      System.out.println("Error in command line input: Need to specify seismic input filename");
      help();
      System.exit(0);
    }
    if( filenameOut == null ) {
      System.out.println("Error in command line input: Need to specify image output filename");
      help();
      System.exit(0);
    }
//    if( propertiesFilename == null ) {
//      System.out.println("Error in command line input: Need to specify properties filename");
//     System.exit(0);
//    }
    
    String libName = System.mapLibraryName( "csJNIlib" );
    try {
      System.load( path + "/" + libName );
    }
    catch( java.lang.UnsatisfiedLinkError e ) {
      System.err.println( e.toString() + "\n" + "java.library.path = " + System.getProperty( "java.library.path" ) + "\n" +
                          " - PlotImage will not run." );
      System.exit( -1 );
    }

    PlotImage plot = null;
    try {
//      System.out.println("new plotimage: " + filenameIn + "\n" + filenameOut + "\n" + propertiesFilename );
      plot = new PlotImage( filenameIn, filenameOut, propertiesFilename );
    }
    catch( Exception e ) {
      System.err.println( "FATAL ERROR occurred in PlotImage. System message: " + e.getMessage() );
      System.exit(-1);
    }
    
    try {
      plot.generatePlot( width, height, evenTraceSpacing, fontsize, headerName, minTime, maxTime,
          convert2db, colorBarOption, colorBarSize, colorBarAnnOption, colorBarAnnSize, maxDecimals, hdrDateType,
          wellPathOverlayAttr, wellPathAttrList );
      plot.closeFiles();
    }
    catch( Exception e ) {
      System.err.println( "FATAL ERROR occurred when generating plot. System message: " + e.getMessage() );
      System.err.println( e.getMessage() );
      System.exit(-1);
    }

    System.exit( 0 );
  }
}
