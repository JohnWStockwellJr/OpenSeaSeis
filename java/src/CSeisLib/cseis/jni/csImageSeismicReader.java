/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

package cseis.jni;

import cseis.seis.csHeader;
import cseis.seis.csHeaderDef;
import cseis.seis.csSeismicTrace;
import cseis.seis.csTraceBuffer;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Iterator;
import javax.imageio.ImageIO;
import javax.imageio.ImageReader;
import javax.imageio.stream.FileImageInputStream;

/**
 *
 */
public class csImageSeismicReader extends csVirtualSeismicReader {
  
  /**
   * 
   * @param filename
   * @throws IOException 
   */
  public csImageSeismicReader( String filename ) throws IOException {
    mySampleInt = 1.0f;
    myCurrentTraceIndex = 0;
    myPeekHeaderIndex = 0;
    mySelectionNotifier = null;
//    myVerticalDomain = verticalDomain;

    FileImageInputStream input = null;
    try {
      input = new FileImageInputStream( new File(filename) );
    }
    catch( FileNotFoundException e1) {
      e1.printStackTrace();
      throw e1;
    }
    catch (IOException e2) {
      e2.printStackTrace();
      throw e2;
    }
    Iterator<ImageReader> iter = ImageIO.getImageReadersByFormatName("JPG");
    if( iter.hasNext() ) {
      ImageReader reader = iter.next();
     
      reader.setInput( input );
      BufferedImage image = null;
      try {
        image = reader.read(0);
      }
      catch (IOException e) {
        e.printStackTrace();
        throw e;
      }
      int numTraces   = image.getWidth();
      int numSamples  = image.getHeight();
      int numHeaders  = 1;
      myTraceBuffer = new csTraceBuffer( numSamples, numHeaders );
      for( int itrc = 0; itrc < numTraces; itrc++ ) {
        csHeader[] headerValues = new csHeader[numHeaders];
        headerValues[0] = new csHeader(itrc+1);
        float[] imageDataArray = new float[numSamples];
        for( int isamp = 0; isamp < numSamples; isamp++ ) {
          imageDataArray[isamp] = (float)image.getRGB( itrc, isamp );
        }
        myTraceBuffer.addTrace( new csSeismicTrace( imageDataArray, headerValues, itrc+1 ) );
      }
      myTraceHeaderDef = new csHeaderDef[numHeaders];
      myTraceHeaderDef[0] = new csHeaderDef("trace","Trace number",csHeaderDef.TYPE_INT);
    }
  }
}

