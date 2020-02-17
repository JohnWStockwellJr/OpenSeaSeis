/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.jni;

import cseis.seis.csISeismicReader;

/**
 * Utilities related to seismic readers.
 * @author Bjorn Olofsson
 */
public class csSeismicReaderUtils {
  public static csISeismicReader createSeismicReader( csISeismicReader refReader, String filename ) throws Exception {
    if( refReader instanceof csNativeSeismicReader ) {
      return new csNativeSeismicReader( filename );
    }
    else if( refReader instanceof csNativeSegyReader ) {
      return new csNativeSegyReader( filename );
    }
    else if( refReader instanceof csNativeSegdReader ) {
      return new csNativeSegdReader( filename );
    }
    else if( refReader instanceof csASCIIReader ) {
      return new csASCIIReader( filename );
    }
    else if( refReader instanceof csNativeRSFReader ) {
      return new csNativeRSFReader( filename );
    }
    else {
//    else if( refReader instanceof csVirtualSeismicReader ) {
//      return new csVirtualSeismicReader();
//    }
      return null;
    }
  }
}

