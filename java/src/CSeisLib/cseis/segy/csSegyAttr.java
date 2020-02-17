/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.segy;

import cseis.jni.csNativeSegyReader;

/**
 * Attributes defining how SEGY and SU files are read in.
 * @author 2013 Felipe Punto
 */
public class csSegyAttr {
  public java.nio.ByteOrder endianSEGYHdr;
  public java.nio.ByteOrder endianSEGYData;
  public java.nio.ByteOrder endianSUHdr;
  public java.nio.ByteOrder endianSUData;
  public int hdrMap;
  public boolean autoScaleCoord;
  public boolean isCustomMap;

  public csSegyAttr() {
    endianSEGYHdr  = java.nio.ByteOrder.BIG_ENDIAN;
    endianSEGYData = java.nio.ByteOrder.BIG_ENDIAN;
    endianSUHdr    = java.nio.ByteOrder.LITTLE_ENDIAN;
    endianSUData   = java.nio.ByteOrder.LITTLE_ENDIAN;;
    hdrMap         = csNativeSegyReader.HDR_MAP_STANDARD;
    autoScaleCoord = true;
    isCustomMap    = false;
  }
  public csSegyAttr( csSegyAttr attr ) {
    endianSEGYHdr  = attr.endianSEGYHdr;
    endianSEGYData = attr.endianSEGYData;
    endianSUHdr    = attr.endianSUHdr;
    endianSUData   = attr.endianSUData;
    hdrMap         = attr.hdrMap;
    autoScaleCoord = attr.autoScaleCoord;
    isCustomMap    = attr.isCustomMap;
  }
}

