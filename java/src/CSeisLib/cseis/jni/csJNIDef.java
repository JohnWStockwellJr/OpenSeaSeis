/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.jni;

/**
 * Definitons that span across JNI interface between C++ and Java classes.
 * @author Bjorn Olofsson
 */
public class csJNIDef {
  /// Data domain types
  // Never change these numbers, for backward compatibility
  public static final int DOMAIN_UNKNOWN = -1;
  /// Time domain
  public static final int DOMAIN_XT = 110;
  /// Frequency-trace domain
  public static final int DOMAIN_FX = 111;
  /// Frequency-K domain
  public static final int DOMAIN_FK = 112;
  /// Time series-K domain
  public static final int DOMAIN_KT = 113;
  /// Depth domain
  public static final int DOMAIN_XD = 114;

  /// Vertical domain
  public static final int VERTICAL_DOMAIN_TIME  = 1001;
  public static final int VERTICAL_DOMAIN_DEPTH = 1002;
  public static final int VERTICAL_DOMAIN_FREQ  = 1003;

  public static final int SORT_NONE = 0;
  public static final int SORT_INCREASING = 1;
  public static final int SORT_DECREASING = 2;

  public static final int SIMPLE_SORT = 101;
  public static final int TREE_SORT   = 201;

  public static final int TYPE_UNKNOWN = 255;
//  public static final int TYPE_EMPTY   = 0;
  public static final int TYPE_INT     = 1;  // 32bit
  public static final int TYPE_FLOAT   = 2;  // 32bit
  public static final int TYPE_DOUBLE  = 3;  // 64bit
  public static final int TYPE_CHAR    = 4;  // 8bit
  public static final int TYPE_STRING  = 5;
  public static final int TYPE_LONG    = 6;  // 64bit integer
//  public static final int TYPE_SHORT   = 11;  // 16bit signed short
//  public static final int TYPE_USHORT  = 12;  // 16bit unsigned short
  
  public static int convertVerticalDomain_fromJNI( int jniDomain ) {
    switch( jniDomain ) {
      case VERTICAL_DOMAIN_TIME:
        return cseis.general.csUnits.DOMAIN_TIME;
      case VERTICAL_DOMAIN_DEPTH:
        return cseis.general.csUnits.DOMAIN_DEPTH;
      case VERTICAL_DOMAIN_FREQ:
        return cseis.general.csUnits.DOMAIN_FREQ;
      default:
        return cseis.general.csUnits.DOMAIN_TIME;
    }
  }
  public static int convertVerticalDomain_toJNI( int domain ) {
    switch( domain ) {
    case cseis.general.csUnits.DOMAIN_TIME:
      return VERTICAL_DOMAIN_TIME;
    case cseis.general.csUnits.DOMAIN_DEPTH:
      return VERTICAL_DOMAIN_DEPTH;
    case cseis.general.csUnits.DOMAIN_FREQ:
      return VERTICAL_DOMAIN_FREQ;
    default:
      return VERTICAL_DOMAIN_TIME;
    }
  }
}

