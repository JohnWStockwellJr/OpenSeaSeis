/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.general;

/**
 * Provides unit definitions and static methods for unit conversions.
 * @author Bjorn Olofsson
 */
public class csUnits {
  public final static int DOMAIN_TIME  = 101;
  public final static int DOMAIN_DEPTH = 102;
  public final static int DOMAIN_FREQ  = 103;
  
  public final static int UNKNOWN_UNIT      = -1;
  public final static int METER             = 1;
  public final static int FEET              = 2;
  public final static int METERS_PER_SECOND = 3;
  public final static int MICROSECONDS_PER_FOOT = 4;
  public final static int ONE_PER_SECOND    = 5;
  public final static int SECONDS           = 6;
  public final static int MILLISECONDS      = 7;
  public final static int SECONDS_PER_METER = 8;
  public final static int KNOTS             = 9;
  public final static int NAUTICAL_MILE     = 10;
  public final static int INCH              = 11;
  public final static int KILOMETER         = 12;
  public final static int HERTZ             = 13;

  public final static String TEXT_UNKNOWN_UNIT         = "N/A";
  public final static String TEXT_METER                = "m";
  public final static String TEXT_FEET                 = "ft";
  public final static String TEXT_METERS_PER_SECOND    = "m/s";
  public final static String TEXT_MICROSECONDS_PER_FOOT[] = { "us/f", "us/ft" };
  public final static String TEXT_ONE_PER_SECOND    = "1/s";
  public final static String TEXT_SECONDS           = "s";
  public final static String TEXT_MILLISECONDS      = "ms";
  public final static String TEXT_SECONDS_PER_METER = "s/m";
  public final static String TEXT_KNOTS             = "knots";
  public final static String TEXT_NAUTICAL_MILE     = "nautical mile";
  public final static String TEXT_INCH              = "in";
  public final static String TEXT_KILOMETER         = "km";
  public final static String TEXT_HERTZ             = "Hz";
  
  public final static String TEXT_UNIT[] = {
      TEXT_UNKNOWN_UNIT,
      TEXT_METER,
      TEXT_FEET,
      TEXT_METERS_PER_SECOND,
      TEXT_MICROSECONDS_PER_FOOT[1],
      TEXT_ONE_PER_SECOND,
      TEXT_SECONDS,
      TEXT_MILLISECONDS,
      TEXT_SECONDS_PER_METER,
      TEXT_KNOTS,
      TEXT_NAUTICAL_MILE,
      TEXT_INCH,
      TEXT_KILOMETER,
      TEXT_HERTZ
  };

  public final static String TEXT_FEET_LONG  = "feet";
  public final static String TEXT_METER_LONG = "meter";
  public final static String TEXT_INCH_LONG = "inch";

  public final static String TEXT_FEET_SHORT  = "f";
  public final static double FEET_TO_METER = 0.3048;
  public final static double INCH_TO_CM    = 2.54;
  public final static double ONE_NAUTICAL_MILE_IN_KM = 1.852;     // [km]
  public final static double ONE_KNOT_IN_KM_PER_HOUR = ONE_NAUTICAL_MILE_IN_KM;

  public static int parseUnit( String unitString ) {
    if( unitString.equalsIgnoreCase(TEXT_METER) )
      return METER;
    if( unitString.equalsIgnoreCase(TEXT_METER_LONG) )
      return METER;
    if( unitString.equalsIgnoreCase(TEXT_FEET) )
      return FEET;
    if( unitString.equalsIgnoreCase(TEXT_FEET_LONG) )
      return FEET;
    if( unitString.equalsIgnoreCase(TEXT_FEET_SHORT) )
      return FEET;
    if( unitString.equalsIgnoreCase(TEXT_METERS_PER_SECOND) )
      return METERS_PER_SECOND;
    if( unitString.equalsIgnoreCase(TEXT_SECONDS) )
      return SECONDS;
    if( unitString.equalsIgnoreCase(TEXT_MILLISECONDS) )
      return MILLISECONDS;
    if( unitString.equalsIgnoreCase(TEXT_KILOMETER) )
      return KILOMETER;
    if( unitString.equalsIgnoreCase(TEXT_HERTZ) )
      return HERTZ;
    if( unitString.equalsIgnoreCase(TEXT_INCH) || unitString.equalsIgnoreCase(TEXT_INCH_LONG) )
      return INCH;
    for( int i = 0; i < TEXT_MICROSECONDS_PER_FOOT.length; i++ ) {
      if( unitString.equalsIgnoreCase(TEXT_MICROSECONDS_PER_FOOT[i]) )
        return MICROSECONDS_PER_FOOT;
    }
    if( unitString.equalsIgnoreCase(TEXT_SECONDS_PER_METER) )
      return SECONDS_PER_METER;
    return UNKNOWN_UNIT;
  }
  /**
   * Unit conversions
   */
  public static double meter2feet( double valueInMeter ) {
    return valueInMeter / FEET_TO_METER;
  }
  public static double feet2meter( double valueInFeet ) {
    return valueInFeet * FEET_TO_METER;
  }
  public static double cm2inch( double valueInCM ) {
    return valueInCM / INCH_TO_CM;
  }
  public static double inch2cm( double valueInInch ) {
    return valueInInch * INCH_TO_CM;
  }
  /**
   * Unit conversions
   */
  public static double knots2kmPerHour( double valueInKnots ) {
    return valueInKnots * ONE_KNOT_IN_KM_PER_HOUR;
  }
  public static double kmPerHour2knots( double valueInKmPerHour ) {
    return valueInKmPerHour / ONE_KNOT_IN_KM_PER_HOUR;
  }
  /**
   * Convert microseconds per foot to meters per second
   *
   * @param valueInUsPerFt double
   * @return double
   */
  public static double us_ft_2_m_s( double valueInUsPerFt ) {
    return( 1000000.0 * FEET_TO_METER / valueInUsPerFt );
  }
  /**
   * Convert array of values into seconds
   * Input array values must be one of the following units:
   * METERS_PER_SECOND, SECONDS_PER_METER, MICROSECONDS_PER_FOOT
   *
   * @param depthValues double[]
   *      Input depths
   * @param values double[]
   *      Input array
   * @param unitIN int
   *      Unit of input array
   * @return double[]
   *      Output array with values in seconds
   */
  public static double[] convertToSeconds( double[] depthValues, double[] values, int unitIN ) {
    double[] valuesSeconds = new double[values.length];
    valuesSeconds[0] = 0.0;
    switch( unitIN ) {
      case METERS_PER_SECOND :
        for( int i = 1; i < values.length; i++ ) {
          valuesSeconds[i] = ( depthValues[i] - depthValues[i-1] )/values[i-1] + valuesSeconds[i-1];
        }
        break;
      case SECONDS_PER_METER :
        for( int i = 1; i < values.length; i++ ) {
          valuesSeconds[i] = ( depthValues[i] - depthValues[i-1] )*values[i-1] + valuesSeconds[i-1];
        }
        break;
      case MICROSECONDS_PER_FOOT :
        for( int i = 1; i < values.length; i++ ) {
          valuesSeconds[i] = ( depthValues[i] - depthValues[i-1] )/us_ft_2_m_s(values[i-1]) + valuesSeconds[i-1];
        }
        break;
      default:
        return null; // Program bug
    }
    return valuesSeconds;
  }

  /**
   * Convert between feet and meter
   *
   * @param value double
   *       Input value
   * @param valueUnit int
   *       Unit of input value
   * @return double
   *       Output value
   */
  public static double convertFeetMeter( double value, int valueUnit ) {
    switch( valueUnit ) {
      case FEET:
        return feet2meter(value);
      case METER:
        return meter2feet(value);
    }
    // Exception
    return value;
  }
}


