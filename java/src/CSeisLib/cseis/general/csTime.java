/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.general;

/**
 * Simple time representation
 * @author Bjorn Olofsson
 */
public class csTime {
  public int days;
  public int hours;
  public int minutes;
  public int seconds;

  public csTime() {
    days    = 0;
    hours   = 0;
    minutes = 0;
    seconds = 0;
  }
  public String toString() {
    return "";
  }
  public static csTime convertHours( double hoursOld ) {
    csTime time = new csTime();
    
    time.days    = (int)( hoursOld / 24.0 );
    time.hours   = (int)( hoursOld - (24.0*time.days) );
    time.minutes = (int)( (hoursOld - (24.0*time.days + time.hours) ) * 60.0 );
    time.seconds = (int)( ( 60.0*(hoursOld - (24.0*time.days + time.hours)) - (double)time.minutes ) * 60.0 );
    
    return time;
  }
}



