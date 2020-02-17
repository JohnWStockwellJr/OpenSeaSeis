/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.general;

import java.util.TimerTask;

/**
 * Timed task to increment something..
 * @author 2013 Felipe Punto
 */
public class csTimerTaskIncrement extends TimerTask {
  private csITimerTaskIncrementListener myListener;
  public csTimerTaskIncrement( csITimerTaskIncrementListener listener ) {
    myListener = listener;
  }
  @Override
  public void run() {
    myListener.increment();
  }
}

