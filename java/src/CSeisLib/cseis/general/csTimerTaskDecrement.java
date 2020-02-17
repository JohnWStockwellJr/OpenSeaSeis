/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.general;

import java.util.TimerTask;

/**
 * Timed task to decrement something..
 * @author 2013 Felipe Punto
 */
public class csTimerTaskDecrement extends TimerTask {
  private csITimerTaskDecrementListener myListener;
  public csTimerTaskDecrement( csITimerTaskDecrementListener listener ) {
    myListener = listener;
  }
  @Override
  public void run() {
    myListener.decrement();
  }
}

