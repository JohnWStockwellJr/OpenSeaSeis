/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.tools;

/**
 * Listener interface which is used to communicate snap shot events.
 * 
 * @author 2011 Bjorn Olofsson
 */
public interface csISnapShotListener {
  /**
   * Select a snap shot.
   * @param itemIndex  Index of snapshot item to be selected.
   * @param forceSelection  true if selection shall be 'forced' in any case.
   */
  public void selectSnapShot( int itemIndex, boolean forceSelection );
  /**
   * Delete a snap shot.
   * @param itemIndex  Index of snap shot item.
   */
  public void deleteSnapShot( int itemIndex );
  /**
   * Move snap shot.
   * @param fromItemIndex  Index of current location of snap shot item.
   * @param toItemIndex    Index of new location of snap shot item.
   */
  public void moveSnapShot( int fromItemIndex, int toItemIndex );
}


