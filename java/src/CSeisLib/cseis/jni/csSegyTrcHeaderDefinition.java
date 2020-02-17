/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.jni;

import java.util.*;

/**
 * SEGY trace header definition
 * @author Bjorn Olofsson
 */
public class csSegyTrcHeaderDefinition {
  private ArrayList<Integer> myByteLoc;
  private ArrayList<Integer> myByteSize;
  private ArrayList<String> myName;
  private ArrayList<String> myDesc;
  private String myTitle;

  /**
   * Constructor
   *
   * @param title Identifier/title of trace header definition
   */
  public csSegyTrcHeaderDefinition( String title ) {
    myTitle = title;
    myByteLoc  = new ArrayList<Integer>();
    myByteSize = new ArrayList<Integer>();
    myName     = new ArrayList<String>();
    myDesc     = new ArrayList<String>();
  }
  /**
   * Set identifying name/title of SEGY trace header definition
   * @param title Identifier/title of SEGY trace header definition
   */
  public void setTitle( String title ) {
    myTitle = title;
  }
  /**
   *
   * @return Indentifying name/title of SEGY trace header definition
   */
  public String title() {
    return myTitle;
  }
  /**
   * Add trace header
   *
   * @param name     Name of trace header
   * @param desc     Description text
   * @param byteLoc  Byte location
   * @param byteSize Byte size
   */
  public void add( String name, String desc, int byteLoc, int byteSize ) {
    int index = getInsertIndex( byteLoc );
    myByteLoc.add( index, byteLoc );
    myByteSize.add( index, byteSize );
    myName.add( index, name );
    myDesc.add( index, desc );
  }
  private int getInsertIndex( int byteLoc ) {
    int index1 = 0;
    int index2 = myByteLoc.size()-1;
    if( index2 < index1 ) return 0;
    if( byteLoc <= myByteLoc.get(index1) ) return index1;
    if( byteLoc >= myByteLoc.get(index2) ) return myByteLoc.size();
    int index = (index1 + index2) / 2;
    while( index2 > index1+1 ) {
      if( byteLoc < myByteLoc.get(index) ) {
        index2 = index;
        index = (index1 + index2) / 2;
      }
      else if( byteLoc > myByteLoc.get(index) ) {
        index1 = index;
        index = (index1 + index2 + 1) / 2;
      }
      else {  // byteLoc == myByteLoc.get(index)
        return index;
      }
    }
    if( byteLoc > myByteLoc.get(index) ) index += 1;
    return index;
  }
  public void clear() {
    myByteLoc.clear();
    myByteSize.clear();
    myName.clear();
    myDesc.clear();
  }
  public int numHeaders() {
    return myName.size();
  }
  public int byteLoc( int index ) {
    return myByteLoc.get(index);
  }
  public int byteSize( int index ) {
    return myByteSize.get(index);
  }
  public String name( int index ) {
    return myName.get(index);
  }
  public String desc( int index ) {
    return myDesc.get(index);
  }
}


