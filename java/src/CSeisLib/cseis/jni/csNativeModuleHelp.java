/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.jni;

/**
 * JNI interface to Module help
 * @author Bjorn Olofsson
 */

public class csNativeModuleHelp {
  private native String native_moduleHtmlHelp( String moduleName );
  private native String native_moduleHtmlExample( String moduleName );
  private native String native_moduleExample( String moduleName );
  private native String native_moduleHtmlListing();
  private native String native_standardHeaderHtmlListing();

  public csNativeModuleHelp() {
  }
  public String moduleHtmlHelp( String moduleName ) {
    return native_moduleHtmlHelp( moduleName );
  }
  public String moduleHtmlExample( String moduleName ) {
    return native_moduleHtmlExample( moduleName );
  }
  public String moduleExample( String moduleName ) {
    return native_moduleExample( moduleName );
  }
  public String moduleHtmlListing() {
    return native_moduleHtmlListing();
  }
  public String standardHeaderHtmlListing() {
//    return "NOTHING";
    return native_standardHeaderHtmlListing();
  }
}


