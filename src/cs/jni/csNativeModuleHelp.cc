/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_jni_csNativeModuleHelp.h"
#include "csHelp.h"
#include "csParamDef.h"
#include "csException.h"
#include "csMethodRetriever.h"
#include <string>

/*
 * Class:     cseis_jni_csNativeModuleHelp
 * Method:    native_moduleHtmlHelp
 * Signature: (Ljava/lang/String;)Ljava/lang/String;
 */
JNIEXPORT jstring JNICALL Java_cseis_jni_csNativeModuleHelp_native_1moduleHtmlHelp
(JNIEnv *env, jobject obj, jstring moduleName_in ) {
  char const* moduleName = (env)->GetStringUTFChars( moduleName_in, NULL );

  cseis_system::csParamDef pdef;
  try {
    cseis_system::csMethodRetriever::getParamMethod( moduleName )( &pdef );
  }
  catch( cseis_geolib::csException& e ) {
    return (env)->NewStringUTF(e.getMessage());
  }
  catch( ... ) {
    return (env)->NewStringUTF("NOT FOUND");
  }

  cseis_system::csHelp helpObj;
  std::string text = "";
  helpObj.moduleHtmlHelp( pdef, text );

  return (env)->NewStringUTF(text.c_str());
}

/*
 * Class:     cseis_jni_csNativeModuleHelp
 * Method:    native_moduleHtmlExample
 * Signature: (Ljava/lang/String;)Ljava/lang/String;
 */
JNIEXPORT jstring JNICALL Java_cseis_jni_csNativeModuleHelp_native_1moduleHtmlExample
(JNIEnv *env, jobject obj, jstring moduleName_in ) {
  char const* moduleName = (env)->GetStringUTFChars( moduleName_in, NULL );

  cseis_system::csParamDef pdef;
  try {
    cseis_system::csMethodRetriever::getParamMethod( moduleName )( &pdef );
  }
  catch( cseis_geolib::csException& e ) {
    return (env)->NewStringUTF(e.getMessage());
  }
  catch( ... ) {
    return (env)->NewStringUTF("NOT FOUND");
  }

  cseis_system::csHelp helpObj;
  std::string text = "";
  helpObj.moduleHtmlExample( true, pdef, text );

  return (env)->NewStringUTF(text.c_str());
}
/*
 * Class:     cseis_jni_csNativeModuleHelp
 * Method:    native_moduleExample
 * Signature: (Ljava/lang/String;)Ljava/lang/String;
 */
JNIEXPORT jstring JNICALL Java_cseis_jni_csNativeModuleHelp_native_1moduleExample
(JNIEnv *env, jobject obj, jstring moduleName_in ) {
  char const* moduleName = (env)->GetStringUTFChars( moduleName_in, NULL );

  cseis_system::csParamDef pdef;
  try {
    cseis_system::csMethodRetriever::getParamMethod( moduleName )( &pdef );
  }
  catch( cseis_geolib::csException& e ) {
    return (env)->NewStringUTF(e.getMessage());
  }
  catch( ... ) {
    return (env)->NewStringUTF("NOT FOUND");
  }

  cseis_system::csHelp helpObj;
  std::string text = "";
  helpObj.moduleExample( false, pdef, text );

  return (env)->NewStringUTF(text.c_str());
}

/*
 * Class:     cseis_jni_csNativeModuleHelp
 * Method:    native_moduleHtmlListing
 * Signature: ()Ljava/lang/String;
 */
JNIEXPORT jstring JNICALL Java_cseis_jni_csNativeModuleHelp_native_1moduleHtmlListing
(JNIEnv *env, jobject obj) {

  try {
    cseis_system::csHelp helpObj;
    std::string text = "";
    helpObj.moduleHtmlListing( text);
    return (env)->NewStringUTF(text.c_str());
  }
  catch( cseis_geolib::csException& e ) {
    return (env)->NewStringUTF(e.getMessage());
  }
  catch( ... ) {
    return (env)->NewStringUTF("NOT FOUND");
  }
}

/*
 * Class:     cseis_jni_csNativeModuleHelp
 * Method:    native_standardHeaderHtmlListing
 * Signature: ()Ljava/lang/String;
 */
JNIEXPORT jstring JNICALL Java_cseis_jni_csNativeModuleHelp_native_1standardHeaderHtmlListing
(JNIEnv *env, jobject obj) {

  try {
    cseis_system::csHelp helpObj;
    std::string text = "";
    helpObj.standardHeaderHtmlListing( text);
    return (env)->NewStringUTF(text.c_str());
  }
  catch( cseis_geolib::csException& e ) {
    return (env)->NewStringUTF(e.getMessage());
  }
  catch( ... ) {
    return (env)->NewStringUTF("NOT FOUND");
  }
}


