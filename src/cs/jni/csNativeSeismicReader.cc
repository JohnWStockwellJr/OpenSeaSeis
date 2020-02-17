/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_jni_csNativeSeismicReader.h"
#include "csGeneralSeismicReader.h"
#include "csSeismicIOConfig.h"
#include "csException.h"
#include "csFlexHeader.h"
#include "csFlexNumber.h"
#include "geolib_defines.h"
#include <string>

/**
 * JNI Native Seismic Reader
 *
 * @author Bjorn Olofsson
 * @date 2007
 */

namespace cseis_jni_csNativeSeismicReader {
  // Global cache of method ids:
  static jmethodID id_setSeismicTraceFloatHeader;
  static jmethodID id_setSeismicTraceDoubleHeader;
  static jmethodID id_setSeismicTraceIntHeader;
  static jmethodID id_setSeismicTraceLongHeader;
  static jmethodID id_setSeismicTraceStringHeader;

  static jmethodID id_setHeaderFloatHeader;
  static jmethodID id_setHeaderDoubleHeader;
  static jmethodID id_setHeaderIntHeader;
  static jmethodID id_setHeaderLongHeader;
}

using namespace cseis_geolib;
using namespace cseis_io;
using namespace cseis_jni_csNativeSeismicReader;

/*
 * Class:     cseis_jni_csNativeSeismicReader
 * Method:    native_createInstance
 * Signature: (Ljava/lang/String;I)J
 */
JNIEXPORT jlong JNICALL Java_cseis_jni_csNativeSeismicReader_native_1createInstance
(JNIEnv *env, jobject obj, jstring filename_in, jint numTracesBuffer )
{
  char const* filename = (env)->GetStringUTFChars( filename_in, NULL );

  csGeneralSeismicReader* reader = NULL;

  try {
    reader = new csGeneralSeismicReader( filename, true, numTracesBuffer );
  }
  catch( csException& e ) {
    fprintf( stderr, "\nError when opening SeaSeis file: %s\n\n", e.getMessage() );
    fflush(stderr);
    jclass newExcCls = (env)->FindClass( "java/lang/Exception");
    if( newExcCls == NULL ) {
      return 0;
    }
    else {
      (env)->ThrowNew( newExcCls, e.getMessage() );
      return 0;
    }
  }

  try {
    bool success = reader->readFileHeader();
    if( !success ) {
      printf("Unknown error occurred when reading file header of SeaSeis input file.\n");
      fflush(stdout);
      return 0;
    }
  }
  catch( csException& e ) {
    fprintf( stderr, "Error when opening SeaSeis file: %s\n", e.getMessage() );
    fflush(stderr);
    jclass newExcCls = (env)->FindClass( "java/lang/Exception");
    if( newExcCls == NULL ) {
      return 0;
    }
    else {
      (env)->ThrowNew( newExcCls, e.getMessage() );
      return 0;
    }
  }

  // Get csSeismicTrace class and method IDs
  jclass class_csSeismicTrace  = env->FindClass("cseis/seis/csSeismicTrace");
  id_setSeismicTraceFloatHeader  = env->GetMethodID(class_csSeismicTrace,"setFloatHeader","(IF)V");
  id_setSeismicTraceDoubleHeader = env->GetMethodID(class_csSeismicTrace,"setDoubleHeader","(ID)V");
  id_setSeismicTraceIntHeader    = env->GetMethodID(class_csSeismicTrace,"setIntHeader","(II)V");
  id_setSeismicTraceLongHeader   = env->GetMethodID(class_csSeismicTrace,"setLongHeader","(IJ)V");
  id_setSeismicTraceStringHeader = env->GetMethodID(class_csSeismicTrace,"setStringHeader","(ILjava/lang/String;)V");

  jclass class_csHeader  = env->FindClass("cseis/seis/csHeader");
  id_setHeaderFloatHeader  = env->GetMethodID(class_csHeader,"setValue","(F)V");
  id_setHeaderDoubleHeader = env->GetMethodID(class_csHeader,"setValue","(D)V");
  id_setHeaderIntHeader    = env->GetMethodID(class_csHeader,"setValue","(I)V");
  id_setHeaderLongHeader   = env->GetMethodID(class_csHeader,"setValue","(J)V");

  return ( reinterpret_cast<jlong>( reader ) );
}

/*
 * Class:     cseis_jni_csNativeSeismicReader
 * Method:    native_closeFile
 * Signature: (J)V
 */
JNIEXPORT void JNICALL Java_cseis_jni_csNativeSeismicReader_native_1closeFile
(JNIEnv *env, jobject obj, jlong ptr_in)
{
  csGeneralSeismicReader* ptr = reinterpret_cast<csGeneralSeismicReader*>(ptr_in);
  ptr->closeFile();
}

/*
 * Class:     cseis_jni_csNativeSeismicReader
 * Method:    native_numTraces
 * Signature: (J)I
 */
JNIEXPORT jint JNICALL Java_cseis_jni_csNativeSeismicReader_native_1numTraces
(JNIEnv *env, jobject obj, jlong ptr_in )
{
  csGeneralSeismicReader* ptr = reinterpret_cast<csGeneralSeismicReader*>(ptr_in);
  csSeismicIOConfig const* config = ptr->getSeismicConfig();
  return( config->numTraces );
}

/*
 * Class:     cseis_jni_csNativeSeismicReader
 * Method:    native_numSamples
 * Signature: (J)I
 */
JNIEXPORT jint JNICALL Java_cseis_jni_csNativeSeismicReader_native_1numSamples
(JNIEnv *env, jobject obj, jlong ptr_in )
{
  csGeneralSeismicReader* ptr = reinterpret_cast<csGeneralSeismicReader*>(ptr_in);
  csSeismicIOConfig const* config = ptr->getSeismicConfig();
  return( config->numSamples );
}

/*
 * Class:     cseis_jni_csNativeSeismicReader
 * Method:    native_numHeaders
 * Signature: (J)I
 */
JNIEXPORT jint JNICALL Java_cseis_jni_csNativeSeismicReader_native_1numHeaders
(JNIEnv *env, jobject obj, jlong ptr_in )
{
  csGeneralSeismicReader* ptr = reinterpret_cast<csGeneralSeismicReader*>(ptr_in);
  return( ptr->numTraceHeaders() );
}

/*
 * Class:     cseis_jni_csNativeSeismicReader
 * Method:    native_sampleInt
 * Signature: (J)F
 */
JNIEXPORT jfloat JNICALL Java_cseis_jni_csNativeSeismicReader_native_1sampleInt
(JNIEnv *env, jobject obj, jlong ptr_in )
{
  csGeneralSeismicReader* ptr = reinterpret_cast<csGeneralSeismicReader*>(ptr_in);
  csSeismicIOConfig const* config = ptr->getSeismicConfig();
  return( config->sampleInt );
}

/*
 * Class:     cseis_jni_csNativeSeismicReader
 * Method:    native_hdrIntValue
 * Signature: (JI)I
 */
JNIEXPORT jint JNICALL Java_cseis_jni_csNativeSeismicReader_native_1hdrIntValue
(JNIEnv *env, jobject obj, jlong ptr_in, jint hdrIndex )
{
  csGeneralSeismicReader* ptr = reinterpret_cast<csGeneralSeismicReader*>(ptr_in);
  return( ptr->hdrIntValue( hdrIndex ) );
}

/*
 * Class:     cseis_jni_csNativeSeismicReader
 * Method:    native_hdrFloatValue
 * Signature: (JI)F
 */
JNIEXPORT jfloat JNICALL Java_cseis_jni_csNativeSeismicReader_native_1hdrFloatValue
(JNIEnv *env, jobject obj, jlong ptr_in, jint hdrIndex )
{
  csGeneralSeismicReader* ptr = reinterpret_cast<csGeneralSeismicReader*>(ptr_in);
  return( ptr->hdrFloatValue( hdrIndex ) );
}

/*
 * Class:     cseis_jni_csNativeSeismicReader
 * Method:    native_hdrDoubleValue
 * Signature: (JI)D
 */
JNIEXPORT jdouble JNICALL Java_cseis_jni_csNativeSeismicReader_native_1hdrDoubleValue
(JNIEnv *env, jobject obj, jlong ptr_in, jint hdrIndex )
{
  csGeneralSeismicReader* ptr = reinterpret_cast<csGeneralSeismicReader*>(ptr_in);
  return( ptr->hdrDoubleValue( hdrIndex ) );
}

/*
 * Class:     cseis_jni_csNativeSeismicReader
 * Method:    native_hdrLongValue
 * Signature: (JI)J
 */
JNIEXPORT jlong JNICALL Java_cseis_jni_csNativeSeismicReader_native_1hdrLongValue
(JNIEnv *env, jobject obj, jlong ptr_in, jint hdrIndex )
{
  csGeneralSeismicReader* ptr = reinterpret_cast<csGeneralSeismicReader*>(ptr_in);
  return( ptr->hdrInt64Value( hdrIndex ) );
}

/*
 * Class:     cseis_jni_csNativeSeismicReader
 * Method:    native_hdrStringValue
 * Signature: (JI)Ljava/lang/String;
 */
JNIEXPORT jstring JNICALL Java_cseis_jni_csNativeSeismicReader_native_1hdrStringValue
(JNIEnv *env, jobject obj, jlong ptr_in, jint hdrIndex )
{
  csGeneralSeismicReader* ptr = reinterpret_cast<csGeneralSeismicReader*>(ptr_in);
  return (env)->NewStringUTF( ptr->hdrStringValue( hdrIndex ).c_str() );
}

/*
 * Class:     cseis_jni_csNativeSeismicReader
 * Method:    native_getNextTrace
 * Signature: (J[FLcseis/seis/csSeismicTrace;)Z
 */
JNIEXPORT jboolean JNICALL Java_cseis_jni_csNativeSeismicReader_native_1getNextTrace
(JNIEnv *env, jobject obj, jlong ptr_in, jfloatArray samples_in, jobject trace_in )
{
  csGeneralSeismicReader* ptr = reinterpret_cast<csGeneralSeismicReader*>(ptr_in);

  // Get csSeismicTrace class
  //  jclass class_csSeismicTrace  = env->GetObjectClass(trace_in);
  // jmethodID id_setFloatHeader  = env->GetMethodID(class_csSeismicTrace,"setFloatHeader","(IF)V");
  // jmethodID id_setDoubleHeader = env->GetMethodID(class_csSeismicTrace,"setDoubleHeader","(ID)V");
  // jmethodID id_setIntHeader    = env->GetMethodID(class_csSeismicTrace,"setIntHeader","(II)V");
  // jmethodID id_setLongHeader   = env->GetMethodID(class_csSeismicTrace,"setLongHeader","(IJ)V");
  // jmethodID id_setStringHeader = env->GetMethodID(class_csSeismicTrace,"setStringHeader","(ILjava/lang/String;)V");

  //  fprintf(stderr,"getNextTrace Current start trace index = %d\n", ptr->getCurrentTraceIndex());
  //  fflush(stderr);

  float const* samplesPtr = NULL;

  try {
    samplesPtr = ptr->readTraceReturnPointer();
  }
  catch( csException& e ) {
    printf( "Error when reading new trace: %s\n", e.getMessage() );
    fflush(stderr);
    jclass newExcCls = (env)->FindClass( "java/lang/Exception");
    if( newExcCls == NULL ) {
      return 0;
    }
    else {
      (env)->ThrowNew( newExcCls, e.getMessage() );
      return 0;
    }
  }

  if( samplesPtr == NULL ) {
    return JNI_FALSE;
  }

  csSeismicIOConfig const* config = ptr->getSeismicConfig();
  jsize numSamples = config->numSamples;
  env->SetFloatArrayRegion( samples_in, 0, numSamples, (jfloat*)samplesPtr );

  int numHeaders = ptr->numTraceHeaders();

  for( int ihdr = 0; ihdr < numHeaders; ihdr++ ) {
    type_t type = ptr->headerType(ihdr);
    //    fprintf(stdout," Header int value: %d %d %s\n", ihdr, ptr->hdrIntValue(ihdr), ptr->headerName( ihdr ).c_str() );
    // fflush(stdout);

    switch( type ) {
      case TYPE_DOUBLE:
        env->CallVoidMethod( trace_in, id_setSeismicTraceDoubleHeader, ihdr, ptr->hdrDoubleValue(ihdr) );
        break;
      case TYPE_FLOAT:
        env->CallVoidMethod( trace_in, id_setSeismicTraceFloatHeader, ihdr, ptr->hdrFloatValue(ihdr) );
        break;
      case TYPE_INT:
        env->CallVoidMethod( trace_in, id_setSeismicTraceIntHeader, ihdr, ptr->hdrIntValue(ihdr) );
        break;
      case TYPE_INT64:
        env->CallVoidMethod( trace_in, id_setSeismicTraceLongHeader, ihdr, ptr->hdrInt64Value(ihdr) );
        break;
      case TYPE_STRING:
        //
        env->CallVoidMethod( trace_in, id_setSeismicTraceStringHeader, ihdr, (env)->NewStringUTF( ptr->hdrStringValue(ihdr).c_str() ) );
        break;
      default: //
        env->CallVoidMethod( trace_in, id_setSeismicTraceIntHeader, ihdr, ptr->hdrIntValue(ihdr) );
        break;
    }
  }

  return JNI_TRUE;
}

/*
 * Class:     cseis_jni_csNativeSeismicReader
 * Method:    native_moveToTrace
 * Signature: (JII)Z
 */
JNIEXPORT jboolean JNICALL Java_cseis_jni_csNativeSeismicReader_native_1moveToTrace
(JNIEnv *env, jobject obj, jlong ptr_in, jint firstTraceIndex, jint numTracesToRead )
{
  csGeneralSeismicReader* ptr = reinterpret_cast<csGeneralSeismicReader*>(ptr_in);

  try {
    //    fprintf(stderr,"JNI Move to trace %d %d\n",firstTraceIndex, numTracesToRead );
    if( ptr->moveToTrace( firstTraceIndex, numTracesToRead ) ) {
      return JNI_TRUE;
    }
  }
  catch( csException& e ) {
    fprintf( stderr, "Error when moving to new position in SeaSeis file: %s\n", e.getMessage() );
    fflush(stderr);
    jclass newExcCls = (env)->FindClass( "java/lang/Exception");
    if( newExcCls == NULL ) {
      return 0;
    }
    else {
      (env)->ThrowNew( newExcCls, e.getMessage() );
      return 0;
    }
  }
  return JNI_FALSE;

}

/*
 * Class:     cseis_jni_csNativeSeismicReader
 * Method:    native_headerName
 * Signature: (JI)Ljava/lang/String;
 */
JNIEXPORT jstring JNICALL Java_cseis_jni_csNativeSeismicReader_native_1headerName
(JNIEnv *env, jobject obj, jlong ptr_in, jint hdrIndex )
{
  csGeneralSeismicReader* ptr = reinterpret_cast<csGeneralSeismicReader*>(ptr_in);
  return (env)->NewStringUTF(ptr->headerName( hdrIndex ).c_str());
}

/*
 * Class:     cseis_jni_csNativeSeismicReader
 * Method:    native_headerDesc
 * Signature: (JI)Ljava/lang/String;
 */
JNIEXPORT jstring JNICALL Java_cseis_jni_csNativeSeismicReader_native_1headerDesc
(JNIEnv *env, jobject obj, jlong ptr_in, jint hdrIndex )
{
  csGeneralSeismicReader* ptr = reinterpret_cast<csGeneralSeismicReader*>(ptr_in);
  return (env)->NewStringUTF(ptr->headerDesc( hdrIndex ).c_str() );
}

/*
 * Class:     cseis_jni_csNativeSeismicReader
 * Method:    native_headerType
 * Signature: (JI)I
 */
JNIEXPORT jint JNICALL Java_cseis_jni_csNativeSeismicReader_native_1headerType
(JNIEnv *env, jobject obj, jlong ptr_in, jint hdrIndex )
{
  csGeneralSeismicReader* ptr = reinterpret_cast<csGeneralSeismicReader*>(ptr_in);
  int headerType = (int)ptr->headerType( hdrIndex );
  return( headerType );
}

/*
 * Class:     cseis_jni_csNativeSeismicReader
 * Method:    native_isFrequencyDomain
 * Signature: (J)Z
 */
JNIEXPORT jboolean JNICALL Java_cseis_jni_csNativeSeismicReader_native_1isFrequencyDomain
(JNIEnv *env, jobject obj, jlong ptr_in ) {
  csGeneralSeismicReader* ptr = reinterpret_cast<csGeneralSeismicReader*>(ptr_in);
  int domain = ptr->getSeismicConfig()->domain;
  bool isFreq = (domain == DOMAIN_FX || domain == DOMAIN_FK);
  if( isFreq ) {
    return JNI_TRUE;
  }
  else {
    return JNI_FALSE;
  }
}

/*
 * Class:     cseis_jni_csNativeSeismicReader
 * Method:    native_setHeaderToPeek
 * Signature: (JLjava/lang/String;)V
 */
JNIEXPORT void JNICALL Java_cseis_jni_csNativeSeismicReader_native_1setHeaderToPeek
(JNIEnv *env, jobject obj, jlong ptr_in, jstring headerName_in ) {
  csGeneralSeismicReader* ptr = reinterpret_cast<csGeneralSeismicReader*>(ptr_in);
  char const* headerName = (env)->GetStringUTFChars( headerName_in, NULL );
  std::string headerNameString( headerName );
  ptr->setHeaderToPeek( headerNameString );
}

/*
 * Class:     cseis_jni_csNativeSeismicReader
 * Method:    native_peekHeaderValue
 * Signature: (JILjava/lang/Number;)V
 */
JNIEXPORT void JNICALL Java_cseis_jni_csNativeSeismicReader_native_1peekHeaderValue
(JNIEnv *env, jobject obj, jlong ptr_in, jint traceIndex_in, jobject value_in ) {
  csGeneralSeismicReader* ptr = reinterpret_cast<csGeneralSeismicReader*>(ptr_in);

  cseis_geolib::csFlexHeader flexValue;
  bool success = false;
  std::string message = "";
  try {
    success = ptr->peekHeaderValue( &flexValue, traceIndex_in );
  }
  catch( csException& e ) {
    message = std::string(e.getMessage());
    //    fprintf(stderr,"PEEK ERROR: %s\n", message.c_str());
    //    fflush(stderr);
  }
  if( !success ) {
    //    fprintf(stderr,"Peek operation failed... %d / %d\n", (int)traceIndex_in, ptr->numTraces());
    //   fflush(stderr);
    jclass newExcCls = (env)->FindClass( "java/lang/Exception");
    if( newExcCls == NULL ) {
      return;
    }
    else {
      if( message.length() > 0 ) {
        (env)->ThrowNew( newExcCls, message.c_str() );
      }
      else {
        (env)->ThrowNew( newExcCls, "Peek operation failed. Last trace..?" );
      }
      return;
    }
  }

  type_t type = flexValue.type();
  switch( type ) {
  case TYPE_DOUBLE:
    env->CallVoidMethod( value_in, id_setHeaderDoubleHeader, flexValue.doubleValue() );
    break;
  case TYPE_FLOAT:
    env->CallVoidMethod( value_in, id_setHeaderFloatHeader, flexValue.floatValue() );
    break;
  case TYPE_INT:
    env->CallVoidMethod( value_in, id_setHeaderIntHeader, flexValue.intValue() );
    break;
  case TYPE_INT64:
    env->CallVoidMethod( value_in, id_setHeaderLongHeader, flexValue.int64Value() );
    break;
  default: //
    env->CallVoidMethod( value_in, id_setHeaderIntHeader, flexValue.intValue() );
    break;
  }

}

/*
 * Class:     cseis_jni_csNativeSeismicReader
 * Method:    native_freeInstance
 * Signature: (J)V
 */
JNIEXPORT void JNICALL Java_cseis_jni_csNativeSeismicReader_native_1freeInstance
(JNIEnv *env, jobject obj, jlong ptr_in) {
  csGeneralSeismicReader* ptr = reinterpret_cast<csGeneralSeismicReader*>(ptr_in);
  delete ptr;
}

/*
 * Class:     cseis_jni_csNativeSeismicReader
 * Method:    native_verticalDomain
 * Signature: (J)I
 */
JNIEXPORT jint JNICALL Java_cseis_jni_csNativeSeismicReader_native_1verticalDomain
(JNIEnv *env, jobject obj, jlong ptr_in) {
  csGeneralSeismicReader* ptr = reinterpret_cast<csGeneralSeismicReader*>(ptr_in);

  int domain = ptr->getSeismicConfig()->domain;
  if( domain == cseis_geolib::DOMAIN_XT || domain == cseis_geolib::DOMAIN_KT ) {
    return cseis_geolib::VERTICAL_DOMAIN_TIME;
  }
  else if( domain == cseis_geolib::DOMAIN_XD ) {
    return cseis_geolib::VERTICAL_DOMAIN_DEPTH;
  }
  else if( domain == cseis_geolib::DOMAIN_FX || domain == cseis_geolib::DOMAIN_FK ) {
    return cseis_geolib::VERTICAL_DOMAIN_FREQ;
  }
  return(  cseis_geolib::VERTICAL_DOMAIN_TIME );
}

/*
 * Class:     cseis_jni_csNativeSeismicReader
 * Method:    native_setSelection
 * Signature: (JLjava/lang/String;Ljava/lang/String;IILcseis/jni/csITraceHeaderScanNotifier;)Z
 */
JNIEXPORT jboolean JNICALL Java_cseis_jni_csNativeSeismicReader_native_1setSelection
(JNIEnv *env, jobject obj, jlong ptr_in, jstring headerValueSelectionText_in, jstring headerName_in, jint sortOrder, jint sortMethod, jobject notifier_in )
{
  csGeneralSeismicReader* reader = reinterpret_cast<csGeneralSeismicReader*>(ptr_in);
  char const* hdrValueSelectionText = (env)->GetStringUTFChars( headerValueSelectionText_in, NULL );
  char const* headerName = (env)->GetStringUTFChars( headerName_in, NULL );

  jclass class_csITraceHeaderScanNotifier  = env->GetObjectClass(notifier_in);
  jmethodID id_notify   = env->GetMethodID(class_csITraceHeaderScanNotifier,"traceHeaderScanNotify","(I)V");
  jmethodID id_continue = env->GetMethodID(class_csITraceHeaderScanNotifier,"traceHeaderScanContinue","()Z");
  int numTracesToRead = 100;

  try {
    reader->setSelectionStep1( hdrValueSelectionText, headerName, sortOrder, sortMethod );
    int traceIndex = 0;
    while( !reader->setSelectionStep2( numTracesToRead ) ) {
      traceIndex += numTracesToRead;
      env->CallVoidMethod( notifier_in, id_notify, traceIndex );
      if( env->CallBooleanMethod( notifier_in, id_continue ) == JNI_FALSE ) break;
    }
    if( !reader->setSelectionStep3( ) ) return JNI_FALSE;
  }
  catch( csException& e ) {
    fprintf( stderr, "\nError when setting trace selection/sorting: %s\n\n", e.getMessage() );
    fflush(stderr);
    jclass newExcCls = (env)->FindClass( "java/lang/Exception");
    if( newExcCls == NULL ) {
      return JNI_FALSE;
    }
    else {
      (env)->ThrowNew( newExcCls, e.getMessage() );
      return JNI_FALSE;
    }
  }
  return JNI_TRUE;
}

/*
 * Class:     cseis_jni_csNativeSeismicReader
 * Method:    native_getNumSelectedTraces
 * Signature: (J)I
 */
JNIEXPORT jint JNICALL Java_cseis_jni_csNativeSeismicReader_native_1getNumSelectedTraces
(JNIEnv *env, jobject obj, jlong ptr_in ) {
  csGeneralSeismicReader* reader = reinterpret_cast<csGeneralSeismicReader*>(ptr_in);
  return reader->getNumSelectedTraces();
}


/*
 * Class:     cseis_jni_csNativeSeismicReader
 * Method:    native_getSelectedValue
 * Signature: (JILcseis/seis/csHeader;)V
 */
JNIEXPORT void JNICALL Java_cseis_jni_csNativeSeismicReader_native_1getSelectedValue
(JNIEnv *env, jobject obj, jlong ptr_in, jint selectedTraceIndex, jobject hdrValue_in ) {
  csGeneralSeismicReader* reader = reinterpret_cast<csGeneralSeismicReader*>(ptr_in);

  // Get csHeader class
  //  jclass class_csHeader  = env->GetObjectClass(hdrValue_in);
  //  jmethodID id_setFloatHeader  = env->GetMethodID(class_csHeader,"setValue","(F)V");
  //  jmethodID id_setDoubleHeader = env->GetMethodID(class_csHeader,"setValue","(D)V");
  // jmethodID id_setIntHeader    = env->GetMethodID(class_csHeader,"setValue","(I)V");
  // jmethodID id_setLongHeader   = env->GetMethodID(class_csHeader,"setValue","(J)V");

  cseis_geolib::csFlexNumber const* flexValue = reader->getSelectedValue( selectedTraceIndex );
  
  type_t type = flexValue->type();
  switch( type ) {
  case TYPE_DOUBLE:
    env->CallVoidMethod( hdrValue_in, id_setHeaderDoubleHeader, flexValue->doubleValue() );
    break;
  case TYPE_FLOAT:
    env->CallVoidMethod( hdrValue_in, id_setHeaderFloatHeader, flexValue->floatValue() );
    break;
  case TYPE_INT:
    env->CallVoidMethod( hdrValue_in, id_setHeaderIntHeader, flexValue->intValue() );
    break;
  case TYPE_INT64:
    env->CallVoidMethod( hdrValue_in, id_setHeaderLongHeader, flexValue->int64Value() );
    break;
  default: //
    env->CallVoidMethod( hdrValue_in, id_setHeaderIntHeader, flexValue->intValue() );
    break;
  }
}

/*
 * Class:     cseis_jni_csNativeSeismicReader
 * Method:    native_getSelectedValues
 * Signature: (JLcseis/jni/csSelectedHeaderBundle;)V
 */
JNIEXPORT void JNICALL Java_cseis_jni_csNativeSeismicReader_native_1getSelectedValues
(JNIEnv *env, jobject obj, jlong ptr_in, jobject hdrBundle_in ) {
  csGeneralSeismicReader* reader = reinterpret_cast<csGeneralSeismicReader*>(ptr_in);

  jclass class_csSelectedHeaderBundle = env->GetObjectClass( hdrBundle_in );

  jmethodID id_setFloatHeader  = env->GetMethodID(class_csSelectedHeaderBundle,"setFloatValue","(IF)V");
  jmethodID id_setDoubleHeader = env->GetMethodID(class_csSelectedHeaderBundle,"setDoubleValue","(ID)V");
  jmethodID id_setIntHeader    = env->GetMethodID(class_csSelectedHeaderBundle,"setIntValue","(II)V");
  jmethodID id_setLongHeader   = env->GetMethodID(class_csSelectedHeaderBundle,"setLongValue","(IJ)V");
  jmethodID id_setTraceIndex   = env->GetMethodID(class_csSelectedHeaderBundle,"setSortedTraceIndex","(II)V");

  if( id_setFloatHeader <= 0  || id_setDoubleHeader <= 0  || id_setIntHeader <= 0  || id_setLongHeader <= 0 ) {
    fprintf(stderr,"Error in JNI layer, csNativeSeismicReader: getSelectedValues. Cannot find Java methods");
    return;
  }

  for( int itrc = 0; itrc < reader->getNumSelectedTraces(); itrc++ ) {
    cseis_geolib::csFlexNumber const* flexValue = reader->getSelectedValue( itrc );
    env->CallVoidMethod( hdrBundle_in, id_setTraceIndex, itrc, reader->getSelectedIndex(itrc) );
    type_t type = flexValue->type();
    switch( type ) {
    case TYPE_DOUBLE:
      env->CallVoidMethod( hdrBundle_in, id_setDoubleHeader, itrc, flexValue->doubleValue() );
      break;
    case TYPE_FLOAT:
      env->CallVoidMethod( hdrBundle_in, id_setFloatHeader, itrc, flexValue->floatValue() );
      break;
    case TYPE_INT:
      env->CallVoidMethod( hdrBundle_in, id_setIntHeader, itrc, flexValue->intValue() );
      break;
    case TYPE_INT64:
      env->CallVoidMethod( hdrBundle_in, id_setLongHeader, itrc, flexValue->int64Value() );
      break;
    default: //
      env->CallVoidMethod( hdrBundle_in, id_setIntHeader, itrc, flexValue->intValue() );
      break;
    }
  }

}
