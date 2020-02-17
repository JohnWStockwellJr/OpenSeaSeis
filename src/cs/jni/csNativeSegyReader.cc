/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_jni_csNativeSegyReader.h"
#include "csStandardHeaders.h"
#include "csSegyTraceHeader.h"
#include "csSegyHdrMap.h"
#include "csSegyReader.h"
#include "csSegyHeaderInfo.h"
#include "csSegyBinHeader.h"
#include "csException.h"
#include "csVector.h"
#include "csFlexHeader.h"
#include "csFlexNumber.h"
#include "geolib_defines.h"

/**
 * JNI Native Segy Reader
 *
 * @author Bjorn Olofsson
 * @date 2007
 */

namespace cseis_jni_csNativeSegyReader {
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
using namespace cseis_jni_csNativeSegyReader;


/*
 * Class:     cseis_jni_csNativeSegyReader
 * Method:    native_createInstance
 * Signature: (Ljava/lang/String;IIZZZ)J
 */
JNIEXPORT jlong JNICALL Java_cseis_jni_csNativeSegyReader_native_1createInstance
(JNIEnv *env, jobject obj,
  jstring filename_in,
  jint nTracesBuffer,
  jint segyHeaderMap,
  jboolean reverseByteOrderData_in,
  jboolean reverseByteOrderHdr_in,
  jboolean autoscale_hdrs_in )
{
  char const* filename = (env)->GetStringUTFChars( filename_in, NULL );

// -----------------------------------------


  csSegyReader::SegyReaderConfig config;
  config.enableRandomAccess = true;
  
  if( reverseByteOrderData_in == JNI_TRUE ) {
    config.reverseByteOrderData = true;
  }
  else {
    config.reverseByteOrderData = false;
  }
  if( reverseByteOrderHdr_in == JNI_TRUE ) {
    config.reverseByteOrderHdr = true;
  }
  else {
    config.reverseByteOrderHdr = false;
  }
  if( autoscale_hdrs_in == JNI_TRUE ) {
    config.autoscaleHdrs = true;
  }
  else {
    config.autoscaleHdrs = false;
  }
  config.segyHeaderMapping = segyHeaderMap;
  config.numTracesBuffer   = nTracesBuffer;
  // Force SU format when SU header map is selected:
  config.isSUFormat        = (segyHeaderMap == csSegyHdrMap::SEGY_SU || segyHeaderMap == csSegyHdrMap::SEGY_SU_BOTH);

  csSegyReader* reader = NULL;
  try {
    reader = new csSegyReader( filename, config );
  }
  catch( csException& e ) {
    fprintf( stderr,"Error when opening SEGY file: %s\n", e.getMessage() );
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
    reader->initialize();
  }
  catch( csException& e ) {
    reader = NULL;
    fprintf(stderr,"Error when initializing SEGY reader object.\nSystem message: %s", e.getMessage() );
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
 * Class:     cseis_jni_csNativeSegyReader
 * Method:    native_closeFile
 * Signature: (J)V
 */
JNIEXPORT void JNICALL Java_cseis_jni_csNativeSegyReader_native_1closeFile
(JNIEnv *env, jobject obj, jlong ptr_in )
{
  csSegyReader* ptr = reinterpret_cast<csSegyReader*>(ptr_in);
  ptr->closeFile();
}

/*
 * Class:     cseis_jni_csNativeSegyReader
 * Method:    native_sampleByteSize
 * Signature: (J)I
 */
JNIEXPORT jint JNICALL Java_cseis_jni_csNativeSegyReader_native_1sampleByteSize
(JNIEnv *env, jobject obj, jlong ptr_in )
{
  csSegyReader* ptr = reinterpret_cast<csSegyReader*>(ptr_in);
  return ptr->sampleByteSize();
}

/*
 * Class:     cseis_jni_csNativeSegyReader
 * Method:    native_numTraces
 * Signature: (J)I
 */
JNIEXPORT jint JNICALL Java_cseis_jni_csNativeSegyReader_native_1numTraces
(JNIEnv *env, jobject obj, jlong ptr_in )
{
  csSegyReader* ptr = reinterpret_cast<csSegyReader*>(ptr_in);
  return( ptr->numTraces() );
}

/*
 * Class:     cseis_jni_csNativeSegyReader
 * Method:    native_numSamples
 * Signature: (J)I
 */
JNIEXPORT jint JNICALL Java_cseis_jni_csNativeSegyReader_native_1numSamples
(JNIEnv *env, jobject obj, jlong ptr_in )
{
  csSegyReader* ptr = reinterpret_cast<csSegyReader*>(ptr_in);
  return( ptr->numSamples() );
}

/*
 * Class:     cseis_jni_csNativeSegyReader
 * Method:    native_numHeaders
 * Signature: (J)I
 */
JNIEXPORT jint JNICALL Java_cseis_jni_csNativeSegyReader_native_1numHeaders
(JNIEnv *env, jobject obj, jlong ptr_in )
{
  csSegyReader* ptr = reinterpret_cast<csSegyReader*>(ptr_in);
  return ptr->numTraceHeaders();
}

/*
 * Class:     cseis_jni_csNativeSegyReader
 * Method:    native_sampleInt
 * Signature: (J)F
 */
JNIEXPORT jfloat JNICALL Java_cseis_jni_csNativeSegyReader_native_1sampleInt
(JNIEnv *env, jobject obj, jlong ptr_in )
{
  csSegyReader* ptr = reinterpret_cast<csSegyReader*>(ptr_in);
  return( ptr->sampleIntMS() );
}

/*
 * Class:     cseis_jni_csNativeSegyReader
 * Method:    native_binHeader
 * Signature: (JLcseis/seis/csSegyBinHeader;)V
 */
JNIEXPORT void JNICALL Java_cseis_jni_csNativeSegyReader_native_1binHeader
(JNIEnv *env, jobject obj, jlong ptr_in, jobject binHdr_in )
{
  csSegyReader* ptr = reinterpret_cast<csSegyReader*>(ptr_in);
  cseis_geolib::csSegyBinHeader const* binHdr = ptr->binHdr();

  // Get csSeismicTrace class
  jclass class_binHdr = env->GetObjectClass(binHdr_in);
  env->SetIntField( binHdr_in, env->GetFieldID(class_binHdr,"jobID","I"), binHdr->jobID );
  env->SetIntField( binHdr_in, env->GetFieldID(class_binHdr,"lineNum","I"), binHdr->lineNum );
  env->SetIntField( binHdr_in, env->GetFieldID(class_binHdr,"reelNum","I"), binHdr->reelNum );

  env->SetShortField( binHdr_in, env->GetFieldID(class_binHdr,"numTraces","S"), (jshort)binHdr->numTraces );
  env->SetShortField( binHdr_in, env->GetFieldID(class_binHdr,"numAuxTraces","S"), (jshort)binHdr->numAuxTraces );
  env->SetShortField( binHdr_in, env->GetFieldID(class_binHdr,"sampleIntUS","S"), (jshort)binHdr->sampleIntUS );
  env->SetShortField( binHdr_in, env->GetFieldID(class_binHdr,"sampleIntOrigUS","S"), (jshort)binHdr->sampleIntOrigUS );
  env->SetShortField( binHdr_in, env->GetFieldID(class_binHdr,"numSamples","S"), (jshort)binHdr->numSamples );
  env->SetShortField( binHdr_in, env->GetFieldID(class_binHdr,"numSamplesOrig","S"), (jshort)binHdr->numSamplesOrig );

  env->SetShortField( binHdr_in, env->GetFieldID(class_binHdr,"dataSampleFormat","S"), binHdr->dataSampleFormat );
  env->SetShortField( binHdr_in, env->GetFieldID(class_binHdr,"fold","S"), binHdr->fold );
  env->SetShortField( binHdr_in, env->GetFieldID(class_binHdr,"sortCode","S"), binHdr->sortCode );
  env->SetShortField( binHdr_in, env->GetFieldID(class_binHdr,"vertSumCode","S"), binHdr->vertSumCode );
  env->SetShortField( binHdr_in, env->GetFieldID(class_binHdr,"sweepFreqStart","S"), binHdr->sweepFreqStart );
  env->SetShortField( binHdr_in, env->GetFieldID(class_binHdr,"sweepFreqEnd","S"), binHdr->sweepFreqEnd );
  env->SetShortField( binHdr_in, env->GetFieldID(class_binHdr,"taperType","S"), binHdr->taperType );
  env->SetShortField( binHdr_in, env->GetFieldID(class_binHdr,"correlatedTraces","S"), binHdr->correlatedTraces );
  env->SetShortField( binHdr_in, env->GetFieldID(class_binHdr,"gainRecovered","S"), binHdr->gainRecovered );
  env->SetShortField( binHdr_in, env->GetFieldID(class_binHdr,"ampRecoveryMethod","S"), binHdr->ampRecoveryMethod );
  env->SetShortField( binHdr_in, env->GetFieldID(class_binHdr,"unitSystem","S"), binHdr->unitSystem );
  env->SetShortField( binHdr_in, env->GetFieldID(class_binHdr,"polarity","S"), binHdr->polarity );
  env->SetShortField( binHdr_in, env->GetFieldID(class_binHdr,"vibPolarityCode","S"), binHdr->vibPolarityCode );
  env->SetShortField( binHdr_in, env->GetFieldID(class_binHdr,"revisionNum","S"), binHdr->revisionNum );
  env->SetShortField( binHdr_in, env->GetFieldID(class_binHdr,"fixedTraceLengthFlag","S"), binHdr->fixedTraceLengthFlag );
  env->SetShortField( binHdr_in, env->GetFieldID(class_binHdr,"numExtendedBlocks","S"), binHdr->numExtendedBlocks );
}

/*
 * Class:     cseis_jni_csNativeSegyReader
 * Method:    native_charHeader
 * Signature: (J)Ljava/lang/String;
 */
JNIEXPORT jstring JNICALL Java_cseis_jni_csNativeSegyReader_native_1charHeader
(JNIEnv *env, jobject obj, jlong ptr_in )
{
  csSegyReader* ptr = reinterpret_cast<csSegyReader*>(ptr_in);
  char const* charHdr = ptr->charHdrBlock();
  return (env)->NewStringUTF(charHdr);
}

/*
 * Class:     cseis_jni_csNativeSegyReader
 * Method:    native_hdrIntValue
 * Signature: (JI)I
 */
JNIEXPORT jint JNICALL Java_cseis_jni_csNativeSegyReader_native_1hdrIntValue
(JNIEnv *env, jobject obj, jlong ptr_in, jint hdrIndex )
{
  csSegyReader* ptr = reinterpret_cast<csSegyReader*>(ptr_in);
  return( ptr->getTraceHeader()->intValue( hdrIndex ) );
}

/*
 * Class:     cseis_jni_csNativeSegyReader
 * Method:    native_hdrFloatValue
 * Signature: (JI)F
 */
JNIEXPORT jfloat JNICALL Java_cseis_jni_csNativeSegyReader_native_1hdrFloatValue
(JNIEnv *env, jobject obj, jlong ptr_in, jint hdrIndex )
{
  csSegyReader* ptr = reinterpret_cast<csSegyReader*>(ptr_in);
  return( ptr->getTraceHeader()->floatValue( hdrIndex ) );
}

/*
 * Class:     cseis_jni_csNativeSegyReader
 * Method:    native_hdrDoubleValue
 * Signature: (JI)D
 */
JNIEXPORT jdouble JNICALL Java_cseis_jni_csNativeSegyReader_native_1hdrDoubleValue
(JNIEnv *env, jobject obj, jlong ptr_in, jint hdrIndex )
{
  csSegyReader* ptr = reinterpret_cast<csSegyReader*>(ptr_in);
  return( ptr->getTraceHeader()->doubleValue( hdrIndex ) );
}

/*
 * Class:     cseis_jni_csNativeSegyReader
 * Method:    native_moveToTrace
 * Signature: (JII)Z
 */
JNIEXPORT jboolean JNICALL Java_cseis_jni_csNativeSegyReader_native_1moveToTrace
(JNIEnv *env, jobject obj, jlong ptr_in, jint firstTraceIndex, jint numTracesToRead )
{
  csSegyReader* ptr = reinterpret_cast<csSegyReader*>(ptr_in);
  try {
    if( ptr->moveToTrace( firstTraceIndex, numTracesToRead ) ) {
      return JNI_TRUE;
    }
  }
  catch( csException& e ) {
    printf( "Error when moving to new position in SEGY file: %s\n", e.getMessage() );
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
 * Class:     cseis_jni_csNativeSegyReader
 * Method:    native_getNextTrace
 * Signature: (J[FLcseis/seis/csSeismicTrace;)Z
 */
JNIEXPORT jboolean JNICALL Java_cseis_jni_csNativeSegyReader_native_1getNextTrace
(JNIEnv *env, jobject obj, jlong ptr_in, jfloatArray samples_in, jobject trace_in )
{
  csSegyReader* ptr = reinterpret_cast<csSegyReader*>(ptr_in);

  float const* samplesPtr = NULL;
  try {
    samplesPtr = ptr->getNextTracePointer();
  }
  catch( csException& e ) {
    printf( "Error when reading from SEGY file: %s\n", e.getMessage() );
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

  jsize numSamples = ptr->binHdr()->numSamples;
  env->SetFloatArrayRegion( samples_in, 0, numSamples, (jfloat*)samplesPtr );

  int numHeaders = ptr->numTraceHeaders();

  csSegyTraceHeader const* trcHdr = ptr->getTraceHeader();

  for( int ihdr = 0; ihdr < numHeaders; ihdr++ ) {
    type_t type = ptr->header(ihdr)->outType;
    switch( type ) {
    case cseis_geolib::TYPE_DOUBLE:
      env->CallVoidMethod( trace_in, id_setSeismicTraceDoubleHeader, ihdr, trcHdr->doubleValue(ihdr) );
      break;
    case cseis_geolib::TYPE_FLOAT:
      env->CallVoidMethod( trace_in, id_setSeismicTraceFloatHeader, ihdr, trcHdr->floatValue(ihdr) );
      break;
    case cseis_geolib::TYPE_STRING:
      env->CallVoidMethod( trace_in, id_setSeismicTraceStringHeader, ihdr, (env)->NewStringUTF( trcHdr->stringValue(ihdr).c_str() ) );
      break;
     default: //       case csNumber::INT:
       env->CallVoidMethod( trace_in, id_setSeismicTraceIntHeader, ihdr, trcHdr->intValue(ihdr) );
       break;
    }
  }

  return JNI_TRUE;
}

/*
 * Class:     cseis_jni_csNativeSegyReader
 * Method:    native_headerName
 * Signature: (JI)Ljava/lang/String;
 */
JNIEXPORT jstring JNICALL Java_cseis_jni_csNativeSegyReader_native_1headerName
(JNIEnv *env, jobject obj, jlong ptr_in, jint headerIndex )
{
  csSegyReader* ptr = reinterpret_cast<csSegyReader*>(ptr_in);
  char const* name = ptr->header( headerIndex )->name.c_str();
  return (env)->NewStringUTF(name);
}

/*
 * Class:     cseis_jni_csNativeSegyReader
 * Method:    native_headerDesc
 * Signature: (JI)Ljava/lang/String;
 */
JNIEXPORT jstring JNICALL Java_cseis_jni_csNativeSegyReader_native_1headerDesc
(JNIEnv *env, jobject obj, jlong ptr_in, jint headerIndex )
{
  csSegyReader* ptr = reinterpret_cast<csSegyReader*>(ptr_in);
  char const* desc = ptr->header( headerIndex )->description.c_str();
  return (env)->NewStringUTF(desc);
}

/*
 * Class:     cseis_jni_csNativeSegyReader
 * Method:    native_headerType
 * Signature: (JI)I
 */
JNIEXPORT jint JNICALL Java_cseis_jni_csNativeSegyReader_native_1headerType
(JNIEnv *env, jobject obj, jlong ptr_in, jint headerIndex )
{
  csSegyReader* ptr = reinterpret_cast<csSegyReader*>(ptr_in);
  int headerType = (int)ptr->header( headerIndex )->outType;
  return headerType;
}

/*
 * Class:     cseis_jni_csNativeSegyReader
 * Method:    native_freeInstance
 * Signature: (J)V
 */
JNIEXPORT void JNICALL Java_cseis_jni_csNativeSegyReader_native_1freeInstance
(JNIEnv *env, jobject obj, jlong ptr_in) {
  csSegyReader* ptr = reinterpret_cast<csSegyReader*>(ptr_in);
  delete ptr;
}


/*
 * Class:     cseis_jni_csNativeSegyReader
 * Method:    native_trcHeaderMap
 * Signature: (ILcseis/jni/csSegyTrcHeaderDefinition;)V
 */
JNIEXPORT void JNICALL Java_cseis_jni_csNativeSegyReader_native_1trcHeaderMap
(JNIEnv *env, jclass cls, jint hdrMap, jobject hdrDef_in)
{
  csSegyHdrMap map( hdrMap );

  jclass class_csSegyHdrDef = env->GetObjectClass(hdrDef_in);
  jmethodID id_add = env->GetMethodID(class_csSegyHdrDef,"add","(Ljava/lang/String;Ljava/lang/String;II)V");


  int numHeaders = map.numHeaders();
  for( int ihdr = 0; ihdr < numHeaders; ihdr++ ) {
    cseis_geolib::csSegyHeaderInfo const* info = map.header( ihdr );
    env->CallVoidMethod( hdrDef_in, id_add,
         (env)->NewStringUTF(info->name.c_str()),
         (env)->NewStringUTF(info->description.c_str()),
         info->byteLoc, info->byteSize );
  }
}

/*
 * Class:     cseis_jni_csNativeSegyReader
 * Method:    native_setHeaderToPeek
 * Signature: (JLjava/lang/String;)V
 */
JNIEXPORT void JNICALL Java_cseis_jni_csNativeSegyReader_native_1setHeaderToPeek
(JNIEnv *env, jobject obj, jlong ptr_in, jstring headerName_in ) {
  csSegyReader* ptr = reinterpret_cast<csSegyReader*>(ptr_in);
  char const* headerName = (env)->GetStringUTFChars( headerName_in, NULL );
  std::string headerNameString( headerName );
  ptr->setHeaderToPeek( headerNameString );
}


/*
 * Class:     cseis_jni_csNativeSegyReader
 * Method:    native_peekHeaderValue
 * Signature: (JILcseis/seis/csHeader;)V
 */
JNIEXPORT void JNICALL Java_cseis_jni_csNativeSegyReader_native_1peekHeaderValue
(JNIEnv *env, jobject obj, jlong ptr_in, jint traceIndex_in, jobject value_in ) {
  csSegyReader* ptr = reinterpret_cast<csSegyReader*>(ptr_in);

  cseis_geolib::csFlexHeader flexValue;
  bool success = false;
  std::string message = "";
  try {
    success = ptr->peekHeaderValue( &flexValue, traceIndex_in );
  }
  catch( csException& e ) {
    message = std::string(e.getMessage());
    //    fprintf(stderr,"PEEK ERROR: %s\n", message.c_str());
  }
  if( !success ) {
    //    fprintf(stderr,"Peek operation failed... %d / %d\n", traceIndex_in, ptr->numTraces());
    //    fflush(stderr);
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
 * Class:     cseis_jni_csNativeSegyReader
 * Method:    native_setSelection
 * Signature: (JLjava/lang/String;Ljava/lang/String;IILcseis/jni/csITraceHeaderScanNotifier;)Z
 */
JNIEXPORT jboolean JNICALL Java_cseis_jni_csNativeSegyReader_native_1setSelection
(JNIEnv *env, jobject obj, jlong ptr_in, jstring headerValueSelectionText_in, jstring headerName_in, jint sortOrder, jint sortMethod, jobject notifier_in )
{
  csSegyReader* reader = reinterpret_cast<csSegyReader*>(ptr_in);
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
    //    fprintf( stderr, "\nError when setting trace selection/sorting: %s\n\n", e.getMessage() );
    //    fflush(stderr);
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
 * Class:     cseis_jni_csNativeSegyReader
 * Method:    native_getNumSelectedTraces
 * Signature: (J)I
 */
JNIEXPORT jint JNICALL Java_cseis_jni_csNativeSegyReader_native_1getNumSelectedTraces
(JNIEnv *env, jobject obj, jlong ptr_in ) {
  csSegyReader* reader = reinterpret_cast<csSegyReader*>(ptr_in);
  return reader->getNumSelectedTraces();
}

/*
 * Class:     cseis_jni_csNativeSegyReader
 * Method:    native_getSelectedValues
 * Signature: (JLcseis/jni/csSelectedHeaderBundle;)V
 */
JNIEXPORT void JNICALL Java_cseis_jni_csNativeSegyReader_native_1getSelectedValues
(JNIEnv *env, jobject obj, jlong ptr_in, jobject hdrBundle_in ) {
  csSegyReader* reader = reinterpret_cast<csSegyReader*>(ptr_in);

  jclass class_csSelectedHeaderBundle = env->GetObjectClass( hdrBundle_in );

  jmethodID id_setFloatHeader  = env->GetMethodID(class_csSelectedHeaderBundle,"setFloatValue","(IF)V");
  jmethodID id_setDoubleHeader = env->GetMethodID(class_csSelectedHeaderBundle,"setDoubleValue","(ID)V");
  jmethodID id_setIntHeader    = env->GetMethodID(class_csSelectedHeaderBundle,"setIntValue","(II)V");
  jmethodID id_setLongHeader   = env->GetMethodID(class_csSelectedHeaderBundle,"setLongValue","(IJ)V");
  jmethodID id_setTraceIndex   = env->GetMethodID(class_csSelectedHeaderBundle,"setSortedTraceIndex","(II)V");

  if( id_setFloatHeader <= 0  || id_setDoubleHeader <= 0  || id_setIntHeader <= 0  || id_setLongHeader <= 0 ) {
    //    fprintf(stderr,"Error in JNI layer, csNativeSeismicReader: getSelectedValues. Cannot find Java methods");
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
