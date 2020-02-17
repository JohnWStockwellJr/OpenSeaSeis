/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_jni_csNativeSegdReader.h"
#include "csSegdReader.h"
#include "csSegdHdrValues.h"
#include "csSegdDefines.h"
#include "csException.h"
#include "csVector.h"
#include <string>

using namespace std;

/*
 * Class:     cseis_jni_csNativeSegdReader
 * Method:    native_createInstance
 * Signature: (Ljava/lang/String;ZZZII)J
 */
JNIEXPORT jlong JNICALL Java_cseis_jni_csNativeSegdReader_native_1createInstance
(JNIEnv *env, jobject obj,
 jstring filename_in,
 jboolean readAuxTraces_in,
 jboolean isDebug_in,
 jboolean thisIsRev0_in,
 jint navInterfaceID_in,
 jint navSystemID_in )
{
  char const* filename = (env)->GetStringUTFChars( filename_in, NULL );

  cseis_segd::csSegdReader::configuration config;
  config.readAuxTraces     = true;
  config.isDebug           = false;
  config.navInterfaceID    = cseis_segd::UNKNOWN;
  config.numSamplesAddOne  = false;
  config.thisIsRev0        = false;
  config.navSystemID       = cseis_segd::NAV_HEADER_NONE;

  if( readAuxTraces_in == JNI_TRUE ) {
    config.readAuxTraces = true;
  }
  else {
    config.readAuxTraces = false;
  }
  if( isDebug_in == JNI_TRUE ) {
    config.isDebug = true;
  }
  else {
    config.isDebug = false;
  }
  if( thisIsRev0_in == JNI_TRUE ) {
    config.thisIsRev0 = true;
  }
  else {
    config.thisIsRev0 = false;
  }

  config.navSystemID    = navSystemID_in;
  config.navInterfaceID = navInterfaceID_in;

  cseis_segd::csSegdReader* reader = NULL;

  bool isNumSamplesError = false;
  int iteration = 0;

  while( iteration < 2 ) {
    iteration += 1;
    isNumSamplesError = false;
    if( iteration > 1 ) {
      delete reader;
      config.numSamplesAddOne = true;
    }
    reader = new cseis_segd::csSegdReader();
    reader->setConfiguration( config );

    try {
      reader->open( filename );
      reader->readNewRecordHeaders();
      cseis_segd::commonRecordHeaderStruct comRecordHdr;
      bool isSuccess = reader->readNextRecord( comRecordHdr );
      if( !isSuccess ) {
        return JNI_FALSE;
      }
      else {
      }
    }
    catch( cseis_geolib::csException& e ) {
      fprintf(stderr,"Error when opening SEGD file '%s'.\nSystem message: %s\n", filename, e.getMessage() );
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
    catch( string text ) {
      fprintf(stderr,"Error when opening SEGD file '%s'.\nSystem message: %s\n", filename, text.c_str() );
      fflush(stderr);
      jclass newExcCls = (env)->FindClass( "java/lang/Exception");
      if( newExcCls == NULL ) {
        return 0;
      }
      else {
        (env)->ThrowNew( newExcCls, text.c_str() );
        return 0;
      }
    }
    catch( ... ) {
      fprintf(stderr,"Unknown error occurred when opening SEGD file '%s'.\n", filename );
      fflush(stderr);
      jclass newExcCls = (env)->FindClass( "java/lang/Exception");
      if( newExcCls == NULL ) {
        return 0;
      }
      else {
        (env)->ThrowNew( newExcCls, "Unknown error occurred when opening SEGD file" );
        return 0;
      }
    }

    cseis_segd::commonChanSetStruct info;
    cseis_segd::commonFileHeaderStruct const* comFileHdr = reader->getCommonFileHeaders();
    float sampleInt = (float)comFileHdr->sampleInt_us / 1000.0;
    int numSamples = comFileHdr->numSamples;

    for( int i = 0; i < reader->numChanSets(); i++ ) {
      reader->retrieveChanSetInfo( i, info );
      if( info.numChannels == 0 ) continue;
      if( info.numSamples > numSamples ) {
        fprintf(stderr,"Number of samples for channel set #%d (=%d) exceeds number of samples specified in SEG-D general header or trace header extension (=%d).\nThis is not supported yet. Please read in a single channel set at a time (can only be done in a SeaSeis batch job).\n",
                i+1, info.numSamples, numSamples);
        fflush(stderr);
        jclass newExcCls = (env)->FindClass( "java/lang/Exception");
        if( newExcCls == NULL ) {
          return 0;
        }
        else {
          (env)->ThrowNew( newExcCls, "Number of samples for different channel sets differ. This is not supported yet.\nPlease read in a single channel set at a time (this can only be done in a SeaSeis batch job)." );
          return 0;
        }
      }
      else if( info.numSamples != numSamples ) {
        fprintf(stderr,"Number of samples for channel set #%d (=%d) differs from number of samples specified in SEG-D general header or trace header extension (=%d).\n   ...trying to add 1 to number of samples computed from general header.\n",
                i+1, info.numSamples, numSamples);
        isNumSamplesError = true;
      }
      if( comFileHdr->sampleInt_us != info.sampleInt_us ) {
        fprintf(stderr,"Sample interval for channel set #%d (=%.4fms) differs from sample interval specified in SEG-D general header (=%.4fms).\n",
                i+1, (float)(info.sampleInt_us/1000.0), sampleInt);
      }
    }
  } // END while

  if( isNumSamplesError ) {
    jclass newExcCls = (env)->FindClass( "java/lang/Exception");
    if( newExcCls == NULL ) {
      return 0;
    }
    else {
      (env)->ThrowNew( newExcCls, "Number of samples for in channel sets differ from number of samples specified in SEG-D general header or trace header extension.\nThis file cannot be read in with SeaView. Try to read it using batch module INPUT_SEGD in SeaSeis." );
      return 0;
    }
  }

  return ( reinterpret_cast<jlong>( reader ) );
}

/*
 * Class:     cseis_jni_csNativeSegdReader
 * Method:    native_closeFile
 * Signature: (J)V
 */
JNIEXPORT void JNICALL Java_cseis_jni_csNativeSegdReader_native_1closeFile
(JNIEnv *env, jobject obj, jlong ptr_in )
{
  cseis_segd::csSegdReader* ptr = reinterpret_cast<cseis_segd::csSegdReader*>(ptr_in);
  ptr->closeFile();
}

/*
 * Class:     cseis_jni_csNativeSegdReader
 * Method:    native_numTraces
 * Signature: (J)I
 */
JNIEXPORT jint JNICALL Java_cseis_jni_csNativeSegdReader_native_1numTraces
(JNIEnv *env, jobject obj, jlong ptr_in )
{
  cseis_segd::csSegdReader* ptr = reinterpret_cast<cseis_segd::csSegdReader*>(ptr_in);
  return( ptr->numTraces() );
}

/*
 * Class:     cseis_jni_csNativeSegdReader
 * Method:    native_numSamples
 * Signature: (J)I
 */
JNIEXPORT jint JNICALL Java_cseis_jni_csNativeSegdReader_native_1numSamples
(JNIEnv *env, jobject obj, jlong ptr_in )
{
  cseis_segd::csSegdReader* ptr = reinterpret_cast<cseis_segd::csSegdReader*>(ptr_in);
  return( ptr->getCommonFileHeaders()->numSamples );
}

/*
 * Class:     cseis_jni_csNativeSegdReader
 * Method:    native_numHeaders
 * Signature: (J)I
 */
JNIEXPORT jint JNICALL Java_cseis_jni_csNativeSegdReader_native_1numHeaders
(JNIEnv *env, jobject obj, jlong ptr_in )
{
  cseis_segd::csSegdReader* ptr = reinterpret_cast<cseis_segd::csSegdReader*>(ptr_in);
  return( ptr->hdrValues()->numHeaders() );
}

/*
 * Class:     cseis_jni_csNativeSegdReader
 * Method:    native_sampleInt
 * Signature: (J)F
 */
JNIEXPORT jfloat JNICALL Java_cseis_jni_csNativeSegdReader_native_1sampleInt
(JNIEnv *env, jobject obj, jlong ptr_in )
{
  cseis_segd::csSegdReader* ptr = reinterpret_cast<cseis_segd::csSegdReader*>(ptr_in);
  return( (jfloat)(ptr->getCommonFileHeaders()->sampleInt_us / 1000.0) );
}

/*
 * Class:     cseis_jni_csNativeSegdReader
 * Method:    native_hdrIntValue
 * Signature: (JI)I
 */
JNIEXPORT jint JNICALL Java_cseis_jni_csNativeSegdReader_native_1hdrIntValue
(JNIEnv *env, jobject obj, jlong ptr_in, jint hdrIndex )
{
  cseis_segd::csSegdReader* ptr = reinterpret_cast<cseis_segd::csSegdReader*>(ptr_in);
  return( ptr->hdrValues()->intValue(hdrIndex) );
}

/*
 * Class:     cseis_jni_csNativeSegdReader
 * Method:    native_hdrFloatValue
 * Signature: (JI)F
 */
JNIEXPORT jfloat JNICALL Java_cseis_jni_csNativeSegdReader_native_1hdrFloatValue
(JNIEnv *env, jobject obj, jlong ptr_in, jint hdrIndex )
{
  cseis_segd::csSegdReader* ptr = reinterpret_cast<cseis_segd::csSegdReader*>(ptr_in);
  return( ptr->hdrValues()->floatValue(hdrIndex) );
}

/*
 * Class:     cseis_jni_csNativeSegdReader
 * Method:    native_hdrDoubleValue
 * Signature: (JI)D
 */
JNIEXPORT jdouble JNICALL Java_cseis_jni_csNativeSegdReader_native_1hdrDoubleValue
(JNIEnv *env, jobject obj, jlong ptr_in, jint hdrIndex )
{
  cseis_segd::csSegdReader* ptr = reinterpret_cast<cseis_segd::csSegdReader*>(ptr_in);
  return( ptr->hdrValues()->doubleValue(hdrIndex) );
}

/*
 * Class:     cseis_jni_csNativeSegdReader
 * Method:    native_getNextTrace
 * Signature: (J[FLcseis/seis/csSeismicTrace;)Z
 */
JNIEXPORT jboolean JNICALL Java_cseis_jni_csNativeSegdReader_native_1getNextTrace
(JNIEnv *env, jobject obj, jlong ptr_in, jfloatArray samples_in, jobject trace_in )
{
  cseis_segd::csSegdReader* ptr = reinterpret_cast<cseis_segd::csSegdReader*>(ptr_in);

  // Get csSeismicTrace class
  jclass class_csSeismicTrace  = env->GetObjectClass(trace_in);
  jmethodID id_setFloatHeader  = env->GetMethodID(class_csSeismicTrace,"setFloatHeader","(IF)V");
  jmethodID id_setDoubleHeader = env->GetMethodID(class_csSeismicTrace,"setDoubleHeader","(ID)V");
  jmethodID id_setIntHeader    = env->GetMethodID(class_csSeismicTrace,"setIntHeader","(II)V");
  jmethodID id_setStringHeader = env->GetMethodID(class_csSeismicTrace,"setStringHeader","(ILjava/lang/String;)V");


  float const* samplesPtr = NULL;

  cseis_segd::commonTraceHeaderStruct comTrcHdr;
  while( (samplesPtr = ptr->getNextTracePointer(comTrcHdr)) == NULL ) {
    try {
      cseis_segd::commonRecordHeaderStruct comRecordHdr;
      bool isSuccess = ptr->readNextRecord( comRecordHdr );
      if( !isSuccess ) {
        return JNI_FALSE;
      }
      else {
        //        vars->recordCounter += 1;
        // vars->traceCounter = 0;
      }
    }
    catch( cseis_geolib::csException& exc ) {
      jclass newExcCls = (env)->FindClass( "java/lang/Exception");
      if( newExcCls != NULL ) {
        (env)->ThrowNew( newExcCls, exc.getMessage() );
      }
      return JNI_FALSE;
    }
    catch( ... ) {
      jclass newExcCls = (env)->FindClass( "java/lang/Exception");
      if( newExcCls != NULL ) {
        (env)->ThrowNew( newExcCls, "\nUnknown exception occurred during SEGD read operation...\nTerminated...\n" );
      }
      return JNI_FALSE;
    }
  }
  if( samplesPtr == NULL ) {
    return JNI_FALSE;
  }

  jsize numSamples = ptr->getCommonFileHeaders()->numSamples;
  env->SetFloatArrayRegion( samples_in, 0, numSamples, (jfloat*)samplesPtr );
  int numHeaders = ptr->hdrValues()->numHeaders();

  cseis_segd::csSegdHdrValues const* trcHdr = ptr->hdrValues();

  for( int ihdr = 0; ihdr < numHeaders; ihdr++ ) {
    cseis_geolib::type_t type = trcHdr->headerType(ihdr);
    switch( type ) {
    case cseis_geolib::TYPE_DOUBLE:
      env->CallVoidMethod( trace_in, id_setDoubleHeader, ihdr, trcHdr->doubleValue(ihdr) );
      break;
    case cseis_geolib::TYPE_FLOAT:
      env->CallVoidMethod( trace_in, id_setFloatHeader, ihdr, trcHdr->floatValue(ihdr) );
      break;
    case cseis_geolib::TYPE_STRING:
      env->CallVoidMethod( trace_in, id_setStringHeader, ihdr, (env)->NewStringUTF( trcHdr->stringValue(ihdr).c_str() ) );
      break;
    default: //       case csNumber::INT:
      env->CallVoidMethod( trace_in, id_setIntHeader, ihdr, trcHdr->intValue(ihdr) );
      break;
    }
  }
  return JNI_TRUE;
}


/*
 * Class:     cseis_jni_csNativeSegdReader
 * Method:    native_moveToTrace
 * Signature: (JII)Z
 */
JNIEXPORT jboolean JNICALL Java_cseis_jni_csNativeSegdReader_native_1moveToTrace
(JNIEnv *env, jobject obj, jlong ptr_in, jint firstTraceIndex, jint numTracesToRead )
{
  //  cseis_segd::csSegdReader* ptr = reinterpret_cast<cseis_segd::csSegdReader*>(ptr_in);
  /*
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
  */
  return JNI_FALSE;
}


/*
 * Class:     cseis_jni_csNativeSegdReader
 * Method:    native_headerName
 * Signature: (JI)Ljava/lang/String;
 */
JNIEXPORT jstring JNICALL Java_cseis_jni_csNativeSegdReader_native_1headerName
(JNIEnv *env, jobject obj, jlong ptr_in, jint hdrIndex )
{
  cseis_segd::csSegdReader* ptr = reinterpret_cast<cseis_segd::csSegdReader*>(ptr_in);
  char const* name = ptr->hdrValues()->headerName( hdrIndex );
  return (env)->NewStringUTF(name);
}

/*
 * Class:     cseis_jni_csNativeSegdReader
 * Method:    native_headerDesc
 * Signature: (JI)Ljava/lang/String;
 */
JNIEXPORT jstring JNICALL Java_cseis_jni_csNativeSegdReader_native_1headerDesc
(JNIEnv *env, jobject obj, jlong ptr_in, jint hdrIndex )
{
  cseis_segd::csSegdReader* ptr = reinterpret_cast<cseis_segd::csSegdReader*>(ptr_in);
  char const* desc = ptr->hdrValues()->headerDesc( hdrIndex );
  return (env)->NewStringUTF(desc);
}

/*
 * Class:     cseis_jni_csNativeSegdReader
 * Method:    native_headerType
 * Signature: (JI)I
 */
JNIEXPORT jint JNICALL Java_cseis_jni_csNativeSegdReader_native_1headerType
(JNIEnv *env, jobject obj, jlong ptr_in, jint hdrIndex )
{
  cseis_segd::csSegdReader* ptr = reinterpret_cast<cseis_segd::csSegdReader*>(ptr_in);
  return( ptr->hdrValues()->headerType(hdrIndex) );
}

/*
 * Class:     cseis_jni_csNativeSegdReader
 * Method:    native_freeInstance
 * Signature: (J)V
 */
JNIEXPORT void JNICALL Java_cseis_jni_csNativeSegdReader_native_1freeInstance
(JNIEnv *env, jobject obj, jlong ptr_in) {
  cseis_segd::csSegdReader* ptr = reinterpret_cast<cseis_segd::csSegdReader*>(ptr_in);
  delete ptr;
}

