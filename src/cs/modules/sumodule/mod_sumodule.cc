/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;


#ifdef PLATFORM_WINDOWS
void init_mod_sumodule_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  writer->error("SUMODULE is not supported on Windows platform");
}
void exec_mod_sumodule_(
                        csTraceGather* traceGather,
                        int* port,
                        int* numTrcToKeep,
                        csExecPhaseEnv* env,
                        csLogWriter* writer ) { }

void cleanup_mod_sumodule_( csExecPhaseEnv* env, csLogWriter* writer ) {
  // Nothing
}

// UNIX systems...:
#else

#include "csSegyTraceHeader.h"
#include "csSegyHdrMap.h"
#include "geolib_endian.h"
#include "geolib_string_utils.h"
#include <cstring>
#include <sstream>

extern "C" {
#include <unistd.h>
#include <fcntl.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <pthread.h>
#include <sys/types.h>
#include <errno.h>
#include <signal.h>
#include <poll.h>
}

void setHeaders( csSegyTraceHeader* segyTrcHdr, csSegyHdrMap* segyHdrMap, int* hdrIndexSegy, type_t* hdrTypeSegy,
                 int numSamples, float sampleInt,  csTraceHeaderDef* hdef, int scalarPolarity );
void mini_sleep( int millisec );
void* threadFunction( void *arg );
pthread_mutex_t mutex1 = PTHREAD_MUTEX_INITIALIZER;
int is_writeable(int fd);
pthread_mutex_t lock_x;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: SUMODULE
 *
 * @author Bjorn Olofsson
 */
namespace mod_sumodule {
  class ThreadParam {
  public:
    ThreadParam() {}
    void setIOFlag( int newIOFlag ) {
      pthread_mutex_lock( &mutex1 );
      io_flag = newIOFlag;
      // fprintf(stdout,"Set IO flag to %d,   %d\n", newIOFlag, fd_c2p);
      pthread_mutex_unlock( &mutex1 );
    }
    int getIOFlag() const {
      return io_flag;
    }
    FILE* file_c2p;
    int fd_c2p;
    int sizeBuffer;
    byte_t* buffer;
    int counter;
    cseis_geolib::csVector<cseis_geolib::byte_t*> bufferList;
    
  private:
    int io_flag;
  };

  struct VariableStruct {
    int*  hdrIndexSegy;
    type_t* hdrTypeSegy;
    pid_t pid;

    int fd_p2c;
    int fd_c2p;
    FILE* file_p2c;
    FILE* file_c2p;
    pthread_t* thread;
    ThreadParam* threadParam;

    bool swapEndian;

    int sizeRead;

    int traceCounter;
    int traceCounterRead;
    csSegyHdrMap* segyHdrMap;
    csSegyTraceHeader* segyTrcHdrWrite;
    csSegyTraceHeader* segyTrcHdrRead;

    //    byte_t* bufferIn;
    byte_t* bufferOut;
    int totalTraceSizeIn;
    int totalTraceSizeOut;

    std::string suModulePath;
    int outputCounter;
    int numSamplesIn;
    float sampleIntIn;
  };
  
  static int const READING  = -1;
  static int const BUFFER_READY = 1;
  static int const FINISHED = 0;
  static int COUNTER = 0;
}
using mod_sumodule::VariableStruct;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_sumodule_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  csTraceHeaderDef* hdef  = env->headerDef;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );
  
  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );

  vars->swapEndian = false;
  vars->hdrIndexSegy = NULL;
  vars->hdrTypeSegy  = NULL;
  vars->pid    = -1;
  vars->fd_p2c   = -1;
  vars->fd_c2p   = -1;
  vars->file_c2p  = NULL;
  vars->file_p2c  = NULL;
  vars->thread     = NULL;
  vars->threadParam = NULL;

  vars->sizeRead = 0;

  vars->traceCounter = 0;
  vars->traceCounterRead = 0;
  vars->segyHdrMap = NULL;
  vars->segyTrcHdrRead  = NULL;
  vars->segyTrcHdrWrite = NULL;
  vars->hdrIndexSegy = NULL;
  vars->hdrTypeSegy  = NULL;
  //  vars->bufferIn       = NULL;
  vars->bufferOut      = NULL;
  vars->totalTraceSizeIn  = 0;
  vars->totalTraceSizeOut = 0;
  vars->outputCounter = 0;
  vars->numSamplesIn = shdr->numSamples;
  vars->sampleIntIn = shdr->sampleInt;

  //--------------------------------------------------
  //
  std::string suBinDir("/opt/su/bin/");
  if( param->exists("cwproot") ) {
    param->getString("cwproot", &suBinDir);
    suBinDir = suBinDir + "/bin/";
  }
  else {
    char* cwpRoot = getenv("CWPROOT");
    if( strlen(cwpRoot) <= 0 ) {
      writer->error("Environment variable CWPROOT not set. Please specify SU root directory with user parameter 'cwproot'");
    }
    string cwpRoot_str(cwpRoot);
    suBinDir = cwpRoot_str + "/bin/";
  }
  //--------------------------------------------------
  //
  if( param->exists("endian_format") ) {
    std::string text;
    param->getString("endian_format",&text);
    if( !text.compare("big") ) {
      vars->swapEndian = isPlatformLittleEndian();
    }
    if( !text.compare("little") ) {
      vars->swapEndian = !isPlatformLittleEndian();
    }
    else {
      writer->error("Unknown option: %s", text.c_str());
    }
  }

  //--------------------------------------------------
  // Read in SU command & argument list
  //
  std::string suCommand;
  csVector<std::string> tokenList;
  param->getString("command", &suCommand);
  tokenize( suCommand.c_str(), tokenList, true );
  if( tokenList.size() < 2 ) {
    writer->error("No arguments given in SU command '%s', Arguments must be separated by spaces.", suCommand.c_str());
  }
  int numArguments = tokenList.size();  // Count module name as argument, needed to pass as argument for execv..!?
  char** argumentList = new char*[numArguments+1];
  for( int i = 0; i < numArguments; i++ ) {
    std::string token = tokenList.at(i);
    argumentList[i] = new char[token.size()+1]; // Add one element to store NULL pointer, also needed for execv to work..
    sprintf(argumentList[i],"%s",token.c_str());
    argumentList[i][token.size()] = '\0';
  }
  argumentList[numArguments] = NULL; // Argument list must be terminated by a NULL pointer for execv to work
  std::string suModulePath = suBinDir + tokenList.at(0);
  vars->suModulePath = suModulePath;

  for( int i = 1; i < numArguments; i++ ) {
    if( !strncmp(argumentList[i],"dt",2) ) {
      shdr->sampleInt = atof(&(argumentList[i][3])) * 1000.0;  // Convert from [s] to [ms]
      writer->line("Output sample interval: %f  (%s)\n", shdr->sampleInt, argumentList[i] );
    }
    if( !strncmp(argumentList[i],"nt",2) ) {
      shdr->numSamples = atoi(&(argumentList[i][3]));
      writer->line("Output number of samples: %d  (%s)\n", shdr->numSamples, argumentList[i] );
    }
  }

  //--------------------------------------------------
  // Set up child process
  //
  std::ostringstream pidStream;
  pidStream << getpid();    // Get process ID
  std::string tmpfile_template = string("tmp" + pidStream.str());
  if( param->exists("temp_dir") ) {
    string dirname_temp;
    param->getString( "temp_dir", &dirname_temp );
    tmpfile_template = string(dirname_temp + "/" + tmpfile_template);
  }
  else {
    tmpfile_template = string("/tmp/" + tmpfile_template);
  }
  char counter_str[6];
  sprintf(counter_str,"_%04d",mod_sumodule::COUNTER++);
  counter_str[5] = '\0';
  //  fprintf(stdout,"%s\n", counter_str);
  std::string filename_p2c(tmpfile_template + counter_str + ".p2c");
  std::string filename_c2p(tmpfile_template + counter_str + ".c2p");

  writer->line("Create/use FIFO files %s and %s", filename_p2c.c_str(), filename_c2p.c_str() );

  int status = mknod(filename_p2c.c_str(), S_IRUSR | S_IWUSR | S_IFIFO, 0);
  if( (status == -1) && (errno != EEXIST) ) {
    fprintf(stdout,"Error creating the named pipe\n");
    exit(-1);
  }
  status = mknod(filename_c2p.c_str(), S_IRUSR | S_IWUSR | S_IFIFO, 0);
  if( (status == -1) && (errno != EEXIST) ) {
    fprintf(stdout,"Error creating the named pipe\n");
    exit(-1);
  }

  vars->pid = fork();
  if( vars->pid == -1 ) {
    fprintf(stdout,"Unable to fork new process\n");
  }

  //--------------------------------------------------------------------------------
  //
  if( vars->pid == 0 ) { // Child process
    vars->fd_p2c = open( filename_p2c.c_str(), O_RDONLY );
    if( vars->fd_p2c < 0 ) { 
      fprintf(stdout,"Child: Error opening the named pipe: %d %d '%s'\n", vars->fd_p2c, errno, strerror(errno));
      exit(-1);
    }
    //    fprintf(stdout,"Child: Error code: %d %s\n", vars->fd_p2c, strerror(errno));
    vars->fd_c2p = open( filename_c2p.c_str(), O_WRONLY );
    if( vars->fd_c2p < 0 ) {
      fprintf(stdout,"Child: Error opening the named pipe: %d %d '%s'\n", vars->fd_c2p, errno, strerror(errno));
      exit(-1);
    }
    //    fprintf(stdout,"Child: Error code: %d %s\n", vars->fd_c2p, strerror(errno));

    dup2(vars->fd_p2c,fileno(stdin));    // Redirect standard input(0) to child 'read pipe'
    dup2(vars->fd_c2p,fileno(stdout));  // Redirect standard output(1) to child 'write pipe'
    close(vars->fd_p2c);
    close(vars->fd_c2p);

    // execv replaces child process with an external one
    if( execv( suModulePath.c_str(), argumentList ) < 0 )  { 
      fprintf(stdout,"External SU process: \n %s", suModulePath.c_str());
      for( int i = 1; i < numArguments; i++ ) { // Skip first element which is the SU module name
        fprintf(stdout," %s", argumentList[i]);
      }
      fprintf(stdout,"\n");
      fprintf(stdout,"External SU process failed, see above.\n");
      kill( getppid(), 9 );  // Kill parent process
      exit(-1);
    }
    // Child process is done. Will not continue from here on
    kill( getppid() , 9 );  // Kill parent process
    writer->error("Failed to run external process");
  }
  //--------------------------------------------------------------------------------
  //
  else { // Parent process
    vars->fd_p2c = open( filename_p2c.c_str(), O_WRONLY );
    if( vars->fd_p2c < 0 ) {
      fprintf(stdout,"Parent: Error opening the named pipe: %d %d '%s'  '%s'\n",
              vars->fd_p2c, errno, strerror(errno), filename_p2c.c_str());
      exit(-1);
    }
    if( edef->isDebug() )  fprintf(stdout,"Parent: Error code: %d %s\n", vars->fd_p2c, strerror(errno));
    //    fd_c2p = open( filename_c2p.c_str(), O_RDONLY );
    vars->file_c2p = fopen( filename_c2p.c_str(), "r");
    vars->fd_c2p = fileno( vars->file_c2p );
    if( vars->fd_c2p < 0 ) {
      fprintf(stdout,"Parent: Error opening the named pipe: %d %d '%s'\n", vars->fd_c2p, errno, strerror(errno));
      exit(-1);
    }
    if( edef->isDebug() ) fprintf(stdout,"Parent: Error code: %d %s\n", vars->fd_c2p, strerror(errno));

    writer->write("External SU process:\n  %s", suModulePath.c_str());
    for( int i = 1; i < numArguments; i++ ) { // Skip first element which is the SU module name
      writer->write(" %s", argumentList[i]);
    }
    writer->write("\n");
  }

  for( int i = 0; i < numArguments; i++ ) {
    delete [] argumentList[i];
  }
  delete [] argumentList;

  //--------------------------------------------------
  // Create objects, allocate memory
  //
  vars->segyHdrMap = new csSegyHdrMap(csSegyHdrMap::SEGY_SU,true);
  vars->segyTrcHdrRead  = new cseis_geolib::csSegyTraceHeader(vars->segyHdrMap);
  vars->segyTrcHdrWrite = new cseis_geolib::csSegyTraceHeader(vars->segyHdrMap);
  int numHeaders = vars->segyTrcHdrWrite->numHeaders();
  vars->hdrIndexSegy = new int[numHeaders];
  vars->hdrTypeSegy  = new type_t[numHeaders];

  vars->totalTraceSizeIn  = 240 + vars->numSamplesIn*4;
  vars->totalTraceSizeOut = 240 + shdr->numSamples*4;
  //  vars->bufferIn  = new byte_t[vars->totalTraceSize];
  vars->bufferOut = new byte_t[vars->totalTraceSizeIn];
  for( int i = 0; i < vars->totalTraceSizeIn; i++ ) {
    //    vars->bufferIn[i] = 0;
    vars->bufferOut[i] = 0;
  }

  setHeaders( vars->segyTrcHdrWrite, vars->segyHdrMap, vars->hdrIndexSegy, vars->hdrTypeSegy, vars->numSamplesIn, vars->sampleIntIn, hdef, 1 );
  setHeaders( vars->segyTrcHdrRead, vars->segyHdrMap, vars->hdrIndexSegy, vars->hdrTypeSegy, shdr->numSamples, shdr->sampleInt, hdef, -1 );

  vars->threadParam = new mod_sumodule::ThreadParam();
  vars->threadParam->file_c2p = vars->file_c2p;
  vars->threadParam->fd_c2p   = vars->fd_p2c;
  vars->threadParam->sizeBuffer = vars->totalTraceSizeOut;
  vars->threadParam->buffer   = new byte_t[vars->totalTraceSizeOut];
  vars->threadParam->setIOFlag( mod_sumodule::READING );
  vars->threadParam->counter = 0;

  vars->thread = new pthread_t();
  int thread_stat = pthread_create( vars->thread, NULL, threadFunction, vars->threadParam );
  if( thread_stat < 0 ) {
    fprintf(stdout,"Error when creating thread\n");
    exit(-1);
  }

}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_sumodule_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;
  csSuperHeader const* shdr = env->superHeader;
  csTraceHeaderDef const* hdef  = env->headerDef;


  // TEMP
  fd_set wio;
  FD_ZERO(&wio);
  FD_SET(vars->fd_p2c, &wio);   
  // TEMP

  if( edef->isDebug() )  fprintf(stdout,"Trace counter %d, pid = %d\n", vars->traceCounter, vars->pid);

  if( traceGather->numTraces() != 0 ) {
    vars->traceCounter += 1;
    csTrace* trace = traceGather->trace(0);

    // Prepare output trace header to be passed on to SU command
    csTraceHeader* trcHdr = trace->getTraceHeader();
    int nHeadersSegy = vars->segyTrcHdrWrite->numHeaders();
    for( int ihdr = 0; ihdr < nHeadersSegy; ihdr++ ) {
      int hdrIdOut = vars->hdrIndexSegy[ihdr];
      if( hdrIdOut < 0 ) continue;
      if( vars->hdrTypeSegy[ihdr] == TYPE_FLOAT ) {
        vars->segyTrcHdrWrite->setFloatValue( ihdr, trcHdr->floatValue( hdrIdOut )  );
      }
      else if( vars->hdrTypeSegy[ihdr] == TYPE_DOUBLE ) {
        vars->segyTrcHdrWrite->setDoubleValue( ihdr, trcHdr->doubleValue(hdrIdOut) );
      }
      else if( vars->hdrTypeSegy[ihdr] == TYPE_INT ) {
        vars->segyTrcHdrWrite->setIntValue( ihdr, trcHdr->intValue(hdrIdOut) );
      }
      else if( vars->hdrTypeSegy[ihdr] == TYPE_INT64 ) {
        vars->segyTrcHdrWrite->setIntValue( ihdr, (int)trcHdr->int64Value(hdrIdOut) );
      }
      else {
        vars->segyTrcHdrWrite->setStringValue( ihdr, trcHdr->stringValue(hdrIdOut) );
      }
    }
    vars->segyTrcHdrWrite->writeHeaderValues(vars->bufferOut,vars->swapEndian,true);
    float* samples = trace->getTraceSamples();
    memcpy( &vars->bufferOut[240], (byte_t*)samples, vars->numSamplesIn*sizeof(float));

    if( vars->swapEndian ) swapEndian4( (char*)(&vars->bufferOut[240]), vars->numSamplesIn*4 );
    traceGather->freeTrace(0);

    if( vars->threadParam->getIOFlag() != mod_sumodule::READING ) {
      vars->traceCounterRead += vars->threadParam->bufferList.size();
      if( edef->isDebug() ) {
        fprintf(stdout,"Parent: Read one more trace: %d of %d, %d %s\n", vars->traceCounterRead, vars->traceCounter, vars->fd_p2c, vars->suModulePath.c_str() );
      }
      pthread_mutex_lock(&lock_x);
      for( int i = 0; i < vars->threadParam->bufferList.size(); i++ ) {
        vars->outputCounter += 1;
        if( edef->isDebug() ) fprintf(stdout,"Output trace #%d\n", vars->outputCounter);
        csTrace* newTrace = traceGather->createTrace( hdef, shdr->numSamples );
        csTraceHeader* newTrcHdr = newTrace->getTraceHeader();
        float* newSamples = newTrace->getTraceSamples();
        byte_t* bufferPtr = vars->threadParam->bufferList.at(i);
        memcpy(newSamples,&bufferPtr[240],shdr->numSamples*sizeof(float));
        vars->segyTrcHdrRead->readHeaderValues( bufferPtr, vars->swapEndian, true );
        delete [] bufferPtr;
        int nHeaders = vars->segyTrcHdrRead->numHeaders();
        for( int ihdr = 0; ihdr < nHeaders; ihdr++ ) {
          int hdrIdOut = vars->hdrIndexSegy[ihdr];
          if( hdrIdOut < 0 ) continue;
          switch( vars->hdrTypeSegy[ihdr] ) {
          case TYPE_FLOAT:
            newTrcHdr->setFloatValue( hdrIdOut, vars->segyTrcHdrRead->floatValue(ihdr) );
            break;
          case TYPE_DOUBLE:
            newTrcHdr->setDoubleValue( hdrIdOut, vars->segyTrcHdrRead->doubleValue(ihdr) );
            break;
          case TYPE_INT:
            newTrcHdr->setIntValue( hdrIdOut, vars->segyTrcHdrRead->intValue(ihdr) );
            break;
          case TYPE_INT64:
            newTrcHdr->setInt64Value( hdrIdOut, (csInt64_t)vars->segyTrcHdrRead->intValue(ihdr) );
            break;
          case TYPE_STRING:
            newTrcHdr->setStringValue( hdrIdOut, vars->segyTrcHdrRead->stringValue(ihdr) );
            break;
          }
        }
      }
      vars->threadParam->bufferList.clear();
      pthread_mutex_unlock(&lock_x);
      // !!! TEMP !!!
      // if( vars->traceCounterRead == vars->traceCounter ) close(vars->fd_c2p);
      if( vars->threadParam->getIOFlag() != mod_sumodule::FINISHED ) vars->threadParam->setIOFlag( mod_sumodule::READING );
    }
    //    else {
    //   traceGather->freeTrace(0);
    // }

    struct timeval tv1;
    tv1.tv_sec  = 0;
    tv1.tv_usec = 10;
    if( edef->isDebug() ) fprintf(stdout,"START select %d...\n", vars->fd_p2c);
    int retVal1 = select( vars->fd_p2c+1, NULL, &wio, NULL, NULL );
    if( edef->isDebug() ) fprintf(stdout,"END   select %d...\n", vars->fd_p2c);
    if( retVal1 < 0 ) {
      fprintf(stdout,"Error in write select\n");
      exit(-1);
    }
    // Continue checking until p2c file descriptor is ready to write
    while( !FD_ISSET(vars->fd_p2c, &wio) ) {
      if( edef->isDebug() ) fprintf(stdout,"Parent: BAD Write 'select' returned %d\n", retVal1);
      retVal1 = select( vars->fd_p2c+1, NULL, &wio, NULL, &tv1 );
    }

    if( !is_writeable(vars->fd_p2c) ) {
      fprintf(stdout,"IS NOT WRITABLE\n");
      close(vars->fd_p2c);
      close(vars->fd_c2p);
      exit(-1);
    }
    if( edef->isDebug() ) fprintf(stdout,"START write %d...\n", vars->fd_p2c);
    int sizeWrite = (int)write( vars->fd_p2c, vars->bufferOut, vars->totalTraceSizeIn );
    if( sizeWrite == -1 ) {
      writer->error("Parent process write error\n");
    }
    if( edef->isDebug() ) fprintf(stdout,"Written data %d to pipe...\n", vars->traceCounter);
    // TEMP!!
    //    mini_sleep(10);  // Help external process empty output pipe... There's got to be a safer and faster way of doing this!

    edef->setTracesAreWaiting();
    if( edef->isDebug() ) fprintf(stdout,"Parent: Wrote %d elements. Total trace size: %d\n", sizeWrite, vars->totalTraceSizeIn);
  }

  //********************************************************************************
  //

  if( edef->isDebug() ) fprintf(stdout, "LAST CALL? %d, num traces: %d, %d %s\n", edef->isLastCall(), traceGather->numTraces(), vars->fd_p2c, vars->suModulePath.c_str() );
  if( edef->isLastCall() ) {
    if( edef->isDebug() ) fprintf(stdout, "*****************\n*****************\n******************\n****************\n");
    //    fflush(vars->file_p2c); // Has no effect
    close(vars->fd_p2c);
    if( edef->isDebug() ) fprintf(stdout, "Entering MAIN loop\n");
    while( vars->threadParam->getIOFlag() != mod_sumodule::FINISHED || vars->threadParam->bufferList.size() > 0 ) {
      if( vars->threadParam->getIOFlag() == mod_sumodule::BUFFER_READY || (vars->threadParam->getIOFlag() == mod_sumodule::FINISHED && vars->threadParam->bufferList.size() > 0) ) {
        vars->traceCounterRead += vars->threadParam->bufferList.size();
        if( edef->isDebug() ) fprintf(stdout,"Parent: Read one more trace: %d of %d\n", vars->traceCounterRead, vars->traceCounter );
        // !!! TEMP !!!
        //  if( vars->traceCounterRead == vars->traceCounter ) close(vars->fd_c2p);

        pthread_mutex_lock(&lock_x);
        for( int i = 0; i < vars->threadParam->bufferList.size(); i++ ) {
          vars->outputCounter += 1;
          if( edef->isDebug() ) fprintf(stdout,"LAST Output trace #%d\n", vars->outputCounter);
          csTrace* newTrace = traceGather->createTrace( hdef, shdr->numSamples );
          csTraceHeader* trcHdr = newTrace->getTraceHeader();
          float* samples = newTrace->getTraceSamples();
          byte_t* bufferPtr = vars->threadParam->bufferList.at(i);
          memcpy(samples,&bufferPtr[240],shdr->numSamples*sizeof(float));
          vars->segyTrcHdrRead->readHeaderValues( bufferPtr, vars->swapEndian, true );
          delete [] bufferPtr;

          int nHeaders = vars->segyTrcHdrRead->numHeaders();
          for( int ihdr = 0; ihdr < nHeaders; ihdr++ ) {
            int hdrIdOut = vars->hdrIndexSegy[ihdr];
            if( hdrIdOut < 0 ) continue;
            switch( vars->hdrTypeSegy[ihdr] ) {
            case TYPE_FLOAT:
              trcHdr->setFloatValue( hdrIdOut, vars->segyTrcHdrRead->floatValue(ihdr) );
              break;
            case TYPE_DOUBLE:
              trcHdr->setDoubleValue( hdrIdOut, vars->segyTrcHdrRead->doubleValue(ihdr) );
              break;
            case TYPE_INT:
              trcHdr->setIntValue( hdrIdOut, vars->segyTrcHdrRead->intValue(ihdr) );
              break;
            case TYPE_INT64:
              trcHdr->setInt64Value( hdrIdOut, (csInt64_t)vars->segyTrcHdrRead->intValue(ihdr) );
              break;
            case TYPE_STRING:
              trcHdr->setStringValue( hdrIdOut, vars->segyTrcHdrRead->stringValue(ihdr) );
              break;
            }
          }
        }
        vars->threadParam->bufferList.clear();
        pthread_mutex_unlock(&lock_x);
        if( vars->threadParam->getIOFlag() != mod_sumodule::FINISHED) vars->threadParam->setIOFlag( mod_sumodule::READING );
        // !!! TEMP !!!
        //if( vars->traceCounterRead >= 10 ) {
        // close(vars->fd_p2c);
        //vars->threadParam->getIOFlag() = mod_sumodule::FINISHED;
        //}
        // else {
        //}
      }
      else {
        if( edef->isDebug() ) fprintf(stdout,"--Main process: Start sleep %d, %d   %s\n", vars->threadParam->getIOFlag(), vars->fd_p2c, vars->suModulePath.c_str());
        mini_sleep(10);
        if( edef->isDebug() ) fprintf(stdout,"--Main process: Stop sleep %d, %d\n", vars->threadParam->getIOFlag(), vars->fd_p2c);
      }
    }
    // wait for thread to finish before continuing
    pthread_join( *vars->thread, NULL );
  }

  return;
}

void setHeaders( csSegyTraceHeader* segyTrcHdr, csSegyHdrMap* segyHdrMap, int* hdrIndexSegy, type_t* hdrTypeSegy,
                 int numSamples, float sampleInt, csTraceHeaderDef* hdef, int scalarPolarity ) {

  // Add mandatory trace headers, if defined in user specified, pre-defined, trace header map
  int num_id = 0;
  int id_nsamp        = num_id++;
  int id_sampint_us   = num_id++;
  int id_scalar_coord = num_id++;
  int id_scalar_elev  = num_id++;
  int id_scalar_stat  = num_id++;
  int id_fold         = num_id++;
  int id_fold_vert    = num_id++;
  int id_data_type    = num_id++;
  int id_trc_type     = num_id++;
  int id_unit_coord   = num_id++;
  int id_time_code    = num_id++;;
  csVector<string> mandatoryHdrNames(num_id);
  csVector<int>    mandatoryHdrValues(num_id);
  mandatoryHdrNames.allocate( num_id, "" );
  mandatoryHdrValues.allocate( num_id, 0 );
  // Set default values for mandatory SEGY output headers:
  mandatoryHdrNames.set("nsamp",       id_nsamp);        mandatoryHdrValues.set( numSamples, id_nsamp );
  mandatoryHdrNames.set("sampint_us",  id_sampint_us);   mandatoryHdrValues.set( (int)(sampleInt*1000), id_sampint_us );
  mandatoryHdrNames.set("scalar_coord", id_scalar_coord); mandatoryHdrValues.set( scalarPolarity * 100, id_scalar_coord );
  mandatoryHdrNames.set("scalar_elev", id_scalar_elev);  mandatoryHdrValues.set( scalarPolarity * 10, id_scalar_elev );
  mandatoryHdrNames.set("scalar_stat", id_scalar_stat);  mandatoryHdrValues.set( scalarPolarity * 1, id_scalar_stat );
  mandatoryHdrNames.set("fold",        id_fold);         mandatoryHdrValues.set( 1, id_fold );
  mandatoryHdrNames.set("fold_vert",   id_fold_vert);    mandatoryHdrValues.set( 1, id_fold_vert );
  mandatoryHdrNames.set("data_type",   id_data_type);    mandatoryHdrValues.set( 1, id_data_type );
  mandatoryHdrNames.set("trc_type",    id_trc_type);     mandatoryHdrValues.set( 1, id_trc_type );
  mandatoryHdrNames.set("unit_coord",  id_unit_coord);   mandatoryHdrValues.set( 1, id_unit_coord );
  mandatoryHdrNames.set("time_code",   id_time_code);    mandatoryHdrValues.set( 4, id_time_code );
  // Check if mandatory headers are defined in user selected pre-set header map
  // Each header that is not defined will be removed from the list of mandatory headers
  int hdrIndex = -1;
  for( int ihdr = mandatoryHdrNames.size()-1; ihdr >= 0; ihdr-- ) {
    string name = mandatoryHdrNames.at(ihdr);
    if( !segyHdrMap->contains(name,&hdrIndex) ) {
      mandatoryHdrNames.remove( ihdr );
      mandatoryHdrValues.remove( ihdr );
    }
  }

  int numHeaders = segyTrcHdr->numHeaders();
  for( int ihdr = 0; ihdr < numHeaders; ihdr++ ) {
    string const& name = segyTrcHdr->headerName(ihdr);
    //    writer->line(" #%2d  %s", ihdr+1, name.c_str());
    if( hdef->headerExists(name) ) {
      hdrIndexSegy[ihdr] = hdef->headerIndex( name );
      hdrTypeSegy[ihdr]  = hdef->headerType( name );
    }
    else {  // Must be mandatory header...
      hdrIndexSegy[ihdr] = -1;
    }
  }

  for( int ihdr = 0; ihdr < mandatoryHdrNames.size(); ihdr++ ) {
    string name = mandatoryHdrNames.at(ihdr);
    int hdrIndex = segyTrcHdr->headerIndex(name);
    //    if( hdrIndex < 0 ) writer->error("Program bug! Wrong mandatory header index....");
    if( hdef->headerExists(name) ) {
      hdrIndexSegy[hdrIndex] = -1;
    }
    segyTrcHdr->setIntValue( hdrIndex, mandatoryHdrValues.at(ihdr) );
  }

}

void mini_sleep( int millisec ) {
  struct timespec req={0},rem={0};
  time_t sec = (int)(millisec/1000);
  millisec = (int)(millisec-(sec*1000));
  req.tv_sec  = sec;
  req.tv_nsec = millisec*1000000L;
  nanosleep(&req,&rem);
}

//--------------------------------------------------------------------------------
//
void* threadFunction( void *arg ) {
  mod_sumodule::ThreadParam* ptr = (mod_sumodule::ThreadParam*)arg;

  // fprintf(stdout,"*** Started thread...\n");
  ptr->setIOFlag( mod_sumodule::READING );
  while( ptr->getIOFlag() != mod_sumodule::FINISHED ) {
    // fprintf(stdout,"Thread process starting... %d, %d\n", ptr->getIOFlag(), ptr->fd_c2p );
    //    if( ptr->getIOFlag() == mod_sumodule::READING ) {
    int sizeRead = (int)fread( ptr->buffer, 1, ptr->sizeBuffer, ptr->file_c2p );
    if( sizeRead <= 0 ) {
      // fprintf(stdout,"Thread process read error\n");
      ptr->setIOFlag( mod_sumodule::FINISHED );
      return NULL;
    }
    ptr->counter += 1;
    // fprintf(stdout,"Read in from SU pipe trace #%d\n", ptr->counter);
    pthread_mutex_lock( &lock_x );
    ptr->bufferList.insertEnd( ptr->buffer );
    ptr->buffer = new byte_t[ptr->sizeBuffer];
    ptr->setIOFlag( mod_sumodule::BUFFER_READY );
    pthread_mutex_unlock( &lock_x );
    // Set thread into sleep mode until main thread has copied the current buffer..?
    // --> doesn't work, can create gridlock
    // instead, copy buffer into vector for later pickup
    //    }
    // else {
    // fprintf(stdout,"Thread process: Start sleep... %d  %d\n", ptr->getIOFlag(), ptr->fd_c2p);
    // mini_sleep(10);
    //    }
    // fprintf(stdout,"Thread process running... %d\n", ptr->getIOFlag());
  }
  // fprintf(stdout,"*** Finished thread...\n");

  return NULL;
}

int is_writeable(int fd) {
  struct pollfd p;
  int ret;

  p.fd = fd;
  p.events = POLLOUT;

  ret = poll(&p, 1, 0);

  if (ret < 0) {
    fprintf(stdout, "poll failed\n");
    exit(EXIT_FAILURE);
  }
  // fprintf(stdout,"Checked poll OK\n");
  return p.revents & POLLOUT;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_sumodule_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_sumodule::VariableStruct* vars = reinterpret_cast<mod_sumodule::VariableStruct*>( env->execPhaseDef->variables() );
  //  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars->fd_c2p >= 0 )  close(vars->fd_p2c);
  if( vars->fd_p2c >= 0 )  close(vars->fd_c2p);
  if( vars->fd_p2c >= 0 )  close(vars->fd_c2p);
  if( vars->segyHdrMap != NULL ) {
    delete vars->segyHdrMap;
    vars->segyHdrMap = NULL;
  }
  if( vars->segyTrcHdrRead != NULL ) {
    delete vars->segyTrcHdrRead;
    vars->segyTrcHdrRead = NULL;
  }
  if( vars->segyTrcHdrWrite != NULL ) {
    delete vars->segyTrcHdrWrite;
    vars->segyTrcHdrWrite = NULL;
  }
  if( vars->hdrIndexSegy != NULL ) {
    delete [] vars->hdrIndexSegy;
    vars->hdrIndexSegy = NULL;
  }
  if( vars->hdrTypeSegy != NULL ) {
    delete [] vars->hdrTypeSegy;
    vars->hdrTypeSegy = NULL;
  }
  if( vars->bufferOut != NULL ) {
    delete [] vars->bufferOut;
    vars->bufferOut = NULL;
  }
  //    if( vars->bufferIn != NULL ) {
  //   delete [] vars->bufferIn;
  //  vars->bufferIn = NULL;
  // }
  delete vars; vars = NULL;
  // Wait until child process has finished, and clean up zombie processes...
  int status = 0;
  waitpid(-1, &status, WNOHANG); // clean up any children
}

#endif
// END: LINUX SYSTEM

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_sumodule_( csParamDef* pdef ) {
  pdef->setModule( "SUMODULE", "Generic wrapper for Seismic Unix (SU) module (EXPERIMENTAL)",
                   "Assumes available SU installation. This module is not guaranteed to work for all SU modules. Only one instance of $SUMODULE can be used in each Seaseis flow" );
  pdef->setVersion( 1, 0 );

  pdef->addParam( "command", "SU command to run", NUM_VALUES_FIXED, "This version only supports running a single SU module");
  pdef->addValue( "\"sugain scale=2.0\"", VALTYPE_STRING, "SU command, including space separated argument list" );

  pdef->addParam( "endian_format", "Endian format of data passed to SU", NUM_VALUES_FIXED, "Depends on whether local SU installation was compiled with BIG or LITTLE endian format" );
  pdef->addValue( "big", VALTYPE_OPTION );
  pdef->addOption( "big", "Pass data in big endian format." );
  pdef->addOption( "little", "Pass data in little endian format." );

  pdef->addParam( "cwproot", "SU root directory", NUM_VALUES_FIXED, "By default, the environment variable $CWPROOT will be queried" );
  pdef->addValue( "", VALTYPE_STRING, "SU root directory, full path name" );
}

//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_sumodule_( csExecPhaseEnv* env, csLogWriter* writer ) {
  //  mod_sumodule::VariableStruct* vars = reinterpret_cast<mod_sumodule::VariableStruct*>( env->execPhaseDef->variables() );
  //  csExecPhaseDef* edef = env->execPhaseDef;
  //  csSuperHeader const* shdr = env->superHeader;
  //  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

extern "C" void _params_mod_sumodule_( csParamDef* pdef ) {
  params_mod_sumodule_( pdef );
}
extern "C" void _init_mod_sumodule_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_sumodule_( param, env, writer );
}
extern "C" bool _start_exec_mod_sumodule_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_sumodule_( env, writer );
}
extern "C" void _exec_mod_sumodule_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_sumodule_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_sumodule_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_sumodule_( env, writer );
}
