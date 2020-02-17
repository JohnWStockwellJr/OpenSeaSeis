/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csRSFHeader.h"
#include "geolib_defines.h"
#include <cstring>

using namespace std;
using namespace cseis_io;

csRSFHeader::csRSFHeader() {
  n1 = 0;
  n2 = 0;
  n3 = 0;
  n4 = 0;
  n5 = 0;
  o1 = 0; o2 = 0; o3 = 0; o4 = 0; o5 = 0;
  e1 = 0; e2 = 0; e3 = 0; e4 = 0; e5 = 0;
  d1 = 0; d2 = 0; d3 = 0; d4 = 0; d5=0;
  world_x1 = 0; world_x2 = 0; world_x3 = 0;
  world_y1 = 0; world_y2 = 0; world_y3 = 0;
  il1 = il2 = il3 = 0;
  xl1 = xl2 = xl3 = 0;
  ild = xld = 0;
  unit1 = csRSFHeader::SAMPLE_UNIT_UNKNOWN;
  unit2 = csRSFHeader::SAMPLE_UNIT_UNKNOWN;
  unit3 = csRSFHeader::SAMPLE_UNIT_UNKNOWN;
  unit4 = csRSFHeader::SAMPLE_UNIT_UNKNOWN;
  unit5 = csRSFHeader::SAMPLE_UNIT_UNKNOWN;
  sou_x = sou_y = sou_z = 0;
  domain1 = cseis_geolib::DOMAIN_XT;
  esize = 4;
  filename_bin = "";
  filename_bin_full_path = "";
  data_format = DATA_FORMAT_UNKNOWN;
}
void csRSFHeader::set( csRSFHeader const& hdr ) {
  n1 = hdr.n1;
  n2 = hdr.n2;
  n3 = hdr.n3;
  n4 = hdr.n4;
  n5 = hdr.n5;
  
  o1 = hdr.o1;
  o2 = hdr.o2;
  o3 = hdr.o3;
  o4 = hdr.o4;
  o5 = hdr.o5;
  
  e1 = hdr.e1;
  e2 = hdr.e2;
  e3 = hdr.e3;
  e4 = hdr.e4;
  e5 = hdr.e5;

  d1 = hdr.d1;
  d2 = hdr.d2;
  d3 = hdr.d3;
  d4 = hdr.d4;
  d5 = hdr.d5;

  world_x1 = hdr.world_x1;
  world_x2 = hdr.world_x2;
  world_x3 = hdr.world_x3;

  world_y1 = hdr.world_y1;
  world_y2 = hdr.world_y2;
  world_y3 = hdr.world_y3;

  il1 = hdr.il1;
  il2 = hdr.il2;
  il3 = hdr.il3;

  xl1 = hdr.xl1;
  xl2 = hdr.xl2;
  xl3 = hdr.xl3;

  ild = hdr.ild;
  xld = hdr.xld;

  sou_x = hdr.sou_x;
  sou_y = hdr.sou_y;
  sou_z = hdr.sou_z;

  unit1 = hdr.unit1;
  unit2 = hdr.unit2;
  unit3 = hdr.unit3;
  unit4 = hdr.unit4;
  unit5 = hdr.unit5;

  domain1 = hdr.domain1;

  data_format = hdr.data_format;
  esize = hdr.esize;
  filename_bin = hdr.filename_bin;
  filename_bin_full_path = hdr.filename_bin_full_path;
}
csRSFHeader::~csRSFHeader() {
}
bool csRSFHeader::setField( char const* name, char const* value ) {
  //  fprintf(stderr,"Set field: '%s' '%s'\n", name, value);
  bool found = false;
  if( !strcmp(name,"in") ) {
    int length = (int)strlen(value);
    if( length > 0 ) {
      if( value[0] == '\"' ) {
        value = &value[1];
        length -= 1;
      }
      if( value[length-1] == '\"' ) {
        length -= 1;
      }
      char* tmpStr = new char[length+1];
      memcpy(tmpStr,value,length);
      tmpStr[length] = '\0';
      filename_bin = tmpStr;
      delete [] tmpStr;
      found = true;
    }
  }
  else if( !strcmp(name,"data_format") ) {
    if( !strncmp(value,"\"native_float\"",14) || !strncmp(value,"native_float",12) ) {
      data_format = DATA_FORMAT_FLOAT;
      found = true;
    }
  }
  else if( !strcmp(name,"esize") ) {
    esize = atoi(value); found = true;
  }
  else if( !strcmp(name,"n1") ) {
    n1 = atoi(value); found = true;
  }
  else if( !strcmp(name,"n2") ) {
    n2 = atoi(value); found = true;
  }
  else if( !strcmp(name,"n3") ) {
    n3 = atoi(value); found = true;
  }
  else if( !strcmp(name,"n4") ) {
    n4 = atoi(value); found = true;
  }
  else if( !strcmp(name,"n5") ) {
    n5 = atoi(value); found = true;
  }
  else if( !strcmp(name,"o1") ) {
    o1 = atof(value); found = true;
  }
  else if( !strcmp(name,"o2") ) {
    o2 = atof(value); found = true;
  }
  else if( !strcmp(name,"o3") ) {
    o3 = atof(value); found = true;
  }
  else if( !strcmp(name,"o4") ) {
    o4 = atof(value); found = true;
  }
  else if( !strcmp(name,"o5") ) {
    o5 = atof(value); found = true;
  }
  else if( !strcmp(name,"e1") ) {
    e1 = atof(value); found = true;
  }
  else if( !strcmp(name,"e2") ) {
    e2 = atof(value); found = true;
  }
  else if( !strcmp(name,"e3") ) {
    e3 = atof(value); found = true;
  }
  else if( !strcmp(name,"e4") ) {
    e4 = atof(value); found = true;
  }
  else if( !strcmp(name,"e5") ) {
    e5 = atof(value); found = true;
  }
  else if( !strcmp(name,"d1") ) {
    d1 = atof(value); found = true;
  }
  else if( !strcmp(name,"d2") ) {
    d2 = atof(value); found = true;
  }
  else if( !strcmp(name,"d3") ) {
    d3 = atof(value); found = true;
  }
  else if( !strcmp(name,"d4") ) {
    d4 = atof(value); found = true;
  }
  else if( !strcmp(name,"d5") ) {
    d5 = atof(value); found = true;
  }
  else if( !strcmp(name,"world_x1") ) {
    world_x1 = atof(value); found = true;
  }
  else if( !strcmp(name,"world_x2") ) {
    world_x2 = atof(value); found = true;
  }
  else if( !strcmp(name,"world_x3") ) {
    world_x3 = atof(value); found = true;
  }
  else if( !strcmp(name,"world_y1") ) {
    world_y1 = atof(value); found = true;
  }
  else if( !strcmp(name,"world_y2") ) {
    world_y2 = atof(value); found = true;
  }
  else if( !strcmp(name,"world_y3") ) {
    world_y3 = atof(value); found = true;
  }
  else if( !strcmp(name,"src_x") ) {
    sou_x = atof(value); found = true;
  }
  else if( !strcmp(name,"src_y") ) {
    sou_y = atof(value); found = true;
  }
  else if( !strcmp(name,"src_z") ) {
    sou_z = atof(value); found = true;
  }
  else if( !strcmp(name,"il1") ) {
    il1 = atof(value); found = true;
  }
  else if( !strcmp(name,"il2") ) {
    il2 = atof(value); found = true;
  }
  else if( !strcmp(name,"il3") ) {
    il3 = atof(value); found = true;
  }
  else if( !strcmp(name,"xl1") ) {
    xl1 = atof(value); found = true;
  }
  else if( !strcmp(name,"xl2") ) {
    xl2 = atof(value); found = true;
  }
  else if( !strcmp(name,"xl3") ) {
    xl3 = atof(value); found = true;
  }
  else if( !strcmp(name,"ild") ) {
    ild = atof(value); found = true;
  }
  else if( !strcmp(name,"xld") ) {
    xld = atof(value); found = true;
  }
  else if( !strncmp(name,"unit",4) ) {
    found = true;
    int unit = csRSFHeader::readSampleUnit(value);
    if( !strcmp(name,"unit1") ) {
      unit1 = unit;
    }
    else if( !strcmp(name,"unit2") ) {
      unit2 = unit;
    }
    else if( !strcmp(name,"unit3") ) {
      unit3 = unit;
    }
    else if( !strcmp(name,"unit4") ) {
      unit4 = unit;
    }
    else if( !strcmp(name,"unit5") ) {
      unit5 = unit;
    }
  }
  else if( !strncmp(name,"label1",6) ) {
    found = true;
    if( !strcmp(value,"depth") ) {
      domain1 = cseis_geolib::DOMAIN_XD;
    }
  }
  return found;
}
void csRSFHeader::dump() {
  dump( stdout );
}
void csRSFHeader::dump( FILE* stream ) {
  dump( stream, false );
}
void csRSFHeader::dump( FILE* stream, bool outputGrid ) {
  fprintf(stream,"in=\"%s\"\n", filename_bin.c_str());
  if( data_format == csRSFHeader::DATA_FORMAT_UNKNOWN ) {
    fprintf(stream,"data_format=\"none\"\n");
  }
  else {
    fprintf(stream,"data_format=\"native_float\"\n");
  }
  fprintf(stream,"esize=%d\n", esize);
  
  fprintf(stream,"n1=%d\n",n1);
  fprintf(stream,"d1=%f\n",d1);
  fprintf(stream,"o1=%f\n",o1);
  fprintf(stream,"e1=%f\n",e1);
  if( unit1 != csRSFHeader::SAMPLE_UNIT_UNKNOWN ) {
    fprintf(stream,"unit1=%d\n",unit1);
  }

  fprintf(stream,"n2=%d\n",n2);
  fprintf(stream,"d2=%f\n",d2);
  fprintf(stream,"o2=%f\n",o2);
  fprintf(stream,"e2=%f\n",e2);
  if( unit2 != csRSFHeader::SAMPLE_UNIT_UNKNOWN ) {
    fprintf(stream,"unit2=%d\n",unit2);
  }

  fprintf(stream,"n3=%d\n",n3);
  fprintf(stream,"d3=%f\n",d3);
  fprintf(stream,"o3=%f\n",o3);
  fprintf(stream,"e3=%f\n",e3);
  if( unit3 != csRSFHeader::SAMPLE_UNIT_UNKNOWN ) {
    fprintf(stream,"unit3=%d\n",unit3);
  }

  if( n4 > 0 ) {
    fprintf(stream,"n4=%d\n",n4);
    fprintf(stream,"d4=%f\n",d4);
    fprintf(stream,"o4=%f\n",o4);
    fprintf(stream,"e4=%f\n",e4);
    if( unit4 != csRSFHeader::SAMPLE_UNIT_UNKNOWN ) {
      fprintf(stream,"unit4=%d\n",unit4);
    }
  }

  if( n5 > 0 ) {
    fprintf(stream,"n5=%d\n",n5);
    fprintf(stream,"d5=%f\n",d5);
    fprintf(stream,"o5=%f\n",o5);
    fprintf(stream,"e5=%f\n",e5);
    if( unit5 != csRSFHeader::SAMPLE_UNIT_UNKNOWN ) {
      fprintf(stream,"unit5=%d\n",unit5);
    }
  }


  if( outputGrid ) {
    fprintf(stream,"xl1=%f\n",xl1);
    fprintf(stream,"il1=%f\n",il1);

    fprintf(stream,"xl2=%f\n",xl2);
    fprintf(stream,"il2=%f\n",il2);

    fprintf(stream,"xl3=%f\n",xl3);
    fprintf(stream,"il3=%f\n",il3);

    fprintf(stream,"xld=%f\n",xld);
    fprintf(stream,"ild=%f\n",ild);

    fprintf(stream,"world_x1=%f\n",world_x1);
    fprintf(stream,"world_y1=%f\n",world_y1);

    fprintf(stream,"world_x2=%f\n",world_x2);
    fprintf(stream,"world_y2=%f\n",world_y2);

    fprintf(stream,"world_x3=%f\n",world_x3);
    fprintf(stream,"world_y3=%f\n",world_y3);
  }
}
bool csRSFHeader::readRSFHeaderFile( char const* filename_rsf, csRSFHeader& hdr ) {
  FILE* fin_rsf = fopen(filename_rsf,"r");
  if( fin_rsf == NULL ) {
    fprintf(stderr,"Could not open input file: %s\n", filename_rsf);
    return false;
  }
  char buffer[132];
  while( fgets(buffer,132,fin_rsf) != NULL ) {
    char* cPtr = strtok(buffer,"=");
    if( cPtr != NULL ) {
      char* name = cPtr;
      cPtr = strtok(NULL,"=");
      if( cPtr != NULL ) {
        char* value = cPtr;
        hdr.setField( name, value );
      }
    }
  }

  // Set binary file name with full path
  int counter = strlen(filename_rsf)-1;
  while( counter > 0 ) {
    if( filename_rsf[counter] == '/' ) {
      break;
    }
    counter -= 1;
  }
  if( counter == 0 ) {
    counter = -1;
  }

  int sizeFile = hdr.filename_bin.length();
  int sizePath = counter + 1;
  char* tmpStr = new char[sizePath + sizeFile + 1];
  memcpy(tmpStr,filename_rsf,sizePath);
  memcpy(&tmpStr[sizePath],hdr.filename_bin.c_str(),sizeFile);
  tmpStr[sizePath+sizeFile] = '\0';
  hdr.filename_bin_full_path = tmpStr;
  delete [] tmpStr;

  fclose(fin_rsf);
  return true;
}
void csRSFHeader::setFilename( char const* theFilename_bin, char const* theFilename_bin_full_path ) {
  filename_bin = theFilename_bin;
  filename_bin = theFilename_bin_full_path;
}

int csRSFHeader::readSampleUnit( std::string text ) {
  if( !text.compare("s") ) {
    return csRSFHeader::SAMPLE_UNIT_S;
  }
  else if( !text.compare("ms") ) {
    return csRSFHeader::SAMPLE_UNIT_MS;
  }
  else if( !text.compare("m") ) {
    return csRSFHeader::SAMPLE_UNIT_M;
  }
  else if( !text.compare("km") ) {
    return csRSFHeader::SAMPLE_UNIT_KM;
  }
  else if( !text.compare("hz") || !text.compare("Hz") ) {
    return csRSFHeader::SAMPLE_UNIT_HZ;
  }
  else {
    return csRSFHeader::SAMPLE_UNIT_UNKNOWN;
  }
}
