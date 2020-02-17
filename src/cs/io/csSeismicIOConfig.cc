/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csSeismicIOConfig.h"
#include "geolib_defines.h"
#include "csHeaderInfo.h"

using namespace cseis_io;

csSeismicIOConfig::csSeismicIOConfig() {
  numSamples = 0;
  sampleInt  = 0;
  domain     = 0;
  fftDataType = 0;
  numSamplesXT = 0;
  sampleIntXT = 0;
  grid_orig_x  = 0.0;
  grid_orig_y  = 0.0;
  grid_orig_il = 0;
  grid_orig_xl = 0;
  grid_binsize_il = 0.0;
  grid_binsize_xl = 0.0;
  grid_azim_il = 0.0;
  grid_azim_xl = 0.0;
}

csSeismicIOConfig::~csSeismicIOConfig() {
  for( int ihdr = 0; ihdr < myTrcHdrList.size(); ++ihdr ) {
    delete myTrcHdrList.at(ihdr);
  }
  myTrcHdrList.clear();
}
int csSeismicIOConfig::numTrcHeaders() const {
  return myTrcHdrList.size();
}

void csSeismicIOConfig::addHeader( cseis_geolib::csHeaderInfo const* info ) {
  addHeader( info->type, info->name, info->description, info->nElements );
}
void csSeismicIOConfig::addHeader( cseis_geolib::type_t type, std::string const& name, std::string const& description, int nElements ) {
  cseis_geolib::csHeaderInfo* newInfo = new cseis_geolib::csHeaderInfo( type, name, description, nElements );
  myTrcHdrList.insertEnd( newInfo );
}

cseis_geolib::csHeaderInfo const* csSeismicIOConfig::headerInfo( int headerIndex ) const {
  return myTrcHdrList.at( headerIndex );
}

cseis_geolib::csHeaderInfo const* csSeismicIOConfig::headerInfo( std::string const& headerName ) const {
  for (int idx = 0 ; idx < myTrcHdrList.size(); ++idx ) {
    if (myTrcHdrList.at(idx)->name == headerName ) { 
      return myTrcHdrList.at(idx);
    }   
  }

  return NULL;
}

int csSeismicIOConfig::headerIndex( std::string const& headerName ) const {
  for (int idx = 0 ; idx < myTrcHdrList.size(); ++idx ) {
    if (myTrcHdrList.at(idx)->name == headerName ) { 
      return idx;
    }
  }

  return -1;
}


