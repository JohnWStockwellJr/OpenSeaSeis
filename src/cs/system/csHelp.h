/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_HELP_H
#define CS_HELP_H

#include <string>
#include <cstdio>

namespace cseis_geolib {
  template <typename T> class csVector;
}
namespace cseis_system {

  class csParamDef;

class csHelp {
 public:
  csHelp();
  csHelp( std::FILE* stream );
  ~csHelp();
  bool moduleTextHelp( std::string moduleNameIn );
  bool moduleTextHelp( std::string moduleNameIn, std::string versionString );
  bool fullHtmlHelp();
  void standardHeaderHelp();
  void moduleHtmlHelp( csParamDef const& pdef, FILE* stream ) const;
  void moduleHtmlHelp( csParamDef const& pdef, std::string& text ) const;
  void moduleHtmlExample( bool printHeader, csParamDef const& pdef, std::string& text ) const;
  void moduleExample( bool printHeader, csParamDef const& pdef, std::string& text ) const;
  void moduleHtmlListing( std::string& text ) const;
  void standardHeaderHtmlListing( std::string& text ) const;
 private:
  void moduleTextHelpShort( cseis_system::csParamDef const& pdef );
  void moduleTextHelpFull( cseis_system::csParamDef const& pdef );
  void moduleTextHelp( bool printHeader, cseis_system::csParamDef const& pdef ) const;
  void moduleHtmlHelp( csParamDef const& pdef ) const;
  void moduleHtmlExample( bool printHeader, csParamDef const& pdef ) const;

  std::FILE* myStream;
};

} // end namespace

#endif

