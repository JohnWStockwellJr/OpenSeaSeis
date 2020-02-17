/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csHelp.h"
#include "csVector.h"
#include <iostream>
#include <string>
#include <cstring>
#include "geolib_string_utils.h"
#include "csParamDef.h"
#include "csException.h"
#include "csMethodRetriever.h"
#include "csStandardHeaders.h"
#include "csHeaderInfo.h"
#include "csGeolibUtils.h"
#include "csSegyHdrMap.h"

using namespace cseis_system;

csHelp::csHelp() 
{
  myStream = stdout;
}

csHelp::csHelp( FILE* stream ) 
{
  myStream = stream;
}

csHelp::~csHelp() 
{
}

void csHelp::moduleTextHelpFull( csParamDef const& pdef ) 
{
  std::fprintf(myStream,"\n#****************************************************\n");
  csParamDescription const* valueDescriptor;
  csParamDescription const* optionDescriptor;
  csParamDescription const* const module = pdef.module();

  if( module == NULL ) {
    fprintf(myStream,"#  ---Module name not defined--- \n");
    return;
  }

  // Special for SU modules: Print out self doc (= original SU doc)
  if( strlen(module->name()) == 2 && !strncmp(module->name(),"SU",2) ) {
    //    fprintf(myStream," *** This seems to be an SU module ***\n");
    fprintf(myStream,"%s\n",pdef.docString().c_str());
    cseis_geolib::csSegyHdrMap suHdrMap( cseis_geolib::csSegyHdrMap::SEGY_SU_ONLY, true );
    fprintf(myStream,"\n **** SU trace header mapping ****\n");
    suHdrMap.dump( myStream );
  }

  fprintf(myStream,"# %s (%s)\n", module->desc(), pdef.versionString().c_str() );
  if( module->descExtra() ) fprintf(myStream,"# %s\n", module->descExtra());
  fprintf(myStream,"#\n# Parameter definitions (default values in <>) :\n");
  int nParams = pdef.numParameters();
  for( int ip = 0; ip < nParams; ip++ ) {
    if( pdef.param(ip)->desc() ) fprintf(myStream, "# %-20s : %s.\n", pdef.param(ip)->name(), pdef.param(ip)->desc() );
    else fprintf(myStream, "# %-20s\n", pdef.param(ip)->name() );
    if( pdef.param(ip)->descExtra() ) fprintf(myStream, "# %-20s : %s.\n", "", pdef.param(ip)->descExtra() );

    //    fprintf(myStream,"# Value definitions:");
    int nValues = pdef.numValues(ip);
    for( int iv = 0; iv < nValues; iv++ ) {
      valueDescriptor = pdef.value(ip,iv);
      if( strlen(valueDescriptor->name()) > 0 ) {
        fprintf(myStream, "# <%s>", valueDescriptor->name() );
      }
      else {
        fprintf(myStream, "# <?>" );
      }
      int fillLength = 18 - std::max( (int)strlen(valueDescriptor->name()), 1 );
      for( int fill = 0; fill < fillLength; fill++ ) {
        fprintf(myStream, " ");
      }

      int nOptions = pdef.numOptions(ip,iv);
      if( valueDescriptor->desc() ) {
        fprintf(myStream, " : %s.\n", valueDescriptor->desc() );
        if( nOptions > 0 ) fprintf(myStream, "#   %-18s   Options: ", "" );
      }
      else if( nOptions > 0 ) {
        fprintf(myStream, " : Options: " );
      }
      else {
        fprintf(myStream, "\n");
      }
      for( int io = 0; io < nOptions; io++ ) {
        optionDescriptor = pdef.option(ip,iv,io);
        fprintf(myStream, "%s (%s)", optionDescriptor->name(), optionDescriptor->desc() );
        if( optionDescriptor->descExtra() ) fprintf(myStream, " : %s.", optionDescriptor->descExtra() );
        if( io < nOptions-1 ) fprintf(myStream," / ");
      }
      if( nOptions > 0 ) fprintf(myStream,"\n");
      if( valueDescriptor->descExtra() ) fprintf(myStream, "#   %-18s : %s.\n", "", valueDescriptor->descExtra() );
    }

  }
  fprintf(myStream,"#\n");
  moduleTextHelp( false, pdef );
}

//--------------------------------------------------------------------------------
//
void csHelp::moduleTextHelpShort( csParamDef const& pdef ) 
{
  fprintf(myStream,"\n#****************************************************\n");
  csParamDescription const* const module = pdef.module();

  if( module == NULL ) {
    fprintf(myStream,"#  ---Module name not defined--- \n");
    return;
  }

  fprintf(myStream,"# %s\n", module->desc());
  if( module->descExtra() ) fprintf(myStream,"# %s\n", module->descExtra());
  int nParams = pdef.numParameters();
  for( int ip = 0; ip < nParams; ip++ ) {
    if( pdef.param(ip)->desc() ) fprintf(myStream, "# %-20s : %s.\n", pdef.param(ip)->name(), pdef.param(ip)->desc() );
    else fprintf(myStream, "# %-20s\n", pdef.param(ip)->name() );
    if( pdef.param(ip)->descExtra() ) fprintf(myStream, "# %-20s   %s.\n", "", pdef.param(ip)->descExtra() );
  }

  moduleTextHelp( false, pdef );
}
//--------------------------------------------------------------------------------
//
void csHelp::moduleTextHelp( bool printHeader, csParamDef const& pdef ) const 
{
  csParamDescription const* const module = pdef.module();
  if( printHeader ) {
    fprintf(myStream,"# %s (%s)\n", module->desc(), pdef.versionString().c_str() );
    if( module->descExtra() ) fprintf(myStream,"# %s\n", module->descExtra());
  }
  int nParams = pdef.numParameters();
  csParamDescription const* valueDescriptor;
  fprintf(myStream,"$%s\n", module->name() );
  for( int ip = 0; ip < nParams; ip++ ) {
    fprintf(myStream, " %-20s ", pdef.param(ip)->name() );
    int nValues = pdef.numValues(ip);
    for( int iv = 0; iv < nValues; iv++ ) {
      valueDescriptor = pdef.value(ip,iv);
      if( strlen(valueDescriptor->name()) != 0 ) fprintf(myStream, "%s  ", valueDescriptor->name() );
      else fprintf(myStream, "?  " );
    }
    fprintf(myStream, "\n" );
  }    
}
//--------------------------------------------------------------------------------
//

bool csHelp::moduleTextHelp( std::string moduleNameIn ) 
{
  return moduleTextHelp( moduleNameIn, "" );
}

bool csHelp::moduleTextHelp( std::string moduleNameIn, std::string versionString ) 
{
  std::string moduleName   = cseis_geolib::toUpperCase(moduleNameIn);
  int numModules      = csMethodRetriever::getNumStandardModules();
  std::string const* names = csMethodRetriever::getStandardModuleNames();
  csParamDef pdef;

  try {
    pdef.clear();
    MParamPtr funcPtr = NULL;
    if( versionString.length() == 0 ) {
      funcPtr = csMethodRetriever::getParamMethod( moduleNameIn );
    }
    else {
      funcPtr = csMethodRetriever::getParamMethod( moduleNameIn, versionString );
     }
    if( funcPtr != NULL ) {
      funcPtr( &pdef );
      moduleTextHelpFull( pdef );
      return true;
    }
  }
  catch( cseis_geolib::csException& exc ) {
    fprintf(myStream,"%s", exc.getMessage() );
    if( versionString.length() > 0 ) return false;
    fprintf(myStream,"...try to find module name starting with string '%s' in standard module list", moduleName.c_str() );
  }
  try {
    // Not found, is this ALL? otherwise try again to match first N letters
    if( !moduleName.compare(",") ) {
      for( int i = 0; i < numModules; i++ ) {
        pdef.clear();
        if( versionString.length() == 0 ) {
          csMethodRetriever::getParamMethod( names[i] )( &pdef );
        }
        else {
          csMethodRetriever::getParamMethod( names[i], versionString )( &pdef );
        }
        moduleTextHelpFull( pdef );
      }
      return true;
    }
    else if( !moduleName.compare(".") ) {
      for( int i = 0; i < numModules; i++ ) {
        pdef.clear();
        if( versionString.length() == 0 ) {
          csMethodRetriever::getParamMethod( names[i] )( &pdef );
        }
        else {
          csMethodRetriever::getParamMethod( names[i], versionString )( &pdef );
        }
        moduleTextHelpShort( pdef );
      }
      return true;
    }
    else {
      int length = moduleName.length();
      for( int i = 0; i < numModules; i++ ) {
        if( (int)(names[i].length()) >= length ) {
          if( !moduleName.compare( names[i].substr(0,length) ) ) {
            pdef.clear();
            if( versionString.length() == 0 ) {
              csMethodRetriever::getParamMethod( names[i] )( &pdef );
            }
            else {
              csMethodRetriever::getParamMethod( names[i], versionString )( &pdef );
            }
            moduleTextHelpFull( pdef );
            return true;
          }
        }
      }
    }
  }
  catch( cseis_geolib::csException& exc ) {
    fprintf(myStream,"%s", exc.getMessage() );
    return false;
  }

  fprintf(myStream,"\n  ...no help found\n\n");
  return false;
}
//----------------------------------------------------------------
//
//
bool csHelp::fullHtmlHelp() 
{
  int numModules = csMethodRetriever::getNumStandardModules();
  std::string const* names = csMethodRetriever::getStandardModuleNames();

  fprintf(myStream,"<html>\n");
  fprintf(myStream,"<body>\n");
  fprintf(myStream,"<h1>SeaSeis Help Page</h1>\n\n");
  fprintf(myStream, " <i>SeaSeis version %s, build date %s</i>\n", CSEIS_VERSION, __DATE__);

  //---------------------------------------------------------------------------------
  //
  fprintf(myStream,"<h2>Module List</h2>\n\n");
  fprintf(myStream,"<DL><DD>\n");
  fprintf(myStream,"<table border=\"0\" cellpadding=\"2\" width=\"100%%\">\n");
  fprintf(myStream,"<COLGROUP> <COL width=\"10%%\"> <COL width=\"90%%\">\n");
  for( int i = 0; i < numModules; i++ ) {
    csParamDef pdef;
    try {
      csMethodRetriever::getParamMethod( names[i] )( &pdef );
    }
    catch( cseis_geolib::csException& exc ) {
      fprintf(stderr,"Module '%s' not found...\n", names[i].c_str());
      continue;
    }
    csParamDescription const* const module = pdef.module();
    fprintf(myStream,"<tr>\n");
    fprintf(myStream,"<td> <a href=\"#%s\"><b>%s</b></a> </td>\n", module->name(), module->name() );
    fprintf(myStream,"<td> <font face=\"arial,helvetica,sans-serif\"> %s </font> </td>\n", module->desc());
    fprintf(myStream,"</tr>\n");
  }
  fprintf(myStream,"</table>\n");
  fprintf(myStream,"</DL>\n");

  //---------------------------------------------------------------------------------
  //
  fprintf(myStream,"<h2>Module Description</h2>\n\n");
  for( int i = 0; i < numModules; i++ ) {
    csParamDef pdef;
    try {
      csMethodRetriever::getParamMethod( names[i] )( &pdef );
    }
    catch( ... ) {
      continue;
    }
    moduleHtmlHelp( pdef );
  }

  //---------------------------------------------------------------------------------
  //
  fprintf(myStream,"<h2>Standard Header Definitions</h2>\n\n");
  fprintf(myStream,"<font face=\"arial,helvetica,sans-serif\">List of trace headers defined in SeaSeis.</font>\n");

  fprintf(myStream,"<DL><DD>\n");
  fprintf(myStream,"<table border=\"0\" cellpadding=\"2\" width=\"100%%\">\n");
  fprintf(myStream,"<COLGROUP> <COL width=\"4%%\"> <COL width=\"10%%\"> <COL width=\"6%%\"> <COL width=\"80%%\">\n");

  cseis_geolib::csVector<cseis_geolib::csHeaderInfo const*> hdrList(40);
  cseis_geolib::csStandardHeaders::getAll( &hdrList );
  for( int ihdr = 0; ihdr < hdrList.size(); ihdr++ ) {
    cseis_geolib::csHeaderInfo const* info = hdrList.at(ihdr);
    fprintf(myStream,"<tr>\n");
    fprintf(myStream,"<td><font face=\"arial,helvetica,sans-serif\">#%d</font></td>\n", ihdr+1 );
    fprintf(myStream,"<td><b>%s</b></td>\n", info->name.c_str() );
    fprintf(myStream,"<td>%s</td>\n", cseis_geolib::csGeolibUtils::typeText(info->type) );
    fprintf(myStream,"<td><font face=\"arial,helvetica,sans-serif\">%s</font></td>\n", info->description.c_str() );
    fprintf(myStream,"</tr>\n");
  }
  fprintf(myStream,"</table>\n");
  fprintf(myStream,"</DL>\n");

  fprintf(myStream,"</body>\n");
  fprintf(myStream,"</html>\n");

  return true;
}
//--------------------------------------------------------------------------------
//
void csHelp::moduleHtmlListing( std::string& text ) const {
  int numModules = csMethodRetriever::getNumStandardModules();
  std::string const* names = csMethodRetriever::getStandardModuleNames();
  //---------------------------------------------------------------------------------
  //
  text.append("<h2>Module List</h2>\n\n");
  text.append("<DD>\n");
  text.append("<table border=\"0\" cellpadding=\"2\" width=\"100%%\">\n");
  text.append("<COLGROUP> <COL width=\"5%%\"> <COL width=\"95%%\">\n");
  for( int i = 0; i < numModules; i++ ) {
    csParamDef pdef;
    try {
      csMethodRetriever::getParamMethod( names[i] )( &pdef );
    }
    catch( cseis_geolib::csException& exc ) {
      //      fprintf(stderr,"Module '%s' not found...\n", names[i].c_str());
      continue;
    }
    csParamDescription const* const module = pdef.module();
    text.append("<tr>\n");
    text = text + "<td> <a href=\"#" + module->name() + "\"><b>" + module->name() + "</b></a> </td>\n";
    text = text + "<td> <font face=\"arial,helvetica,sans-serif\"> " + module->desc() + " </font> </td>\n";
    text.append("</tr>\n");
  }
  text.append("</table>\n");
  text.append("\n");
}
//--------------------------------------------------------------------------------
//
void csHelp::standardHeaderHtmlListing( std::string& text ) const {
  char buffer[1024];
  //---------------------------------------------------------------------------------
  //
  text.append("<h2>Standard Header Definitions</h2>\n\n");

  text.append("<DD>\n");
  text.append("<table border=\"0\" cellpadding=\"2\" width=\"100%%\">\n");
  text.append("<COLGROUP> <COL width=\"1%%\"> <COL width=\"5%%\"> <COL width=\"6%%\"> <COL width=\"88%%\">\n");

  cseis_geolib::csVector<cseis_geolib::csHeaderInfo const*> hdrList(40);
  cseis_geolib::csStandardHeaders::getAll( &hdrList );
  for( int ihdr = 0; ihdr < hdrList.size(); ihdr++ ) {
    cseis_geolib::csHeaderInfo const* info = hdrList.at(ihdr);
    text.append("<tr>\n");
    buffer[0] = '\0'; sprintf(buffer,"<td><font face=\"arial,helvetica,sans-serif\">#%d</font></td>\n", ihdr+1 );
    text.append(buffer);
    text = text + "<td><b>" + info->name.c_str() + "</b></td>\n";
    text = text + "<td>" + cseis_geolib::csGeolibUtils::typeText(info->type) + "</td>\n";
    text = text + "<td><font face=\"arial,helvetica,sans-serif\">" + info->description.c_str() + "</font></td>\n";
    text.append("</tr>\n");
  }
  text.append("</table>\n");
  text.append("\n");
}
//--------------------------------------------------------------------------------
//
//
void csHelp::moduleHtmlExample( bool printHeader, csParamDef const& pdef ) const 
{
  csParamDescription const* const module = pdef.module();
  
  if( printHeader ) {
    fprintf(myStream,"<tr><td colspan=3><font face=\"Courier New,serif\" size=-1># %s</font></td></tr>\n", module->desc());
    //    if( module->descExtra() ) fprintf(myStream,"<tr><td colspan=3><font face=\"Courier New,serif\" size=-1> # %s</font></td></tr>\n", module->descExtra());
  }
  int nParams = pdef.numParameters();
  csParamDescription const* valueDescriptor;
  fprintf(myStream,"<tr><td colspan=2> <font face=\"Courier New,serif\" size=-1><b>$%s</b></font></td></tr>\n", module->name() );

  for( int ip = 0; ip < nParams; ip++ ) {
    fprintf(myStream,"<tr>\n");
    fprintf(myStream, "<td></td>\n");
    fprintf(myStream, "<td> <font face=\"Courier New,serif\" size=-1> %s </font> </td>\n", pdef.param(ip)->name() );
    fprintf(myStream, "<td> <font face=\"Courier New,serif\" size=-1>");
    int nValues = pdef.numValues(ip);
    for( int iv = 0; iv < nValues; iv++ ) {
      valueDescriptor = pdef.value(ip,iv);
      if( strlen(valueDescriptor->name()) != 0 ) fprintf(myStream, "%s  ", valueDescriptor->name() );
      else fprintf(myStream, "?  " );
    }
    fprintf(myStream, "</font></td>\n");
  }    

}
void csHelp::moduleHtmlExample( bool printHeader, csParamDef const& pdef, std::string& text ) const 
{
  csParamDescription const* const module = pdef.module();
  if( printHeader ) {
    text = text + "<tr><td colspan=3><font face=\"Courier New,serif\" size=-1># " + module->desc() + "</font></td></tr>\n";
  }
  int nParams = pdef.numParameters();
  csParamDescription const* valueDescriptor;
  text = text + "<tr><td colspan=2> <font face=\"Courier New,serif\" size=-1><b>$" + module->name() + "</b></font></td></tr>\n";

  for( int ip = 0; ip < nParams; ip++ ) {
    text.append("<tr>\n");
    text.append("<td></td>\n");
    text = text + "<td> <font face=\"Courier New,serif\" size=-1> " + pdef.param(ip)->name() + " </font> </td>\n";
    text.append("<td> <font face=\"Courier New,serif\" size=-1>");
    int nValues = pdef.numValues(ip);
    for( int iv = 0; iv < nValues; iv++ ) {
      valueDescriptor = pdef.value(ip,iv);
      if( strlen(valueDescriptor->name()) != 0 ) {
        text.append(valueDescriptor->name());
        text += "  ";
      }
      else {
        text.append("?  ");
      }
    }
    text.append("</font></td>\n");
  }
}
void csHelp::moduleExample( bool printHeader, csParamDef const& pdef, std::string& text ) const 
{
  csParamDescription const* const module = pdef.module();
  if( printHeader ) {
    text = text + "# " + module->desc() + "\n";
  }
  int nParams = pdef.numParameters();
  csParamDescription const* valueDescriptor;
  text = text + "$" + module->name() + "\n";

  int maxNumLetters = 0;
  for( int ip = 0; ip < nParams; ip++ ) {
    maxNumLetters = std::max( maxNumLetters, (int)strlen( pdef.param(ip)->name() ) );
  }
  maxNumLetters += 2;

  for( int ip = 0; ip < nParams; ip++ ) {
    text = text + " " + pdef.param(ip)->name() + " ";
    for( int il = (int)strlen( pdef.param(ip)->name() ) + 2; il < maxNumLetters; il++ ) {
      text += " ";
    }
    int nValues = pdef.numValues(ip);
    for( int iv = 0; iv < nValues; iv++ ) {
      valueDescriptor = pdef.value(ip,iv);
      if( strlen(valueDescriptor->name()) != 0 ) {
        text.append(valueDescriptor->name());
        text += "  ";
      }
      else {
        text.append("?  ");
      }
    }
    text.append("\n");
  }
}
//--------------------------------------------------------------------------------
//
//
void csHelp::moduleHtmlHelp( csParamDef const& pdef ) const 
{
  moduleHtmlHelp( pdef, myStream );
}
void csHelp::moduleHtmlHelp( csParamDef const& pdef, FILE* stream ) const 
{
  csParamDescription const* valueDescriptor;
  csParamDescription const* optionDescriptor;
  csParamDescription const* const module = pdef.module();

  if( module == NULL ) {
    fprintf(stream,"#  ---Module name not defined--- \n");
    return;
  }

  fprintf(stream,"<table border=\"0\" cellpadding=\"5\" bgcolor=\"#222222\" width=\"100%%\">\n");
  fprintf(stream,"<td bgcolor=\"#eeeeff\">\n");

  fprintf(stdout,"<a name=\"%s\"></a>\n", module->name());
  fprintf(stream,"<h2>%s (%s)</h2>\n", module->name(), pdef.versionString().c_str());
  fprintf(stdout,"<DL><DD><font face=\"arial,helvetica,sans-serif\"> \n");
  fprintf(stream,"%s.", module->desc());
  if( module->descExtra() ) fprintf(stream," %s.", module->descExtra());

  std::string selfDoc = pdef.docString();
  if( selfDoc.length() > 0 ) {
    fprintf(stream,"<DD><p>\n%s\n<p>", selfDoc.c_str());
  }
  fprintf(stream,"</font>\n");

  int nParams = pdef.numParameters();
  for( int ip = 0; ip < nParams; ip++ ) {
    int nValues = pdef.numValues(ip);
    if( pdef.numValues(ip) == 1 && pdef.param(ip)->type() == cseis_system::NUM_VALUES_VARIABLE ) {
      fprintf(stream, "<p> <b>%s</b> <i>(List)</i></p>\n", pdef.param(ip)->name() );
    }
    else {
      fprintf(stream, "<p> <b>%s</b> </p>\n", pdef.param(ip)->name() );
    }
    //if( nValues == 1 ) fprintf(stream, "<p> <b>%s</b>   <i>(takes 1 value)</i></p>\n", pdef.param(ip)->name() );
    //else fprintf(stream, "<p> <b>%s</b>   <i>(takes %d values)</i></p>\n", pdef.param(ip)->name(), nValues );
    fprintf(stream,"<font face=\"arial,helvetica,sans-serif\">");
    if( pdef.param(ip)->descExtra() ) fprintf(stream, "%s<br>%s\n", pdef.param(ip)->desc(), pdef.param(ip)->descExtra() );
    else if( pdef.param(ip)->desc() ) fprintf(stream, "%s\n", pdef.param(ip)->desc() );
    fprintf(stream,"</font>");

    fprintf(stream,"<table border=\"1\" cellpadding=\"5\" bgcolor=\"#eeeeaa\" width=\"100%%\">\n");
    fprintf(stream,"<COLGROUP> <COL width=\"10%%\"> <COL width=\"90%%\">\n");

    for( int iv = 0; iv < nValues; iv++ ) {
      fprintf(stream,"<tr>\n");
      valueDescriptor = pdef.value(ip,iv);
      if( strlen(valueDescriptor->name()) > 0 ) {
        fprintf(stream,"<td>%s</td>\n", valueDescriptor->name() );
      }
      else {
//??    fprintf(stream,"<td>?</td>\n", valueDescriptor->name() );
        fprintf(stream,"<td>?</td>\n");        
      }
      fprintf(stream, "<td>" );
      int nOptions = pdef.numOptions(ip,iv);
      if( nOptions > 0 ) {
        fprintf(stream,"<table border=\"0\" cellpadding=\"2\" align=\"left\" valign=\"top\">\n");
        for( int io = 0; io < nOptions; io++ ) {
          fprintf(stream,"<tr>\n");
          optionDescriptor = pdef.option(ip,iv,io);
          fprintf(stream, "<td valign=top><b>%s</b></td> <td valign=top><font face=\"arial,helvetica,sans-serif\">\n", optionDescriptor->name() );
          if( optionDescriptor->descExtra() ) fprintf(stream,"%s.\n%s.", optionDescriptor->desc(), optionDescriptor->descExtra() );
          else fprintf(stream,"%s.", optionDescriptor->desc() );
          fprintf(stream,"</font></td></tr>\n");
        }
        fprintf(stream,"</table>\n");
      }
      else if( valueDescriptor->desc() ) {
        fprintf(stream, "<font face=\"arial,helvetica,sans-serif\">%s.</font>\n", valueDescriptor->desc() );
      }
      if( valueDescriptor->descExtra() ) {
        fprintf(stream, "<font face=\"arial,helvetica,sans-serif\"> \n%s.</font>\n", valueDescriptor->descExtra() );
      }
      fprintf(stream, "</td>\n" );
      fprintf(stream,"</tr>\n");
    }
    fprintf(stream,"</td>\n");
    fprintf(stream,"</table>\n");
  }

  fprintf(stdout,"<p><h3><font face=\"arial,helvetica,sans-serif\"><b>Example:</b></font></p>\n<DD>\n");
  fprintf(stream,"<table border=\"0\" bgcolor=\"#ffffff\"width=\"100%%\">\n");
  fprintf(stream,"<COLGROUP> <COL width=\"1%%\"> <COL width=\"9%%\"> <COL width=\"90%%\">\n");
  moduleHtmlExample( true, pdef );
  fprintf(stdout,"\n");

  fprintf(stream,"</table>\n");

  fprintf(stdout,"</DL>\n");

  fprintf(stream,"</table>\n");
}

//--------------------------------------------------------------------------------

void csHelp::moduleHtmlHelp( csParamDef const& pdef, std::string& text ) const
{
  csParamDescription const* valueDescriptor;
  csParamDescription const* optionDescriptor;
  csParamDescription const* const module = pdef.module();

  if( module == NULL ) {
    text.append("#  ---Module name not defined--- \n");
    return;
  }
  text.append("<table border=\"0\" cellpadding=\"5\" bgcolor=\"#222222\" width=\"100%%\">\n");
  text.append("<td bgcolor=\"#eeeeff\">\n");

  text = text + "<h2>" + module->name() + " (" + pdef.versionString().c_str() + ")</h2>\n";
  text = text + module->desc() + ".";
  if( module->descExtra() ) {
    text = text + module->descExtra() + ".";
  }
  std::string selfDoc = pdef.docString();
  if( selfDoc.length() > 0 ) {
    text = text + "<DD><p>\n" + selfDoc.c_str() + "\n<p>";
  }
  text.append("</font>\n");

  int nParams = pdef.numParameters();
  for( int ip = 0; ip < nParams; ip++ ) {
    int nValues = pdef.numValues(ip);
    if( pdef.numValues(ip) == 1 && pdef.param(ip)->type() == cseis_system::NUM_VALUES_VARIABLE ) {
      text = text + "<p> <b>" + pdef.param(ip)->name() + "</b> <i>(List)</i></p>\n";
    }
    else {
      text = text + "<p> <b>" + pdef.param(ip)->name() + "</b> </p>\n";
    }
    text.append("<font face=\"arial,helvetica,sans-serif\">");
    if( pdef.param(ip)->descExtra() ) {
      text = text + pdef.param(ip)->desc() + "<br>" + pdef.param(ip)->descExtra() + "\n";
    }
    else if( pdef.param(ip)->desc() ) {
      text = text + pdef.param(ip)->desc() + "\n";
    }
    text.append("</font>");

    text.append("<table border=\"1\" cellpadding=\"5\" bgcolor=\"#eeeeaa\" width=\"100%%\">\n");
    text.append("<COLGROUP> <COL width=\"10%%\"> <COL width=\"90%%\">\n");

    for( int iv = 0; iv < nValues; iv++ ) {
      text.append("<tr>\n");
      valueDescriptor = pdef.value(ip,iv);
      if( strlen(valueDescriptor->name()) > 0 ) {
        text += "<td>";
        text += valueDescriptor->name();
        text += "</td>\n";
      }
      else {
        text.append("<td>?</td>\n");
      }
      text.append("<td>");
      int nOptions = pdef.numOptions(ip,iv);
      if( nOptions > 0 ) {
        text.append("<table border=\"0\" cellpadding=\"2\" align=\"left\" valign=\"top\">\n");
        for( int io = 0; io < nOptions; io++ ) {
          text.append("<tr>\n");
          optionDescriptor = pdef.option(ip,iv,io);
          text += "<td valign=top><b>";
          text += optionDescriptor->name();
          text += "</b></td> <td valign=top><font face=\"arial,helvetica,sans-serif\">\n";
          if( optionDescriptor->descExtra() ) {
            text += optionDescriptor->desc();
            text += ".\n";
            text += optionDescriptor->descExtra();
            text += ".";
          }
          else {
            text += optionDescriptor->desc();
            text += ".";
          }
          text.append("</font></td></tr>\n");
        }
        text.append("</table>\n");
      }
      else if( valueDescriptor->desc() ) {
        text += "<font face=\"arial,helvetica,sans-serif\">";
        text += valueDescriptor->desc();
        text += ".</font>\n";
      }
      if( valueDescriptor->descExtra() ) {
        text += "<font face=\"arial,helvetica,sans-serif\"> \n";
        text +=  valueDescriptor->descExtra();
        text += ".</font>\n";
      }
      text.append("</td>\n");
      text.append("</tr>\n");
    }
    text.append("</td>\n");
    text.append("</table>\n");
  }

  text.append("<table border=\"0\" bgcolor=\"#ffffff\"width=\"100%%\">\n");
  text.append("<COLGROUP> <COL width=\"1%%\"> <COL width=\"9%%\"> <COL width=\"90%%\">\n");
  moduleHtmlExample( true, pdef, text );

  text.append("</table>\n");
  text.append("</table>\n");
}
//--------------------------------------------------------------------------------
//
//
void csHelp::standardHeaderHelp() 
{
  cseis_geolib::csVector<cseis_geolib::csHeaderInfo const*> hdrList(40);
  cseis_geolib::csStandardHeaders::getAll( &hdrList );
  for( int ihdr = 0; ihdr < hdrList.size(); ihdr++ ) {
    cseis_geolib::csHeaderInfo const* info = hdrList.at(ihdr);
    fprintf(myStream,"Trace header #%3d: %-20s %-10s %s\n", ihdr+1, info->name.c_str(), cseis_geolib::csGeolibUtils::typeText(info->type), info->description.c_str() );
  }
}
