/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csSeismicWriter.h"
#include "csFileUtils.h"
#include <sstream>
#include <ctime>

extern "C" {
  #include <unistd.h>
}

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: IMAGE
 *
 * @author Bjorn Olofsson
 * @date   2009
 */
namespace mod_image {
  struct VariableStruct {
    std::string propertiesFilename;
    std::string seismicFilename;
    std::string command;
    bool keepPropertyFile;
    bool isTmpFile;
    cseis_system::csSeismicWriter* writer;
    int nTracesOut;
    bool isFirstCall;
  };

   static const int WIGGLE_FILL_NONE   = 0;
   static const int WIGGLE_FILL_POS    = 1;
   static const int WIGGLE_FILL_NEG    = 2;
   static const int WIGGLE_FILL_POSNEG = 3;
  
   static const int WIGGLE_COLOR_FIXED    = 41;
   static const int WIGGLE_COLOR_VARIABLE = 42;

   static const int SCALE_TYPE_SCALAR = 31;
   static const int SCALE_TYPE_RANGE  = 32;
   static const int SCALE_TYPE_TRACE  = 33;

   static const int TRACE_SCALING_MAXIMUM = 51;
   static const int TRACE_SCALING_AVERAGE = 52;
  
   static const int WIGGLE_TYPE_LINEAR = 11;
   static const int WIGGLE_TYPE_CUBIC  = 12;

   static const float POLARITY_NORMAL    = 1.0f;
   static const float POLARITY_REVERSED  = -1.0f;

   static const int VA_TYPE_2DSPLINE = 21;
   static const int VA_TYPE_VERTICAL = 22;
   static const int VA_TYPE_DISCRETE = 23;

   static const int PLOT_DIR_VERTICAL   = 0;
   static const int PLOT_DIR_HORIZONTAL = 1;
}
using namespace mod_image;

std::string random_name( std::string pre, int size, std::string post );

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_image_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  //  csTraceHeaderDef* hdef = env->headerDef;
  //  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );


#if PLATFORM_WINDOWS
  // For Windows, add double quotation marks around plotimage command line arguments
  // This is because Windows batch files only recognise maximum 9 command line arguments
  vars->command            = "plotimage.bat \"";
#else
  vars->command            = "plotimage.sh";
#endif
  vars->seismicFilename    = "";
  vars->keepPropertyFile   = false;
  vars->isTmpFile          = true;
  vars->writer     = NULL;
  vars->propertiesFilename = "";
  vars->nTracesOut = 0;
  vars->isFirstCall = true;

  int height = 800;
  int width  = 600;

  std::string text;
  stringstream out;
  int value;
  float valueFloat;
  int valueInt;
//------------------------------------------------------------

  std::string filename_image;
  param->getStringAtLine( "filename", &filename_image );
  vars->command.append( " -o " + filename_image );

  // Create random temporary file name
  // Make sure that file name is unique
  std::ostringstream pid;
  pid << getpid();  // Get process ID
  std::string tmpfile_template = string("tmp" + pid.str() + filename_image);
  // Replace letters '/', '\' and ',' in file name by underscores '_'
  for( int i = 0; i < (int)tmpfile_template.length(); i++ ) {
    char c = tmpfile_template[i];
    if( c == '/' || c == '\\' || c == '.' ) {
      //|| c == '\ ' ) {
      tmpfile_template[i] = '_';
    }
  }

  if( param->exists("temp_dir") ) {
    string dirname_temp;
    param->getString( "temp_dir", &dirname_temp );
#ifdef PLATFORM_WINDOWS
    tmpfile_template = string(dirname_temp + "\\" + tmpfile_template);
  }
  else {
    tmpfile_template = string("C:\\Windows\\Temp\\" + tmpfile_template);
#else
    tmpfile_template = string(dirname_temp + "/" + tmpfile_template);
  }
  else {
    tmpfile_template = string("/tmp/" + tmpfile_template);
#endif
  }

  vars->propertiesFilename = random_name( tmpfile_template + "", 10, ".prop" );

  if( edef->isDebug() ) {
    writer->line("Temporary file template: '%s'", tmpfile_template.c_str());
    writer->line("Temporary property file: '%s'", vars->propertiesFilename.c_str());
  }

  if( param->exists("seismic_filename") ) {
    vars->isTmpFile = false;
    param->getStringAtLine( "seismic_filename", &text );
    vars->command.append( " -f " + text );
    FILE* fin = fopen(text.c_str(),"r");
    if( fin == NULL ) writer->error("Input file not found: %s", text.c_str());
    fclose(fin);
  }
  else {
    vars->isTmpFile = true;
    vars->seismicFilename = random_name( tmpfile_template + "", 10, ".cseis" );
    if( !csFileUtils::createDoNotOverwrite( vars->seismicFilename ) ) {
      writer->error("Cannot open SeaSeis output file '%s'", vars->seismicFilename.c_str() );
    }
    vars->command.append( " -f " + vars->seismicFilename );
  }
  //------------------------------------------------
  if( param->exists( "height" ) ) {
    param->getInt( "height", &height );
  }
  out << height;
  vars->command.append( " -height " + out.str() );
  out.str("");
  if( param->exists( "width" ) ) {
    param->getInt( "width", &width );
  }
  out << width;
  vars->command.append( " -width " + out.str() );
  out.str("");
  //------------------------------------------------

  if( param->exists( "fontsize" ) ) {
    param->getInt( "fontsize", &value );
    out << value;
    vars->command.append( " -fontsize " + out.str() );
    out.str("");
  }
  if( param->exists( "window" ) ) {
    param->getFloat( "window", &valueFloat, 0 ); 
    char bufferFloat[80];
    sprintf(bufferFloat,"%.3f",valueFloat);
    // 'out' does not tolerate floats on Windows.. leads to program crash
    out << bufferFloat;
    vars->command.append( " -min_time " + out.str() );
    out.str("");
    param->getFloat( "window", &valueFloat, 1 );
    sprintf(bufferFloat,"%.3f",valueFloat);
    out << bufferFloat;
    vars->command.append( " -max_time " + out.str() );
    out.str("");
  }
  if( param->exists( "trace_spacing" ) ) {
    param->getString( "trace_spacing", &text );
    if( !text.compare("even") ) {
      vars->command.append( " -trace_spacing even " );
    }
    else if( !text.compare("auto") ) {
      vars->command.append( " -trace_spacing auto " );
    }
    else {
      writer->error("Unknown option: %s", text.c_str() );
    }
  }
  if( param->exists( "db" ) ) {
    param->getString( "db", &text );
    if( !text.compare("yes") ) {
      vars->command.append( " -db yes " );
    }
    else if( !text.compare("no") ) {
      vars->command.append( " -db no " );
    }
    else {
      writer->error("Unknown option: %s", text.c_str() );
    }
  }
  if( param->exists( "color_bar" ) ) {
    param->getString( "color_bar", &text, 0 );
    if( text.compare("bottom") && text.compare("right") ) {
      writer->error("Unknown option: %s", text.c_str() );
    }
    if( !text.compare("bottom") ) {
      writer->error("Option not supported yet: %s", text.c_str() );
    }
    text = "right";
    if( param->getNumValues("color_bar" ) > 1 ) {
      param->getInt( "color_bar", &valueInt, 1 );
      out << valueInt;
    }
    vars->command.append( " -color_bar " + text + out.str() );
    out.str("");
  }
  if( param->exists( "color_bar_ann" ) ) {
    param->getString( "color_bar_ann", &text, 0 );
    if( text.compare("simple") && text.compare("fancy") ) {
      writer->error("Unknown option: %s", text.c_str() );
    }
    if( param->getNumValues("color_bar_ann" ) > 1 ) {
      param->getInt( "color_bar_ann", &valueInt, 1 );
      out << valueInt;
    }
    vars->command.append( " -color_bar_ann " + text + out.str() );
    out.str("");
  }
  if( param->exists( "header" ) ) {
    param->getString( "header", &text );
    vars->command.append( " -header " + text );
  }
  if( param->exists( "properties" ) ) {
    param->getString( "properties", &text );
    if( !text.compare("remove") ) {
      vars->keepPropertyFile = false;
    }
    else if( !text.compare("keep") ) {
      vars->keepPropertyFile = true;
    }
    else {
      writer->error("Unknown option: %s", text.c_str() );
    }
  }

  //================================================================================
  // Read in properties and save properties file
  //
  /*
  FILE* test = fopen( vars->propertiesFilename.c_str(), "r" );
  while( test != NULL ) {
    writer->line("Info: Random property filename '%s' already exists. Change random name...", vars->propertiesFilename.c_str() );
    fclose( test );
    vars->propertiesFilename.append("x");
    test = fopen( vars->propertiesFilename.c_str(), "r" );
  }
  */
  FILE* foutProp = fopen( vars->propertiesFilename.c_str(), "w" );

  vars->command.append( " -p " + vars->propertiesFilename );

  if( param->exists( "title" ) ) {
    param->getString( "title", &text );
    fprintf(foutProp,"title:%s\n", text.c_str());
  }
  if( param->exists( "title_vert_axis" ) ) {
    param->getString( "title_vert_axis", &text );
    fprintf(foutProp,"titleVertAxis:%s\n", text.c_str());
  }
  if( param->exists("disp_scalar")  ) {
    param->getFloat( "disp_scalar", &valueFloat );
    fprintf(foutProp,"dispScalar:%e\n", valueFloat);
  }
  if( param->exists("min_value")  ) {
    param->getFloat( "min_value", &valueFloat );
    fprintf(foutProp,"minValue:%e\n", valueFloat );
  }
  if( param->exists("max_value")  ) {
    param->getFloat( "max_value", &valueFloat );
    fprintf(foutProp,"maxValue:%e\n", valueFloat );
  }
  if( param->exists("time_major_inc")  ) {
    if( !param->exists("time_minor_inc") ) writer->error("Both time_major_inc and time_minor_inc need to be specified, or none (automatic setting)."); 
    fprintf(foutProp,"isTimeLinesAuto:false\n");
    param->getFloat( "time_major_inc", &valueFloat );
    fprintf(foutProp,"timeLineMajorInc:%f\n", valueFloat );
  }
  if( param->exists("time_minor_inc")  ) {
    if( !param->exists("time_major_inc") ) writer->error("Both time_major_inc and time_minor_inc need to be specified, or none (automatic setting)"); 
    param->getFloat( "time_minor_inc", &valueFloat );
    fprintf(foutProp,"timeLineMinorInc:%f\n", valueFloat );
  }
  if( param->exists("time_max_decimals")  ) {
    param->getInt( "time_max_decimals", &value );
    out << value;
    vars->command.append( " -max_decimals " + out.str() );
    out.str("");
  }
  if( param->exists("hdr_date_format")  ) {
    param->getString( "hdr_date_format", &text );
    if( !text.compare("none") ) {
      // Nothing
    }
    else if( !text.compare("hours") ) {
      vars->command.append( " -hdr_date_type h" );
    }
    else if( !text.compare("days") ) {
      vars->command.append( " -hdr_date_type d" );
    }
    else if( !text.compare("day_month") ) {
      vars->command.append( " -hdr_date_type m" );
    }
    else if( !text.compare("day_month_year") ) {
      vars->command.append( " -hdr_date_type y" );
    }
    else {
      writer->error("Unknown option: ", text.c_str());
    }
  }

  if( param->exists("trace_clip")  ) {
    param->getFloat( "trace_clip", &valueFloat );
    if( valueFloat == 0.0 ) {
      valueFloat = 1.0;
      fprintf(foutProp,"doTraceClipping:false\n" );
    }
    else {
      fprintf(foutProp,"doTraceClipping:true\n" );
    }
    fprintf(foutProp,"traceClip:%f\n", valueFloat );
    
  }
  if( param->exists("polarity")  ) {
    param->getString( "polarity", &text );
    if( !text.compare("normal") ) {
      fprintf(foutProp,"polarity:1.0\n");
    }
    else if( !text.compare("reverse") ) {
      fprintf(foutProp,"polarity:-1.0\n");
    }
    else {
      writer->error("Unknown option: %s", text.c_str());
    }
  }

  if( param->exists("scale_type")  ) {
    param->getString( "scale_type", &text );
    //    fprintf(":\n", value);
    if( !text.compare("scalar") ) {
      fprintf(foutProp,"scaleType:%d\n", SCALE_TYPE_SCALAR);
    }
    else if( !text.compare("range") ) {
      fprintf(foutProp,"scaleType:%d\n", SCALE_TYPE_RANGE);
    }
    else if( !text.compare("full_trace") ) {
      fprintf(foutProp,"scaleType:%d\n", SCALE_TYPE_TRACE);
    }
    else {
      writer->error("Unknown option: %s", text.c_str());
    }
  }
  if( param->exists("wiggle")  ) {
    param->getString( "wiggle", &text );
    if( !text.compare("none") ) {
      fprintf(foutProp,"showWiggle:false\n");
    }
    else {
      fprintf(foutProp,"showWiggle:true\n");
      if( !text.compare("linear") ) {
        fprintf(foutProp,"wiggleType:%d\n", WIGGLE_TYPE_LINEAR);
      }
      else if( !text.compare("cubic") ) {
        fprintf(foutProp,"wiggleType:%d\n", WIGGLE_TYPE_CUBIC);
      }
      else {
        writer->error("Unknown option: %s", text.c_str());
      }
    }
  }
  if( param->exists("vi_type")  ) {
    param->getString( "vi_type", &text );
    if( !text.compare("none") ) {
      fprintf(foutProp,"isVIDisplay:false\n");
      fprintf(foutProp,"viType:%d\n", VA_TYPE_DISCRETE);
    }
    else {
      fprintf(foutProp,"isVIDisplay:true\n");
      if( !text.compare("discrete") ) {
        fprintf(foutProp,"viType:%d\n", VA_TYPE_DISCRETE);
      }
      else if( !text.compare("vertical") ) {
        fprintf(foutProp,"viType:%d\n", VA_TYPE_VERTICAL);
      }
      else if( !text.compare("spline") ) {
        fprintf(foutProp,"viType:%d\n", VA_TYPE_2DSPLINE);
      }
      else {
        writer->error("Unknown option: %s", text.c_str());
      }
    }
  }
  if( param->exists("show_zero_lines")  ) {
    param->getString( "show_zero_lines", &text );
    if( !text.compare("yes") ) {
      fprintf(foutProp,"showZeroLines:true\n");
    }
    else if( !text.compare("no") ) {
      fprintf(foutProp,"showZeroLines:false\n");
    }
    else {
      writer->error("Unknown option: %s", text.c_str());
    }
  }
  if( param->exists("show_time_lines")  ) {
    param->getString( "show_time_lines", &text );
    if( !text.compare("yes") ) {
      fprintf(foutProp,"showTimeLines:true\n");
    }
    else if( !text.compare("no") ) {
      fprintf(foutProp,"showTimeLines:false\n");
    }
    else {
      writer->error("Unknown option: %s", text.c_str());
    }
  }
  if( param->exists("pos_fill")  ) {
    param->getString( "pos_fill", &text );
    if( !text.compare("yes") ) {
      fprintf(foutProp,"isPosFill:true\n");
    }
    else if( !text.compare("no") ) {
      fprintf(foutProp,"isPosFill:false\n");
    }
    else {
      writer->error("Unknown option: %s", text.c_str());
    }
  }
  if( param->exists("neg_fill")  ) {
    param->getString( "neg_fill", &text );
    if( !text.compare("yes") ) {
      fprintf(foutProp,"isNegFill:true\n");
    }
    else if( !text.compare("no") ) {
      fprintf(foutProp,"isNegFill:false\n");
    }
    else {
      writer->error("Unknown option: %s", text.c_str());
    }
  }
  if( param->exists("log_scale")  ) {
    param->getString( "log_scale", &text );
    if( !text.compare("yes") ) {
      fprintf(foutProp,"isLogScale:true\n");
    }
    else if( !text.compare("no") ) {
      fprintf(foutProp,"isLogScale:false\n");
    }
    else {
      writer->error("Unknown option: %s", text.c_str());
    }
  }
  if( param->exists("var_color")  ) {
    param->getString( "var_color", &text );
    if( !text.compare("yes") ) {
      fprintf(foutProp,"isVariableColor:true\n");
    }
    else if( !text.compare("no") ) {
      fprintf(foutProp,"isVariableColor:false\n");
    }
    else {
      writer->error("Unknown option: %s", text.c_str());
    }
  }
  if( param->exists("plot_direction")  ) {
    param->getString( "plot_direction", &text );
    if( !text.compare("vertical") ) {
      fprintf(foutProp,"plotDirection:%d\n", PLOT_DIR_VERTICAL);
    }
    else if( !text.compare("horizontal") ) {
      fprintf(foutProp,"plotDirection:%d\n", PLOT_DIR_HORIZONTAL);
    }
    else {
      writer->error("Unknown option: %s", text.c_str());
    }
  }
  if( param->exists("vi_color_map")  ) {
    param->getString( "vi_color_map", &text );
    fprintf(foutProp,"viColorMap:%s\n", text.c_str());
  }
  if( param->exists("wiggle_color_map")  ) {
    param->getString( "wiggle_color_map", &text );
    fprintf(foutProp,"wiggleColorMap:%s\n", text.c_str());
  }
  if( param->exists("custom_vi_color_map")  ) {
    if( param->exists("vi_color_map")  ) {
      writer->warning("Specify either parameter vi_color_map OR custom_vi_color_map, not both. custom_vi_color_map will be ignored");
    }
    else {
      param->getString( "custom_vi_color_map", &text );
      fprintf(foutProp,"viColorMap:%s\n", text.c_str());
    }
  }
  if( param->exists("custom_wiggle_color_map")  ) {
    if( param->exists("wiggle_color_map")  ) {
      writer->warning("Specify either parameter wiggle_color_map OR custom_wiggle_color_map, not both. custom_wiggle_color_map will be ignored");
    }
    else {
      param->getString( "custom_wiggle_color_map", &text );
      fprintf(foutProp,"wiggleColorMap:%s\n", text.c_str());
    }
  }

  fclose(foutProp);

#if PLATFORM_WINDOWS
  vars->command.append( "\"" );  // Add double quotation marks around plotimage command line arguments
#endif

  if( edef->isDebug() ) fprintf(stderr,"Command to execute: '%s'\n", vars->command.c_str() );
  writer->line("\nCommand to execute:\n    %s\n", vars->command.c_str() );
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_image_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  //  csExecPhaseDef* edef = env->execPhaseDef;
  csSuperHeader const* shdr = env->superHeader;
  csTraceHeaderDef const* hdef = env->headerDef;
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );

  csTrace* trace = traceGather->trace(0);


  if( vars->isFirstCall ) {
    vars->isFirstCall = false;
    try {
      vars->writer = new csSeismicWriter( vars->seismicFilename, 20 );
    }
    catch( csException& exc ) {
      writer->error("Error occurred when opening SeaSeis file. System message:\n%s", exc.getMessage() );
    }
  }


  if( vars->isTmpFile ) {
    float* samples = trace->getTraceSamples();
    char const* hdrValueBlock = trace->getTraceHeader()->getTraceHeaderValueBlock();

    try {
      bool success = true;
      if( vars->nTracesOut == 0 ) {
        success = vars->writer->writeFileHeader( shdr, hdef );
      }
      if( success ) success = vars->writer->writeTrace( samples, hdrValueBlock );
      if( !success ) {
        writer->error("Unknown error occurred when writing to SeaSeis file.\n" );
      }
    }
    catch( csException& exc ) {
      writer->error("Error occurred when writing to SeaSeis file. System message:\n%s", exc.getMessage() );
    }
    vars->nTracesOut += 1;
  }


  return;
}
//********************************************************************************
// Parameter definition
//
//
//********************************************************************************
void params_mod_image_( csParamDef* pdef ) {
  pdef->setModule( "IMAGE", "Create image file from seismic display.." );

  pdef->addParam( "filename", "Output image file name", NUM_VALUES_FIXED,
                  "Supported extensions/formats: jpg/jpeg, png, bmp, gif (or upper case letters)" );
  pdef->addValue( "", VALTYPE_STRING, "Output image file name" );

  pdef->addParam( "seismic_filename", "Input seismic file name", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Input seismic file name" );

  pdef->addParam( "temp_dir", "Temporary file directory for storing temporary files", NUM_VALUES_FIXED,
                  "Default is /tmp for Linux/UNIX system, and C:/Windows/Temp for Windows system");
  pdef->addValue( "", VALTYPE_STRING, "Directory name" );

  pdef->addParam( "properties", "Action to perform on temporary 'properties file' that is auto-generated to communicate display settings between SeaSeis and PlotImage", NUM_VALUES_FIXED );
  pdef->addValue( "remove", VALTYPE_OPTION );
  pdef->addOption( "remove", "Remove temporary 'properties file' after module finished");
  pdef->addOption( "keep", "Keep temporary 'properties file', do not remove");

  pdef->addParam( "width", "Image width in pixels", NUM_VALUES_FIXED );
  pdef->addValue( "600", VALTYPE_NUMBER, "Image width in pixels" );

  pdef->addParam( "height", "Image height in pixels", NUM_VALUES_FIXED );
  pdef->addValue( "800", VALTYPE_NUMBER, "Image height in pixels" );

  pdef->addParam( "trace_spacing", "Trace spacing", NUM_VALUES_FIXED );
  pdef->addValue( "even", VALTYPE_OPTION );
  pdef->addOption( "even", "Ensure even trace spacing (traces/pixel = integer value)",
                   "Generally, this will lead to an adjustment in the image width" );
  pdef->addOption( "auto", "Automatic trace spacing dictated by image width" );

  pdef->addParam( "header", "Trace header for trace annotation", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Trace header name" );

  pdef->addParam( "hdr_date_format", "Date format", NUM_VALUES_FIXED, "Only applicable for integer type trace headers containing seconds since 1970-1-1" );
  pdef->addValue( "none", VALTYPE_OPTION );
  pdef->addOption( "none", "Do not convert trace header value to date" );
  pdef->addOption( "hours", "Display time in number of hours" );
  pdef->addOption( "days", "Display time in number of days" );
  pdef->addOption( "day_month", "Display time as day/month" );
  pdef->addOption( "day_month_year", "Display time as day/month/year" );

  pdef->addParam( "title", "Title text", NUM_VALUES_FIXED,
                  "To print trace header value from first data trace in title text, use the following syntax: %headerName%" );
  pdef->addValue( "", VALTYPE_STRING, "Title text" );

  pdef->addParam( "fontsize", "Font size in pixels", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "Font size in pixels" );

  pdef->addParam( "window", "Window to display", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "Minimum time/frequency to display [ms] or [Hz]" );
  pdef->addValue( "", VALTYPE_NUMBER, "Maximum time/frequency to display [ms] or [Hz]" );

  pdef->addParam( "db", "Convert to dB?", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "no", "Do not convert to dB");
  pdef->addOption( "yes", "Convert (absolute) amplitudes to dB" );

  pdef->addParam( "color_bar", "Select this option to show the colour bar (Variable intensity/VI colour map only)", NUM_VALUES_FIXED,
                  "Min/max values are defined by given min/max range, parameters 'min_value', and 'max_value'" );
  pdef->addValue( "right", VALTYPE_OPTION );
  pdef->addOption( "bottom", "Display colour bar at bottom of image");
  pdef->addOption( "right", "Display colour bar on the right hand side of image" );
  pdef->addValue( "30", VALTYPE_NUMBER, "Size of colour bar" );

  pdef->addParam( "color_bar_ann", "Colour bar annotation", NUM_VALUES_FIXED );
  pdef->addValue( "simple", VALTYPE_OPTION );
  pdef->addOption( "simple", "Simple annotation");
  pdef->addOption( "fancy", "Fancy annotation" );
  pdef->addValue( "0", VALTYPE_NUMBER, "Size of annotation area" );

  //---------------------------------------------------------------------
  //
  //
  pdef->addParam( "disp_scalar", "Display scalar. Data is scaled by this value before plotting", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "" );

  //  pdef->addParam( "full_trace_scalar", "Display scalar for 'full trace' scaling", NUM_VALUES_FIXED );
  //  pdef->addValue( "", VALTYPE_NUMBER, "" );

  pdef->addParam( "scale_type", "Type of scaling", NUM_VALUES_FIXED );
  pdef->addValue( "scalar", VALTYPE_OPTION );
  pdef->addOption( "scalar", "Use simple scalar to scale data before plotting" );
  pdef->addOption( "range", "Apply range (min/max) for plotting" );
  pdef->addOption( "full_trace", "Apply full trace equalisation scalar before plotting" );

  pdef->addParam( "min_value", "Minimum value for 'range' scaling", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "" );

  pdef->addParam( "max_value", "Maximum value for 'range' scaling", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "" );

  pdef->addParam( "wiggle", "Wiggle trace display", NUM_VALUES_FIXED );
  pdef->addValue( "linear", VALTYPE_OPTION );
  pdef->addOption( "linear", "Apply linear interpolation between samples" );
  pdef->addOption( "cubic", "Apply cubic interpolation between samples" );
  pdef->addOption( "none", "Do not plot wiggle trace" );

  pdef->addParam( "vi_type", "Variable intensity plot", NUM_VALUES_FIXED );
  pdef->addValue( "none", VALTYPE_OPTION );
  pdef->addOption( "spline", "2D spline interpolation" );
  pdef->addOption( "vertical", "Vertical interpolation" );
  pdef->addOption( "discrete", "Discrete samples" );
  pdef->addOption( "none", "No variable intensity plot" );

  pdef->addParam( "pos_fill", "Fill positive wiggle?", NUM_VALUES_FIXED );
  pdef->addValue( "yes", VALTYPE_OPTION );
  pdef->addOption( "yes", "Fill positive wiggle" );
  pdef->addOption( "no", "Do not fill positive wiggle" );

  pdef->addParam( "neg_fill", "Fill negative wiggle?", NUM_VALUES_FIXED );
  pdef->addValue( "yes", VALTYPE_OPTION );
  pdef->addOption( "yes", "Fill negative wiggle" );
  pdef->addOption( "no", "Do not fill negative wiggle" );

  pdef->addParam( "var_color", "Use variable colour for wiggle fill?", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "yes", "Use variable colour for wiggle fill" );
  pdef->addOption( "no", "Use constant colour for wiggle fill" );

  pdef->addParam( "show_zero_lines", "Show zero lines (amplitude = 0) for each trace?", NUM_VALUES_FIXED );
  pdef->addValue( "yes", VALTYPE_OPTION );
  pdef->addOption( "yes", "Show zero line" );
  pdef->addOption( "no", "Do not show zero line" );

  pdef->addParam( "show_time_lines", "Show time lines?", NUM_VALUES_FIXED );
  pdef->addValue( "yes", VALTYPE_OPTION );
  pdef->addOption( "yes", "Show time lines" );
  pdef->addOption( "no", "Do not show time lines" );

  pdef->addParam( "time_major_inc", "Time axis annotation: Major line increment [ms/Hz]", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "" );

  pdef->addParam( "time_minor_inc", "Time axis annotation: Minor line increment [ms/Hz]", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "" );

  pdef->addParam( "time_max_decimals", "Maximum number of decimal places to display for time/frequency axis annotation", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "" );

  pdef->addParam( "plot_direction", "Plot direction: Vertical or horizontal", NUM_VALUES_FIXED );
  pdef->addValue( "vertical", VALTYPE_OPTION );
  pdef->addOption( "vertical", "Plot data vertically (normal)." );
  pdef->addOption( "horizontal", "Plot data horizontally (90deg transposed)." );

  pdef->addParam( "wiggle_color_map", "Colour map for variable colour wiggle fill", NUM_VALUES_FIXED );
  pdef->addValue( "blue_white_red", VALTYPE_OPTION );
  pdef->addOption( "gray_w2b", "Gray scale, white-black" );
  pdef->addOption( "gray_b2w", "Gray scale, black-white" );
  pdef->addOption( "gray_bwb", "Gray scale, black-white-black" );
  pdef->addOption( "gray_wbw", "Gray scale, white-black-white" );
  pdef->addOption( "blue_white_red", "Blue white red" );
  pdef->addOption( "blue_white_red2", "Blue white red (2)" );
  pdef->addOption( "black_white_orange", "Black white orange" );
  pdef->addOption( "black_white_red", "Black white red" );
  pdef->addOption( "rainbow", "Rainbow colour scale" );
  pdef->addOption( "rainbow_black", "Rainbow colour scale with black" );
  pdef->addOption( "rainbow_2", "Rainbow(2) colour scale" );
  pdef->addOption( "rainbow_mirror", "Rainbow colour scale, mirrored around zero" );
  pdef->addOption( "cold_warm", "Cold to warm colour" );
  pdef->addOption( "brown", "Brown" );
  pdef->addOption( "default", "Default color map (good for interpretation)" );

  pdef->addParam( "vi_color_map", "Colour map for variable intensity plot", NUM_VALUES_FIXED );
  pdef->addValue( "gray_w2b", VALTYPE_OPTION );
  pdef->addOption( "gray_w2b", "Gray scale, white-black" );
  pdef->addOption( "gray_b2w", "Gray scale, black-white" );
  pdef->addOption( "gray_bwb", "Gray scale, black-white-black" );
  pdef->addOption( "gray_wbw", "Gray scale, white-black-white" );
  pdef->addOption( "blue_white_red", "Blue white red" );
  pdef->addOption( "blue_white_red2", "Blue white red (2)" );
  pdef->addOption( "black_white_orange", "Black white orange" );
  pdef->addOption( "black_white_red", "Black white red" );
  pdef->addOption( "rainbow", "Rainbow colour scale" );
  pdef->addOption( "rainbow_black", "Rainbow colour scale with black" );
  pdef->addOption( "rainbow_2", "Rainbow(2) colour scale" );
  pdef->addOption( "rainbow_mirror", "Rainbow colour scale, mirrored around zero" );
  pdef->addOption( "cold_warm", "Cold to warm colour" );
  pdef->addOption( "brown", "Brown" );
  pdef->addOption( "default", "Default color map (good for interpretation)" );

  pdef->addParam( "custom_wiggle_color_map", "Custom colour map for variable colour wiggle fill", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Exact name of custom color map." );

  pdef->addParam( "custom_vi_color_map", "Custom colour map for variable intensity plot", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Exact name of custom color map." );

  pdef->addParam( "trace_clip", "Number of traces where seismic wiggle is clipped", NUM_VALUES_FIXED,
                  "Set to 0 to avoid any clipping" );
  pdef->addValue( "", VALTYPE_NUMBER, "" );

  pdef->addParam( "polarity", "Plot polarity convention", NUM_VALUES_FIXED );
  pdef->addValue( "normal", VALTYPE_OPTION );
  pdef->addOption( "normal", "Plot normal polarity (negative value are on left hand side)" );
  pdef->addOption( "reverse", "Plot reverse polarity (negative values are on right hand side)" );

  pdef->addParam( "log_scale", "Plot is log scale?", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "no", "Normal plot" );
  pdef->addOption( "yes", "Plot as log scale (time/frequency axis)" );

  pdef->addParam( "title_vert_axis", "Title text for vertical axis", NUM_VALUES_FIXED, "If not specified, default text will be used" );
  pdef->addValue( "", VALTYPE_STRING );

  /*  pdef->addParam( "override", "Override domain specified in superheader.", NUM_VALUES_VARIABLE );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "no", "Acknowledge domain found in super header" );
  pdef->addOption( "fx", "Override domain found in super header. Set to FK." );
  pdef->addOption( "fx", "Override domain found in super header. Set to FX." );
  pdef->addOption( "xt", "Override domain found in super header. Set to XT (time)." );
  pdef->addOption( "xd", "Override domain found in super header. Set to XD (depth)." );
  pdef->addValue( "amp_phase", VALTYPE_OPTION, "Override data type (in case of FX/FK)" );
  pdef->addOption( "amp_phase", "Set FFT data type to amplitude/phase spectrum" );
  pdef->addOption( "amp", "Set FFT data type to amplitude spectrum" );
  pdef->addOption( "psd", "Set FFT data type to psd spectrum" );
  pdef->addOption( "complex", "Set FFT data type to complex spectrum (only applies to FK)" );
  pdef->addOption( "real_imag", "Set FFT data type to real/imaginary spectrum" );
  */
}


//************************************************************************************************
// Start exec phase
//
//*************************************************************************************************
bool start_exec_mod_image_( csExecPhaseEnv* env, csLogWriter* writer ) {
//  mod_image::VariableStruct* vars = reinterpret_cast<mod_image::VariableStruct*>( env->execPhaseDef->variables() );
//  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;
//  csTraceHeaderDef const* hdef = env->headerDef;
  return true;
}

//************************************************************************************************
// Cleanup phase
//
//*************************************************************************************************
void cleanup_mod_image_( csExecPhaseEnv* env, csLogWriter* writer ) {
  mod_image::VariableStruct* vars = reinterpret_cast<mod_image::VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;
  if( vars->writer != NULL ) {
    delete vars->writer;
    vars->writer = NULL;
  }

  int returnFlag = 0;

  if( edef->isDebug() ) fprintf(stderr,"Executing command: '%s'\n", vars->command.c_str() );
  writer->line("\n...executing command:\n    %s\n", vars->command.c_str() );
  returnFlag = system( vars->command.c_str() );
  if( returnFlag != 0 ) {
    writer->error("$IMAGE: Error occurred when executing system command: '%s'", vars->command.c_str());
  }
#if PLATFORM_WINDOWS
  std::string command_rm = "del ";
#else
  std::string command_rm = "rm -f ";
#endif

  if( !vars->keepPropertyFile ) {
    std::string command1 = command_rm + vars->propertiesFilename;
    returnFlag = system( command1.c_str() );
    if( returnFlag != 0 ) {
      writer->error("$IMAGE: Error occurred when executing system command: '%s'", command1.c_str());
    }
  }
  if( vars->isTmpFile ) {
    std::string command2 = command_rm + vars->seismicFilename;
    returnFlag = system( command2.c_str() );
    if( returnFlag != 0 ) {
      writer->error("$IMAGE: Error occurred when executing system command: '%s'", command2.c_str());
    }
  }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_image_( csParamDef* pdef ) {
  params_mod_image_( pdef );
}
extern "C" void _init_mod_image_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_image_( param, env, writer );
}
extern "C" bool _start_exec_mod_image_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_image_( env, writer );
}
extern "C" void _exec_mod_image_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_image_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_image_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_image_( env, writer );
}


std::string random_name( std::string pre, int size, std::string post ) {
  std::stringstream str;
  srand( time(NULL) );
  int counter = 0;
  while( counter < size ) {
    char c = rand() % 26 + 'a';
    str << c;
    counter += 1;
  }

  // Check if temporary file already exists
  string text = pre + str.str() + post;
  
  FILE* test = fopen( text.c_str(), "r" );
  while( test != NULL ) {
    fclose( test );
    str << 'x';
    text = pre + str.str() + post;    
    test = fopen( text.c_str(), "r" );
  }
  
  return text;
}

