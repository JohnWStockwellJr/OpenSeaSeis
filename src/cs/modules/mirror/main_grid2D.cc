/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <string>
#include <cstring>
#include "geolib_string_utils.h"
#include "csPoint3D.h"
#include "csVector.h"
#include "geolib_string_utils.h"

using namespace std;
using namespace cseis_geolib;

void interpolate2D_regular_grid( int numPointsIn,
                                 csPoint3D const* pin,
                                 csPoint3D p1,
                                 double xinc,
                                 double yinc,
                                 int nxOut,
                                 int nyOut,
                                 double* zout,
                                 double no_value,
                                 bool fillHoles,
                                 bool smoothGrid,
                                 int smoothSize );

static const double NO_VALUE = std::numeric_limits<double>::max();

int main(int argc, char **argv) {

  double xinc = 0.0;
  double yinc = 0.0;
  bool set_xinc = false;
  bool set_yinc = false;
  bool fillHoles = true;
  bool smoothGrid = false;
  int smoothSize = 0;
  bool isXOnly = false;

  if( argc == 1 ) {
    fprintf(stderr," Type '%s -h' for help.\n", argv[0]);
    return(-1) ;
  }
  int iArg = 1;   // Start with 1 to skip name of program executable
  while( iArg < argc ) {
    if( argv[iArg][0] == '-' ) {
      if( strlen(argv[iArg]) < 2 ) {
        fprintf(stderr," Syntax error in command line: Unknown option: '%s'\n", argv[iArg]);
        fprintf(stderr," Type '%s -h' for help.\n", argv[0]);
        return(-1);
      }
      char option = argv[iArg][1];
      if ( option == 'h' ) {
        fprintf(stderr," 2D grid regularization.\n" );
        fprintf(stderr," Usage:  %s -x <bin size X> -y <bin size Y> [-n] [-h] < <filein> > <fileout>\n", argv[0]);
        fprintf(stderr," -x <bin size X> :  Bin size in X direction\n");
        fprintf(stderr," -y <bin size Y> :  Bin size in Y direction\n");
        fprintf(stderr," -n              :  No interpolation of empty bins\n");
        fprintf(stderr,"                 :  Grid bins that do not contain any input data points will not be interpolated\n");
        fprintf(stderr," -s <num bins>   :  Smooth grid after interpolation, use n bins in each direction\n");
        fprintf(stderr," -h              :  Print this page\n");
        fprintf(stderr," -x_only         :  Input file does not contain a Y coordinate, only X & Z\n");
        fprintf(stderr," <filein>        :  Input ASCII file. Format: X Y Z\n");
        fprintf(stderr," <fileout>       :  Output ASCII file. Format: X Y Z\n");
        return(-1);
      }
      else if ( option == 'n' ) {
        fillHoles = false;
        ++iArg;
      }
      else if ( option == 'x' ) {
        if( !strcmp(argv[iArg],"-x_only") ) {
          iArg++;
          isXOnly = true;
          yinc = 1.0;
          set_yinc = true;
        }
        else {
          iArg++;
          if( iArg == argc ) {
            fprintf(stderr,"Missing argument for option -%c\n", option);
            return(-1);
          }
          xinc = atof(argv[iArg]);
          set_xinc = true;
          ++iArg;
        }
      }
      else if ( option == 'y' ) {
        iArg++;
        if( iArg == argc ) {
          fprintf(stderr,"Missing argument for option -%c\n", option);
          return(-1);
        }
        yinc = atof(argv[iArg]);
        set_yinc = true;
        ++iArg;
      }
      else if ( option == 's' ) {
        iArg++;
        if( iArg == argc ) {
          fprintf(stderr,"Missing argument for option -%c\n", option);
          return(-1);
        }
        smoothSize = atoi(argv[iArg]);
        smoothGrid = true;
        ++iArg;
      }
      else {
        fprintf(stderr," Syntax error in command line: Unknown option: '%s'\n", argv[iArg]);
        fprintf(stderr," Type '%s -h' for help.\n", argv[0]);
        return(-1) ;
      }
    }
    else {
      fprintf(stderr," Syntax error in command line:\n Expected command line option (starting with '-'), found '%s'\n", argv[iArg]);
      fprintf(stderr," Type '%s -h' for help.\n", argv[0]);
      return(-1) ;
    }
  }

  if( !set_xinc || !set_yinc ) {
    fprintf(stderr," Grid bin size not set in input argument list.\n");
    fprintf(stderr," Type '%s -h' for help.\n", argv[0]);
    return(-1);
  }
  else if( xinc <= 0.0 || yinc <= 0.0 ) {
    fprintf(stderr," Incorrect grid bin size: %f %f\n", xinc, yinc);
    fprintf(stderr," Type '%s -h' for help.\n", argv[0]);
    return(-1);
  }

  fprintf(stderr,"Generate regular grid with bin size:   %f x %f\n", xinc, yinc);

  //
  //--------------------------------------------------------------------------------
  //

  char buffer[1024];
//  FILE* file = fopen("file.in","r");
  csVector<csPoint3D> pointList(1000);
  csVector<std::string> tokenList;

  double xmin = std::numeric_limits<double>::max();
  double ymin = std::numeric_limits<double>::max();
  double xmax = -(std::numeric_limits<double>::max()-1);
  double ymax = -(std::numeric_limits<double>::max()-1);

  int indexZ = 2;
  if( isXOnly ) indexZ = 1;

  while( fgets( buffer, 1024, stdin ) != NULL ) {
    tokenList.clear();
    tokenize( buffer, tokenList );
    double x = atof(tokenList.at(0).c_str());
    double y = 1.0;
    if( !isXOnly ) {
      y = atof(tokenList.at(1).c_str());
    }
    double z = atof(tokenList.at(indexZ).c_str());
    if( x < xmin ) xmin = x;
    if( x > xmax ) xmax = x;
    if( y < ymin ) ymin = y;
    if( y > ymax ) ymax = y;
    pointList.insertEnd( csPoint3D(x,y,z) );
  }

  int numPoints = pointList.size();

  csPoint3D p1(xmin,ymin,0);
  csPoint3D p2(xmax,ymax,0);

  int nx = (int)( (p2.x-p1.x)/xinc ) + 1;
  int ny = (int)( (p2.y-p1.y)/yinc ) + 1;
  p2.x = p1.x + nx*xinc;
  p2.y = p1.y + ny*yinc;

  double* zout = new double[nx*ny];

  interpolate2D_regular_grid( numPoints, pointList.toArray(),
                              p1, xinc, yinc, nx, ny, zout, NO_VALUE, fillHoles,
                              smoothGrid, smoothSize );

  for( int ix = 0; ix < nx; ix++ ) {
    double xout = ix*xinc + p1.x;
    for( int iy = 0; iy < ny; iy++ ) {
      double yout = iy*yinc + p1.y;
      int bin = ix*ny + iy;
      if( zout[bin] != NO_VALUE ) {
        if( !isXOnly ) fprintf(stdout,"%f %f %f\n", xout, yout, zout[bin]);
        else fprintf(stdout,"%f %f\n", xout, zout[bin]);
      }
    }
  }

}

