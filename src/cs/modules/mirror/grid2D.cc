/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <string>
#include "csPoint3D.h"

using namespace std;
using namespace cseis_geolib;

/**
 * Interpolate sparse/irregular 2D map to regular 2D grid
 * Use nearest neighbour falling into each bin.
 * If more than one point falls into the bin, use weighted sum.
 * Grid bins where no input points fall into can be linearly interpolated.
 * No extrapolation is done outside the area covered my input points.
 *
 * @param numPointsIn  Number of points in input arrays
 * @param pin          Input point (Z coordinate is the value that shall be interpolated)
 * @param p1           Location of first grid bin (XY coordinates only)
 * @param xinc         Bin size in X direction
 * @param yinc         Bin size in Y direction
 * @param nxOut        Number of output bins in X dimension
 * @param nyOut        Number of output bins in Y dimension
 * @param zout         Output array, interpolated values on regular grid
 * @param no_value     Value that indicates an invalid value, or 'no value'
 * @param fillHoles    true if bins that do not contain an input point shall be interpolated
 */
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
                                 int smoothSize )
{
  double* weights = new double[nxOut*nyOut];

  // Initialise output grid. Set Z value = no_value to indicate that this bin has not been assigned yet
  for( int ix = 0; ix < nxOut; ix++ ) {
    for( int iy = 0; iy < nyOut; iy++ ) {
      int binIndex = ix*nyOut + iy;
      zout[binIndex]   = no_value;
      weights[binIndex] = 0.0;
    }
  }

  // Scan through all existing points, compute the bin into which each point falls
  for( int i = 0; i < numPointsIn; i++ ) {
    double x = pin[i].x;
    double y = pin[i].y;
    double z = pin[i].z;
    int binx = (int)( (x-p1.x)/xinc + 0.5 );
    int biny = (int)( (y-p1.y)/yinc + 0.5 );
    if( binx < 0 || binx >= nxOut || biny < 0 || biny >= nyOut ) continue;
    int bin = binx*nyOut + biny;
    double xout = p1.x + binx*xinc;
    double yout = p1.y + biny*yinc;
    double distance = sqrt( pow(x-xout,2) + pow(y-yout,2) );
    if( distance < 1e-6 ) {
      distance = 1e-6;
    }
    double w = 1.0/distance;
    weights[bin] += w;
    if( zout[bin] == no_value ) zout[bin] = 0.0;
    zout[bin]   += w*z;
  }

  // Compute mean value for existing bins
  for( int ix = 0; ix < nxOut; ix++ ) {
//    double xout = ix*xinc + p1.x;
    for( int iy = 0; iy < nyOut; iy++ ) {
//      double yout = iy*yinc + p1.y;
      int bin = ix*nyOut + iy;
      if( zout[bin] != no_value ) {
        zout[bin] /= weights[bin];
        weights[bin] = 0.0;
      }
    }
  }

  // Interpolate empty bins
  if( fillHoles ) {
    int maxIterations = 3;
    int iteration = 0;
    bool assignedPoint = true;
    while( iteration < maxIterations && assignedPoint ) {
      assignedPoint = false;
      iteration += 1;

      // Interpolate empty bins along Y axis
      for( int ix = 1; ix < nxOut-1; ix++ ) {
        //      double xout = ix*xinc + p1.x;
        int iy_start = 0;
        int iy_end = nyOut-1;
        while( iy_start < nyOut && zout[ix*nyOut+iy_start] == no_value ) iy_start += 1;
        while( iy_end > iy_start && zout[ix*nyOut+iy_end] == no_value ) iy_end -= 1;
        
        for( int iy = iy_start; iy < iy_end; iy++ ) {
          int bin = ix*nyOut + iy;
          if( zout[bin] == no_value ) {
            double z_y1   = zout[ix*nyOut + iy-1];
            double z      = z_y1;
            double weight = 1.0;
            int iy2 = iy+1;
            while( zout[ix*nyOut + iy2] == no_value ) iy2 += 1;
            //        iy2 = iy+1;
            double z_y2 = zout[ix*nyOut + iy2];
            double w = 1/(double)(iy2-iy);
            if( zout[ix*nyOut + iy2] != no_value ) {
              z += z_y2 * w;
              weight += w;
            }
            double z_x1 = zout[(ix-1)*nyOut + iy];
            if( z_x1 != no_value ) {
              z += z_x1;
              weight += 1.0;
            }
            double z_x2 = zout[(ix+1)*nyOut + iy];
            if( z_x2 != no_value ) {
              z += z_x2;
              weight += 1.0;
            }
            zout[bin] = z/weight;
            weights[bin] = no_value;   // Signal that this bin has been re-interpolated
            assignedPoint = true;
          }
        }
      }
      
      // Interpolate empty bins along X axis
      //      for( int iy = 1; iy < nyOut-1; iy++ ) {
      for( int iy = 0; iy < nyOut; iy++ ) {
        //      double yout = iy*yinc + p1.y;
        int ix_start = 0;
        int ix_end = nxOut-1;
        while( ix_start < nxOut && ( zout[ix_start*nyOut+iy] == no_value || weights[ix_start*nyOut+iy] == no_value ) ) ix_start += 1;
        while( ix_end > ix_start && ( zout[ix_end*nyOut+iy] == no_value || weights[ix_end*nyOut+iy] == no_value ) ) ix_end -= 1;
        
        for( int ix = ix_start; ix < ix_end; ix++ ) {
          int bin = ix*nyOut + iy;
          if( zout[bin] == no_value || weights[bin] == no_value ) {
            double z_x1   = zout[(ix-1)*nyOut + iy];
            double z      = z_x1;
            double weight = 1.0;
            int ix2 = ix+1;
            while( zout[ix2*nyOut + iy] == no_value || weights[ix2*nyOut + iy] == no_value ) ix2 += 1;
            double z_x2 = zout[ix2*nyOut + iy];
            double w = 1/(double)(ix2-ix);
            z += z_x2 * w;
            weight += w;
            double z_y1 = zout[ix*nyOut + iy-1];
            double z_y2 = zout[ix*nyOut + iy+1];
            
            if( z_y1 != no_value && weights[ix*nyOut + iy-1] != no_value ) {
              z += z_y1;
              weight += 1.0;
            }
            if( z_y2 != no_value && weights[ix*nyOut + iy+1] != no_value ) {
              z += z_y2;
              weight += 1.0;
            }
            if( weights[bin] == no_value ) {
              zout[bin] = 0.5 * ( zout[bin] + z/weight );
            }
            else {
              zout[bin] = z/weight;
            }
            assignedPoint = true;
          }
        }
      }
      
      for( int ix = 0; ix < nxOut; ix++ ) {
        for( int iy = 0; iy < nyOut; iy++ ) {
          int bin = ix*nyOut + iy;
          if( zout[bin] != no_value ) {
            weights[bin] = 0.0;
          }
        }
      }
      
    } // END while iteration
  }


  if( smoothGrid ) {
    for( int ix = 0; ix < nxOut; ix++ ) {
      for( int iy = 0; iy < nyOut; iy++ ) {
        int bin = ix*nyOut + iy;

        int counter = 0;
        double sum = 0.0;
        for( int ix2 = ix-smoothSize; ix2 <= ix+smoothSize; ix2++ ) {
          if( ix2 < 0 || ix2 >= nxOut ) continue;
          for( int iy2 = iy-smoothSize; iy2 <= iy+smoothSize; iy2++ ) {
            if( iy2 < 0 || iy2 >= nyOut ) continue;
            int bin2 = ix2*nyOut + iy2;
            if( zout[bin2] != no_value ) {
              sum = sum + zout[bin2];
              counter += 1;
            }
          }
        }
        if( counter > 0 ) {
          weights[bin] = sum/(double)counter;  // Temporary holding location
          //          weights[bin] = 0.0;
        }
        else {
          weights[bin] = no_value;
        }
      }
    }
    for( int ix = 0; ix < nxOut; ix++ ) {
      for( int iy = 0; iy < nyOut; iy++ ) {
        int bin = ix*nyOut + iy;
        zout[bin] = weights[bin];
      }
    }
  }


  delete [] weights;
}

