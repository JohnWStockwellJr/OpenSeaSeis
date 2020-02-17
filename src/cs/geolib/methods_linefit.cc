/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

/*-------------------------------------------------------------------
 * Fit 3D line to data points
 *
 * Computes 3D unit direction vector that best fits the data points in the given time window.
 * The solution is based upon the minimimum distance of data points to the direction vector.
 *
 * Note that the computed direction vector may be positive or negative.
 * The direction vector, when fitted to the data points, may not go through the origin (0,0,0), but may be offset to the origin.
 *
 * Suggested improvements:
 * - Compute perpendicular distance vector from origin
 * - Compute polarity and scale output vector accordingly.
 *   Requires additional time window
 * - Allocate matrices outside of this module...
 *
 *-------------------------------------------------------------------*/

#include <cmath>
#include <cstdio>
#include "geolib_mem.h"
#include "geolib_math.h"
#include "geolib_methods.h"

// These should be declared in header file!
#define LINEFIT_3D    11
#define LINEFIT_LS    12
#define LINEFIT_SVD   13

namespace cseis_geolib {

  bool linefit_3d( float* xSamples,
                  float* ySamples,
                  float* zSamples,
                  int firstSample,
                  int lastSample,
                  int force_origin,
                  float* vec_out )
  {
    const int NCOLS = 3;
    const int ID_X = 0;
    const int ID_Y = 1;
    const int ID_Z = 2;

    float norm;
    float** mat_a;
    float** mat_u;
    float** mat_v;
    float*  vec_w;
    float*  vec_work;
    float x_mean, y_mean, z_mean;
    float max_value;
    int max_col_index;
    int nSamples;
    int isamp, icol;
    bool ret = false;

    //------------------------------------------------------------

    nSamples = lastSample - firstSample + 1;

    //  for( isamp = firstSample; isamp <= lastSample; isamp++ ) {
    //   fprintf(stdout,"%f  %f  %d SAMPXY\n", xSamples[isamp], ySamples[isamp], isamp);
    //   fprintf(stdout,"%f  %f  %d SAMPXZ\n", xSamples[isamp], zSamples[isamp], isamp);
    //   fprintf(stdout,"%f  %f  %d SAMPYZ\n", ySamples[isamp], zSamples[isamp], isamp);
    //  }

    mat_a = allocate_matrix( nSamples, NCOLS );
    mat_u = allocate_matrix( nSamples, NCOLS );
    mat_v = allocate_matrix( NCOLS, NCOLS );
    vec_w    = allocate_vector( nSamples );
    vec_work = allocate_vector( nSamples );

    if( !mat_a || !mat_u || !mat_v || !vec_w || !vec_work) {
      return false;
    }

    //-------------------------------
    // Compute average xyz values, fill inversion matrix with reduced data point values
    // x - x_mean, y - y_mean, z - z_mean

    x_mean = y_mean = z_mean = 0.0;
    // Unless line is forced through origin, compute centroid point by averaging x,y,z:
    if( !force_origin ) {
      for( isamp = firstSample; isamp <= lastSample; isamp++ ) {
        x_mean += xSamples[isamp];
        y_mean += ySamples[isamp];
        z_mean += zSamples[isamp];
      }
      x_mean /= (float)nSamples;
      y_mean /= (float)nSamples;
      z_mean /= (float)nSamples;
    }
    for( isamp = 0; isamp < nSamples; isamp++ ) {
      mat_a[isamp][ID_X] = xSamples[firstSample+isamp] - x_mean;
      mat_a[isamp][ID_Y] = ySamples[firstSample+isamp] - y_mean;
      mat_a[isamp][ID_Z] = zSamples[firstSample+isamp] - z_mean;
      for( icol = 0; icol < NCOLS; icol++ ) {
        mat_u[isamp][icol] = mat_a[isamp][icol];
      }
    }

    //------------------------------
    // Decompose 'inversion' matrix
    //
    ret = svd_decomposition( mat_u, nSamples, NCOLS, vec_w, mat_v, vec_work );
    if( !ret ) {
      return false;
    }

    // Find maximum singular value, corresponding singular vector is the unit direction vector
    max_value = -1000.0;
    max_col_index = 0;
    for( icol = 0; icol < NCOLS; icol++ ) {
      if( vec_w[icol] > max_value ) {
        max_value = vec_w[icol];
        max_col_index = icol;
      }
    }

    // Assign and normalise direction vector
    // Note: Normalisation should be redundant because we are already dealing with unit vectors
    norm = 0.0;
    for( icol = 0 ; icol < NCOLS; icol++ ) {
      vec_out[icol] = mat_v[icol][max_col_index];
      norm += CS_SQR(vec_out[icol]);
    }
    // Work-around: If only one sample is given, SVD will have computed zero output
    // This really only applies to testing
    if( nSamples == 1 ) {
      vec_out[0] = xSamples[0];
      vec_out[1] = ySamples[0];
      vec_out[2] = zSamples[0];
      norm = 0.0;
      for( icol = 0 ; icol < NCOLS; icol++ ) {
        norm += CS_SQR(vec_out[icol]);
      }
    }

    norm = sqrt(norm);
    if( norm == 0.0 ) norm = 1.0;
    for( icol = 0 ; icol < NCOLS; icol++ ) {
      vec_out[icol] /= norm;
    }

    // fprintf(stdout,"----------firstSample: %d  nSamples: %d  norm= %f-------------\n", firstSample, nSamples, norm);
    //  fprintf(stdout,"Mean values: %12.6f %12.6f %12.6f\n", x_mean, y_mean, z_mean );

    //  printf("VEC_W  %10.4f %10.4f %10.4f   %d\n", vec_w[0],vec_w[1],vec_w[2], max_col_index);
    // printf("  MAT_V1  %10.4f %10.4f %10.4f   \n", mat_v[0][0], mat_v[1][0], mat_v[2][0]);
    // printf("  MAT_V2  %10.4f %10.4f %10.4f   \n", mat_v[0][1], mat_v[1][1], mat_v[2][1]);
    // printf("  MAT_V3  %10.4f %10.4f %10.4f   \n", mat_v[0][2], mat_v[1][2], mat_v[2][2]);
    //  fprintf(stdout,"Singular values: %12.6f %12.6f %12.6f\n", vec_w[0], vec_w[1], vec_w[2]);
    //  fprintf(stdout,"Min singular vector:  %12.6f  %12.6f  %12.6f\n", mat_v[0][min_col_index], mat_v[1][min_col_index], mat_v[2][min_col_index]);
    //  orig[1]-mat_v[1][max_col_index], orig[2]-mat_v[2][max_col_index]);
    //  fprintf(stdout,"Original vector:  %12.6f  %12.6f  %12.6f\n", xSamples[495], ySamples[495], zSamples[495] );
    //  fprintf(stdout,"Singular vector:  %12.6f  %12.6f  %12.6f\n", mat_v[0][max_col_index], mat_v[1][max_col_index], mat_v[2][max_col_index]);

    free_matrix( mat_a );
    free_matrix( mat_u );
    free_matrix( mat_v );
    free_vector( vec_w );
    free_vector( vec_work );
    return true;
  }


  bool linefit_3d_all( float* xSamples,
                      float* ySamples,
                      float* zSamples,
                      int firstSample,
                      int lastSample,
                      int force_origin,
                      float* vec_out,
                      float* vec_axes, // Size of major, minor, and second minor axes of best fit spheroid
                      float* vec_minor ) // 3D vector of minor spheroid axis
  {
    const int NCOLS = 3;
    const int ID_X = 0;
    const int ID_Y = 1;
    const int ID_Z = 2;

    float norm;
    float** mat_a;
    float** mat_u;
    float** mat_v;
    float*  vec_w;
    float*  vec_work;
    float x_mean, y_mean, z_mean;
    int max_col_index, min_col_index, med_col_index;
    int nSamples;
    int isamp, icol;
    bool ret = false;

    //------------------------------------------------------------

    nSamples = lastSample - firstSample + 1;

    //  for( isamp = firstSample; isamp <= lastSample; isamp++ ) {
    //   fprintf(stdout,"%f  %f  %d SAMPXY\n", xSamples[isamp], ySamples[isamp], isamp);
    //   fprintf(stdout,"%f  %f  %d SAMPXZ\n", xSamples[isamp], zSamples[isamp], isamp);
    //   fprintf(stdout,"%f  %f  %d SAMPYZ\n", ySamples[isamp], zSamples[isamp], isamp);
    //  }

    mat_a = allocate_matrix( nSamples, NCOLS );
    mat_u = allocate_matrix( nSamples, NCOLS );
    mat_v = allocate_matrix( NCOLS, NCOLS );
    vec_w    = allocate_vector( nSamples );
    vec_work = allocate_vector( nSamples );

    if( !mat_a || !mat_u || !mat_v || !vec_w || !vec_work) {
      return false;
    }

    //-------------------------------
    // Compute average xyz values, fill inversion matrix with reduced data point values
    // x - x_mean, y - y_mean, z - z_mean

    x_mean = y_mean = z_mean = 0.0;
    // Unless line is forced through origin, compute centroid point by averaging x,y,z:
    if( !force_origin ) {
      for( isamp = firstSample; isamp <= lastSample; isamp++ ) {
        x_mean += xSamples[isamp];
        y_mean += ySamples[isamp];
        z_mean += zSamples[isamp];
      }
      x_mean /= (float)nSamples;
      y_mean /= (float)nSamples;
      z_mean /= (float)nSamples;
    }
    for( isamp = 0; isamp < nSamples; isamp++ ) {
      mat_a[isamp][ID_X] = xSamples[firstSample+isamp] - x_mean;
      mat_a[isamp][ID_Y] = ySamples[firstSample+isamp] - y_mean;
      mat_a[isamp][ID_Z] = zSamples[firstSample+isamp] - z_mean;
      for( icol = 0; icol < NCOLS; icol++ ) {
        mat_u[isamp][icol] = mat_a[isamp][icol];
      }
    }

    //------------------------------
    // Decompose 'inversion' matrix
    //
    ret = svd_decomposition( mat_u, nSamples, NCOLS, vec_w, mat_v, vec_work );
    if( !ret ) {
      return false;
    }

    // Find min/medium/maximum singular value, corresponding singular vector is the unit direction vector

    // Assume NCOLS = 3
    if( vec_w[0] >= vec_w[1] ) {
      if( vec_w[1] >= vec_w[2] ) {
        max_col_index = 0;
        med_col_index = 1;
        min_col_index = 2;
      }
      else if( vec_w[2] >= vec_w[0] ) {
        max_col_index = 2;
        med_col_index = 0;
        min_col_index = 1;
      }
      else {
        max_col_index = 0;
        med_col_index = 2;
        min_col_index = 1;
      }
    }
    else {  // vec_w[1] >= vec_w[0]
      if( vec_w[2] >= vec_w[1] ) {
        max_col_index = 2;
        med_col_index = 1;
        min_col_index = 0;
      }
      else if( vec_w[0] >= vec_w[2] ) {
        max_col_index = 1;
        med_col_index = 0;
        min_col_index = 2;
      }
      else {
        max_col_index = 1;
        med_col_index = 2;
        min_col_index = 0;
      }
    }

    // Assign and normalise direction vector (normalisation should be redundant)
    norm = 0.0;
    for( icol = 0 ; icol < NCOLS; icol++ ) {
      vec_out[icol] = mat_v[icol][max_col_index] * vec_w[max_col_index];
      norm += CS_SQR(vec_out[icol]);
    }

    norm = sqrt(norm);
    if( norm == 0.0 ) norm = 1.0;
    for( icol = 0 ; icol < NCOLS; icol++ ) {
      vec_out[icol] /= norm;
    }

    vec_axes[0] = vec_w[max_col_index];// / norm;
    vec_axes[1] = vec_w[med_col_index];// / norm;
    vec_axes[2] = vec_w[min_col_index];// / norm;
    for( icol = 0 ; icol < NCOLS; icol++ ) {
      vec_minor[icol] = mat_v[icol][med_col_index] * vec_w[med_col_index] / norm;
      //    vec_third[icol]  = mat_v[icol][min_col_index] * vec_w[min_col_index] / norm;
    }

    // fprintf(stdout,"----------firstSample: %d  nSamples: %d  norm= %f-------------\n", firstSample, nSamples, norm);
    //  fprintf(stdout,"Mean values: %12.6f %12.6f %12.6f\n", x_mean, y_mean, z_mean );

    //  printf("VEC_W  %10.4f %10.4f %10.4f   %d\n", vec_w[0],vec_w[1],vec_w[2], max_col_index);
    // printf("  MAT_V1  %10.4f %10.4f %10.4f   \n", mat_v[0][0], mat_v[1][0], mat_v[2][0]);
    // printf("  MAT_V2  %10.4f %10.4f %10.4f   \n", mat_v[0][1], mat_v[1][1], mat_v[2][1]);
    // printf("  MAT_V3  %10.4f %10.4f %10.4f   \n", mat_v[0][2], mat_v[1][2], mat_v[2][2]);
    //  fprintf(stdout,"Singular values: %12.6f %12.6f %12.6f\n", vec_w[0], vec_w[1], vec_w[2]);
    //  fprintf(stdout,"Min singular vector:  %12.6f  %12.6f  %12.6f\n", mat_v[0][min_col_index], mat_v[1][min_col_index], mat_v[2][min_col_index]);
    //  orig[1]-mat_v[1][max_col_index], orig[2]-mat_v[2][max_col_index]);
    //  fprintf(stdout,"Original vector:  %12.6f  %12.6f  %12.6f\n", xSamples[495], ySamples[495], zSamples[495] );
    //  fprintf(stdout,"Singular vector:  %12.6f  %12.6f  %12.6f\n", mat_v[0][max_col_index], mat_v[1][max_col_index], mat_v[2][max_col_index]);

    free_matrix( mat_a );
    free_matrix( mat_u );
    free_matrix( mat_v );
    free_vector( vec_w );
    free_vector( vec_work );
    return true;
  }

//---------------------------------------------------------------------------------------------
  
  bool linefit_3d_all( float* xSamples,
                      float* ySamples,
                      float* zSamples,
                      int firstSample,
                      int lastSample,
                      int force_origin,
                      float* vec_out,
                      float* vec_axes )  // Size of major, minor, and second minor axes of best fit spheroid
  {
    const int NCOLS = 3;
    const int ID_X = 0;
    const int ID_Y = 1;
    const int ID_Z = 2;

    float norm;
    float** mat_a;
    float** mat_u;
    float** mat_v;
    float*  vec_w;
    float*  vec_work;
    float x_mean, y_mean, z_mean;
    int max_col_index, min_col_index, med_col_index;
    int nSamples;
    int isamp, icol;
    bool ret = false;

    //------------------------------------------------------------

    nSamples = lastSample - firstSample + 1;

    //  for( isamp = firstSample; isamp <= lastSample; isamp++ ) {
    //   fprintf(stdout,"%f  %f  %d SAMPXY\n", xSamples[isamp], ySamples[isamp], isamp);
    //   fprintf(stdout,"%f  %f  %d SAMPXZ\n", xSamples[isamp], zSamples[isamp], isamp);
    //   fprintf(stdout,"%f  %f  %d SAMPYZ\n", ySamples[isamp], zSamples[isamp], isamp);
    //  }

    mat_a = allocate_matrix( nSamples, NCOLS );
    mat_u = allocate_matrix( nSamples, NCOLS );
    mat_v = allocate_matrix( NCOLS, NCOLS );
    vec_w    = allocate_vector( nSamples );
    vec_work = allocate_vector( nSamples );

    if( !mat_a || !mat_u || !mat_v || !vec_w || !vec_work) {
      return false;
    }

    //-------------------------------
    // Compute average xyz values, fill inversion matrix with reduced data point values
    // x - x_mean, y - y_mean, z - z_mean

    x_mean = y_mean = z_mean = 0.0;
    // Unless line is forced through origin, compute centroid point by averaging x,y,z:
    if( !force_origin ) {
      for( isamp = firstSample; isamp <= lastSample; isamp++ ) {
        x_mean += xSamples[isamp];
        y_mean += ySamples[isamp];
        z_mean += zSamples[isamp];
      }
      x_mean /= (float)nSamples;
      y_mean /= (float)nSamples;
      z_mean /= (float)nSamples;
    }
    for( isamp = 0; isamp < nSamples; isamp++ ) {
      mat_a[isamp][ID_X] = xSamples[firstSample+isamp] - x_mean;
      mat_a[isamp][ID_Y] = ySamples[firstSample+isamp] - y_mean;
      mat_a[isamp][ID_Z] = zSamples[firstSample+isamp] - z_mean;
      for( icol = 0; icol < NCOLS; icol++ ) {
        mat_u[isamp][icol] = mat_a[isamp][icol];
      }
    }

    //------------------------------
    // Decompose 'inversion' matrix
    //
    ret = svd_decomposition( mat_u, nSamples, NCOLS, vec_w, mat_v, vec_work );
    if( !ret ) {
      return false;
    }

    // Find min/medium/maximum singular value, corresponding singular vector is the unit direction vector

    // Assume NCOLS = 3
    if( vec_w[0] >= vec_w[1] ) {
      if( vec_w[1] >= vec_w[2] ) {
        max_col_index = 0;
        med_col_index = 1;
        min_col_index = 2;
      }
      else if( vec_w[2] >= vec_w[0] ) {
        max_col_index = 2;
        med_col_index = 0;
        min_col_index = 1;
      }
      else {
        max_col_index = 0;
        med_col_index = 2;
        min_col_index = 1;
      }
    }
    else {  // vec_w[1] >= vec_w[0]
      if( vec_w[2] >= vec_w[1] ) {
        max_col_index = 2;
        med_col_index = 1;
        min_col_index = 0;
      }
      else if( vec_w[0] >= vec_w[2] ) {
        max_col_index = 1;
        med_col_index = 0;
        min_col_index = 2;
      }
      else {
        max_col_index = 1;
        med_col_index = 2;
        min_col_index = 0;
      }
    }

    // Assign and normalise direction vector (normalisation should be redundant)
    norm = 0.0;
    for( icol = 0 ; icol < NCOLS; icol++ ) {
      vec_out[icol] = mat_v[icol][max_col_index] * vec_w[max_col_index];
      norm += CS_SQR(vec_out[icol]);
    }

    norm = sqrt(norm);
    if( norm == 0.0 ) norm = 1.0;
    for( icol = 0 ; icol < NCOLS; icol++ ) {
      vec_out[icol] /= norm;
    }

    vec_axes[0] = vec_w[max_col_index];// / norm;
    vec_axes[1] = vec_w[med_col_index];// / norm;
    vec_axes[2] = vec_w[min_col_index];// / norm;
    //  for( icol = 0 ; icol < NCOLS; icol++ ) {
    //    vec_second[icol] = mat_v[icol][med_col_index] * vec_w[med_col_index] / norm;
    //    vec_third[icol]  = mat_v[icol][min_col_index] * vec_w[min_col_index] / norm;
    //  }

    // fprintf(stdout,"----------firstSample: %d  nSamples: %d  norm= %f-------------\n", firstSample, nSamples, norm);
    //  fprintf(stdout,"Mean values: %12.6f %12.6f %12.6f\n", x_mean, y_mean, z_mean );

    //  printf("VEC_W  %10.4f %10.4f %10.4f   %d\n", vec_w[0],vec_w[1],vec_w[2], max_col_index);
    // printf("  MAT_V1  %10.4f %10.4f %10.4f   \n", mat_v[0][0], mat_v[1][0], mat_v[2][0]);
    // printf("  MAT_V2  %10.4f %10.4f %10.4f   \n", mat_v[0][1], mat_v[1][1], mat_v[2][1]);
    // printf("  MAT_V3  %10.4f %10.4f %10.4f   \n", mat_v[0][2], mat_v[1][2], mat_v[2][2]);
    //  fprintf(stdout,"Singular values: %12.6f %12.6f %12.6f\n", vec_w[0], vec_w[1], vec_w[2]);
    //  fprintf(stdout,"Min singular vector:  %12.6f  %12.6f  %12.6f\n", mat_v[0][min_col_index], mat_v[1][min_col_index], mat_v[2][min_col_index]);
    //  orig[1]-mat_v[1][max_col_index], orig[2]-mat_v[2][max_col_index]);
    //  fprintf(stdout,"Original vector:  %12.6f  %12.6f  %12.6f\n", xSamples[495], ySamples[495], zSamples[495] );
    //  fprintf(stdout,"Singular vector:  %12.6f  %12.6f  %12.6f\n", mat_v[0][max_col_index], mat_v[1][max_col_index], mat_v[2][max_col_index]);

    free_matrix( mat_a );
    free_matrix( mat_u );
    free_matrix( mat_v );
    free_vector( vec_w );
    free_vector( vec_work );
    return true;
  }

  //---------------------------------------------------------------------------------------------
  //
  //
  bool linefit_3d_2step( float* xSamples,
                        float* ySamples,
                        float* zSamples,
                        int nSamples,
                        int method,
                        int force_origin,
                        float* vec_out )
  {
    const int NCOLS = 3;

    float norm;
    float** mat_a;
    float** mat_u;
    float** mat_v;
    float*  vec_w;
    float*  vec_work;
    float*  vec_b;
    float*  vec_x1;
    float*  vec_x2;
    float x_mean, y_mean, z_mean;
    int isamp, icol;
    float ssxx, ssxy, ssxy2;
    bool ret = false;

    //------------------------------------------------------------

    mat_a = allocate_matrix( nSamples, 2 );
    mat_u = allocate_matrix( nSamples, 2 );
    mat_v = allocate_matrix( 2, 2 );
    vec_w    = allocate_vector( nSamples );
    vec_work = allocate_vector( nSamples );
    vec_b    = allocate_vector( nSamples );
    vec_x1   = allocate_vector( 2 );
    vec_x2   = allocate_vector( 2 );

    if( !mat_a || !mat_u || !mat_v || !vec_w || !vec_work || !vec_b || !vec_x1 || !vec_x2 ) {
      return false;
    }

    //-------------------------------
    // Compute average xyz values
    //

    x_mean = y_mean = z_mean = 0.0;
    if( !force_origin ) {
      for( isamp = 0; isamp < nSamples; isamp++ ) {
        x_mean += xSamples[isamp];
        y_mean += ySamples[isamp];
        z_mean += zSamples[isamp];
      }
      x_mean /= (float)nSamples;
      y_mean /= (float)nSamples;
      z_mean /= (float)nSamples;
    }

    if( method == LINEFIT_LS ) {
      //--------------------------------------------------------------------------------
      // (1) Ordinary least-square's method
      //
      ssxy = 0.0;
      ssxy2 = 0.0;
      ssxx = 0.0;
      for( isamp = 0; isamp < nSamples; isamp++ ) {
        ssxy += xSamples[isamp]*zSamples[isamp];
        ssxy2 += xSamples[isamp]*ySamples[isamp];
        ssxx += xSamples[isamp]*xSamples[isamp];
      }
      //  printf("y(x): a b  %10.4f  %10.4f\n", ssxy/ssxx, y_mean - ssxy/ssxx*x_mean);
      if( fabs(ssxx) > 1.0e-09 ) {
        ssxx -= nSamples*x_mean*x_mean;
        ssxy -= nSamples*x_mean*z_mean;
        ssxy2 -= nSamples*x_mean*y_mean;
        vec_out[0] = 1.0;
        vec_out[1] = vec_out[0] * (y_mean - ssxy/ssxx*x_mean) + ssxy/ssxx;
        vec_out[2] = vec_out[0] * (z_mean - ssxy2/ssxx*x_mean) + ssxy2/ssxx;
      }
      else {
        // X component is zero, try to solve for line in zy plane
        ssxy = 0.0;
        ssxy2 = 0.0;
        ssxx = 0.0;
        for( isamp = 0; isamp < nSamples; isamp++ ) {
          ssxy  += ySamples[isamp]*zSamples[isamp];
          ssxx  += ySamples[isamp]*ySamples[isamp];
        }
        if( fabs(ssxx) > 1.0e-09 ) {
          ssxx  -= nSamples*y_mean*y_mean;
          ssxy  -= nSamples*y_mean*x_mean;
          ssxy2 -= nSamples*y_mean*z_mean;
          vec_out[0] = 0.0;
          vec_out[1] = 1.0;
          vec_out[2] = vec_out[1] * (z_mean - ssxy/ssxx*y_mean) + ssxy/ssxx;
        }
        else {
          // X and Y components are both zero --> pathological case
          vec_out[0] = 0.0;
          vec_out[1] = 0.0;
          vec_out[2] = 1.0;
        }
      }
    }
    else { // if( method == LINEFIT_SVD ) {
      //--------------------------------------------------------------------------------
      // (2) Solve using SVD
      //
      // (a) Solve z(x) = ax + b
      //
      for( isamp = 0; isamp < nSamples; isamp++ ) {
        vec_b[isamp] = zSamples[isamp] - z_mean;
      }
      for( isamp = 0; isamp < nSamples; isamp++ ) {
        mat_u[isamp][0] = xSamples[isamp] - x_mean;
        mat_u[isamp][1] = 1.0;
        if( force_origin )
          mat_u[isamp][1] = 0.0;  // Force line to go through origin (y = a*x + b, b=0)
      }
      ret = svd_decomposition( mat_u, nSamples, 2, vec_w, mat_v, vec_work );
      if( !ret ) {
        return false;
      }
      //if( vec_w[0] > vec_w[1] ) {
      //  if( vec_w[1]/vec_w[0] < 0.0001 ) vec_w[1] = 0.0;
      // }
      //else if( vec_w[0]/vec_w[1] < 0.0001 ) vec_w[0] = 0.0;
      svd_linsolve( mat_u, vec_w, mat_v, nSamples, 2, vec_b, vec_x2, vec_work );
      //--------------------------------------------------------------------------------
      // (b) Solve y(x) = cx + d
      //
      for( isamp = 0; isamp < nSamples; isamp++ ) {
        vec_b[isamp] = ySamples[isamp] - y_mean;
      }
      for( isamp = 0; isamp < nSamples; isamp++ ) {
        mat_u[isamp][0] = xSamples[isamp] - x_mean;
        mat_u[isamp][1] = 1.0;
        if( force_origin )
          mat_u[isamp][1] = 0.0;  // Force line to go through origin (y = a*x + b, b=0)
      }
      if( !svd_decomposition( mat_u, nSamples, 2, vec_w, mat_v, vec_work ) ) {
        return false;
      }
      svd_linsolve( mat_u, vec_w, mat_v, nSamples, 2, vec_b, vec_x1, vec_work );
      //--------------------------------------------------------------------------------
      vec_out[0] = 1.0;
      vec_out[1] = vec_out[0] * vec_x1[0] + vec_x1[1];
      vec_out[2] = vec_out[0] * vec_x2[0] + vec_x2[1];

    } // END solve by SVD


    if( force_origin )
      printf("----------------FORCE--------------------\n");

    // Assign and normalise direction vector (normalisation should be redundant)
    norm = 0.0;
    for( icol = 0 ; icol < NCOLS; icol++ ) {
      norm += CS_SQR(vec_out[icol]);
    }
    norm = sqrt(norm);
    for( icol = 0 ; icol < NCOLS; icol++ ) {
      vec_out[icol] /= norm;
    }

    //  fprintf(stdout,"Mean values: %12.6f %12.6f %12.6f\n", x_mean, y_mean, z_mean );

    free_matrix( mat_a );
    free_matrix( mat_u );
    free_matrix( mat_v );
    free_vector( vec_w );
    free_vector( vec_work );
    free_vector( vec_x1 );
    free_vector( vec_x2 );
    free_vector( vec_b );
    return false;
  }


} // namespace

