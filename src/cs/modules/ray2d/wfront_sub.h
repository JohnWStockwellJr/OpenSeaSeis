/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef WFRONT_SUB_H
#define WFRONT_SUB_H

#include <cstdio>
#include "csMatrixFStyle.h"

extern "C" {
  void model_( int& n_int, int& nboreholes, int& numIntAll, int* npoints_int, float* x_int, float* z_int, int* iii,
               int* nx_grid, int* nz_grid, float* x_grid, float* z_grid, float* v_grid,
               int* index_grid,
               float* v_top, float* v_bottom, float* vpvs_ratio, int& f_out,
               int& flag_smooth,
               float* b_int, float* c_int, float* d_int,
               float& left_border, float& right_border,
               int* iii2, int& error,
               int& MAXPOINTS_INT, int& MAXP_XGRID, int& MAXP_ZGRID );

  void get_z_( float& xin, float* x_int, float* z_int, float* d_int, float*c_int, float* b_int, int& npoints,
               float& zout, int& actp,
               int& MAXPOINTS_INT );

  void velocity_( float& x_source, float& z_source, int& nx_grid, int& nz_grid,
                  float* x_grid, float* z_grid, float* v_grid,
                  int& veltype, float& vpvs_ratio, int& lay_source,
                  float* vel,
                  int& MAXP_XGRID, int& MAXP_ZGRID );

  void velonly_( float& x0, float& z0, int& nx_grid, int& nz_grid,
                 float* x_grid, float* z_grid, float* v_grid,
                 int& veltype, float& vpvs_ratio, int& actlayer,
                 float& vel,
                 int& MAXP_XGRID, int& MAXP_ZGRID);
 
  void sortraycode_( int* numCodeSteps, int& numCodes, int* codeStepPlot, int& f_out,
                     int& int_source, int& lay_source, int* code, int* hold, int& flag_dup, int& error,
                     int& maxNumCodeSteps);

  void check_sou_loc_( float* x_int, float* z_int, float* d_int, float* c_int, float* b_int, int* numPointsInt,
                       float& x_source, float& z_source, float& dzdx_source,
                       int& int_source, int& lay_source, int& error, int& f_out,
                       int& numInt, int& maxPointsInt );

  void wfsub_( int& numInt, int* numPointsInt, float* x_int, float* z_int, int* iii, float* rho1, float* rho2, float* vpvs_ratio,
               int* nx_grid, int* nz_grid, float* x_grid, float* z_grid, float* v_grid,
               float* b_int, float* c_int, float* d_int,
               float& left_border, float& right_border, int* iii2,
               int& nray_source, int& lay_source, int& int_source, float& x_source, float& z_source, float& t_source, float& dzdx_source,
               float& anglein_rad, float& angleout_rad,
               float& dtray, float& spread_max, float& angle_max,
               int& flag_smooth, int& flag_compounds, int& flag_surface, int& flag_stack, int& flag_2point,
               int& nreclines, int* int_rec, int* lay_rec, int* nrec, float* xz_rec, int& numBoreholes,
               int& step_spread, int& step_wfronts, int& step1_wfronts, int& last_wfronts,
               int* code, int* numCodeSteps, int& numCodes, int* codeStepPlot, int* hold,
               int& f_out, int& f_wfronts, int& f_xrays, int& f_rayout, int& f_rayouttmp, int& f_timeout,
               int& maxNumReceivers, int& maxNumCodeSteps,
               int& maxPointsXGrid, int& maxPointsZGrid, int& maxAngles2Point,
               int& numIntAll, int& maxPointsInt, int& maxNumArrivals,
               float* time, float* amplitude, float* phase, float* angle, int& comp_out, int& numTimeCodes, float* ccp );
}

int wfront_check_params(
          int n_int,
          cseis_geolib::csMatrixFStyle<float> const& x_int,
          int* npoints_int,
          int* index_grid,
          cseis_geolib::csMatrixFStyle<float> const& x_grid,
          cseis_geolib::csMatrixFStyle<float> const& z_grid,
          int* nx_grid,
          int* nz_grid,
          int lay_source, int int_source,
          float anglein, float angleout,
          int nray_source,
          int nreclines,
          int* lay_rec,
          int* int_rec,
          int f_wfronts,
          FILE* f_out,
          int step_wfronts,
          int ncodes,
          cseis_geolib::csMatrixFStyle<int> const& code,
          int* numCodeSteps,
          int flag_2point,
          int flag_stack );

bool addRickerWavelet( float freq_hz, float dampingFactor, float phaseShift_rad, float sampleInt_s, float timeSamp1_s,
                       float eventTime_s, float eventPhase_rad, float eventAmplitude, float* samples, int nSamples );

#endif

