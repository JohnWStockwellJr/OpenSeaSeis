/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include <cmath>
#include <cstdlib>
#include <string>
#include "geolib_math.h"

namespace {
  int N_ITERATIONS = 30;
}

//using namespace cseis;
namespace cseis_geolib {

  float pythag(float a, float b)
    {
      float absa,absb;
      absa = fabs(a);
      absb = fabs(b);
      if( absa > absb ) return absa*sqrt( 1.0 + CS_SQR(absb/absa) );
      else return (absb == 0.0 ? 0.0 : absb*sqrt( 1.0 + CS_SQR(absa/absb) ) );
    }

  /*--------------------------------------------------------------------------------
   * Singular value decomposition, numerical recipes
   *
   * Return false if maximum number of iterations is reached
   *
   */

  bool svd_decomposition( float** mat_a, int nrows, int ncols, float vec_w[], float** mat_v, float vec_work[] ) {
    int flag,i,its,j,jj,k,val,nm;
    float anorm,c,f,g,h,s,scale,x,y,z;
    nm = 0;
  
    val = 0;
    g = scale = anorm = 0.0;
    for( i = 0; i <= ncols-1; i++ ) {
      val = i+1;
      vec_work[i]=scale*g;
      g=s=scale=0.0;
      if (i <= nrows-1) {
        for( k = i; k <= nrows-1; k++ ) scale += fabs(mat_a[k][i]);
        if( scale ) {
          for( k = i; k <= nrows-1; k++ ) {
            mat_a[k][i] /= scale;
            s += mat_a[k][i]*mat_a[k][i];
          }
          f=mat_a[i][i];
          g = -SPECIAL_SIGN(sqrt(s),f);
          h=f*g-s;
          mat_a[i][i]=f-g;
          for (j=val;j<=ncols-1;j++) {
            for (s=0.0,k=i;k<=nrows-1;k++) s += mat_a[k][i]*mat_a[k][j];
            f=s/h;
            for (k=i;k<=nrows-1;k++) mat_a[k][j] += f*mat_a[k][i];
          }
          for (k=i;k<=nrows-1;k++) mat_a[k][i] *= scale;
        }
      }
      vec_w[i]=scale *g;
      g=s=scale=0.0;
      if (i <= nrows-1 && i != ncols-1) {
        for (k=val;k<=ncols-1;k++) scale += fabs(mat_a[i][k]);
        if (scale) {
          for (k=val;k<=ncols-1;k++) {
            mat_a[i][k] /= scale;
            s += mat_a[i][k]*mat_a[i][k];
          }
          f=mat_a[i][val];
          g = -SPECIAL_SIGN(sqrt(s),f);
          h=f*g-s;
          mat_a[i][val]=f-g;
          for( k = val; k <= ncols-1; k++ ) vec_work[k] = mat_a[i][k]/h;
          for( j = val; j <= nrows-1; j++ ) {
            for( s = 0.0, k = val; k <= ncols-1; k++ ) s += mat_a[j][k]*mat_a[i][k];
            for( k = val; k <= ncols-1; k++ ) mat_a[j][k] += s*vec_work[k];
          }
          for( k = val; k <= ncols-1; k++ ) mat_a[i][k] *= scale;
        }
      }
      anorm = std::max(anorm,(float)(fabs(vec_w[i])+fabs(vec_work[i])));
    }
    for( i = ncols-1; i >= 0; i-- ) {
      if (i < ncols-1) {
        if (g) {
          for( j = val; j <= ncols-1; j++ )
            mat_v[j][i]=(mat_a[i][j]/mat_a[i][val])/g;
          for( j = val; j <= ncols-1; j++ ) {
            for( s = 0.0,k = val; k <= ncols-1; k++ ) s += mat_a[i][k]*mat_v[k][j];
            for( k = val; k <= ncols-1; k++ ) mat_v[k][j] += s*mat_v[k][i];
          }
        }
        for( j = val; j <= ncols-1; j++ ) mat_v[i][j]=mat_v[j][i]=0.0;
      }
      mat_v[i][i]=1.0;
      g=vec_work[i];
      val=i;
    }
    for( i = std::min(nrows-1,ncols-1); i >= 0; i-- ) {
      val=i+1;
      g=vec_w[i];
      for( j = val; j <= ncols-1; j++ ) mat_a[i][j]=0.0;
      if (g) {
        g=1.0/g;
        for( j = val; j <= ncols-1; j++ ) {
          for( s = 0.0, k = val; k <= nrows-1; k++ ) s += mat_a[k][i]*mat_a[k][j];
          f=(s/mat_a[i][i])*g;
          for( k = i; k <= nrows-1; k++ ) mat_a[k][j] += f*mat_a[k][i];
        }
        for( j = i; j <= nrows-1; j++ ) mat_a[j][i] *= g;
      } else for( j = i; j <= nrows-1; j++ ) mat_a[j][i]=0.0;
      ++mat_a[i][i];
    }
    for( k = ncols-1; k >= 0; k-- ) {
      for( its = 1; its <= N_ITERATIONS; its++ ) {
        flag=1;
        for( val = k; val >= 0; val-- ) {
          nm=val-1;
          if ((float)(fabs(vec_work[val])+anorm) == anorm) {
            flag=0;
            break;
          }
          if ((float)(fabs(vec_w[nm])+anorm) == anorm) break;
        }
        if (flag) {
          c=0.0;
          s=1.0;
          for( i = val; i <= k; i++ ) {
            f=s*vec_work[i];
            vec_work[i]=c*vec_work[i];
            if ((float)(fabs(f)+anorm) == anorm) break;
            g=vec_w[i];
            h=pythag(f,g);
            vec_w[i]=h;
            h=1.0/h;
            c=g*h;
            s = -f*h;
            for( j = 0; j <= nrows-1; j++ ) {
              y=mat_a[j][nm];
              z=mat_a[j][i];
              mat_a[j][nm]=y*c+z*s;
              mat_a[j][i]=z*c-y*s;
            }
          }
        }
        z=vec_w[k];
        if (val == k) {
          if (z < 0.0) {
            vec_w[k] = -z;
            for( j = 0; j <= ncols-1; j++ ) mat_v[j][k] = -mat_v[j][k];
          }
          break;
        }
        if (its == N_ITERATIONS) return false;
        x=vec_w[val];
        nm=k-1;
        y=vec_w[nm];
        g=vec_work[nm];
        h=vec_work[k];
        f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
        g=pythag(f,1.0);
        f=((x-z)*(x+z)+h*((y/(f+SPECIAL_SIGN(g,f)))-h))/x;
        c=s=1.0;
        for( j = val; j <= nm; j++ ) {
          i=j+1;
          g=vec_work[i];
          y=vec_w[i];
          h=s*g;
          g=c*g;
          z=pythag(f,h);
          vec_work[j]=z;
          c=f/z;
          s=h/z;
          f=x*c+g*s;
          g = g*c-x*s;
          h=y*s;
          y *= c;
          for( jj = 0; jj <= ncols-1; jj++ ) {
            x=mat_v[jj][j];
            z=mat_v[jj][i];
            mat_v[jj][j]=x*c+z*s;
            mat_v[jj][i]=z*c-x*s;
          }
          z=pythag(f,h);
          vec_w[j]=z;
          if (z) {
            z=1.0/z;
            c=f*z;
            s=h*z;
          }
          f=c*g+s*y;
          x=c*y-s*g;
          for( jj = 0; jj <= nrows-1; jj++ ) {
            y=mat_a[jj][j];
            z=mat_a[jj][i];
            mat_a[jj][j]=y*c+z*s;
            mat_a[jj][i]=z*c-y*s;
          }
        }
        vec_work[val]=0.0;
        vec_work[k]=f;
        vec_w[k]=x;
      }
    }
    return true;
  }


  //********************************************************************************
  //********************************************************************************
  //********************************************************************************

  double pythag(double a, double b)
    {
      double absa,absb;
      absa = fabs(a);
      absb = fabs(b);
      if( absa > absb ) return absa*sqrt( 1.0 + CS_SQR(absb/absa) );
      else return (absb == 0.0 ? 0.0 : absb*sqrt( 1.0 + CS_SQR(absa/absb) ) );
    }

  /*--------------------------------------------------------------------------------
   * Singular value decomposition, numerical recipes
   *
   * Return false if maximum number of iterations is reached
   *
   */

  bool svd_decomposition( double** mat_a, int nrows, int ncols, double vec_w[], double** mat_v, double vec_work[] ) {
    int flag,i,its,j,jj,k,val,nm;
    double anorm,c,f,g,h,s,scale,x,y,z;
    nm = 0;
  
    val = 0;
    g = scale = anorm = 0.0;
    for( i = 0; i <= ncols-1; i++ ) {
      val = i+1;
      vec_work[i]=scale*g;
      g=s=scale=0.0;
      if (i <= nrows-1) {
        for( k = i; k <= nrows-1; k++ ) scale += fabs(mat_a[k][i]);
        if( scale ) {
          for( k = i; k <= nrows-1; k++ ) {
            mat_a[k][i] /= scale;
            s += mat_a[k][i]*mat_a[k][i];
          }
          f=mat_a[i][i];
          g = -SPECIAL_SIGN(sqrt(s),f);
          h=f*g-s;
          mat_a[i][i]=f-g;
          for (j=val;j<=ncols-1;j++) {
            for (s=0.0,k=i;k<=nrows-1;k++) s += mat_a[k][i]*mat_a[k][j];
            f=s/h;
            for (k=i;k<=nrows-1;k++) mat_a[k][j] += f*mat_a[k][i];
          }
          for (k=i;k<=nrows-1;k++) mat_a[k][i] *= scale;
        }
      }
      vec_w[i]=scale *g;
      g=s=scale=0.0;
      if (i <= nrows-1 && i != ncols-1) {
        for (k=val;k<=ncols-1;k++) scale += fabs(mat_a[i][k]);
        if (scale) {
          for (k=val;k<=ncols-1;k++) {
            mat_a[i][k] /= scale;
            s += mat_a[i][k]*mat_a[i][k];
          }
          f=mat_a[i][val];
          g = -SPECIAL_SIGN(sqrt(s),f);
          h=f*g-s;
          mat_a[i][val]=f-g;
          for( k = val; k <= ncols-1; k++ ) vec_work[k] = mat_a[i][k]/h;
          for( j = val; j <= nrows-1; j++ ) {
            for( s = 0.0, k = val; k <= ncols-1; k++ ) s += mat_a[j][k]*mat_a[i][k];
            for( k = val; k <= ncols-1; k++ ) mat_a[j][k] += s*vec_work[k];
          }
          for( k = val; k <= ncols-1; k++ ) mat_a[i][k] *= scale;
        }
      }
      anorm = std::max(anorm,(fabs(vec_w[i])+fabs(vec_work[i])));
    }
    for( i = ncols-1; i >= 0; i-- ) {
      if (i < ncols-1) {
        if (g) {
          for( j = val; j <= ncols-1; j++ )
            mat_v[j][i]=(mat_a[i][j]/mat_a[i][val])/g;
          for( j = val; j <= ncols-1; j++ ) {
            for( s = 0.0,k = val; k <= ncols-1; k++ ) s += mat_a[i][k]*mat_v[k][j];
            for( k = val; k <= ncols-1; k++ ) mat_v[k][j] += s*mat_v[k][i];
          }
        }
        for( j = val; j <= ncols-1; j++ ) mat_v[i][j]=mat_v[j][i]=0.0;
      }
      mat_v[i][i]=1.0;
      g=vec_work[i];
      val=i;
    }
    for( i = IMIN(nrows-1,ncols-1); i >= 0; i-- ) {
      val=i+1;
      g=vec_w[i];
      for( j = val; j <= ncols-1; j++ ) mat_a[i][j]=0.0;
      if (g) {
        g=1.0/g;
        for( j = val; j <= ncols-1; j++ ) {
          for( s = 0.0, k = val; k <= nrows-1; k++ ) s += mat_a[k][i]*mat_a[k][j];
          f=(s/mat_a[i][i])*g;
          for( k = i; k <= nrows-1; k++ ) mat_a[k][j] += f*mat_a[k][i];
        }
        for( j = i; j <= nrows-1; j++ ) mat_a[j][i] *= g;
      } else for( j = i; j <= nrows-1; j++ ) mat_a[j][i]=0.0;
      ++mat_a[i][i];
    }
    for( k = ncols-1; k >= 0; k-- ) {
      for( its = 1; its <= N_ITERATIONS; its++ ) {
        flag=1;
        for( val = k; val >= 0; val-- ) {
          nm=val-1;
          if ((double)(fabs(vec_work[val])+anorm) == anorm) {
            flag=0;
            break;
          }
          if( nm < 0 ) return false;   // This can occur if there are inconsistencies in the input data
          if ((double)(fabs(vec_w[nm])+anorm) == anorm) break;
        }
        if (flag) {
          c=0.0;
          s=1.0;
          for( i = val; i <= k; i++ ) {
            f=s*vec_work[i];
            vec_work[i]=c*vec_work[i];
            if ((double)(fabs(f)+anorm) == anorm) break;
            g=vec_w[i];
            h=pythag(f,g);
            vec_w[i]=h;
            h=1.0/h;
            c=g*h;
            s = -f*h;
            for( j = 0; j <= nrows-1; j++ ) {
              y=mat_a[j][nm];
              z=mat_a[j][i];
              mat_a[j][nm]=y*c+z*s;
              mat_a[j][i]=z*c-y*s;
            }
          }
        }
        z=vec_w[k];
        if (val == k) {
          if (z < 0.0) {
            vec_w[k] = -z;
            for( j = 0; j <= ncols-1; j++ ) mat_v[j][k] = -mat_v[j][k];
          }
          break;
        }
        if (its == N_ITERATIONS) return false;
        x=vec_w[val];
        nm=k-1;
        y=vec_w[nm];
        g=vec_work[nm];
        h=vec_work[k];
        f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
        g=pythag(f,1.0);
        f=((x-z)*(x+z)+h*((y/(f+SPECIAL_SIGN(g,f)))-h))/x;
        c=s=1.0;
        for( j = val; j <= nm; j++ ) {
          i=j+1;
          g=vec_work[i];
          y=vec_w[i];
          h=s*g;
          g=c*g;
          z=pythag(f,h);
          vec_work[j]=z;
          c=f/z;
          s=h/z;
          f=x*c+g*s;
          g = g*c-x*s;
          h=y*s;
          y *= c;
          for( jj = 0; jj <= ncols-1; jj++ ) {
            x=mat_v[jj][j];
            z=mat_v[jj][i];
            mat_v[jj][j]=x*c+z*s;
            mat_v[jj][i]=z*c-x*s;
          }
          z=pythag(f,h);
          vec_w[j]=z;
          if (z) {
            z=1.0/z;
            c=f*z;
            s=h*z;
          }
          f=c*g+s*y;
          x=c*y-s*g;
          for( jj = 0; jj <= nrows-1; jj++ ) {
            y=mat_a[jj][j];
            z=mat_a[jj][i];
            mat_a[jj][j]=y*c+z*s;
            mat_a[jj][i]=z*c-y*s;
          }
        }
        vec_work[val]=0.0;
        vec_work[k]=f;
        vec_w[k]=x;
      }
    }
    return true;
  }

} // namespace cseis

