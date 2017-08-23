//
//    File: matrix3.h
//
//    (C) 2000-2008 Helmut Cantzler
//
//    Licensed under the terms of the Lesser General Public License.
//

#ifndef _MATRIX3_H
#define _MATRIX3_H

#include <stdio.h>
#include <stdarg.h>
#include <math.h>

#define MSIZE 3

template <class T>
class Matrix3
{
public:
  Matrix3() {
    for (int i=0; i < MSIZE; i++)
      for (int j=0; j < MSIZE; j++)
	data[i][j]=0.0;
  }
  Matrix3(const T d[MSIZE][MSIZE]) {
    for (int i=0; i < MSIZE; i++)
      for (int j=0; j < MSIZE; j++)
	data[i][j]=d[i][j];
  }
  Matrix3(const T *d) {
    for (int i=0; i < MSIZE; i++)
      for (int j=0; j < MSIZE; j++)
	data[i][j]=*d++;
  }

  /*
  Matrix3(T first, ...) {
    vaList ap;
    vaStart(ap, first);
    for (int i=0; i < MSIZE; i++)
      for (int j=0; j < MSIZE; j++)
	data[i][j]=vaArg(ap, T);
    vaEnd(ap);
  }
  */

  // display data
  void display(void) const {
    for (int i=0; i < MSIZE; i++)
      {
	for (int j=0; j < MSIZE; j++)
	  printf("%f ", data[i][j]);
	printf("\n");
      }
    printf("\n");
  }

  // make identity matrix
  void setIdentity(void) {
    for (int i=0; i < MSIZE; i++)
      for (int j=0; j < MSIZE; j++)
	data[i][j]= i == j ? 1.0 : 0.0;
  }

  // make rotation matrix
  void setRotationFromAngles(float x, float y, float z) {
    data[0][0] = cos(y)*cos(z);
    data[0][1] = -cos(y)*sin(z); 
    data[0][2] = sin(y);
    data[1][0] = sin(x)*sin(y)*cos(z) + cos(x)*sin(z);
    data[1][1] = -sin(x)*sin(y)*sin(z) + cos(x)*cos(z);
    data[1][2] = -sin(x)*cos(y);
    data[2][0] = -cos(x)*sin(y)*cos(z) + sin(x)*sin(z);
    data[2][1] = cos(x)*sin(y)*sin(z) + sin(x)*cos(z);
    data[2][2] = cos(x)*cos(y);
  }

  // make rotation matrix
  void setRotationFromQuaternion(float x, float y, float z, float w) {
    float xx = x*x, xy = x*y, xz = x*z, xw = x*w;
    float yy = y*y, yz = y*z, yw = y*w;
    float zz = z*z, zw = z*w;

    data[0][0] = 1.0 - 2.0 * ( yy + zz );
    data[0][1] =       2.0 * ( xy - zw );
    data[0][2] =       2.0 * ( xz + yw );
    data[1][0] =       2.0 * ( xy + zw );
    data[1][1] = 1.0 - 2.0 * ( xx + zz );
    data[1][2] =       2.0 * ( yz - xw );
    data[2][0] =       2.0 * ( xz - yw );
    data[2][1] =       2.0 * ( yz + xw );
    data[2][2] = 1.0 - 2.0 * ( xx + yy );
  }

  // get rotation angles from rotation matrix
  void getRotation(float *x, float *y, float *z) const {
    *x = atan( (-get(1,2)) / get(2,2) );
    *y = asin( get(0,2) );
    *z = atan( (-get(0,1)) / get(0,0) );
  }

  // add matrix
  void add(const Matrix3<T> *d) {
    for (int i=0; i < MSIZE; i++)
      for (int j=0; j < MSIZE; j++)
	data[i][j]+=d->data[i][j];
  }

  // multiply with scalar
  void mul(T s) {
    for (int i=0; i < MSIZE; i++)
      for (int j=0; j < MSIZE; j++)
	data[i][j]*=s;
  }

  // multiply with vector
  void mul(const T *vector, T *result) const {
    for (int i=0; i < MSIZE; i++)
      {
	result[i]=0;
	for (int j=0; j < MSIZE; j++)
	  result[i]+=vector[j]*data[i][j];
      }
  }

  // multiply with matrix
  void mul(const Matrix3<T> *d) {
    T tmp[MSIZE][MSIZE];
    for (int i=0; i < MSIZE; i++)
      for (int j=0; j < MSIZE; j++)
      {
	tmp[i][j] = 0.0;
	for (int l=0; l < MSIZE; l++)
	  tmp[i][j] += data[i][l] * d->data[l][j];
      }
    for (int i=0; i < MSIZE; i++)
      for (int j=0; j < MSIZE; j++)
	data[i][j] = tmp[i][j];
  }

  T get(int s1, int s2) const {
    return data[s1][s2];
  }

private:
  T data[MSIZE][MSIZE];
};

#endif
