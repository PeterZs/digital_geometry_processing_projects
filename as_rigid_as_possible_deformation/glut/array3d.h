//
//    File: array3d.h
//
//    (C) 2000-2008 Helmut Cantzler
//
//    Licensed under the terms of the Lesser General Public License.
//

#ifndef _ARRAY3D_H
#define _ARRAY3D_H

template <class T>
class Array3D
{
 public:

  Array3D()
    {
      sizeX=sizeY=sizeZ=0;
      data=NULL;
    }

  Array3D(int sx, int sy, int sz)
    {
      sizeX=sx; sizeY=sy; sizeZ=sz;
      data= new T [sizeX * sizeY * sizeZ];
    }

  ~Array3D()
    {
      delete [] data;
    }

  const T* getData(void) const
    {
      return data;
    }

  T get(int x, int y, int z) const
    {
      return data[(y * sizeX + x) * sizeZ + z];
    } 

  void set(int x, int y, int z, T d)
    {
      data[(y * sizeX + x) * sizeZ + z] = d;
    }

  int getSizeX() const
    {
      return sizeX;
    }

  int getSizeY() const
    {
      return sizeY;
    }

  int getSizeZ() const
    {
      return sizeZ;
    }

 private:
  T *data;
  int sizeX, sizeY, sizeZ;
};

#endif
