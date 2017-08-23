//
//    File: edge.h
//
//    (C) 2000-2008 Helmut Cantzler
//
//    Licensed under the terms of the Lesser General Public License.
//

#ifndef _EDGE_H
#define _EDGE_H

#include <list>

#include "triangle.h"
#include "vertex.h"
#include "mathvector.h"

using namespace std;

class Vertex;
class Triangle;

class Edge
{
 public:
  Edge(Vertex *v1, Vertex *v2);
  ~Edge();

  void set(const MathVector *v1, const MathVector *v2);
  void calcProperties(void);
  // void addTriangle(Triangle* triangle);
  // void deleteTriangle(Triangle* triangle);
  int changeVertex(const Vertex *oldV, Vertex *newV);

  int onSameEdge(const Edge *e, float tolerance) const;
  int equal(const Edge *e) const;
  int equal(const Vertex *v1, const Vertex *v2) const;

  float length(void) const;
  const float* centroid(void) const;
  const MathVector* mathCentroid(void) const;
  const MathVector* mathOrientation(void) const;

  Edge *makeFirstTriangle(const Triangle*);

  Vertex *vertices[2];
  Triangle *triangles[2];
  //list<Triangle*> triangles;

  unsigned int nrTri, number, name;

 private:
  float edgeLength;
  MathVector edgeCentroid, edgeOrientation;
};

#endif

