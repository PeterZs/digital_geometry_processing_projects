//
//    File: triangle.h
//
//    (C) 2000-2008 Helmut Cantzler
//
//    Licensed under the terms of the Lesser General Public License.
//

#ifndef _TRIANGLE_H
#define _TRIANGLE_H

#include <vector>
#include "edge.h"
#include "vertex.h"
#include "mathvector.h"
#include "misc.h"

using namespace std;

enum TriangleType{ TRIANGLE, POLYGON_TRIANGLE, SENTINEL_TRIANGLE };

class Triangle
{
public:
  Triangle(Vertex *v1, Vertex *v2, Vertex *v3, TriangleType type = TRIANGLE);
  //Triangle(Edge *e1, Edge *e2, Edge *e3);
  ~Triangle();

  //  void set(MathVector *v1, MathVector *v2, MathVector *v3);
  void calcProperties(void);
  void setVertices(const MathVector *v1, const MathVector *v2, 
		   const MathVector *v3);
  int changeVertex(const Vertex *oldV, Vertex *newV);
  void moveCentroid(const MathVector *v);
  void negateNormal(void);

  void setTextCoordinates(int i, float s, float t);
  void setTextCoordinates(int i, const pair<float, float> *texCoord);
  float getTextT(int i) const;
  float getTextS(int i) const;

  const float* centroid(void) const;
  const MathVector* mathCentroid(void) const;
  const float* floatNormal(void) const;
  const MathVector* mathNormal(void) const;

  float size(void) const;
  float perimeter(void) const;
  float distanceToOrigin(void) const;
  float angle(const Vertex *v) const;

  bool equal(const Triangle *tri) const;
  bool onSamePlane(const Triangle *tri, float distanceTolerance,
		  float orientationTolerance) const;
  bool isFromPolygon(void) const;

  vector<Triangle*> neighbors(void) const;
  vector<Triangle*> neighborsOnPlane(void) const;
  int getSurroundingPlane(void) const;

  Vertex* nextVertex(Vertex *v);
  Vertex* prevVertex(Vertex *v);
  static Triangle* getSentinel() {return &sentinel;}

  Edge *edges[3];
  Vertex *vertices[3];
  unsigned int plane, name;

private:
  MathVector surfaceNormal, surfaceCentroid;
  float surfaceSize, surfacePerimeter, surfaceDistanceToOrigin;
  bool fromPolygon;
  TriangleType type;
  static Triangle sentinel;
};

#endif
