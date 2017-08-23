//
//    File: vertex.h
//
//    (C) 2000-2008 Helmut Cantzler
//
//    Licensed under the terms of the Lesser General Public License.
//

#ifndef _VERTEX_H
#define _VERTEX_H

#include <list>
#include <utility>

#include "triangle.h"
#include "edge.h"
#include "matrix3.h"
#include "mathvector.h"

using namespace std;

class Triangle;
class Edge;

class Vertex
{
public:
  Vertex();
  Vertex(const Vertex &v);
  Vertex(float x, float y, float z);
  ~Vertex();

  void addTriangle(Triangle* triangle);
  void deleteTriangle(Triangle* triangle);
  void addEdge(Edge* edge);
  void deleteEdge(Edge* edge);

  const float& x(void) const{return data.v[0];}
  const float& y(void) const{return data.v[1];}
  const float& z(void) const{return data.v[2];}
  float& x(void) {return data.v[0];}
  float& y(void) {return data.v[1];}
  float& z(void) {return data.v[2];}

  const float* floatData(void) const;
  const MathVector* mathData(void) const;

  const float* floatNormal(void) const;
  void addNormal(const MathVector *norm);
  void resetNormal();

  float length(void) const;
  void set(const Vertex *v);
  void set(const MathVector *v);
  void move(const MathVector *v);
  void move(float x, float y, float z);
  void scale(float x, float y, float z);
  void scale(float scale);
  void rotate(const Matrix3<float> *r);
  int equal(const Vertex *v) const;

  void setTextCoordinates(float s, float t);
  void setTextCoordinates(const pair<float, float>* texCoord);
  float getTextT(void);
  float getTextS(void);

  list<Triangle*>* getTriangles(void);
  list<Edge*>* getEdges(void);

  static float distance(const Vertex *v1, const Vertex *v2);

  int operator<(const Vertex &v) const;

  unsigned int name, number;

private:
  MathVector *normal, data;
  float textS, textT;
   int nrTri, nrEdg;
  list<Edge*> edges;
  list<Triangle*> triangles;
};

#endif
