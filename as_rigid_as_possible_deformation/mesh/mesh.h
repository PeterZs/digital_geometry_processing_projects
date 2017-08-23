//
//    File: mesh.h
//
//    (C) 2000-2008 Helmut Cantzler
//
//    Licensed under the terms of the Lesser General Public License.
//

#ifndef _MESH_H
#define _MESH_H

#include <list>
#include <vector>
#include <map>
#include <unordered_map>
#include <set>

#include <stdio.h>
#include <string.h>
#include <cassert>

#include "triangle.h"
#include "edge.h"
#include "vertex.h"

#ifdef WIN32
   #define strcasecmp stricmp
#endif

using namespace std;


class Mesh
{
public:
  Mesh();
  virtual ~Mesh();
  
  void clear(void);

  void setMesh(Mesh *mesh);
  void addTriangle(Triangle *t);
  void addEdge(Edge *e);
  void addVertex(Vertex *v);
  void addMesh(Mesh *mesh);
  void addTriangles(list<Triangle*> *newTriangles,
		    const char *textureName = "", int clearMap = TRUE);
  void addEdges(list<Edge*> *newEdges, int clearMap = TRUE);
  void addVertices(list<Vertex*> *newVertices);

  void remove(Vertex *v);
  void remove(Edge *e);
  void remove(Triangle *t);

  Vertex* getVertex(unsigned int name) const;
  Edge* getEdge(unsigned int name) const;
  Triangle* getTriangle(unsigned int name) const;

  list<Vertex*>* getVertices(void) const;
  list<Edge*>* getEdges(void) const;
  list<Edge*>* getEdges(list<Edge*> *es) const;
  list<Triangle*>* getTriangles(void) const;

  void setName(const char *n);
  void setPath(const char *p);
  char* getName() const;
  char* getPath() const;

  void writePoints(FILE *f) const;
  void writeGtsPoints(FILE *f) const;

  void createEdges(void);

  int numberOfVertices(void) const;
  int numberOfEdges(void) const;
  int numberOfTriangles(void) const;
  float getXMin(void) const;
  float getXMax(void) const;
  float getYMin(void) const;
  float getYMax(void) const;
  float getZMin(void) const;
  float getZMax(void) const;
  float averageTriangleSize(void) const;
  float getMaxVertexLength(void) const;
  MathVector getCentroid(void) const;
  void setMinMaxValues(const float* = NULL);

  void move(const MathVector *v);
  void moveToCentre(void);
  void scale(float scale);
  void scaleIntoNormalSphere(void);
  void calcOriginalCoordinates(const Vertex *v, Vertex *org) const;
  MathVector getModelCentroid(void) const;
  float getModelScale(void) const;
  void scaleAccordingToReferenceMesh(const Mesh *mesh);

  void negateSurfaceNormals(void);
  void removeDoublePoints(void);

  Vertex* findClosedPoint(const Vertex *v) const;

  // --------------------------- Selection -------------------

  void clearSelection(int i=0);
  void selectVertex(Vertex* v, unsigned int i);
  set<Vertex*>* getSelectedVertices(const unsigned short i);
  float* getSelectedCenter(const unsigned short i);

  int n_select() const {return NSELECT;}


protected:
  static constexpr unsigned short  NSELECT = 9;

  // helper function for createEdges
  typedef  map< pair<Vertex*,Vertex*>, Edge* > EdgeMap;

  Edge* getEdge(EdgeMap *edgeMap, Vertex *v1, Vertex *v2, int& res);

  // data
  list<Triangle*> *triangles;
  list<Edge*> *edges;
  list<Vertex*> *vertices;
  set<Vertex*> selectedVertices[NSELECT];
  float selectedCenter[NSELECT][3];
  int verNr, triNr, edgeNr;
  char *fileName, *path;

  // a map to avoid adding vertices more than once
  // in addTriangles() & addEdges()
  map<Vertex*, Vertex*> vertexMap;

  MathVector modelCentroid;
  float modelScale;

  float xMin, xMax, yMin, yMax, zMin, zMax;
};

#endif
