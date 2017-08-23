//
//    File: edge.cc
//
//    (C) 2000-2008 Helmut Cantzler
//
//    Licensed under the terms of the Lesser General Public License.
//

#include "triangle.h"
#include "edge.h"

Edge::Edge(Vertex *v1, Vertex *v2)
{
  vertices[0]=v1;   v1->addEdge(this);
  vertices[1]=v2;   v2->addEdge(this);

  nrTri=0;

  triangles[0] = Triangle::getSentinel();
  triangles[1] = Triangle::getSentinel();

  calcProperties();
}

Edge::~Edge()
{
  for (int i=0; i < 2; i++)
    vertices[i]->deleteEdge(this);
}

void Edge::set(const MathVector *v1, const MathVector *v2)
{
  vertices[0]->set(v1);
  vertices[1]->set(v2);

  calcProperties();
}

void Edge::calcProperties(void)
{
  edgeCentroid.setZero();
  for (int v=0; v < 2; v++)
    edgeCentroid.add(vertices[v]->mathData());
  edgeCentroid.scale(1.0/2.0);

  // connection between the vertices
  MathVector::sub(vertices[0]->mathData(), vertices[1]->mathData(), 
		  &edgeOrientation);
  // length of the connection
  edgeLength = edgeOrientation.length();
  // connection => orientation
  edgeOrientation.normalize();
}

int Edge::equal(const Vertex *v1, const Vertex *v2) const
{
  if ((vertices[0] == v1 && vertices[1] == v2) ||
      (vertices[0] == v2 && vertices[1] == v1))
    return TRUE;
  else
    return FALSE;
}

int Edge::equal(const Edge *e) const
{
  if (vertices[0]->equal(e->vertices[0]) && vertices[1]->equal(e->vertices[1]))
    return TRUE;
  if (vertices[0]->equal(e->vertices[1]) && vertices[1]->equal(e->vertices[0]))
    return TRUE;

  return FALSE;
}

int Edge::changeVertex(const Vertex *oldV, Vertex *newV)
{
  int i=0;

  while (i < 2 && vertices[i] != oldV)
    i++;

  if (i != 2)
    {
      vertices[i]=newV;
      return 1;
    }

  return 0;
}

int Edge::onSameEdge(const Edge *e, float tolerance) const
{
  // Check if both points from e are on the edge

  MathVector point, ori;
  float distance[2];

  //  if (MathVector::angle2(&edgeOrientation, e->mathOrientation())
  //      >= GRAD2RAD(70.0))
  //    return 0;

  for (int i=0; i < 2; i++)
    {
      ori=edgeOrientation;
      point.set(e->vertices[i]->mathData());
      point.sub(vertices[0]->mathData());
      ori.scale(MathVector::dotProduct(&point, &ori));
      point.sub(&ori);
      distance[i]=point.length();

    }  

  return distance[0] <= tolerance && distance[1] <= tolerance;
}

//void Edge::(Triangle* triangle)
//{
//  triangles.push_back(triangle);
//  nrTri++;
//}
//
//void Edge::deleteTriangle(Triangle* triangle)
//{
//  triangles.remove(triangle);
//  nrTri--;
//}

const MathVector* Edge::mathOrientation(void) const
{
  return &edgeOrientation;
}

float Edge::length(void) const
{
  return edgeLength;
}

const float* Edge::centroid(void) const
{
  return edgeCentroid.v;
}

const MathVector* Edge::mathCentroid(void) const
{
  return &edgeCentroid;
}

//static void swap(void *p1, void *p2)
//{
//	void *tmp;
//	tmp = p1;
//	p1 = p2;
//	p2 = tmp;
//}
Edge* Edge::makeFirstTriangle(const Triangle* tri)
{
	if(tri == triangles[0]) return this;
	else if(tri == triangles[1])
	{
		swap(triangles[0], triangles[1]);
		swap(vertices[0], vertices[1]);
		edgeOrientation.negation();
		return this;
	}
	else return NULL;
}
