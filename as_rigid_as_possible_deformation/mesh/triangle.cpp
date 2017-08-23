//
//    File: triangle.cc
//
//    (C) 2000-2008 Helmut Cantzler
//
//    Licensed under the terms of the Lesser General Public License.
//

#include "triangle.h"

#include <cassert>

Triangle Triangle::sentinel(NULL, NULL, NULL, SENTINEL_TRIANGLE);

Triangle::Triangle(Vertex *v1, Vertex *v2, Vertex *v3, TriangleType type_in)
{
	int v;

	type = type_in;
	if(type == SENTINEL_TRIANGLE)
	{
		fromPolygon=false;
		plane=0;
		vertices[0]=vertices[1]=vertices[2]=NULL;
		edges[0]=edges[1]=edges[2]=NULL;
		name = -1;
	}
	else
	{
		plane=0;
		fromPolygon = type == POLYGON_TRIANGLE ? true : false;

		vertices[0]=v1;
		vertices[1]=v2;
		vertices[2]=v3;

		edges[0]=edges[1]=edges[2]=NULL;

		calcProperties();

		for (v=0; v < 3; v++)
		{
			vertices[v]->addTriangle(this);
			vertices[v]->addNormal(&surfaceNormal);
		}
	}
}

//Triangle::Triangle(Edge *e1, Edge *e2, Edge *e3)
//{
//  int i;
//
//  plane=0;
//  fromPolygon = false;
//
//  edges[0]=e1;
//  edges[1]=e2;
//  edges[2]=e3;
//
//  for (i=0; i < 3; i++)
//    edges[i]->addTriangle(this);
//
//  if ( (edges[0]->vertices[0] != edges[1]->vertices[0]) &&
//       (edges[0]->vertices[0] != edges[1]->vertices[1]) )
//  {
//      vertices[0]=edges[0]->vertices[0];
//      vertices[1]=edges[0]->vertices[1];
//  }
//  else
//  {
//      vertices[0]=edges[0]->vertices[1];
//      vertices[1]=edges[0]->vertices[0];
//  }
//
//  if ( (edges[1]->vertices[0] != edges[0]->vertices[0]) &&
//       (edges[1]->vertices[0] != edges[0]->vertices[1]) )
//  {
//      vertices[2]=edges[1]->vertices[0];
//  }
//  else
//  {
//      vertices[2]=edges[1]->vertices[1];
//  }
//
//  calcProperties();
//
//  for (i=0; i < 3; i++)
//    {
//      vertices[i]->addTriangle(this);
//      vertices[i]->addNormal(&surfaceNormal);
//    }
//}

void Triangle::setVertices(const MathVector *v1, const MathVector *v2,
			   const MathVector *v3)
{
  vertices[0]->set(v1);
  vertices[1]->set(v2);
  vertices[2]->set(v3);

  calcProperties();
}

void Triangle::calcProperties(void)
{
  // calculate surface centroid and register triangle at the vertices
  surfaceCentroid.setZero();
  for (int v=0; v < 3; v++)
    surfaceCentroid.add(vertices[v]->mathData());
  surfaceCentroid.scale(1.0/3.0);

  // calculate the vectors between the vertices
  MathVector c[3];

  MathVector::sub(vertices[0]->mathData(), vertices[1]->mathData(), &c[0]);
  MathVector::sub(vertices[0]->mathData(), vertices[2]->mathData(), &c[1]);
  MathVector::sub(vertices[2]->mathData(), vertices[1]->mathData(), &c[2]);

  // calculate surface normal
  MathVector::crossProduct(&c[0], &c[1], &surfaceNormal);
  surfaceNormal.normalize();

  // calculate distance to the origin from the centroid in direction
  // of the surface normal
  surfaceDistanceToOrigin = MathVector::dotProduct(&surfaceCentroid,
						   &surfaceNormal);
  // calculate surface perimeter
  surfacePerimeter = c[0].length() + c[1].length() + c[2].length();

  // calculate surface size
  float s=surfacePerimeter/2.0;
  surfaceSize=sqrt(s * (s-c[0].length()) * (s-c[1].length()) * (s-c[2].length()));
}

Triangle::~Triangle()
{
	if(type != SENTINEL_TRIANGLE)
		for (int i=0; i < 3; i++)
		{
			vertices[i]->deleteTriangle(this);
			// if (edges[i] != NULL) 	edges[i]->deleteTriangle(this);
		}
}

void Triangle::moveCentroid(const MathVector *v)
{
  surfaceCentroid.add(v);
}

void Triangle::negateNormal(void)
{
  Vertex *dummy;

  dummy=vertices[1];
  vertices[1]=vertices[2];
  vertices[2]=dummy;

  surfaceNormal.negation();
  surfaceDistanceToOrigin = -surfaceDistanceToOrigin;
}

const float* Triangle::centroid(void) const
{
  return surfaceCentroid.v;
}

const MathVector* Triangle::mathCentroid(void) const
{
  return &surfaceCentroid;
}

const float* Triangle::floatNormal(void) const
{
  return surfaceNormal.v;
}

const MathVector* Triangle::mathNormal(void) const
{
  return &surfaceNormal;
}

float Triangle::size(void) const
{
  return surfaceSize;
}

float Triangle::perimeter(void) const
{
  return surfacePerimeter;
}

float Triangle::distanceToOrigin(void) const
{
  return surfaceDistanceToOrigin;
}

float Triangle::angle(const Vertex *v) const
{
  MathVector c[3];
  int i, k;

  // find the two vectors between the vertices
  for (i=k=0; i < 3; i++)
    if (v != vertices[i])
      MathVector::sub(v->mathData(), vertices[i]->mathData(), &c[k++]);

  // calculate the angle between the edges
  return k == 2 ? MathVector::angle90(&c[0], &c[1]) : 0;
}

vector<Triangle*> Triangle::neighbors(void) const
{
  vector<Triangle*> nei;

  for (int i=0; i < 3; i++)
  {
	  if( edges[i]->makeFirstTriangle(this) == NULL )
	  {
		  printf("Fatal error. \n");
		  assert(0); throw;
	  }
	  nei.push_back(edges[i]->triangles[1]);
 }

  return nei;
}

vector<Triangle*> Triangle::neighborsOnPlane(void) const
{
  vector<Triangle*> nei = neighbors();
  vector<Triangle*>::iterator it;

  for (it=nei.begin(); it != nei.end(); )
    if ((*it)->plane == plane)
      it++;
    else
      it=nei.erase(it);
  
  return nei;
}

int Triangle::getSurroundingPlane(void) const
{
  vector<Triangle*> nei = neighbors();
  int i, k, size = (int) nei.size();

  if (size == 1)
    return nei[0]->plane;

  // look for two which are the same
  for (i=0; i < size-1; i++)
    for (k=i+1; k < size; k++)
      if (nei[i]->plane == nei[k]->plane)
	return nei[i]->plane;

  return 0;
}

int Triangle::changeVertex(const Vertex *oldV, Vertex *newV)
{
  int i=0;

  while (i < 3 && vertices[i] != oldV)
    i++;

  if (i != 3)
    {
      vertices[i]=newV;
      return 1;
    }

  return 0;
}

void Triangle::setTextCoordinates(int i, float s, float t)
{
  vertices[i]->setTextCoordinates(s, t); 
}

void Triangle::setTextCoordinates(int i, const pair<float, float>* texCoord)
{
  vertices[i]->setTextCoordinates(texCoord); 
}

float Triangle::getTextT(int i) const
{
  return vertices[i]->getTextT();
}

float Triangle::getTextS(int i) const
{
  return vertices[i]->getTextS();
}

bool Triangle::equal(const Triangle *tri) const
{
  return vertices[0] == tri->vertices[0] &&
    vertices[1] == tri->vertices[1] && vertices[2] == tri->vertices[2];
}

bool Triangle::onSamePlane(const Triangle *tri, float distanceTolerance,
			  float orientationTolerance) const
{
  MathVector normal2 = tri->surfaceNormal;
  float distance2 = MathVector::dotProduct(&tri->surfaceCentroid,
					   &surfaceNormal);

  // triangle facing in different directions?
  if (MathVector::angle90(&surfaceNormal, &normal2) > 
      GRAD2RAD(orientationTolerance))
    return 0;

  // optimisation ?!?!
  //
  //  if (MathVector::angle90(&surfaceNormal, &normal2) > GRAD2RAD(20.0))
  //
  // because the normals are already normalised
  //  if (acos(MathVector::dotProduct(&surfaceNormal, &normal2)) > GRAD2RAD(20.0))
  //
  // one function call less
  //  if (acos(MathVector::dotProduct(&surfaceNormal, &normal2)) > GRAD2RAD(20.0))
  //
  // move cos on other side
  //  if (MathVector::dotProduct(&surfaceNormal, &normal2) < cos(GRAD2RAD(20.0)))

  float d = surfaceDistanceToOrigin - distance2;

  // the centroids are on the same plane?
  return fabs(d)*100 < distanceTolerance;
}

bool Triangle::isFromPolygon(void) const
{
  return fromPolygon;
}

Vertex* Triangle::nextVertex(Vertex *v)
{
	if (v == vertices[0]) return vertices[1];
	if (v == vertices[1]) return vertices[2];
	if (v == vertices[2]) return vertices[0];
	else return NULL;
}
Vertex* Triangle::prevVertex(Vertex *v)
{
	if (v == vertices[0]) return vertices[2];
	if (v == vertices[1]) return vertices[0];
	if (v == vertices[2]) return vertices[1];
	else return NULL;
}
