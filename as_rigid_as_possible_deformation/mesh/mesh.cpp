//
//    File: mesh.cc
//
//    (C) 2000-2008 Helmut Cantzler
//
//    Licensed under the terms of the Lesser General Public License.
//

#include "mesh.h"

#ifndef WIN32
    #include <values.h>
    #include <limits.h>
#else
    #include <float.h>
    #define MAXFLOAT FLT_MAX
#endif

Mesh::Mesh()
{
  edgeNr = verNr = triNr = 0;
  xMin = xMax = yMin = yMax = zMin = zMax = 0;

  modelCentroid.setZero();
  modelScale = 1.0;

  triangles = new list<Triangle*>;
  vertices = new list<Vertex*>;
  edges = new list<Edge*>;

  path = NULL;
  fileName = NULL;
}

Mesh::~Mesh()
{
  clear();

  delete triangles;
  delete edges;
  delete vertices;

  if(path) delete[] path;
  if(fileName) delete[] fileName;
}

void Mesh::clear(void)
{
  list<Triangle*>::iterator it;
  list<Edge*>::iterator ie;
  list<Vertex*>::iterator iv;

  for (it=triangles->begin(); it != triangles->end(); it++)
    delete *it;
  for (ie=edges->begin(); ie != edges->end(); ie++)
    delete *ie;
  for (iv=vertices->begin(); iv != vertices->end(); iv++)
    delete *iv;

  triangles->clear();
  edges->clear();
  vertices->clear();

  edgeNr = verNr = triNr = 0;
  xMin = xMax = yMin = yMax = zMin = zMax = 0;
}

int Mesh::numberOfVertices(void) const
{
  return verNr;
}

int Mesh::numberOfEdges(void) const
{
  return edgeNr;
}

int Mesh::numberOfTriangles(void) const
{
  return triNr;
}

float Mesh::getXMin(void) const
{
  return xMin;
}

float Mesh::getXMax(void) const
{
  return xMax;
}

float Mesh::getYMin(void) const
{
  return yMin;
}

float Mesh::getYMax(void) const
{
  return yMax;
}

float Mesh::getZMin(void) const
{
  return zMin;
}

float Mesh::getZMax(void) const
{
  return zMax;
}

void Mesh::addTriangle(Triangle *t)
{
  triangles->push_back(t);
  t->name = triNr++;
}

void Mesh::addEdge(Edge *e)
{
  edges->push_back(e);
  e->name = edgeNr++;
}

void Mesh::addVertex(Vertex *v)
{
  vertices->push_back(v);
  v->name = verNr++;
}

void Mesh::remove(Vertex *v)
{

  vertices->remove(v);
  delete v;
}

void Mesh::remove(Edge *e)
{

  edges->remove(e);

  delete e;
}

void Mesh::remove(Triangle *t)
{

  triangles->remove(t);

  delete t;
}

void Mesh::addVertices(list<Vertex*> *newVertices)
{
  list<Vertex*>::iterator iv;
  Vertex *ver;


  for (iv=newVertices->begin(); iv != newVertices->end(); iv++)
    {
      ver = new Vertex( (*iv)->x(), (*iv)->y(), (*iv)->z() );
      addVertex(ver);
    }
}

void Mesh::addEdges(list<Edge*> *newEdges, int clearMap)
{
  // mapping from old to new vertices
  map<Vertex*, Vertex*>::iterator iv;
  list<Edge*>::iterator ie;
  Vertex *ver, *v[2];
  Edge *edg;
  int i;

  for (ie=newEdges->begin(); ie != newEdges->end(); ie++)
    {
      // Get or create the two new vertices
      for (i=0; i < 2; i++)
	{
	  ver=(*ie)->vertices[i];
	  // avoid using vertices more than once
	  iv=vertexMap.find(ver);
	  if (iv == vertexMap.end())
	    {
	      // Create one new vertex
	      v[i] = new Vertex( ver->x(), ver->y(), ver->z() );
	      addVertex(v[i]);
	      // old vertex => new vertex
	      vertexMap[ver] = v[i];
	    }
	  else
	    // We found the vertex in our list
	    v[i]=(*iv).second;
	}
	    
      // Create one new edge with the two vertices
      edg = new Edge(v[0], v[1]);
      addEdge(edg);
    }


  if (clearMap)
    vertexMap.clear();
}

void Mesh::addTriangles(list<Triangle*> *newTriangles, 
			 const char *textureName, int clearMap)
{
  // mapping from old to new vertices
  map<Vertex*, Vertex*>::iterator iv;
  list<Triangle*>::iterator it;
  Vertex *ver, *v[3];
  Triangle *tri;
  int i;

  for (it=newTriangles->begin(); it != newTriangles->end(); it++)
    {
      // Get or create the three new vertices
      for (i=0; i < 3; i++)
	{
	  ver=(*it)->vertices[i];
	  // avoid using vertices more than once
	  iv=vertexMap.find(ver);
	  if (iv == vertexMap.end())
	    {
	      // Create one new vertex
	      v[i] = new Vertex( ver->x(), ver->y(), ver->z() );
	      addVertex(v[i]);
	      // old vertex => new vertex
	      vertexMap[ver] = v[i];
	    }
	  else
	    // We found the vertex in our list
	    v[i]=(*iv).second;
	}
      
      // Create one new triangle with the three vertices
      tri = new Triangle(v[0], v[1], v[2]);
      for (i=0; i < 3; i++)
	tri->setTextCoordinates( i, (*it)->getTextS(i),
				 (*it)->getTextT(i) );
      addTriangle(tri);
    }


  if (clearMap)
    vertexMap.clear();
}

void Mesh::addMesh(Mesh *mesh)
{
	if (mesh != NULL)
	{
		if (mesh->numberOfTriangles() > 0)
			addTriangles(mesh->triangles);
		else if (mesh->numberOfEdges() > 0)
			addEdges(mesh->getEdges());
		else if (mesh->numberOfVertices() > 0)
			addVertices(mesh->getVertices());
	}
}

void Mesh::setMesh(Mesh *mesh)
{
  modelCentroid = mesh->modelCentroid;
  modelScale = mesh->modelScale;

  clear();

  addMesh(mesh);
  setMinMaxValues();
}

Vertex* Mesh::getVertex(unsigned int name) const
{
  list<Vertex*>::iterator iv;

  for (iv=vertices->begin(); iv != vertices->end(); iv++)
    if ((*iv)->name == name)
      return *iv;

  return NULL;
}

Edge* Mesh::getEdge(unsigned int name) const
{
  list<Edge*>::iterator ie;

  for (ie=edges->begin(); ie != edges->end(); ie++)
    if ((*ie)->name == name)
      return *ie;

  return NULL;
}

Triangle* Mesh::getTriangle(unsigned int name) const
{
  list<Triangle*>::iterator it;

  for (it=triangles->begin(); it != triangles->end(); it++)
    if ((*it)->name == name)
      return *it;

  return NULL;
}

void Mesh::setName(const char *n)
{
  if(fileName) delete[] fileName;
  fileName = new char[strlen(n)+1];
  strcpy(fileName, n);
}

void Mesh::setPath(const char *p)
{
  if(path) delete[] path;
  path = new char[strlen(p)+1];
  strcpy(path, p);
}

char* Mesh::getName() const
{
  return fileName;
}

char* Mesh::getPath() const
{
  return path;
}

list<Triangle*>* Mesh::getTriangles(void) const
{
  return triangles;
}

list<Vertex*>* Mesh::getVertices(void) const
{
  return vertices;
}

list<Edge*>* Mesh::getEdges(void) const
{
  return edges;
}


list<Edge*>* Mesh::getEdges(list<Edge*> *es) const
{
  list<Edge*> *edgeList = new list<Edge*>;
  list<Edge*>::iterator ie1, ie2;

  for (ie1=es->begin(); ie1 != es->end(); ie1++)
    for (ie2=edges->begin(); ie2 != edges->end(); ie2++)
      if ((*ie2)->equal(*ie1))
	edgeList->push_back(*ie2);

  return edgeList;
}

Edge* Mesh::getEdge(EdgeMap *edgeMap, Vertex *v1, Vertex *v2, int& res)
{
  EdgeMap::iterator ei;

  // search in map (hash table)
  // first try to get the edge
  ei=edgeMap->find( pair<Vertex*,Vertex*> (v1, v2) );
  res = 0;

  // second try to get the edge (reverse order)
  if (ei == edgeMap->end())
  {
    ei=edgeMap->find( pair<Vertex*,Vertex*> (v2, v1) );
    res = 1;
  }

  if (ei == edgeMap->end())
  {
	  // create new edge
	  Edge *e = new Edge(v1, v2);
	  addEdge(e);

	  // new hash entry
	  (*edgeMap)[ pair<Vertex*,Vertex*> (v1, v2) ] = e;
	  res = 2;
	  return e;
  }
  else  return (*ei).second;
}

void Mesh::createEdges(void)
{
  list<Triangle*>::iterator it;
  EdgeMap edgeMap;
  Edge *e;

  for (it=triangles->begin(); it != triangles->end(); it++)
  {
	  Triangle *tri = *it;
	  int res;

	  for (int i = 0 ; i < 3 ; i++)
	  {
		  const int j = (i + 1) % 3;

		  e = getEdge(&edgeMap, tri->vertices[i], tri->vertices[j], res);
		  if (res == 0)
		  {
			  printf("Triangle %d, verts %d %d also belong to \n", tri->name, tri->vertices[i]->name, tri->vertices[j]->name);
			  printf("Triangle %d or %d.  \n", e->triangles[0]->name, e->triangles[1]->name);
			  assert(0 && "Non-manifold mesh.");
			  throw;
		  }
		  else if (res == 1)
		  {
			  assert(e->triangles[0] != Triangle::getSentinel());
			  assert(e->triangles[1] == Triangle::getSentinel());
			  e->triangles[1] = tri;
		  }
		  else if (res == 2)
		  {
			  assert(e->triangles[0] == Triangle::getSentinel());
			  assert(e->triangles[1] == Triangle::getSentinel());
			  e->triangles[0] = tri;
		  }
		  tri->edges[i] = e;
	  }

    }
}

float Mesh::averageTriangleSize(void) const
{
  list<Triangle*>::iterator it;
  float size=0.0;

  for (it=triangles->begin(); it != triangles->end(); it++)
    size+=(*it)->perimeter();

  return size/triNr;
}

void Mesh::setMinMaxValues(const float *override)
{
	if(override)
	{
		xMin = override[0];
		yMin = override[1];
		zMin = override[2];
		xMax = override[3];
		yMax = override[4];
		zMax = override[5];
		return;
	}

  list<Vertex*>::iterator iv;

  xMin = xMax = yMin = yMax = zMin = zMax = 0;

  // find smalles & largest vertex in all directions (x, y, z)
  for (iv=vertices->begin(); iv != vertices->end(); iv++)
    {
      if (xMin > (*iv)->x())
	xMin = (*iv)->x();
      if (xMax < (*iv)->x())
	xMax = (*iv)->x();
      if (yMin > (*iv)->y())
	yMin = (*iv)->y();
      if (yMax < (*iv)->y())
	yMax = (*iv)->y();
      if (zMin > (*iv)->z())
	zMin = (*iv)->z();
      if (zMax < (*iv)->z())
	zMax = (*iv)->z();
    }
}

void Mesh::moveToCentre(void)
{
  modelCentroid = getCentroid();
  modelCentroid.negation();
  move(&modelCentroid);
  modelCentroid.negation();
  //printf("Centriod is at %f %f %f \n", modelCentroid.v[0], modelCentroid.v[1], modelCentroid.v[2]);
}

MathVector Mesh::getCentroid(void) const
{
  list<Vertex*>::iterator iv;
  MathVector c;

  // find centroid
  for (iv=vertices->begin(); iv != vertices->end(); iv++)
  {
    c.add((*iv)->mathData());
  }
  c.scale(1.0/verNr);

  return c;
}

void Mesh::move(const MathVector *v)
{
  list<Triangle*>::iterator it;
  list<Vertex*>::iterator iv;

  // Move vertices and triangle centroids
  for (iv=vertices->begin(); iv != vertices->end(); iv++)
    (*iv)->move(v);
  for (it=triangles->begin(); it != triangles->end(); it++)
    (*it)->moveCentroid(v);
}

void Mesh::scaleIntoNormalSphere(void)
{
  // Scaling mesh into +/- 1 sphere
  modelScale = sqrt(getMaxVertexLength())/2;
  scale( modelScale );
  modelScale = 1.0 / modelScale;
  setMinMaxValues();
}

float Mesh::getMaxVertexLength(void) const
{
  list<Vertex*>::iterator iv;
  float length=0.0;

  // Finding the length
  for (iv=vertices->begin(); iv != vertices->end(); iv++)
    length = MAX( float, (*iv)->length(), length );

  return length;
}

void Mesh::scale(float scale)
{
  list<Triangle*>::iterator it;
  list<Edge*>::iterator ie;
  list<Vertex*>::iterator iv;

  // Scaling vertices
  for (iv=vertices->begin(); iv != vertices->end(); iv++)
    (*iv)->scale(scale);
  // Recalculate triangle centroid, size, perimeter and so on
  for (it=triangles->begin(); it != triangles->end(); it++)
    (*it)->calcProperties();
  // Recalculate edge centroid, length and so on
  for (ie=edges->begin(); ie != edges->end(); ie++)
    (*ie)->calcProperties();
}

void Mesh::calcOriginalCoordinates(const Vertex *v, Vertex *org) const
{
  org->set(v);
  org->scale(modelScale);
  org->move(&modelCentroid);
}

MathVector Mesh::getModelCentroid(void) const
{
  return modelCentroid;
}

float Mesh::getModelScale(void) const
{
  return modelScale;
}

void Mesh::scaleAccordingToReferenceMesh(const Mesh *mesh)
{
  modelCentroid = mesh->getModelCentroid();
  modelScale = 1.0/mesh->getModelScale();
  modelCentroid.negation();

  move(&modelCentroid);
  scale(modelScale);

  modelScale = 1.0/modelScale;
  modelCentroid.negation();
  setMinMaxValues();
}

void Mesh::clearSelection(int i)
{
	assert( (i >= 0) && (i < NSELECT) );
	selectedVertices[i].clear();
	selectedCenter[i][0] = 0.;
	selectedCenter[i][1] = 0.;
	selectedCenter[i][2] = 0.;
}

void Mesh::selectVertex(Vertex* v, unsigned int i)
{
	  if(!((i>=0) && (i<NSELECT)))
	  {
		  printf("ONLY %d SELECTION SETS AVAILABLE. \n", NSELECT);
		  assert(0); throw;
	  }
	  if (v != NULL)
	  {
		  bool prev_selected =false;

		  for (unsigned short j = 0 ; j < i ; j++)
		  {
			  if(selectedVertices[j].count(v))
			  {
				  prev_selected = true;
				  break;
			  }
		  }

		  if(!prev_selected)
		  {
			  for (int dim=0 ; dim < 3 ; dim++)
				  selectedCenter[i][dim] =
						  float(selectedVertices[i].size() * selectedCenter[i][dim] + v->floatData()[dim]) /
						  float(selectedVertices[i].size() + 1);
			  selectedVertices[i].insert(v);
		  }
	  }
}

set<Vertex*>* Mesh::getSelectedVertices(const unsigned short i)
{
	  if(!((i>=0) && (i<NSELECT)))
	  {
		  printf("ONLY %d SELECTION SETS AVAILABLE. \n", NSELECT);
		  assert(0); throw;
	  }
	  return &selectedVertices[i];
}

float* Mesh::getSelectedCenter(const unsigned short i)
{
	  if(!((i>=0) && (i<NSELECT)))
	  {
		  printf("ONLY %d SELECTION SETS AVAILABLE. \n", NSELECT);
		  assert(0); throw;
	  }
	  return selectedCenter[i];

}


void Mesh::negateSurfaceNormals(void)
{
  list<Triangle*>::iterator it;

  // negate all triangle surface normals
  for (it=triangles->begin(); it != triangles->end(); it++)
    (*it)->negateNormal();
}

void Mesh::removeDoublePoints(void)
{
  // the MathVector represents the x, y and z coordinates of the vertex
  map<MathVector, Vertex*> vertexMap;
  map<MathVector, Vertex*>::iterator iMap;
  list<Vertex*>::iterator iList;
  list<Triangle*>::iterator it;
  list<Triangle*> *verTriangles;
  list<Edge*>::iterator ie;
  list<Edge*> *verEdges;

  for (iList=vertices->begin(); iList != vertices->end();)
    {
      // Check if we have the vertex already
      iMap=vertexMap.find( *(*iList)->mathData() );
      if (iMap == vertexMap.end())
	{
	  // Create a new entry
	  vertexMap[ *(*iList)->mathData() ] = *iList;
	  iList++;
	}
      else
	{
	  // change all references to the vertex
	  verTriangles=(*iList)->getTriangles();
	  for (it=verTriangles->begin(); it != verTriangles->end(); it++)
	    (*it)->changeVertex(*iList, (*iMap).second);
	  verEdges=(*iList)->getEdges();
	  for (ie=verEdges->begin(); ie != verEdges->end(); ie++)
	    (*ie)->changeVertex(*iList, (*iMap).second);


	  // delete vertex and remove from list
	  delete *iList;
	  iList=vertices->erase(iList);
	}
    }

  verNr= (int) vertices->size();
}

Vertex* Mesh::findClosedPoint(const Vertex *v) const
{
  list<Vertex*>::iterator iv;
  MathVector con;
  Vertex *best = NULL;
  float bestLength;

  bestLength=MAXFLOAT;
  for (iv=vertices->begin(); iv != vertices->end(); iv++)
    {
      MathVector::sub((*iv)->mathData(), v->mathData(), &con);
      if (con.length() < bestLength)
	{
	  bestLength=con.length();
	  best=*iv;
	}
    }
 
  return best;
}

void Mesh::writeGtsPoints(FILE *f) const
{
  int n;
  list<Vertex*>::iterator iv;

  // header
  fprintf(f,"%d 0 0\n", verNr);

  // vertex::number == vertex::name ???

  n=0;   // vertices
  for (iv=vertices->begin(); iv != vertices->end(); iv++)
    {
      fprintf(f,"%f %f %f\n", (*iv)->x(), (*iv)->y(), (*iv)->z());
      n++;
      (*iv)->number=n;
    }
}

void Mesh::writePoints(FILE *f) const
{
  list<Vertex*>::iterator iv;

  for (iv=vertices->begin(); iv != vertices->end(); iv++)
    fprintf(f,"%f %f %f\n", (*iv)->x(), (*iv)->y(), (*iv)->z());
}
