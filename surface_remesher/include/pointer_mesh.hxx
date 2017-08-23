/*
 * pointer_mesh.hxx
 *
 *  Created on: Nov 29, 2016
 *      Author: hooshi
 *
 *  A duplicate of the cartel mesh that uses pointers instead of
 *  indices.
 */

#ifndef INCLUDE_POINTER_MESH_HXX_
#define INCLUDE_POINTER_MESH_HXX_


#include <cmath>
#include <map>
#include <set>
#include <iostream>
#include <memory>
#include <fstream>
#include <list>

#include <armadillo>

#include "indexed_list.hxx"


/****************************************************************************\
                                Typedefs
\****************************************************************************/
typedef arma::Col<double>::fixed<3> Vector3;
typedef arma::Mat<double>::fixed<3,3> Mat3;
struct Edge;
struct Vertex;
struct Face;
class Mesh;

/****************************************************************************\
                                Half Edge
\****************************************************************************/

struct HalfEdge: public IndexedItem<HalfEdge>
{
    HalfEdge *next;
    HalfEdge *twin;
    Edge *edge;
    Vertex *vert;
    Face *face;

    HalfEdge(HalfEdge* next_in=NULL,
    		HalfEdge* twin_in=NULL,
			Edge* edge_in=NULL,
			Vertex* vertex_in=NULL,
			Face* face_in=NULL):
    	next(next_in),
		twin(twin_in),
		edge(edge_in),
		vert(vertex_in),
		face(face_in)
    {}
};

/****************************************************************************\
                               Vertex
\****************************************************************************/

struct Vertex: public IndexedItem<Vertex>
{
    HalfEdge *he;
    Vector3 xyz;
    Vector3 normal;

    // Fidelity and projection
    Vector3 abcr;
    Face *facer;

};

struct BdryVertex: public IndexedItem<BdryVertex>
{
	Vertex *vertex;
};

/****************************************************************************\
                                Edge
\****************************************************************************/

struct Edge: public IndexedItem<Edge>
{
	HalfEdge *he;
};

struct BdryEdge: public IndexedItem<BdryEdge>
{
	Edge *edge;
};

/****************************************************************************\
                                Face
\****************************************************************************/
struct Face: public IndexedItem<Face>
{
	HalfEdge *he;
    bool xy2abc(const Vector3 &xy, Vector3& abc) const;
    void abc2xy(const Vector3 &abc, Vector3& xyz) const;

    // struct pinfo_cmp
    // {
    //	bool operator() ( const std::pair<Mesh*, Face*>& p1, const std::pair<Mesh*, Face*>& p2) const
    //	{return std::less<Mesh*>()(p1.first, p2.first);}
    // };
    //std::set< std::pair<Mesh*, Face*>, pinfo_cmp> pinfo;
    std::list< std::pair<Mesh*, Face*> > pinfo;
};

/****************************************************************************\
                                 Mesh
\****************************************************************************/
class Mesh
{
    /****************************************************************************\
                                Constructor and Destructor
    \****************************************************************************/

public:

    Mesh();
    ~Mesh();
    void clear_memory();
    void init( const std::vector<double>& xyz, const std::vector<int>& tris);
    void verify() const;

    /****************************************************************************\
                                      Iterators
    \****************************************************************************/

public:

    // Half edges around a vertex pointing towards it
    class ring_iterator
    {
    private:
        HalfEdge* _cur;
        HalfEdge* _end;
    public:
        ring_iterator():_cur(NULL), _end(NULL){}
        bool init(Vertex*);
        bool reset_boundary();
        bool advance();
        HalfEdge* half_edge() {return _cur;}
    };

    /****************************************************************************\
                                Mesh Modifications
    \****************************************************************************/

    HalfEdge* find_connecting_hedge( Vertex* vfrom, Vertex* vto );
    int count_connecting_hedges(Vertex* vfrom, Vertex* vto );
    bool collapse_edge(Edge*, Vertex*);
    bool split_edge(Edge*, Vertex*);
    bool flip_edge(Edge*);

    bool locate_point(const Vector3 &xy, Face *&face, Vector3& abc) const;

    /****************************************************************************\
            Constructing Boundary Info ( for subdivision or dual creation)
    \****************************************************************************/
    void get_boundary_info(IndexedList<BdryVertex>&, IndexedList<BdryEdge>&, const int two_way_map_index);

    /****************************************************************************\
                                     Queries
    \****************************************************************************/
    static Face* hole() {return &_hole;}

    /****************************************************************************\
                                     Members
    \****************************************************************************/

private:
    static Face _hole;

public:

    IndexedList<Vertex> verts;
    IndexedList<Edge> edges;
    IndexedList<Face> faces;
    IndexedList<HalfEdge> half_edges;

};


#endif /* INCLUDE_POINTER_MESH_HXX_ */
