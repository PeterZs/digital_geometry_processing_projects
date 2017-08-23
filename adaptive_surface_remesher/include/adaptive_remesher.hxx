/*
 * adaptive_remesher.hxx
 *
 *  Created on: Dec 1, 2016
 *      Author: hooshi
 */

#ifndef INCLUDE_ADAPTIVE_REMESHER_HXX_
#define INCLUDE_ADAPTIVE_REMESHER_HXX_

#include <list>

#include "GetPot"

#include "pointer_mesh.hxx"


class Projector;

/************************************************************************\
                              Remesher
\************************************************************************/

class AdaptiveRemesher
{
public:
	enum AreaRemeshMode{CONSTANT, CURVATURE, REFINE00};

	Mesh *mesh() {return &_mesh;}
	const AreaRemeshMode& arearemesh_mode() const{return _arearemesh_mode;}
	AreaRemeshMode& arearemesh_mode(){return _arearemesh_mode;}
	const int& n_tverts() const{return _target_verts;}
	int& n_tverts(){return _target_verts;}
	Projector* proj() const {return _proj.get();}

	//--------------------------------------- constructor
	AdaptiveRemesher(Mesh& mesh_in, GetPot &getpot, AreaRemeshMode=CONSTANT);

	//--------------------------------------- quality increasing functions
	void init();
	//
	void remesh_flipping();
	void random_flipping();
	//
	bool remesh_areabased(Vertex*, Mesh &helper);
	void remesh_areabased();
	void remesh_areabased(const int n_step, const int n_area);
	//
	bool remesh_smooth_laplacian(Vertex*, Mesh &helper);
	void remesh_smooth_laplacian(const int n_step = 1);
	//
	bool is_obtuse(HalfEdge *he) const;
	bool is_obtuse(Edge *edge) const;
	void remesh_obtuse_angles();
	//
	bool remesh_splitedges();
	bool remesh_collapseedges();

	//--------------------------------------- internal guys
	bool should_flip(Edge *) const;
	bool can_flip(Edge *) const;
	//
	double vertex_areabased_mu(Vertex*) const;
	//
	bool tri_fidelity(Vertex *vertices[], double &Qtri, double &Qvert) const;
	double tri_quality(Face *) const;
	Vertex* split_edge(Edge *);
	//
	Vertex* collapse_edge(Edge *) ;

	//--------------------------------------- Mapping Stuff
	void geodesic_edge(Edge*  , Vector3 []) const;
	void geodesic_vert(Vertex*, Mesh& ) const;


private:
	Mesh &_mesh;
	AreaRemeshMode _arearemesh_mode;
	int _target_verts;
	double _ctheta_tri, _ctheta_vert;
	std::unique_ptr<Projector> _proj;
};

/************************************************************************\
                              Projector
\************************************************************************/

// Abstract
class Projector
{
protected:
	Mesh &_mesh;
	Projector(Mesh &mesh): _mesh(mesh){}
public:
   virtual void copy_vert_data(const Vertex* vfrom, Vertex* vto) const
    {
    	vto->xyz = vfrom->xyz;
    	vto->normal = vfrom->normal;
    	vto->facer = vfrom->facer;
    	vto->abcr = vfrom->abcr;
    }

	virtual void find_point(Face*, const Vector3& abc, Vertex *out) =0;
	virtual void init()=0;
	virtual ~Projector(){}
	static Projector* build(GetPot &getpot, Mesh&);

};

// Plane
class PlaneProjector: public Projector
{
public:
	void init();
	PlaneProjector(Mesh &mesh):Projector(mesh){}
	void find_point(Face* face, const Vector3& abc, Vertex *out);
};

// Sphere
class SphereProjector: public Projector
{
public:
	SphereProjector(Mesh &mesh, const Vector3 scale_in);

	Vector3 scale;
	Mat3 trans, i_trans, n_trans;
	void init();
	void find_point(Face*, const Vector3& abc, Vertex *out);
};

// Projection on current mesh!
class CurrentProjector: public Projector
{
public:
	CurrentProjector(Mesh&);

	void init();
	void find_point(Face*, const Vector3& abc, Vertex *out);
};

// General Case: Patch
class PatchProjector: public Projector
{
public:
	Mesh mesh_ref;
	std::list<Mesh*> patches;

	PatchProjector(Mesh& mesh_in);
	~PatchProjector();
	void init();
	void find_point(Face*, const Vector3& abc, Vertex *out);
	bool create_patch(Face*, Mesh *& patch, Face *pfaces[3]);

	int bfs_search(std::set<Face*>& members, Face *source, Face *targ1, Face *targ2, int = -1);
};

#endif /* INCLUDE_ADAPTIVE_REMESHER_HXX_ */
