/*
 * pointer_mesh.cxx
 *
 *  Created on: Nov 30, 2016
 *      Author: hooshi
 */

#include <algorithm>

#include "pointer_mesh.hxx"

/****************************************************************************\
                                     Face
\****************************************************************************/

// Only planar intersection
bool Face::xy2abc(const Vector3 &xy, Vector3& abc) const
{
	Vertex *v[3] = {he->vert, he->next->vert, he->next->next->vert};
	const Vector3 e1 = v[1]->xyz - v[0]->xyz;
	const Vector3 e2 = v[2]->xyz - v[0]->xyz;
	const Vector3 e3 = v[0]->xyz - xy;
	const Vector3 e4 = v[1]->xyz - xy;
	const Vector3 e5 = v[2]->xyz - xy;

	const Vector3 Atotv = arma::cross(e1, e2);
	const Vector3 n_tri = arma::normalise(Atotv);
	const double Atot = arma::norm(Atotv);
	const Vector3 A0v = arma::cross(e4, e5);
	const Vector3 A1v = arma::cross(e5, e3);
	const Vector3 A2v = arma::cross(e3, e4);

	abc(0) = arma::dot(A0v, n_tri) / Atot;
	abc(1) = arma::dot(A1v, n_tri) / Atot;
	abc(2) = arma::dot(A2v, n_tri) / Atot;

	if( fabs( 1 - abc(2) - abc(1) - abc(0) ) > 1e-10  )
	{
		printf("i_face: %d  \n", this->idx());
		std::cout << "abc: " << abc.t();
		std::cout << "xy:  " << xy.t() ;
		std::cout << "r1: " << v[0]->xyz.t();
		std::cout << "r2: " << v[1]->xyz.t();
		std::cout << "r3: " << v[2]->xyz.t();
		return false;
	}

	return true;
}

void Face::abc2xy(const Vector3 &abc, Vector3& xyz) const
{
	Vertex *v[3] = {he->vert, he->next->vert, he->next->next->vert};
	assert( fabs( 1 - abc(2) - abc(1) - abc(0) ) < 1e-10 );
	xyz = abc(0) * v[0]->xyz + abc(1) * v[1]->xyz + abc(2) * v[2]->xyz;
}

/****************************************************************************\
                                Constructor and Destructor
\****************************************************************************/
Face Mesh::_hole;


Mesh::Mesh(){}

Mesh::~Mesh() {this->clear_memory();}

void Mesh::clear_memory()
{
	verts.clear_memory();
	edges.clear_memory();
	half_edges.clear_memory();
	faces.clear_memory();
}

void Mesh::init( const std::vector<double>& xyz, const std::vector<int>& tris)
{
	assert( tris.size() % 3 == 0 );
	assert( xyz.size() % 3 == 0 );

	typedef std::map< std::pair<Vertex *, Vertex*>, HalfEdge* > EdgeMap;
	this->clear_memory();
	const int n_tris = int(tris.size() / 3);
	const int n_verts = int(xyz.size() / 3);
	EdgeMap edge_map;


	/*
	 * Add the vertices
	 */
	for (int idx_vert = 0 ; idx_vert < n_verts ; idx_vert++)
	{
		Vertex *v = new Vertex;
		v->xyz(0) = xyz[3*idx_vert+0];
		v->xyz(1) = xyz[3*idx_vert+1];
		v->xyz(2) = xyz[3*idx_vert+2];
		verts.add_member(v);
	}
	assert(n_verts == verts.n_mems());

	/*
	 * Add the faces and half edges.
	 */
	for( int tri = 0 ; tri < n_tris; tri++)
	{
		// get the vertices
		Vertex *fverts[3] =
		{
				verts.get_member_ptr(tris[3*tri+0]),
				verts.get_member_ptr(tris[3*tri+1]),
				verts.get_member_ptr(tris[3*tri+2])
		};

		// create half edges and the faces
		Face *face = new Face; 	faces.add_member(face);

		HalfEdge* he[3];
		he[0] = new HalfEdge( NULL, NULL, NULL, fverts[0], face); half_edges.add_member(he[0]);
		he[1] = new HalfEdge( NULL, NULL, NULL, fverts[1], face); half_edges.add_member(he[1]);
		he[2] = new HalfEdge( NULL, NULL, NULL, fverts[2], face); half_edges.add_member(he[2]);

		he[0]->next = he[1];
		he[1]->next = he[2];
		he[2]->next = he[0];
		face->he = he[0];
		fverts[0]->he = he[0];
		fverts[1]->he = he[1];
		fverts[2]->he = he[2];

		// get the twin for each of the half edges
		for (int i = 0 ; i < 3 ; i++)
		{
			int j = (i+1)%3;

			EdgeMap::iterator it = edge_map.find(std::make_pair( fverts[j], fverts[i] ) );
			if( (it != edge_map.end()))
			{
				 assert(it->first.first == fverts[j]);
				 assert(it->first.second == fverts[i]);
				 he[i]->twin = it->second;
				 it->second->twin = he[i];
				 edge_map.erase( it );
			}
			else
			{
				edge_map.insert( std::make_pair( std::make_pair(fverts[i], fverts[j]), he[i] ) );
			}
		}

	} // End of triangles
	assert(n_tris == faces.n_mems());


	/*
	 * Add half edges for the holes
	 */
	for( auto it = edge_map.begin() ; it != edge_map.end() ; )
	{
		HalfEdge* he_interior = it->second;
		HalfEdge* he_bdry = new HalfEdge(NULL, he_interior, NULL, it->first.second, hole());
		he_interior->twin = he_bdry;
		half_edges.add_member(he_bdry);
		it = edge_map.erase( it );
	} // End of residue iterators
	if(edge_map.size() != 0)
	{
		std::cerr << "[Error] Corrupted mesh, n_mems left: " << edge_map.size() << std::endl;
		throw;
	}

	/*
	 * Set he->next for hole half edges
	 */
	for (int he_idx = 0 ; he_idx < half_edges.n_mems() ; he_idx++)
	{
		HalfEdge *he = half_edges.get_member_ptr(he_idx);
		if(he->face != hole())
		{
			assert(he->next);
			continue;
		}

		assert(!he->next);
		HalfEdge *he_next = he->twin;
		while(he_next->face != hole())
		{
			he_next = he_next->next->next->twin;
		}
		he->next = he_next;
	} // End of hald edges

	/*
	 * Create the edges
	 */
	for (int he_idx = 0 ; he_idx < half_edges.n_mems() ; he_idx++)
	{
		HalfEdge *he = half_edges.get_member_ptr(he_idx);

		// Only do this once
		if(he->idx() > he->twin->idx())
		{
			assert(he->edge);
			continue;
		}

		Edge *edge = new Edge; edges.add_member(edge);

		he->edge = edge;
		he->twin->edge = edge;
		// Point to boundary first
		edge->he = (he->twin->face == hole() ? he->twin : he);

	} // End of hald edges

}

/****************************************************************************\
                                      Iterators
\****************************************************************************/

bool Mesh::ring_iterator::init(Vertex* v)
{
    if( !v->he ) return false;
    _cur = _end = v->he->twin; // point towards v
    return true;
}

bool Mesh::ring_iterator::reset_boundary()
{
	HalfEdge *snapshot = _cur;
    do
    {
        // start from the boundary edge located to the right of the vertex.
        if ( (_cur->face == Mesh::hole()) && (_cur->next->face == Mesh::hole()) )
        {
            _end = _cur;
            return true;
        }
        _cur = _cur->next->twin;
    } while (_cur != snapshot);
    return false;
}


bool Mesh::ring_iterator::advance()
{
    _cur = _cur->next->twin;
    return (_cur != _end);
}

/****************************************************************************\
                            Mesh Modifications
\****************************************************************************/

HalfEdge* Mesh::find_connecting_hedge( Vertex* vfrom, Vertex* vto )
{
    ring_iterator it;
    if( !it.init(vfrom) ) return NULL;

    do
    {
        if( it.half_edge()->vert == vto ) return it.half_edge()->twin;
    }while( it.advance() );

    return NULL;
}

int Mesh::count_connecting_hedges(Vertex* vfrom, Vertex* vto )
{
	int ans=0;
    ring_iterator it;
    if( !it.init(vfrom) ) return ans;

    do
    {
        if( it.half_edge()->vert == vto ) ans++;
    }while( it.advance() );

    return ans;
}

bool Mesh::collapse_edge(Edge* edge, Vertex* vertex)
{

	/*
	 * Cannot collapse boundary edges.
	 */
	if(  edge->he->face == hole() ) return false;

    /*
      Record the id's before starting the process.
    */
	HalfEdge *heBase = edge->he;
	HalfEdge *heTwin = heBase->twin;

    HalfEdge *heBorder[4];
    heBorder[0] = heBase->next->twin;
    heBorder[1] = heBase->next->next->twin;
    heBorder[2] = heTwin->next->next->twin;
    heBorder[3] = heTwin->next->twin;
    Edge *edgeBorder[4];
    edgeBorder[0] = heBorder[0]->edge;
    edgeBorder[1] = heBorder[1]->edge;
    edgeBorder[2] = heBorder[2]->edge;
    edgeBorder[3] = heBorder[3]->edge;

    /*
     * Cannot collapse cases where border edges are on the boundary.
     * Need to add more code.
     */
    if( heBorder[0]->face==hole() || heBorder[1]->face==hole() ||
    		heBorder[2]->face==hole() || heBorder[3]->face==hole() ) return false;

    /*
      Check that the collapse does not screw the data structure
    */
    if( heBorder[1]->next->twin->next == heBorder[0] ) return false;
    if( heBorder[2]->next->twin->next == heBorder[3] ) return false;
    if( this->count_connecting_hedges(heBorder[3]->vert, heBorder[0]->vert) > 1 ) return false;
    if( this->count_connecting_hedges(heTwin->vert, heBase->vert) > 1 ) return false;

    // Capture the indices of things (2 faces & 6 half-edges) we want
    // to delete.
    Face* fToDelete[] = { heBase->face, heTwin->face };
    HalfEdge* heToDelete[] = { heBase, heBase->next, heBase->next->next, heTwin, heTwin->next, heTwin->next->next };
    Vertex* vertToDelete[] = {heBase->vert, heTwin->vert};
    Edge* edgeToDelete[] = {edge, edgeBorder[1], edgeBorder[3]};

    // We can't be deleting border edges!
    if( std::find( heBorder, heBorder + 4, heToDelete[0] ) != heBorder + 4 ) return false;
    if( std::find( heBorder, heBorder + 4, heToDelete[1] ) != heBorder + 4 ) return false;
    if( std::find( heBorder, heBorder + 4, heToDelete[2] ) != heBorder + 4 ) return false;
    if( std::find( heBorder, heBorder + 4, heToDelete[3] ) != heBorder + 4 ) return false;
    if( std::find( heBorder, heBorder + 4, heToDelete[4] ) != heBorder + 4 ) return false;
    if( std::find( heBorder, heBorder + 4, heToDelete[5] ) != heBorder + 4 ) return false;

    /*
      Adjust the connectivities:
      half edge twin.
      half edge vertex.
      vertex to half edge.
      -- No: face to half edge!
    */

    // Add the new vertex and all associated data.
    verts.add_member(vertex);
    Vertex* modVerts[] = { heBase->next->next->vert, vertex, heTwin->next->next->vert };

    // Half edge to vertex
    HalfEdge *heIt = heBase->next->twin->next;
    HalfEdge * heEnd = heTwin;
    for( ; heIt != heEnd; heIt = heIt->twin->next )
    {
        assert( heIt->vert == heTwin->vert );
        heIt->vert = vertex;
    }
    heIt = heTwin->next->twin->next;
    heEnd = heBase;
    for( ; heIt != heEnd; heIt = heIt->twin->next)
    {
        assert( heIt->vert == heBase->vert );
        heIt->vert = vertex;
    }

    // Vertex to half edge
    modVerts[0]->he = heBorder[0];
    modVerts[1]->he = heBorder[1];
    modVerts[2]->he = heBorder[3];

    // Half edge twin.
    heBorder[0]->twin = heBorder[1];
    heBorder[1]->twin = heBorder[0];
    heBorder[2]->twin = heBorder[3];
    heBorder[3]->twin = heBorder[2];

    // Edge
   heBorder[0]->edge = heBorder[1]->edge = edgeBorder[0];
   heBorder[2]->edge = heBorder[3]->edge = edgeBorder[2];
   edgeBorder[0]->he = heBorder[0];
   edgeBorder[2]->he = heBorder[2];

    /*
      Now update the add/remove data structure
    */

    // verts
    verts.remove_member(vertToDelete[0]);
    verts.remove_member(vertToDelete[1]);

    // faces
    faces.remove_member(fToDelete[0]);
    faces.remove_member(fToDelete[1]);

    // he's
    half_edges.remove_member(heToDelete[0]);
    half_edges.remove_member(heToDelete[1]);
    half_edges.remove_member(heToDelete[2]);
    half_edges.remove_member(heToDelete[3]);
    half_edges.remove_member(heToDelete[4]);
    half_edges.remove_member(heToDelete[5]);

    // edges
    edges.remove_member(edgeToDelete[0]);
    edges.remove_member(edgeToDelete[1]);
    edges.remove_member(edgeToDelete[2]);

    return true;
}

bool Mesh::split_edge(Edge* edge_in, Vertex* vert_in)
{

	HalfEdge *he[13], *he2t, *he6t;
	Edge *e[9];
	Vertex *v[6];
	Face *f[5];

	/*
	 * Boundary
	 */
	if(edge_in->he->face == hole())
	{

		/*
		 * Get all the current entities
		 */
		// half edges
		he[1] = edge_in->he; he[2] = he[1]->next;
		he[4] = he[1]->twin; he[5] = he[4]->next; he[6] = he[5]->next;
		he6t = he[6]->twin;

		// vertices
		v[1] = he[5]->vert;
		v[2] = he[4]->vert;
		v[4] = he[6]->vert;
		v[5] = vert_in;

		// faces
		f[2] = he[4]->face;

		// edges
		e[1] = he[4]->edge;
		e[4] = he[5]->edge;
		e[5] = he[6]->edge;

		/*
		 * Create all the new entities
		 */
		{
			int i;
			i = 7; he[i] = new HalfEdge();half_edges.add_member(he[i]);
			i = 10; he[i] = new HalfEdge();half_edges.add_member(he[i]);
			i = 11; he[i] = new HalfEdge();half_edges.add_member(he[i]);
			i = 12; he[i] = new HalfEdge();half_edges.add_member(he[i]);
			i = 6 ; e[i] = new Edge();  edges.add_member(e[i]);
			i = 8 ; e[i] = new Edge();  edges.add_member(e[i]);
			i = 3 ; f[i] = new Face(); faces.add_member(f[i]);
			verts.add_member(v[5]);
		}


		/*
		 * Fix the connectivity
		 */

		// ---------------------------- half_edges

		he[4]->vert = he[11]->vert = he[7]->vert = v[5];
		he[10]->vert = v[2];
		he[12]->vert = v[4];

		he[6]->edge = he[11]->edge = e[8];
		he[10]->edge = he[7]->edge = e[6];
		he[12]->edge = e[5];

		he[7]->face = Mesh::hole();
		he[10]->face = he[11]->face = he[12]->face = f[3];

		he[1]->next = he[7]; he[7]->next = he[2];
		he[10]->next = he[11]; he[11]->next = he[12]; he[12]->next = he[10];

		he6t->twin = he[12]; he[12]->twin = he6t;
		he[7]->twin = he[10]; he[10]->twin = he[7];
		he[11]->twin = he[6]; he[6]->twin = he[11];

		// ---------------------------- faces
		f[3]->he = he[10];

		// ---------------------------- verts
		v[5]->he = he[7];
		v[2]->he = he[10];

		// ---------------------------- edges
		e[5]->he = he6t;
		e[6]->he = he[7];
		e[8]->he = he[6];

	} // end of if( edge is boundary )

	/*
	 * Internal
	 */
	else
	{

	/*
	 * Get all the current entities
	 */
	// half edges
	he[1] = edge_in->he; he[2] = he[1]->next; he[3] = he[2]->next;
	he[4] = he[1]->twin; he[5] = he[4]->next; he[6] = he[5]->next;
	he2t = he[2]->twin;  he6t = he[6]->twin;

	// vertices
	v[1] = he[1]->vert;
	v[2] = he[4]->vert;
	v[3] = he[3]->vert;
	v[4] = he[6]->vert;
	v[5] = vert_in;

	// faces
	f[1] = he[1]->face;
	f[2] = he[4]->face;

	// edges
	e[1] = he[1]->edge;
	e[2] = he[2]->edge;
	e[3] = he[3]->edge;
	e[4] = he[5]->edge;
	e[5] = he[6]->edge;

	/*
	 * Create all the new entities
	 */
	for (int i = 7 ; i <= 12 ; i++)
	{
		he[i] = new HalfEdge();
		half_edges.add_member(he[i]);
	}
	for (int i = 6 ; i <= 8 ; i++)
	{
		 e[i] = new Edge();
		 edges.add_member(e[i]);
	}
	for (int i = 3 ; i <= 4 ; i++)
	{
		f[i] = new Face();
		faces.add_member(f[i]);
	}
	verts.add_member(v[5]);

	/*
	 * Fix the connectivity
	 */

	// ---------------------------- half_edges

	he[4]->vert = he[2]->vert = he[7]->vert = he[11]->vert = v[5];
	he[8]->vert = he[10]->vert = v[2];
	he[12]->vert = v[4];
	he[9]->vert  = v[3];

	he[6]->edge = he[11]->edge = e[8];
	he[2]->edge = he[9]->edge = e[7];
	he[10]->edge = he[7]->edge = e[6];
	he[8]->edge = e[2];
	he[12]->edge = e[5];

	he[7]->face = he[8]->face = he[9]->face = f[4];
	he[10]->face = he[11]->face = he[12]->face = f[3];

	he[7]->next = he[8]; he[8]->next = he[9]; he[9]->next = he[7];
	he[10]->next = he[11]; he[11]->next = he[12]; he[12]->next = he[10];

	he2t->twin = he[8]; he[8]->twin = he2t;
	he6t->twin = he[12]; he[12]->twin = he6t;
	he[7]->twin = he[10]; he[10]->twin = he[7];
	he[11]->twin = he[6]; he[6]->twin = he[11];
	he[9]->twin = he[2]; he[2]->twin = he[9];

	// ---------------------------- faces
	f[4]->he = he[7];
	f[3]->he = he[10];

	// ---------------------------- verts
	v[5]->he = he[7];
	v[2]->he = he[8];

	// ---------------------------- edges
	e[2]->he = he2t;
	e[5]->he = he6t;
	e[6]->he = he[7];
	e[7]->he = he[2];
	e[8]->he = he[6];

	} // End of if( edge is internal )

	/*
	 * All done
	 */
	return true;
}

bool Mesh::flip_edge(Edge* edge)
{

	HalfEdge *he1, *he2, *he3, *he4, *he5, *he6;
	Face *f1, *f2;
	Vertex *v1, *v2, *v3, *v4;

	// Half edges
	he1 = edge->he;
	he2 = he1->next;
	he3 = he2->next;
	he4 = he1->twin;
	he5 = he4->next;
	he6 = he5->next;

	// Vertices
	v1 = he1->vert;
	v2 = he4->vert;
	v3 = he3->vert;
	v4 = he6->vert;

	// Faces
	f1 = he1->face;
	f2 = he4->face;

	// Fix connectivity for faces
	f1->he = he5;
	f2->he = he2;

	// Fix connectivity for vertices
	v1->he = he5;
	v2->he = he2;

	// Fix connectivity for half edges
	he1->next = he3; he3->next = he5; he5->next = he1;
	he4->next = he6; he6->next = he2; he2->next = he4;
	he1->vert = v4;
	he4->vert = v3;
	he3->face = f1;
	he5->face = f1;
	he2->face = f2;
	he6->face = f2;


	return true;
}

void Mesh::verify() const
{
	int i, iEnd;

	/*
	 * Check faces
	 */
    for( i = 0, iEnd = faces.n_mems(); i < iEnd; ++i )
    {
    	int c = 0;
       	Face *face = faces.get_member_ptr(i);
    	HalfEdge* it = face->he;
        assert( it->next != face->he );
        do
        {
            assert( it->face == face );
            assert( half_edges.is_member(it->next) );
            assert( half_edges.is_member(it->twin) );
            assert( verts.is_member(it->vert) );
            assert( ( it->twin->face == hole() ) || it->next->twin->face != it->face );
            it = it->next;
            c++;
        }while( it != face->he);
        if(c != 3) printf("[Error] c = %d \n", c);
        assert(c==3);
    }

    /*
     * Check vertices
     */
    for( i = 0, iEnd = verts.n_mems(); i < iEnd; ++i )
    {
    	Vertex *v = verts.get_member_ptr(i);
        assert( half_edges.is_member(v->he) );
        HalfEdge* it = v->he;
        assert( it->vert == v );
        assert( it->face != hole() );
    }

    /*
     * Check half edges
     */
    for( i = 0, iEnd = half_edges.n_mems(); i < iEnd; ++i )
    {
        const HalfEdge* it = half_edges.get_member_ptr(i);
        assert( verts.is_member(it->vert) );
        assert( ( it->face == hole() ) || faces.is_member(it->face) );

        assert( half_edges.is_member(it->next) );
        assert( it->next != it );
        assert( it->next->face == it->face );
        assert( it->next->vert != it->vert );

        assert( half_edges.is_member(it->twin) );
        assert( it->twin->twin == it );
        assert( it->twin->face != it->face );
        assert( it->twin->vert == it->next->vert );
    }

    /*
     * Check edges
     */
    assert(edges.n_mems() == half_edges.n_mems()/2);
    for(i = 0 ; i < edges.n_mems() ; i++)
    {
    	Edge *e = edges.get_member_ptr(i);
    	HalfEdge *he = e->he;

    	assert(he->edge == e);
    	assert(he->twin->edge == e);
    	assert(he->twin->face != hole());
    }
}

// Only works for 2d mesh
bool Mesh::locate_point(const Vector3 &xy, Face *&face, Vector3& abc) const
{
	assert( faces.is_member(face) );

	bool ierr;
	ierr = face->xy2abc(xy, abc);
	if(!ierr) return false;
	Face *next[] = {face->he->next->twin->face, face->he->next->next->twin->face, face->he->twin->face};

	// debug
	// std::cout << "checking face " << face->idx() << " abc are " << abc.t();

	if( (abc(0) > -1e-10) && (abc(1) > -1e-10) && (abc(2) > -1e-10) )
	{
			return true;
	}
	else if( (abc(0) < 0) && (next[0] != hole()) )
	{
		face = next[0];
		return this->locate_point(xy, face, abc);
	}
	else if( (abc(1) < 0) && (next[1] != hole()) )
	{
		face = next[1];
		return this->locate_point(xy, face, abc);
	}
	else if( (abc(2) < 0) && (next[2] != hole()) )
	{
		face = next[2];
		return this->locate_point(xy, face, abc);
	}

	return false;
}

void Mesh::get_boundary_info(IndexedList<BdryVertex>& boundary_verts, IndexedList<BdryEdge>& boundary_edges, const int two_way_map_index)
{

	boundary_verts.clear_memory();
	boundary_edges.clear_memory();

	//Loop over all the edges and create the edge to bdry_edge mapping.
    for (int i=0 ; i < edges.n_mems() ; i++)
    {
        Edge *edge= edges.get_member_ptr(i);
        if(two_way_map_index >= 0)
        {
        	if((int)edge->misc.size() <= two_way_map_index) edge->misc.resize(two_way_map_index+1);
        	edge->misc[two_way_map_index] = (void*)NULL;
        }

        // Only if on the boundary
        if(edge->he->face != hole()) continue;
        BdryEdge *bedge = new BdryEdge();
        bedge->edge = edge;
        boundary_edges.add_member(bedge);

        // Only if two_way_map_index is set
        if(two_way_map_index >= 0)
        	edge->misc[two_way_map_index] = static_cast<void*>(bedge);
    }

    //Loop over all vertices
    ring_iterator ringit;

    for (int i=0 ; i < verts.n_mems() ; i++)
    {
    	Vertex *vert= verts.get_member_ptr(i);
    	if(two_way_map_index >= 0)
    	{
    		if((int)vert->misc.size() <= two_way_map_index) vert->misc.resize(two_way_map_index+1);
    		vert->misc[two_way_map_index] = (void*)NULL;
    	}

    	// Only if on the boundary
    	ringit.init(vert);
		if(!ringit.reset_boundary()) continue;
		BdryVertex *bvert = new BdryVertex();
		bvert->vertex = vert;
		boundary_verts.add_member(bvert);

		// Only if two_way_map_index is set
		if(two_way_map_index >= 0)
			vert->misc[two_way_map_index] = static_cast<void*>(bvert);
    }

} // All done
