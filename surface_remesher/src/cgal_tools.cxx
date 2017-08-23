/*
 * cgal_tools.cxx
 *
 *  Created on: Dec 15, 2016
 *      Author: shayan
 *
 * id: http://cgal-discuss.949826.n4.nabble.com/Getting-facet-indexes-from-polyhedron-3-td4553195.html
 * builder: cgal-discuss.949826.n4.nabble.com/Getting-facet-indexes-from-polyhedron-3-td4553195.html
 *
 */


#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Parameterization_polyhedron_adaptor_3.h>
#include <CGAL/parameterize.h>

#include <sstream>
#include <iostream>
#include <algorithm>
#include <fstream>

#include "GetPot"
#include "armadillo"

#include "cgal_tools.hxx"
#include "mesh_io.hxx"
#include "pointer_mesh.hxx"

using namespace std;

typedef CGAL::Simple_cartesian<double>     CGAL_Kernel;
typedef CGAL::Polyhedron_3<CGAL_Kernel, CGAL::Polyhedron_items_with_id_3>  CGAL_Polyhedron;
typedef CGAL_Polyhedron::HalfedgeDS   CGAL_HalfedgeDS;


// A modifier converting a Mesh mesh to CGAL Polyhedron_3
template <class HDS>
class CGAL_Mesh2Poly_Object : public CGAL::Modifier_base<HDS>
{
public:
	CGAL_Mesh2Poly_Object(Mesh &mesh) : m_mesh(mesh) {}

    void operator()(HDS& hds)
    {

        // Postcondition: `hds' is a valid polyhedral surface.
        CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true);

        // begin surface
        B.begin_surface(m_mesh.verts.n_mems(), m_mesh.faces.n_mems());

        // vertices
        typedef typename HDS::Vertex::Point CGAL_Vertex;
        for ( int i = 0 ; i < m_mesh.verts.n_mems() ; i++ )
        {
        	Vertex *vert = m_mesh.verts.get_member_ptr(i);
        	typename HDS::Vertex_handle vh =
        			B.add_vertex(CGAL_Vertex(vert->xyz(0), vert->xyz(1), vert->xyz(2)));
        	vh->id() = i;
        }

        // triangles
        for ( int i = 0 ; i < m_mesh.faces.n_mems() ; i++ )
        {
        	Face *face = m_mesh.faces.get_member_ptr(i);
        	typename HDS::Face_handle fh = B.begin_facet();
            B.add_vertex_to_facet(face->he->vert->idx());
            B.add_vertex_to_facet(face->he->next->vert->idx());
            B.add_vertex_to_facet(face->he->next->next->vert->idx());
            B.end_facet();
            fh->id() = i;
        }

        // end surface
        B.end_surface();
    }

private:
    Mesh& m_mesh;
};


void CGAL_mesh2polyhedron(Mesh &mesh, CGAL_Polyhedron &P)
{
	CGAL_Mesh2Poly_Object<CGAL_HalfedgeDS> builder(mesh);
    P.delegate(builder);
}

void CGAL_polyhedron2mesh(CGAL_Polyhedron &P, Mesh &mesh)
{
	std::vector<double> coords;
	std::vector<int> conn;

	for(CGAL_Polyhedron::Vertex_iterator it = P.vertices_begin(); it != P.vertices_end() ; ++it)
	{
		coords.push_back( it->point().x() );
		coords.push_back( it->point().y() );
		coords.push_back( it->point().z() );
	}

	for(CGAL_Polyhedron::Face_iterator it = P.facets_begin(); it != P.facets_end() ; ++it)
	{
		CGAL_Polyhedron::Halfedge_around_facet_circulator circ = it->facet_begin();
		do
		{
			conn.push_back( int(circ->vertex()->id()) );
		} while ( ++circ != it->facet_begin());
	}

	mesh.init(coords, conn);
}

bool CGAL_parametrize(Mesh &mesh, CGAL_Polyhedron& P, std::string& errmsg)
{

	    // adaptor
	    typedef CGAL::Parameterization_polyhedron_adaptor_3<CGAL_Polyhedron>
	    Parameterization_polyhedron_adaptor;
	    Parameterization_polyhedron_adaptor mesh_adaptor(P);

	    // parametrizer
	    typedef CGAL::Parameterizer_traits_3<Parameterization_polyhedron_adaptor> Parameterizer;
	    Parameterizer::Error_code err = CGAL::parameterize(mesh_adaptor);

	    switch(err)
	    {
	    case Parameterizer::OK:
	    	break;
	    case Parameterizer::ERROR_EMPTY_MESH:
	    	errmsg = "empty mesh.";
	    	return false;
	    	break;
	    case Parameterizer::ERROR_NON_TRIANGULAR_MESH:
	    	errmsg = "non triangular mesh.";
	    	return false;
	    	break;
	    case Parameterizer::ERROR_NO_TOPOLOGICAL_DISC:
	    	errmsg = "mesh is not topological to a disk.";
	    	return false;
	    	break;
	    case Parameterizer::ERROR_BORDER_TOO_SHORT:
	    	errmsg = "border is too short.";
	    	return false;
	    	break;
	    default:
	    	errmsg = Parameterizer::get_error_message(err);
	    	return false;
	    	break;
	    };

	    // Set U and V
	    CGAL_Polyhedron::Vertex_const_iterator pVertex;
	    for (pVertex = P.vertices_begin(); pVertex != P.vertices_end();	pVertex++)
	    {
	    	// (u,v) pair is stored in any halfedge
	    	const double u = mesh_adaptor.info(pVertex->halfedge())->uv().x();
	    	const double v = mesh_adaptor.info(pVertex->halfedge())->uv().y();
	    	Vertex *vert = mesh.verts.get_member_ptr( int( pVertex->id() ) );
	    	vert->xyz(0) = u;
	    	vert->xyz(1) = v;
	    	vert->xyz(2) = 0;
	    }

	    return true;
} // End CGAL_parametrize

bool CGAL_circular_map( Mesh& mesh, std::string& errmsg)
{
	// Create a CGAL mesh
	CGAL_Polyhedron P;
	CGAL_mesh2polyhedron(mesh, P);

	// Parametrize
	bool scss;
	scss = CGAL_parametrize(mesh, P, errmsg);

	return scss;
}


int CGAL_map_test(int argc, char *argv[])
{
	GetPot getpot(argc,argv);
	Mesh mesh, m2;
	MeshIO io(mesh), io2(m2);
	AdaptiveRemesher arm(mesh,getpot);
	Mesh::ring_iterator it;
	stringstream ss;
	PatchProjector proj(mesh);

	// Read the mesh
	string iname = getpot.follow("meshes/square.off", "-i");
	io.read_auto(iname);

	// Create a CGAL mesh
	CGAL_Polyhedron P;
	CGAL_mesh2polyhedron(mesh, P);

    // Create the second mesh
    CGAL_polyhedron2mesh(P, m2);
    io2.write_vtk("cgal2.vtk");

    // test
    for (int i = 0 ; i < m2.faces.n_mems() ; i++)
    {
    	Face *f0 = mesh.faces.get_member_ptr(i);
    	Face *f1 = m2.faces.get_member_ptr(i);
    	Vertex *fv0[3] = {f0->he->vert, f0->he->next->vert, f0->he->next->next->vert};
    	Vertex *fv1[3] = {f0->he->vert, f1->he->next->vert, f1->he->next->next->vert};

    	if(i == 0)
    	{
    		std::cout << fv0[0]->idx() << " " << fv0[1]->idx() << " " << fv0[2]->idx() << std::endl;
    		std::cout << fv1[0]->idx() << " " << fv1[1]->idx() << " " << fv1[2]->idx() << std::endl;
    	}

    	assert( fv0[0]->idx() == fv1[0]->idx() );
    	assert( fv0[1]->idx() == fv1[1]->idx() );
    	assert( fv0[2]->idx() == fv1[2]->idx() );
    }

    std::string msg;
    bool scss;
    scss = CGAL_parametrize(m2, P, msg);
    if(!scss) std::cout << "[Error] parametrizer: " << msg << std::endl;
    io2.write_vtk("cgal3.vtk");

	return 0;

}


