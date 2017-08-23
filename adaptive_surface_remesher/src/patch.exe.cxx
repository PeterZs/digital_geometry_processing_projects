/*
 * patch.exe.cxx
 *
 * Test creation of patches
 * id: http://cgal-discuss.949826.n4.nabble.com/Getting-facet-indexes-from-polyhedron-3-td4553195.html
 * builder: cgal-discuss.949826.n4.nabble.com/Getting-facet-indexes-from-polyhedron-3-td4553195.html
 */


#include <sstream>
#include <iostream>
#include <algorithm>
#include <fstream>

#include "cgal_tools.hxx"
#include "mesh_io.hxx"
#include "adaptive_remesher.hxx"

#include "GetPot"
#include "armadillo"

using namespace std;


/*
 * Test one round of splitting (edge split)
 */
int test1(int argc, char *argv[])
{
	GetPot getpot(argc,argv);
	Mesh mesh;
	MeshIO io(mesh);
	AdaptiveRemesher arm(mesh,getpot);
	Mesh::ring_iterator it;
	stringstream ss;
	PatchProjector proj(mesh);
	FILE *fl;

	// Read the mesh
	string iname = getpot.follow("meshes/square.off", "-i");
	io.read_auto(iname);

	/*
	 * Initialize the projection object
	 */
	proj.init();

	/*
	 * Write the mesh
	 */
	fl = MeshIO::open_file("output1.vtk", ".vtk", "w");
	io.write_vtk(fl);
	fclose(fl);

	/*
	 * Create a patch
	 */
	Face face;
	HalfEdge hes[3];
	Vertex verts[3];

	//int targets[] = {650, 675, 830};
	//int targets[] = {723, 839, 840};
	// int targets[] = {1700, 1708, 1709};
	// int targets[] = {1679, 1708, 1709}; // cow
	// int targets[] = {9936, 5254, 123}; // camel
	// int targets[] = {9935, 9935, 9934}; // camel
	int targets[] = {1822, 1823, 1815};

	verts[0].facer = mesh.faces.get_member_ptr(targets[0]); verts[0].abcr = "0.25 0.25 0.5";
	verts[1].facer = mesh.faces.get_member_ptr(targets[1]); verts[1].abcr = "0.25 0.5 0.25";
	verts[2].facer = mesh.faces.get_member_ptr(targets[2]); verts[2].abcr = "0.5 0.25 0.25";

	hes[0].next = &hes[1]; hes[1].next = &hes[2]; hes[2].next = &hes[0];
	hes[0].vert = &verts[0]; hes[1].vert = &verts[1]; hes[2].vert = &verts[2];

	face.he = &hes[0];

	Mesh *patch; Face *pfaces[3];
	proj.create_patch(&face, patch, pfaces);
	// MeshIO(*patch).write_vtk("patch0.vtk");

	bool scss; string errmsg;
	scss = CGAL_circular_map(*patch, errmsg);
	if(!scss) std::cout << "[Error] parametrizer" << errmsg;
	// MeshIO(*patch).write_vtk("patch1.vtk");

	return 0;
}

int main(int argc, char *argv[])
{

	test1(argc, argv);
	//test_CGAL_in(argc, argv);
	return 0;
}

