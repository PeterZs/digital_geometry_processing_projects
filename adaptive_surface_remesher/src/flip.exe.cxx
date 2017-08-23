/*
 * exec_subdiv.cxx
 *
 * This executable reads a mesh, subdivides for a number of times and
 * writes each subdivided mesh for visualization.
 */

#include <sstream>

#include "parameters.hxx"
#include "mesh_io.hxx"
#include "adaptive_remesher.hxx"

#include "GetPot"
#include "armadillo"

#include <iostream>

using std::cout;
using std::endl;
using std::string;
using std::stringstream;


int main(int argc, char *argv[])
{
	GetPot getpot(argc, argv);
	Mesh mesh;
	MeshIO io(mesh);
	AdaptiveRemesher arm(mesh, getpot);

//	io.read_auto("meshes/test/test.off");
//	io.write_vtk("output1.vtk");
//	io.print_info();
//
//	mesh.flip_edge(mesh.find_connecting_hedge( mesh.verts.get_member_ptr(1), mesh.verts.get_member_ptr(3) )->edge);
//	io.write_vtk("output2.vtk");
//	io.print_info();


//	io.read_auto("meshes/test/big.off");
//	io.write_vtk("output1.vtk");
//
//	mesh.flip_edge(mesh.find_connecting_hedge( mesh.verts.get_member_ptr(4), mesh.verts.get_member_ptr(12) )->edge);
//	io.write_vtk("output2.vtk");
//	io.print_info();
//
//	mesh.flip_edge(mesh.find_connecting_hedge( mesh.verts.get_member_ptr(4), mesh.verts.get_member_ptr(8) )->edge);
//	io.write_vtk("output3.vtk");
//	cout << endl;
//	io.print_info();
//
//	arm.remesh_flipping();
//	io.write_vtk("output4.vtk");


	std::string iname = getpot.follow("meshes/square.off", "-i");

	io.read_auto(iname);
	io.write_vtk("output1.vtk");

	//arm.random_flipping();
	//io.write_vtk("output2.vtk");

	arm.remesh_flipping();
	io.write_vtk("output3.vtk");

	return 0;
}
