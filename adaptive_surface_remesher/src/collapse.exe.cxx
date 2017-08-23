/*
 * collapse.exe.cxx
 *
 * This executable reads a mesh, subdivides for a number of times and
 * writes each subdivided mesh for visualization.
 */

#include <sstream>
#include <iostream>
#include <algorithm>


#include "parameters.hxx"
#include "mesh_io.hxx"
#include "adaptive_remesher.hxx"

#include "GetPot"
#include "armadillo"


using std::cout;
using std::endl;
using std::string;
using std::stringstream;



///*
// * Some helper functions
// */
//static void collapse_data(Mesh &m, const int iv1, const int iv2, Edge *&eout, Vertex *&vout)
//{
//	Vertex *v1 = m.verts.get_member_ptr(iv1);
//	Vertex *v2 = m.verts.get_member_ptr(iv2);
//	vout  = new Vertex;
//	eout = m.find_connecting_hedge(v1, v2)->edge;
//	vout->xyz = (v1->xyz + v2->xyz) / 2. ;
//}
//
//static void area_remesh(Mesh& m,
//		const int n_step = 10, const int n_area = 3, const std::string name = "output", int file = 0)
//{
//
//	AdaptiveRemesher arm(m);
//	MeshIO io(m);
//	stringstream ss;
//	FILE *fl;
//
//	for (int i = 0 ; i < n_step ; i++)
//	{
//		for (int j = 0 ; j < n_area ; j++)
//		{
//			arm.remesh_areabased();
//			ss.str(""); ss << "output" << file << ".vtk";
//			fl = MeshIO::open_file(ss.str(), ".vtk", "w");
//			io.write_vtk(fl);
//			io.write_vtk_fidelity(fl, arm);
//			fclose(fl);
//			file++;
//		}
//
//		arm.remesh_flipping();
//		ss.str(""); ss << "output" << file << ".vtk";
//		fl = MeshIO::open_file(ss.str(), ".vtk", "w");
//		io.write_vtk(fl);
//		io.write_vtk_fidelity(fl, arm);
//		fclose(fl);
//		file++;
//	}
//}
//
///*
// * Test collapse
// */
//int test1(int argc, char *argv[])
//{
//	Mesh mesh;
//	MeshIO io(mesh);
//	AdaptiveRemesher arm(mesh);
//	Mesh::ring_iterator it;
//	FILE *fl;
//
//
//	io.read_auto("meshes/test/big.off");
//	fl = MeshIO::open_file("output1.vtk", ".vtk", "w");
//	io.write_vtk(fl);
//	io.write_vtk_fidelity(fl, arm);
//	fclose(fl);
//
//
//	Vertex *vmid;
//	Edge *e;
//	std::vector<int> idx;
//
//	collapse_data(mesh, 8, 5, e, vmid); 	mesh.verify();
//	mesh.collapse_edge(e,vmid);
//	fl = MeshIO::open_file("output2.vtk", ".vtk", "w");
//	io.write_vtk(fl);
//	io.write_vtk_fidelity(fl, arm);
//	fclose(fl);
//
//
//	collapse_data(mesh, 5, 9, e, vmid);     mesh.verify();
//	mesh.collapse_edge(e,vmid);
//	fl = MeshIO::open_file("output3.vtk", ".vtk", "w");
//	io.write_vtk(fl);
//	io.write_vtk_fidelity(fl, arm);
//	fclose(fl);
//
//	collapse_data(mesh, 5, 0, e, vmid);    mesh.verify();
//	mesh.split_edge(e,vmid);
//	fl = MeshIO::open_file("output4.vtk", ".vtk", "w");
//	io.write_vtk(fl);
//	io.write_vtk_fidelity(fl, arm);
//	fclose(fl);
//
//	arm.remesh_flipping();                  mesh.verify();
//	fl = MeshIO::open_file("output5.vtk", ".vtk", "w");
//	io.write_vtk(fl);
//	io.write_vtk_fidelity(fl, arm);
//	fclose(fl);
//
//	area_remesh(mesh, 10, 3, "output", 6);
//	mesh.verify();
//
//	return 0;
//}


/*
 * Test one round of splitting (edge split)
 */
int test2(int argc, char *argv[])
{
	GetPot getpot(argc, argv);
	Mesh mesh;
	MeshIO io(mesh);
	AdaptiveRemesher arm(mesh, getpot);
	Mesh::ring_iterator it;
	FILE *fl;
	stringstream ss;
	int file;

	arm.arearemesh_mode() = AdaptiveRemesher::CONSTANT;

	string iname = getpot.follow("meshes/test/big_ref.off", "-i");
	io.read_auto(iname);
	arm.init();

	fl = MeshIO::open_file("output1.vtk", ".vtk", "w");
	io.write_vtk(fl);
	io.write_vtk_data(fl, arm, "tri_quality");
	fclose(fl);

	file = 2;

	for (int i = 0 ; i < 1 ; i++)
	{
		arm.remesh_splitedges();
		ss.str(""); ss << "output" << file << ".vtk";
		fl = MeshIO::open_file(ss.str(), ".vtk", "w"); file++;
		io.write_vtk(fl);
		io.write_vtk_data(fl, arm, "tri_quality");
		fclose(fl);

		arm.remesh_areabased(1, 2);
		ss.str(""); ss << "output" << file;
		fl = MeshIO::open_file(ss.str(), ".vtk", "w"); file++;
		io.write_vtk(fl);
		io.write_vtk_data(fl, arm, "tri_quality");
		fclose(fl);
	}

	arm.remesh_areabased(10, 3);
	ss.str(""); ss << "output" << file;
	fl = MeshIO::open_file(ss.str(), ".vtk", "w"); file++;
	io.write_vtk(fl);
	io.write_vtk_data(fl, arm, "tri_quality");
	fclose(fl);

	for (int i = 0 ; i < 1 ; i++)
	{
		arm.remesh_collapseedges();
		ss.str(""); ss << "output" << file << ".vtk";
		fl = MeshIO::open_file(ss.str(), ".vtk", "w"); file++;
		io.write_vtk(fl);
		io.write_vtk_data(fl, arm, "tri_quality");
		fclose(fl);

		arm.remesh_areabased(1, 2);
		ss.str(""); ss << "output" << file;
		fl = MeshIO::open_file(ss.str(), ".vtk", "w"); file++;
		io.write_vtk(fl);
		io.write_vtk_data(fl, arm, "tri_quality");
		fclose(fl);
	}


	arm.remesh_areabased(10, 3);
	ss.str(""); ss << "output" << file;
	fl = MeshIO::open_file(ss.str(), ".vtk", "w"); file++;
	io.write_vtk(fl);
	io.write_vtk_data(fl, arm, "tri_quality");
	fclose(fl);


	return 0;
}

int main(int argc, char *argv[])
{

	test2(argc, argv);
	return 0;
}
