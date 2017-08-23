/*
 * exec_subdiv.cxx
 *
 * This executable reads a mesh, subdivides for a number of times and
 * writes each subdivided mesh for visualization.
 */

#include <sstream>
#include <iostream>

#include "parameters.hxx"
#include "mesh_io.hxx"
#include "adaptive_remesher.hxx"

#include "GetPot"
#include "armadillo"


using std::cout;
using std::endl;
using std::string;
using std::stringstream;



int test1(int argc, char *argv[])
{
	GetPot getpot(argc, argv);
	Mesh mesh, m2d;
	MeshIO io(mesh);
	AdaptiveRemesher arm(mesh, getpot);
	string name = getpot.follow("output", "-o");
	stringstream ss;
	PatchProjector *proj = dynamic_cast<PatchProjector*>( arm.proj() );
	int i , j = 0;

	io.read_auto("meshes/test/big.off");
	ss.str(""); ss << name << j << ".vtk"; j++;	io.write_vtk(ss.str());

	arm.init();

	Vector3 xyz;

//	i=5; arm.remesh_areabased(mesh.verts.get_member_ptr(i),  m2d);
//	ss.str(""); ss << name << j << ".vtk"; j++;	io.write_vtk(ss.str());
//	if( proj ) cout << "total n_patches " << proj->patches.size() << endl;
//
//	i=12; arm.remesh_areabased(mesh.verts.get_member_ptr(i),  m2d);
//	ss.str(""); ss << name << j << ".vtk"; j++;	io.write_vtk(ss.str());
//	if( proj ) cout << "total n_patches " << proj->patches.size() << endl;

	i=8; arm.remesh_areabased(mesh.verts.get_member_ptr(i),  m2d);
	ss.str(""); ss << name << j << ".vtk"; j++;	io.write_vtk(ss.str());
	if( proj ) cout << "total n_patches " << proj->patches.size() << endl;

	i=9; arm.remesh_areabased(mesh.verts.get_member_ptr(i),  m2d);
	ss.str(""); ss << name << j << ".vtk"; j++;	io.write_vtk(ss.str());
	if( proj ) cout << "total n_patches " << proj->patches.size() << endl;

	i=8; arm.remesh_areabased(mesh.verts.get_member_ptr(i),  m2d);
	ss.str(""); ss << name << j << ".vtk"; j++;	io.write_vtk(ss.str());
	if( proj ) cout << "total n_patches " << proj->patches.size() << endl;

	i=9; arm.remesh_areabased(mesh.verts.get_member_ptr(i),  m2d);
	ss.str(""); ss << name << j << ".vtk"; j++;	io.write_vtk(ss.str());
	if( proj ) cout << "total n_patches " << proj->patches.size() << endl;

	if( proj )
	{
		cout << "total n_patches " << proj->patches.size() << endl;

		for (auto itp = proj->patches.begin() ; itp != proj->patches.end() ; ++itp)
		{
			ss.str(""); ss << name << j << ".vtk"; j++;	MeshIO(**itp).write_vtk(ss.str());
		}
	}

	return 0;
}

int test2(int argc, char *argv[])
{
	GetPot getpot(argc, argv);
	Mesh mesh, m2d;
	MeshIO io(mesh);
	AdaptiveRemesher arm(mesh, getpot);
	stringstream ss;


	string iname = getpot.follow("meshes/test/big.off", "-i");
	io.read_auto(iname);
	io.write_vtk("output1.vtk");

	int i;
	for(i = 0 ; i < mesh.verts.n_mems() ; i++)
	{
		Vertex *v = mesh.verts.get_member_ptr(i);
		if(!arm.remesh_areabased(v , m2d)) continue;

		ss.str("");
		ss << "output" << i+2 << ".vtk";
		io.write_vtk(ss.str());
	}

	arm.remesh_flipping();
	ss.str("");
	ss << "output" << i+2 << ".vtk";
	io.write_vtk(ss.str());


	return 0;
}

int test3(int argc, char *argv[])
{
	GetPot getpot(argc, argv);
	Mesh mesh;
	MeshIO io(mesh);
	AdaptiveRemesher arm(mesh,getpot);
	stringstream ss;

	string iname = getpot.follow("meshes/test/big.off", "-i");
	io.read_auto(iname);
	io.write_vtk("output0.vtk");

	for (int i = 0 ; i < 10 ; i++)
	{
		arm.remesh_areabased();
		ss.str(""); ss << "output" << 2*i+1 << ".vtk"; io.write_vtk(ss.str());

		arm.remesh_flipping();
		ss.str(""); ss << "output" << 2*i+2 << ".vtk"; io.write_vtk(ss.str());
	}

	return 0;
}

int test4(int argc, char *argv[])
{
	GetPot getpot(argc, argv);
	Mesh mesh;
	MeshIO io(mesh);
	AdaptiveRemesher arm(mesh, getpot);
	stringstream ss;

	string iname = getpot.follow("meshes/test/big.off", "-i");
	io.read_auto(iname);
	io.write_vtk("output0.vtk");

	const int n_step = 10;
	const int n_area = 3;
	int file = 1;

	arm.arearemesh_mode() = AdaptiveRemesher::CONSTANT;

	for (int i = 0 ; i < n_step ; i++)
	{
		for (int j = 0 ; j < n_area ; j++)
		{
			arm.remesh_areabased();
			ss.str(""); ss << "output" << file << ".vtk"; io.write_vtk(ss.str()); file++;
		}

		arm.remesh_flipping();
		ss.str(""); ss << "output" << file << ".vtk"; io.write_vtk(ss.str()); file++;
	}

	return 0;
}

int main(int argc, char *argv[])
{

	test1(argc, argv);
	 //test2(argc, argv);
	// test3(argc, argv);
	// test4(argc, argv);

	return 0;
}
