/*
 * exec_subdiv.cxx
 *
 * This executable reads a mesh, subdivides for a number of times and
 * writes each subdivided mesh for visualization.
 */

#include <sstream>

#include "parameters.hxx"
#include "mesh_io.hxx"
#include "indexed_list.hxx"
#include "GetPot"

using std::cout;
using std::endl;
using std::string;
using std::stringstream;

/*
 * Test the IndexedList and IndexedItem classes.
 */
class Dummy: public IndexedItem<Dummy>
{
	public:
	double v;
	Dummy(double v_in):v(v_in){}
};

int test_index(int argc, char *argv[])
{
	IndexedList<Dummy> list;
	list.add_member(new Dummy(2.));
	list.add_member(new Dummy(3.));
	list.add_member(new Dummy(4.));

	for (int i = 0 ; i < list.n_mems() ; i++)
	{
		Dummy *ptr = list.get_member_ptr(i);
		int idx = list.get_member_idx(ptr);

		printf("%d %d %lf \n", i, idx, ptr->v);
	}

	list.remove_member(1); // 2 4
	list.add_member(new Dummy(8.)); // 2 4 8
	list.add_member(new Dummy(9.)); // 2 4 8 9
	list.remove_member(0); // 9 4 8
	list.add_member(new Dummy(12)); // 9 4 8 12

	printf("N new members: %d \n", list.n_mems());
	for (int i = 0 ; i < list.n_mems() ; i++)
	{
		Dummy *ptr = list.get_member_ptr(i);
		int idx = list.get_member_idx(ptr);

		printf("%d %d %lf \n", i, idx, ptr->v);
	}


	list.clear_memory();
	return 0;
}

/*
 * Test reading and writing a cartel mesh.
 */
int test_io(int argc, char *argv[])
{
	// getpot reads the command line
	GetPot getpot(argc, argv);
	// The mesh
	Mesh mesh;
	AdaptiveRemesher arm(mesh, getpot);
	MeshIO io(mesh);

	// Input file name
	const string fname_in = getpot.follow("DUMMY", "-i");

	// Read the mesh
	io.read_auto(fname_in);
	arm.init();

	// write the mesh
	stringstream ss;
	ss.str("");
	ss << "output";
	FILE *fl = MeshIO::open_file(ss.str(), ".vtk", "w");
	io.write_vtk(fl);
	io.write_vtk_data(fl, arm, "vert_normals");
	io.write_vtk_data(fl, arm, "tri_quality");
	io.write_vtk_data(fl, arm, "tri_fidelity");
	fclose(fl);
	if(getpot.search("-p")) MeshIO(mesh).print_info();

	ss << "_dual";
	io.write_vtk_dual(ss.str());

	return 0;
}

/*
 * Test the ring iterator
 */
int test_ringiterator(int, char *[])
{
	Mesh mesh;
	MeshIO io(mesh);

	io.read_auto("meshes/test/big.off");
	io.write_vtk("output.vtk");

	Mesh::ring_iterator it;
	Vertex *v5  = mesh.verts.get_member_ptr(5);
	Vertex *v4 = mesh.verts.get_member_ptr(4);

	// Test inside
	cout << "init 5: " << it.init(v5) << endl;
	do
	{
		cout << it.half_edge()->vert->idx() << " ";
	}while(it.advance());
	cout << endl;
	cout << "reset boundary 5: " << it.reset_boundary() << endl;
	do
	{
		cout << it.half_edge()->vert->idx() << " ";
	}while(it.advance());
	cout << endl;

	// Test outside
	cout << "init 4: " << it.init(v4) << endl;
	do
	{
		cout << it.half_edge()->vert->idx() << " ";
	}while(it.advance());
	cout << endl;
	cout << "reset boundary 4: " << it.reset_boundary() << endl;
	do
	{
		cout << it.half_edge()->vert->idx() << " ";
	}while(it.advance());
	cout << endl;

	return 0;
}

int test_flip(int, char *[])
{
	Mesh mesh;
	MeshIO io(mesh);

	io.read_auto("meshes/test/big.off");
	io.write_vtk("output1.vtk");

	mesh.flip_edge(mesh.find_connecting_hedge( mesh.verts.get_member_ptr(4), mesh.verts.get_member_ptr(12) )->edge);
	io.write_vtk("output2.vtk");

	mesh.flip_edge(mesh.find_connecting_hedge( mesh.verts.get_member_ptr(8), mesh.verts.get_member_ptr(11) )->edge);
	io.write_vtk("output3.vtk");

	mesh.flip_edge(mesh.find_connecting_hedge( mesh.verts.get_member_ptr(9), mesh.verts.get_member_ptr(7) )->edge);
	io.write_vtk("output4.vtk");

	mesh.flip_edge(mesh.find_connecting_hedge( mesh.verts.get_member_ptr(4), mesh.verts.get_member_ptr(5) )->edge);
	io.write_vtk("output5.vtk");

	return 0;
}

int test_locator(int, char*[])
{
	Mesh mesh;
	MeshIO io(mesh);
	io.read_off("meshes/test/big.off");
	io.write_vtk("output.vtk");

	// Let's create a point
	Face *f_in = mesh.faces.get_member_ptr(20);
	Face *f_out = mesh.faces.get_member_ptr(1);

	Vector3 abc_in, abc_out, xy_in, xy_out;
	abc_in(0) = 0.1; abc_in(1) = 0.2; abc_in(2)=0.7;
	f_in->abc2xy(abc_in, xy_in);
	mesh.locate_point(xy_in, f_out, abc_out);
	f_out->abc2xy(abc_out, xy_out);

	std::cout << " input is face " << f_in->idx() << " abc " << abc_in.t();
	std::cout << " xy is " << xy_in.t();
	std::cout << " found point in face " << f_out->idx() << " abc " << abc_out.t() ;
	std::cout << " xy is " << xy_out.t();

	return 0;
}

int test_geodesic(int argc , char *argv[])
{
	GetPot getpot(argc, argv);
	Mesh mesh, m2;
	MeshIO io(mesh), io2(m2);

	AdaptiveRemesher arm(mesh, getpot);

	io.read_off("meshes/test/big.off");
	mesh.flip_edge(mesh.find_connecting_hedge( mesh.verts.get_member_ptr(0), mesh.verts.get_member_ptr(5) )->edge);
	io.write_vtk("output.vtk");

	arm.geodesic_vert(mesh.verts.get_member_ptr(5), m2);
	io2.write_vtk("geodesic.vtk");

	return 0;
}

int main(int argc, char *argv[])
{
	//test_index(argc, argv);
	test_io(argc, argv);
	//test_ringiterator(argc, argv);
	//test_flip(argc, argv);
	//test_locator(argc,argv);
	//test_geodesic(argc,argv);
	return 0;
}
