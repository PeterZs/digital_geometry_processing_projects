/*
 * split.exe.cxx
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



/*
 * Some helper functions
 */
//static void split_data(Mesh &m, const int iv1, const int iv2, Edge *&eout, Vertex *&vout)
//{
//	Vertex *v1 = m.verts.get_member_ptr(iv1);
//	Vertex *v2 = m.verts.get_member_ptr(iv2);
//	vout  = new Vertex;
//	eout = m.find_connecting_hedge(v1, v2)->edge;
//	vout->xyz = (v1->xyz + v2->xyz) / 2. ;
//}

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
//			// io.write_vtk_fidelity(fl, arm);
//			fclose(fl);
//			file++;
//		}
//
//		arm.remesh_flipping();
//		ss.str(""); ss << "output" << file << ".vtk";
//		fl = MeshIO::open_file(ss.str(), ".vtk", "w");
//		io.write_vtk(fl);
//		// io.write_vtk_fidelity(fl, arm);
//		fclose(fl);
//		file++;
//	}
//}

/*
 * Test insertion (edge split)
 */
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
//	split_data(mesh, 0, 5, e, vmid);
//	cout << "splitting " << e->he->vert->idx() << " " << e->he->twin->vert->idx() << endl;
//	mesh.split_edge(e,vmid); idx.push_back(vmid->idx());
//	fl = MeshIO::open_file("output2.vtk", ".vtk", "w");
//	io.write_vtk(fl);
//	io.write_vtk_fidelity(fl, arm);
//	fclose(fl);
//
//
//	split_data(mesh, 4, 18, e, vmid);
//	mesh.split_edge(e,vmid); idx.push_back(vmid->idx());
//	fl = MeshIO::open_file("output3.vtk", ".vtk", "w");
//	io.write_vtk(fl);
//	//io.write_vtk_fidelity(fl, arm);
//	fclose(fl);
//
//	split_data(mesh, 0, 19, e, vmid);
//	mesh.split_edge(e,vmid); idx.push_back(vmid->idx());
//	fl = MeshIO::open_file("output4.vtk", ".vtk", "w");
//	io.write_vtk(fl);
//	//io.write_vtk_fidelity(fl, arm);
//	fclose(fl);
//
//	split_data(mesh, 5, 19, e, vmid);
//	mesh.split_edge(e,vmid); idx.push_back(vmid->idx());
//	fl = MeshIO::open_file("output5.vtk", ".vtk", "w");
//	io.write_vtk(fl);
//	io.write_vtk_fidelity(fl, arm);
//	fclose(fl);
//
//
//	area_remesh(mesh, 10, 3, "output", 4);
//
//	return 0;
//}

/*
 * Test index_list's misc and sort!
 */
//int test2(int argc, char *argv[])
//{
//
//	// Test the sort thing
//	IndexedList<Edge> edges;
//	double val[20];
//	Edge *ie[20];
//
//	for (int i = 0 ; i < 20 ; i++)
//	{
//		val[i] = sin( i ) * 1e3;
//		ie[i] = new Edge();
//		edges.add_member(ie[i]);
//		ie[i]->misc.push_back(&val[i]);
//		printf("v %5.0lf -> %4d \n", val[i], ie[i]->idx());
//	}
//	printf("\n");
//
//	std::sort(ie, ie+20, Edge::cmp_functor<double, 0>());
//
//	for (int i = 0 ; i < 20 ; i++)
//	{
//		ie[i]->misc.push_back(&val[i]);
//		printf("v %5.0lf -> %4d \n", *static_cast<double*>(ie[i]->misc[0]), ie[i]->idx());
//	}
//
//
//	return 0;
//}


/*
 * Test one round of splitting (edge split)
 */
int test3(int argc, char *argv[])
{
	GetPot getpot(argc,argv);
	Mesh mesh, m2;
	MeshIO io(mesh);
	AdaptiveRemesher arm(mesh,getpot);
	Mesh::ring_iterator it;
	FILE *fl;
	stringstream ss;
	int file = 0;

	string iname = getpot.follow("meshes/test/big.off", "-i");
	arm.arearemesh_mode() = AdaptiveRemesher::CONSTANT;

	io.read_auto(iname);
	arm.init();

        // write the file
	ss.str(""); ss << "output" << file;
	fl = MeshIO::open_file(ss.str(), ".vtk", "w"); file++;
	io.write_vtk(fl);
	io.write_vtk_data(fl, arm, "tri_quality");
	fclose(fl);

        // One round splitting
        arm.n_tverts() = 1e6;
        for (int i = 0 ; i < 0 ; i++)
        {

            // Split
            arm.remesh_splitedges();
            ss.str(""); ss << "output" << file;
            fl = MeshIO::open_file(ss.str(), ".vtk", "w"); file++;
            io.write_vtk(fl);
            io.write_vtk_data(fl, arm, "tri_quality");
            fclose(fl);

            // Area based
            for (int j = 0 ; j < 3 ; j++)
            {
                for (int k = 0 ; k < mesh.verts.n_mems() ; k++)
                {
                    Vertex *vert = mesh.verts.get_member_ptr(k);

                    if( arm.remesh_areabased(vert, m2) )
                    {
                        ss.str(""); ss << "output" << file;
                        fl = MeshIO::open_file(ss.str(), ".vtk", "w"); file++;
                        io.write_vtk(fl);
                        io.write_vtk_data(fl, arm, "tri_quality");
                        fclose(fl);
                    }
                }
            
                arm.remesh_flipping();
                ss.str(""); ss << "output" << file;
                fl = MeshIO::open_file(ss.str(), ".vtk", "w"); file++;
                io.write_vtk(fl);
                io.write_vtk_data(fl, arm, "tri_quality");
                fclose(fl);
            }
            
        }
        
        // 3 area based remeshing
        for (int i = 0  ; i < 1 ; i++)
        {
            for (int j = 0 ; j < 3 ; j++)
            {
                arm.remesh_areabased();
                ss.str(""); ss << "output" << file;
                fl = MeshIO::open_file(ss.str(), ".vtk", "w"); file++;
                io.write_vtk(fl);
                io.write_vtk_data(fl, arm, "tri_quality");
                fclose(fl);
      
            }
            
            arm.remesh_flipping();
            ss.str(""); ss << "output" << file;
            fl = MeshIO::open_file(ss.str(), ".vtk", "w"); file++;
            io.write_vtk(fl);
            io.write_vtk_data(fl, arm, "tri_quality");
            fclose(fl);
        }

        // 5 smoothing
        for (int i = 0  ; i < 2 ; i++)
        {
            arm.remesh_smooth_laplacian(1);
            ss.str(""); ss << "output" << file;
            fl = MeshIO::open_file(ss.str(), ".vtk", "w"); file++;
            io.write_vtk(fl);
            io.write_vtk_data(fl, arm, "tri_quality");
            fclose(fl);
        }

	return 0;
}

int main(int argc, char *argv[])
{

	// test1(argc, argv);
	 test3(argc, argv);
	return 0;
}
