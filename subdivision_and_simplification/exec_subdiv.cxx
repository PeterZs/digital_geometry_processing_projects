/*
 * exec_subdiv.cxx
 *
 * This executable reads a mesh, subdivides for a number of times and
 * writes each subdivided mesh for visualization.
 */

#include <sstream>

#include "edit_mesh.hxx"
#include "mesh_io.hxx"

// A Cartel always needs its Drugs. Pot in this case! :)
#include "GetPot"

using namespace hooshi;
using std::cout;
using std::endl;
using std::string;
using std::stringstream;

void print_help()
{
	cout << "The usage is: " << endl;
	cout << "./exec_subdiv <options>" << endl;
	cout << "Where the options are (The marked ones are mandatory): " << endl;
	cout << "-i* <input file name> " << endl;
	cout << "-o <output file name> " << endl;
	cout << "-n <number of subdivision steps> " << endl;
	cout << "-m <0 or 1 or 2> where 0->loop, 1->sqrt3, 2->butterfly" << endl;
	cout << "--vtk: also write vtk files." << endl;
	cout << "-h, --help: print this message " << endl;	

}

int main(int argc, char *argv[])
{

	// getpot reads the command line
	GetPot getpot(argc, argv);
	// The mesh
	EditMesh mesh;

	// Check if we should print help
	if(getpot.search(2, "--help", "-h"))
	{
		print_help();
		return 0;
	}

	// Options
	const uint n_subdiv = getpot.follow(1, "-n");
	const uint method = getpot.follow(0, "-m");
	const bool should_write_vtk = getpot.search("--vtk");	
	if( ! getpot.search("-i") )
	{
		print_help();
		return 0 ;
	}
	const string fname_in = getpot.next("NA");
	const string fname_out_default = fname_in.substr(0, fname_in.length()-4);
	string fname_out = getpot.follow(fname_out_default.c_str(), "-o");

	
	// Read the mesh
	MeshIO(mesh).read_auto(fname_in);
					 
	// Do the subdivions
	stringstream ss;

	// 0'th step
	ss.str("");
	ss << fname_out << "_" << 0 ;
	MeshIO(mesh).write_obj(ss.str());
	if (should_write_vtk) MeshIO(mesh).write_vtk(ss.str());
	
	for (uint i=0 ; i < n_subdiv ; i++)
	{

		std::cout << "subdividing step " << i+1 << " ... ";
		
		// subdivide
		switch(method)
		{
		case 0:
			mesh.subdivide_loop();
			break;
		case 1:
			mesh.subdivide_sqrt3();
			break;
		case 2:
			mesh.subdivide_butterfly();
			break;
		default:
			std::cout << "Invalid method " << method << "." << std::endl;
			assert(0);
			throw;
		}

		std::cout << " done." << std::endl;

		// write the mesh
		ss.str("");
		ss << fname_out << "_" << i + 1;
		MeshIO(mesh).write_obj(ss.str());
		if(should_write_vtk) MeshIO(mesh).write_vtk(ss.str());
	}


	return 0;
}
