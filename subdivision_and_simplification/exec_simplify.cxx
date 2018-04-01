/*
 * exec_simplify.cxx
 *
 * This executable reads a mesh, simplifies it to have a certain
 * number of vertices and then writes the answer. To demostrate the
 * ability to retrieve old meshes, there is also the option to go back
 * and write the old meshes as well.
 */


#include <sstream>

#include "edit_mesh.hxx"
#include "mesh_io.hxx"

// A Cartel always needs its Drugs. Pot in this case! :)
#include "GetPot"

using namespace std;
using namespace hooshi;

void print_help()
{
	cout << "The usage is: " << endl;
	cout << "./exec_simplify <options>" << endl;
	cout << "Where the options are (The marked ones are mandatory): " << endl;
	cout << "-i* <input file name> " << endl;
	cout << "-o <output file name> " << endl;
        cout << "-n X1 Y1 X2 Y2 ... where Xi is target #removed_vertex and Yi number of files to be written by then" << endl;
	cout << "-m <0 or 1> where 1->vertex decimation, 2->quartic edge collapse" << endl;
        cout << "-nb <# files to write as proof> show ability to restore original mesh" << endl;
	cout << "-h, --help: print this message " << endl;	
}

int main(int argc, char *argv[])
{

    try
    {
        
	// getpot reads the command line
	GetPot getpot(argc, argv);
	// The mesh
	EditMesh mesh;
        MeshIO io(mesh);
        stringstream ss;

	// Check if we should print help
	if(getpot.search(2, "--help", "-h"))
	{
		print_help();
		return 0;
	}

	// Options
	const int method = getpot.follow(1, "-m");
	const int n_backward = getpot.follow(0, "-nb");	
	if( ! getpot.search("-i") )
	{
            printf("* Input file must be specified! *\n");
            printf("* Use exec_simplify --help! *\n");
            return 0 ;
	}
	const string fname_in = getpot.next("NA");
	const string fname_out = getpot.follow("out", "-o");
        const vector<string> nout_str = getpot.nominus_followers("-n");

        // Parse the nout_str
        vector<int> nvr, nfl;
        if(nout_str.size() % 2 != 0) throw "-n does not have even size";
        for(int i = 0 ;  i < (int)nout_str.size() ; i+=2)
        {
            int tmp;
            
            ss.str(""); ss.clear(); ss<<nout_str[i]; ss>>tmp; nvr.push_back(tmp);
            if(ss.bad()) throw "-n is not well formatted";

            ss.str(""); ss.clear(); ss<<nout_str[i+1]; ss>>tmp; nfl.push_back(tmp);
            if(ss.bad()) throw "-n is not well formatted";
        }
        for (int i = 0 ; i <  (int)nfl.size() ; i++)
        {
            printf("step %d: remove total of %d verts, write %d files \n", i, nvr[i], nfl[i]);
        }

        // Read the mesh
        bool wrote_this_step =false;
        io.read_auto(fname_in);

        // Write the first step
        ss.str(""); ss.clear();
        ss<<fname_out<<"_m"<<method<<"_"<<0<<".vtk";
        io.write_vtk(ss.str());
        wrote_this_step = true;

        //Init simplification
        mesh.init_simplification(method);
        
        // Start simplification and write the results
        int isim = 0;
        bool ressim;
        for(int step=0 ; step < (int)nvr.size() ; step++)
        {
            // Figure out how many files have to be written
            if( (step>0) && (nvr[step] < nvr[step-1]) ) throw "-n arguments must be increasing";
            const int nprev  = (step == 0 ? 0 : nvr[step-1]);
            const int n_write = (nvr[step] - nprev) / nfl[step];

            // Start simplification
            for (; isim < nvr[step] ; isim++)
            {
                ressim = mesh.simplify();
                wrote_this_step = false;
                if(!ressim)
                {
                    printf("could not do more than %d removals\n", (int)isim);
                    break;
                }
                if( ( (isim-nprev) % n_write == 0 ) && (!wrote_this_step) )
                {
                    printf("Removed %d vertices, writing file ... \n", isim);
                    ss.str(""); ss.clear();
                    ss<<fname_out<<"_m"<<method<<"_"<<isim<<".vtk";
                    io.write_vtk(ss.str());
                    wrote_this_step = true; 
                }
#ifndef NDEBUG        
                mesh.verify();
#endif
            }

            // Write the last step if we have not already done it.
            if(!wrote_this_step)
            {
                ss.str(""); ss.clear();
                ss<<fname_out<<"_m"<<method<<"_"<<isim<<".vtk";
                io.write_vtk(ss.str());
                wrote_this_step = true; 
            }

            // Check if we cannot simplify further
            if(!ressim) break;
        }

        // Go back to the initial mesh and write it
        const int n_write = isim / n_backward;
        const int iinit = isim+1;
        if( (method==1) && (n_backward > 0) )
        {
            for (; isim >= 0 ; isim--)
            {
                if(! mesh.restore_last_simplification_step()) break;
                wrote_this_step = false;
                if( (iinit-isim) % n_write == 0  )
                {
                    ss.str(""); ss.clear();
                    ss<<fname_out<<"_m"<<method<<"_"<<2*iinit-isim<<".vtk";
                    io.write_vtk(ss.str());
                    wrote_this_step = true; 
                }
            }
            if(!wrote_this_step)
            {
                ss.str(""); ss.clear();
                ss<<fname_out<<"_m"<<method<<"_"<<2*iinit-isim<<".vtk";
                io.write_vtk(ss.str());
                wrote_this_step = true; 
            }
        }


        // Finalize simplification
        mesh.finalize_simplification();
    }
    catch(const char* msg)
    {
        printf("execption occured: %s \n", msg);
    }
	return 0;
}

