/*
 * split.exe.cxx
 *
 * This executable reads a mesh, subdivides for a number of times and
 * writes each subdivided mesh for visualization.
 */

#include <sstream>
#include <iostream>
#include <algorithm>
#include <ctime>

#include "parameters.hxx"
#include "mesh_io.hxx"
#include "adaptive_remesher.hxx"

#include "GetPot"
#include "armadillo"

using namespace std;


#define WRITE_FILE(i, name, io, arm)                            \
    do{                                                         \
        stringstream ss; FILE *fl;                              \
        ss.str(""); ss << name << i;                            \
        fl = MeshIO::open_file(ss.str(), ".vtk", "w"); i++;     \
        io.write_vtk(fl);                                       \
        io.write_vtk_data(fl, arm, "tri_quality");              \
        io.write_vtk_data(fl, arm, "tri_fidelity");             \
        io.write_vtk_data(fl, arm, "tri_obtuse");               \
        fclose(fl);                                             \
    }while(0)

#define WRITE_SIZE(arm, mesh)                                           \
    do{                                                                 \
        cout << "[log] n_vertices: " << mesh.verts.n_mems() << ", ";    \
        PatchProjector *proj = dynamic_cast<PatchProjector*>(arm.proj()); \
        if (proj) cout << "n_patches, " << proj->patches.size();        \
        cout << endl;                                                   \
    }while(0)

#define MEASURE_TIME_BEGIN(t_init, t_total)     \
    do{                                         \
        t_init = clock();                       \
    }while(0)

#define MEASURE_TIME_END(t_init, t_total)                     \
    do{                                                         \
        clock_t t_end = clock();                                \
        t_total += (double)(t_end - t_init) / CLOCKS_PER_SEC;   \
    }while(0)

vector<double> MEASURE_QUALITY(AdaptiveRemesher& arm)
{
    vector<double> quality(6,0);
    for (int i = 0 ; i < arm.mesh()->faces.n_mems() ; i++ )
    {
        Face *face = arm.mesh()->faces.get_member_ptr(i);
        const double qq = arm.tri_quality(face);
        if (qq < 10) quality[0]++;
        else if (qq < 20) quality[1]++;
        else if (qq < 30) quality[2]++;
        else if (qq < 40) quality[3]++;
        else if (qq < 50) quality[4]++;
        else if (qq < 60) quality[5]++;
    }

    for (int j = 0 ; j < 6 ; j++)
        quality[j] = quality[j] / arm.mesh()->faces.n_mems() * 100;

    return quality;
}

/*
 * The main function
 */
int test1(int argc, char *argv[])
{
    GetPot getpot(argc,argv);
    Mesh mesh;
    MeshIO io(mesh);
    AdaptiveRemesher arm(mesh,getpot);
    Mesh::ring_iterator it;
    int file =1;
    double time_split = 0,
        time_collapse = 0,
        time_areabased = 0,
        time_laplacian = 0;
    clock_t time_init;
    vector<double> quality_init, quality_final;
    int n_verts_begin, n_verts_end;

    // Read the mesh
    string iname = getpot.follow("meshes/square.off", "-i");
    string oname = getpot.follow("output", "-o");
    int n_target_verts = getpot("n_verts", -1);
    int remove_obtuse = getpot("r_obtuse", 0);
    arm.arearemesh_mode() = AdaptiveRemesher::CONSTANT;
    io.read_auto(iname);
    quality_init = MEASURE_QUALITY(arm);

    /*
     * Initialize the projection object
     */
    arm.init();

    // Write the mesh
    n_verts_begin = mesh.verts.n_mems();
    WRITE_FILE(file, oname, io, arm);
    WRITE_SIZE(arm, mesh);

    /*
     * collapse or insertion
     */
    if (n_target_verts > 0)
    {
        bool scss;
        int mode; /* 0-> collapse, 1->insertion */

        arm.n_tverts() = n_target_verts;
        if ( mesh.verts.n_mems() > n_target_verts) mode = 0;
        else mode = 1;

        do
        {
                    
            if( mode == 0 )
            {
                MEASURE_TIME_BEGIN(time_init, time_collapse);
                scss = arm.remesh_collapseedges();
                MEASURE_TIME_END(time_init, time_collapse);
            }
            else
            {
                MEASURE_TIME_BEGIN(time_init, time_split); 
                scss = arm.remesh_splitedges();
                MEASURE_TIME_END(time_init, time_split);
            }
            WRITE_FILE(file, oname, io, arm);
            WRITE_SIZE(arm, mesh);

            MEASURE_TIME_BEGIN(time_init, time_areabased); 
            arm.remesh_areabased(1, 3);
            MEASURE_TIME_END(time_init, time_areabased);
            WRITE_FILE(file, oname, io, arm);
        }while( scss );
    }

    /*
     * Insertion
     */
    if(remove_obtuse)
    {
        MEASURE_TIME_BEGIN(time_init, time_split); 
        arm.remesh_obtuse_angles();
        MEASURE_TIME_END(time_init, time_split);
        WRITE_FILE(file, oname, io, arm);
        WRITE_SIZE(arm, mesh);
    }

    /*
     * Area based smoothing and flipping
     */
    MEASURE_TIME_BEGIN(time_init, time_areabased);
    arm.remesh_areabased(10, 3);
    MEASURE_TIME_END(time_init, time_areabased);
    WRITE_FILE(file, oname, io, arm);
    WRITE_SIZE(arm, mesh);

    /*
     * laplacian smoothing
     */
    MEASURE_TIME_BEGIN(time_init, time_laplacian);
    arm.remesh_smooth_laplacian(10);
    MEASURE_TIME_END(time_init, time_laplacian);
    WRITE_FILE(file, oname, io, arm);
    WRITE_SIZE(arm, mesh);

    /*
     * save the statistics
     */
    string statfilename = oname + ".stats";
    FILE *statsfl =  fopen(statfilename.c_str(), "w");

    // NUMBER OF VERTICES
    n_verts_end = mesh.verts.n_mems();
    fprintf(statsfl, "# N_VERTS_B: %d  \n", n_verts_begin);
    fprintf(statsfl, "# N_VERTS_E: %d  \n", n_verts_end);
    

    
    // TIME
    double time_total = time_areabased + time_collapse + time_split + time_laplacian;
    fprintf(stdout, "# INSERTION: %lf s \n", time_split);
    fprintf(stdout, "# COLLAPSE: %lf s \n", time_collapse);
    fprintf(stdout, "# AREABASED: %lf s \n", time_areabased);
    fprintf(stdout, "# LAPLACIAN: %lf s \n", time_laplacian);
    fprintf(stdout, "# TOTAL: %lf s \n", time_total);
    fprintf(statsfl, "# INSERTION: %lf s \n", time_split);
    fprintf(statsfl, "# COLLAPSE: %lf s \n", time_collapse);
    fprintf(statsfl, "# AREABASED: %lf s \n", time_areabased);
    fprintf(statsfl, "# LAPLACIAN: %lf s \n", time_laplacian);
    fprintf(statsfl, "# TOTAL: %lf s \n", time_total);

    // MESH QUALITY
    quality_final = MEASURE_QUALITY(arm);
    fprintf(statsfl, "# %10s %10s %10s %10s %10s %10s \n",
            "0 - 10", "10 - 20", "20 - 30",
            "30 - 40", "40 - 50", "50 - 60"); 
    fprintf(statsfl, "%10.2lf %10.2lf %10.2lf %10.2lf %10.2lf %10.2lf \n",
            quality_init[0], quality_init[1], quality_init[2],
            quality_init[3], quality_init[4], quality_init[5]);
    fprintf(statsfl, "%10.2lf %10.2lf %10.2lf %10.2lf %10.2lf %10.2lf \n",
            quality_final[0], quality_final[1], quality_final[2],
            quality_final[3], quality_final[4], quality_final[5]);

    fprintf(stdout, "%10s %10s %10s %10s %10s %10s \n",
            "0 - 10", "10 - 20", "20 - 30",
            "30 - 40", "40 - 50", "50 - 60"); 
    fprintf(stdout, "%10.2lf %10.2lf %10.2lf %10.2lf %10.2lf %10.2lf \n",
            quality_init[0], quality_init[1], quality_init[2],
            quality_init[3], quality_init[4], quality_init[5]);
    fprintf(stdout, "%10.2lf %10.2lf %10.2lf %10.2lf %10.2lf %10.2lf \n",
            quality_final[0], quality_final[1], quality_final[2],
            quality_final[3], quality_final[4], quality_final[5]);

    fclose(statsfl);
    return 0;
}

int main(int argc, char *argv[])
{

    test1(argc, argv);
    return 0;
}
