/*
 * exec_test1.cxx
 *
 * A bunch of tests. Also my playground!
 */

#include <iostream>
#include <cstdio>
#include <sstream>

#include "edit_mesh.hxx"
#include "mesh_io.hxx"

using namespace std;
using namespace hooshi;

void make_mesh1(EditMesh &m1)
{
      
    vector<size_t> conn=
    	{
            0, 1, 6,
            1, 2, 6,
            2, 3, 6,
            6, 3, 4,
            6, 4, 5,
            6, 5, 0,
    	};

    const double ang = M_PI / 3.;
    vector<double> xyz   =
    	{
            cos(0), sin(0), 0,
            cos(ang), sin(ang), 0,
            cos(2*ang), sin(2*ang), 0,
            cos(3*ang), sin(3*ang), 0,
            cos(4*ang), sin(4*ang), 0,
            cos(5*ang), sin(5*ang), 0,
            0, 0, 0,
    	};	


        m1.init(xyz, conn);

}

void make_mesh2(EditMesh &m1)
{
      
    vector<size_t> conn=
    	{
            0, 1, 6,
            1, 2, 6,
            2, 3, 6,
            6, 3, 4,
            6, 4, 5,
            6, 5, 0,

            7,8,1,
            8,9,2,
            9,10,3,
            10,11,4,
            5,11,12,
            7,0,12,

            7,1,0,
            8,2,1,
            2,9,3,
            3,10,4,
            4,11,5,
            5,12,0
    	};

    const double ang = M_PI / 3.;
    vector<double> xyz   =
    	{
            cos(0), sin(0), 0,
            cos(ang), sin(ang), 0,
            cos(2*ang), sin(2*ang), 0,
            cos(3*ang), sin(3*ang), 0,
            cos(4*ang), sin(4*ang), 0,
            cos(5*ang), sin(5*ang), 0,
            0, 0, 0,
            2*cos(ang/2.), 2*sin(ang/2.), 0,
            2*cos(1.5*ang),2* sin(1.5*ang), 0,
            2*cos(2.5*ang),2* sin(2.5*ang), 0,
            2*cos(3.5*ang),2* sin(3.5*ang), 0,
            2*cos(4.5*ang),2* sin(4.5*ang), 0,
            2*cos(5.5*ang),2* sin(5.5*ang), 0            
    	};	


        m1.init(xyz, conn);

}

void test1()
{

    EditMesh m1;

    vector<size_t> conn=
    	{
    		0, 1, 6,
    		1, 2, 6,
    		2, 3, 6,
    		6, 3, 4,
    		6, 4, 5,
    		6, 5, 0,
    	};

    const double ang = M_PI / 3.;
    vector<double> xyz   =
    	{
    		cos(0), sin(0), 0,
    		cos(ang), sin(ang), 0,
    		cos(2*ang), sin(2*ang), 0,
    		cos(3*ang), sin(3*ang), 0,
    		cos(4*ang), sin(4*ang), 0,
    		cos(5*ang), sin(5*ang), 0,
    		0, 0, 0	
    	};	
    m1.init(xyz, conn);
    //MeshIO(m1).read_auto("pyramid_0.obj");
	
    FILE *t1 = fopen("t1.txt", "w");
    FILE *t2 = fopen("t2.txt", "w");
	
    m1.init_simplification(1);
    m1.print_info(t1);
    MeshIO(m1).write_vtk("test1_0.vtk");
	
    m1.simplify_by_removing_vertex(6, {0,1,2,0,2,3,0,3,5,5,3,4});
    //m1.simplify_by_removing_vertex(4, {3,0,1,3,1,2});
    m1.verify();
    MeshIO(m1).write_vtk("test1_1.vtk");
	
    m1.restore_last_simplification_step();
    m1.print_info(t2);
    m1.verify();
    MeshIO(m1).write_vtk("test1_2.vtk");

    fclose(t1);
    fclose(t2);
}


void test2()
{
    EditMesh m1;
    //MeshIO(m1).read_auto("../meshes/p1_inputs/closed/cube.obj");
    //MeshIO(m1).read_auto("../meshes/p1_inputs/closed/sphere.obj");
    //MeshIO(m1).read_auto("../meshes/p1_inputs/closed/camel_mc.off");
    //MeshIO(m1).read_auto("../Cartel/Mesh/cow1.obj");
    //MeshIO(m1).read_auto("../Cartel/Mesh/cow1.obj");    
    MeshIO(m1).read_auto("../../meshes/p1_inputs/closed/sphere.obj");
    //MeshIO(m1).read_auto("../meshes/test.obj");
    // MeshIO(m1).read_auto("../meshes/nicolo/439.off");    
  
    // vector<size_t> conn=
    // 	{
    //         0, 1, 6,
    //         1, 2, 6,
    //         2, 3, 6,
    //         6, 3, 4,
    //         6, 4, 5,
    //         6, 5, 0,
    // 	};

    // const double ang = M_PI / 3.;
    // vector<double> xyz   =
    // 	{
    //         cos(0), sin(0), 0,
    //         cos(ang), sin(ang), 0,
    //         cos(2*ang), sin(2*ang), 0,
    //         cos(3*ang), sin(3*ang), 0,
    //         cos(4*ang), sin(4*ang), 0,
    //         cos(5*ang), sin(5*ang), 0,
    //         0, 0, 0	
    // 	};	
    // m1.init(xyz, conn);
    

	
    m1.init_simplification(1);
    MeshIO(m1).write_vtk("test2_0.vtk");

    // std::vector<size_t> triverts;
    const int nstep = 1000;
    stringstream ss;
    int i;
    for (i = 0 ; i < nstep ; i++)
    {
        //  m1.analyze_vertex_for_removal(m1.get_vert_size()-1, triverts);
        // m1.simplify_by_removing_vertex(m1.get_vert_size()-1, triverts);
        bool res = m1.simplify(); if(!res){ printf("could not do further\n");  break; assert(res);}
        m1.verify();
        ss.str(""); ss << "test2_" << i+1 << ".vtk";
        if(i>900 ) MeshIO(m1).write_vtk(ss.str());
    }
    ss.str(""); ss << "test2_" << i+1 << ".vtk";
    MeshIO(m1).write_vtk(ss.str());

    // for (int i = 0 ; i < nstep ; i++)
    // {
    //     m1.restore_last_simplification_step();
    //     m1.verify();
    //     ss.str(""); ss << "test2_" << i+nstep+1 << ".vtk";
    //     if(i%10 == 0 ) MeshIO(m1).write_vtk(ss.str()); 
    // }
  
}

void test3(uint method, const char* name)
{
    EditMesh m1;
    //MeshIO(m1).read_auto("../meshes/p1_inputs/closed/cube.obj");
    // MeshIO(m1).read_auto("../../meshes/p1_inputs/closed/sphere.obj");
    //MeshIO(m1).read_auto("../../meshes/p1_inputs/closed/camel_mc.off");
    //MeshIO(m1).read_auto("../../meshes/cartel/cow2.obj");
//    MeshIO(m1).read_auto("../../meshes/cartel/camel.obj");
//        MeshIO(m1).read_auto("../../meshes/p1_inputs/closed/sphere.obj");

    //MeshIO(m1).read_auto("../meshes/p1_inputs/closed/bunny_mc.off");
    //MeshIO(m1).read_auto("../meshes/test.obj");
    //MeshIO(m1).read_auto("../../meshes/raptor/178_raptor.off");
    //MeshIO(m1).read_auto("../../meshes/cartel/horse.obj");
    
    //make_mesh2(m1);
    MeshIO(m1).read_auto(name);
    


    stringstream ss;
    m1.init_simplification(method);
    ss << "test3_m" << method << "_" << 0 << ".vtk";
    MeshIO(m1).write_vtk(ss.str());

    //m1.simplify_by_collapsing_edge(m1.find_twin(6, 2)->twin);
    // m1.simplify();
    // m1.verify();
    // MeshIO(m1).write_vtk("test3_1.vtk");
    
    // std::vector<size_t> triverts;
    const int nstep = 1e6;
    const int nwrite = 10;

    int i;
    for (i = 0 ; i < nstep ; i++)
    {
        //  m1.analyze_vertex_for_removal(m1.get_vert_size()-1, triverts);
        // m1.simplify_by_removing_vertex(m1.get_vert_size()-1, triverts);
        bool res = m1.simplify(); if(!res){ printf("could not do more than %d\n", i);  break; assert(res);}
        ss.str(""); ss << "test3_m" << method << "_" << i+1 << ".vtk";
        if((i%nwrite ==0) ) MeshIO(m1).write_vtk(ss.str());
        // else if(m1.get_face_size()==994) MeshIO(m1).write_vtk(ss.str());
        // else if(m1.get_face_size()==532) MeshIO(m1).write_vtk(ss.str());
        // else if(m1.get_face_size()==248) MeshIO(m1).write_vtk(ss.str());
        // else if(m1.get_face_size()==64) MeshIO(m1).write_vtk(ss.str());
#ifndef NDEBUG        
        m1.verify();
#endif
    }
    ss.str(""); ss << "test3_m" << method << "_" << i+1 << ".vtk";
    if((i%nwrite !=0) ) MeshIO(m1).write_vtk(ss.str());

    // if(method==1)
    // for (int i = 0 ; i < nstep ; i++)
    // {
    //     if(! m1.restore_last_simplification_step()) break;
    //     //m1.verify();
    //     ss.str(""); ss << "test3_m" << method << "_" << i+nstep+1  << ".vtk";
    //     if(i%10 == 0 ) MeshIO(m1).write_vtk(ss.str()); 
    // }

    m1.finalize_simplification();
}

int main(int argc, char *argv[])
{
    try
    {
    //test1();
    //test2();
    assert(argc==3);
    test3(stoi(argv[1]), argv[2]);
    }
    catch(const char* msg)
    {
        printf("Error: %s\n", msg);
    }
    return 0;
    
}
