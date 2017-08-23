#include <cstdio>
#include <fstream>
#include <cstring>
#include <iostream>
#include <unordered_map>
#include <set>

#include "mesh_io.hxx"

using std::unordered_map;
using std::fstream;
using std::string;
using std::vector;

#ifndef MIN
#define MIN(x,y) ( (x) < (y) ? (x) : (y) )
#endif

FILE* MeshIO::open_file(const string& fname, const string& format, const string& mode)
{
    // open the file
    FILE *fl;
    string final_fname;
	
    if( fname.find(format) == fname.length()-format.length()) 
        final_fname = fname; 
    else
        final_fname = fname + format;
	
    fl = fopen(final_fname.c_str(), mode.c_str());
    
    if(!fl)
    {
        std::cout << "[Error] Could not find the file " << final_fname << "." << std::endl;
        throw "Bad input file";
    }

    return fl;
}

void MeshIO::read_auto(const string& fname)
{
    if( fname.find(".obj") == fname.length()-4) this->read_obj(fname); 
    else if( fname.find(".off") == fname.length()-4) this->read_off(fname);
    else
    {
        std::cout << "Unrecognized format " << fname  << std::endl;
        throw "Unsupported input file format";
    }
}

void MeshIO::read_obj(const string& fname)
{
    // open the file
    FILE *fl = open_file(fname, ".obj", "r");
	
    /*
     * read the file
     */
    vector<double> vcoord;
    vector<int> ftovert;

    // Read each line
    int face_data_mode = -1;
    for(;;)
    {

        int c;
        c = fgetc(fl);
    
        if(c == EOF )
        {
            break;
        }
        else if (c == (int)'v')
        {
            c = fgetc(fl);
            if (c == (int)' ')
            {
                double x,y,z;
                const int success=fscanf(fl, "%lf %lf %lf", &x, &y, &z);assert(success==3);
                vcoord.push_back(x);
                vcoord.push_back(y);
                vcoord.push_back(z);
            }
        }
        else if (c == (int) 'f')
        {

            // read the face data in a single line
            char line[1000];		
            char * success = fgets(line, 1000, fl); assert(success == line);

            // find how the verts are written in a line
            // if it is not already found.
            if( face_data_mode < 0)
            {
                // read the first word in the string
                char word[100];
                sscanf(line, "%s", word);
                //printf("word: %s \n", word);
                const int n_char = strlen(word);
                int n_slash = 0;

                // count the number of slashes in there
                // and decide the format based on that.
                for (int i = 0 ; i < n_char ; i++)
                {
                    if(word[i] == '/')
                    {
                        if(word[i+1] == '/')
                        {
                            face_data_mode = 3;
                            break;
                        }
                        else
                        {
                            n_slash++;
                        }
                    }
                }

                if(n_slash == 2) face_data_mode = 4;
                else if(n_slash == 1) face_data_mode = 2;
                else if(face_data_mode!=3) face_data_mode = 1;

                std::cout << "MeshIO::read_obj: mode= " << face_data_mode << std::endl;
            }
    
		
            // read the vertex numbers
            int n_read, v0, v1, v2;
            switch(face_data_mode)
            {
            case 1:
                n_read = sscanf(line, "%d %d %d", &v0, &v1, &v2);
                break;
            case 2:
                n_read = sscanf(line, "%d%*c%*d %d%*c%*d %d%*c%*d", &v0, &v1, &v2);
                break;
            case 3:
                n_read = sscanf(line, "%d%*c%*c%*d %d%*c%*c%*d %d%*c%*c%*d ", &v0, &v1, &v2);
                break;
            case 4:
                n_read = sscanf(line, "%d%*c%*d%*c%*d %d%*c%*d%*c%*d %d%*c%*d%*c%*d", &v0, &v1, &v2);
                break;
            default:
                break;
            }
            if(n_read!=3) printf("mode: %d , line: %s, n_read: %d\n", face_data_mode, line, n_read);
            // printf("%d %d %d \n" , v0, v1, v2);
            assert(n_read==3);
			
            // shove the data in the corespoding array
            ftovert.push_back(v0-1);
            ftovert.push_back(v1-1);
            ftovert.push_back(v2-1);
        }
        else if ( (c == (int)'\n') || (c == (int)' ') )
        {
            continue;
        }
        else // if ( c == (int)'#' )
        {
            char line[1000];
            char *success = fgets(line, 1000, fl); assert(success == line);
        }
    }
	
    // close the file
    fclose(fl);

    // Remove the unused vertices
    const int n_uverts = vcoord.size();
    std::vector<bool> is_vert_used(n_uverts, false);
	
    for(int face=0 ; face < int(ftovert.size()) ; face++)
    {
        is_vert_used[ftovert[face]]=true;
    }

    std::vector<int> uvert_to_vert(n_uverts);
    int n_verts = 0;
    for(int uv = 0 ; uv < n_uverts ; uv++)
    {
        if(is_vert_used[uv])
        {
            uvert_to_vert[uv] = n_verts;
            n_verts++;
        }
        else
        {
            uvert_to_vert[uv] = -1;
        }
    }

    for(int i=0 ; i < int(ftovert.size()) ; i++)
    {
        ftovert[i]=uvert_to_vert[ftovert[i]];
    }

    std::vector<double> vc2(n_verts*3, -1);
    for(int i= 0 ; i < n_uverts ; i++)
    {
        if(is_vert_used[i])
        {
            vc2[uvert_to_vert[i]*3 + 0] = vcoord[i*3 + 0];
            vc2[uvert_to_vert[i]*3 + 1] = vcoord[i*3 + 1];
            vc2[uvert_to_vert[i]*3 + 2] = vcoord[i*3 + 2];
        }
    }
	
    mesh().init(vc2, ftovert);
}

static void skip_comment(FILE *fl)
{

    assert(fl);
	
    for(;;)
    {
        int c;
        c = fgetc(fl);
    
        if(c == EOF )
        {
            assert(0 && "Unexpected EOF");
            break;
        }
        else if ( ( (c >= (int)'0') && (c <= (int)'9') ) ||
                  (c == (int)'.')  ||
                  (c == (int)'-')	)
        {
            int success=fseek(fl, -1, SEEK_CUR); assert(success==0);
            break;
        }
    }
}

void MeshIO::read_off(const std::string& fname)
{

    // open the file
    FILE *fl = open_file(fname, ".off", "r");

    /*
     * read the file
     */
    int n_verts, n_faces;
    vector<double> vcoord;
    vector<int> ftovert;
    int success;

    // Read the number of verts and faces and edges
    skip_comment(fl);
    success=fscanf(fl, "%d %d %*d", &n_verts, &n_faces); assert(success==2);
    vcoord.reserve(n_verts*3);
    ftovert.reserve(n_faces*3);

    // Read vertex coordinates
    for (int v = 0 ; v < n_verts ; v++)
    {
        skip_comment(fl);
        double x,y,z;
        success=fscanf(fl, "%lf %lf %lf", &x, &y, &z); assert(success==3);
        vcoord.push_back(x); vcoord.push_back(y); vcoord.push_back(z);
    }

    // Read face vertices
    for (int f = 0 ; f < n_faces ; f++)
    {
        skip_comment(fl);
        int nv, v0,v1,v2;
        success=fscanf(fl, "%d %d %d %d", &nv, &v0, &v1, &v2); assert(success==4);
        assert(nv==3 && "Only triangular faces are supported");
        ftovert.push_back(v0);	ftovert.push_back(v1); ftovert.push_back(v2);
    }
    
    // close the file
    fclose(fl);

    // init the mesh
    mesh().init(vcoord, ftovert);
}


void MeshIO::write_off(const string& fname) const
{
	// open the file
	FILE *fl = open_file(fname, ".off", "w");

	// write the number of verts and faces and edges
	fprintf(fl, "%d %d %d", mesh().verts.n_mems(), mesh().faces.n_mems(), 0);

	// Read vertex coordinates
	for (int v = 0 ; v < mesh().verts.n_mems() ; v++)
	{
		Vertex *vert = mesh().verts.get_member_ptr(v);
		fprintf(fl, "%lf %lf %lf \n", vert->xyz(0), vert->xyz(1), vert->xyz(2));
	}
	fprintf(fl, "\n");

	// Read face vertices
	for (int f = 0 ; f < mesh().faces.n_mems() ; f++)
	{
		Face *face = mesh().faces.get_member_ptr(f);
		const int v1 = face->he->vert->idx();
		const int v2 = face->he->next->vert->idx();
		const int v3 = face->he->next->next->vert->idx();
		fprintf(fl, "%d %d %d %d \n", 3, v1, v2, v3);
	}
	fprintf(fl, "\n");

	// close the file
	fclose(fl);
}

void MeshIO::write_vtk(const string& fname)
{
	FILE *fl = open_file(fname, ".vtk", "w");
	this->write_vtk(fl);
	fclose(fl);
}

void MeshIO::write_vtk(FILE* fl)
{
	assert(fl);

	/*
	 * Write the vtk file.
	 */

	// write the header
	fprintf(fl, "# vtk DataFile Version 2.0\n");
	fprintf(fl, "Cartel output mesh\n");
	fprintf(fl, "ASCII\n");
	fprintf(fl, "DATASET UNSTRUCTURED_GRID\n");
	fprintf(fl, "\n");

	// write the vertices
	fprintf(fl, "POINTS %d float\n", mesh().verts.n_mems());
	for (int vidx = 0 ; vidx < mesh().verts.n_mems() ; vidx++)
	{
		const Vector3 &xyz = mesh().verts.get_member_ptr(vidx)->xyz;
		fprintf(fl, "%e %e %e \n", xyz(0), xyz(1), xyz(2));
	}
	fprintf(fl, "\n");

	// write the faces
	fprintf(fl, "CELLS %d %d \n", mesh().faces.n_mems(), mesh().faces.n_mems()*4);
	for (int f = 0 ; f < mesh().faces.n_mems() ; f++)
	{
		Face *face = mesh().faces.get_member_ptr(f);
		HalfEdge *he = face->he;
		int verts[] =
		{
				mesh().verts.get_member_idx(he->vert),
				mesh().verts.get_member_idx(he->next->vert),
				mesh().verts.get_member_idx(he->next->next->vert),
		};
		fprintf(fl, "3 %d %d %d \n", verts[0], verts[1], verts[2]);
	}
	fprintf(fl, "\n");

	// write the face types
	fprintf(fl, "CELL_TYPES %d \n", mesh().faces.n_mems());
	for (int f = 0 ; f < mesh().faces.n_mems() ; f++)
	{
		fprintf(fl, "5\n");
	}
	fprintf(fl, "\n");

	_cell_data_header = _vert_data_header = false;
}

void MeshIO::write_vtk_vert_header(FILE* fl)
{
	if(!_vert_data_header)
	{
		fprintf(fl, "POINT_DATA %d \n", mesh().verts.n_mems());
		_vert_data_header = true;
	}
}

void MeshIO::write_vtk_cell_header(FILE* fl)
{
	if(!_cell_data_header)
	{
		fprintf(fl, "CELL_DATA %d \n", mesh().faces.n_mems());
		_cell_data_header = true;
	}
}

void MeshIO::write_vtk_data(FILE* fl, AdaptiveRemesher& arm, std::string objective)
{
	if(arm.mesh() != &this->mesh())
	{
		printf("[Error] write_vtk_data needs a remesher defined for the same mesh! \n");
		printf("arm.mesh: %p, io.mesh: %p \n", static_cast<void*>(arm.mesh()), static_cast<void*>(&this->mesh()));
		throw;
	}

	// Write the face fidelity for each face
	if(objective == "vert_normals")
	{
		this->write_vtk_vert_header(fl);
		fprintf(fl, "NORMALS vert_normals float\n");

		Vector3 normal;
		for(int i = 0 ; i < mesh().verts.n_mems() ; i++)
		{
			Vertex *vert = mesh().verts.get_member_ptr(i);
			fprintf(fl, "%lf %lf %lf \n", double(vert->normal(0)), double(vert->normal(1)), double(vert->normal(2)));
		}
	}
	else if(objective == "tri_quality")
	{
		this->write_vtk_cell_header(fl);
		fprintf(fl, "SCALARS tri_quality float 1 \n");
		fprintf(fl, "LOOKUP_TABLE default \n");

		for(int i = 0 ; i < mesh().faces.n_mems() ; i++)
		{
			Face *face = mesh().faces.get_member_ptr(i);
			const double error = arm.tri_quality(face);
			fprintf(fl, "%lf \n", error);
		}
	}
	else if(objective == "tri_fidelity")
	{
		this->write_vtk_cell_header(fl);
		fprintf(fl, "SCALARS tri_fidelity float 3 \n");
		fprintf(fl, "LOOKUP_TABLE default \n");

		for(int i = 0 ; i < mesh().faces.n_mems() ; i++)
		{
			Face *face = mesh().faces.get_member_ptr(i);
			Vertex *verts[3] = {face->he->vert, face->he->next->vert, face->he->next->next->vert};
			double Etri=0. , Evert = 1.;
			bool is_okay = arm.tri_fidelity(verts, Etri, Evert);
			fprintf(fl, "%lf %lf %.0lf \n", Etri, Evert, (is_okay ? 1. : 0.) );
		}
	}
	else if(objective == "tri_obtuse")
	{
		this->write_vtk_cell_header(fl);
		fprintf(fl, "SCALARS tri_obtuse int 1 \n");
		fprintf(fl, "LOOKUP_TABLE default \n");

		for(int i = 0 ; i < mesh().faces.n_mems() ; i++)
		{
			Face *face = mesh().faces.get_member_ptr(i);
			int obtuse = 0;
			if ( arm.is_obtuse(face->he) ) obtuse = 1;
			if ( arm.is_obtuse(face->he->next) ) obtuse = 1;
			if ( arm.is_obtuse(face->he->next->next) ) obtuse = 1;
			fprintf(fl, "%d \n", obtuse );
		}
	}
	else
	{
		printf("[Error] write_vtk_data:: objective %s is invalid. \n", objective.c_str());
		throw;
	}
}

void MeshIO::write_vtk_data(FILE* fl, std::vector<int>& data, std::string name)
{
	this->write_vtk_cell_header(fl);
	fprintf(fl, "SCALARS tri_%s int 1 \n", name.c_str());
	fprintf(fl, "LOOKUP_TABLE default \n");

	for(int i = 0 ; i < mesh().faces.n_mems() ; i++)
	{
		fprintf(fl, "%d \n", data[i] );
	}

}

void MeshIO::write_vtk_dual(const std::string& fname)
{
	FILE *fl = open_file(fname, ".vtk", "w");

	/*
	 * Create the boundary info
	 * Vertex numbering would be
	 * BdryVerts - Edges - Cells
	 */
	IndexedList<BdryVertex> boundary_verts;
	IndexedList<BdryEdge>   boundary_edges;
	mesh().get_boundary_info(boundary_verts, boundary_edges, 0);

	/*
	 * Write the vtk header.
	 */
	fprintf(fl, "# vtk DataFile Version 2.0\n");
	fprintf(fl, "Cartel output dual mesh\n");
	fprintf(fl, "ASCII\n");
	fprintf(fl, "DATASET UNSTRUCTURED_GRID\n");
	fprintf(fl, "\n");

	/*
	 * Write the vertices.
	 * BdryVerts - Edges - Cells
	 */
	const int n_new_verts = boundary_verts.n_mems() + mesh().edges.n_mems() +  mesh().faces.n_mems();
	const int n_new_faces = mesh().verts.n_mems();
	const int n_bverts = boundary_verts.n_mems();
	const int n_edges =  mesh().edges.n_mems();
	// const int n_faces =  mesh().faces.n_mems();

	fprintf(fl, "POINTS %d float\n", n_new_verts);

	for (int i = 0 ; i < boundary_verts.n_mems() ; i++)
	{
		const Vector3 &xyz = boundary_verts.get_member_ptr(i)->vertex->xyz;
		fprintf(fl, "%e %e %e \n", xyz(0), xyz(1), xyz(2));
	}

	for (int i = 0 ; i < mesh().edges.n_mems() ; i++)
	{
		Edge *edge =  mesh().edges.get_member_ptr(i);
		Vector3 xyz = (edge->he->vert->xyz + edge->he->twin->vert->xyz) / 2.;
		fprintf(fl, "%e %e %e \n", xyz(0), xyz(1), xyz(2));
	}

	for (int i = 0 ; i < mesh().faces.n_mems() ; i++)
	{
		HalfEdge *he = mesh().faces.get_member_ptr(i)->he;
		Vector3 xyz = (he->vert->xyz + he->next->vert->xyz + he->next->next->vert->xyz)/3.;
		fprintf(fl, "%e %e %e \n", xyz(0), xyz(1), xyz(2));
	}

	fprintf(fl, "\n");

	/*
	 * Create the dual cells
	 */
	std::vector< std::vector<int> > dual_points(n_new_faces);
	int size_conn_list;

	for (int vidx = 0 ; vidx < mesh().verts.n_mems() ; vidx++)
	{
		Vertex *vert = mesh().verts.get_member_ptr(vidx);
		Mesh::ring_iterator ringit;
		ringit.init(vert);

		// Inner vertex
		if(!ringit.reset_boundary())
		{
			do
			{
				Face *facep = ringit.half_edge()->twin->face;
				Edge *edge = ringit.half_edge()->edge;
				dual_points[vidx].push_back( n_bverts + n_edges + facep->idx());
				dual_points[vidx].push_back(  n_bverts + edge->idx() );

			}while(ringit.advance());
		} // Done with internal vertex
		// Boundary vertex
		else
		{
			// First boundary edge
			Face *face_beg = ringit.half_edge()->twin->face;
			Edge *edge_beg = ringit.half_edge()->edge;
			dual_points[vidx].push_back( n_bverts + n_edges + face_beg->idx());
			dual_points[vidx].push_back(  n_bverts + edge_beg->idx() );
			ringit.advance();

			// Second boundary edge
			dual_points[vidx].push_back( static_cast<BdryVertex*>(vert->misc[0])->idx() );
			dual_points[vidx].push_back(  n_bverts + ringit.half_edge()->edge->idx() );
			ringit.advance();

			// Other edges
			do
			{
				Face *facep = ringit.half_edge()->twin->face;
				Edge *edge = ringit.half_edge()->edge;
				dual_points[vidx].push_back( n_bverts + n_edges + facep->idx());
				dual_points[vidx].push_back(  n_bverts + edge->idx() );

			}while(ringit.advance());

		} // Done with bdry vertex
		size_conn_list += dual_points[vidx].size();
	}

	/*
	 * Write the dual cells.
	 * VTK polygon 7
	 */
	fprintf(fl, "CELLS %d %d \n", n_new_faces, n_new_faces + size_conn_list);
	for (auto it = dual_points.begin() ;it != dual_points.end() ; ++it)
	{
		fprintf(fl, "%d ", int(it->size()));
		for (auto it2 = it->begin() ; it2 != it->end() ; ++it2 )
			fprintf(fl, "%d ", *it2);
		fprintf(fl, "\n");
	}
	fprintf(fl, "\n");

	// write the face types
	fprintf(fl, "CELL_TYPES %d \n", n_new_faces);
	for (int f = 0 ; f < n_new_faces ; f++)
	{
		fprintf(fl, "7\n");
	}
	fprintf(fl, "\n");


	// close the file
	fclose(fl);
}



void MeshIO::print_info(FILE *fl) const
{
    assert(fl);

    fprintf(fl , "%8s %8s %8s %8s %8s %8s %8s \n","#HE", "vbegin", "vend", "face", "twin", "next", "edge"); ;
    for (int i=0 ; i < mesh().half_edges.n_mems(); i++)
    {
    	HalfEdge *he = mesh().half_edges.get_member_ptr(i);
        const int twin = he->twin->idx();
        const int vbeg =  he->vert->idx();
        const int vend = he->twin->vert->idx();
        const int face = he->face->idx();
        const int next = he->next->idx();
        const int edge = he->edge->idx();

        fprintf(fl , "%8d %8d %8d %8d %8d %8d %8d(%8d)\n", i, vbeg, vend, face, twin, next, edge, he->edge->he->idx());
    }

}

