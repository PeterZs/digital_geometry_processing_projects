#include <cstdio>
#include <fstream>
#include <cstring>
#include <iostream>
#include <unordered_map>
#include <set>

#include "mesh_io.h"

using std::unordered_map;
using std::fstream;
using std::string;
using std::vector;
using std::size_t;

#ifndef MIN
#define MIN(x,y) ( (x) < (y) ? (x) : (y) )
#endif

const size_t HOLE_INDEX = static_cast<size_t>(-1);

static FILE* open_file(const string& fname, const string& format, const string& mode)
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
        std::cout << "Could not find the file " << final_fname << "." << std::endl;
        // Do whatever you can to stop the process !!
        assert(0);
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
        std::cout << "Unrecognized format " << fname << "." << std::endl;
        assert(0);
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
    vector<float> vcoord;
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
            ftovert.push_back(std::size_t(v0-1));
            ftovert.push_back(std::size_t(v1-1));
            ftovert.push_back(std::size_t(v2-1));
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
    const size_t n_uverts = vcoord.size();
    std::vector<bool> is_vert_used(n_uverts, false);

    for(size_t face=0 ; face < ftovert.size() ; face++)
    {
        is_vert_used[ftovert[face]]=true;
    }

    std::vector<size_t> uvert_to_vert(n_uverts);
    size_t n_verts = 0;
    for(size_t uv = 0 ; uv < n_uverts ; uv++)
    {
        if(is_vert_used[uv])
        {
            uvert_to_vert[uv] = n_verts;
            n_verts++;
        }
        else
        {
            uvert_to_vert[uv] = HOLE_INDEX;
        }
    }

    for(size_t i=0 ; i < ftovert.size() ; i++)
    {
        ftovert[i]=uvert_to_vert[ftovert[i]];
    }

    std::vector<double> vc2(n_verts*3, -1);
    for(size_t i= 0 ; i < n_uverts ; i++)
    {
        if(is_vert_used[i])
        {
            vc2[uvert_to_vert[i]*3 + 0] = vcoord[i*3 + 0];
            vc2[uvert_to_vert[i]*3 + 1] = vcoord[i*3 + 1];
            vc2[uvert_to_vert[i]*3 + 2] = vcoord[i*3 + 2];
        }
    }

    // Create the mesh
    vector<Vertex*> vertIndex;
    vertIndex.reserve(vc2.size()/3);

    for (int iv = 0 ; iv < (int)vc2.size() ; iv+=3)
    {
    	vertIndex.push_back( new Vertex((float)vc2[iv], (float)vc2[iv+1], (float)vc2[iv+2] ) );
    	_mesh->addVertex(vertIndex.back());
    }
    for(int it = 0 ; it < (int)ftovert.size() ; it+=3)
    {
    	_mesh->addTriangle(new Triangle(vertIndex[ftovert[it]], vertIndex[ftovert[it+1]], vertIndex[ftovert[it+2]]));
    }
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
    int success;
    vector<Vertex*> vertIndex;

    // Read the number of verts and faces and edges
    skip_comment(fl);
    success=fscanf(fl, "%d %d %*d", &n_verts, &n_faces); assert(success=2);
    vertIndex.reserve(n_verts);

    // Read vertex coordinates
    for (int iv = 0 ; iv < n_verts ; iv++)
    {
        skip_comment(fl);
        Vertex *v = new Vertex;
        vertIndex.push_back(v);
        success=fscanf(fl, "%f %f %f", &v->x(), &v->y(), &v->z()); assert(success==3);
        _mesh->addVertex(v);
    }

    // Read face vertices
    for (int f = 0 ; f < n_faces ; f++)
    {
        skip_comment(fl);
        int nv, v0,v1,v2;
        success=fscanf(fl, "%d %d %d %d", &nv, &v0, &v1, &v2); assert(success==4);
        assert(nv==3 && "Only triangular faces are supported");
        _mesh->addTriangle(new Triangle(vertIndex[v0], vertIndex[v1], vertIndex[v2]));
    }

    // close the file
    fclose(fl);
}


//void MeshIO::write_obj(const string& fname) const
//{
//    string final_fname;
//    if( fname.find(".obj") == fname.length()-4)
//        final_fname = fname;
//    else
//        final_fname = fname + ".obj";
//
//    // open the file
//    fstream fs;
//    fs.open (final_fname, fstream::out);
//    if(!fs.is_open())
//    {
//        std::cout << "Could not find " << final_fname << "." << std::endl;
//        assert(0);
//        throw;
//    }
//
//    // Write the mesh to file
//    _mesh.write_to_obj_stream(fs);
//
//    // close the file
//    fs.close();
//}


void MeshIO::write_vtk(const string& fname) const
{
	// open the file
	FILE *fl = open_file(fname, ".vtk", "w");

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
	fprintf(fl, "POINTS %d float\n", (int)_mesh->numberOfVertices());
	for (auto it = _mesh->getVertices()->begin() ; it != _mesh->getVertices()->end() ; ++it)
	{
		Vertex *v = *it;
		fprintf(fl, "%e %e %e \n", v->x(), v->y(), v->z());
	}
	fprintf(fl, "\n");

	// write the faces
	fprintf(fl, "CELLS %d %d \n", _mesh->numberOfTriangles(), _mesh->numberOfTriangles()*4);
	for (auto it = _mesh->getTriangles()->begin() ; it != _mesh->getTriangles()->end() ; ++it)
	{
		Triangle *tri = *it;
		fprintf(fl, "3 %d %d %d \n", tri->vertices[0]->name, tri->vertices[1]->name, tri->vertices[2]->name);
	}
	fprintf(fl, "\n");

	// write the face types
	fprintf(fl, "CELL_TYPES %d \n", (int)_mesh->numberOfTriangles());
	for (int f = 0 ; f < _mesh->numberOfTriangles() ; f++)
	{
		fprintf(fl, "5\n");
	}

	// close the file
	fclose(fl);
}

void MeshIO::printf_info()
{
	printf("n_edge: %d \n", _mesh->numberOfEdges());

	for (auto ie = _mesh->getEdges()->begin() ; ie != _mesh->getEdges()->end() ; ie++)
	{
		const Edge *e = *ie;
		printf("edge(v,v,t,t): %5d %5d %5d %5d %5d \n",
				 e->name, e->vertices[0]->name, e->vertices[1]->name, e->triangles[0]->name, e->triangles[1]->name);
	}

	for (auto it = _mesh->getTriangles()->begin() ; it != _mesh->getTriangles()->end() ; it++)
	{
		Triangle *t = *it;
		t->edges[0]->makeFirstTriangle(t);
		t->edges[1]->makeFirstTriangle(t);
		t->edges[2]->makeFirstTriangle(t);
		printf("tri: %5d %5d %5d %5d \n",
				t->name, t->edges[0]->triangles[1]->name, t->edges[1]->triangles[1]->name, t->edges[2]->triangles[1]->name);
	}
}
