#include <cstdio>
#include <fstream>
#include <cstring>
#include <iostream>
#include <unordered_map>
#include <set>

#include "mesh_io.hxx"

namespace hooshi {
  
using std::unordered_map;
using std::fstream;
using std::string;
using std::vector;
using std::size_t;

#ifndef MIN
#define MIN(x,y) ( (x) < (y) ? (x) : (y) )
#endif

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
    vector<double> vcoord;
    vector<std::size_t> ftovert;

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
	
    _mesh.init(vc2, ftovert);
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
    vector<std::size_t> ftovert;
    int success;

    // Read the number of verts and faces and edges
    skip_comment(fl);
    success=fscanf(fl, "%d %d %*d", &n_verts, &n_faces); assert(success=2);
    vcoord.reserve(size_t(n_verts)*3);
    ftovert.reserve(size_t(n_faces)*3);

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
        ftovert.push_back(size_t(v0)); 	ftovert.push_back(size_t(v1)); ftovert.push_back(size_t(v2));
    }
    
    // close the file
    fclose(fl);

    // init the mesh
    _mesh.init(vcoord, ftovert);
}


void MeshIO::write_obj(const string& fname) const
{

    string final_fname;
    if( fname.find(".obj") == fname.length()-4) 
        final_fname = fname; 
    else
        final_fname = fname + ".obj";

    // open the file
    fstream fs;
    fs.open (final_fname, fstream::out);
    if(!fs.is_open())
    {
        std::cout << "Could not find " << final_fname << "." << std::endl;
        assert(0);
        throw;
    }

    // Write the mesh to file
    _mesh.write_to_obj_stream(fs);

    // close the file
    fs.close();
}

static size_t omit_inactive(vector<bool>& is_active, unordered_map<size_t, size_t>& map)
{
    typedef std::pair<size_t, size_t> Pair;         
    size_t i,j;
       
    for (i=0, j=0 ; i < is_active.size() ; i++)
    {
        if(is_active[i])
        {
            map.insert(Pair(i,j));
            j++;
        }
    }

    return j;
}


void MeshIO::write_vtk(const string& fname) const
{
    // If the mesh does not have any inactive members
    if (!_mesh._is_simplification_in_progress)
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
        fprintf(fl, "POINTS %d float\n", (int)_mesh.get_vert_size());
        for (auto it = _mesh._verts.begin() ; it != _mesh._verts.end() ; ++it)
        {
            const Eigen::Vector3d &v = *it;
            fprintf(fl, "%e %e %e \n", v.x(), v.y(), v.z());
        }
        fprintf(fl, "\n");
	
        // write the faces
        fprintf(fl, "CELLS %d %d \n", (int)_mesh.get_face_size(), (int)_mesh.get_face_size()*4);
        std::size_t verts[3];
        for (uint f = 0 ; f < _mesh.get_face_size() ; f++)
        {
            _mesh.getIndicesForFace(f, verts);
            fprintf(fl, "3 %d %d %d \n", (int)verts[0], (int)verts[1], (int)verts[2]);
        }
        fprintf(fl, "\n");
	
        // write the face types
        fprintf(fl, "CELL_TYPES %d \n", (int)_mesh.get_face_size());
        for (uint f = 0 ; f < _mesh.get_face_size() ; f++)
        {
            fprintf(fl, "5\n");
        }

        // close the file
        fclose(fl);

    }

    else
    {
        // Use an unordered map to omit inactive entities.
    	unordered_map<size_t, size_t> v_map;
        unordered_map<size_t, size_t> f_map;

        const size_t nv = omit_inactive(_mesh._is_vert_active, v_map);
        const size_t nf = omit_inactive(_mesh._is_face_active, f_map);
        assert(nv == _mesh.get_vert_size());
        assert(nf == _mesh.get_face_size());
	
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
        fprintf(fl, "POINTS %d float\n", (int)nv);
        for (size_t i=0 ; i != _mesh._verts.size() ; ++i)
        {
            if(_mesh._is_vert_active[i])
            {
                // const size_t j = v_map[i];
                const Eigen::Vector3d &v = _mesh._verts[i];
                fprintf(fl, "%e %e %e \n", v.x(), v.y(), v.z());
            }
        }
        fprintf(fl, "\n");
	
        // write the faces
        fprintf(fl, "CELLS %d %d \n", (int)nf, (int)nf*4);
        std::size_t verts[3];
        for (uint f = 0 ; f < _mesh._face_to_he.size() ; f++)
        {
            if(_mesh._is_face_active[f])
            {
                // printf("%d \n", (int)_mesh._face_to_he[f]);
                _mesh.getIndicesForFace(f, verts);
                fprintf(fl, "3 %d %d %d \n", (int)v_map[verts[0]], (int)v_map[verts[1]], (int)v_map[verts[2]]);
            }
        }
        fprintf(fl, "\n");
	
        // write the face types
        fprintf(fl, "CELL_TYPES %d \n", (int)nf);
        for (uint f = 0 ; f < _mesh.get_face_size() ; f++)
        {
            fprintf(fl, "5\n");
        }
        fprintf(fl, "\n");

        // Unique vertex numbers
        fprintf(fl, "POINT_DATA %d \n", (int)_mesh.get_vert_size());
        fprintf(fl, "SCALARS unique_id int 1 \n");
        fprintf(fl, "LOOKUP_TABLE default \n");
        for (size_t i=0 ; i != _mesh._verts.size() ; ++i)
        {
            if(_mesh._is_vert_active[i])  fprintf(fl, "%d \n", (int)i);
        }
        fprintf(fl, "\n");
                    
        // write the priorities
        if(_mesh._is_simplification_in_progress == 1)
        {
            //  Method1: write the data for each vertex
            // fprintf(fl, "SCALARS priority float 1 \n");
            // fprintf(fl, "LOOKUP_TABLE default \n");
            // for (size_t i=0 ; i != _mesh._verts.size() ; ++i)
            // {
            //     if(_mesh._is_vert_active[i])
            //     {
            //         auto it = _mesh._priadd_vd[i];
            //         if(it == _mesh._prival_vd.end()) fprintf(fl, "-1000 \n");
            //         else fprintf(fl, "%f \n", (float)it->first);
            //     }
            // }
            // fprintf(fl, "\n");

            // Method2: write only the 10 biggest ones
            fprintf(fl, "CELL_DATA %d \n", (int)_mesh.get_face_size());
            fprintf(fl, "SCALARS priority int 1 \n");
            fprintf(fl, "LOOKUP_TABLE default \n");
            
            std::unordered_map<size_t, int> firstten;
            int value = 10;
            for (auto it=_mesh._prival_vd.begin(); it!=_mesh._prival_vd.end() ; ++it)
            {
                // Add the value to all faces
                vface_iterator vfiter;
                _mesh.init_iterator(vfiter, it->second);
                do
                {
                    firstten.insert(std::make_pair(_mesh.deref_iterator(vfiter), value));
                } while(_mesh.advance_iterator(vfiter));

               
                // if(value==10) printf("first vert is %d \n", (int)it->second);
                value--;
                if(value == -1 ) break;
            }
            
            for (size_t i=0 ; i != _mesh._face_to_he.size() ; ++i)
            {
                auto fans = firstten.find(i);
                if(fans != firstten.end())  fprintf(fl, "%d \n", fans->second );
                else if(_mesh._is_face_active[i]) fprintf(fl, "-1 \n" );               
            }

        }
        else if(_mesh._is_simplification_in_progress == 2)
        {
            //  Method1: write the data for each vertex
            // fprintf(fl, "SCALARS vertex_error float 1 \n");
            // fprintf(fl, "LOOKUP_TABLE default \n");
            // for (size_t i=0 ; i != _mesh._verts.size() ; ++i)
            // {
            //     if(_mesh._is_vert_active[i])
            //     {
            //         Eigen::Vector4d xyz(_mesh._verts[i].x(), _mesh._verts[i].y(), _mesh._verts[i].z(), 1);
            //         const double deviation = (xyz.transpose() * _mesh._vertexQ[i] * xyz)(0);
            //         fprintf(fl, "%f \n", (float)deviation);
            //     }
            // }
            // fprintf(fl, "\n");

            // Method2: write only the 10 biggest ones
            fprintf(fl, "CELL_DATA %d \n", (int)_mesh.get_face_size());
            fprintf(fl, "SCALARS priority int 1 \n");
            fprintf(fl, "LOOKUP_TABLE default \n");
            
            std::unordered_map<size_t, int> firstten;
            int value = 10;
            for (auto it=_mesh._prival_ec.begin(); it!=_mesh._prival_ec.end() ; ++it)
            {
                EditMesh::PriorityEC *p = *it;

                if( (p->he[0] == HOLE_INDEX) || (p->he[1] == HOLE_INDEX) ) continue;
                
                // Add the value to all faces
                vface_iterator vfiter;
                _mesh.init_iterator(vfiter, _mesh._he_data[p->he[0]].vert);
                do
                {
                    firstten.insert(std::make_pair(_mesh.deref_iterator(vfiter), value));
                } while(_mesh.advance_iterator(vfiter));
                _mesh.init_iterator(vfiter, _mesh._he_data[p->he[1]].vert);
                do
                {
                    firstten.insert(std::make_pair(_mesh.deref_iterator(vfiter), value));
                } while(_mesh.advance_iterator(vfiter));

               
                // if(value==10) printf("first vert is %d \n", (int)it->second);
                value--;
                if(value == -1 ) break;
            }
            
            for (size_t i=0 ; i != _mesh._face_to_he.size() ; ++i)
            {
                auto fans = firstten.find(i);
                if(fans != firstten.end())  fprintf(fl, "%d \n", fans->second );
                else if(_mesh._is_face_active[i]) fprintf(fl, "-1 \n" );               
            }

        }
    
        
        // close the file
        fclose(fl);

    }

}

} // End of hooshi
