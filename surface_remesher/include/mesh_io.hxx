/*
 * mesh_io.hxx
 *
 * A variety of functionality for reading and writing meshes in different
 * formats.
 */

#ifndef MESH_IO_HXX
#define MESH_IO_HXX

#include <string>

#include "pointer_mesh.hxx"
#include "adaptive_remesher.hxx"

// EditMesh: All functions are members of this class.
class MeshIO
{
	// pointer to the mesh that we are working on.
	Mesh &_mesh;
	bool _cell_data_header;
	bool _vert_data_header;

	void write_vtk_vert_header(FILE*) ;
	void write_vtk_cell_header(FILE*) ;

public:
	// Trivial constructor
	MeshIO(Mesh &mesh_in): _mesh(mesh_in), _cell_data_header(false), _vert_data_header(false)  {}
	Mesh& mesh() {return _mesh;}
	const Mesh& mesh() const {return _mesh;}

	// Functions to read a mesh.
	void read_auto(const std::string&);
	void read_obj(const std::string&);
	void read_off(const std::string&);


	// Functions to write a mesh
	void write_vtk(const std::string&);
	void write_vtk(FILE*);
	void write_vtk_data(FILE*, AdaptiveRemesher&, std::string objective);
	void write_vtk_data(FILE*, std::vector<int>&, std::string name);
	void write_off(const std::string&) const;
	void write_vtk_dual(const std::string&);
	void print_info(FILE *fl = stdout) const;

	// Open a file
	static FILE* open_file(const std::string& fname, const std::string& format, const std::string& mode);
};

#endif
