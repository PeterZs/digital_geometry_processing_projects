/*
 * mesh_io.hxx
 *
 * A variety of functionality for reading and writing meshes in different
 * formats.
 */

#ifndef HOOSHI_MESH_IO_HXX
#define HOOSHI_MESH_IO_HXX

#include <string>

#include "edit_mesh.hxx"

namespace hooshi {
  
// EditMesh: All functions are members of this class.
class MeshIO
{
	// pointer to the mesh that we are working on.
	EditMesh &_mesh;

public:
	// Trivial constructor
	MeshIO(EditMesh &mesh_in): _mesh(mesh_in)  {}

	// Functions to read a mesh.
	void read_auto(const std::string&);
	void read_obj(const std::string&);
	void read_off(const std::string&);

	// Functions to write a mesh
	void write_vtk(const std::string&) const;
	void write_obj(const std::string&) const;
};

} // End of hooshi

#endif

  
