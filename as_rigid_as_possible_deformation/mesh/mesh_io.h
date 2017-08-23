/*
 * mesh_io.hpp
 *
 * A variety of functionality for reading and writing meshes in different
 * formats.
 */

#ifndef MESH_IO_HPP
#define MESH_IO_HPP

#include <string>

#include "mesh.h"

// EditMesh: All functions are members of this class.
class MeshIO
{
	// pointer to the mesh that we are working on.
	Mesh *_mesh;

public:
	// Trivial constructor
	MeshIO(Mesh *mesh_in): _mesh(mesh_in)  {}

	// Functions to read a mesh.
	void read_auto(const std::string&);
	void read_obj(const std::string&);
	void read_off(const std::string&);

	// Functions to write a mesh
	void write_vtk(const std::string&) const;
	//void write_obj(const std::string&) const;
	void printf_info();
};

#endif
