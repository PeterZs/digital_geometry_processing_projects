/*
 * make_sphere.cxx
 *
 * Reads the reference sphere mesh and modifies it.
 */

#include <string>
#include "GetPot"

#include "adaptive_remesher.hxx"
#include "mesh_io.hxx"

using namespace std;

int main(int argc, char *argv[])
{
	GetPot getpot(argc, argv);

	string iname = getpot.follow("unknown.off", "-i");
	string oname = getpot.follow("unknown.off", "-o");

	Mesh mesh;
	MeshIO io(mesh);
	AdaptiveRemesher arm(mesh, getpot);
	const SphereProjector *proj = dynamic_cast<SphereProjector const*>(arm.proj());

	cout << "scale: " << proj->scale.t() ;
	io.read_auto(iname);

	/*
	 * FInd the center of the sphere
	 */
	Vertex *v;

	Vector3 cent(arma::fill::zeros);
	for (int iv = 0 ; iv < mesh.verts.n_mems() ; iv++)
	{
		v = mesh.verts.get_member_ptr(iv);
		cent += v->xyz;
	}
	cent /= mesh.verts.n_mems();
	cout << "center: " << cent.t() ;

	v = mesh.verts.get_member_ptr(0);
	const double rad = arma::norm(v->xyz - cent);
	cout << "radius: " << rad << endl;

	for (int iv = 0 ; iv < mesh.verts.n_mems() ; iv++)
	{
		v = mesh.verts.get_member_ptr(iv);
		v->xyz = (v->xyz - cent) / rad;
		v->xyz = proj->trans * v->xyz ;
	}

	// write the mesh and the normals
	FILE *fl;
	fl = MeshIO::open_file(oname, ".vtk", "w");
	io.write_vtk(fl);
	io.write_vtk_data(fl, arm, "vert_normals");
	fclose(fl);
	io.write_off(oname);

	return 0;
}
