/*
 * arap.h
 *
 * Declaration of the Deformation class.
 */

#ifndef _ARAP_H
#define _ARAP_H

#include <map>

#include "mesh.h"

// Petsc linear solver code!
#include <petscksp.h>
// Armadillo library.
#include <armadillo>

// Basically Armadillo handles small matrices and svd.
// Petsc handles huge matrices and linear solvers.

/*
 * A wrapper for the PETSc Vector class.
 * Contains both the vector and a pointer to its data.
 */
struct PetscVecWrapper
{
	// Petsc Vector.
	Vec vec;
	// Pointer to vector data.
	double *v;
	PetscVecWrapper():vec(NULL), v(NULL){}
};

// Used functionality from the armadillo library.
typedef	arma::vec::fixed<3> vec3d;
typedef arma::mat::fixed<3,3> matrix3d;

/*
 * Deformation.
 *
 * Deformation takes a mesh as input. And deforms it using the ARAP method.
 */
class Deformation
{
private:

	// Mesh to deform
	Mesh *_mesh;
	// Rotation matrices for each vertex
	arma::field<matrix3d> _rot;
	// Original location of vertices
	arma::field<vec3d>    _xyz0;
	// Current guess for the location of the vertices.
	arma::field<vec3d>    _xyz;
	// Right hand side of the linear part of the ARAP system.
	arma::field<vec3d>    _rhs;
	// A vector to have random access to each vertex, cause mesh
	// class only stores them in a list.
	vector<Vertex*>       _vertIndex;
	// Weight of each edge.
	vector<double>        _edgew;
	// Matrix A, in the ARAP method.
	Mat                   _A;
	// Matrix A^T * A in the ARAP method.
	Mat 				  _ATA;
	// The vector to get the solution of Ax=b from petsc.
	PetscVecWrapper       _X;
	// Right hand side: b.
	PetscVecWrapper       _B;
	// A^T * b.
	PetscVecWrapper       _ATB;
	// Petsc Linear solver.
	KSP                   _ksp;

	// Vertex number and location of the anchor vertices to be
	// imposed as constraints.
	map<int, vec3d >   _anchors;

	// Becomes true if the anchor vertices have changed, so that
	// the matrix has to be assembled and refactored again.
	bool _should_refactor;

	// The weight for imposing anchor vertices.
	double _beta;
	// Tolerance for checking convergence of the nonlinear solver.
	double _tol;
	// Should we use uniform weights for edges?
	bool _should_use_uniform_wgt;
	// Maximum number of allowed nonlinear iterations.
	int _n_max_iter;
	// Level of verbosity.
	int _verbosity;
	// How many deformations we have made so far.
	int _n_total_deformations;

	// Assemble the left hand side
	PetscErrorCode assembleMat();
	// Assemble the right hand side
	PetscErrorCode assembleRHS();
	// Find the rotation matrices again
	PetscErrorCode updateRots();

public:
	// Petsc error code for error checking.
	PetscErrorCode ierr;

	// Constructor and destructor.
	Deformation(Mesh *mesh_in);
	PetscErrorCode destroy();
	~Deformation();

	// THis is just for testing.
	vector<Vertex*>*  getVertIndex(){return &_vertIndex;}

	// Reset the position of vertices in the mesh to their original.
	PetscErrorCode resetMeshPositions();

	// Deform the mesh
	// Input: anchorverts: index of anchors.
	//        newpos:      location of each anchor in order.
	PetscErrorCode deform(vector<int> anchorverts, vector<float> newpos);
};

/*
 * Functions to move data from Petsc to armadillo and vice versa.
 */

inline void mergeDim(const int dim, PetscVecWrapper in, arma::field<vec3d>& ans)
{
	assert( (dim == 0) || (dim == 1) || (dim == 2) );

	for (int iv = 0 ; iv < (int)ans.size() ; iv++)
	{
		ans(iv)(dim) = in.v[iv];
	}
}

inline void extractDim(const int dim, arma::field<vec3d>& in, PetscVecWrapper ans)
{
	assert( (dim == 0) || (dim == 1) || (dim == 2) );

	for (int iv = 0 ; iv < (int)in.size() ; iv++)
	{
		ans.v[iv] = in(iv)(dim);
	}
}


#endif /*_ARAP_H*/
