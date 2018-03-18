/*
 * arap.cpp
 *
 * Function definitions for the Deformation class.
 */


#include "arap.h"
#include "mesh_io.h"
#include <string>
#include <sstream>

#define LOGME(x,y)         do { if(_verbosity >= x) printf(y)        }while(0)
#define LOGME1(x,y,z)      do { if(_verbosity >= x) printf(y ,z);    }while(0)
#define LOGME2(x,y,z,w)    do { if(_verbosity >= x) printf(y, z, w); }while(0)

Deformation::Deformation(Mesh *mesh_in):
  _mesh(mesh_in),
  _rot(_mesh->numberOfVertices()),
  _xyz0(_mesh->numberOfVertices()),
  _xyz(_mesh->numberOfVertices()),
  _rhs(_mesh->numberOfVertices()),
  _A(NULL), _ATA(NULL), _ksp(NULL),
  _should_refactor(true),
  _beta(10),
  _tol(5e-3),
  _should_use_uniform_wgt(false),
  _n_max_iter(50),
  _verbosity(2),
  _n_total_deformations(0)
{

  assert(_mesh->numberOfEdges() > 0);
  const int n_vert = _mesh->numberOfVertices();

  /*
   * Read the options from command line
   */
  PetscBool flg;
  ierr=PetscOptionsGetInt(NULL, NULL, "-arap_max_iter", &_n_max_iter, NULL);CHKERRV(ierr);
  ierr=PetscOptionsGetReal(NULL, NULL, "-arap_beta", &_beta, NULL);CHKERRV(ierr);
  ierr=PetscOptionsGetReal(NULL, NULL, "-arap_tol", &_tol, NULL);CHKERRV(ierr);
  flg = (_should_use_uniform_wgt ? PETSC_TRUE : PETSC_FALSE);CHKERRV(ierr);
  ierr=PetscOptionsGetBool(NULL, NULL, "-arap_uni_wgt", &flg, NULL);CHKERRV(ierr);
  _should_use_uniform_wgt = (bool)flg;
  ierr=PetscOptionsGetInt(NULL, NULL, "-arap_v", &_verbosity, NULL);CHKERRV(ierr);


  /*
   *   Fill xyz, field and vertIndex
   */
  _vertIndex.reserve(n_vert);
  for (auto iv = _mesh->getVertices()->begin() ; iv != _mesh->getVertices()->end() ; ++iv)
    {
      Vertex *v = *iv;
      _vertIndex.push_back(v);
      _rot(v->name).eye();
      _xyz0(v->name)(0) = v->x();
      _xyz0(v->name)(1) = v->y();
      _xyz0(v->name)(2) = v->z();
    }

  /*
   * Find the edge weights
   */
  _edgew.reserve(_mesh->numberOfEdges());
  if(_should_use_uniform_wgt)
    {
      _edgew.resize(_mesh->numberOfEdges(), 1.);
    }
  else
    {
      for(auto ie = _mesh->getEdges()->begin(); ie != _mesh->getEdges()->end() ; ++ie)
	{
	  Edge *e = *ie;
	  _edgew.push_back(0);

	  for(int it=0 ; it < 2 ; it++)
	    {
	      Triangle *t = e->triangles[it];
	      if(t == Triangle::getSentinel()) continue;

	      Vertex *v0 = e->vertices[0];
	      Vertex *v1 = t->nextVertex(v0);
	      Vertex *v2 = t->nextVertex(v1);

	      if(v1 == e->vertices[1])
		{
		  swap(v2, v0);
		  swap(v0, v1);
		}
	      else assert (v2 == e->vertices[1]);

	      vec3d e12 = _xyz0[v2->name] - _xyz0[v1->name];
	      vec3d e10 = _xyz0[v0->name] - _xyz0[v1->name];
	      const double cos_ang = arma::norm_dot(e12, e10);
	      double ang;
	      const double tol = 1e-6;

	      if      ( (cos_ang > 1)  && (cos_ang < 1+tol ) )  ang = 0;
	      else if ( (cos_ang < -1) && (cos_ang > -1-tol) ) ang = M_PI;
	      else ang = acos(cos_ang);
	      double weight = 0.5 * 1/tan(ang);
	      assert( finite(weight) );

	      if(weight < 0.1)
		{
		  LOGME2(1, "Negative weight v: %d, %d \n", v0->name, v2->name);
		  weight = 0.1;
		}
	      _edgew.back() += weight;
	    }
	  //printf("edge_back: %d %d %lf \n", e->vertices[0]->name, e->vertices[1]->name, _edgew.back());
	} // End of for over edges
    } // End of else(use_uniform_weigt)

  /*
   * Create x and b.
   */
  ierr=VecCreateSeq(PETSC_COMM_SELF, n_vert, &_X.vec);CHKERRV(ierr);
  ierr=VecDuplicate(_X.vec, &_ATB.vec);CHKERRV(ierr);
  ierr=VecGetArray(_X.vec, &_X.v);CHKERRV(ierr);
  ierr=VecGetArray(_ATB.vec, &_ATB.v);CHKERRV(ierr);

  /*
   * Create the KSP
   */
  ierr=KSPCreate(PETSC_COMM_SELF, &_ksp);CHKERRV(ierr);


}

Deformation::~Deformation(){}

PetscErrorCode
Deformation::destroy()
{
  if(_A) {ierr=MatDestroy(&_A); CHKERRQ(ierr);}
  if(_ATA) {ierr=MatDestroy(&_ATA); CHKERRQ(ierr);}
  if(_ksp) {ierr=KSPDestroy(&_ksp);CHKERRQ(ierr);}
  if(_X.vec)
    {
      ierr=VecRestoreArray(_X.vec, &_X.v);CHKERRQ(ierr);
      ierr=VecDestroy(&_X.vec);CHKERRQ(ierr);
    }
  if(_B.vec)
    {
      ierr=VecRestoreArray(_B.vec, &_B.v);CHKERRQ(ierr);
      ierr=VecDestroy(&_B.vec);CHKERRQ(ierr);
    }
  if(_ATB.vec)
    {
      ierr=VecRestoreArray(_ATB.vec, &_ATB.v);CHKERRQ(ierr);
      ierr=VecDestroy(&_ATB.vec);CHKERRQ(ierr);
    }
  return 0;
}

PetscErrorCode
Deformation::assembleMat()
{

  if(!_should_refactor) return 0;

  /*
   * Create A and B.
   */
  if(_A) {ierr = MatDestroy(&_A);CHKERRQ(ierr);}
  if(_ATA) {ierr = MatDestroy(&_ATA);CHKERRQ(ierr);}
  if(_B.vec)
    {
      ierr=VecRestoreArray(_B.vec, &_B.v);CHKERRQ(ierr);
      ierr=VecDestroy(&_B.vec);CHKERRQ(ierr);
    }

  vector<int> n_nzero;

  const int n_col = _mesh->numberOfVertices();
  const int n_row = _mesh->numberOfVertices()+_anchors.size();
  n_nzero.reserve(n_row);
  for (int iv=0 ; iv < n_col; iv++)
    {
      n_nzero.push_back( _vertIndex[iv]->getEdges()->size()+1 );
    }
  for (int ia=0 ; ia < (int)_anchors.size() ; ia++)
    {
      n_nzero.push_back(1);
    }
  ierr=MatCreateSeqAIJ(PETSC_COMM_SELF, n_row, n_col, -1, n_nzero.data(), &_A); CHKERRQ(ierr);
  ierr=VecCreateSeq(PETSC_COMM_SELF, n_row, &_B.vec);CHKERRQ(ierr);
  ierr=VecGetArray(_B.vec, &_B.v);CHKERRQ(ierr);


  /*
   * Assemble the rigidity
   */
  int idx[2];
  double lmat[4];

  //ierr=MatZeroEntries(_A);CHKERRQ(ierr);

  // Assemble over edges
  for(auto ie = _mesh->getEdges()->begin(); ie != _mesh->getEdges()->end() ; ++ie)
    {
      Edge *e = *ie;
      idx[0] = e->vertices[0]->name; idx[1] = e->vertices[1]->name;
      lmat[0] =  _edgew[e->name];  lmat[1] = -_edgew[e->name];
      lmat[2] = -_edgew[e->name];  lmat[3] =  _edgew[e->name];
      ierr= MatSetValues(_A, 2, idx, 2, idx, lmat, ADD_VALUES);CHKERRQ(ierr);
    }

  /*
   * Assemble boundary condition constraints.
   */
  assert(_anchors.size());

  int irow = n_col;
  for (auto ian = _anchors.begin() ; ian != _anchors.end() ; ian++)
    {
      const int icol = ian->first;
      ierr= MatSetValue(_A, irow, icol, 1.*_beta, ADD_VALUES);CHKERRQ(ierr);
      irow++;
    }

  ierr=MatAssemblyBegin(_A, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr=MatAssemblyEnd(_A, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  //ierr=MatView(_A, PETSC_VIEWER_DRAW_SELF);CHKERRV(ierr);

  /*
   * Create the ATA
   */
  ierr=MatTransposeMatMult(_A,_A,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&_ATA);CHKERRQ(ierr);

  //ierr=MatView(_ATA, PETSC_VIEWER_DRAW_SELF);CHKERRV(ierr);

  /*
   * Set KSP to solve A
   */
  PC pc;
  ierr=KSPSetOperators(_ksp, _ATA, _ATA);CHKERRQ(ierr);

  // default options
  ierr=KSPSetType(_ksp, "preonly");CHKERRQ(ierr);
  ierr=KSPGetPC(_ksp, &pc);CHKERRQ(ierr);
  ierr=PCSetType(pc, "cholesky");CHKERRQ(ierr);

  // override default options from command line
  ierr=KSPSetFromOptions(_ksp);CHKERRQ(ierr);

  /*
   * Set the flag.
   */
  _should_refactor = false;

  return 0;
}

PetscErrorCode
Deformation::assembleRHS()
{
  const int n_col = _mesh->numberOfVertices();
  const int n_row = _mesh->numberOfVertices()+_anchors.size();

  _rhs.set_size(n_row);
  for (int iv=0 ; iv < _mesh->numberOfVertices(); iv++)
    {
      _rhs(iv).zeros();
    }

  for(auto ie = _mesh->getEdges()->begin(); ie != _mesh->getEdges()->end() ; ++ie)
    {
      Edge *e = *ie;
      const int iv0 = e->vertices[0]->name;
      const int iv1 = e->vertices[1]->name;
      vec3d val = 0.5 * _edgew[e->name] * (_rot[iv0] + _rot[iv1]) * (_xyz0[iv1] - _xyz0[iv0]);
      _rhs[iv1] += val;
      _rhs[iv0] -= val;
    }

  // Set the boundary conditions
  assert(_anchors.size());

  int irow = n_col;
  for (auto ian = _anchors.begin() ; ian != _anchors.end() ; ian++)
    {
      //const int icol = ian->first;
      vec3d &val = ian->second;
      _rhs(irow) = val * _beta;
      irow++;
    }

  return 0;
}

PetscErrorCode
Deformation::updateRots()
{
  // Init the S matrix as zeros
  arma::field<matrix3d> smat(_mesh->numberOfVertices());
  matrix3d value,U,V;
  vec3d S;

  for (auto iv = _mesh->getVertices()->begin() ; iv != _mesh->getVertices()->end() ; ++iv)
    {
      smat((*iv)->name).zeros();
    }

  // Assemble them over each edge
  for(auto ie = _mesh->getEdges()->begin(); ie != _mesh->getEdges()->end() ; ++ie)
    {
      Edge *e = *ie;
      Vertex *vi = e->vertices[0];
      Vertex *vj = e->vertices[1];
      vec3d eij = _xyz0[vi->name] - _xyz0[vj->name];
      vec3d eijpr = _xyz[vi->name] - _xyz[vj->name];

      value = _edgew[e->name] * eij * eijpr.t();
      smat(vi->name) += value;
      smat(vj->name) += value;
    }

  for (auto iv = _mesh->getVertices()->begin() ; iv != _mesh->getVertices()->end() ; ++iv)
    {
      Vertex *v = *iv;
      ierr=!arma::svd(U, S, V, smat(v->name), "std");
      if(ierr)
	{
	  SETERRQ1(PETSC_COMM_SELF, 1, "SVD failed on vertex %d", v->name);
	}
      else if (arma::det(U) * arma::det(V) < 0 )
	{
	  U.col(2) *= -1;
	}
      _rot(v->name) = V * U.t();
      //		_rot(v->name).print();
    }

  return 0;
}


PetscErrorCode
Deformation::deform(vector<int> anchorverts, vector<float> newpos)
{

  // Return if nothing is selected
  if(anchorverts.size() == 0) return 0;

  /*
   * Preprocess anchors.
   */

  // See if anchors are changed.
  bool changed = false;
  if (anchorverts.size() != _anchors.size()) changed = true;
  else
    {
      for (auto it = anchorverts.begin() ; it!= anchorverts.end() ; ++it)
	{
	  if (_anchors.find(*it) == _anchors.end())
	    {
	      changed= true;
	      break;
	    }
	}
    }

  // Update the anchor locations.
  {
    _anchors.clear();
    int j;
    vector<int>::iterator it;
    for (it = anchorverts.begin(), j=0 ; it!= anchorverts.end() ; ++it, ++j)
      {
	vec3d pos;
	pos(0)=newpos[3*j+0];
	pos(1)=newpos[3*j+1];
	pos(2)=newpos[3*j+2];
	_anchors.insert( make_pair(*it, pos ));
      }
  }

  _should_refactor = _should_refactor || changed ;

  /*
   * Get the initial guess from the mesh.
   */
  for (auto iv = _mesh->getVertices()->begin() ; iv != _mesh->getVertices()->end() ; ++iv)
    {
      Vertex *v = *iv;
      _xyz(v->name)(0) = v->x();
      _xyz(v->name)(1) = v->y();
      _xyz(v->name)(2) = v->z();
    }


  /*
   * Extrace the RHS and solve.
   */
  int iiter;
  double res;
  vec3d delta;
  stringstream ss;
  for (iiter = 0 ; iiter < _n_max_iter ; iiter++)
    {
      // Reset residual
      res = -1e6;

      // Write a file
      if(_verbosity == 2)
	{
	  ss.str("");
	  ss << "svd_" << _n_total_deformations << "_" << iiter << ".vtk";
	  MeshIO(_mesh).write_vtk(ss.str());
	}

      // Assemble the matrix and RHS
      ierr=this->assembleMat();CHKERRQ(ierr);
      ierr=this->assembleRHS();CHKERRQ(ierr);

      // Solve the linear system of equations
      for (int dim=0 ; dim < 3 ; dim++)
	{
	  extractDim(dim,_rhs, _B);
	  ierr=VecAssemblyBegin(_B.vec);CHKERRQ(ierr);
	  ierr=VecAssemblyEnd(_B.vec);CHKERRQ(ierr);
	  ierr=MatMultTranspose(_A,_B.vec, _ATB.vec);CHKERRQ(ierr);
	  ierr=KSPSolve(_ksp, _ATB.vec, _X.vec);CHKERRQ(ierr);
	  // Find the change in the positions and return if we have converged
	  for (auto iv = _mesh->getVertices()->begin() ; iv != _mesh->getVertices()->end() ; ++iv)
	    {
	      Vertex *v = *iv;
	      res = std::max(res, std::abs( _X.v[v->name] -_xyz(v->name)(dim) ) );
	    }
	  mergeDim(dim, _X, _xyz);
	}
      LOGME1(2, "Iteration solution change : %lf \n", res);

      // Find the new rotation matrices.
      ierr=this->updateRots();CHKERRQ(ierr);


      // Propagate the solution to the mesh.
      // We do it in every step so we can write the mesh.
      if(_verbosity == 2)
	{
	  for (auto iv = _mesh->getVertices()->begin() ; iv != _mesh->getVertices()->end() ; ++iv)
	    {
	      Vertex *v = *iv;
	      v->x()  = _xyz(v->name)(0) ;
	      v->y()  = _xyz(v->name)(1) ;
	      v->z()  = _xyz(v->name)(2) ;
	    }
	}

      if(res < _tol)	break;
    }

  // Log convergance
  if (iiter < _n_max_iter) LOGME2(1, "ARAP converged, n_step= %d, final res = %lf \n", iiter+1, res);
  else LOGME2(1, "ARAP did not converge, n_step= %d, final res = %lf \n", iiter, res);

  // Write a last file
  if(_verbosity == 2)
    {
      ss.str("");
      ss << "svd_" << _n_total_deformations << "_" << iiter << ".vtk";
      MeshIO(_mesh).write_vtk(ss.str());
    }
  else if(_verbosity == 1)
    {
      ss.str("");
      ss << "svd_" << _n_total_deformations  << ".vtk";
      MeshIO(_mesh).write_vtk(ss.str());
    }

  /*
   * Propagate the solution to the mesh
   */
  // update positions
  for (auto iv = _mesh->getVertices()->begin() ; iv != _mesh->getVertices()->end() ; ++iv)
    {
      Vertex *v = *iv;
      v->x()  = _xyz(v->name)(0) ;
      v->y()  = _xyz(v->name)(1) ;
      v->z()  = _xyz(v->name)(2) ;
      v->resetNormal();
    }
  // update normals.
  for (auto it = _mesh->getTriangles()->begin() ; it != _mesh->getTriangles()->end() ; ++it)
    {
      Triangle *t = *it;
      t->calcProperties();
      for (int tv=0; tv < 3; tv++) t->vertices[tv]->addNormal(t->mathNormal());
    }


  _n_total_deformations++;
  return 0;
}

PetscErrorCode
Deformation::resetMeshPositions()
{
  for (auto iv = _mesh->getVertices()->begin() ; iv != _mesh->getVertices()->end() ; ++iv)
    {
      Vertex *v = *iv;
      v->x()  = _xyz0(v->name)(0) ;
      v->y()  = _xyz0(v->name)(1) ;
      v->z()  = _xyz0(v->name)(2) ;
      v->resetNormal();
    }
  for (auto it = _mesh->getTriangles()->begin() ; it != _mesh->getTriangles()->end() ; ++it)
    {
      Triangle *t = *it;
      t->calcProperties();
      for (int tv=0; tv < 3; tv++) t->vertices[tv]->addNormal(t->mathNormal());
    }

  return 0;
}
