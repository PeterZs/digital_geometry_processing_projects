#include <stack>
#include <sstream>
#include <cstdlib>
#include <iostream>
#include <chrono>
#include <random>
#include <algorithm>

#include "adaptive_remesher.hxx"
#include "mesh_io.hxx"
#include "predicates.hxx"

static double acos_tol(const double val, const double tol = 1e-6)
{
	if( (val < -1-tol) || (val > 1+tol) )
	{
		printf("[Error] Bad value for acos %lf \n", val);
		assert(0); throw;
	}

	if(val < -1) return M_PI;
	else if(val > 1) return 0;
	else return acos(val);
}



void AdaptiveRemesher::init()
{
	_proj->init();
}

AdaptiveRemesher::AdaptiveRemesher(Mesh& mesh_in, GetPot& getpot, AreaRemeshMode mode_in):
_mesh(mesh_in),
_arearemesh_mode(mode_in),
_target_verts(mesh_in.verts.n_mems()),
_ctheta_tri( -1 ),
_ctheta_vert( -1 ),
_proj(Projector::build(getpot, mesh_in))
{
	const double ttri = getpot("ttri", 20.);
	const double tvert = getpot("tvert", 20.);

	_ctheta_tri = cos(M_PI / 180. * ttri);
	_ctheta_vert = cos(M_PI / 180. * tvert);
	printf("[log] adaptive remesher created, ttri = %.2lf, tvert=%lf \n" , ttri, tvert);
}

void AdaptiveRemesher::remesh_flipping()
{
	int n_flip = 0;
	//std::stringstream ss;
	//ss << "flip" << n_flip << ".vtk";
	//MeshIO(*mesh()).write_vtk(ss.str());

	std::stack<Edge*> cands;

	for(int iedge=0 ; iedge < mesh()->edges.n_mems() ; iedge++)
	{
		Edge *edge = mesh()->edges.get_member_ptr(iedge);
		edge->is_marked = this->should_flip(edge);
		if(edge->is_marked) cands.push(edge);
	}

	while(!cands.empty())
	{
		Edge *edge = cands.top();
		cands.pop();
		edge->is_marked = false;

		// printf("Inspecting %d %d \n", edge->he->vert->idx(), edge->he->twin->vert->idx());
		if( this->should_flip(edge) )
		{
			n_flip++;
			//ss.str("");
			// printf("Flipping %d %d \n", edge->he->vert->idx(), edge->he->twin->vert->idx());
			mesh()->flip_edge(edge);
			//ss << "flip" << n_flip << ".vtk";
			//MeshIO(*mesh()).write_vtk(ss.str());

			Edge *edge2;
			edge2 = edge->he->next->edge;               if(!edge2->is_marked) cands.push(edge2);
			edge2 = edge->he->next->next->edge;         if(!edge2->is_marked) cands.push(edge2);
			edge2 = edge->he->twin->next->edge;         if(!edge2->is_marked) cands.push(edge2);
			edge2 = edge->he->twin->next->next->edge;   if(!edge2->is_marked) cands.push(edge2);

		} // End of if(should_flip)
	} // End of while(stack is not empty)
	std::cout << "[Log] Total " << n_flip << " delaunay flips. " << std::endl;
}

void AdaptiveRemesher::random_flipping()
{

	//int n_flip = 0;
	//std::stringstream ss;
	//ss << "rflip" << n_flip << ".vtk";
	//MeshIO(*mesh()).write_vtk(ss.str());

	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::minstd_rand0 generator (seed);

	const int n_edge = mesh()->edges.n_mems();

	for (int i = 0 ; i < n_edge ; i++)
	{
		// get a random edge
		const int e = generator() % n_edge;
		Edge *edge = mesh()->edges.get_member_ptr(e);

		if(this->can_flip(edge))
		{
			// n_flip ++;
			mesh()->flip_edge(edge);
			//ss.str("");
			//ss << "rflip" << n_flip << ".vtk";
			//MeshIO(*mesh()).write_vtk(ss.str());
		}
	}

}

bool AdaptiveRemesher::remesh_areabased(Vertex* v0, Mesh& m2d)
{
	/*
	 * Cannot do anything about the boundary
	 */
	Mesh::ring_iterator it;
	it.init(v0);
	if(it.reset_boundary()) return false;

	/*
	 * Map the neighbourhood to 2d
	 */
	this->geodesic_vert(v0, m2d);

	/*
	 * Find the new location of the vertex in 2d
	 */
	std::vector<double> mu, a, b, c, A;
	double mutot(0), Atot(0), mufirst;
	int n_neigh(0);

	// Find all the weights
	Vertex *v02d = m2d.verts.get_member_ptr(0);
	it.init(v02d);
	do
	{
		Vertex *vi = it.half_edge()->vert;
		Vertex *vii = it.half_edge()->next->next->vert;

		n_neigh++;
		a.push_back(vi->xyz(1)-vii->xyz(1));
		b.push_back(vii->xyz(0)-vi->xyz(0));
		c.push_back(vi->xyz(0)*vii->xyz(1) - vi->xyz(1)*vii->xyz(0));
		mu.push_back( 1. /*this->vertex_areabased_mu(vi) */);
		A.push_back( a.back() * v02d->xyz(0) + b.back() * v02d->xyz(1) + c.back() );

		mutot += mu.back();
		Atot += A.back();

	}while(it.advance());

	// Find the correct mu
	mufirst = mu[0];
	for (int i = 0 ; i < n_neigh-1 ; i++)
	{
		mu[i] = 0.5*(mu[i]+mu[i+1]) / mutot;
		//printf("%lf \n", mu[i]);
	}
	//printf("\n");
	mu.back() = 0.5 * ( mu.back() + mufirst ) / mutot;

	// Form the system to solve
	arma::mat LHS(n_neigh, 2), Q, R;
	arma::vec RHS(n_neigh), QtRHS, sol;

	for (int i = 0 ; i < n_neigh ; i++)
	{
		LHS(i,0) = a[i]; LHS(i,1) = b[i];
		RHS(i) = -c[i] + mu[i] * Atot;
	}

	// Solve using QR
	if(!arma::qr(Q, R, LHS)) return false;
	QtRHS = Q.t() * RHS;
	arma::solve(sol, R.submat(0,0,1,1), QtRHS.subvec(0,1) );

	/*
	 * Map back to 3d
 	 */
	bool ierr;
	Vertex vertnew;
	Vector3 abc, sol3;
	Face *f2d = m2d.faces.get_member_ptr(0), *f3d;

	sol3(0)=sol(0); sol3(1)=sol(1); sol3(2)=0;
	ierr = m2d.locate_point(sol3, f2d, abc);
	// debug
	if(!ierr)
	{
		return false;
		printf("ivert: %d \n", v0->idx());
		MeshIO(*this->mesh()).write_vtk("geodesic_master.vtk");
		MeshIO(m2d).write_vtk("geodesic0.vtk");
		m2d.verts.get_member_ptr(0)->xyz = sol3;
		MeshIO(m2d).write_vtk("geodesic1.vtk");
		assert(ierr);
	}

	f3d = static_cast<Face*>(f2d->misc[0]);
	proj()->find_point(f3d,abc,&vertnew);
	//	MeshIO(m2d).write_vtk("areabased1.vtk");
	//	v02d->xyz = sol3;
	//	MeshIO(m2d).write_vtk("areabased2.vtk");


	/*
	 * Check that the new location adheres to the error metrics
	 */
	Vertex *cands[3];
	double Q1, Q2;
	it.init(v0);
	do
	{
		cands[0]=it.half_edge()->vert;
		cands[1]=&vertnew;
		cands[2]=it.half_edge()->next->next->vert;
		if(!this->tri_fidelity(cands, Q1, Q2)) return false;
	}while(it.advance());
	proj()->copy_vert_data(&vertnew, v0);

	return true;
}

void AdaptiveRemesher::remesh_areabased()
{
	Mesh m2d;
	for(int i = 0 ; i < mesh()->verts.n_mems() ; i++)
	{
		Vertex *v = mesh()->verts.get_member_ptr(i);
		this->remesh_areabased(v , m2d);
	}
}

void AdaptiveRemesher::remesh_areabased(const int n_step, const int n_area)
{
	for (int step = 0 ; step < n_step ; step++)
	{
		for (int area = 0 ; area < n_area ; area++)
		{
			this->remesh_areabased();
		}
		this->remesh_flipping();
	}
}

bool AdaptiveRemesher::remesh_smooth_laplacian(Vertex* v0, Mesh &m2d)
{
	/*
	 * Cannot do anything about the boundary
	 */
	Mesh::ring_iterator it;
	it.init(v0);
	if(it.reset_boundary()) return false;

	/*
	 * Map the neighbourhood to 2d
	 */
	this->geodesic_vert(v0, m2d);
	// MeshIO(m2d).write_vtk("geodesic.vtk");

	/*
	 * Find the new location in 3d
	 */
	it.init(m2d.verts.get_member_ptr(0));
	Vector3 xyz2dnew = arma::fill::zeros;
	int n_neigh = 0;
	do
	{
		xyz2dnew += it.half_edge()->vert->xyz;
		n_neigh++;
	}while(it.advance());
	xyz2dnew /= n_neigh;

	/*
	 * Map back to 3d
	 */
	bool ierr;
	Vertex vertnew;
	Vector3 abc;
	Face *f2d = m2d.faces.get_member_ptr(0), *f3d;

	ierr = m2d.locate_point(xyz2dnew, f2d, abc);
	if(!ierr)
	{
		return false;
		printf("ivert: %d \n", v0->idx());
		MeshIO(*this->mesh()).write_vtk("geodesic_master.vtk");
		MeshIO(m2d).write_vtk("geodesic0.vtk");
		m2d.verts.get_member_ptr(0)->xyz = xyz2dnew;
		MeshIO(m2d).write_vtk("geodesic1.vtk");
		assert(ierr);
	}

	f3d = static_cast<Face*>(f2d->misc[0]);
	proj()->find_point(f3d,abc,&vertnew);

	/*
	 * Check that the new location adheres to the error metrics
	 */
	Vertex *cands[3];
	double Q1, Q2;
	it.init(v0);
	do
	{
		cands[0]=it.half_edge()->vert;
		cands[1]=&vertnew;
		cands[2]=it.half_edge()->next->next->vert;
		if(!this->tri_fidelity(cands, Q1, Q2)) return false;
	}while(it.advance());
	proj()->copy_vert_data(&vertnew, v0);

	return true;
}

void AdaptiveRemesher::remesh_smooth_laplacian(const int n_step)
{
	Mesh m2d;
	for (int  i = 0 ; i < n_step ; i++)
	{
		for(int i = 0 ; i < mesh()->verts.n_mems() ; i++)
		{
			Vertex *v = mesh()->verts.get_member_ptr(i);
			this->remesh_smooth_laplacian(v , m2d);
		}
	}
}


bool AdaptiveRemesher::is_obtuse(HalfEdge *he) const
{
	//double Q1, Q2;
	Vertex *v[] = {he->next->next->vert, he->vert, he->next->vert};

	if(he->face == Mesh::hole()) return false;
	// if( !this->tri_fidelity(v, Q1, Q2) ) return false;

	const Vector3 e01 = arma::normalise( v[1]->xyz - v[0]->xyz );
	const Vector3 e02 = arma::normalise( v[2]->xyz - v[0]->xyz );
	const double angle = acos_tol( arma::dot(e01, e02) );
	return angle > M_PI/2.;
}

bool AdaptiveRemesher::is_obtuse(Edge *edge) const
{
	return is_obtuse(edge->he) || is_obtuse(edge->he->twin);
}

void AdaptiveRemesher::remesh_obtuse_angles()
{
	int n_insert = 0;
	std::stack<Edge*> cands;

	for(int iedge=0 ; iedge < mesh()->edges.n_mems() ; iedge++)
	{
		Edge *edge = mesh()->edges.get_member_ptr(iedge);
		edge->is_marked = is_obtuse(edge);
		if(edge->is_marked) cands.push(edge);
	}

	Vertex *v;
	Edge  *edge;
	//Mesh::ring_iterator it; Edge *edge2;

	while(!cands.empty())
	{
		edge = cands.top();
		cands.pop();
		edge->is_marked = false;

		if( is_obtuse(edge) )
		{
			v = this->split_edge(edge);
			if(!v) continue;
			n_insert++;
		} // End of if(should_flip)
	} // End of while(stack is not empty)
	std::cout << "[Log] Inserted total " << n_insert << " verts to remove obtuse angles. " << std::endl;
}

bool AdaptiveRemesher::remesh_splitedges()
{
	int n_inserted = 0;
	std::vector<double> pool_error( mesh()->edges.n_mems(), 0 );
	std::vector<Edge *> stack(mesh()->edges.n_mems(), NULL);
	Mesh::ring_iterator itring;

	/*
	 * Find the error for each edge
	 * Split the triangles that have the least quality.
	 */
	for (int it = 0 ; it < mesh()->faces.n_mems() ; it++)
	{
		Face *face = mesh()->faces.get_member_ptr(it);
		HalfEdge *he;

		const double error = 60. - this->tri_quality(face);
		he = face->he; pool_error[ he->edge->idx() ] += error*error;
		he = he->next; pool_error[ he->edge->idx() ] += error*error;
		he = he->next; pool_error[ he->edge->idx() ] += error*error;

		face->is_marked = false;
	}
	Mesh::hole()->is_marked = false;

	/*
	 * set the reference right for each edge
	 */
	for (int ie = 0 ; ie < mesh()->edges.n_mems() ; ie++)
	{
		Edge *edge = mesh()->edges.get_member_ptr(ie);
		edge->misc.resize(1);
		edge->misc[0] = &pool_error[ie];
		stack[ie] = edge;
		if(edge->he->face == Mesh::hole()) pool_error[ie]+=1000;
	}

	/*
	 * Sort the edges
	 */
	std::sort(stack.begin(), stack.end(), Edge::cmp_functor<double, 0>());

	/*
	 * Start splitting edges
	 */
	for (int i = stack.size()-1 ; i >= 0 ; i--)
	{
		Edge *edge = stack[i];

		// Split the edge if neighbours are not marked
		if ( edge->he->face->is_marked || edge->he->twin->face->is_marked ) continue;

		// Check the insertion to be good.
		const double rat = 1;
		HalfEdge *he= edge->he;
		if(he->face != Mesh::hole())
		{
			const double l1 = arma::norm(he->vert->xyz - he->next->vert->xyz); he=he->next;
			const double l2 = arma::norm(he->vert->xyz - he->next->vert->xyz); he=he->next;
			const double l3 = arma::norm(he->vert->xyz - he->next->vert->xyz); he=he->next;
			if( (l1 < l2/rat) || (l1 < l3/rat) ) continue;
		}
		he = he->twin;
		{
			const double l1 = arma::norm(he->vert->xyz - he->next->vert->xyz); he=he->next;
			const double l2 = arma::norm(he->vert->xyz - he->next->vert->xyz); he=he->next;
			const double l3 = arma::norm(he->vert->xyz - he->next->vert->xyz); he=he->next;
			if( (l1 < l2/rat) || (l1 < l3/rat) ) continue;
		}

		// Find the location of the new vertex, midpoint for now!!
		Vertex *vertnew = this->split_edge(edge);
		if(!vertnew) continue;
		itring.init(vertnew);
		do
		{
			itring.half_edge()->face->is_marked = true;
		}while(itring.advance());
		Mesh::hole()->is_marked = false;
		n_inserted++;

		// Check if we have inserted enough
		if( mesh()->verts.n_mems() > n_tverts() )
		{
			std::cout << "[log] inserted " << n_inserted
					<< " points, finished with inserting." << std::endl;
			return false;
		}

	} // End of all candidates

	std::cout << "[log] inserted " << n_inserted
			<< " points, will insert another round." << std::endl;
	return true;
} // End of split edges

bool AdaptiveRemesher::remesh_collapseedges()
{
	int n_collapsed = 0;
	std::vector<double> pool_error( mesh()->edges.n_mems(), 0 );
	std::vector<Edge *> stack(mesh()->edges.n_mems(), NULL);
	Vertex *vert;
	Mesh::ring_iterator itring;

	/*
	 * Find the error for each edge
	 * Unmark all triangles
	 */
	for (int it = 0 ; it < mesh()->faces.n_mems() ; it++)
	{
		Face *face = mesh()->faces.get_member_ptr(it);
		HalfEdge *he;

		const double error = this->tri_quality(face);
		he = face->he; pool_error[ he->edge->idx() ] = std::max(error, pool_error[ he->edge->idx() ] );
		he = he->next; pool_error[ he->edge->idx() ] = std::max(error, pool_error[ he->edge->idx() ] );
		he = he->next; pool_error[ he->edge->idx() ] = std::max(error, pool_error[ he->edge->idx() ] );

		face->is_marked = false;
	}
	Mesh::hole()->is_marked = false;

	/*
	 * set the reference right for each edge
	 */
	for (int ie = 0 ; ie < mesh()->edges.n_mems() ; ie++)
	{
		Edge *edge = mesh()->edges.get_member_ptr(ie);
		stack[ie] = edge;
		edge->misc.resize(2);
		edge->misc[0] = &pool_error[ie];
	}

	/*
	 * Sort the edges
	 */
	std::sort(stack.begin(), stack.end(), Edge::cmp_functor<double, 0>());
	for (int ie = 0 ; ie < mesh()->edges.n_mems() ; ie++)
	{
		stack[ie]->misc[1] = &stack[ie];
	}

	/*
	 * Start collapsing edges
	 */
	for (int i = stack.size()-1 ; i >= 0 ; i--)
	{
		Edge *edge = stack[i];

		// If the edge still exists and is not on the boundary
		// if(!edge) printf("prev deleted \n");
		// else printf("%d \n", edge->idx());
		if ( (!edge) || (edge->he->face == Mesh::hole())
				|| (edge->he->face->is_marked) || (edge->he->twin->face->is_marked) ) continue;

		// Get the location of deleted faces in the stack to
		// set them to NULL after deleting
		Edge **stack_e[5];
		stack_e[0] = static_cast<Edge**>(edge->he->edge->misc[1]);
		stack_e[1] = static_cast<Edge**>(edge->he->next->edge->misc[1]);
		stack_e[2] = static_cast<Edge**>(edge->he->next->next->edge->misc[1]);
		stack_e[3] = static_cast<Edge**>(edge->he->twin->next->edge->misc[1]);
		stack_e[4] = static_cast<Edge**>(edge->he->twin->next->next->edge->misc[1]);

		assert( stack_e[0][0] == edge);
		assert( stack_e[1][0] == edge->he->next->edge);
		assert( stack_e[2][0] == edge->he->next->next->edge);
		assert( stack_e[3][0] == edge->he->twin->next->edge);
		assert( stack_e[4][0] == edge->he->twin->next->next->edge);


		// Find the location of the new vertex
		// midpoint for now!!
		vert = this->collapse_edge(edge);
		if( !vert ) continue;
		itring.init(vert);
		do
		{
			itring.half_edge()->face->is_marked = true;
		}while(itring.advance());
		Mesh::hole()->is_marked = false;
		stack_e[0][0] = stack_e[1][0] = stack_e[2][0] =
				stack_e[3][0] = stack_e[4][0] = static_cast<Edge*>(NULL);
		n_collapsed++;

		// Check if we have collapsed enough
		if( mesh()->verts.n_mems() < n_tverts() )
		{
			std::cout << "[log] Collapsed " << n_collapsed
					<< " edges, finished with collapsing." << std::endl;
			return false;
		}

	} // End of all candidates

	if(n_collapsed > 0){
		std::cout << "[log] Collapsed " << n_collapsed
				<< " edges, will collapse another round." << std::endl;
		return true;
	}else{
		std::cout << "[log] Collapsed " << n_collapsed
				<< " edges, finished with collapsing." << std::endl;
		return false;
	}
}

bool AdaptiveRemesher::should_flip(Edge * edge) const
{
	if(!this->can_flip(edge)) return false;

	// Map the edge to two dimesions
	Vector3 xyz2d[4];
	this->geodesic_edge(edge, xyz2d);
	const double pred = incircle(xyz2d[0].memptr(), xyz2d[1].memptr(), xyz2d[2].memptr(), xyz2d[3].memptr());
	if (pred > 0) return true;
	else return false;

	/* Alternatively use length
	Vertex *v[4];
	v[0] = edge->he->vert;
	v[1] = edge->he->twin->vert;
	v[2] = edge->he->next->next->vert;
	v[3] = edge->he->twin->next->next->vert;
	const Vector3 enow = v[0]->xyz - v[1]->xyz;
	const Vector3 enew = v[2]->xyz - v[3]->xyz;
	const double lnow = arma::norm(enow);
	const double lnew = arma::norm(enew);

	if( this->can_flip(edge) &&  (lnew < lnow) ) return true;
	else return false;
	*/
}

bool AdaptiveRemesher::can_flip(Edge * edge) const
{
	if(edge->he->face == Mesh::hole()) return false;

	Vertex *v[4];
	v[0] = edge->he->vert;
	v[1] = edge->he->twin->vert;
	v[2] = edge->he->next->next->vert;
	v[3] = edge->he->twin->next->next->vert;

	// See if the fidelity is preseverted
	Vertex *tris_new[3];
	double Qtri, Qvert;
	tris_new[0]=v[3]; 	tris_new[1]=v[1]; 	tris_new[2]=v[2];
	if( !tri_fidelity(tris_new, Qtri, Qvert) ) return false;
	tris_new[0]=v[0]; 	tris_new[1]=v[3]; 	tris_new[2]=v[2];
	if( !tri_fidelity(tris_new, Qtri, Qvert) ) return false;

	return true;

	/* Compare normals, in three dimenstions */

	/* Orientation, only for two dimensions
	Vertex *v1 = edge->he->vert;
	Vertex *v2 = edge->he->twin->vert;
	Vertex *v3 = edge->he->next->next->vert;
	Vertex *v4 = edge->he->twin->next->next->vert;
	const double pred1 = orient2d(v3->xyz.memptr(), v4->xyz.memptr(), v1->xyz.memptr());
	const double pred2 = orient2d(v3->xyz.memptr(), v4->xyz.memptr(), v2->xyz.memptr());
	return ( (pred1 > 0) && (pred2 < 0)  )  || ( (pred2 > 0) && (pred1 < 0)  );
	*/
}

double AdaptiveRemesher::vertex_areabased_mu(Vertex* v) const
{
	switch(arearemesh_mode())
	{
	case CONSTANT:
	{
		return 1. ;
	}
	case CURVATURE:
	{
		assert(0);
		throw;
		return -1;
	}
	case REFINE00:
	{
		// std::cout <<  pow(arma::norm(v->xyz)+0.1 , 20) << std::endl;
		Vector3 cent ;
		cent(0)=cent(1)=1; cent(2)=0;
		double r = arma::norm(v->xyz-cent );
		if(r < 0.05 ) r = 0.05;
		//if(r<0.005) return 1;
		//else if(r<0.3) return r*r*r*r*r*r;
		//else return r*r*r;
		// return r*r*r*r*r*r + 0.01;
		return r ;

	}
	} // end switch

	// die a horrible death
	return 5000;
}


bool  AdaptiveRemesher::tri_fidelity(Vertex *v[], double &Qtri, double &Qvert) const
{
	// Get vertex normals
	Vector3 n_tri, e01, e02;

	e01 = v[1]->xyz - v[0]->xyz;
	e02 = v[2]->xyz - v[0]->xyz;
	n_tri = arma::normalise( arma::cross(e01, e02) );

	/*
	 * Triangle normal vs. vertex normals
	 */
	Qtri = -1.;
	Qtri = std::max( Qtri, double( arma::dot(n_tri, v[0]->normal) ) );
	Qtri = std::max( Qtri, double( arma::dot(n_tri, v[1]->normal) ) );
	Qtri = std::max( Qtri, double( arma::dot(n_tri, v[2]->normal) ) );

	/*
	 * Vertex normals compared to one another
	 */
	Qvert = -1.;
	Qvert = std::max( Qvert, double( arma::dot(v[0]->normal, v[1]->normal) ) );
	Qvert = std::max( Qvert, double( arma::dot(v[1]->normal, v[2]->normal) ) );
	Qvert = std::max( Qvert, double( arma::dot(v[0]->normal, v[2]->normal) ) );

	return (Qtri > _ctheta_tri) && (Qvert > _ctheta_vert) ;
}

double AdaptiveRemesher::tri_quality(Face *f) const
{
	Vertex *v1 = f->he->vert;
	Vertex *v2 = f->he->next->vert;
	Vertex *v3 = f->he->next->next->vert;

	/* For now define the error as 60 - min_angle */
	double min_angle = M_PI;
	Vector3 e12  = v2->xyz - v1->xyz; e12 /= arma::norm(e12);
	Vector3 e23  = v3->xyz - v2->xyz; e23 /= arma::norm(e23);
	Vector3 e31  = v1->xyz - v3->xyz; e31 /= arma::norm(e31);

	const double ang1 = acos_tol(arma::dot(-e12, e23));
	const double ang2 = acos_tol(arma::dot(-e23, e31));
	const double ang3 = acos_tol(arma::dot(-e31, e12));
	min_angle = std::min(min_angle, ang1);
	min_angle = std::min(min_angle, ang2);
	min_angle = std::min(min_angle, ang3);

	return 180./M_PI * min_angle;
}

Vertex* AdaptiveRemesher::split_edge(Edge * edge)
{

	/*
	 * Find the new location
	 */
	HalfEdge *he = edge->he->twin;
	Vertex *v = new Vertex;
	Vector3 abc;
	if(he->face->he == he)                  abc = "0.5 0.5 0";
	else if(he->face->he == he->next)       abc = "0.5 0 0.5";
	else if(he->face->he == he->next->next) abc = "0 0.5 0.5";
	else
	{
		printf("[Error]");
		throw;
	}
	proj()->find_point(he->face, abc, v);

	/*
	 * Call fidelity
	 * not needed for now
	 */
//	Vertex *vs[] = {he->vert, he->next->vert, he->next->next->vert, he->twin->next->next->vert};
//	Vertex *cands[]= {v, NULL, NULL};
//	double Q1, Q2;
//	cands[1]=vs[1]; cands[2]=vs[2]; if(!this->tri_fidelity(cands, Q1, Q2)){  delete v; return NULL; }
//	cands[1]=vs[2]; cands[2]=vs[0]; if(!this->tri_fidelity(cands, Q1, Q2)){  delete v; return NULL; }
//	cands[1]=vs[0]; cands[2]=vs[3]; if(!this->tri_fidelity(cands, Q1, Q2)){  delete v; return NULL; }
//	cands[1]=vs[3]; cands[2]=vs[1]; if(!this->tri_fidelity(cands, Q1, Q2)){  delete v; return NULL; }

	/*
	 * check if we can insert this
	 */
	if( mesh()->split_edge(edge, v) ) return v;
	else
	{
		delete v;
		return NULL;
	}
}

Vertex* AdaptiveRemesher::collapse_edge(Edge * edge)
{
	HalfEdge *he;
	Mesh::ring_iterator it;
	const double rat = 2;

	he= edge->he;
	if(he->face != Mesh::hole())
	{
		const double l1 = arma::norm(he->vert->xyz - he->next->vert->xyz); he=he->next;
		const double l2 = arma::norm(he->vert->xyz - he->next->vert->xyz); he=he->next;
		const double l3 = arma::norm(he->vert->xyz - he->next->vert->xyz); he=he->next;
		if( (l1 > l2*rat) || (l1 > l3*rat) ) return NULL;
	}
	he = he->twin;
	{
		const double l1 = arma::norm(he->vert->xyz - he->next->vert->xyz); he=he->next;
		const double l2 = arma::norm(he->vert->xyz - he->next->vert->xyz); he=he->next;
		const double l3 = arma::norm(he->vert->xyz - he->next->vert->xyz); he=he->next;
		if( (l1 > l2*rat) || (l1 > l3*rat) ) return NULL;
	}

	/*
	 * Find the new location
	 */
	Vertex *v = new Vertex;
	Vector3 abc;
	if(he->face->he == he)                  abc = "0.5 0.5 0";
	else if(he->face->he == he->next)       abc = "0.5 0 0.5";
	else if(he->face->he == he->next->next) abc = "0 0.5 0.5";
	else
	{
		printf("[Error]");
		throw;
	}
	proj()->find_point(he->face, abc, v);

	/*
	 * Call fidelity
	 * not needed for now
	 */
	Vertex *cands[]= {v, NULL, NULL};
	double Q1, Q2;

	it.init(edge->he->vert);
	do
	{
		if( (it.half_edge()->face == edge->he->face) || (it.half_edge()->face == edge->he->twin->face) ) continue;
		cands[1] = it.half_edge()->next->next->vert;
		cands[2] = it.half_edge()->vert;
		if(!this->tri_fidelity(cands, Q1, Q2)){ delete v; return NULL; }
	}while(it.advance());
	it.init(edge->he->twin->vert);
	do
	{
		if( (it.half_edge()->face == edge->he->face) || (it.half_edge()->face == edge->he->twin->face) ) continue;
		cands[1] = it.half_edge()->next->next->vert;
		cands[2] = it.half_edge()->vert;
		if(!this->tri_fidelity(cands, Q1, Q2)){ delete v; return NULL; }
	}while(it.advance());

	/*
	 * check if we can insert this
	 */
	if( mesh()->collapse_edge(edge, v) ) return v;
	else
	{
		delete v;
		return NULL;
	}
}

void AdaptiveRemesher::geodesic_edge(Edge* edge, Vector3 xyz_m[4]) const
{
	// Three dimensional coordinates
	Vertex *v0 = edge->he->vert;
	Vertex *v1 = edge->he->twin->vert;
	Vertex *v2 = edge->he->next->next->vert;
	Vertex *v3 = edge->he->twin->next->next->vert;

	// Find the angles and distances
	Vector3 e01 = v1->xyz - v0->xyz;
	Vector3 e02 = v2->xyz - v0->xyz;
	Vector3 e03 = v3->xyz - v0->xyz;

	const double l0 = arma::norm(e02); e02 /= l0;
	const double l1 = arma::norm(e03); e03 /= l1;
	const double l2 = arma::norm(e01); e01 /= l2;

	const double c0 = arma::dot(e01, e02);
	const double c1 = arma::dot(e01, e03);
	const double s0 = sqrt(1. - c0*c0); assert(finite(s0));
	const double s1 = sqrt(1. - c1*c1); assert(finite(s1));

	xyz_m[0](0) = 0;     xyz_m[0](1) = 0;      xyz_m[0](2) = 0;
	xyz_m[1](0) = l2;    xyz_m[1](1) = 0;      xyz_m[1](2) = 0;
	xyz_m[2](0) = l0*c0; xyz_m[2](1) = l0*s0;  xyz_m[2](2) = 0;
	xyz_m[3](0) = l1*c1; xyz_m[3](1) = -l1*s1; xyz_m[3](2) = 0;
	//xyz_m[0] = v0->xyz;
	//xyz_m[1] = v1->xyz;
	//xyz_m[2] = v2->xyz;
	//xyz_m[3] = v3->xyz;
}

void AdaptiveRemesher::geodesic_vert(Vertex* v, Mesh& mout) const
{
	mout.clear_memory();

	std::vector<double> len, theta;
	std::vector<HalfEdge*> hes3d;
	double theta_tot, alpha;
	int ivprev, ivcurr;
	Mesh::ring_iterator it;
	Vertex *vcurr, *vprev;
	Vector3 ecurr, eprev;


	/*
	 * Get the half edges
	 */
	it.init(v);
	do
	{
		hes3d.push_back(it.half_edge());
	}while(it.advance());

	/*
	 * Create theta and len
	 */
	theta_tot = 0;
	for(int i = 0 ; i < int(hes3d.size()) ; i++)
	{
		HalfEdge *he = hes3d[int(hes3d.size()) - i - 1];
		vcurr = he->vert;
		vprev = he->next->next->vert;
		ecurr = vcurr->xyz - v->xyz;
		eprev = arma::normalise(vprev->xyz - v->xyz);
		len.push_back(arma::norm(ecurr));
		ecurr /= len.back();
		theta.push_back( acos_tol(arma::dot(ecurr, eprev)) );
		theta_tot += theta.back();

		// printf("vcurr %d, vprev %d \n", vcurr->idx(), vprev->idx() );

	}while(it.advance());
	alpha =  2. * M_PI / theta_tot;

	// printf("alpha is %lf \n", alpha);

	/*
	 * Create the mesh
	 */
	std::vector<int> conn;
	std::vector<double> coord(3, 0.);
	theta_tot = 0;
	for (int i = 0 ; i < int(hes3d.size()) ; i++)
	{
		HalfEdge *he = hes3d[int(hes3d.size()) - i - 1];
		theta_tot += alpha * theta[i];
		ivprev = 1 +  ( i == 0 ? int(len.size()) - 1 : i -1);
		ivcurr = 1 + i;


		coord.push_back( len[i] * cos(theta_tot) );
		coord.push_back( len[i] * sin(theta_tot) );
		coord.push_back(0.);

		if( he->face->he == he)
		{
			conn.push_back(ivcurr);
			conn.push_back(0);
			conn.push_back(ivprev);
		}
		else if( he->face->he == he->next )
		{
			conn.push_back(0);
			conn.push_back(ivprev);
			conn.push_back(ivcurr);
		}
		else
		{
			assert( he->face->he == he->next->next );
			conn.push_back(ivprev);
			conn.push_back(ivcurr);
			conn.push_back(0);
		}
	}
	mout.init(coord, conn);

	/*
	 * Set the backward mapping
	 */
	for (int iface = 0 ; iface < mout.faces.n_mems() ; iface++)
	{
		Face *face2d = mout.faces.get_member_ptr(iface);
		face2d->misc.resize(1);
		face2d->misc[0] = hes3d[mout.faces.n_mems() - iface - 1]->face ;

		// printf(" Face2d %d maps to face3d %d \n", face2d->idx(), static_cast<Face*>(face2d->misc[0])->idx() );
	}

}


/************************************************************************\
                              Projection
\************************************************************************/
Projector* Projector::build(GetPot &getpot, Mesh &mesh)
{
	 std::string type = getpot("proj_type", "patch");

	 if(type == "plane")
	 {
		 printf("[log] projector: plane. \n");
		 return new PlaneProjector(mesh);
	 }
	 if(type == "sphere")
	 {
		 printf("[log] projector: sphere. \n");
		 std::string scale_str = getpot("proj_scale", "1.0 1.0 1.0");
		 Vector3 scale_vec(scale_str);
		 return new SphereProjector(mesh, scale_vec);
	 }
	 if(type == "current")
	 {
		 printf("[log] projector: current mesh. \n");
		 return new CurrentProjector(mesh);
	 }
	 if(type == "patch")
	 {
		 printf("[log] projector: overlapping patches. \n");
		 return new PatchProjector(mesh);
	 }

	 else
	 {
		 printf("[Error] projector %s not found. \n", type.c_str());
		 assert(0); throw;
	 }
}

void PlaneProjector::find_point(Face* face, const Vector3& abc, Vertex *out)
{
	face->abc2xy(abc, out->xyz);
	out->normal="0 0 1";
}

void PlaneProjector::init()
{
	/*
	 * Find the normals at each vertex
	 */
	for (int i = 0 ; i < _mesh.verts.n_mems() ; i++)
	{
		Vertex *v0 = _mesh.verts.get_member_ptr(i);
		v0->normal = "1 0 0";
	}
}

SphereProjector::SphereProjector(Mesh &mesh_in, const Vector3 scale_in):
		Projector(mesh_in), scale(scale_in)
{
	trans.zeros();
	trans(0,0) = scale(0);
	trans(1,1) = scale(1);
	trans(2,2) = scale(2);
	//
	if(arma::det(trans) < 1e-5)
	{
		printf("[Error] uninvertible transformation matrix for sphere. \n");
		std::cout << trans << std::endl;
		assert(0); throw;
	}
	i_trans = trans.i();
	//
	n_trans = i_trans.t();

}

void SphereProjector::init()
{
	/* find the normals */
	for (int iv = 0 ; iv < _mesh.verts.n_mems() ; iv++)
	{
		Vertex *v = _mesh.verts.get_member_ptr(iv);
		v->normal = i_trans * v->xyz;
		arma::normalise(v->normal);
		v->normal = n_trans * v->normal;
		arma::normalise(v->normal);
	}
}

void SphereProjector::find_point(Face* face, const Vector3& abc, Vertex *out)
{
	face->abc2xy(abc, out->xyz);
	// Project the vector to reference coordinate system
	out->xyz = i_trans * out->xyz;
	// Project the sphere surface
	arma::normalise(out->xyz);
	out->normal = out->xyz;
	// Bring back the point to the physical space
	out->xyz = trans * out->xyz;
	out->normal = n_trans * out->normal;
	arma::normalise(out->normal);
}

CurrentProjector::CurrentProjector(Mesh& mesh_in):Projector(mesh_in)
{}

void CurrentProjector::init()
{
	/*
	 * Find the normals at each vertex
	 */
	Mesh::ring_iterator it;
	Vertex *v0,*v1,*v2;
	for (int i = 0 ; i < _mesh.verts.n_mems() ; i++)
	{
		v0 = _mesh.verts.get_member_ptr(i);
		v0->normal.zeros();
		it.init(v0);
		do
		{
			v1 = it.half_edge()->next->next->vert;
			v2 = it.half_edge()->vert;
			v0->normal += arma::cross((v1->xyz - v0->xyz), (v2->xyz - v0->xyz));
		}while(it.advance());
		v0->normal = arma::normalise( v0->normal );
	}
}


void CurrentProjector::find_point(Face* face, const Vector3& abc, Vertex *out)
{
	face->abc2xy(abc, out->xyz);
	Vertex *fvs[] = {face->he->vert, face->he->next->vert, face->he->next->next->vert};
	out->normal = arma::normalise( abc(0)*fvs[0]->normal + abc(1)*fvs[1]->normal + abc(2)*fvs[2]->normal );

}

