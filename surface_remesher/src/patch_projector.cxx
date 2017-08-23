/*
 * patch_projector.cxx
 *
 *  Created on: Dec 14, 2016
 *      Author: hooshi
 */

#include <queue>
#include <stack>
#include <sstream>
#include <map>
#include <iostream>

#include "cgal_tools.hxx"
#include "adaptive_remesher.hxx"
#include "mesh_io.hxx"

using namespace std;

#define CHECK_ABC(abc) assert( fabs( 1 - abc(2) - abc(1) - abc(0) ) < 1e-10 )

PatchProjector::PatchProjector(Mesh& mesh_in): Projector(mesh_in){}
PatchProjector::~PatchProjector()
{
	for(auto it = patches.begin() ; it != patches.end() ; ++it)
	{
		delete *it;
	}
}

void PatchProjector::init()
{
	/*
	 * Create the reference mesh
	 */
	std::vector<int> conn;
	std::vector<double> coords;

	for (int iv = 0 ; iv < _mesh.verts.n_mems() ; iv++)
	{
		Vertex *v = _mesh.verts.get_member_ptr(iv);
		coords.push_back( v->xyz(0) );
		coords.push_back( v->xyz(1) );
		coords.push_back( v->xyz(2) );
	}

	for (int it = 0 ; it < _mesh.faces.n_mems() ; it++)
	{
		Face *t = _mesh.faces.get_member_ptr(it);
		conn.push_back( t->he->vert->idx() );
		conn.push_back( t->he->next->vert->idx() );
		conn.push_back( t->he->next->next->vert->idx() );
	}

	mesh_ref.init(coords, conn);
	assert(mesh_ref.faces.n_mems() == _mesh.faces.n_mems());
	assert(mesh_ref.verts.n_mems() == _mesh.verts.n_mems());
	assert(mesh_ref.half_edges.n_mems() == _mesh.half_edges.n_mems());
	assert(mesh_ref.edges.n_mems() == _mesh.edges.n_mems());

	/*
	 * Find the normal at each point on the reference mesh and
	 * current mesh.
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
		mesh_ref.verts.get_member_ptr(i)->normal = v0->normal;
	}

	/*
	 * Set the reference data for the current mesh vertex
	 */
	for (int it = 0 ; it < _mesh.faces.n_mems() ; it++)
	{
		Face *tr = mesh_ref.faces.get_member_ptr(it);
		Vertex *vr[3] = {tr->he->vert, tr->he->next->vert, tr->he->next->next->vert};

		Face *tm = _mesh.faces.get_member_ptr(it);
		Vertex *vm[3] = {tm->he->vert, tm->he->next->vert, tm->he->next->next->vert};

		assert( vr[0]->idx() == vm[0]->idx() );
		assert( vr[1]->idx() == vm[1]->idx() );
		assert( vr[2]->idx() == vm[2]->idx() );

		vm[0]->facer = vm[1]->facer = vm[2]->facer = tr;
		vm[0]->abcr = "1 0 0";
		vm[1]->abcr = "0 1 0";
		vm[2]->abcr = "0 0 1";
	}

}


void PatchProjector::find_point(Face* face, const Vector3& abc, Vertex *out)
{
	static int n_visit=0;

	Vertex *fvs[] = {face->he->vert, face->he->next->vert, face->he->next->next->vert};
	Mesh *patch;
	Face *pfaces[3] = {NULL, NULL, NULL};
	std::list< std::pair<Mesh*, Face*> >::iterator pit[3];

	/*
	 * Special case all three points are on the same facer;
	 */
	if( (fvs[0]->facer == fvs[1]->facer) && (fvs[1]->facer == fvs[2]->facer) )
	{
		out->facer = fvs[0]->facer;
		Vertex *rfvs[] = {out->facer->he->vert, out->facer->he->next->vert, out->facer->he->next->next->vert};

		out->abcr = abc(0) * fvs[0]->abcr + abc(1) * fvs[1]->abcr + abc(2) * fvs[2]->abcr;
		CHECK_ABC(out->abcr);
		out->facer->abc2xy(out->abcr, out->xyz);
		out->normal = out->abcr(0)*rfvs[0]->normal + out->abcr(1)*rfvs[1]->normal  + out->abcr(2)*rfvs[2]->normal;
		//printf("found from past. \n");
		return ;
	}

	/*
	 *  Search for a common patch
	 */
	stringstream ss;
	n_visit++;
	// printf("n_visit: %d \n", n_visit);
	for (pit[0] = fvs[0]->facer->pinfo.begin() ; pit[0] != fvs[0]->facer->pinfo.end() ; pit[0]++)
	{
		patch = pit[0]->first;
		pfaces[0] = pit[0]->second;

		// search first patch
		//pit[1] = fvs[1]->facer->pinfo.find(*pit[0]);
		//pit[2] = fvs[2]->facer->pinfo.find(*pit[0]);
		//if( ( pit[1] != fvs[1]->facer->pinfo.end() ) && ( pit[2] != fvs[2]->facer->pinfo.end() ) )
		//	{
		//		pfaces[1] = pit[1]->second; assert(pit[1]->first == patch);
		//		pfaces[2] = pit[2]->second; assert(pit[2]->first == patch);
		//      printf("found patch \n");
		//		break;
		//		printf("found patch \n");
		//	}

		// search first patch
		// MeshIO(*pit[0]->first).write_vtk("hello.vtk");
		bool success = false;
		for( pit[1] = fvs[1]->facer->pinfo.begin() ; pit[1] != fvs[1]->facer->pinfo.end() ; pit[1]++ )
		{
			if(pit[0]->first == pit[1]->first)
			{
				for( pit[2] = fvs[2]->facer->pinfo.begin() ; pit[2] != fvs[2]->facer->pinfo.end() ; pit[2]++ )
				{
					if(pit[0]->first == pit[2]->first)
					{
						pfaces[1] = pit[1]->second; assert(pit[1]->first == patch);
						pfaces[2] = pit[2]->second; assert(pit[2]->first == patch);
						// printf("found patch \n");
						success = true;
						break;
					}
				}
				if(success) break;
			}
		} // End of search over 1 and 2
		if(success) break;

	} // End of total search

	/*
	 * Create the patch if it was not found.
	 */
	if (!pfaces[1])
	{
		assert(!pfaces[2]);

		//printf("not found patch. \n");
		this->create_patch(face, patch, pfaces);
	}

	/*
	 * some checks
	 */
	assert( patch->faces.is_member(pfaces[0]) );
	assert( patch->faces.is_member(pfaces[1]) );
	assert( patch->faces.is_member(pfaces[2]) );

	/*
	 * find the location of the point in the patch
	 */
	bool scss;
	Vector3 pxyz[3], pxyztarget;
	Face *pfacetarget;
	pfaces[0]->abc2xy(fvs[0]->abcr, pxyz[0]);
	pfaces[1]->abc2xy(fvs[1]->abcr, pxyz[1]);
	pfaces[2]->abc2xy(fvs[2]->abcr, pxyz[2]);
	pxyztarget = abc(0)*pxyz[0] + abc(1)*pxyz[1] + abc(2)*pxyz[2];
	pfacetarget = pfaces[0];
	scss = patch->locate_point(pxyztarget, pfacetarget, out->abcr);
	if( !scss )
	{
		MeshIO(*patch).write_vtk("badpatch.vtk");
		printf("[Error] could not find point in patch: ");
		std::cout << pxyztarget.t();
		printf("[Error] reference triangles are: %d %d %d \n",
				pfaces[0]->idx(), pfaces[1]->idx(), pfaces[2]->idx());
		printf("[Error] patch members: %d \n", patch->faces.n_mems());
		throw;
	}
	CHECK_ABC(out->abcr);

	/*
	 * Find the location of the point in the original mesh
	 */
	out->facer = static_cast<Face*>(pfacetarget->misc[0]);
	Vertex *rfv[]={out->facer->he->vert, out->facer->he->next->vert, out->facer->he->next->next->vert};
	out->facer->abc2xy(out->abcr, out->xyz);
	out->normal = out->abcr(0)*rfv[0]->normal + out->abcr(1)*rfv[1]->normal + out->abcr(2)*rfv[2]->normal;
}


static void bfs_inspect(std::set<Face*>& members, std::set<Face*>& bfspool, std::queue<Face*>& bfsfront, Face *source)
{

	Face *nextface[3] = {source->he->twin->face, source->he->next->twin->face, source->he->next->next->twin->face};
	bool topology;

	for (int i = 0 ; i < 3 ; i++)
	{
		// printf("next face: %d \n", nextface[i]->idx() );
		if(nextface[i] != Mesh::hole() && ( !bfspool.count(nextface[i]) ))
		{
			/*
			 * Check topology
			 */
			// find number of edges inside
			Face *neighs[] = {nextface[i]->he->twin->face, nextface[i]->he->next->twin->face, nextface[i]->he->next->next->twin->face};
			Vertex *verts[] = {nextface[i]->he->next->next->vert, nextface[i]->he->vert, nextface[i]->he->next->vert};
			int cnt[] = { int(members.count(neighs[0])), int(members.count(neighs[1])), int(members.count(neighs[2])) };
			Vertex *suspect = NULL;
			Mesh::ring_iterator it;
			topology = true;

			if (cnt[0] + cnt[1] + cnt[2] != 2)
			{
				if(cnt[0]) suspect = verts[0];
				else if(cnt[1]) suspect = verts[1];
				else if(cnt[2]) suspect = verts[2];
				else
				{
					fprintf(stderr, "[Error] topology check in BFS search failed. \n");
					assert(0);
					throw;
				}

				it.init(suspect);
				do
				{
					if( members.count( it.half_edge()->face ) )
					{
						topology = false;
						// printf("topology failure ! \n");
						break;
					}
				}while(it.advance());
			}

			/*
			 * Add to front, pool and members
			 */
			if( topology )
			{
				bfspool.insert(nextface[i]);
				members.insert(nextface[i]);
				bfsfront.push(nextface[i]);
			}
			//	printf("Added to pool \n");
		}
		else
		{
			//	printf("Already in pool\n");
		}
	} // End of for over neighbours
}

int PatchProjector::bfs_search(std::set<Face*>& members, Face *source, Face *targ1, Face *targ2, int target_num)
{
	std::queue<Face*> bfsfront;
	std::set<Face*> bfspool;
	Face *nextface;

	members.insert( source );
	bfspool.insert( source );
	bfsfront.push( source );

	/*
	 * Do the bfs search
	 */
	//int i =0 ;
	while( !bfsfront.empty() )
	{
		//i ++;
		nextface = bfsfront.front(); 	bfsfront.pop();
		bfs_inspect(members, bfspool, bfsfront, nextface);

		/*
		 * debug
		std::vector<int> bfs1(mesh_ref.faces.n_mems(), 0);
		FILE *fl;
		MeshIO io(mesh_ref);
		stringstream ss;
		ss.str("");
		ss << "bfsstep_" << i << ".vtk";
		fl = fopen(ss.str().c_str(), "w");
		io.write_vtk(fl);
		for (auto it = bfspool.begin() ; it != bfspool.end() ; ++it)
		{
			bfs1[(*it)->idx()] = 1;
		}
		bfs1[source->idx()] = bfs1[targ1->idx()] = bfs1[targ2->idx()] = 2;
		io.write_vtk_data(fl, bfs1, "bfs");
		fclose(fl);
		*/

		// ending criteria
		if ( (target_num < 0) && bfspool.count(targ1) && bfspool.count(targ2) ) break;
		else if ( (target_num > 0) && (int(bfspool.size()) > target_num) ) break;
	}

	return int(bfspool.size());
}

bool PatchProjector::create_patch(Face* mface, Mesh *& patch, Face *pfaces[3])
{

	/*
	 * Get the verts of the face
	 */
	Vertex *mverts[3] = {mface->he->vert, mface->he->next->vert, mface->he->next->next->vert};

	/*
	 * Create a mesh
	 */
	patch = new Mesh;
	patches.push_back(patch);

	/*
	 * BFS data
	 */
	std::set<Face*> members;
	std::queue< Face* > ears;

	/*
	 * Do the bfs search
	 */

	int n_mem1;
	n_mem1 = this->bfs_search(members, mverts[0]->facer, mverts[1]->facer, mverts[2]->facer, -1);
	this->bfs_search(members, mverts[1]->facer, mverts[2]->facer, mverts[0]->facer, n_mem1 );
	this->bfs_search(members, mverts[2]->facer, mverts[0]->facer, mverts[1]->facer, n_mem1 );

	//if( members.size() > 500 )
	//{
	//	std::vector<int> bfs1(mesh_ref.faces.n_mems(), 0);
	//	FILE *fl;
	//	MeshIO io(mesh_ref);
	//	fl = fopen("bfs0.vtk", "w");
	//	io.write_vtk(fl);
	//	for (auto it = members.begin() ; it != members.end() ; ++it)
	//	{
	//		bfs1[(*it)->idx()] = 1;
	//	}
	//	bfs1[mverts[0]->facer->idx()] = bfs1[mverts[1]->facer->idx()] = bfs1[mverts[2]->facer->idx()] = 2;
	//	io.write_vtk_data(fl, bfs1, "bfs1");
	//	printf("[Error] bfs search failed, see bfs.vtk! \n");
	//	throw;
	//}

	/*
	 * Trim the edges
	 */
	for (auto it = members.begin() ; it != members.end() ; ++it)
	{
		Face *fc = *it;
		assert(fc != Mesh::hole());
		Face *neighs[] = {fc->he->twin->face, fc->he->next->twin->face, fc->he->next->next->twin->face};
		int cnt[] = { int(members.count(neighs[0])), int(members.count(neighs[1])), int(members.count(neighs[2])) };

		if( cnt[0] + cnt[1] + cnt[2] == 1) ears.push(fc);
	}

	while(!ears.empty())
	{
		Face *fc = ears.front(); ears.pop();
		assert(fc != Mesh::hole());
		Face *neighs[] = {fc->he->twin->face, fc->he->next->twin->face, fc->he->next->next->twin->face};

		auto itear = members.find(fc);
		int cnt[] = { int(members.count(neighs[0])), int(members.count(neighs[1])), int(members.count(neighs[2])) };

		if( (fc != mverts[0]->facer) && (fc != mverts[1]->facer) && (fc != mverts[2]->facer)
				&& (itear != members.end() ) && (cnt[0] + cnt[1] + cnt[2] == 1) )
		{
			members.erase(itear);
			if(neighs[0] != Mesh::hole()) ears.push(neighs[0]);
			if(neighs[1] != Mesh::hole()) ears.push(neighs[1]);
			if(neighs[2] != Mesh::hole()) ears.push(neighs[2]);
		}
	}


	/*
	bfs1.resize(0);
	bfs1.resize(mesh_ref.faces.n_mems(), 0);
	for (auto it = members.begin() ; it != members.end() ; ++it)
	{
		bfs1[(*it)->idx()] = 1;
	}
	bfs1[mverts[0]->facer->idx()] = bfs1[mverts[1]->facer->idx()] = bfs1[mverts[2]->facer->idx()] = 2;
	io.write_vtk_data(fl, bfs1, "bfs2");
	fclose(fl);
	*/

	/*
	 * Create a two dimensional mesh
	 */
	vector<int> conn;
	vector<double> coords;
	map<int, int> face_gtol, vert_gtol;
	int n_vertex = 0, n_face = 0;

	// create the maps
	for (auto it = members.begin() ; it != members.end() ; ++it)
	{
		Face *fc = *it;
		Vertex *verts[] = {fc->he->vert, fc->he->next->vert, fc->he->next->next->vert};
		if (!vert_gtol.count(verts[0]->idx())) { vert_gtol.insert( make_pair(verts[0]->idx(), n_vertex) ); n_vertex++; }
		if (!vert_gtol.count(verts[1]->idx())) { vert_gtol.insert( make_pair(verts[1]->idx(), n_vertex) ); n_vertex++; }
		if (!vert_gtol.count(verts[2]->idx())) { vert_gtol.insert( make_pair(verts[2]->idx(), n_vertex) ); n_vertex++; }
		face_gtol.insert( make_pair(fc->idx(), n_face) ); n_face++;
	}

	// create vertex coordinates
	coords.resize( vert_gtol.size() * 3 );
	for (auto it = vert_gtol.begin(); it != vert_gtol.end(); ++it)
	{
		Vertex *v = mesh_ref.verts.get_member_ptr(it->first);
		coords[it->second*3 + 0] =  v->xyz(0);
		coords[it->second*3 + 1] = v->xyz(1) ;
		coords[it->second*3 + 2] = v->xyz(2) ;
	}

	// create the connectivity
	conn.resize( face_gtol.size() * 3);
	for (auto it = face_gtol.begin(); it != face_gtol.end(); ++it)
	{
		Face *t = mesh_ref.faces.get_member_ptr(it->first);
		conn[3*it->second+0] = vert_gtol[ t->he->vert->idx() ] ;
		conn[3*it->second+1] = vert_gtol[ t->he->next->vert->idx() ] ;
		conn[3*it->second+2] = vert_gtol[ t->he->next->next->vert->idx()  ];
	}

	// create the mesh
	patch->init(coords, conn);


	/*
	 * Set the mapping right
	 */
	// face to face
	for (auto it = face_gtol.begin(); it != face_gtol.end(); ++it)
	{
		Face *fref = mesh_ref.faces.get_member_ptr(it->first);
		Face *fpatch = patch->faces.get_member_ptr(it->second);

		fpatch->misc.push_back(fref);
		fref->pinfo.push_back( make_pair(patch, fpatch) );

		assert ( patch->faces.is_member(fpatch) );
		assert ( mesh_ref.faces.is_member( static_cast<Face*>(fpatch->misc[0]) ) );
		assert( static_cast<Face*>(fpatch->misc[0])->pinfo.back().second == fpatch );
		assert( static_cast<Face*>(fpatch->misc[0])->pinfo.back().first == patch );
	}

	pfaces[0] = patch->faces.get_member_ptr( face_gtol[ mverts[0]->facer->idx() ] );
	pfaces[1] = patch->faces.get_member_ptr( face_gtol[ mverts[1]->facer->idx() ] );
	pfaces[2] = patch->faces.get_member_ptr( face_gtol[ mverts[2]->facer->idx() ] );

	assert( static_cast<Face*>(pfaces[0]->misc[0])->pinfo.back().second == pfaces[0] );
	assert( static_cast<Face*>(pfaces[1]->misc[0])->pinfo.back().second == pfaces[1] );
	assert( static_cast<Face*>(pfaces[2]->misc[0])->pinfo.back().second == pfaces[2] );
	assert( static_cast<Face*>(pfaces[0]->misc[0])->pinfo.back().first == patch );
	assert( static_cast<Face*>(pfaces[1]->misc[0])->pinfo.back().first == patch );
	assert( static_cast<Face*>(pfaces[2]->misc[0])->pinfo.back().first == patch );
	assert ( patch->faces.is_member(pfaces[0]) );
	assert ( patch->faces.is_member(pfaces[1]) );
	assert ( patch->faces.is_member(pfaces[2]) );

	//	std::vector<int> bfs2(patch->faces.n_mems(), 0);
	//	MeshIO io2(*patch);
	//	FILE *fl2 = fopen("patch0.vtk", "w");
	//	io2.write_vtk(fl2);
	//	for (auto it = members.begin() ; it != members.end() ; ++it)
	//	{
	//		bfs2[ face_gtol [(*it)->idx()] ] = 1;
	//	}
	//	bfs2[ face_gtol[ mverts[0]->facer->idx()] ] = bfs2[ face_gtol[ mverts[1]->facer->idx()] ] = bfs2[ face_gtol[ mverts[2]->facer->idx()] ] = 2;
	//	io2.write_vtk_data(fl2, bfs2, "bfs");

	/*
	 * Map the patch into a unit circle.
	 */
	string errmsg; bool scss;
	scss = CGAL_circular_map(*patch, errmsg);
	if(!scss) std::cerr << "[Mapping Error] " << errmsg << std::endl;

	//	fl2 = fopen("patch1.vtk", "w");
	//	io2.write_vtk(fl2);
	//	io2.write_vtk_data(fl2, bfs2, "bfs");

	return scss;
}

