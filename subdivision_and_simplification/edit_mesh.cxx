/* Copyright (c) Darcy Harisson, Russell Gillette
 * April 2014
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use, copy,
 * modify, merge, publish, distribute, sublicense, and/or sell copies
 * of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include "edit_mesh.hxx"
#include "mesh_io.hxx"

#include <Eigen/Geometry>
#include <cmath>

#include <map>
#include <set>
#include <iostream>
#include <memory>
#include <fstream>

namespace hooshi {
  
  using std::multimap;
  using std::map;
  using std::vector;
  using std::cout;
  using std::endl;
  using std::size_t;
  using std::pair;
  using std::make_pair;

  //#if defined(NDEBUG) && defined(ALWAYS_ASSERT)
  //#undef NDEBUG
  //#endif
  //#include <cassert>

  EditMesh *loadEditMeshFromFile(std::string file_name);

  namespace detail{
    inline void init( HalfEdge& he, std::size_t next, std::size_t twin, std::size_t vert, std::size_t face )
    {
      he.next = next;
      he.twin = twin;
      he.vert = vert;
      he.face = face;
    }

    void delete_face( std::vector<std::size_t>& faceData, std::vector<HalfEdge>& heData, std::size_t f )
    {
      assert( f < faceData.size() );

      // In order to delete the face properly, we need to move a face from the end of the list to overwrite 'f'. Then we need to update the 
      // indices stored in the moved face's half-edges.
      faceData[f] = faceData.back();
      faceData.pop_back();

      if( f != faceData.size() ){
	//std::clog << "Reindexed face " << faceData.size() << " to " << f << std::endl;

	std::size_t he = faceData[f];
	do {
	  assert( heData[he].face == faceData.size() );

	  heData[he].face = f;
	  he = heData[he].next;
	} while( he != faceData[f] );
      }
    }

    template <int N>
    void delete_faces( std::vector<std::size_t>& faceData, std::vector<HalfEdge>& heData, std::size_t (&fToDelete)[N] )
    {
      // Sort the faces by decreasing index so that we can safely delete them all without causing any of them to be accidentally re-indexed (which
      // cause 'fToDelete' to contain invalid indices). This also chooses the optimal deletion order to minimize re-indexing.
      std::sort( fToDelete, fToDelete + N, std::greater<std::size_t>() );
      for( std::size_t i = 0; i < N; ++i )
	detail::delete_face( faceData, heData, fToDelete[i] );
    }
  }

  void init_adjacency( std::size_t numVertices, const std::vector<std::size_t>& faces, std::vector< HalfEdge >& _he_data, std::vector< std::size_t >& _face_to_he, std::vector< std::size_t >& _vert_to_he ){
    typedef std::map< std::pair<std::size_t, std::size_t>, std::size_t > edge_map_type;
	
    assert( faces.size() % 3 == 0 && "Invalid data specified for faces. Must have 3 vertex indices per face." );

    edge_map_type edgeMap; // Use a temporary map to find edge pairs.

    _he_data.reserve( faces.size() ); // Assume there are 3 edges per face.
    _face_to_he.resize( faces.size() / 3 );
    _vert_to_he.resize( numVertices, HOLE_INDEX ); // Init with HOLE_INDEX since a vert might be floating w/ no faces.

    for( std::size_t i = 0, iEnd = faces.size(); i < iEnd; i+=3 ){
      std::size_t f[] = { faces[i], faces[i+1], faces[i+2] };
      std::size_t fIndex = i / 3;

      // The index of the first (of three) half-edges associated with the current face.
      std::size_t heIndex = _he_data.size();

      HalfEdge he[3];
      detail::init( he[0], heIndex+1, HOLE_INDEX, f[0], fIndex );
      detail::init( he[1], heIndex+2, HOLE_INDEX, f[1], fIndex );
      detail::init( he[2], heIndex, HOLE_INDEX, f[2], fIndex );
#ifdef USE_PREV
      he[0].prev = heIndex+2;
      he[1].prev = heIndex;
      he[2].prev = heIndex+1;
#endif
			
      // These will be set each time a vertex is referenced, but that's fine. The last assignment will stick.
      _face_to_he[ fIndex ] = heIndex;
      _vert_to_he[ f[0] ] = heIndex;
      _vert_to_he[ f[1] ] = heIndex+1;
      _vert_to_he[ f[2] ] = heIndex+2;

      edge_map_type::iterator it;

      it = edgeMap.lower_bound( std::make_pair( f[0], f[1] ) );
      if( it != edgeMap.end() && it->first.first == f[0] && it->first.second == f[1] ){
	_he_data[it->second].twin = heIndex;
	he[0].twin = it->second;
	edgeMap.erase( it );
      } else {
	he[0].twin = HOLE_INDEX;
	edgeMap.insert( it, std::make_pair( std::make_pair( f[1], f[0] ), heIndex ) ); // NOTE: Reversed order since we are matching opposite HalfEdge.
      }

      it = edgeMap.lower_bound( std::make_pair( f[1], f[2] ) );
      if( it != edgeMap.end() && it->first.first == f[1] && it->first.second == f[2] ){
	_he_data[it->second].twin = heIndex+1;
	he[1].twin = it->second;
	edgeMap.erase( it );
      } else {
	he[1].twin = HOLE_INDEX;
	edgeMap.insert( it, std::make_pair( std::make_pair( f[2], f[1] ), heIndex+1 ) ); // NOTE: Reversed order since we are matching opposite HalfEdge.
      }

      it = edgeMap.lower_bound( std::make_pair( f[2], f[0] ) );
      if( it != edgeMap.end() && it->first.first == f[2] && it->first.second == f[0] ){
	_he_data[it->second].twin = heIndex+2;
	he[2].twin = it->second;
	edgeMap.erase( it );
      } else {
	he[2].twin = HOLE_INDEX;
	edgeMap.insert( it, std::make_pair( std::make_pair( f[0], f[2] ), heIndex+2 ) ); // NOTE: Reversed order since we are matching opposite HalfEdge.
      }

      _he_data.push_back( he[0] );
      _he_data.push_back( he[1] );
      _he_data.push_back( he[2] );
    }

    // Keep track of the last edge we processed so we can hook up HalfEdge::prev as we go.
    // const std::size_t prev = HOLE_INDEX;

    // Add half-edges for any holes. Any edges still in the map are holes.
    edge_map_type::iterator it = edgeMap.begin();
    while( it != edgeMap.end() ){
      HalfEdge he;
      detail::init( he, HOLE_INDEX, it->second, it->first.first, HOLE_INDEX );
#ifdef USE_PREV
      he.prev = prev;
      prev = _he_data.size(); // Size is the index of the HalfEdge we are about to push into the list.
#endif

      _he_data[he.twin].twin = _he_data.size();
      _he_data.push_back( he );

      // const std::size_t curVert = it->first.first;
      std::size_t nextVert = it->first.second; // We are about to erase this information, so store it to use later.

      edgeMap.erase( it ); // We are done with this edge now.

      HalfEdge* twinPrev = &_he_data[_he_data[_he_data[he.twin].next].next];
      while( twinPrev->twin != HOLE_INDEX && _he_data[twinPrev->twin].face != HOLE_INDEX ){
	assert( _he_data[twinPrev->next].vert == nextVert );
	assert( _he_data[twinPrev->twin].vert == nextVert );
	twinPrev = &_he_data[_he_data[_he_data[twinPrev->twin].next].next];
      }

      if( twinPrev->twin == HOLE_INDEX ){
	// We haven't processed the next edge in the loop yet. Let's do so now so we can assume the index of the next half-edge.
	_he_data.back().next = _he_data.size();
	it = edgeMap.find( std::make_pair( nextVert, twinPrev->vert ) );
				
	assert( it != edgeMap.end() );
      }else{
	assert( _he_data[twinPrev->twin].vert == nextVert );
	assert( _he_data[twinPrev->twin].face == HOLE_INDEX );

	// We already processed this edge and have a valid index for the next HalfEdge.
	_he_data.back().next = twinPrev->twin;
#ifdef USE_PREV
	_he_data[ twinPrev->twin ].prev = prev; // Complete the loop
	prev = HOLE_INDEX;
#endif
	it = edgeMap.begin(); // Arbitrarily pick the next edge in the list.
      }
    }

    assert( edgeMap.empty() );
  }

  void EditMesh::init( const std::vector<double>& xyzPositions, const std::vector<std::size_t>& triangleVerts, const bool is_reinit ){

    // first clear the mesh.
    if(is_reinit) prepare_for_reinit();
    else clear();
	
    assert( xyzPositions.size() % 3 == 0 && "Invalid vertex positions for EditMesh::init(). Must have 3 values per-vertex." );
    assert( triangleVerts.size() % 3 == 0 && "Invalid face data for EditMesh::init(). Must have 3 vertex indices per face." );

    //_verts.resize( Eigen::NoChange, xyzPositions.size() / 3 );
    _verts.resize( xyzPositions.size() / 3 );
    _is_vert_selected.resize( xyzPositions.size() / 3 );

    // The Eigen matrix has the same format as the incoming vector so we can straight copy it.
    // HACK: This is pretty sketchy and relies on Eigen::Vector3d having the same layout as a double[3] and nothing extra or fancy alignment.
    std::copy( xyzPositions.begin(), xyzPositions.end(), _verts.front().data() );

    init_adjacency( xyzPositions.size() / 3, triangleVerts, _he_data, _face_to_he, _vert_to_he );
  }

  EditMesh* EditMesh::clone()
  {
    EditMesh* m = new EditMesh();
    m->_verts = _verts;
    m->_he_data = _he_data;
    m->_face_to_he = _face_to_he;
    m->_vert_to_he = _vert_to_he;

    return m;
  }

  HalfEdge* EditMesh::find_twin( std::size_t vFrom, std::size_t vTo ){
    vvert_iterator it;
    if( !this->init_iterator( it, vFrom ) )
      return NULL;
	
    do{
      if( this->deref_iterator( it ) == vTo )
	return const_cast<HalfEdge*>( it.m_cur ); // Gross. This is just laziness.
    }while( this->advance_iterator( it ) );

    return NULL;
  }

  uint EditMesh::count_edge( std::size_t vFrom, std::size_t vTo )
  {
    uint ans = 0;
    vvert_iterator it;
    if( !this->init_iterator( it, vFrom ) )
      return ans;
	
    do
      {
        if( this->deref_iterator( it ) == vTo )  ans++;
      }while( this->advance_iterator( it ) );

    return ans;
  }

  HalfEdge* EditMesh::find_edge( std::size_t vFrom, std::size_t vTo ){
    if( HalfEdge* he = this->find_twin( vFrom, vTo ) )
      return &_he_data[ he->twin ];
    return NULL;
  }

  bool EditMesh::flip_edge( HalfEdge &he ){
    HalfEdge &twin = _he_data[ he.twin ];

    if (he.face == HOLE_INDEX ||
        twin.face == HOLE_INDEX)
      return false;

    std::size_t he_tri[3];
    std::size_t twin_tri[3];

    // prep: gather half edge indices in
    // the order they should be after flip
    he_tri[0]   = he.next;
    twin_tri[0] = twin.next;
    he_tri[1]   = twin.twin;
    twin_tri[1] = he.twin;
    he_tri[2]   = _he_data[ twin_tri[0] ].next;
    twin_tri[2] = _he_data[ he_tri[0] ].next;

    if( _he_data[ he_tri[2] ].vert == _he_data[ twin_tri[2] ].vert )
      return false;

    // step 1: ensure he's verts don't point to
    // either HalfEdge (does not break mesh)
    _vert_to_he[ he.vert ] = twin_tri[0];
    _vert_to_he[ twin.vert ] = he_tri[0];

    // step 2: set the he's vert to new originating vert
    he.vert = _he_data[ twin_tri[2] ].vert;
    twin.vert = _he_data[ he_tri[2] ].vert;
    
    // step 3: ensure the faces point to one
    // of the half edges connected to them
    _face_to_he[ he.face ] = he_tri[0];
    _face_to_he[ twin.face ] = twin_tri[0];

    // step 4: fix two edges that will point
    // to the wrong face
    _he_data[he_tri[2]].face = he.face;
    _he_data[twin_tri[2]].face = twin.face;

    // step 5: ensure half edges point to
    // each other
    for( int i=0; i<3; ++i ) {
      _he_data[ he_tri[i] ].next = he_tri[(i+1)%3];
      _he_data[ twin_tri[i] ].next = twin_tri[(i+1)%3];
    }

    return true;
  }

  // IMPORTANT: Given a collection of half-edges to delete (ex. When
  // removing a face we need to kill 2, 4, or 6 half-edges) they must be
  // deleting in decreasing index order!
  void EditMesh::delete_HalfEdge_impl( std::size_t he )
  {
    assert( (_he_data[he].vert >= _vert_to_he.size() || _vert_to_he[_he_data[he].vert] != he) && "Deleting this HalfEdge leaves a dangling link from a vertex. Must handle this first" );
		
    // Move a HalfEdge from the end overtop of the HalfEdge we are
    // deleting, then update the indices of linked HalfEdges.
    _he_data[he] = _he_data.back();
    _he_data.pop_back();

    // We may have just deleted the item at the end of the list, so we
    // have nothing to update since the indices didn't change.
    if( he != _he_data.size() ){
      const HalfEdge& heMoved = _he_data[he];

      // If the moved HalfEdge was the arbitrary HalfEdge linked to
      // the vertex, update it.
      if( _vert_to_he[heMoved.vert] == _he_data.size() )
	_vert_to_he[heMoved.vert] = he;

      // If the moved HalfEdge was the arbitrary HalfEdge linked to
      // the face, update it.
      if( heMoved.face != HOLE_INDEX && _face_to_he[heMoved.face] == _he_data.size() )
	_face_to_he[heMoved.face] = he;

      assert( heMoved.twin < _he_data.size() );
      assert( _he_data[heMoved.twin].twin == _he_data.size() );
      _he_data[heMoved.twin].twin = he;

      // NOTE: If we are deleting a bundle of HalfEdges, then by
      //       definition we must call delete_HalfEdge() in
      //       decreasing order of indices. That prevents me from
      //       having to worry about moving a partially destroyed
      //       HalfEdge into the 'he' position.

#ifdef USE_PREV
      assert( _he_data[heMoved.prev].next == _he_data.size() );
      _he_data[heMoved.prev].next = he;

      assert( _he_data[heMoved.next].prev == _he_data.size() );
      _he_data[heMoved.next].prev = he;
#else
      // Have to loop around the face until we find the HalfEdge
      // using 'heMoved' as its 'next' entry, then update it.
      std::size_t hePrev = heMoved.next;
      while( _he_data[hePrev].next != _he_data.size() )
	hePrev = _he_data[hePrev].next;

      assert( _he_data[hePrev].next == _he_data.size() );
      _he_data[hePrev].next = he;
#endif
    }
  }

  template <std::size_t N>
  void EditMesh::delete_HalfEdges_impl( std::size_t (&heToDelete)[N] ){
    std::sort( heToDelete, heToDelete + N, std::greater<std::size_t>() );
    for( std::size_t i = 0; i < N; ++i )
      this->delete_HalfEdge_impl( heToDelete[i] );
  }

  bool g_debug = false;

  std::size_t EditMesh::collapse_edge( std::size_t he ){
    assert( he < _he_data.size() );
    assert( _he_data[he].face != HOLE_INDEX && _he_data[_he_data[he].twin].face != HOLE_INDEX && "Cannot collapse a boundary edge" );

    const HalfEdge& heBase = _he_data[he];
    const HalfEdge& heTwin = _he_data[heBase.twin];

    // We are going to delete the faces on either side of the chosen
    // edge, so we need to delete 3 HalfEdges and patch up the twin
    // links on the 4 bordering edges.
    std::size_t heBorder[4];
    heBorder[0] = _he_data[ heBase.next ].twin;
    heBorder[1] = _he_data[ _he_data[ heBase.next ].next ].twin;
    heBorder[2] = _he_data[ _he_data[ heTwin.next ].next ].twin;
    heBorder[3] = _he_data[ heTwin.next ].twin;

    // TODO: Relax this assertion. We should be able to collapse a
    // spike jutting into a hole.
    assert( ( _he_data[ heBorder[0] ].face != HOLE_INDEX || _he_data[ heBorder[1] ].face != HOLE_INDEX ) && "Cannot collapse an edge on a face with holes on either side." );
    assert( ( _he_data[ heBorder[2] ].face != HOLE_INDEX || _he_data[ heBorder[3] ].face != HOLE_INDEX ) && "Cannot collapse an edge on a face with holes on either side." );

    // Check if we can actually collapse. This checks for a degree 3
    // vertex at the vertices not on the edge we are collapsing.
    if( _he_data[ _he_data[ _he_data[ heBorder[1] ].next ].twin ].next == heBorder[0] )
      return HOLE_INDEX;
    if( _he_data[ _he_data[ _he_data[ heBorder[2] ].next ].twin ].next == heBorder[3] )
      return HOLE_INDEX;

    // Capture the indices of things (2 faces & 6 half-edges) we want
    // to delete.
    std::size_t fToDelete[] = { heBase.face, heTwin.face };
    std::size_t heToDelete[] = { he, heBase.next, _he_data[ heBase.next ].next, heBase.twin, heTwin.next, _he_data[ heTwin.next ].next };
	
#ifndef NDEBUG
    // We can't be deleting border edges!
    for( auto i : heToDelete ){
      if( std::find( heBorder, heBorder + 4, i ) != heBorder + 4 )
	return HOLE_INDEX;	
      //assert( std::find( heBorder, heBorder + 4, i ) == heBorder + 4 );
    }

    if( g_debug ){
      std::vector< std::set<std::size_t> > verts( 3 );

      verts[0].insert( heBase.vert );
      verts[0].insert( heTwin.vert );

      for( size_t i = 1; i < verts.size(); ++i ){
	for( auto v : verts[i-1] ){
	  vvert_iterator it;
	  this->init_iterator( it, v );
	  do{
	    verts[i].insert( this->deref_iterator( it ) );
	  }while( this->advance_iterator( it ) );
	}
      }

      std::vector<std::size_t> orderedVerts( verts.back().begin(), verts.back().end() );
      std::set<std::size_t> faces;

      std::vector< double > vpos;
      std::vector< std::size_t > finds;

      for( auto v : orderedVerts ){
	vpos.push_back( _verts[v].x() ); vpos.push_back( _verts[v].y() ); vpos.push_back( _verts[v].z() );
	//std::clog << "m.add_vert( " << _verts[v].x() << ", " << _verts[v].y() << ", " << _verts[v].z() << " );" << std::endl;
      }

      // Visit the 1-ring
      for( auto v : verts[1] ){
	vface_iterator it;
	this->init_iterator( it, v );
	do{
	  if( this->deref_iterator( it ) != HOLE_INDEX && faces.find( this->deref_iterator( it ) ) == faces.end() ){
	    faces.insert( this->deref_iterator( it ) );

	    fvert_iterator itFace;
	    this->init_iterator( itFace, this->deref_iterator( it ) );

	    std::size_t f[3];
	    std::size_t i = 0;
	    do{
	      f[i++] = std::find( orderedVerts.begin(), orderedVerts.end(), this->deref_iterator( itFace ) ) - orderedVerts.begin();
	    }while( this->advance_iterator( itFace ) );

	    finds.push_back( f[0] ); finds.push_back( f[1] ); finds.push_back( f[2] );
	    //std::clog << "m.add_face( " << f[0] << ", " << f[1] << ", " << f[2] << " );" << std::endl;
	  }	
	}while( this->advance_iterator( it ) );
      }

      std::size_t base = std::find( orderedVerts.begin(), orderedVerts.end(), heBase.vert ) - orderedVerts.begin();
      std::size_t twin = std::find( orderedVerts.begin(), orderedVerts.end(), heTwin.vert ) - orderedVerts.begin();
      std::clog << "m.collapse_edge( " << base << ", " << twin << " );" << std::endl;

      EditMesh m;
      m.init( vpos, finds );
      std::ofstream fout( "debug.obj" );
      m.write_to_obj_stream( fout );
      fout.close();
    }
#endif

    // We may also need to fix the vertex->HalfEdge link for the
    // verts using these faces. There are technically 4, but we only
    // update the 3 that are not going to be deleted.
    std::size_t verts[] = { this->prev( heBase ).vert, heBase.vert, this->prev( heTwin ).vert };

    // Move the base vertex (arbitrarily) to the middle of the
    // edge. Could leave it where it is, or do something fancier too.
    _verts[heBase.vert] = 0.5 * ( _verts[heBase.vert] + _verts[heTwin.vert] ); 

    // Adjust all the twin's 1-ring to link to the vertex we are not
    // going to delete.
    std::size_t heIt = this->twin(this->next(heBase)).next;
    std::size_t heEnd = heBase.twin;
    for( ; heIt != heEnd; heIt = this->twin( _he_data[heIt] ).next ){
      assert( _he_data[heIt].vert == heTwin.vert );
		
      // Associate to the other vertex now, so we can delete this one.
      _he_data[heIt].vert = heBase.vert;
    }

    // Fix the vert associations if required, picking a non-hole face.
    if( _vert_to_he[ verts[0] ] == _he_data[ heBorder[1] ].twin )
      _vert_to_he[ verts[0] ] = (_he_data[ heBorder[0] ].face != HOLE_INDEX) ? heBorder[0] : _he_data[ heBorder[1] ].next;
    if( _vert_to_he[ verts[1] ] == he || _vert_to_he[ verts[1] ] == heTwin.next )
      _vert_to_he[ verts[1] ] = (_he_data[ heBorder[1] ].face != HOLE_INDEX) ? heBorder[1] : heBorder[2];
    if( _vert_to_he[ verts[2] ] == _he_data[ heBorder[2] ].twin )
      _vert_to_he[ verts[2] ] = (_he_data[ heBorder[3] ].face != HOLE_INDEX) ? heBorder[3] : _he_data[ heBorder[2] ].next;

    // "Delete" the other vertex
    _vert_to_he[heTwin.vert] = HOLE_INDEX;

    // Collapse the two triangles bordering our chosen half-edge by
    // connecting the opposite edges together.
    _he_data[ heBorder[0] ].twin = heBorder[1];
    _he_data[ heBorder[1] ].twin = heBorder[0];
    _he_data[ heBorder[2] ].twin = heBorder[3];
    _he_data[ heBorder[3] ].twin = heBorder[2];

    // Have to delete the faces in the proper order.
    if( fToDelete[0] < fToDelete[1] )
      std::swap( fToDelete[0], fToDelete[1] );

    this->delete_HalfEdges_impl( heToDelete );
    detail::delete_face( _face_to_he, _he_data, fToDelete[0] );
    detail::delete_face( _face_to_he, _he_data, fToDelete[1] );

    return verts[1];
  }

  std::size_t EditMesh::add_face( std::size_t (&v)[3] ){
    std::size_t faceIndex = _face_to_he.size();
    std::size_t heIndex = _he_data.size();

    // Find the half-edges on the hole face we are filling. We must either:
    //  1. Find no half-edges, if all vertices are unconnected from the mesh.
    //  3. Find one half-edge, if one of the vertices is not connected to the existing mesh.
    //  4. Find two half-edges, if we are adding a triangle inside of a polygonal hole.
    //  2. Find three half-edges, if we are filling an existing triangular hole.
    HalfEdge* he[] = { 
      this->find_edge( v[0], v[1] ), 
      this->find_edge( v[1], v[2] ), 
      this->find_edge( v[2], v[0] ) };
	
    // Find the first half-edge we need to modify. This is an edge
    std::size_t base = HOLE_INDEX;
    for( std::size_t i = 0; i < 3 && base == HOLE_INDEX; ++i ){
      if( he[i] ){
	assert( he[i]->face == HOLE_INDEX && "Non-manifold mesh detected. Cannot connect to an edge which already has two incident faces (ie. One side must be a hole)" );
	if( !he[(i+2)%3] )
	  base = i;
      }
    }

    if( base == HOLE_INDEX ){
      // This triangle is not connected to any others, or we
      // completely filled a triangular hole.
      if( he[0] /*|| he[1] || he[2]*/ ){
	assert( he[0] && he[1] && he[2] );
	assert( he[0]->face == HOLE_INDEX && he[1]->face == HOLE_INDEX && he[2]->face == HOLE_INDEX );
	assert( &_he_data[ he[0]->next ] == he[1] && &_he_data[ he[1]->next ] == he[2] && &_he_data[ he[2]->next ] == he[0] );
			
	// Update the face index of the triangular hole to convert
	// it to a face.
	he[0]->face = he[1]->face = he[2]->face = faceIndex;
	_face_to_he.push_back( he[2]->next );
      }else{
	assert( !he[0] && !he[1] && !he[2] );
	assert( _vert_to_he[v[0]] == HOLE_INDEX && _vert_to_he[v[1]] == HOLE_INDEX && _vert_to_he[v[2]] == HOLE_INDEX && "Non-manifold mesh detected. Cannot have two hole faces incident on a vertex." );

	// Make 3 new half-edges for the triangle, and 3 new
	// half-edges for the hole outside of the triangle.
	HalfEdge newHe[6];
	detail::init( newHe[0], heIndex+1, heIndex+5, v[0], faceIndex );
	detail::init( newHe[1], heIndex+2, heIndex+4, v[1], faceIndex );
	detail::init( newHe[2], heIndex  , heIndex+3, v[2], faceIndex );
	detail::init( newHe[3], heIndex+4, heIndex+2, v[0], HOLE_INDEX );
	detail::init( newHe[4], heIndex+5, heIndex+1, v[2], HOLE_INDEX );
	detail::init( newHe[5], heIndex+3, heIndex  , v[1], HOLE_INDEX );
#ifdef USE_PREV
	newHe[0].prev = heIndex+2;
	newHe[1].prev = heIndex;
	newHe[2].prev = heIndex+1;

	newHe[3].prev = heIndex+5;
	newHe[4].prev = heIndex+3;
	newHe[5].prev = heIndex+4;
#endif

	_vert_to_he[ v[0] ] = heIndex;
	_vert_to_he[ v[1] ] = heIndex+1;
	_vert_to_he[ v[2] ] = heIndex+2;

	_face_to_he.push_back( heIndex );
	_he_data.push_back( newHe[0] );
	_he_data.push_back( newHe[1] );
	_he_data.push_back( newHe[2] );
	_he_data.push_back( newHe[3] );
	_he_data.push_back( newHe[4] );
	_he_data.push_back( newHe[5] );
      }
    }else{
      std::size_t next = (base+1)%3, prev = (base+2)%3;
      std::size_t baseIndex = static_cast<std::size_t>( he[base] - &_he_data.front() );

      assert( !he[prev] );

      if( he[next] ){
	// We have two edges to steal from the hole, and we need to add two new half-edges
	HalfEdge newHe[2];
	detail::init( newHe[0], baseIndex, heIndex+1, v[prev], faceIndex );
	detail::init( newHe[1], he[next]->next, heIndex, v[base], HOLE_INDEX );

#ifdef USE_PREV
	newHe[0].prev = he[base]->next;
	newHe[1].prev = he[base]->prev;

	_he_data[ he[base]->prev ].next = heIndex + 1;
	_he_data[ he[next]->next ].prev = heIndex + 1;

	he[next]->next = heIndex;
	he[base]->prev = heIndex;
#else
	// Have to find the previous HalfEdge in the polygonal
	// hole so we can point it to the new half-edge in the
	// hole.
	HalfEdge* hePrev = &_he_data[ he[next]->next ];
	while( &_he_data[hePrev->next] != he[base] ){
	  hePrev = &_he_data[hePrev->next];
	  assert( hePrev != he[next] ); // To catch weirdness.
	}
	assert( &_he_data[hePrev->next] == he[base] );
			
	hePrev->next = heIndex + 1;
	he[next]->next = heIndex;
#endif

	// Update the face indices of the half-edges to indicate
	// they are in a triangle now.
	he[base]->face = he[next]->face = faceIndex;

	_face_to_he.push_back( heIndex );
	_he_data.push_back( newHe[0] );
	_he_data.push_back( newHe[1] );
      }else{
	assert( _vert_to_he[ v[prev] ] == HOLE_INDEX && "Non-manifold mesh detected. Cannot have two hole faces incident on a vertex." );

	// We have one edge to steal from the hole, and we need to
	// add four new half-edges.
	HalfEdge newHe[4];
	detail::init( newHe[0], baseIndex, heIndex+2, v[prev], faceIndex );
	detail::init( newHe[1], heIndex  , heIndex+3, v[next], faceIndex );
	detail::init( newHe[2], heIndex+3, heIndex  , v[base], HOLE_INDEX );
	detail::init( newHe[3], he[base]->next, heIndex+1, v[prev], HOLE_INDEX );

#ifdef USE_PREV
	newHe[0].prev = heIndex+1;
	newHe[1].prev = baseIndex;
	newHe[2].prev = he[base]->prev;
	newHe[3].prev = heIndex+2;

	_he_data[ he[base]->prev ].next = heIndex+2;
	_he_data[ he[base]->next ].prev = heIndex+3;

	he[base]->prev = heIndex;
	he[base]->next = heIndex+1;
#else
	// Have to find the previous HalfEdge in the polyognal
	// hole so we can point it to the new half-edge in the
	// hole.
	HalfEdge* hePrev = &_he_data[ he[base]->next ];
	while( &_he_data[hePrev->next] != he[base] ){
	  hePrev = &_he_data[hePrev->next];
	  assert( hePrev != he[next] ); // To catch weirdness.
	}
	assert( &_he_data[hePrev->next] == he[base] );
			
	hePrev->next = heIndex+2;
	he[base]->next = heIndex+1;
#endif

	// Update the face indices of the half-edges to indicate
	// they are in a triangle now.
	he[base]->face = faceIndex;

	_vert_to_he[v[prev]] = heIndex;
	_face_to_he.push_back( heIndex );
	_he_data.push_back( newHe[0] );
	_he_data.push_back( newHe[1] );
	_he_data.push_back( newHe[2] );
	_he_data.push_back( newHe[3] );
      }
    }

    return faceIndex;
  }

  void EditMesh::delete_face( std::size_t f ){
    assert( f < _face_to_he.size() );

    // We can assume that this face has 3 half-edges.
    std::size_t heIndices[3];
    heIndices[0] = _face_to_he[f];
	
    HalfEdge* he[3];
    he[0] = &_he_data[heIndices[0]];
    he[1] = &_he_data[he[0]->next];
    he[2] = &_he_data[he[1]->next];

    heIndices[1] = he[0]->next;
    heIndices[2] = he[1]->next;
	
    assert( he[0]->face == f && he[1]->face == f && he[2]->face == f );
    assert( he[2]->next == _face_to_he[f] );

    // Search for an edge that has a neighbor, but its prev edge
    // doesn't. This is a canonical place to construct the algorithm
    // from.
    std::size_t base = HOLE_INDEX;
    for( std::size_t i = 0; i < 3 && base == HOLE_INDEX; ++i ){
      if( _he_data[he[i]->twin].face != HOLE_INDEX && _he_data[he[(i+2)%3]->twin].face == HOLE_INDEX )
	base = i;
    }

    if( base == HOLE_INDEX ){
      if( _he_data[he[0]->twin].face == HOLE_INDEX ){
	// This is a lone triangle, so delete its half-edges and
	// the exterior hole surrounding it too.

	// TODO: Remove the floating vertices? Currently we are leaving them.
	_vert_to_he[he[0]->vert] = HOLE_INDEX;
	_vert_to_he[he[1]->vert] = HOLE_INDEX;
	_vert_to_he[he[2]->vert] = HOLE_INDEX;

	// Delete all of the edges (both inside & outside
	// half-edges). Must do this last since indices can change
	// arbitrarily when deleting.
	std::size_t toDelete[] = { 
	  heIndices[0], heIndices[1], heIndices[2], 
	  he[0]->twin, he[1]->twin, he[2]->twin 
	};
			
	this->delete_HalfEdges_impl( toDelete );
	detail::delete_face( _face_to_he, _he_data, f );
      }else{
	// This is an interior triangle. Only have to change the
	// face_index to HOLE_INDEX for these edges.

	// Adjust any vertex references to new edges in non-hole faces.
	if( _vert_to_he[he[0]->vert] == heIndices[0] )
	  _vert_to_he[he[0]->vert] = he[2]->twin;
	if( _vert_to_he[he[1]->vert] == heIndices[1] )
	  _vert_to_he[he[1]->vert] = he[0]->twin;
	if( _vert_to_he[he[2]->vert] == heIndices[2] )
	  _vert_to_he[he[2]->vert] = he[1]->twin;

	// Flag all these half-edges as being a hole now.
	he[0]->face = he[1]->face = he[2]->face = HOLE_INDEX;
	detail::delete_face( _face_to_he, _he_data, f );
      }
    }else{
      std::rotate( he, he+base, he+3 );
      std::rotate( heIndices, heIndices+base, heIndices+3 );
      assert( _he_data[he[0]->twin].face != HOLE_INDEX );
      assert( _he_data[he[2]->twin].face == HOLE_INDEX );

      if( _he_data[he[1]->twin].face != HOLE_INDEX ){
	// We have one edge to remove, and a hole to connect to.
#ifdef USE_PREV
	he[1]->next = _he_data[he[2]->twin].next;
	he[0]->prev = _he_data[he[2]->twin].prev;
	_he_data[he[1]->next].prev = heIndices[1];
	_he_data[he[0]->prev].next = heIndices[0];
#else
	he[1]->next = _he_data[he[2]->twin].next;

	std::size_t hePrev = he[1]->next;
	while( _he_data[hePrev].next != he[2]->twin )
	  hePrev = _he_data[hePrev].next;

	assert( _he_data[hePrev].next == he[2]->twin );
	_he_data[hePrev].next = heIndices[0];
#endif

	assert( _he_data[ _vert_to_he[ he[0]->vert ] ].face != HOLE_INDEX );
	assert( _he_data[ _vert_to_he[ he[1]->vert ] ].face != HOLE_INDEX );
	assert( _he_data[ _vert_to_he[ he[2]->vert ] ].face != HOLE_INDEX );

	// We may need to update the vertices if they referenced
	// the edges we are deleting. Choose new half-edges that
	// are inside non-hole triangles.
	if( _vert_to_he[he[0]->vert] == heIndices[0] )
	  _vert_to_he[he[0]->vert] = _he_data[he[0]->twin].next;
	if( _vert_to_he[he[1]->vert] == heIndices[1] )
	  _vert_to_he[he[1]->vert] = he[0]->twin;
	if( _vert_to_he[he[2]->vert] == heIndices[2] )
	  _vert_to_he[he[2]->vert] = he[1]->twin;

	assert( _he_data[ _vert_to_he[ he[0]->vert ] ].face != HOLE_INDEX );
	assert( _he_data[ _vert_to_he[ he[1]->vert ] ].face != HOLE_INDEX );
	assert( _he_data[ _vert_to_he[ he[2]->vert ] ].face != HOLE_INDEX );

	he[0]->face = he[1]->face = HOLE_INDEX;

	std::size_t toDelete[] = { heIndices[2], he[2]->twin };

	// Delete the edges and face. Must do this last since
	// indices can change arbitrarily when deleting.
	this->delete_HalfEdges_impl( toDelete );
	detail::delete_face( _face_to_he, _he_data, f );
      }else{
	// We have two edges to remove, a vertex that will become
	// floating, and a hole to connect to.
#ifdef USE_PREV
	he[0]->next = _he_data[he[1]->twin].next;
	he[0]->prev = _he_data[he[2]->twin].prev;
	_he_data[he[0]->next].prev = heIndices[0];
	_he_data[he[0]->prev].next = heIndices[0];
#else
	he[0]->next = _he_data[he[1]->twin].next;

	std::size_t hePrev = he[0]->next;
	while( _he_data[hePrev].next != he[2]->twin )
	  hePrev = _he_data[hePrev].next;

	assert( _he_data[hePrev].next == he[2]->twin );
	_he_data[hePrev].next = heIndices[0];
#endif

	// We may need to update the vertices if they referenced
	// the edges we are deleting. Choose new half-edges that
	// are inside non-hole triangles.
	if( _vert_to_he[he[1]->vert] == heIndices[1] )
	  _vert_to_he[he[1]->vert] = he[0]->twin;
	if( _vert_to_he[he[0]->vert] == he[2]->twin || _vert_to_he[he[0]->vert] == heIndices[0] )
	  _vert_to_he[he[0]->vert] = _he_data[he[0]->twin].next;
	_vert_to_he[he[2]->vert] = HOLE_INDEX;

	// Update the face association of the one half-edge we are
	// keeping (it joins the hole).
	he[0]->face = HOLE_INDEX;
			
	// Delete the edges and face. Must do this last since
	// indices can change arbitrarily when deleting.
	std::size_t toDelete[] = { 
	  heIndices[1], heIndices[2], 
	  he[1]->twin, he[2]->twin 
	};

	this->delete_HalfEdges_impl( toDelete );
	detail::delete_face( _face_to_he, _he_data, f );
      }
    }
  }

  std::size_t EditMesh::split_face_center( std::size_t f, std::size_t (*pOutFaceIndices)[3] ){
    assert( f < _face_to_he.size() );
    assert( _face_to_he[f] < _he_data.size() && _he_data[_face_to_he[f]].face == f );

    std::size_t he[3];
    he[0] = _face_to_he[f];
    he[1] = _he_data[he[0]].next;
    he[2] = _he_data[he[1]].next;

    assert( _he_data[he[2]].next == he[0] );
    assert( _he_data[he[0]].vert < _verts.size() && _he_data[he[1]].vert < _verts.size() && _he_data[he[2]].vert < _verts.size() );
	
    // New vert at face center
    Eigen::Vector3d newVert = ( _verts[ _he_data[ he[0] ].vert ] + _verts[ _he_data[ he[1] ].vert ] + _verts[ _he_data[ he[2] ].vert ] ) / 3.0;

    std::size_t newVertIndex = _verts.size();
    _verts.push_back( newVert );

    // Each half-edge gets associated to a new face, and we add 6
    // half-edges from the old vertices to the new.
    std::size_t newHeIndex = _he_data.size();
    std::size_t newFaceIndex = _face_to_he.size();

    if( pOutFaceIndices ){
      (*pOutFaceIndices)[0] = f;
      (*pOutFaceIndices)[1] = newFaceIndex;
      (*pOutFaceIndices)[2] = newFaceIndex+1;
    }

    // Create six new half-edges connecting the center vertex to the
    // old triangle corners.
    HalfEdge newHe[6];
    detail::init( newHe[0], newHeIndex+1, newHeIndex+3, _he_data[he[1]].vert, f );
    detail::init( newHe[1], he[0]       , newHeIndex+4, newVertIndex        , f );
    detail::init( newHe[2], newHeIndex+3, newHeIndex+5, _he_data[he[2]].vert, newFaceIndex );
    detail::init( newHe[3], he[1]       , newHeIndex  , newVertIndex        , newFaceIndex );
    detail::init( newHe[4], newHeIndex+5, newHeIndex+1, _he_data[he[0]].vert, newFaceIndex+1 );
    detail::init( newHe[5], he[2]       , newHeIndex+2, newVertIndex        , newFaceIndex+1 );

    // Connect the old half-edges to the new ones, and update their
    //face association.  _he_data[he[0]].face = f;
    _he_data[he[0]].next = newHeIndex;
    _he_data[he[1]].face = newFaceIndex;
    _he_data[he[1]].next = newHeIndex+2;
    _he_data[he[2]].face = newFaceIndex+1;
    _he_data[he[2]].next = newHeIndex+4;

#ifdef USE_PREV
    newHe[0].prev = he[0];
    newHe[1].prev = newHeIndex;
    newHe[2].prev = he[1];
    newHe[3].prev = newHeIndex+2;
    newHe[4].prev = he[2];
    newHe[5].prev = newHeIndex+4;

    _he_data[he[0]].prev = newHeIndex+1;
    _he_data[he[1]].prev = newHeIndex+3;
    _he_data[he[2]].prev = newHeIndex+5;
#endif

    _vert_to_he.push_back( newHeIndex+3 ); // Arbitrary from 1, 3 & 5
    _face_to_he[f] = he[0];
    _face_to_he.push_back( he[1] );
    _face_to_he.push_back( he[2] );
    _he_data.push_back( newHe[0] );
    _he_data.push_back( newHe[1] );
    _he_data.push_back( newHe[2] );
    _he_data.push_back( newHe[3] );
    _he_data.push_back( newHe[4] );
    _he_data.push_back( newHe[5] );

    return newVertIndex;
  }

  void EditMesh::split_boundary_edge( std::size_t heToSplit, std::size_t (*pOutVertIndices)[2], std::size_t (*pOutFaceIndices)[3] ){
    assert( heToSplit < _he_data.size() );
    assert( _he_data[heToSplit].face != HOLE_INDEX && _he_data[_he_data[heToSplit].twin].face == HOLE_INDEX );

    HalfEdge& heBase = _he_data[heToSplit];
    HalfEdge& heNext = _he_data[heBase.next];
    HalfEdge& hePrev = _he_data[heNext.next];
    HalfEdge& heTwin = _he_data[heBase.twin];

    Eigen::Vector3d newVert1 = _verts[ heBase.vert ] + ( _verts[ heNext.vert ] - _verts[ heBase.vert ] ) / 3.0;
    Eigen::Vector3d newVert2 = _verts[ heBase.vert ] + ( _verts[ heNext.vert ] - _verts[ heBase.vert ] ) * (2.0 / 3.0);

    // Construct 2 new faces and 8 new HalfEdges connecting the new
    // verts to the off-edge vert.
    std::size_t newVertIndex = _verts.size();
    std::size_t newHeIndex = _he_data.size();
    std::size_t newFaceIndex = _face_to_he.size();

    if( pOutVertIndices ){
      (*pOutVertIndices)[0] = newVertIndex;
      (*pOutVertIndices)[1] = newVertIndex+1;
    }

    if( pOutFaceIndices ){
      (*pOutFaceIndices)[0] = heBase.face;
      (*pOutFaceIndices)[1] = newFaceIndex;
      (*pOutFaceIndices)[2] = newFaceIndex+1;
    }

    HalfEdge newHe[8];
    detail::init( newHe[0], heNext.next , newHeIndex+1, newVertIndex   , heBase.face );
    detail::init( newHe[1], newHeIndex+4, newHeIndex  , hePrev.vert    , newFaceIndex );
    detail::init( newHe[2], newHeIndex+1, newHeIndex+3, newVertIndex+1 , newFaceIndex );
    detail::init( newHe[3], newHeIndex+5, newHeIndex+2, hePrev.vert    , newFaceIndex+1 );
    detail::init( newHe[4], newHeIndex+2, newHeIndex+6, newVertIndex   , newFaceIndex );
    detail::init( newHe[5], heBase.next , heBase.twin , newVertIndex+1 , newFaceIndex+1 );
    detail::init( newHe[6], newHeIndex+7, newHeIndex+4, newVertIndex+1 , HOLE_INDEX );
    detail::init( newHe[7], heTwin.next , heToSplit   , newVertIndex   , HOLE_INDEX );

#ifdef USE_PREV
    newHe[0].prev = heToSplit;
    newHe[1].prev = newHeIndex+2;
    newHe[2].prev = newHeIndex+4;
    newHe[3].prev = heBase.next;
    newHe[4].prev = newHeIndex+1;
    newHe[5].prev = newHeIndex+3;
    newHe[6].prev = heBase.twin;
    newHe[7].prev = newHeIndex+6;

    heNext.prev = newHeIndex+5;
    _he_data[heTwin.next].prev = newHeIndex+7
#endif

      heBase.next = newHeIndex;
    heBase.twin = newHeIndex+7;
    heTwin.next = newHeIndex+6;
    heTwin.twin = newHeIndex+5;
    heNext.next = newHeIndex+3;
    heNext.face = newFaceIndex+1;

    _verts.push_back( newVert1 );
    _verts.push_back( newVert2 );
    _vert_to_he.push_back( newHeIndex+4 );
    _vert_to_he.push_back( newHeIndex+5 );
    _face_to_he[heBase.face] = heToSplit;
    _face_to_he.push_back( newHeIndex+4 );
    _face_to_he.push_back( newHeIndex+5 );
    _he_data.push_back( newHe[0] );
    _he_data.push_back( newHe[1] );
    _he_data.push_back( newHe[2] );
    _he_data.push_back( newHe[3] );
    _he_data.push_back( newHe[4] );
    _he_data.push_back( newHe[5] );
    _he_data.push_back( newHe[6] );
    _he_data.push_back( newHe[7] );
  }

  double EditMesh::get_cotan_weight( const vvert_iterator& it ) const {
    const HalfEdge *itCur = it.m_cur;
    const HalfEdge *itTwin = &_he_data[ itCur->twin ];

    double result = 0;

    Eigen::Vector3d a = this->get_vertex( itTwin->vert );
    Eigen::Vector3d b = this->get_vertex( itCur->vert );

    assert( ( itCur->face != HOLE_INDEX || itTwin->face != HOLE_INDEX ) && "Invalid mesh: edge with no face on either side" );

    if( itCur->face != HOLE_INDEX ){
      Eigen::Vector3d c = this->get_vertex( this->prev( *itCur ).vert );
      Eigen::Vector3d e0 = b - c;
      Eigen::Vector3d e1 = a - c;

      // We use the dot product and norm of the cross product to get cos and sin respectively. cotan = cos / sin
      result += static_cast<double>( e0.dot(e1) ) / e0.cross(e1).norm();
    }

    if( itTwin->face != HOLE_INDEX ){
      Eigen::Vector3d c = this->get_vertex( this->prev( *itTwin ).vert );
      Eigen::Vector3d e0 = a - c;
      Eigen::Vector3d e1 = b - c;

      result += static_cast<double>( e0.dot(e1) ) / e0.cross(e1).norm();
    }
	
    return result;
  }

  double EditMesh::get_mean_value_weight( const vvert_iterator& it ) const {
    const HalfEdge *itCur = it.m_cur;
    const HalfEdge *itTwin = &_he_data[ itCur->twin ];

    double result = 0;

    Eigen::Vector3d a = this->get_vertex( itTwin->vert );
    Eigen::Vector3d b = this->get_vertex( itCur->vert );
	
    Eigen::Vector3d e0 = (b - a);
    double eLen = e0.norm();

    e0 /= eLen;

    assert( ( itCur->face != HOLE_INDEX || itTwin->face != HOLE_INDEX ) && "Invalid mesh: edge with no face on either side" );

    if( itCur->face != HOLE_INDEX ){
      Eigen::Vector3d c = this->get_vertex( this->prev( *itCur ).vert );
      Eigen::Vector3d e1 = (c - a).normalized();

      result += std::tan( 0.5 * std::acos( e0.dot(e1) ) );
    }

    if( itTwin->face != HOLE_INDEX ){
      Eigen::Vector3d c = this->get_vertex( this->prev( *itTwin ).vert );
      Eigen::Vector3d e1 = (c - a).normalized();

      result += std::tan( 0.5 * std::acos( e0.dot(e1) ) );
    }
	
    return result / eLen;
  }

  Eigen::Vector3d EditMesh::get_normal( const vface_iterator& it ) const {
    Eigen::Vector3d a = this->get_vertex( it.m_next->vert );
    Eigen::Vector3d b = this->get_vertex( it.m_cur->vert );
    Eigen::Vector3d c = this->get_vertex( _he_data[ it.m_cur->next ].vert );
	
    return ( a - c ).cross( b - c ).normalized();
  }

  Eigen::Vector4d EditMesh::get_plane( const vface_iterator& it ) const {
    Eigen::Vector3d a = this->get_vertex( it.m_next->vert );
    Eigen::Vector3d b = this->get_vertex( it.m_cur->vert );
    Eigen::Vector3d c = this->get_vertex( _he_data[ it.m_cur->next ].vert );
    Eigen::Vector3d n = ( a - c ).cross( b - c ).normalized();

    // Plane equation Ax + By + Cz + D = 0 -> n.dot( [x,y,z] ) - n.dot( c ) = 0
    return Eigen::Vector4d( n.x(), n.y(), n.z(), -n.dot( c ) );
  }

  Eigen::Vector3d EditMesh::get_vnormal( std::size_t vertex ) const {
    vvert_iterator vit;
    init_iterator(vit, vertex);

    Eigen::Vector3d normal;
    normal.setZero();

    Eigen::Vector3d center = this->get_vertex( vertex );
    Eigen::Vector3d vec_prev;
    Eigen::Vector3d vec_curr = this->get_vertex( deref_iterator(vit) ) - center;

    advance_iterator(vit);
    vit.m_end = vit.m_cur;
    do {
      vec_prev = vec_curr;
      vec_curr = this->get_vertex( deref_iterator(vit) ) - center;

      if (_he_data[vit.m_cur->twin].face != HOLE_INDEX)
	normal += vec_curr.cross(vec_prev);
    } while (advance_iterator(vit));

    return normal.normalized();
  }

  Eigen::Vector3d EditMesh::get_fnormal( std::size_t face ) const {
    const HalfEdge *e1 = &_he_data[ _face_to_he[ face ] ];
    const HalfEdge *e2 = &this->next(*e1);
    const HalfEdge *e3 = &this->next(*e2);

    Eigen::Vector3d a = this->get_vertex( e1->vert );
    Eigen::Vector3d b = this->get_vertex( e2->vert );
    Eigen::Vector3d c = this->get_vertex( e3->vert );

    Eigen::Vector3d normal = (b - a).cross(c - a);

    return normal.normalized();
  }

  void EditMesh::getIndicesForFace(size_t tri_index, size_t indicesForFace[3]) const {
    fvert_iterator fvit;
    init_iterator( fvit, tri_index );

    for( size_t i = 0; i < 3; i++ ) {
      indicesForFace[i] = deref_iterator(fvit);
      advance_iterator( fvit );
    }
  }

  Eigen::Vector3d EditMesh::getFaceMidpoint(size_t tri_index) {
    fvert_iterator fvit;
    init_iterator( fvit, tri_index );

    Eigen::Vector3d average = Eigen::Vector3d::Zero();
    int count = 0;
    for( size_t i = 0; i < 3; i++ ) {
      size_t index = deref_iterator(fvit);
      if (index != HOLE_INDEX) {
	count++;
	average += get_vertex(index);
      }
      advance_iterator( fvit );
    }
    return average / count;
  }

  void EditMesh::get_draw_data( float *verts, int *indices ) const {

    /* get each vertex only once. This is good for efficiency
     * but results in bad looking meshes due to each vertex 
     * having a fixed normal

     for( std::size_t i = 0, iEnd = _face_to_he.size(); i < iEnd; i++ ){
     const HalfEdge* he = &_he_data[ _face_to_he[i] ];

     for( int j = 0; j < 3; j++){
     indices[3*i +j] = he->vert;
     he = &this->next(*he);
     }
     }

     for( std::size_t i = 0, iEnd = _vert_to_he.size(); i < iEnd; i++ ){
     Eigen::Vector3d vert = this->get_vertex( i );
     for( int j = 0; j < 3; j++)
     verts[3*i+j] = (float) vert[j];
     }
    */

    // for each face
    int continuous_id = 0;
    for( std::size_t i = 0, iEnd = _face_to_he.size(); i < iEnd; i++ )
      {
	if( is_simplification_in_progress() && (!_is_face_active[i]) ) continue;
	const HalfEdge* he = &_he_data[ _face_to_he[i] ];

	// for each vertex of the face
	for( int j = 0; j < 3; j++){
	  Eigen::Vector3d vert = get_vertex(he->vert);
	  if(indices) indices[3*continuous_id+j] = 3*continuous_id+j;

	  // for each component of the vertex
	  for( int k = 0; k < 3; k++){
	    // this is where we convert from the double-precision of the data structure to float for the graphics card
	    verts[3*(3*continuous_id+j) + k] = static_cast<float>(vert[k]);
	  }
	  he = &this->next(*he);
	}
	continuous_id++;
      }
  }

  void EditMesh::get_draw_normals( float *normals ) const {

    /* this finds the averaged vertex normals which results in
     * poor looking meshes when they are not smooth

     for( std::size_t i = 0, iEnd = _vert_to_he.size(); i < iEnd; i++ ){
     Eigen::Vector3d normal = this->get_normal( i );
     for( int j = 0; j < 3; j++)
     normals[3*i+j] = (float) normal[j];
     }
    */

    int continuous_id = 0;
    for( std::size_t f = 0, iEnd = _face_to_he.size(); f < iEnd; ++f )
      {
	if( is_simplification_in_progress() && (!_is_face_active[f]) ) continue;

	Eigen::Vector3d fnormal = this->get_fnormal( f );
	//HalfEdge he = _he_data[ _face_to_he[f] ];

	// for each vertex of the face
	for( int j = 0; j < 3; j++){
	  //Eigen::Vector3d vnormal = this->get_vertex(he.vert);

	  // for each component of the vertex
	  for( int k = 0; k < 3; k++){
	    // this is where we convert from the double-precision of the data structure to float for the graphics card
	    normals[3*(3*continuous_id+j) + k] = static_cast<float>(fnormal[k]);
	  }
	  //he = next(he);
	}
	continuous_id++;
      }
  }

  void EditMesh::get_draw_selection( int *selection ) const
  {
    int continuous_id = 0;
    for( std::size_t i = 0, iEnd = _face_to_he.size(); i < iEnd; i++ )
      {
	if( is_simplification_in_progress() && (!_is_face_active[i]) ) continue;
	size_t verts[3];
	getIndicesForFace(i, verts);

	for( int j = 0; j < 3; j++)
	  {
	    selection[3* continuous_id+j] = _is_vert_selected[verts[j]];
	  }
	continuous_id++;
      }
  }

  // call instead of init to test edge flip
  // easiest way is hacking it into mesh constructor
  void EditMesh::test_flip() {
    std::vector<double> xyz;
    std::vector<std::size_t> faces;

    // four verts
    xyz.push_back(-1); xyz.push_back(0); xyz.push_back(0);
    xyz.push_back(0); xyz.push_back(1); xyz.push_back(1);
    xyz.push_back(0); xyz.push_back(1); xyz.push_back(-1);
    xyz.push_back(1); xyz.push_back(0); xyz.push_back(0);

    // two triangles
    faces.push_back(0); faces.push_back(1);
    faces.push_back(2); faces.push_back(2);
    faces.push_back(1); faces.push_back(3);

    this->init(xyz, faces);

    HalfEdge *he = &_he_data[0];
    for( size_t i = 0; i < _he_data.size(); ++i ) {
      he = &_he_data[i];
      if (he->face != HOLE_INDEX &&
	  _he_data[ he->twin ].face != HOLE_INDEX )
	break;
    }
    flip_edge(*he);
    this->_edit_count++;
  }

  void EditMesh::write_to_obj_stream( std::ostream& stream ) const {
    for( auto& v : _verts )
      stream << "v " << v.x() << ' ' << v.y() << ' ' << v.z() << std::endl;
    stream << std::endl;
    for( std::size_t i = 0, iEnd = _face_to_he.size(); i < iEnd; ++i ){
      fvert_iterator it;
      this->init_iterator( it, i );
      stream << "f ";
      bool isFirst = true;
      do {
	if( !isFirst )
	  stream << ' ';
	isFirst = false;
	stream << this->deref_iterator( it )+1;
      } while ( this->advance_iterator( it ) );
      stream << std::endl;
    }
  }

  void EditMesh::verify() const {
    for( std::size_t i = 0, iEnd = _face_to_he.size(); i < iEnd; ++i ){

      if( _is_simplification_in_progress && (!_is_face_active[i]) ) continue;
        
      std::size_t c = 0;
		
      const HalfEdge* it = &_he_data[ _face_to_he[i] ];
      assert( it->next != _face_to_he[i] );
      while( it->next != _face_to_he[i] ){
	assert( it->face == i );
	assert( it->next != HOLE_INDEX && it->twin != HOLE_INDEX && it->vert < _vert_to_he.size() );
	assert( ( _he_data[ it->twin ].face == HOLE_INDEX || _he_data[ _he_data[it->next].twin ].face != _he_data[ it->twin ].face ) && "Can't have two edges shared between the same faces!" );
	it = &_he_data[it->next];
	assert( ++c < 1000000 ); // This isn't strictly a problem, but probably no face has a million verts in it.
      }
    }

    for( std::size_t i = 0, iEnd = _vert_to_he.size(); i < iEnd; ++i ){

      if( _is_simplification_in_progress && (!_is_vert_active[i]) ) continue;
 
      assert( _vert_to_he[i] == HOLE_INDEX || _vert_to_he[i] < _he_data.size() );
      if( _vert_to_he[i] != HOLE_INDEX ){
	const HalfEdge* it = &_he_data[ _vert_to_he[i] ];
	assert( it->vert == i );
	assert( it->face != HOLE_INDEX && "By convention, vertices should not reference hole faces" );
      }
    }

    for( std::size_t i = 0, iEnd = _he_data.size(); i < iEnd; ++i ){

      if( _is_simplification_in_progress && (!_is_he_active[i]) ) continue;
 
      const HalfEdge* it = &_he_data[i];
      assert( it->vert < _vert_to_he.size() );
      assert( it->face == HOLE_INDEX || it->face < _face_to_he.size() );

      assert( it->next < _he_data.size() );
      assert( it->next != i );
      assert( _he_data[it->next].face == it->face );
      assert( _he_data[it->next].vert != it->vert );

      assert( it->twin < _he_data.size() );
      assert( _he_data[it->twin].twin == i );
      assert( _he_data[it->twin].face != it->face );
      assert( _he_data[it->twin].vert == _he_data[it->next].vert );

#ifdef USE_PREV
      assert( it->prev < _he_data.size() );
      assert( it->prev != i );
      assert( _he_data[it->next].prev == i );
      assert( _he_data[it->prev].next == i );
      assert( _he_data[it->prev].face == it->face );
      assert( _he_data[it->prev].vert != it->vert );
#endif
    }
  }

  void EditMesh::test(){
    EditMesh m1, m2, m3;

    std::vector<double> v;
    std::vector<std::size_t> f;

    v.push_back( 0 ); v.push_back( 0 ); v.push_back( 0 );
    v.push_back( 1 ); v.push_back( 0 ); v.push_back( 0 );
    v.push_back( 0 ); v.push_back( 1 ); v.push_back( 0 );
    v.push_back( 1 ); v.push_back( 1 ); v.push_back( 0 );
    v.push_back( 2 ); v.push_back( 1 ); v.push_back( 0 );
    v.push_back( 1 ); v.push_back( 2 ); v.push_back( 0 );

    f.push_back( 0 ); f.push_back( 1 ); f.push_back( 2 );
    f.push_back( 2 ); f.push_back( 1 ); f.push_back( 3 );
    f.push_back( 3 ); f.push_back( 1 ); f.push_back( 4 );
    f.push_back( 4 ); f.push_back( 5 ); f.push_back( 3 );
    f.push_back( 3 ); f.push_back( 5 ); f.push_back( 2 );

    m1.init( v, f );

    for( std::size_t i = 0, iEnd = v.size(); i < iEnd; i += 3 )
      assert( m2.add_vertex( v[i], v[i+1], v[i+2] ) == i/3 );
    for( std::size_t i = 0, iEnd = f.size(); i < iEnd; i += 3 )
      assert( m2.add_face( *reinterpret_cast<std::size_t(*)[3]>( &f[i] ) ) == i/3 );

    assert( m1.get_face_size() == m2.get_face_size() );
    assert( m1.get_vert_size() == m2.get_vert_size() );

    m1.verify();
    m2.verify();

    m2.delete_face( 0 );
    m2.verify();

    m2.delete_face( 0 ); // Was face 4 originally
    m2.verify();

    m2.delete_face( 0 ); // Was face 3 originally
    m2.verify();

    m2.delete_face( 0 ); // Was face 2 originally
    m2.verify();

    m2.delete_face( 0 ); // Was face 1 originally
    m2.verify();

    assert( m2.get_face_size() == 0 );

    m2.add_face( 0, 1, 2 );
    m2.verify();

    m2.split_face_center( 0 );
    m2.verify();

    assert( m2.get_face_size() == 3 );

    m2.clear();
    m2.add_vertex( 0, 0, 0 );
    m2.add_vertex( 1, 0, 0 );
    m2.add_vertex( 0, 1, 0 );
    m2.add_face( 0, 1, 2 );
    m2.split_boundary_edge( m2.find_twin( 0, 1 )->twin );
    m2.verify();

    assert( m2.get_face_size() == 3 );

    m3.add_vertex( -1, 0, 0 );
    m3.add_vertex( 1, 0, 0 );
    m3.add_vertex( 0, 1, 0 );
    m3.add_vertex( 0, -1, 0 );
    m3.add_vertex( -1, 1, 0 );
    m3.add_vertex( -1, -1, 0 );
    m3.add_vertex( 1, 1, 0 );
    m3.add_vertex( 1, -1, 0 );

    m3.add_face( 0, 1, 2 );
    m3.add_face( 0, 2, 4 );
    m3.add_face( 0, 4, 5 );
    m3.add_face( 0, 5, 3 );
    m3.add_face( 0, 3, 1 );
	
    m3.add_face( 1, 3, 7 );
    m3.add_face( 1, 7, 6 );
    m3.add_face( 1, 6, 2 );
    m3.verify();

    m3.flip_edge( *m3.find_edge( 1, 0 ) );
    m3.verify();
    m3.flip_edge( *m3.find_edge( 2, 3 ) );
    m3.verify();

    std::size_t newVert = m3.collapse_edge( m3.find_edge( 1, 0 )->twin );
    assert( newVert != HOLE_INDEX );
    m3.verify();

    /*m2.clear();
      m2.add_vertex( 0.010744, 0.483695, 0.298761 );
      m2.add_vertex( 0.010538, 0.484281, 0.305409 );
      m2.add_vertex( 0.014906, 0.48369, 0.304997 );
      m2.add_vertex( 0.006473, 0.484811, 0.30548 );
      m2.add_vertex( 0.010333, 0.484867, 0.312038 );
      m2.add_vertex( 0.004998, 0.485704, 0.314376 );
      m2.add_vertex( 0.010129, 0.485783, 0.323883 );
      m2.add_vertex( 0.016209, 0.484307, 0.313866 );
      m2.add_face( 7, 1, 4 );
      m2.add_face( 7, 6, 1 );
      m2.add_face( 6, 3, 1 );
      m2.add_face( 5, 3, 6 );
      m2.add_face( 5, 0, 3 );
      m2.add_face( 3, 2, 4 );
      m2.add_face( 3, 0, 2 );
      m2.add_face( 1, 3, 4 );

      for( std::size_t i = m2.get_face_size(); i > 0; --i ){
      m2.delete_face( i-1 );
      m2.verify();
      }*/

    m1.clear();
    m1.add_vertex( -1, 0, 0 );
    m1.add_vertex( 1, 0, 0 );
    m1.add_vertex( 0, 1, 0 );
    m1.add_vertex( -2, 1, 0 );
    m1.add_vertex( -2, -1, 0 );
    m1.add_vertex( 0, -1, 0 );
    m1.add_vertex( 2, -1, 0 );
    m1.add_vertex( 2, 1, 0 );
    m1.add_vertex( -2, 2, 0 );
    m1.add_vertex( -3, 1, 0 );
    m1.add_vertex( -3, -1, 0 );
    m1.add_vertex( -2, -2, 0 );
    m1.add_vertex( 2, -2, 0 );
    m1.add_vertex( 3, -1, 0 );
    m1.add_vertex( 3, 1, 0 );
    m1.add_vertex( 2, 2, 0 );

    m1.add_face( 0, 1, 2 );
    m1.add_face( 0, 2, 3 );
    m1.add_face( 0, 3, 4 );
    m1.add_face( 0, 4, 5 );
    m1.add_face( 0, 5, 1 );
    m1.add_face( 1, 5, 6 );
    m1.add_face( 1, 6, 7 );
    m1.add_face( 1, 7, 2 );
    m1.add_face( 3, 2, 8 );
    m1.add_face( 3, 8, 9 );
    m1.add_face( 3, 9, 10 );
    m1.add_face( 3, 10, 4 );
    m1.add_face( 4, 10, 11 );
    m1.add_face( 4, 11, 5 );
    m1.add_face( 5, 11, 12 );
    m1.add_face( 5, 12, 6 );
    m1.add_face( 6, 12, 13 );
    m1.add_face( 6, 13, 14 );
    m1.add_face( 6, 14, 7 );
    m1.add_face( 7, 14, 15 );
    m1.add_face( 7, 15, 2 );
    m1.add_face( 2, 15, 8 );

    m1.verify();

    std::vector< std::size_t > faces;
    vface_iterator it;
    if( m1.init_iterator( it, 0 ) ){
      do{
	std::vector< std::size_t >::iterator itInsert = std::lower_bound( faces.begin(), faces.end(), m1.deref_iterator( it ), std::greater<std::size_t>() );
	if( itInsert == faces.end() || *itInsert != m1.deref_iterator( it ) )
	  faces.insert( itInsert, m1.deref_iterator( it ) );
      }while( m1.advance_iterator( it ) );
    }
    if( m1.init_iterator( it, 1 ) ){
      do{
	std::vector< std::size_t >::iterator itInsert = std::lower_bound( faces.begin(), faces.end(), m1.deref_iterator( it ), std::greater<std::size_t>() );
	if( itInsert == faces.end() || *itInsert != m1.deref_iterator( it ) )
	  faces.insert( itInsert, m1.deref_iterator( it ) );
      }while( m1.advance_iterator( it ) );
    }

    std::swap( faces[faces.size()-1], faces[faces.size()-2] );

    for( auto face : faces ){
      m1.delete_face( face );
      m1.verify();
    }
  }

  bool EditMesh::is_safe_addface( std::size_t v1, std::size_t v2, std::size_t v3 ) {

    std::set<std::size_t> vv;
    vv.insert( _vert_to_he[v1] );
    vv.insert( _vert_to_he[v2] );
    vv.insert( _vert_to_he[v3] );


    // if there's one disconnected vertex its safe
    // more than one is not
    if( vv.size() < 3 )
      return false;
    else if( vv.count(HOLE_INDEX) > 0 )
      return true;

    vv.clear();
    vv.insert(v1);
    vv.insert(v2);
    vv.insert(v3);

    int count = 0;
    vvert_iterator vit;
    init_iterator(vit, v1);
    do {
      if (vit.m_cur->vert == v2 ||
	  vit.m_cur->vert == v3) {
	count++;
	break;
      }
    } while( this->advance_iterator(vit) );

    init_iterator(vit, v2);
    do {
      if (vit.m_cur->vert == v1 ||
	  vit.m_cur->vert == v3) {
	count++;
	break;
      }
    } while( this->advance_iterator(vit) );

    // if two of the vertices have a next within the triplet, its safe
    if ( count > 1 )
      return true;

    // else unsafe
    return false;
    //return (v1 != HOLE_INDEX && v2 != HOLE_INDEX) ||
    //       (v2 != HOLE_INDEX && v3 != HOLE_INDEX) ||
    //       (v3 != HOLE_INDEX && v1 != HOLE_INDEX);
  }


  void EditMesh::example() {
    // iterate over all halfedges of the mesh
    for (size_t i = 0; i < _he_data.size(); ++i) {
      // const HalfEdge* he = &_he_data[i];
    }

    // iterate over all faces of the mesh
    for (size_t f = 0; f < _face_to_he.size(); ++f) {
      // one of its halfedges is _face_to_he[i];
    }

    // iterate over all vertices of the mesh
    for (size_t v = 0; v < _vert_to_he.size(); ++v) {
      // one of its halfedges is _vert_to_he[i];
    }

    // iterate over first ring of a vertex
    size_t vFrom = 42;
    vvert_iterator it;
    if( !this->init_iterator( it, vFrom ) ) { /* error, vertex doesn't exist or has no neighbors */ }
    do {
      // size_t vTo = deref_iterator(it);
    } while ( this->advance_iterator(it) );

    // similarly for other iterators (vface, fvert, fface)

    // get the 3 halfedges of a face
    {
      // size_t face = 42;
      // size_t he_idx = _face_to_he[face];
      // const HalfEdge* he1 = &_he_data[he_idx];
      // const HalfEdge* he2 = &next(*he1);
      // const HalfEdge* he3 = &next(*he2);
    }
  }

  /***************************************************************************\
                              Assignment 1
\***************************************************************************/

  // ---------------------------------------------------------------------
  // Helper functions
  // ---------------------------------------------------------------------

  // Count the number of edges, and define a map between
  // the half edges and the edges.
  void EditMesh::init_edge_maps()
  {

    // If the maps are already defined just return.
    if(_is_edge_map_defined) return;

    // Associate an edge number to each half edge and his twin.
    _he_to_edge.resize(_he_data.size());
    std::fill(_he_to_edge.begin(), _he_to_edge.end(), HOLE_INDEX);
    std::size_t n_edge=0;
	
    for (std::size_t i=0 ; i < _he_data.size() ; i++)
      {
        HalfEdge &he= _he_data[i];

        // If the half edge is not already visited, then
        // assign an edge to both it and its twin.
        if (_he_to_edge[i] == HOLE_INDEX)
	  {		
            _he_to_edge[i] = n_edge;
            if(he.twin != HOLE_INDEX) _he_to_edge[he.twin] = n_edge;		
            n_edge++;
	  }
      }

    // Now that we know the number of edges, for each edge find its
    // one of its corresponding half edges.
    // Also count the number of boundary edges.
    _edge_to_he.resize(n_edge);
    size_t n_bdry_edges = 0;

    for (std::size_t ihe=0 ; ihe < _he_data.size() ; ihe++)
      {
        const HalfEdge &he= _he_data[ihe];
        std::size_t iedge = _he_to_edge[ihe];
        assert(iedge != HOLE_INDEX);

        // The edge must point to the internal half edge in case of a
        // boundary edge.
        const bool dont_point_to_me = (he.face == HOLE_INDEX);

        if (!dont_point_to_me)
	  {
            _edge_to_he[iedge] = ihe;
	  }
        else //Also count the number of bdry vertices
	  {
            ++n_bdry_edges;
	  }
      }

    // Now loop over all the edges and create the edge to bdry_edge mapping.
    // If an edge is not on the boundary, its map will store the invalid
    // value of HOLE_INDEX.
    _edge_to_bedge.resize(n_edge);
    _bedge_to_edge.resize(n_bdry_edges);
	
    size_t bdry_edge = 0;
    for (size_t edge=0 ; edge < n_edge ; edge++)
      {
        const size_t ihe1 = _edge_to_he[edge];
        const size_t ihe2 = _he_data[ihe1].twin;

        // If the edge was on the boundary
        assert(_he_data[ihe1].face != HOLE_INDEX);
        if (_he_data[ihe2].face == HOLE_INDEX)
	  {
            _edge_to_bedge[edge] = bdry_edge;
            _bedge_to_edge[bdry_edge] = edge;
            ++bdry_edge;
	  }
        else
	  {
            _edge_to_bedge[edge] = HOLE_INDEX;
	  }
      }


    // Announce that the maps are now defined.
    _is_edge_map_defined = true;
	
  }

  // get the edge numbers for the face
  void EditMesh::get_edges_for_face( const size_t face, std::size_t edges[3] ) const
  {
    assert(_is_edge_map_defined);

    // We must guarantee that the first edge has the first vert
    // that getIndicesForFace returns as a member.
    std::size_t verts[3];
    getIndicesForFace(face, verts);

    // get the half edges starting from the one containing verts[0]
    const std::size_t& ihe1 = _face_to_he[face];
    const std::size_t& ihe2 = _he_data[ihe1].next;
    const std::size_t& ihe3 = _he_data[ihe2].next;

    // find the edge corresponding to each half edge.
    edges[0] = _he_to_edge[ihe1];	
    edges[1] = _he_to_edge[ihe2];
    edges[2] = _he_to_edge[ihe3];

    // make sure that the first half-edge contains point 1
    assert(_he_data[ihe1].vert == verts[0] );
		
  }

  void EditMesh::get_verts_for_edge( const size_t edge, std::size_t verts[2] ) const
  {
    const std::size_t ihe = _edge_to_he[edge];
    // the _edge_to_he map is created in a way that the 
    // half_edge should never be adjacent to a  hole.
    assert(_he_data[ihe].face != HOLE_INDEX);

    verts[0] = _he_data[ihe].vert;
    verts[1] = _he_data[_he_data[ihe].twin].vert;		
  }
  
  // ---------------------------------------------------------------------
  //  Uniform refinement by adding a vertex in the middle of each edge.
  //  Increases the number of faces by a factor of 4.
  // ---------------------------------------------------------------------

  void EditMesh::uniformly_refine_face_based(const bool use_init_adj)
  {
    init_edge_maps();
    const std::size_t n_edge = _edge_to_he.size();
    VecOfVerts newverts(n_edge);

    for(std::size_t e=0 ; e < n_edge ; e++)
      {
        std::size_t verts[2];
        get_verts_for_edge(e, verts);
        newverts[e] = 0.5 * (_verts[verts[0]] + _verts[verts[1]] );
      }

    uniformly_refine_face_based(newverts, use_init_adj);
  }

  // This version uses the init_adjacency to build the adjacency list
  // from scratch.
  void EditMesh::uniformly_refine_face_based(const VecOfVerts& newverts, const bool should_use_init_adjacency)
  {

    // Make sure the size for the new verts are correct
    this->init_edge_maps();
    const std::size_t n_edges = _edge_to_he.size();
    const std::size_t n_verts = get_vert_size();
    const std::size_t n_faces = get_face_size();
    const std::size_t n_hes = _he_data.size();
    const std::size_t n_bdry_edges = _bedge_to_edge.size();
	
    assert(newverts.size() == n_edges);

    if (should_use_init_adjacency)
      {
        // Reinitialize a new mesh.
        std::vector<double> verts;
        std::vector<std::size_t> f2vert;

        // Add all the current points to m2.
        for (auto v=_verts.begin(); v != _verts.end(); v++)
	  {
            verts.push_back(v->x());
            verts.push_back(v->y());
            verts.push_back(v->z());
	  }
        // Add all the new verts
        for (auto v=newverts.begin(); v != newverts.end(); v++)
	  {
            verts.push_back(v->x());
            verts.push_back(v->y());
            verts.push_back(v->z());
	  }

        // Now for each old face, create four new faces.
        for (std::size_t face = 0 ; face < n_faces ; face++)
	  {
            std::size_t verts[3];
            std::size_t edges[3];
            get_edges_for_face(face, edges);
            getIndicesForFace(face, verts);
			
            // three corner points
            const std::size_t &a = verts[0];
            const std::size_t &b = verts[1];
            const std::size_t &c = verts[2];

            // three mid points
            const std::size_t e = n_verts + edges[1];
            const std::size_t f = n_verts + edges[2];
            const std::size_t g = n_verts + edges[0];

            // create the new four faces
            f2vert.push_back(a); f2vert.push_back(g); f2vert.push_back(f);
            f2vert.push_back(g); f2vert.push_back(b); f2vert.push_back(e);
            f2vert.push_back(f); f2vert.push_back(e); f2vert.push_back(c);
            f2vert.push_back(g); f2vert.push_back(e); f2vert.push_back(f);
	  }

        // Create the new mesh.
        this->init(verts, f2vert, true);
      } // End of use_init_adjacency

    // Not using init_adj
    else
      {
        // We will create a new mesh and initialize its connectivity.
        // Then we will swap its connectivity with the current mesh.
        EditMesh mesh;

        /*
         * We know the number of new entities in advance.
         */
        assert( _he_data.size() / 2 == n_edges);
        const size_t n_verts_new = n_verts + n_edges;
        const size_t n_faces_new = 4*n_faces;
        const size_t n_he_new = 2 * _he_data.size() + 6 * n_faces;

        mesh._he_data.resize(n_he_new);
        mesh._vert_to_he.resize(n_verts_new);
        mesh._face_to_he.resize(n_faces_new);

        /*
         * Add the point locations.
         */

        // resize
        mesh._verts.resize(0);
        mesh._verts.reserve(n_verts_new);
		
        // Add all the current points to mesh.
        for (auto v=_verts.begin(); v != _verts.end(); v++)
	  mesh._verts.push_back(*v);
		
        // Add all the new verts
        for (auto v=newverts.begin(); v != newverts.end(); v++)
	  mesh._verts.push_back(*v);

        /*
         * Loop over all faces and find their contribution to the connectivity.
         * The numberings in this part of the code is based on an arbitrary
         * numbering I have drawn in my notebook :) !
         */

        for (size_t face = 0 ; face < n_faces ; face++)
	  {
            /*
             * Get all the involved old and new numberings for book keeping.
             */
	   
            // Old half edges
            const size_t he0 = _face_to_he[face];
            const size_t he1 = _he_data[he0].next;
            const size_t he2 = _he_data[he1].next;
            const size_t he3 = _he_data[he0].twin;
            const size_t he4 = _he_data[he2].twin;
            const size_t he5 = _he_data[he1].twin;

            // Old vertices
            const size_t v0 = _he_data[he0].vert;
            const size_t v1 = _he_data[he1].vert;
            const size_t v2 = _he_data[he2].vert;

            // New faces
            const size_t sf0 = face * 4 + 0;
            const size_t sf1 = face * 4 + 1;
            const size_t sf2 = face * 4 + 2;
            const size_t sf3 = face * 4 + 3;

            // New verts
            const size_t sv0 = v0;
            const size_t sv1 = v1;
            const size_t sv2 = v2;
            const size_t sv3 = n_verts + _he_to_edge[he1];
            const size_t sv4 = n_verts + _he_to_edge[he2];
            const size_t sv5 = n_verts + _he_to_edge[he0];

            // New half edges
            const size_t she0 = he0*2 + 0;
            const size_t she1 = he0*2 + 1;
            const size_t she2 = he1*2 + 0;
            const size_t she3 = he1*2 + 1;
            const size_t she4 = he2*2 + 0;
            const size_t she5 = he2*2 + 1;
            const size_t she6 = n_hes * 2 + face * 6 + 0;
            const size_t she7 = n_hes * 2 + face * 6 + 1;
            const size_t she8 = n_hes * 2 + face * 6 + 2;
            const size_t she9 = n_hes * 2 + face * 6 + 3;
            const size_t she10 = n_hes * 2 + face * 6 + 4;
            const size_t she11 = n_hes * 2 + face * 6 + 5;
            const size_t she12 = he3*2 + 0;
            const size_t she13 = he3*2 + 1;
            const size_t she14 = he4*2 + 0;
            const size_t she15 = he4*2 + 1;
            const size_t she16 = he5*2 + 0;
            const size_t she17 = he5*2 + 1;
	    				
            /*
             * Setup the connectivity.
             */
			
            // Half edges

            // -------------------- she0
            mesh._he_data[she0].next = she9;
            mesh._he_data[she0].twin = she13;
            mesh._he_data[she0].face = sf0;
            mesh._he_data[she0].vert = sv0;

            // she1
            mesh._he_data[she1].next = she2;
            mesh._he_data[she1].twin = she12;
            mesh._he_data[she1].face = sf1;
            mesh._he_data[she1].vert = sv5;			

            // she2
            mesh._he_data[she2].next = she11;
            mesh._he_data[she2].twin = she17;
            mesh._he_data[she2].face = sf1;
            mesh._he_data[she2].vert = sv1;			

            // she3
            mesh._he_data[she3].next = she4;
            mesh._he_data[she3].twin = she16;
            mesh._he_data[she3].face = sf2;
            mesh._he_data[she3].vert = sv3;			

            // she4
            mesh._he_data[she4].next = she10;
            mesh._he_data[she4].twin = she15;
            mesh._he_data[she4].face = sf2;
            mesh._he_data[she4].vert = sv2;			

            // she5
            mesh._he_data[she5].next = she0;
            mesh._he_data[she5].twin = she14;
            mesh._he_data[she5].face = sf0;
            mesh._he_data[she5].vert = sv4;			

            //----------------------- she6
            mesh._he_data[she6].next = she7;
            mesh._he_data[she6].twin = she9;
            mesh._he_data[she6].face = sf3;
            mesh._he_data[she6].vert = sv4;			

            // she7
            mesh._he_data[she7].next = she8;
            mesh._he_data[she7].twin = she11;
            mesh._he_data[she7].face = sf3;
            mesh._he_data[she7].vert = sv5;			

            // she8
            mesh._he_data[she8].next = she6;
            mesh._he_data[she8].twin = she10;
            mesh._he_data[she8].face = sf3;
            mesh._he_data[she8].vert = sv3;			

            // ---------------------- she9
            mesh._he_data[she9].next = she5;
            mesh._he_data[she9].twin = she6;
            mesh._he_data[she9].face = sf0;
            mesh._he_data[she9].vert = sv5;			

            // she10
            mesh._he_data[she10].next = she3;
            mesh._he_data[she10].twin = she8;
            mesh._he_data[she10].face = sf2;
            mesh._he_data[she10].vert = sv4;			

            // she11
            mesh._he_data[she11].next = she1;
            mesh._he_data[she11].twin = she7;
            mesh._he_data[she11].face = sf1;
            mesh._he_data[she11].vert = sv3;			

#ifdef USE_PREV
            mesh._he_data[she0].prev = she5;
            mesh._he_data[she1].prev = she11;
            mesh._he_data[she2].prev = she1;
            mesh._he_data[she3].prev = she10;
            mesh._he_data[she4].prev = she3;
            mesh._he_data[she5].prev = she9;
            mesh._he_data[she6].prev = she8;
            mesh._he_data[she7].prev = she6;
            mesh._he_data[she8].prev = she7;
            mesh._he_data[she9].prev = she0;
            mesh._he_data[she10].prev = she4;
            mesh._he_data[she11].prev = she2;
#endif /*USE_PREV*/
			
            // ----------------------- Faces
            mesh._face_to_he[sf0] = she0;
            mesh._face_to_he[sf1] = she1;
            mesh._face_to_he[sf2] = she3;
            mesh._face_to_he[sf3] = she6;
						
            // Vertices. Some stuff will be done more than once, but it
            // will not cause any problem.
            mesh._vert_to_he[sv0] = she0;
            mesh._vert_to_he[sv1] = she2;
            mesh._vert_to_he[sv2] = she4;
            mesh._vert_to_he[sv3] = she3;
            mesh._vert_to_he[sv4] = she5;
            mesh._vert_to_he[sv5] = she1;		
			
	  } // End of loop over faces
		

        /*
         * Loop over all boundary HE's.
         */
        for (size_t bedge=0 ; bedge<n_bdry_edges ; bedge++)
	  {
            /*
             * Get all the involved old and new numberings for book keeping.
             */

            // Old edge
            const size_t edge = _bedge_to_edge[bedge];
	   
            // Old half edges
            const size_t he1 = _he_data[ _edge_to_he[edge] ].twin;
            const size_t he2 = _he_data[he1].next;
            const size_t he4 = _he_data[he1].twin;

            // Old vertices
            const size_t v0 = _he_data[he1].vert;

            // New verts
            const size_t sv0 = v0;
            const size_t sv2 = n_verts + edge;

            // New half edges
            const size_t she1 = he1*2 + 0;
            const size_t she2 = he1*2 + 1;
            const size_t she3 = he2*2 + 0;
            const size_t she5 = he4*2 + 0;
            const size_t she6 = he4*2 + 1;

            /*
             * Setup the connectivity.
             */
            // she1
            mesh._he_data[she1].next = she2;
            mesh._he_data[she1].twin = she6;
            mesh._he_data[she1].face = HOLE_INDEX;
            mesh._he_data[she1].vert = sv0;

            // she2
            mesh._he_data[she2].next = she3;
            mesh._he_data[she2].twin = she5;
            mesh._he_data[she2].face = HOLE_INDEX;
            mesh._he_data[she2].vert = sv2;

#ifdef USE_PREV
            const size_t he0 = _he_data[he1].prev;
            const size_t he5 = _he_data[he0].twin;
            const size_t she0 = he0*2 + 1;
            mesh._he_data[she1].prev = she0;
            mesh._he_data[she2].prev = she1;
#endif /*USE_PREV*/
			
	  } // End of for over boundary edges

        /*
         * Devour the new data
         */
        _he_data.swap(mesh._he_data);
        _face_to_he.swap(mesh._face_to_he);
        _vert_to_he.swap(mesh._vert_to_he);
        _verts.swap(mesh._verts);
        _is_edge_map_defined = false;
		
      } // End of !should_use_init_adjacency
	
  } // All done

  // ---------------------------------------------------------------------
  //  Uniform refinement by adding a vertex in the middle of each face,
  //  and then flipping each old edge.
  //  Increases the number of faces by a factor of 3.
  // ---------------------------------------------------------------------

  void EditMesh::uniformly_refine_edge_based()
  {
    init_edge_maps();
    VecOfVerts face_verts(get_face_size());
    VecOfVerts bedge_verts(_bedge_to_edge.size()*2);
	

    for(std::size_t f=0 ; f < get_face_size() ; f++)
      {
        std::size_t verts[3];
        getIndicesForFace(f, verts);
        face_verts[f] = 1/3. * (_verts[verts[0]] + _verts[verts[1]] + _verts[verts[2]] );
      }

    for(std::size_t bedge=0 ; bedge < _bedge_to_edge.size() ; bedge++)
      {
        const std::size_t edge = _bedge_to_edge[bedge];
        const std::size_t ihe1 = _edge_to_he[edge];
        const std::size_t ihe2 = _he_data[ihe1].twin;
        const std::size_t v1 = _he_data[ihe1].vert;
        const std::size_t v2 = _he_data[ihe2].vert;

        bedge_verts[bedge*2 + 0] = 2/3. * _verts[v1] + 1/3. * _verts[v2];
        bedge_verts[bedge*2 + 1] = 1/3. * _verts[v1] + 2/3. * _verts[v2];		
      }

    uniformly_refine_edge_based(face_verts, bedge_verts);
  }


  // a helper function
  void EditMesh::find_correct_vert_on_bedge(const std::size_t ihe, std::size_t &ibdry_edge, uint &ivert)
  {

    assert(_he_data[ihe].face != HOLE_INDEX);
    assert(_he_data[_he_data[ihe].twin].face != HOLE_INDEX);
	
    // next he is one the boundary
    std::size_t he2 = _he_data[ihe].next;
    if(_he_data[ _he_data[he2].twin].face == HOLE_INDEX)
      {
        const std::size_t target_edge = _he_to_edge[he2];
        ibdry_edge = _edge_to_bedge[target_edge];
        ivert = 0;
        return;
      }

    // second half edge is on the boundary
    he2 = _he_data[he2].next;
    if(_he_data[ _he_data[he2].twin].face == HOLE_INDEX)
      {
        const std::size_t target_edge = _he_to_edge[he2];
        ibdry_edge = _edge_to_bedge[target_edge];
        ivert = 1;
        return;
      }

    // Must not reach here.
    assert(0);
  }

  void EditMesh::uniformly_refine_edge_based
  (const VecOfVerts& new_face_verts, const VecOfVerts& new_bedge_verts)
  {
    // debug
    // for (uint i = 0 ; i < _edge_to_he.size() ; i++)
    // {
    // 	const std::size_t ihe1 = _edge_to_he[i];
    // 	const std::size_t ihe2 = _he_data[ihe1].twin;
    // 	const std::size_t v1 = _he_data[ihe1].vert;
    // 	const std::size_t v2 = _he_data[ihe2].vert;
    // 	const std::size_t bedge = _edge_to_bedge[i];
    // 	std::cout << "he1:" << ihe1 << std::endl
    // 			  << "he2:" << ihe2 << std::endl
    // 			  << "v1:" << v1 << std::endl
    // 			  << "v2:" << v2 << std::endl
    // 			  << "bedge:" << bedge << std::endl
    // 			  << std::endl;			
    // }

    // Make sure the size for the new verts are correct
    // There should be one vertex per face and two vertices per
    // boundary edge.
    this->init_edge_maps();
    const std::size_t n_edges = _edge_to_he.size();
    const std::size_t n_bdry_edge = _bedge_to_edge.size();
    const std::size_t n_verts = get_vert_size();
    const std::size_t n_faces = get_face_size();

    assert(new_face_verts.size() == n_faces);
    assert(new_bedge_verts.size() == n_bdry_edge*2); 
	
    // Reinitialize a new mesh.
    std::vector<double> verts;
    std::vector<std::size_t> f2vert;

    // Are we in an even or odd refinement step?
    const bool is_even_step = (_n_edge_based_refinement_steps % 2 == 0);
    std::vector<std::size_t> face_to_gface(get_face_size());
    std::size_t n_good_faces = 0;
    
    // Add all the current points to m2.
    for (auto v=_verts.begin(); v != _verts.end(); v++)
      {
        verts.push_back(v->x());
        verts.push_back(v->y());
        verts.push_back(v->z());
      }
    // Add all the new face verts
    size_t iface = 0;
    n_good_faces=0;
    for (auto v=new_face_verts.begin(); v != new_face_verts.end(); v++)
      {
        if( is_even_step || (!this->isBoundaryFace(iface)) )
	  {			
            verts.push_back(v->x());
            verts.push_back(v->y());
            verts.push_back(v->z());
            face_to_gface[iface] = n_good_faces;
            n_good_faces++;
	  }
        iface++;
      }
    // Add all the new boundary edge verts if we are in an odd step.
    if(!is_even_step)
      {
        for (auto v=new_bedge_verts.begin(); v != new_bedge_verts.end(); v++)
	  {
            verts.push_back(v->x());
            verts.push_back(v->y());
            verts.push_back(v->z());
	  }
      }

    // Now for each old edge, create its corresponding new triangles.
    // We have three cases:
    for (std::size_t edge = 0 ; edge < n_edges ; edge++)
      {
        const std::size_t ihe1 = _edge_to_he[edge];
        const std::size_t ihe2 = _he_data[ihe1].twin;
        const std::size_t face1 = _he_data[ihe1].face;
        const std::size_t face2 = _he_data[ihe2].face;
		
        assert(face1 != HOLE_INDEX);

        // Internal edges
        if (face2 != HOLE_INDEX )
	  {

            const bool face1_on_bdry = this->isBoundaryFace(face1);
            const bool face2_on_bdry = this->isBoundaryFace(face2);
			

            // edge end points.
            const std::size_t a = _he_data[ihe1].vert;
            const std::size_t b = _he_data[ihe2].vert;

            // left and right vertices to create two triangles for this
            // edge. To be found.
            std::size_t e;
            std::size_t f;

            // If we are in an even refinement step, or if both adjacent
            // faces are not on the boundary.
            if ( is_even_step || ( (!face1_on_bdry) && (!face2_on_bdry) ))
	      {
                f = n_verts + face_to_gface[ _he_data[ihe1].face ];				
                e = n_verts + face_to_gface[ _he_data[ihe2].face ];
	      }			
            // If both adjacent faces are on the boundary.
            else if (face1_on_bdry && face2_on_bdry)
	      {
                std::size_t bedge;
                uint lv;

                find_correct_vert_on_bedge(ihe1, bedge, lv);
                f = n_verts + n_good_faces + bedge * 2 + lv;
				
                find_correct_vert_on_bedge(ihe2, bedge, lv);
                e = n_verts + n_good_faces + bedge * 2 + lv;				
	      }			
            // If only face1 is on the boundary.
            else if (face1_on_bdry)
	      {
                std::size_t bedge;
                uint lv;

                find_correct_vert_on_bedge(ihe1, bedge, lv);
                f = n_verts + n_good_faces + bedge * 2 + lv;

                e = n_verts + face_to_gface[ _he_data[ihe2].face ];
	      }			
            // If only face2 is on the boundary.
            else if (face2_on_bdry)
	      {
                std::size_t bedge;
                uint lv;

                f = n_verts + face_to_gface[ _he_data[ihe1].face ];
				
                find_correct_vert_on_bedge(ihe2, bedge, lv);
                e = n_verts + n_good_faces + bedge * 2 + lv;				  
	      }			
            // We should not reach here!
            else
	      {
                assert(0 && "Die a horrible death!!");
	      }

            // create the triangles
            f2vert.push_back(a); f2vert.push_back(e); f2vert.push_back(f);
            f2vert.push_back(e); f2vert.push_back(b); f2vert.push_back(f);
	  }

        // Boundary edges
        else
	  {
            std::size_t a,b,c;
            const std::size_t bedge = _edge_to_bedge[edge];

            if (is_even_step)
	      {
                a = _he_data[ihe1].vert;
                b = _he_data[ihe2].vert;
                c = n_verts + face_to_gface[ _he_data[ihe1].face ];
	      }
            else
	      {
                a = n_verts + n_good_faces + 2 * bedge;
                b = n_verts + n_good_faces + 2 * bedge + 1;
                std::size_t tmp;
                tmp = _he_data[ihe1].next;
                tmp = _he_data[tmp].next;
                c = _he_data[tmp].vert;
	      }
			
            // create the triangle
            f2vert.push_back(a); f2vert.push_back(b); f2vert.push_back(c);
	  }
      }
	
    // Create the new mesh.
    this->init(verts, f2vert,true);

    // Increase the counter
    _n_edge_based_refinement_steps++;
  }

  // ---------------------------------------------------------------------
  //  Sqrt3 subdivision - with boundary support
  // ---------------------------------------------------------------------
  void EditMesh::subdivide_sqrt3()
  {

    // To correctly handle the boundaries we need to know if we are in an
    // even or an odd step.
    this->init_edge_maps();
    const bool is_even_step = (_n_edge_based_refinement_steps % 2 == 0);
    const size_t n_boundary_edges = _bedge_to_edge.size();

    /*
     * Find the location of the new verticies that have to be placed
     * in the middle of each face.
     *
     * Note that when we are refining an odd step, the values we find
     * for boundary faces will be ignored which is totally ok.
     */
    VecOfVerts face_verts(get_face_size());

    for(std::size_t f=0 ; f < get_face_size() ; f++)
      {
        std::size_t verts[3];
        getIndicesForFace(f, verts);
        face_verts[f] = 1/3. * (_verts[verts[0]] + _verts[verts[1]] + _verts[verts[2]] );

      }

    /*
     * Find the new location of each old vertex which is not on
     * the boundary.
     */
    VecOfVerts vert_verts(get_vert_size());

    for (std::size_t v=0 ; v < get_vert_size() ; v++)
      {
        vvert_iterator viter;

        if (!this->init_iterator(viter, v))
	  assert(0 && "Lonely vertex.");

        // The next part of the code will take care of
        // boundary vertices, so just skip them.
        if (!this->reset_boundary_iterator(viter))
	  {
            // Find the number of neighbours and their sum.
            Eigen::Vector3d sum_neigh = Eigen::Vector3d::Zero();
            std::size_t n_neigh = 0; 

            do
	      {
                const size_t vn = deref_iterator(viter);
                sum_neigh += _verts[vn];
                n_neigh++;
	      } while (this->advance_iterator(viter));

            // Find the new location of the oldvert
            const double alpha = ( 4 - 2 * cos(2 * M_PI / n_neigh) ) / 9. ;
            vert_verts[v] = (1-alpha) * _verts[v] + alpha / n_neigh * sum_neigh;
	  } // end of interior
      }

    /*
     * Find the location of vertices located on the boundary.
     *
     * Note that we will be dealing with both the bedge_verts and
     * vert_verts arrays.
     * When we are in an even step, the bedge_verts will be ignored 
     * and the vert_verts located on the boundary will retain their
     * previous position.
     * However, in an odd step smoothing must occur according to the 
     * sqrt(3) paper.
     */
    VecOfVerts bedge_verts(n_boundary_edges*2);

    // Even step
    if(is_even_step)
      {
        // Just put the boundary verts in the same place as
        // they were before.
        for (size_t bedge=0 ; bedge<n_boundary_edges ; bedge++)
	  {
            const size_t edge = _bedge_to_edge[bedge];
            const size_t he = _edge_to_he[edge];
            const size_t end_vert = _he_data[ _he_data[he].next ].vert;
            vert_verts[end_vert] = _verts[end_vert];
	  }
      }

    // Odd step
    else
      {
        for (size_t bedge=0 ; bedge<n_boundary_edges ; bedge++)
	  {

            // get all the involved edge, boundary edge and half
            // edge numberings for book keeping.
            const size_t edge = _bedge_to_edge[bedge];
            const size_t he = _he_data[ _edge_to_he[edge] ].twin;
            const size_t heP = _he_data[he].next;
            const size_t edgeP = _he_to_edge[heP];
            const size_t bedgeP = _edge_to_bedge[edgeP];
            const size_t hePP = _he_data[heP].next;

            // get the index handle for all the points in the
            // stencil
            const size_t oR = _he_data[he  ].vert;
            const size_t oC = _he_data[heP ].vert;
            const size_t oL = _he_data[hePP].vert;

            const size_t nR = 2 * bedge + 0;
            const size_t nC = oC;
            const size_t nL = 2 * bedgeP + 1;

            // apply the smoothing.
            const double coeff = 1/double(27);
            bedge_verts[nL] = coeff * (10*_verts[oL] + 16*_verts[oC] + 1 *_verts[oR] );
            vert_verts [nC] = coeff * (4 *_verts[oL] + 19*_verts[oC] + 4 *_verts[oR] );
            bedge_verts[nR] = coeff * (1 *_verts[oL] + 16*_verts[oC] + 10*_verts[oR] );

            // const double coeff = 1/double(3);
            // bedge_verts[nL] = coeff * (1*_verts[oL] + 2*_verts[oC] + 0 *_verts[oR] );
            // vert_verts [nC] = coeff * (0 *_verts[oL] + 3*_verts[oC] + 0 *_verts[oR] );
            // bedge_verts[nR] = coeff * (0 *_verts[oL] + 2*_verts[oC] + 1*_verts[oR] );
	  }	
      }

    // Update the vertex coordinates, and then refine the mesh
    // according to the new locations.
    _verts.swap(vert_verts);	
    uniformly_refine_edge_based(face_verts, bedge_verts);
  } 


  // ---------------------------------------------------------------------
  //  loop subdivision - with boundary support
  // ---------------------------------------------------------------------
  void EditMesh::subdivide_loop(const bool should_use_init_adjacency)
  {
    // Find the location of the new verticies that have to be placed
    // in the middle of each edge.
    init_edge_maps();
    const size_t n_edge = _edge_to_he.size();
    VecOfVerts edgeverts(n_edge);

    for(std::size_t e=0 ; e < n_edge ; e++)
      {

        const std::size_t ihe1 = _edge_to_he[e];
        assert(_he_data[ihe1].face != HOLE_INDEX && "Sanity Check");
        const std::size_t ihe2 = _he_data[ihe1].next;
        const std::size_t ihe3 = _he_data[ihe2].next;
		
        const std::size_t ihe4 = _he_data[ihe1].twin;
		
        const std::size_t l = _he_data[ihe1].vert;
        const std::size_t r = _he_data[ihe2].vert;

        // For internal edges
        if (_he_data[ihe4].face != HOLE_INDEX)
	  {
            const std::size_t ihe5 = _he_data[ihe4].next;
            const std::size_t ihe6 = _he_data[ihe5].next;

            const std::size_t t = _he_data[ihe3].vert;
            const std::size_t b = _he_data[ihe6].vert;
		
            edgeverts[e] =
	      3/8. * (_verts[l] + _verts[r]) + 1/8. * (_verts[t] + _verts[b]);
	  }
        // Boundary edges
        else
	  {
            edgeverts[e] =	1/2. * (_verts[l] + _verts[r]);
	  }
      }

    // Find the new location of each old vertices.
    VecOfVerts vertverts(get_vert_size());

    for (std::size_t v=0 ; v < get_vert_size() ; v++)
      {
        vvert_iterator it;
		
        if (!this->init_iterator(it, v))
	  {
            std::cout <<"vertex: " << v << std::endl;
            assert( 0 && "Lonely vertex");
            vertverts[v] = _verts[v];
            break;
	  }

        // boundary vertex
        if (this->reset_boundary_iterator(it))
	  {
            const size_t l = deref_iterator(it);
            this->advance_iterator(it);
            const size_t r = deref_iterator(it);

            vertverts[v] = 1/8. * (_verts[r] + _verts[l]) + 3/4.*_verts[v]  ;

            // debug
            // std::cout << l << " " << v << "  " << r << std::endl;
			
	  } // End of boundary vertex
		
        // inner vertex
        else
	  {
            // Find the number of neighbours and their sum.
            Eigen::Vector3d sum_neigh = Eigen::Vector3d::Zero();
            std::size_t n_neigh = 0;

            do
	      {
                const size_t vn = deref_iterator(it);
                sum_neigh += _verts[vn];
                n_neigh++;
	      } while (this->advance_iterator(it));

            // Find the new location of the oldvert
            const double c1 = 3 + 2 * cos(2*M_PI/n_neigh);
            const double beta = 5/8. - c1*c1 / 64. ;
            vertverts[v] = (1-beta) * _verts[v] + beta / n_neigh * sum_neigh;
	
	  } // End of inner vertex
		
      }
    _verts.swap(vertverts);

    // Refine the mesh according to the new locations.
    uniformly_refine_face_based(edgeverts, should_use_init_adjacency);	
  }

  // ---------------------------------------------------------------------
  //  Butterfly subdivision
  //  Boundary support not implemented.
  // ---------------------------------------------------------------------
  void EditMesh::subdivide_butterfly()
  {
    // Find the location of the new verticies that have to be placed
    // in the middle of each edge.
    init_edge_maps();
    const size_t n_edge = _edge_to_he.size();
    VecOfVerts edgeverts(n_edge);

    for(std::size_t edge=0 ; edge < n_edge ; edge++)
      {

        /*
         * Traverse through the mesh and find all the needed vertices.
         */
		
        std::size_t iheA = _edge_to_he[edge];
        std::size_t iheB = _he_data[iheA].twin;
        assert(_he_data[iheA].face != HOLE_INDEX && "Boundary not supported");
        assert(_he_data[iheB].face != HOLE_INDEX && "Boundary not supported");

        // Triangle A
        const std::size_t a = _he_data[iheA].vert;

        iheA = _he_data[iheA].next;		
        std::size_t iheE = _he_data[iheA].twin;
        assert(_he_data[iheE].face != HOLE_INDEX && "Boundary not supported");
        const std::size_t b = _he_data[iheA].vert;
		
        iheA = _he_data[iheA].next;
        const std::size_t c = _he_data[iheA].vert;
        std::size_t iheF = _he_data[iheA].twin;
        assert(_he_data[iheF].face != HOLE_INDEX && "Boundary not supported");

        // Triangle F
        iheF = _he_data[iheF].next;
        iheF = _he_data[iheF].next;
        const std::size_t h =  _he_data[iheF].vert;

        // Triangle E
        iheE = _he_data[iheE].next;
        iheE = _he_data[iheE].next;
        const std::size_t g =  _he_data[iheE].vert;

        // Triangle B
        iheB = _he_data[iheB].next;		
        std::size_t iheD = _he_data[iheB].twin;
        assert(_he_data[iheD].face != HOLE_INDEX && "Boundary not supported");
    	
        iheB = _he_data[iheB].next;
        const std::size_t d = _he_data[iheB].vert;
        std::size_t iheC = _he_data[iheB].twin;
        assert(_he_data[iheC].face != HOLE_INDEX && "Boundary not supported");

        // Triangle C
        iheC = _he_data[iheC].next;
        iheC = _he_data[iheC].next;
        const std::size_t f =  _he_data[iheC].vert;

        // Triangle D
        iheD = _he_data[iheD].next;
        iheD = _he_data[iheD].next;
        const std::size_t e =  _he_data[iheD].vert;


        // Find the new coordinates
        edgeverts[edge] =
	  1/16. *
	  (
	   8 * (_verts[a] + _verts[b]) +
	   2 * (_verts[c] + _verts[d]) -
	   1 * (_verts[e] + _verts[f] + _verts[g] + _verts[h])
	   );
      }

    // Refine the mesh according to the new locations.
    // Old vertices are not displaced.
    uniformly_refine_face_based(edgeverts);	
  }


  /***************************************************************************\
                              Assignment 2
\***************************************************************************/
  void EditMesh::print_info(FILE *fl)
  {
    assert(fl);

    fprintf(fl , "%8s %8s %8s %8s %8s %8s \n","#HE", "vbegin", "vend", "face", "twin", "next"); ;
    for (size_t he=0 ; he < _he_data.size(); he++)
      {
        const int twin = _he_data[he].twin;
        const int vbeg = _he_data[he].vert;
        const int vend = _he_data[twin].vert;
        const int face = _he_data[he].face;
        const int next = _he_data[he].next;

        fprintf(fl , "%8d %8d %8d %8d %8d %8d \n", (int)he, vbeg, vend, face, twin, next); ;
      }
  }

  void EditMesh::print_he_verts(const std::size_t he, FILE* fl)
  {
    assert(fl);
    assert(he < _he_data.size());

    const size_t v1 = _he_data[he].vert;
    const size_t twin = _he_data[he].twin;
    const size_t v2 = _he_data[twin].vert;
    printf("verts: %d %d \n", (int)v1, (int)v2);
  }

  void EditMesh::init_simplification(const uint type)
  {
    // Should not call this twice.
    assert(_is_simplification_in_progress == 0);

    // is active arrays
    _is_vert_active.resize(_vert_to_he.size()); 
    _is_face_active.resize(_face_to_he.size());
    _is_he_active.resize(_he_data.size());

    std::fill(_is_he_active.begin(), _is_he_active.end(), true);
    std::fill(_is_face_active.begin(), _is_face_active.end(), true);
    std::fill(_is_vert_active.begin(), _is_vert_active.end(), true);
	
    // added and deleted entities per step
    _n_deleted_verts.resize(0);
    _n_deleted_hes.resize(0);
    _n_deleted_faces.resize(0);
    _n_added_verts.resize(0);
    _n_added_hes.resize(0);
    _n_added_faces.resize(0);
    _id_deleted_verts.resize(0);
    _id_deleted_hes.resize(0);
    _id_deleted_faces.resize(0);

    // Set the active entities
    _n_verts_active = _vert_to_he.size();
    _n_faces_active = _face_to_he.size();
    _n_hes_active = _he_data.size();
 
    // step count
    _n_simplification_steps = 0;

    // Set type of Simplification
    if( (type != 1) && (type != 2) ) throw "only method 1 and 2 are supported";
    _is_simplification_in_progress = type;

    // set up priority
    switch(_is_simplification_in_progress)
      {
      case 1:
	{
	  // Find the priority of each vertex.
	  _priadd_vd.resize(_vert_to_he.size());
	  std::fill(_priadd_vd.begin(), _priadd_vd.end(), _prival_vd.end());
	  for(size_t v=0 ; v < get_vert_size() ; v++) update_priority(v);
	  break;
	}
      case 2:
	{
	  // Find the Q matrix for each vertex
        
	  _vertexQ.resize(_vert_to_he.size(), Eigen::Matrix4d::Zero());

	  // Find K_p for every triangle and sum it up over its vertices.
	  geo::Plane triplane;
	  size_t triverts[3];
	  for (size_t face=0 ; face < get_face_size() ; face++)
	    {
	      this->getIndicesForFace(face, triverts);
	      triplane.n = geo::normal(_verts[triverts[0]], _verts[triverts[1]], _verts[triverts[2]]);
	      triplane.b = _verts[triverts[0]];
	      const double a = triplane.n.x();
	      const double b = triplane.n.y();
	      const double c = triplane.n.z();
	      const double d = -triplane.n.dot(triplane.b);

	      // debug: check the constant are correct.
	      // std::cout << a * _verts[triverts[0]].x() + b * _verts[triverts[0]].y() + c * _verts[triverts[0]].z() + d << std::endl;
            
	      for(uint triv=0 ; triv < 3 ; triv++)
		{
		  Eigen::Matrix4d& Q = _vertexQ[triverts[triv]];
		  Q(0,0) += a*a; Q(0,1) += a*b; Q(0,2) += a*c; Q(0,3) += a*d;
		  Q(1,0) += b*a; Q(1,1) += b*b; Q(1,2) += b*c; Q(1,3) += b*d; 
		  Q(2,0) += c*a; Q(2,1) += c*b; Q(2,2) += c*c; Q(2,3) += c*d; 
		  Q(3,0) += d*a; Q(3,1) += d*b; Q(3,2) += d*c; Q(3,3) += d*d; 
		}
	    }

	  // debug: check that Q's are correct
	  // for (size_t v=0 ; v < _vert_to_he.size() ; v++)
	  // {
	  //     Eigen::Vector4d xyz(_verts[v].x(), _verts[v].y(), _verts[v].z(), 1);
	  //     std::cout << xyz.transpose() * _vertexQ[v] * xyz << std::endl;
	  // }

	  // Update the priority of each half edge In order to prevent
	  // doing it twice and save some time, between the half edge
	  // and the twin do it for the one with a smaller id.
	  _priadd_ec.resize(_he_data.size(), _prival_ec.end());
	  for (size_t halfedge=0 ; halfedge < _he_data.size(); halfedge++)
	    {
	      if(halfedge > _he_data[halfedge].twin) continue;
	      update_priority(halfedge);
	    }
               
	  break;
	}
      default:
        break;
      }

  }

  void EditMesh::finalize_simplification()
  {
    switch (_is_simplification_in_progress)
      {
      case 1: break;
      case 2:
        for (auto it = _prival_ec.begin() ; it != _prival_ec.end() ; ++it)
	  {
            delete *it;
	  }
        _prival_ec.clear();
	_priadd_ec.clear();
        break;
      default: break;
      }
    _is_simplification_in_progress= 0;

  }
  void EditMesh::update_priority(const std::size_t id)
  {
    switch(_is_simplification_in_progress)
      {
      case 1:
	{
	  // Remove the item from the queue
	  auto it = _priadd_vd[id];
	  if(it != _prival_vd.end())
	    {
	      assert(it->second == id);
	      _prival_vd.erase(it);
	    }
	  _priadd_vd[id] = _prival_vd.end();
         
	  // If the item is active, calculate the new value
	  // and then insert it.
	  if(!_is_vert_active[id]) return;
        
	  // If the item is on boundary do not give it
	  // any priority
	  vvert_iterator viter;
	  bool flg;
	  flg = this->init_iterator(viter, id);
	  flg = flg && this->reset_boundary_iterator(viter);
	  if(flg) break;

	  vector<size_t> rv, hv;  get_ring_data(id, rv, hv);
	  //printf("number of guys in the ring: %d \n", (int)rv.size());

	  /*
	    Old Method: angle divided by 2PI
	  */
	  // double value = 0;
	  // for (uint i = 0 ; i < rv.size() ; i++)
	  // {
	  //     const uint j = (i+1) % rv.size();
	  //     Eigen::Vector3d v1 = _verts[rv[i]] - _verts[id];
	  //     Eigen::Vector3d v2 = _verts[rv[j]] - _verts[id];
	  //     const double cosangle =  v1.dot(v2) / sqrt( v1.dot(v1) * v2.dot(v2) );
	  //     assert(cosangle < 1);
	  //     assert(cosangle > -1);
	  //     const double angle = acos(cosangle);
	  //     assert(finite(angle));
	  //     value +=  angle;
	  //     // printf("Angle between %d,%d, and %d is: %lf \n", (int)rv[i], (int)id, (int)rv[j], 180. / M_PI * angle);                 
	  // }
	  // value /= -(2. * M_PI);

	  /*
	    Check for a crease. If there is one the criterion is distance
	    from the crease line. Otherwise, distance from the average plane.
	  */

	  // Count the number of crease
	  vector<size_t> vcrease;
	  for (uint i = 0 ; i < rv.size() ; i++)
	    {
	      const uint j = (i+1) % rv.size();
	      const uint k = (i+2) % rv.size();
	      Eigen::Vector3d n1 = geo::normal(_verts[id], _verts[rv[i]], _verts[rv[j]]);
	      Eigen::Vector3d n2 = geo::normal(_verts[id], _verts[rv[j]], _verts[rv[k]]);
	      const double ang = geo::angle_n(n1, n2);
	      // if(id==269) printf("VERTEX 269: angle %d %d %d is %lf , %lf\n ", rv[i], rv[j], rv[k], 180. /M_PI *ang, n1.dot(n2));
	      if(ang > _min_crease_angle) 
		{
		  vcrease.push_back(rv[j]);
		  //printf("vertex %d has creases: %lf %lf \n", (int)id, 180. / M_PI * ang, n1.dot(n2));
		}
	    }


	  // Give weight according to vertex type
	  switch (vcrease.size())
	    {
	    case 0:
	      {
		// no crease, average plane
		geo::Plane plane;
		geo::avg_plane(_verts, id, rv, plane);
		const double value = std::abs(geo::dist_from_plane(plane, _verts[id])); 
		auto result = _prival_vd.insert(PriorityPairVD(value, id));
		_priadd_vd[id] = result;
		break;
	      }

	    case 2:
	      {
		// edge vertex, distance from crease
		const double value = std::abs(geo::dist_from_line(_verts[id], _verts[vcrease[0]], _verts[vcrease[1]]) );
		auto result = _prival_vd.insert(PriorityPairVD(value, id));
		_priadd_vd[id] = result;
		break;
	      }

	    default:
	      // corner vertex, do not do anything
	      break;
	    } // End of switch(n_crease);
         
        
	  break;
	} // End of case(vertex decimation)
      case 2:
	{

	  /*
	    Do nothing for inactive half edges.  Their priority node
	    will be popped when it reaches the head of the queue in the
	    simplify function.
	  */
	  if(!_is_he_active[id]) return;
	  // debug
	  //printf("%d ", (int)id);
  
	  /*
	    Delete the current priority object inside the multimap.
	  */
        
	  const size_t idtwin = _he_data[id].twin;
	  assert(_is_he_active[idtwin]);

	  std::set<size_t> examiner;
	  if( (_priadd_ec[id] == _priadd_ec[idtwin]) && (_priadd_ec[id] != _prival_ec.end()) )
	    {
	      // Get the iterator and the priority
	      auto it = _priadd_ec[id];

	      // debug
	      // make sure that the half edges in the priority are consistent
	      PriorityEC *priority = *it;
	      examiner.clear();
	      examiner.insert(priority->he[0]);
	      examiner.insert(priority->he[1]);
	      assert(examiner.count(id));
	      assert(examiner.count(idtwin));

	      // Erase the priority and set the new address
	      delete priority;
	      _prival_ec.erase(it);
	      _priadd_ec[id] = _prival_ec.end();
	      _priadd_ec[idtwin] = _prival_ec.end();            
	    }
	  else
	    {
	      const size_t guys[] = {id, idtwin};
	      for (uint i = 0 ; i < 2 ; i++)
		{
		  const size_t id2 = guys[i];
		  const size_t idtwin2 = guys[(i+1) % 2];
                
		  if(_priadd_ec[id2] != _prival_ec.end())
		    {
		      auto it = _priadd_ec[id2];
                    
		      // debug
		      // make sure that the twin is not in here
		      PriorityEC *priority = *it;
		      examiner.clear();
		      examiner.insert(priority->he[0]);
		      examiner.insert(priority->he[1]);
		      assert(examiner.count(id2));
		      assert(!examiner.count(idtwin2));

		      // Remove this guy from the priority
		      const size_t wone = (priority->he[0] == id2 ? 0 : 1);
		      priority->he[wone] = HOLE_INDEX;
		      _priadd_ec[id2] = _prival_ec.end();
		    }
		}
            
	    }
	  assert(_priadd_ec[id]     == _prival_ec.end());
	  assert(_priadd_ec[idtwin] == _prival_ec.end());


	  /*
	    Find the new priority and push it into the container.
	    We have to do the following steps:
	    1- Find the new Q for the edge.
	    2- Solve for vbar.
	    3- find vbar.transpose() Q vbar
	  */
	  const size_t v1 = _he_data[id].vert;
	  const size_t v2 = _he_data[idtwin].vert;
	  PriorityEC *priority = new PriorityEC;
	  priority->he[0] = id;
	  priority->he[1] = idtwin;

	  // set the new Q matrix
	  priority->Q = _vertexQ[v1] + _vertexQ[v2];

	  // Solve for vbar.
	  double RHS[] = {-priority->Q(0,3), -priority->Q(1,3), -priority->Q(2,3)};
	  double LHS[3][3] =
            {
	      {priority->Q(0,0), priority->Q(0,1), priority->Q(0,2)},
	      {priority->Q(1,0), priority->Q(1,1), priority->Q(1,2)},
	      {priority->Q(2,0), priority->Q(2,1), priority->Q(2,2)}
            };
	  double LHSinv[3][3];
	  double result[4] = {0,0,0,1};
	  const bool invertible = geo::Invert3x3(LHS, LHSinv);
        
	  if(invertible)
	    {
	      geo::MultVec(LHSinv, RHS, result);
	      priority->vbar(0) = result[0];
	      priority->vbar(1) = result[1];
	      priority->vbar(2) = result[2];           
	    }
	  else
	    {
	      priority->vbar.array() = (_verts[v1] + _verts[v2]) / 2.;
	      result[0] = priority->vbar(0); 
	      result[1] = priority->vbar(1); 
	      result[2] = priority->vbar(2); 
	    }

	  // Find the new dv
	  Eigen::Vector4d vbar4(result);
	  priority->dv = (vbar4.transpose() * priority->Q * vbar4)(0);

	  // push the priority into the queue
	  auto it = _prival_ec.insert(priority);
	  _priadd_ec[id] = it;
	  _priadd_ec[idtwin] = it;

	  // debug: view the details of the edge
	  // std::cout << "verts: "<< v1 << ", " << v2 << std::endl; 
	  // std::cout << priority->Q << std::endl;
	  // std::cout << "Invertible: " << invertible << std::endl;
	  // std::cout << "Error: " << priority->dv << std::endl;
	  // std::cout << "VBar: " << priority->vbar.transpose() << std::endl;
	  // std::cout << "VNesfe: " << ((_verts[v1] + _verts[v2])/2.).transpose() << std::endl;        
	  // std::cout << std::endl;
        
	  break;
	}
      default:
        break;
      }
  }

  bool EditMesh::restore_last_simplification_step()
  {
    if(_n_simplification_steps==0) return false; 
    
    /*
     * Remember the vertices that are just added.
     */
    vector<size_t> bupverts;
    
    /*
     * POP back anything added in the previous step!
     */

    // verts
    for (uint i = 0 ; i < _n_added_verts.back() ; i++)
      {
        _is_vert_active.pop_back();
        _verts.pop_back();
        _vert_to_he.pop_back();		
      }

    // faces
    for (uint i = 0 ; i < _n_added_faces.back() ; i++)
      {
        _is_face_active.pop_back();
        _face_to_he.pop_back();
      }

    // half edges
    for (uint i = 0 ; i < _n_added_hes.back() ; i++)
      {
        _is_he_active.pop_back();
        _he_data.pop_back();		
      }

    /*
     * Find the new active entities
     */
    _n_verts_active = _n_verts_active - _n_added_verts.back() + _n_deleted_verts.back();
    _n_faces_active = _n_faces_active - _n_added_faces.back() + _n_deleted_faces.back();
    _n_hes_active = _n_hes_active - _n_added_hes.back() + _n_deleted_hes.back();

    /*
     * Set the deleted entities to active and fix their connectivities
     */
    // verts
    for (uint i = 0 ; i < _n_deleted_verts.back() ; i++)
      {
        const size_t vert = _id_deleted_verts.back();
        _id_deleted_verts.pop_back();
        _is_vert_active[vert] = true;
        bupverts.push_back(vert);
      }

    // faces
    for (uint i = 0 ; i < _n_deleted_faces.back() ; i++)
      {
        const size_t face = _id_deleted_faces.back();
        _id_deleted_faces.pop_back();
        _is_face_active[face] = true;
      }

    // half edges
    for (uint i = 0 ; i < _n_deleted_hes.back() ; i++)
      {
        const size_t he1 = _id_deleted_hes.back();
        _id_deleted_hes.pop_back();
        _is_he_active[he1] = true;

        // get the next and previous edges
        // const size_t face = _he_data[he1].face;
        const size_t twin = _he_data[he1].twin;
        const size_t vert = _he_data[he1].vert;

        // set the connectivity
        _he_data[twin].twin = he1;
        _vert_to_he[vert] = he1;

        // should be correct
        // _he_data[he2].face = face;
        // _he_data[he3].face = face;
        // _he_data[he3].next = he1;
        // _he_data[he2].next = he3;
        // _vert_to_he[vert] = he1;
      }


    /*
     * Pop back the n_deleted_entities
     */
    _n_deleted_hes.pop_back();
    _n_deleted_faces.pop_back();
    _n_deleted_verts.pop_back();

    _n_added_hes.pop_back();
    _n_added_faces.pop_back();
    _n_added_verts.pop_back();
	
    // reduce the number of simplification steps
    _n_simplification_steps--;

    /*
     * Find the priorities if modified.
     */
    switch(_is_simplification_in_progress)
      {
      case 1:
	{
	  assert(bupverts.size() == 1);
	  const size_t v = bupverts.back();
	  vector<size_t> rv, hv;
	  get_ring_data(v, rv, hv);

	  update_priority(v);
	  for (uint i=0 ; i < rv.size() ; i++) update_priority(rv[i]);

	  return true;
	  break;
	}
      case 2:
	{
	  // I am currently too lazy to implement this!
	  assert(0 && "Not implemented");
	  throw "Not implemented";
	  return false;    
	  break;
	}
      default: return false; break;
      }

  }

  void EditMesh::get_ring_data(const std::size_t vertex, std::vector<size_t>& ring_verts, std::vector<size_t>& ring_hes)
  {
    ring_verts.resize(0);
    ring_hes.resize(0);
    
    /*
     * Save the surrounding half edges.
     */
    vvert_iterator it;
    if(!this->init_iterator(it, vertex))
      {
        std::cout << "Isolated vertex does not have a ring" << std::endl;
        assert(0); throw "Isolated vertex does not have a ring";
      }
    // else if(this->reset_boundary_iterator(it))
    // {
    //     std::cout << "Boundary vertex does not have a ring" << std::endl;
    //     assert(0); throw;
    // }
    else do
	   {
             ring_hes.push_back(_he_data[it.m_cur->next].next);
             // debug
             // std::cout << "adding half edge " << ring_hes.back() << " to the ring" << std::endl;
	   }while(this->advance_iterator(it));

    /*
     * Add all the surrounding vertices.
     */
    for (int i=(int)ring_hes.size()-1 ; i >= 0  ; i--)
      {
        // Add the half edge to the map
        const size_t he1 = ring_hes[i];
        const size_t gv1 = _he_data[he1].vert;
		
        // Add the vertex to the ring
        ring_verts.push_back(gv1);
      }

  }

  bool EditMesh::analyze_vertex_for_removal(const std::size_t vertex, std::vector<size_t>& tri_verts)
  {
    // assert not boundary
    // If three vertices lie on an existing triangle: return false
    // ...
    // Else: return find_best_triangulation(ring, ring_verts, verts);

    // Find the ring vertices
    vector<size_t> rv, ringhes;
    get_ring_data(vertex, rv, ringhes);    

    /*
      Check for creases!
      If there are one or three creases do not remove this vertex.
      If there are two creases remove
    */
    vector<uint> rvc;
    for (uint i = 0 ; i < rv.size() ; i++)
      {
        const uint j = (i+1) % rv.size();
        const uint k = (i+2) % rv.size();
        Eigen::Vector3d n1 = geo::normal(_verts[vertex], _verts[rv[i]], _verts[rv[j]]);
        Eigen::Vector3d n2 = geo::normal(_verts[vertex], _verts[rv[j]], _verts[rv[k]]);
        const double ang = geo::angle_n(n1, n2);
        // if(id==269) printf("VERTEX 269: angle %d %d %d is %lf , %lf\n ", rv[i], rv[j], rv[k], 180. /M_PI *ang, n1.dot(n2));
        if(ang > _min_crease_angle) 
	  {
            rvc.push_back(j);
            //printf("vertex %d has creases: %lf %lf \n", (int)id, 180. / M_PI * ang, n1.dot(n2));
	  }
      }
    assert( (rvc.size()==2) || (rvc.size()==0) );

    // Find the average plane
    geo::Plane aveplane;
    geo::avg_plane(_verts, vertex, rv, aveplane);
    //debug
    // std::cout << "average plane, b: " << aveplane.b.transpose() << " n: " << aveplane.n.transpose() << std::endl;

    /*
      Find the triangulation
    */
    if(vertex==24) std::cout <<"ncrease: " << rvc.size() << std::endl;
    
    if(rvc.size()==0)
      {
	const bool ans = find_best_hole_triangulation(rv, tri_verts, aveplane);
	return ans;
      }
    else if(rvc.size()==2)
      {
        //std::cout << " removing crease \n ";
        vector<size_t> srv1, srv2, stv1, stv2;

        const uint i = rvc[0];
        const uint j = rvc[1];
        
        srv1.push_back(rv[i]);
        for (uint k = (i+1)%(rv.size()) ; k != j ; k = (k+1)%rv.size())  srv1.push_back(rv[k]);
        srv1.push_back(rv[j]);
        
        srv2.push_back(rv[j]);
        for (uint k = (j+1)%(rv.size()) ; k != i ; k = (k+1)%rv.size())  srv2.push_back(rv[k]);
        srv2.push_back(rv[i]);

        const bool ans1 = find_best_hole_triangulation(srv1, stv1, aveplane);
        const bool ans2 = find_best_hole_triangulation(srv2, stv2, aveplane);

        // if(vertex==24)
        // {
        //     for(uint l=0 ; l < srv1.size() ; l++) cout << srv1[l] << " ";
        //     cout << "\n";
        //     for(uint l=0 ; l < srv2.size() ; l++) cout << srv2[l] << " ";
        //     cout << "\n";
        //     std::cout <<"ans1: " << ans1  << " " << ans2<< std::endl;
        // }
        
        if(ans1 && ans2)
	  {
            tri_verts = stv1;
            for(uint l=0 ; l < stv2.size() ; l++) tri_verts.push_back(stv2[l]);
            return true;
	  }
        
        return false;
      }
    else
      {
        throw "Must either have 2 or 0 creases at a vertex if it has passed the priority queue.";
      }

    return false;
  }

  bool EditMesh::find_best_hole_triangulation(const std::vector<size_t>& rv, std::vector<size_t>& tri_verts, const geo::Plane& avep)
  {
    // Base case of divide and conquer
    if(rv.size() == 3)
      {

        // Make sure the normal does not flip too much
        // auto newnormal = geo::normal(_verts[rv[0]], _verts[rv[1]], _verts[rv[2]]);
        // if( geo::angle_n(newnormal, avep.n) > 60. / 180. * M_PI) return false;

        // somehow check connectivity (tet case)
        
        tri_verts = rv;
        return true;
      }

    /*
      Find all the valid dividing planes and sort them according
      to the min/max ratio
    */
    std::multimap<double, pair<uint, uint> > validdividers;
    bool success;
    double mindist, dividerlength, ratio;
    geo::Plane divp;
    int sameside;
    
    for (uint i = 0 ; i < rv.size()-2 ; i++ )
      {
        const uint upper_bound = (i==0 ? rv.size() - 1 : rv.size() - 0);
        
        for (uint j = i + 2 ; j < upper_bound; j++)
	  {
            // debug
            // printf("connecting %d to %d , i: %d, j: %d size: %d\n", rv[i] , rv[j], i  , j , rv.size()); // continue;

            /*
              Make sure the dividing plane does not create new connectivity
            */
            if (this->find_twin(rv[i], rv[j]) != NULL) continue;

            /*
              Find the dividing plane.
            */
            geo::div_plane(avep, _verts[rv[i]], _verts[rv[j]], divp);

            /*
              Make sure all points are on the same side of the div plane
              Find the maxdist 
	    */
            mindist = 1e6;
            sameside = -2;
            
            // right divided section
            success=true;
            for (uint k = i+1 ; k < j ; k++)
	      {
                const double dist = geo::dist_from_plane(divp, _verts[rv[k]]);
                if (sameside == -2) sameside = (dist > 0 ? 1 : -1);                
                //debug
                // cout << "vert " << rv[k] << " is on side " << sameside << " of " << rv[i] << " and " << rv[j] << std::endl;                    
                if (sameside * dist < 0)
		  {
                    // cout << "faild ... \n";
                    success=false;
                    break;
		  }
                mindist=MIN(mindist, ABS(dist));
	      }
            // cout << "srv1: ";
            // for(uint zz = 0 ; zz < srv1.size() ; zz++)
            // {
            // cout << srv1[zz] << " " ;
            // }
            // cout << endl;
            if(!success) continue;

            // left side
            sameside *= -1;
            success = true;
            for (uint k = (j+1)%(rv.size()) ; k != i ; k = (k+1)%rv.size())
	      {
                const double dist = geo::dist_from_plane(divp, _verts[rv[k]]);
                // cout << "vert " << rv[k] << " is on side " << (dist > 0 ? 1 : -1) << " of " << rv[i] << " and " << rv[j] << std::endl;
                if (sameside * dist < 0)
		  {
                    //   cout << "faild ... \n";
                    success=false;
                    break;
		  }
                mindist=MIN(mindist, ABS(dist));
	      }
            // cout << "srv2: ";
            // for(uint zz = 0 ; zz < srv2.size() ; zz++)
            // {
            //     cout << srv2[zz] << " " ;
            // }
            // cout << endl;
            if(!success) continue;

            /*
              Find the length of the dividing plane
            */
            Eigen::Vector3d p1,p2,p3;
            geo::project_on_plane(avep, _verts[rv[i]], p1);
            geo::project_on_plane(avep, _verts[rv[j]], p2);
            p3 = p2 - p1;
            dividerlength = sqrt(p3.dot(p3));
            ratio = mindist / dividerlength;
            
            if(ratio > 0.1) validdividers.insert(make_pair(-ratio, make_pair(i,j)));
	  }
      } // End of for to find candidate edges


    /*
      Find the solution by recursively trying out all the valid
      dividing planes.
    */
    
    // Loop over all possible dividing planes
    vector<size_t> srv1, srv2, stv1, stv2;
    for (auto it = validdividers.begin() ; it != validdividers.end(); it++)
      {
        const uint i = it->second.first;
        const uint j = it->second.second;
        
        srv1.resize(0);
        srv1.push_back(rv[i]);
        for (uint k = i+1 ; k < j ; k++)  srv1.push_back(rv[k]);
        srv1.push_back(rv[j]);
        // cout << "srv1: ";
        // for(uint zz = 0 ; zz < srv1.size() ; zz++)
        // {
        // cout << srv1[zz] << " " ;
        // }
        // cout << endl;

        srv2.push_back(rv[j]);
        for (uint k = (j+1)%(rv.size()) ; k != i ; k = (k+1)%rv.size())  srv2.push_back(rv[k]);
        srv2.push_back(rv[i]);
        // cout << "srv2: ";
        // for(uint zz = 0 ; zz < srv2.size() ; zz++)
        // {
        //     cout << srv2[zz] << " " ;
        // }
        // cout << endl;

        // Now try to subdivide
        success = find_best_hole_triangulation(srv1, stv1, avep); if(!success) continue;
        success = find_best_hole_triangulation(srv2, stv2, avep); if(!success) continue;
            
        tri_verts = stv1;
        for(uint l=0 ; l < stv2.size() ; l++) tri_verts.push_back(stv2[l]);
        // cout << "FINAL: ";
        // for(uint zz = 0 ; zz < tri_verts.size() ; zz++)
        // {
        //     cout << tri_verts[zz] << " " ;
        // }
        // cout << endl;
        return true;            
      } // End of for over candidate edges

    return false;
 
    // Temporary solution
    // connect all the vertices to a single one
    // tri_verts.resize(0);
    // for (uint i = 1 ; i < (rv.size()-1) ; i++)
    // {
    //     tri_verts.push_back(rv[0]);
    //     tri_verts.push_back(rv[i]);
    //     tri_verts.push_back(rv[i+1]);
    // }
    // return true;
  }

  void EditMesh::simplify_by_removing_vertex(const std::size_t vertex, const std::vector<size_t>& tri_verts)
  {
    // Typedef for a common map that we are going to use
    typedef std::pair< size_t, size_t> Pair;
    typedef std::pair< Pair  , size_t> VEdgePair; 
    typedef std::map < Pair  , size_t> VEdgeMap;

    // The vertices and half edges on the edge, and the
    // map between vertices and half edges.
    VEdgeMap vedgemap;
    std::vector<size_t> ring_verts;
    std::vector<size_t> ring_hes;

    // Number of new entities
    const uint n_tris_new = tri_verts.size()/3;
    uint n_tris_rmv=0;
    uint n_hes_new= 0;

    /*
     * Save all the half edges around this guy, and count the number
     * of triangles to be destroyed.
     */
    vvert_iterator it;
    if(!this->init_iterator(it, vertex))
      {
        std::cout << "Cannot remove isolated vertex" << std::endl;
        assert(0); throw "Cannot remove isolated vertex";
      }
    // else if(this->reset_boundary_iterator(it))
    // {
    //     std::cout << "Cannot remove boundary vertex" << std::endl;
    //     throw;
    //     assert(0); throw;
    // }
    else do
	   {
             n_tris_rmv++;
             ring_hes.push_back(_he_data[it.m_cur->next].next);
             // debug
             // std::cout << "adding half edge " << ring_hes.back() << " to the ring" << std::endl;
	   }while(this->advance_iterator(it));

    /*
     * Add outer half edges on the ring and their twins to the map. The
     * advance_iterator goes in clockwise direction, while we want to
     * save the vertices in counter-clockwise direction. So we have to
     * do inverse iteration.
     */
    for (int i=(int)ring_hes.size()-1 ; i >= 0  ; i--)
      {
        // Add the half edge to the map
        const size_t he1 = ring_hes[i];
        const size_t he2 = _he_data[he1].twin;
        const size_t gv1 = _he_data[he1].vert;
        const size_t gv2 = _he_data[he2].vert;

        // don't do this anymore, i.e., replace the inner hedges as well.
        // auto result = vedgemap.insert( VEdgePair(Pair(gv1, gv2), he1) ); assert(result.second);
        auto result =    vedgemap.insert( VEdgePair(Pair(gv2, gv1), he2) ); assert(result.second);
		
        // Add the vertex to the ring
        ring_verts.push_back(gv1);
        // debug
        // std::cout << "adding vertex " << ring_verts.back() << " to the ring" << std::endl;        
      }

    /*
     * Add the new half edges to the mesh and the map.
     * Do not set their connectivity yet.
     */
    for (uint tri=0 ; tri < n_tris_new ; tri++)
      {
        const size_t gv0 = tri_verts[tri*3+0];
        const size_t gv1 = tri_verts[tri*3+1];
        const size_t gv2 = tri_verts[tri*3+2];
        Pair vpair[3];
        vpair[0] = Pair(gv0, gv1);
        vpair[1] = Pair(gv1, gv2);
        vpair[2] = Pair(gv2, gv0);

        for (uint i = 0 ; i < 3 ; i++)
	  {
            auto it = vedgemap.find(vpair[i]);
            if(it == vedgemap.end())
	      {
                const size_t newhe = _he_data.size();
                _he_data.push_back(HalfEdge());
                n_hes_new++;
                _is_he_active.push_back(true);
                vedgemap.insert( VEdgePair(vpair[i], newhe) );
                //debug
                //std::cout << "created he: " << newhe << std::endl;
	      }
            // Else should never trigger
            // else
            // {
            // 	assert( vedgemap.find( Pair(vpair[i].second, vpair[i].first) ) != vedgemap.end() );
            // 	assert( vedgemap.find( Pair(vpair[i].second, vpair[i].first) )->second == _he_data[it->second].twin );
            // }
	  } // End of for over triangle edges
        // debug
        // std::cout << "created triangle consisting of: " << gv0 << ", " << gv1<< ", "<< gv2 << std::endl;                
      } // End of for over triangles

    /*
     * Count the number of entities to be removed, added and updated.
     */
    const size_t n_hes_rmv = n_tris_rmv * 3;
	
    // Update these numbers

    // Total steps
    _n_simplification_steps++;
	
    // Added
    _n_added_hes.push_back(n_hes_new);
    _n_added_verts.push_back(0);
    _n_added_faces.push_back(n_tris_new);

    // Removed
    _n_deleted_hes.push_back(n_hes_rmv);
    _n_deleted_verts.push_back(1);
    _n_deleted_faces.push_back(n_tris_rmv);

    // Active
    _n_verts_active = _n_verts_active - 1;
    _n_faces_active = _n_faces_active - n_tris_rmv + n_tris_new;
    _n_hes_active =   _n_hes_active   - n_hes_rmv + n_hes_new;

    /*
     * Remove the old stuff.
     */
	
    //Vert
    _is_vert_active[vertex] = false;
    _id_deleted_verts.push_back(vertex);

    // Faces and half edges
    for (uint i=0 ; i < ring_hes.size()  ; i++)
      {
        const size_t he1 = ring_hes[i];
        const size_t he2 = _he_data[he1].next;
        const size_t he3 = _he_data[he2].next;
				
        const size_t face = _he_data[he1].face;

        _is_he_active  [ he1 ] = false;
        _is_he_active  [ he2 ] = false;
        _is_he_active  [ he3 ] = false;
        _id_deleted_hes.push_back  (he1);
        _id_deleted_hes.push_back  (he2);
        _id_deleted_hes.push_back  (he3);

        _is_face_active[face] = false;
        _id_deleted_faces.push_back(face);
      }

    /*
     * Add the new face, and fix the new connectivity.
     *
     * The assertions will ensure that the triangulation is
     * truly a decomposition of the hole.
     */
    for (uint tri=0 ; tri < n_tris_new ; tri++)
      {
        // Get the vertices
        const size_t v1 = tri_verts[tri*3+0];
        const size_t v2 = tri_verts[tri*3+1];
        const size_t v3 = tri_verts[tri*3+2];

        // Get the half edges 
        const auto it1 = vedgemap.find(Pair(v1, v2)); assert(it1 != vedgemap.end());
        const auto it2 = vedgemap.find(Pair(v2, v3)); assert(it2 != vedgemap.end());
        const auto it3 = vedgemap.find(Pair(v3, v1)); assert(it3 != vedgemap.end());
        const auto it4 = vedgemap.find(Pair(v2, v1)); assert(it4 != vedgemap.end());
        const auto it5 = vedgemap.find(Pair(v3, v2)); assert(it5 != vedgemap.end());
        const auto it6 = vedgemap.find(Pair(v1, v3)); assert(it6 != vedgemap.end());

        const size_t he1 = it1->second;
        const size_t he2 = it2->second;
        const size_t he3 = it3->second;
        const size_t he4 = it4->second;
        const size_t he5 = it5->second;
        const size_t he6 = it6->second;

        // Add a new face, and its connectivity
        const size_t newface = _face_to_he.size();
        _face_to_he.push_back(he1);
        _is_face_active.push_back(true);
							  
        // Set the connectivity for the vertices
        _vert_to_he[v1] = he1;
        _vert_to_he[v2] = he2;
        _vert_to_he[v3] = he3;

        // Fix the half edge connectivity

        // The twins
        _he_data[he4].twin = he1;
        _he_data[he5].twin = he2;
        _he_data[he6].twin = he3;

        // The current guys
        _he_data[he1].face = newface;
        _he_data[he1].next = he2;
        _he_data[he1].twin = he4;
        _he_data[he1].vert = v1;

        _he_data[he2].face = newface;
        _he_data[he2].next = he3;
        _he_data[he2].twin = he5;
        _he_data[he2].vert = v2;

        _he_data[he3].face = newface;
        _he_data[he3].next = he1;
        _he_data[he3].twin = he6;
        _he_data[he3].vert = v3;		
      }

    /*
     * Change the priority of the vertices affected
     * No need to change the size of _priadd_vd because
     * the number of vertices has stayed the same.
     */
    update_priority(vertex);
    for( auto it= ring_verts.begin() ; it != ring_verts.end() ; ++it )
      {
        update_priority(*it);
      }
  }

  // Copy of collapse edge with my own modifications
  bool EditMesh::simplify_by_collapsing_edge(const std::size_t he, const Eigen::Vector3d *loc)
  {
    /*
      Record the id's before starting the process.
    */
    const bool _debug = false;

    assert(_is_simplification_in_progress == 2);    
    assert( he < _he_data.size() );
    assert( _is_he_active[he] );
    assert( _is_he_active[_he_data[he].twin] );
    if( (_he_data[he].face == HOLE_INDEX) || (_he_data[_he_data[he].twin].face == HOLE_INDEX) )
      throw "Cannot collapse a boundary edge" ;

    const HalfEdge& heBase = _he_data[he];
    const HalfEdge& heTwin = _he_data[heBase.twin];

    // We are going to delete the faces on either side of the chosen
    // edge, so we need to delete 3 HalfEdges and patch up the twin
    // links on the 4 bordering edges.
    std::size_t heBorder[4];
    heBorder[0] = _he_data[ heBase.next ].twin;
    heBorder[1] = _he_data[ _he_data[ heBase.next ].next ].twin;
    heBorder[2] = _he_data[ _he_data[ heTwin.next ].next ].twin;
    heBorder[3] = _he_data[ heTwin.next ].twin;

    // TODO: Relax this assertion. We should be able to collapse a
    // spike jutting into a hole.
    assert( ( _he_data[ heBorder[0] ].face != HOLE_INDEX || _he_data[ heBorder[1] ].face != HOLE_INDEX ) && "Cannot collapse an edge on a face with holes on either side." );
    assert( ( _he_data[ heBorder[2] ].face != HOLE_INDEX || _he_data[ heBorder[3] ].face != HOLE_INDEX ) && "Cannot collapse an edge on a face with holes on either side." );

    
    /*
      Check that the collapse does not screw the data structure
    */
    if( _he_data[ _he_data[ _he_data[ heBorder[1] ].next ].twin ].next == heBorder[0] ) return false;
    if( _he_data[ _he_data[ _he_data[ heBorder[2] ].next ].twin ].next == heBorder[3] ) return false;
    if(this->count_edge(_he_data[heBorder[3]].vert, _he_data[heBorder[0]].vert) > 1) return false;
    if(this->count_edge(heTwin.vert, heBase.vert) > 1) return false;
        
    
    // if(_n_simplification_steps == 9712)
    // {

    //     size_t bi[3], ti[3];
    //     getIndicesForFace(_he_data[heTwin.twin].face, bi);
    //     getIndicesForFace(_he_data[heBase.twin].face, ti);
        
    //     cout << "base tri: " << bi[0] << " " << bi[1] << " " << bi[2] << "\n";
    //     cout << "twin tri: " << ti[0] << " " << ti[1] << " " << ti[2] << "\n";
            
    //     cout << "base " ; print_he_verts(heTwin.twin);
    //     cout << "bdry0 "; print_he_verts(heBorder[0]);
    //     cout << "bdry1 "; print_he_verts(heBorder[1]);
    //     cout << "bdry2 "; print_he_verts(heBorder[2]);
    //     cout << "bdry3 "; print_he_verts(heBorder[3]);
    //     cout << "next, bdry1 "; print_he_verts(_he_data[ _he_data[ _he_data[ heBorder[1] ].next ].twin ].next);
    //     cout << "next, bdry3 "; print_he_verts(_he_data[ _he_data[ _he_data[ heBorder[2] ].next ].twin ].next);
        
    //     // std::cout << "nextnextnext: " << tmp << " " << _he_data[tmp].vert << " " << _he_data[twin].vert << std::endl;
    //     // std::cout << "border0: " << heBorder[0] << " " << _he_data[heBorder[0]].vert << " " <<  _he_data[btwin].vert << std::endl;
    // }

    // Capture the indices of things (2 faces & 6 half-edges) we want
    // to delete.
    std::size_t fToDelete[] = { heBase.face, heTwin.face };
    std::size_t heToDelete[] = { he, heBase.next, _he_data[ heBase.next ].next, heBase.twin, heTwin.next, _he_data[ heTwin.next ].next };
    std::size_t vertToDelete[] = {heBase.vert, heTwin.vert};
    
    // We can't be deleting border edges!
    for( auto i : heToDelete )
      {
        if( std::find( heBorder, heBorder + 4, i ) != heBorder + 4 )
	  return false;	
        //assert( std::find( heBorder, heBorder + 4, i ) == heBorder + 4 );
      }


    // Write down the edge and the neighbourhood before collapsing.
#ifndef NDEBUG

    if( _debug )
      {
        std::vector< std::set<std::size_t> > verts( 3 );

        verts[0].insert( heBase.vert );
        verts[0].insert( heTwin.vert );

        for( size_t i = 1; i < verts.size(); ++i )
	  {
            for( auto v : verts[i-1] )
	      {
                vvert_iterator it;
                this->init_iterator( it, v );
                do
		  {
                    verts[i].insert( this->deref_iterator( it ) );
		  }while( this->advance_iterator( it ) );
	      }
	  }

        std::vector<std::size_t> orderedVerts( verts.back().begin(), verts.back().end() );
        std::set<std::size_t> faces;

        std::vector< double > vpos;
        std::vector< std::size_t > finds;

        for( auto v : orderedVerts )
	  {
            vpos.push_back( _verts[v].x() ); vpos.push_back( _verts[v].y() ); vpos.push_back( _verts[v].z() );
            //std::clog << "m.add_vert( " << _verts[v].x() << ", " << _verts[v].y() << ", " << _verts[v].z() << " );" << std::endl;
	  }

        // Visit the 1-ring
        for( auto v : verts[1] )
	  {
            vface_iterator it;
            this->init_iterator( it, v );
            do{
	      if( this->deref_iterator( it ) != HOLE_INDEX && faces.find( this->deref_iterator( it ) ) == faces.end() )
                {
		  faces.insert( this->deref_iterator( it ) );

		  fvert_iterator itFace;
		  this->init_iterator( itFace, this->deref_iterator( it ) );

		  std::size_t f[3];
		  std::size_t i = 0;
		  do{
		    f[i++] = std::find( orderedVerts.begin(), orderedVerts.end(), this->deref_iterator( itFace ) ) - orderedVerts.begin();
		  }while( this->advance_iterator( itFace ) );

		  finds.push_back( f[0] ); finds.push_back( f[1] ); finds.push_back( f[2] );
		  //std::clog << "m.add_face( " << f[0] << ", " << f[1] << ", " << f[2] << " );" << std::endl;
                }	
            }while( this->advance_iterator( it ) );
	  }

        std::size_t base = std::find( orderedVerts.begin(), orderedVerts.end(), heBase.vert ) - orderedVerts.begin();
        std::size_t twin = std::find( orderedVerts.begin(), orderedVerts.end(), heTwin.vert ) - orderedVerts.begin();
        std::clog << "m.collapse_edge( " << base << ", " << twin << " );" << std::endl;

        EditMesh m;
        m.init( vpos, finds );
        MeshIO(m).write_vtk("edge_collapse.vtk");
      }
#endif

    
    /*
      Adjust the connectivities:
      half edge twin.
      half edge vertex.
      vertex to half edge.
      -- No: face to half edge!
    */
    
    // Add the new vertex and all associated data.
    const size_t newvert = _vert_to_he.size();
    _is_vert_active.push_back(true);
    _vertexQ.push_back(_vertexQ[heBase.vert] + _vertexQ[heTwin.vert]);
    if(loc) _verts.push_back(*loc);
    else _verts.push_back(0.5*(_verts[heBase.vert] + _verts[heTwin.vert]));
    _vert_to_he.push_back(HOLE_INDEX);

    std::size_t verts[] = { this->prev( heBase ).vert, newvert, this->prev( heTwin ).vert };

    // Half edge to vertex
    std::size_t heIt = this->twin(this->next(heBase)).next;
    std::size_t heEnd = heBase.twin;
    for( ; heIt != heEnd; heIt = this->twin( _he_data[heIt] ).next )
      { 
        assert( _he_data[heIt].vert == heTwin.vert );
        _he_data[heIt].vert = newvert;
      }
    heIt = this->twin(this->next(heTwin)).next;
    heEnd = heTwin.twin;
    for( ; heIt != heEnd; heIt = this->twin( _he_data[heIt] ).next )
      { 
        assert( _he_data[heIt].vert == heBase.vert );
        _he_data[heIt].vert = newvert;
      }

    // Vertex to half edge
    _vert_to_he[ verts[0] ] = (_he_data[ heBorder[0] ].face != HOLE_INDEX) ? heBorder[0] : _he_data[ heBorder[1] ].next;  
    _vert_to_he[ verts[1] ] = (_he_data[ heBorder[1] ].face != HOLE_INDEX) ? heBorder[1] : heBorder[2];
    _vert_to_he[ verts[2] ] = (_he_data[ heBorder[3] ].face != HOLE_INDEX) ? heBorder[3] : _he_data[ heBorder[2] ].next;

    // Half edge twin.
    _he_data[ heBorder[0] ].twin = heBorder[1];
    _he_data[ heBorder[1] ].twin = heBorder[0];
    _he_data[ heBorder[2] ].twin = heBorder[3];
    _he_data[ heBorder[3] ].twin = heBorder[2];

    /*
      Now update the add/remove data structure
    */
    
    // set n_deleted_entities
    _n_deleted_hes.push_back(6);
    _n_deleted_verts.push_back(2);
    _n_deleted_faces.push_back(2);
    
    // set n_added_entities
    _n_added_hes.push_back(0);
    _n_added_verts.push_back(1);
    _n_added_faces.push_back(0);

    // set n_active_entities
    _n_verts_active = _n_verts_active - 1;
    _n_faces_active = _n_faces_active - 2;
    _n_hes_active =   _n_hes_active   - 6;

    // set id_deleted, and is_active[] = false
    
    // verts
    _is_vert_active[vertToDelete[0]] = false;
    _is_vert_active[vertToDelete[1]] = false;
    _id_deleted_verts.push_back(vertToDelete[0]);
    _id_deleted_verts.push_back(vertToDelete[1]);

    // faces
    _is_face_active[fToDelete[0]] = false;
    _is_face_active[fToDelete[1]] = false;
    _id_deleted_faces.push_back(fToDelete[0]);
    _id_deleted_faces.push_back(fToDelete[1]);

    // he's
    for (uint i = 0 ; i < 6 ; i++)
      {
        _is_he_active[heToDelete[i]] = false;
        _id_deleted_hes.push_back(heToDelete[i]);
      }

    /*
      Increase number of simplification steps
    */
    _n_simplification_steps++;

    /*
      Update the priorities of the current half edges
    */
    vvert_iterator viter;
    this->init_iterator(viter, newvert);
    do
      {
        update_priority(viter.m_cur->twin);
      }while(this->advance_iterator(viter));

    return true;
  }

  bool EditMesh::simplify()
  {
    assert(_is_simplification_in_progress);

    switch(_is_simplification_in_progress)
      {
      case 1:
	{
	  vector<size_t> triverts;
	  bool success;
	  for(;;)
	    {
	      auto itbeg = _prival_vd.begin();
            
	      if(itbeg == _prival_vd.end()) return false;
	      assert(_is_vert_active[itbeg->second]);
	      //printf("Removing vertex %d ... ", (int) itbeg->second);

	      success=analyze_vertex_for_removal(itbeg->second, triverts);
	      if(!success)
		{
		  _priadd_vd[itbeg->second] = _prival_vd.end();
		  _prival_vd.erase(itbeg);
		  //printf("unsuccessful. \n");
		  continue;
		}
            
	      simplify_by_removing_vertex(itbeg->second, triverts);
	      //printf("success. \n");
	      return true;
	    }
        
	  break;
	}
      case 2:
	{
	  bool success;

	  // Try to find an edge to collapse
	  for(;;)
	    {
	      auto itbeg = _prival_ec.begin();

	      // Check if there are no more priorities
	      if(itbeg == _prival_ec.end()) return false;

	      // Check for invalid priorities
	      PriorityEC *p = *itbeg;            
	      if( (p->he[0]==HOLE_INDEX) || (p->he[1]==HOLE_INDEX) )
		{
		  if(p->he[0]!=HOLE_INDEX)
		    {
		      _priadd_ec[p->he[0]] = _prival_ec.end();
		      assert(!_is_he_active[p->he[0]]);
		    }
		  if(p->he[1]!=HOLE_INDEX)
		    {
		      _priadd_ec[p->he[1]] = _prival_ec.end();
		      assert(!_is_he_active[p->he[1]]);
		    }
		  delete p;
		  _prival_ec.erase(itbeg);
		  continue;
		}


	      // Try to collapse the edge
	      // printf("Trying to collapse %d %d \n", (int)_he_data[p->he[0]].vert, (int)_he_data[p->he[1]].vert);
	      assert(_is_he_active[p->he[1]]);
	      assert(_is_he_active[p->he[0]]);
	      success=simplify_by_collapsing_edge(p->he[0], &p->vbar);

	      _priadd_ec[p->he[1]] = _prival_ec.end();
	      _priadd_ec[p->he[0]] = _prival_ec.end();
	      delete p;
	      _prival_ec.erase(itbeg);

	      if(success) return true;
	    }
	}
      default:
        return false;
      }

    return false;
  }

  void EditMesh::updateBBox()
  {
    Eigen::AlignedBox3d bb;
    for (auto &v:_verts)
      {
	bb.extend(v);
      }
    bboxMin = bb.min();
    bboxMax = bb.max();
  }

} // End of hooshi
