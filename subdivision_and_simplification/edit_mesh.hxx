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

/*
 * edit_mesh.hxx
 *
 * Header file for the class EditMesh, which encapsulates a manifold
 * triangular mesh.
 */

#ifndef HOOSHI_EDIT_MESH_HXX
#define HOOSHI_EDIT_MESH_HXX

#include <Eigen/Core>

#include <cstdlib>
#include <vector>
#include <set>
#include <map>
#include <memory>

#include "geometry.hxx"

//#define USE_PREV

namespace hooshi {

  const std::size_t HOLE_INDEX = static_cast<std::size_t>( -1 );


  /****************************************************************************\
                                Half Edge
\****************************************************************************/
  struct HalfEdge
  {
    // Index of next half-edge in the loop.
    std::size_t next; 
#ifdef USE_PREV
    // Index of the previous half-edge in the loop. Alternatively: use
    // next->next if you have strictly triangle meshes.
    // On boundaries your life will be slightly more difficult without
    // this. You would need to use vvert_iterator.
    std::size_t prev; 
#endif
    // Index of half-edge that is in other face that shares this
    // edge. Open edges?
    std::size_t twin;
    // Index of vertex at the start of this edge.
    std::size_t vert;
    // Index of face to the "left" of this edge. Can be HOLE_INDEX
    // for boundary half edge
    std::size_t face; 
  };


  /****************************************************************************\
                                 Iterators
\****************************************************************************/

  // Verts around a vert
  class vvert_iterator
  {
  private:
    friend class EditMesh;
    const HalfEdge* m_cur;
    const HalfEdge* m_end;
  };

  // Faces around a vert
  class vface_iterator
  {
  private:
    friend class EditMesh;
    const HalfEdge* m_cur;
    const HalfEdge* m_end;
    const HalfEdge* m_next;
  };

  // Verts of a face
  class fvert_iterator
  {
  private:
    friend class EditMesh;
    const HalfEdge* m_cur;
    const HalfEdge* m_end;
  };

  // Faces around a face
  class fface_iterator
  {
  private:
    friend class EditMesh;
    const HalfEdge* m_cur;
    const HalfEdge* m_end;
  };


  /****************************************************************************\
                                 EditMesh
\****************************************************************************/
  class EditMesh
  {
  public:

    EditMesh();
    EditMesh *clone();

    /**
     * Initialize the mesh from existing data.  \param xyzPositions A
     * list of doubles, storing the vertex data interleaved
     * (ie. X1,Y1,Z1,X2,Y2,Z2,etc.)  \param triangleVerts A list of
     * vertex indices, where each run of 3 defines a triangle in the
     * mesh. (ie. T1.v1, T1.v2, T1.v3, T2.v1, T2.v2, T2.v3, etc.)
     */
    void init( const std::vector<double>& xyzPositions, const std::vector<std::size_t>& triangleVerts, const bool is_reinit = false );

    void clear();
    void prepare_for_reinit();

    std::size_t add_vertex( double x, double y, double z );
    std::size_t add_face( std::size_t v1, std::size_t v2, std::size_t v3 );
    std::size_t add_face( std::size_t (&v)[3] );
	
    void delete_face( std::size_t face );
	
    std::size_t split_face_center( std::size_t f, std::size_t (*pOutFaceIndices)[3] = NULL );

    /**
     * returns the position of a vertex
     * \param i the vertex index (not the he index)
     * \return Vector3d the vertex position values
     */
    Eigen::Vector3d get_vertex( std::size_t vertex ) const;
    Eigen::Vector3d get_vnormal( std::size_t vertex ) const;
    Eigen::Vector3d get_fnormal( std::size_t face ) const;

    void set_vertex( std::size_t i, const Eigen::Vector3d& v );

    void example();
    
    /************************************
     * Mesh modification functions
     ************************************/
  public:
	
    // Create maps from half edge indices to edge and boundary edge indices.
    void init_edge_maps();

    // Prints information about the mesh
    void print_info(FILE* = stdout);
    void print_he_verts(const std::size_t, FILE* = stdout);
    
    // -------------------------   Subdivision ----------------------------
	
    // Uniformly refines the mesh by breaking each face into four
    // faces.
    //
    // Input:
    //
    // new_edge_verts:
    // The location of vertices that have to be placed on each edge.
    // The version that does not take this argument was used for debugging
    // purposes and just finds the coordinates using linear interpolation.
    //
    // use_init_adjacency:
    // If set to true a version of algorithm will be be used which
    // uses the init_adjacency function and is therefore O(nlogn). If
    // set to false a different version will be used which uses the
    // initial adjacency to find the new adjacency and is O(n).
    void uniformly_refine_face_based(const bool use_init_adj);
    void uniformly_refine_face_based(const VecOfVerts& new_edge_verts, const bool should_use_init_adjacency = false);

    // Uniformly refines the mesh by using the same strategy used in
    // sqrt(3) subdivision.
    //
    // Input:
    //
    // new_face_verts:
    // The location of vertices that have to be placed inside each face.
    // The version that does not take this argument was used for debugging
    // purposes and just finds the coordinates using linear interpolation.
    //
    // new_bedge_verts:
    // The location of vertices that have to be placed on each boundary edge.
    // They are only used for odd refinement steps.
    // The version that does not take this argument was used for debugging
    // purposes and just finds the coordinates using linear interpolation.
    //
    // Note: This algorithm always uses init_adjacency.
    void uniformly_refine_edge_based();
    void uniformly_refine_edge_based(const VecOfVerts& new_face_verts, const VecOfVerts& new_bedge_verts);
	
    // Helper function for uniformly_refine_edge_based
    // Input:
    // ihe:        The index of a half edge inside a boundary face.
    // Output:
    // ibdry_edge: The index of the boundary edge that has to be connected
    //             to the half edge for refinement.
    // ivert:      Whether the first or the second vert on the boundary edge
    //             has to be connected to the center of the half edge.
    void find_correct_vert_on_bedge(const std::size_t ihe, std::size_t &ibdry_edge, uint &ivert);	


    // Subdivide the mesh using the sqrt(3) algorithm.
    void subdivide_sqrt3();

    // Subdivide the mesh using the loop algorithm.
    // Input:
    // should_use_init_adjacency: Indicates which algorithm has to be
    // used to find the adjacency list of the refined mesh.
    void subdivide_loop(const bool should_use_init_adjacency  = false);

    // Subdivide the mesh using the butterfly algorithm.
    // Note: Does not support meshes with open boundaries.
    void subdivide_butterfly();

    // -------------------------  Simplification  ----------------------------

 
    // Prepare the mesh to start the simplification process.
    // input: type: 1 for vertex decimation, 2 for quartic edge collapse.
    void init_simplification(const uint type);

    uint  is_simplification_in_progress() const {return  _is_simplification_in_progress;}
    // Free the memory used while simplifying. Not calling this causes
    // memory leak.
    void finalize_simplification();

    // Update the priority of an entity.
    // Input:
    // id: id of vertex in case of vertex decimation,
    //     id of half edge in case of edge collapse
    void update_priority(const std::size_t id);
	
    // Restores the last edge or vertex removal.
    // Currently only works for vertex decimation.
    bool restore_last_simplification_step();

    // Vertex removal specific functions ----------------------------------

    // get the ring of vertices and half edges around a verts
    void get_ring_data(const std::size_t vertex, std::vector<size_t>& ring_verts, std::vector<size_t>& ring_hes);

    // determine whether a vertex can be removed and return the triangles to fill the hole.
    // Input:
    //   vertex: verted id
    // Output:
    //   tri_verts: vertices of the triangles to replace the hole.
    bool analyze_vertex_for_removal(const std::size_t vertex, std::vector<size_t>& tri_verts);

    // Helper function for analyze_vertex_for_removal(). Recursively searches for the best triangulation.
    bool find_best_hole_triangulation(const std::vector<size_t>& ring_verts, std::vector<size_t>& tri_verts, const geo::Plane& avep);

    // Removes a vertex and fills the hole.
    // Input: vertex: the vertex id
    //        tri_verts: id of vertices for triangles to replace the hole.
    void simplify_by_removing_vertex(const std::size_t vertex,const std::vector<size_t>& tri_verts);

    // Removes an edge and replaces it with a vertex
    // Input: he: id of half edge (either would work)
    //        loc: location of vertex to replace half edge ( if set to NULL will use the edge middle point)
    bool simplify_by_collapsing_edge(const std::size_t he, const Eigen::Vector3d* loc = NULL);

    // Simplifies a step based on the type passed to init_simplifaction()
    bool simplify();
	
    

    /*================================================
     * Iterator functions
     *================================================*/
  public:
    /**
     * Initialize an iterator that visits the 1-ring of 'vertex'.
     * \param it The iterator to initialize.
     * \param vertex The index of the vertex to iterate around.
     * \return False if the vertex has no neighbors (ie. a floating vertex)
     */
    bool init_iterator( vvert_iterator& it, std::size_t vertex ) const;
	
    /**
     * move the iterator to a boundary, if one exists, and set its end
     * to the current location
     * \param it The iterator to move
     * \return False if there is no boundary
     */
    bool reset_boundary_iterator( vvert_iterator& it ) const;

    /**
     * \param face index
     * \return true if on boundary false else
     */
    bool isBoundaryFace( std::size_t i ) const;

    /**
     * Advances a vertex iterator to the next vertex in the 1-ring.
     * \param it The iterator to advance
     * \return False if the iterator has completed the loop.
     * It may continue to be used at this point as it will merely restart
     * the loop.
     */
    bool advance_iterator( vvert_iterator& it ) const;
	
    std::size_t deref_iterator( const vvert_iterator& it ) const;
    std::size_t deref_iterator_left_face( const vvert_iterator& it ) const;
    std::size_t deref_iterator_right_face( const vvert_iterator& it ) const;
    std::size_t deref_iterator_left_edge( const vvert_iterator& it ) const;
    std::size_t deref_iterator_right_edge( const vvert_iterator& it ) const;

    Eigen::Vector3d get_vertex( const vvert_iterator& it ) const;
    double get_cotan_weight( const vvert_iterator& it ) const;
    double get_mean_value_weight( const vvert_iterator& it ) const;

    /**
     * Intialize an iterator that visits the faces in the 1-ring of 'vertex'
     */
    bool init_iterator(vface_iterator& it, std::size_t vertex ) const;
    bool advance_iterator( vface_iterator& it ) const;
    std::size_t deref_iterator( const vface_iterator& it ) const;
    Eigen::Vector3d get_normal( const vface_iterator& it ) const;
    Eigen::Vector4d get_plane( const vface_iterator& it ) const;

    bool init_iterator(fvert_iterator& it, std::size_t face ) const;
    bool advance_iterator( fvert_iterator& it ) const;
    std::size_t deref_iterator( const fvert_iterator& it ) const;

    bool init_iterator(fface_iterator& it, std::size_t face ) const;
    bool advance_iterator( fface_iterator& it ) const;
    std::size_t deref_iterator( const fface_iterator& it ) const;

    /*====================================================
     * Helper Functions
     *====================================================*/
  public:
    void getIndicesForFace( size_t tri_index, size_t indicesForFace[3] ) const;
    Eigen::Vector3d getFaceMidpoint(size_t tri_index);
	
    // Extra query functions
    void get_edges_for_face( const size_t face, size_t edges[3] ) const;
    void get_verts_for_edge( const size_t edge, size_t verts[2] ) const;


    /*====================================================
     * Reflection Interface (learn stuff about the mesh)
     *====================================================*/
  public:
    /* get a list of all vertices/face indices in the mesh as floats in
     * contiguous memory. Used for rendering so loss of precision unimportant
     */
    void get_draw_data( float *verts, int *indices ) const;
    void get_draw_normals( float *normals ) const;
    void get_draw_selection( int *selection ) const;
    int  get_edit_count() const;
    void get_face_neighbors(int face_index, size_t neighbors[3]);
    void flag_edited();

    /* selection interface */
    void select_vert( size_t index );
    void deselect_vert( size_t index );
    void deselect_allVerts();
    bool isSelected( size_t index );

    /* get information about internal state */
    std::size_t get_vert_size() const;
    std::size_t get_face_size() const;	

    void test_flip();
    static void test();
    void verify() const;

    void write_to_obj_stream( std::ostream& stream ) const;

  public:
    const HalfEdge& prev( const HalfEdge& cur ) const;
    const HalfEdge& next( const HalfEdge& cur ) const;
    const HalfEdge& twin( const HalfEdge& cur ) const;

    HalfEdge* find_edge( std::size_t vFrom, std::size_t vTo );
    HalfEdge* find_twin( std::size_t vFrom, std::size_t vTo );
    // Semihack for simplification. If this function returns anything
    // bigger than one then you do not have a manifold mesh. Although,
    // your data structure might not have failed yet.
    uint count_edge( std::size_t vFrom, std::size_t vTo ); 

  private:
    std::size_t collapse_edge( std::size_t he ); // TODO REMOVE

    void delete_HalfEdge_impl( std::size_t he );

    template <std::size_t N>
    void delete_HalfEdges_impl( std::size_t (&edges)[N] );
	
    // Splits a boundary edge into 3 by adding 2 vertices on the
    // specified edge, connecting them to the other vertex. Replaces
    // this face with 3 new ones.
    void split_boundary_edge( std::size_t he, std::size_t (*pOutVertIndices)[2] = NULL, std::size_t (*pOutFaceIndices)[3] = NULL );

    /**
     * a simple helper function to determine if adding a particular face
     * will break mesh manifoldness
     * \param input vertices that could become a face
     * \return if this face will break manifoldness
     */
    bool is_safe_addface( std::size_t v1, std::size_t v2, std::size_t v3 );
	
    /**
     * Flip the passed in edge to connect the two opposing vertices
     * Flips both the passed in HalfEdge and its twin
     * \param cur The half edge to flip
     * \return True if successfully flipped, false if edge
     */
    bool flip_edge( HalfEdge &cur );

    /* =========================================
     * Mesh Collision Tests
     * =========================================*/
  public:
    void updateBBox();

    // TODO: update when mesh is changed
    // NOTE: as of right now only used by skeleton animation code
    Eigen::Vector3d bboxMin;
    Eigen::Vector3d bboxMax;
    Eigen::Vector3d bSphereCenter;
    double          bSphereRadius;

    /*****************************************
     * Mesh Variables
     *****************************************/
  private:
    // for mesh read and write.
    friend class MeshIO;
	
    // friend class mesh_adjacency;

    // increment this value every change
    int _edit_count;

    // -------------------------  Vertex Coordinates  -------------------

    // Coordinate of the vertices
    VecOfVerts _verts;

    // -------------------------  HalfEdge data structure ------------

    // All the half-edges that make up the mesh.
    std::vector< HalfEdge > _he_data;
	
    // A mapping from face index to an arbitrary half-edge on its boundary.
    std::vector< std::size_t > _face_to_he;
	
    // A mapping from vertex index to an arbitrary half-edge
    // originating from this vertex. Can be "HOLE_INDEX" for
    // unconnected vertices.
    std::vector< std::size_t > _vert_to_he;


    // -------------------------  Subdivision ------------------------

    // How many times have we refined the mesh.
    // We need to know this in order to refine the boundary correctly
    // for the sqrt3 subdivision.
    uint _n_edge_based_refinement_steps;

    // Mappings for edges.
    bool _is_edge_map_defined;
	
    // _he_to_edge maps each HalfEdge to an edge.
    // _edge_to_he maps each edge to one of its corresponding HalfEdges.
    std::vector< std::size_t > _he_to_edge;
    std::vector< std::size_t > _edge_to_he;

    // Mapping from boundary edges to edges and vice versa.
    // if an edge is not located on the boundary, the corresponding
    // entry in the _edge_to_bedge vector will be HOLE_INDEX.
    std::vector< std::size_t > _edge_to_bedge;
    std::vector< std::size_t > _bedge_to_edge;

    // -------------------------  Simplification  ---------------------

    // Input options
    static constexpr double _min_crease_angle = M_PI / 180 * 85;
    
    // Are we in the middle of simplifying
    // 0 no, 1 vert removal, 2 edge collapse
    uint _is_simplification_in_progress;

    // Current level of simplification
    std::size_t _n_simplification_steps;

    // Is entitiy active?
    std::vector<bool> _is_vert_active;
    std::vector<bool> _is_face_active;
    std::vector<bool> _is_he_active;

    // Entities that were deleted at each step
    // and how many they were.
    std::vector<uint> _n_deleted_verts;
    std::vector<uint> _n_deleted_hes;
    std::vector<uint> _n_deleted_faces;
    std::vector<uint> _n_added_verts;
    std::vector<uint> _n_added_hes;
    std::vector<uint> _n_added_faces;
    std::vector<std::size_t> _id_deleted_verts;
    std::vector<std::size_t> _id_deleted_hes;
    std::vector<std::size_t> _id_deleted_faces;
	
    // Number of active entities that we have to keep track of.
    std::size_t _n_verts_active;
    std::size_t _n_hes_active;
    std::size_t _n_faces_active;

    // Priorities of who to remove -- Vertex Decimation ---------------
    typedef double PriorityVD;
    typedef std::pair<PriorityVD, size_t> PriorityPairVD;
    // Priority queue for vertices.
    typedef std::multimap<PriorityVD, size_t> PriorityContainerVD;
    
    // Priority queue for vertices.
    PriorityContainerVD _prival_vd;
    
    // Each vertex must know which priority object it is associated with.
    // If there is none, then the value would be _prival_vd.end().
    std::vector<PriorityContainerVD::iterator> _priadd_vd;

    // Priorities of who to remove -- Edge Collapse -------------------
    struct PriorityEC
    {
      // Amount of error
      double dv;
      // Location of vertex to replace this edge
      Eigen::Vector3d vbar;
      // The Quartic matrix
      Eigen::Matrix4d Q;
      // Half edges associated with this edge
      size_t he[2];

      // Operator to sort this types to be passed as a comparing
      // functor to the multiset.
      bool operator() (const PriorityEC *p1, const PriorityEC *p2)
      {
	return p1->dv < p2->dv;
      }
    };
    // Priority queue for edges.
    typedef std::multiset<PriorityEC*, PriorityEC> PriorityContainerEC;

    // Priority queue for edges.
    PriorityContainerEC _prival_ec;
    // Each half edge must know which priority object it is associated with.
    // If there is none, then the value would be _prival_ec.end().
    std::vector<PriorityContainerEC::iterator> _priadd_ec;
    // Quartic matrix of each vertex
    std::vector<Eigen::Matrix4d> _vertexQ;
    
    // -------------------------  Selection  ----------------------------

    // if a vertex is selected or not
    std::vector< bool > _is_vert_selected;
	
  };


  inline EditMesh::EditMesh():
    _edit_count(0),
    _n_edge_based_refinement_steps(0),
    _is_edge_map_defined(false),
    _is_simplification_in_progress(0){}

  inline void EditMesh::clear()
  {
    _edit_count = 0;
    _n_edge_based_refinement_steps = 0;

    _he_data.clear();
    _face_to_he.clear();
    _vert_to_he.clear();
    _is_vert_selected.clear();
    _verts.clear();
    _is_edge_map_defined = false;
    _is_simplification_in_progress = 0;
    _he_to_edge.clear();
    _edge_to_he.clear();
  }

  inline void EditMesh::prepare_for_reinit()
  {
    _he_data.resize(0);
    _face_to_he.resize(0);
    _vert_to_he.resize(0);
    _is_vert_selected.resize(0);
    _verts.resize(0);
    _is_edge_map_defined = false;
    _is_simplification_in_progress = 0;
    _he_to_edge.resize(0);
    _edge_to_he.resize(0);
  }

  inline std::size_t EditMesh::add_vertex( double x, double y, double z )
  {
    std::size_t newIndex = _verts.size();
    _verts.emplace_back( x,y,z );
    _vert_to_he.push_back( HOLE_INDEX );
    return newIndex;
  }

  inline Eigen::Vector3d EditMesh::get_vertex( std::size_t vertex ) const
  {
    return _verts[vertex];
  }

  inline void EditMesh::set_vertex( std::size_t i, const Eigen::Vector3d& v )
  {
    _verts[i] = v;
  }

  inline std::size_t EditMesh::add_face( std::size_t v1, std::size_t v2, std::size_t v3 )
  {
    std::size_t v[] = { v1, v2, v3 };
    return this->add_face( v );
  }

  inline bool EditMesh::init_iterator( vvert_iterator& it, std::size_t vertex ) const
  {
    std::size_t vertToHE = _vert_to_he[ vertex ];
			
    if( vertToHE == HOLE_INDEX )
      return false;

    // Store a pointer to the (arbitrary) first half-edge pointing
    // into the specified vertex. This implies
    // _he_data[m_cur->twin].vert == vertex &
    // _he_data[m_cur->next].vert == vertex.  By iterating around the
    // half-edges pointing into 'vertex' we will visit all the
    // vertices in the 1-ring.
    it.m_cur = it.m_end = &_he_data[ _he_data[ vertToHE ].twin ];
    return true;
  }

  inline bool EditMesh::reset_boundary_iterator( vvert_iterator &it ) const
  {
    const HalfEdge *cur = it.m_cur;
    do
      {
        // hooshi: second condition makes sure that we always start
        // from the boundary edge located to the right of the vertex.
        if ( (it.m_cur->face == HOLE_INDEX) &&
             (_he_data[it.m_cur->next].face == HOLE_INDEX) )
	  {
            it.m_end = it.m_cur;
            return true;
	  }
        it.m_cur = &_he_data[ _he_data[ it.m_cur->next ].twin ];
      } while (cur != it.m_cur);

    return false;
  }

  inline bool EditMesh::isBoundaryFace( std::size_t i ) const
  {
    std::size_t he_index = _face_to_he[i];
    const HalfEdge *he = &_he_data[ he_index ];
    while( he->next != he_index ) {
      if( _he_data[he->twin].face == HOLE_INDEX )
	return true;
      he = &_he_data[ he->next ];
    }
    return false;
  }

  inline bool EditMesh::advance_iterator( vvert_iterator& it ) const
  {
    it.m_cur = &_he_data[ _he_data[ it.m_cur->next ].twin ];
    return it.m_cur != it.m_end;
  }

  inline std::size_t EditMesh::deref_iterator( const vvert_iterator& it ) const
  {
    return it.m_cur->vert;
  }

  inline std::size_t EditMesh::deref_iterator_left_face( const vvert_iterator& it ) const
  {
    return _he_data[it.m_cur->twin].face;
  }

  inline std::size_t EditMesh::deref_iterator_right_face( const vvert_iterator& it ) const
  {
    return it.m_cur->face;
  }

  inline std::size_t EditMesh::deref_iterator_left_edge( const vvert_iterator& it ) const
  {
    return it.m_cur->twin;
  }

  inline std::size_t EditMesh::deref_iterator_right_edge( const vvert_iterator& it ) const
  {
    return _he_data[ it.m_cur->twin ].twin;
  }

  inline bool EditMesh::init_iterator( vface_iterator& it, std::size_t vertex ) const
  {
    std::size_t vertToHE = _vert_to_he[ vertex ];
			
    if( vertToHE == HOLE_INDEX )
      return false;

    // Store a pointer to the (arbitrary) first half-edge pointing
    // into the specified vertex. This implies
    // _he_data[m_cur->twin].vert == vertex &
    // _he_data[m_cur->next].vert == vertex.  By iterating around the
    // half-edges pointing into 'vertex' we will visit all the
    // vertices in the 1-ring.
    it.m_cur = it.m_end = &_he_data[ _he_data[ vertToHE ].twin ];
    it.m_next = &_he_data[ _he_data[ it.m_cur->next ].twin ];
    return true;
  }

  inline bool EditMesh::advance_iterator( vface_iterator& it ) const
  {
    it.m_cur = it.m_next;
    it.m_next = &_he_data[ _he_data[ it.m_next->next ].twin ];
    return it.m_cur != it.m_end;
  }

  inline std::size_t EditMesh::deref_iterator( const vface_iterator& it ) const
  {
    return it.m_cur->face;
  }

  inline bool EditMesh::init_iterator(fvert_iterator& it, std::size_t face ) const {
    assert( face < _face_to_he.size() );
    std::size_t faceToHE = _face_to_he[ face ];
    assert( faceToHE < _he_data.size() );

    // Store a pointer to the (arbitrary) first half-edge pointing
    // into the specified vertex. This implies
    // _he_data[m_cur->twin].vert == vertex &
    // _he_data[m_cur->next].vert == vertex.  By iterating around the
    // half-edges pointing into 'vertex' we will visit all the
    // vertices in the 1-ring.
    it.m_cur = it.m_end = &_he_data[ faceToHE ];
    return true;
  }

  inline bool EditMesh::advance_iterator( fvert_iterator& it ) const
  {
    it.m_cur = &_he_data[ it.m_cur->next ];
    return it.m_cur != it.m_end;
  }

  inline std::size_t EditMesh::deref_iterator( const fvert_iterator& it ) const
  {
    return it.m_cur->vert;
  }

  inline bool EditMesh::init_iterator(fface_iterator& it, std::size_t face ) const
  {
    assert( face < _face_to_he.size() );
    std::size_t faceToHE = _face_to_he[ face ];
    assert( faceToHE < _he_data.size() );

    // Store a pointer to the (arbitrary) first half-edge pointing
    // into the specified vertex. This implies
    // _he_data[m_cur->twin].vert == vertex &
    // _he_data[m_cur->next].vert == vertex.  By iterating around the
    // half-edges pointing into 'vertex' we will visit all the
    // vertices in the 1-ring.
    it.m_cur = it.m_end = &_he_data[ faceToHE ];
    return true;
  }

  inline bool EditMesh::advance_iterator( fface_iterator& it ) const
  {
    it.m_cur = &_he_data[ it.m_cur->next ];
    return it.m_cur != it.m_end;
  }

  inline std::size_t EditMesh::deref_iterator( const fface_iterator& it ) const
  {
    return _he_data[ it.m_cur->twin ].face;
  }

  inline const HalfEdge& EditMesh::prev( const HalfEdge& cur ) const
  {
#ifdef USE_PREV
    return _he_data[ cur.prev ];
#else
    return _he_data[ _he_data[ cur.next ].next ];
#endif
  }

  inline const HalfEdge& EditMesh::next( const HalfEdge& cur ) const
  {
    return _he_data[ cur.next ];
  }

  inline const HalfEdge& EditMesh::twin( const HalfEdge& cur ) const
  {
    return _he_data[ cur.twin ];
  }

  inline void EditMesh::select_vert( size_t index )
  {
    _is_vert_selected[index] = true;
    flag_edited();
  }

  inline void EditMesh::deselect_vert( size_t index )
  {
    _is_vert_selected[index] = false;
    flag_edited();
  }

  inline void EditMesh::deselect_allVerts()
  {
    for (size_t i = 0; i < _is_vert_selected.size(); ++i)
      {
        _is_vert_selected[i] = false;
      }
    flag_edited();
  }

  inline bool EditMesh::isSelected( size_t index )
  {
    return _is_vert_selected[index];
  }

  inline std::size_t EditMesh::get_vert_size() const
  {
    if(_is_simplification_in_progress)
      return _n_verts_active;
    else
      return _vert_to_he.size();
  }

  inline std::size_t EditMesh::get_face_size() const
  {
    if(_is_simplification_in_progress)
      return _n_faces_active;
    else
      return _face_to_he.size();
  }

  inline int EditMesh::get_edit_count() const
  {
    return _edit_count;
  }

  inline void EditMesh::get_face_neighbors(int face_index, size_t neighbors[3])
  {
    fface_iterator fit;
    init_iterator(fit, face_index);
    int i = 0;
    do {
      neighbors[i++] = deref_iterator(fit);
    } while(advance_iterator(fit));

    for (; i < 3; i++)
      neighbors[i] = HOLE_INDEX;
  }

  inline void EditMesh::flag_edited()
  {
    ++_edit_count;
  }

} // End of namespace hooshi

#endif /* EDIT_MESH_HXX */
