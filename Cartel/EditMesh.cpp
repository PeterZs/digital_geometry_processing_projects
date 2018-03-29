/* Copyright (c) Darcy Harisson, Russell Gillette
 * April 2014
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated
 * documentation files (the "Software"), to deal in the Software without
 * restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and
 * to permit persons to whom the Software is furnished to do so, subject to the
 * following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
 * NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
 * CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF
 * OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
#define _SCL_SECURE_NO_WARNINGS

#include "EditMesh.h"

#include <Eigen/Geometry>
#include <cmath>

#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <set>

//#if defined(NDEBUG) && defined(ALWAYS_ASSERT)
//#undef NDEBUG
//#endif
//#include <cassert>

EditMesh *loadEditMeshFromFile(std::string file_name);

namespace detail
{
inline void init(half_edge &he, std::size_t next, std::size_t twin,
                 std::size_t vert, std::size_t face)
{
  he.next = next;
  he.twin = twin;
  he.vert = vert;
  he.face = face;
}

void delete_face(std::vector<std::size_t> &faceData,
                 std::vector<half_edge> &heData, std::size_t f)
{
  assert(f < faceData.size());

  // In order to delete the face properly, we need to move a face from the end
  // of the list to overwrite 'f'. Then we need to update the
  // indices stored in the moved face's half-edges.
  faceData[f] = faceData.back();
  faceData.pop_back();

  if (f != faceData.size())
    {
      // std::clog << "Reindexed face " << faceData.size() << " to " << f <<
      // std::endl;

      std::size_t he = faceData[f];
      do
        {
          assert(heData[he].face == faceData.size());

          heData[he].face = f;
          he = heData[he].next;
        }
      while (he != faceData[f]);
    }
}

template <int N>
void delete_faces(std::vector<std::size_t> &faceData,
                  std::vector<half_edge> &heData, std::size_t (&fToDelete)[N])
{
  // Sort the faces by decreasing index so that we can safely delete them all
  // without causing any of them to be accidentally re-indexed (which
  // cause 'fToDelete' to contain invalid indices). This also chooses the
  // optimal deletion order to minimize re-indexing.
  std::sort(fToDelete, fToDelete + N, std::greater<std::size_t>());
  for (std::size_t i = 0; i < N; ++i)
    detail::delete_face(faceData, heData, fToDelete[i]);
}
}

void init_adjacency(std::size_t numVertices,
                    const std::vector<std::size_t> &faces,
                    std::vector<half_edge> &m_heData,
                    std::vector<std::size_t> &m_faceData,
                    std::vector<std::size_t> &m_vertData)
{
  typedef std::map<std::pair<std::size_t, std::size_t>, std::size_t>
      edge_map_type;

  assert(
      faces.size() % 3 == 0 &&
      "Invalid data specified for faces. Must have 3 vertex indices per face.");

  edge_map_type edgeMap; // Use a temporary map to find edge pairs.

  m_heData.reserve(faces.size()); // Assume there are 3 edges per face.
  m_faceData.resize(faces.size() / 3);
  m_vertData.resize(numVertices, HOLE_INDEX); // Init with HOLE_INDEX since a
                                              // vert might be floating w/ no
                                              // faces.

  for (std::size_t i = 0, iEnd = faces.size(); i < iEnd; i += 3)
    {
      std::size_t f[] = {faces[i], faces[i + 1], faces[i + 2]};
      std::size_t fIndex = i / 3;

      // The index of the first (of three) half-edges associated with the
      // current face.
      std::size_t heIndex = m_heData.size();

      half_edge he[3];
      detail::init(he[0], heIndex + 1, HOLE_INDEX, f[0], fIndex);
      detail::init(he[1], heIndex + 2, HOLE_INDEX, f[1], fIndex);
      detail::init(he[2], heIndex, HOLE_INDEX, f[2], fIndex);
#ifdef USE_PREV
      he[0].prev = heIndex + 2;
      he[1].prev = heIndex;
      he[2].prev = heIndex + 1;
#endif

      // These will be set each time a vertex is referenced, but that's fine.
      // The last assignment will stick.
      m_faceData[fIndex] = heIndex;
      m_vertData[f[0]] = heIndex;
      m_vertData[f[1]] = heIndex + 1;
      m_vertData[f[2]] = heIndex + 2;

      edge_map_type::iterator it;

      it = edgeMap.lower_bound(std::make_pair(f[0], f[1]));
      if (it != edgeMap.end() && it->first.first == f[0] &&
          it->first.second == f[1])
        {
          m_heData[it->second].twin = heIndex;
          he[0].twin = it->second;
          edgeMap.erase(it);
        }
      else
        {
          he[0].twin = HOLE_INDEX;
          edgeMap.insert(it, std::make_pair(std::make_pair(f[1], f[0]),
                                            heIndex)); // NOTE: Reversed order
                                                       // since we are matching
                                                       // opposite half_edge.
        }

      it = edgeMap.lower_bound(std::make_pair(f[1], f[2]));
      if (it != edgeMap.end() && it->first.first == f[1] &&
          it->first.second == f[2])
        {
          m_heData[it->second].twin = heIndex + 1;
          he[1].twin = it->second;
          edgeMap.erase(it);
        }
      else
        {
          he[1].twin = HOLE_INDEX;
          edgeMap.insert(it, std::make_pair(std::make_pair(f[2], f[1]),
                                            heIndex + 1)); // NOTE: Reversed
                                                           // order since we are
                                                           // matching opposite
                                                           // half_edge.
        }

      it = edgeMap.lower_bound(std::make_pair(f[2], f[0]));
      if (it != edgeMap.end() && it->first.first == f[2] &&
          it->first.second == f[0])
        {
          m_heData[it->second].twin = heIndex + 2;
          he[2].twin = it->second;
          edgeMap.erase(it);
        }
      else
        {
          he[2].twin = HOLE_INDEX;
          edgeMap.insert(it, std::make_pair(std::make_pair(f[0], f[2]),
                                            heIndex + 2)); // NOTE: Reversed
                                                           // order since we are
                                                           // matching opposite
                                                           // half_edge.
        }

      m_heData.push_back(he[0]);
      m_heData.push_back(he[1]);
      m_heData.push_back(he[2]);
    }

  // Keep track of the last edge we processed so we can hook up half_edge::prev
  // as we go.
  std::size_t prev = HOLE_INDEX;

  // Add half-edges for any holes. Any edges still in the map are holes.
  edge_map_type::iterator it = edgeMap.begin();
  while (it != edgeMap.end())
    {
      half_edge he;
      detail::init(he, HOLE_INDEX, it->second, it->first.first, HOLE_INDEX);
#ifdef USE_PREV
      he.prev = prev;
      prev = m_heData.size(); // Size is the index of the half_edge we are about
                              // to push into the list.
#endif

      m_heData[he.twin].twin = m_heData.size();
      m_heData.push_back(he);

      std::size_t curVert = it->first.first;
      std::size_t nextVert = it->first.second; // We are about to erase this
                                               // information, so store it to
                                               // use later.

      edgeMap.erase(it); // We are done with this edge now.

      half_edge *twinPrev = &m_heData[m_heData[m_heData[he.twin].next].next];
      while (twinPrev->twin != HOLE_INDEX &&
             m_heData[twinPrev->twin].face != HOLE_INDEX)
        {
          assert(m_heData[twinPrev->next].vert == nextVert);
          assert(m_heData[twinPrev->twin].vert == nextVert);
          twinPrev = &m_heData[m_heData[m_heData[twinPrev->twin].next].next];
        }

      if (twinPrev->twin == HOLE_INDEX)
        {
          // We haven't processed the next edge in the loop yet. Let's do so now
          // so we can assume the index of the next half-edge.
          m_heData.back().next = m_heData.size();
          it = edgeMap.find(std::make_pair(nextVert, twinPrev->vert));

          assert(it != edgeMap.end());
        }
      else
        {
          assert(m_heData[twinPrev->twin].vert == nextVert);
          assert(m_heData[twinPrev->twin].face == HOLE_INDEX);

          // We already processed this edge and have a valid index for the next
          // half_edge.
          m_heData.back().next = twinPrev->twin;
#ifdef USE_PREV
          m_heData[twinPrev->twin].prev = prev; // Complete the loop
          prev = HOLE_INDEX;
#endif
          it = edgeMap.begin(); // Arbitrarily pick the next edge in the list.
        }
    }

  assert(edgeMap.empty());
}

void EditMesh::init(const std::vector<double> &xyzPositions,
                    const std::vector<std::size_t> &triangleVerts)
{
  assert(xyzPositions.size() % 3 == 0 && "Invalid vertex positions for "
                                         "EditMesh::init(). Must have 3 values "
                                         "per-vertex.");
  assert(triangleVerts.size() % 3 == 0 && "Invalid face data for "
                                          "EditMesh::init(). Must have 3 "
                                          "vertex indices per face.");

  // m_vertices.resize( Eigen::NoChange, xyzPositions.size() / 3 );
  m_vertices.resize(xyzPositions.size() / 3);
  m_selected.resize(xyzPositions.size() / 3);

  // The Eigen matrix has the same format as the incoming vector so we can
  // straight copy it.
  // HACK: This is pretty sketchy and relies on Eigen::Vector3d having the same
  // layout as a double[3] and nothing extra or fancy alignment.
  std::copy(xyzPositions.begin(), xyzPositions.end(),
            m_vertices.front().data());

  init_adjacency(xyzPositions.size() / 3, triangleVerts, m_heData, m_faceData,
                 m_vertData);
}

half_edge *EditMesh::find_twin(std::size_t vFrom, std::size_t vTo)
{
  vvert_iterator it;
  if (!this->init_iterator(it, vFrom))
    return NULL;

  do
    {
      if (this->deref_iterator(it) == vTo)
        return const_cast<half_edge *>(
            it.m_cur); // Gross. This is just laziness.
    }
  while (this->advance_iterator(it));

  return NULL;
}

half_edge *EditMesh::find_edge(std::size_t vFrom, std::size_t vTo)
{
  if (half_edge *he = this->find_twin(vFrom, vTo))
    return &m_heData[he->twin];
  return NULL;
}

bool EditMesh::flip_edge(half_edge &he)
{
  half_edge &twin = m_heData[he.twin];

  if (he.face == HOLE_INDEX || twin.face == HOLE_INDEX)
    return false;

  std::size_t he_tri[3];
  std::size_t twin_tri[3];

  // prep: gather half edge indices in
  // the order they should be after flip
  he_tri[0] = he.next;
  twin_tri[0] = twin.next;
  he_tri[1] = twin.twin;
  twin_tri[1] = he.twin;
  he_tri[2] = m_heData[twin_tri[0]].next;
  twin_tri[2] = m_heData[he_tri[0]].next;

  if (m_heData[he_tri[2]].vert == m_heData[twin_tri[2]].vert)
    return false;

  // step 1: ensure he's verts don't point to
  // either half_edge (does not break mesh)
  m_vertData[he.vert] = twin_tri[0];
  m_vertData[twin.vert] = he_tri[0];

  // step 2: set the he's vert to new originating vert
  he.vert = m_heData[twin_tri[2]].vert;
  twin.vert = m_heData[he_tri[2]].vert;

  // step 3: ensure the faces point to one
  // of the half edges connected to them
  m_faceData[he.face] = he_tri[0];
  m_faceData[twin.face] = twin_tri[0];

  // step 4: fix two edges that will point
  // to the wrong face
  m_heData[he_tri[2]].face = he.face;
  m_heData[twin_tri[2]].face = twin.face;

  // step 5: ensure half edges point to
  // each other
  for (int i = 0; i < 3; ++i)
    {
      m_heData[he_tri[i]].next = he_tri[(i + 1) % 3];
      m_heData[twin_tri[i]].next = twin_tri[(i + 1) % 3];
    }

  return true;
}

// IMPORTANT: Given a collection of half-edges to delete (ex. When removing a
// face we need to kill 2, 4, or 6 half-edges) they must be deleting in
// decreasing index order!
void EditMesh::delete_half_edge_impl(std::size_t he)
{
  assert((m_heData[he].vert >= m_vertData.size() ||
          m_vertData[m_heData[he].vert] != he) &&
         "Deleting this half_edge leaves a dangling link from a vertex. Must "
         "handle this first");

  // Move a half_edge from the end overtop of the half_edge we are deleting,
  // then update the indices of linked half_edges.
  m_heData[he] = m_heData.back();
  m_heData.pop_back();

  // We may have just deleted the item at the end of the list, so we have
  // nothing to update since the indices didn't change.
  if (he != m_heData.size())
    {
      const half_edge &heMoved = m_heData[he];

      // If the moved half_edge was the arbitrary half_edge linked to the
      // vertex, update it.
      if (m_vertData[heMoved.vert] == m_heData.size())
        m_vertData[heMoved.vert] = he;

      // If the moved half_edge was the arbitrary half_edge linked to the face,
      // update it.
      if (heMoved.face != HOLE_INDEX &&
          m_faceData[heMoved.face] == m_heData.size())
        m_faceData[heMoved.face] = he;

      assert(heMoved.twin < m_heData.size());
      assert(m_heData[heMoved.twin].twin == m_heData.size());
      m_heData[heMoved.twin].twin = he;

// NOTE: If we are deleting a bundle of half_edges, then by definition we must
// call delete_half_edge() in decreasing order of indices. That prevents
//       me from having to worry about moving a partially destroyed half_edge
//       into the 'he' position.

#ifdef USE_PREV
      assert(m_heData[heMoved.prev].next == m_heData.size());
      m_heData[heMoved.prev].next = he;

      assert(m_heData[heMoved.next].prev == m_heData.size());
      m_heData[heMoved.next].prev = he;
#else
      // Have to loop around the face until we find the half_edge using
      // 'heMoved' as its 'next' entry, then update it.
      std::size_t hePrev = heMoved.next;
      while (m_heData[hePrev].next != m_heData.size())
        hePrev = m_heData[hePrev].next;

      assert(m_heData[hePrev].next == m_heData.size());
      m_heData[hePrev].next = he;
#endif
    }
}

template <std::size_t N>
void EditMesh::delete_half_edges_impl(std::size_t (&heToDelete)[N])
{
  std::sort(heToDelete, heToDelete + N, std::greater<std::size_t>());
  for (std::size_t i = 0; i < N; ++i)
    this->delete_half_edge_impl(heToDelete[i]);
}

bool g_debug = false;

std::size_t EditMesh::collapse_edge(std::size_t he)
{
  assert(he < m_heData.size());
  assert(m_heData[he].face != HOLE_INDEX &&
         m_heData[m_heData[he].twin].face != HOLE_INDEX &&
         "Cannot collapse a boundary edge");

  const half_edge &heBase = m_heData[he];
  const half_edge &heTwin = m_heData[heBase.twin];

  // We are going to delete the faces on either side of the chosen edge,
  // so we need to delete 3 half_edges and patch up the twin links on the 4
  // bordering edges.
  std::size_t heBorder[4];
  heBorder[0] = m_heData[heBase.next].twin;
  heBorder[1] = m_heData[m_heData[heBase.next].next].twin;
  heBorder[2] = m_heData[m_heData[heTwin.next].next].twin;
  heBorder[3] = m_heData[heTwin.next].twin;

  // TODO: Relax this assertion. We should be able to collapse a spike jutting
  // into a hole.
  assert((m_heData[heBorder[0]].face != HOLE_INDEX ||
          m_heData[heBorder[1]].face != HOLE_INDEX) &&
         "Cannot collapse an edge on a face with holes on either side.");
  assert((m_heData[heBorder[2]].face != HOLE_INDEX ||
          m_heData[heBorder[3]].face != HOLE_INDEX) &&
         "Cannot collapse an edge on a face with holes on either side.");

  // Check if we can actually collapse. This checks for a degree 3 vertex at the
  // vertices not on the edge we are collapsing.
  if (m_heData[m_heData[m_heData[heBorder[1]].next].twin].next == heBorder[0])
    return HOLE_INDEX;
  if (m_heData[m_heData[m_heData[heBorder[2]].next].twin].next == heBorder[3])
    return HOLE_INDEX;

  // Capture the indices of things (2 faces & 6 half-edges) we want to delete.
  std::size_t fToDelete[] = {heBase.face, heTwin.face};
  std::size_t heToDelete[] = {
      he,          heBase.next, m_heData[heBase.next].next,
      heBase.twin, heTwin.next, m_heData[heTwin.next].next};

#ifndef NDEBUG
  // We can't be deleting border edges!
  for (auto i : heToDelete)
    {
      if (std::find(heBorder, heBorder + 4, i) != heBorder + 4)
        return HOLE_INDEX;
      // assert( std::find( heBorder, heBorder + 4, i ) == heBorder + 4 );
    }

  if (g_debug)
    {
      std::vector<std::set<std::size_t>> verts(3);

      verts[0].insert(heBase.vert);
      verts[0].insert(heTwin.vert);

      for (size_t i = 1; i < verts.size(); ++i)
        {
          for (auto v : verts[i - 1])
            {
              vvert_iterator it;
              this->init_iterator(it, v);
              do
                {
                  verts[i].insert(this->deref_iterator(it));
                }
              while (this->advance_iterator(it));
            }
        }

      std::vector<std::size_t> orderedVerts(verts.back().begin(),
                                            verts.back().end());
      std::set<std::size_t> faces;

      std::vector<double> vpos;
      std::vector<std::size_t> finds;

      for (auto v : orderedVerts)
        {
          vpos.push_back(m_vertices[v].x());
          vpos.push_back(m_vertices[v].y());
          vpos.push_back(m_vertices[v].z());
          // std::clog << "m.add_vert( " << m_vertices[v].x() << ", " <<
          // m_vertices[v].y() << ", " << m_vertices[v].z() << " );" <<
          // std::endl;
        }

      // Visit the 1-ring
      for (auto v : verts[1])
        {
          vface_iterator it;
          this->init_iterator(it, v);
          do
            {
              if (this->deref_iterator(it) != HOLE_INDEX &&
                  faces.find(this->deref_iterator(it)) == faces.end())
                {
                  faces.insert(this->deref_iterator(it));

                  fvert_iterator itFace;
                  this->init_iterator(itFace, this->deref_iterator(it));

                  std::size_t f[3];
                  std::size_t i = 0;
                  do
                    {
                      f[i++] =
                          std::find(orderedVerts.begin(), orderedVerts.end(),
                                    this->deref_iterator(itFace)) -
                          orderedVerts.begin();
                    }
                  while (this->advance_iterator(itFace));

                  finds.push_back(f[0]);
                  finds.push_back(f[1]);
                  finds.push_back(f[2]);
                  // std::clog << "m.add_face( " << f[0] << ", " << f[1] << ", "
                  // << f[2] << " );" << std::endl;
                }
            }
          while (this->advance_iterator(it));
        }

      std::size_t base =
          std::find(orderedVerts.begin(), orderedVerts.end(), heBase.vert) -
          orderedVerts.begin();
      std::size_t twin =
          std::find(orderedVerts.begin(), orderedVerts.end(), heTwin.vert) -
          orderedVerts.begin();
      std::clog << "m.collapse_edge( " << base << ", " << twin << " );"
                << std::endl;

      EditMesh m;
      m.init(vpos, finds);
      std::ofstream fout("debug.obj");
      m.write_to_obj_stream(fout);
      fout.close();
    }
#endif

  // We may also need to fix the vertex->half_edge link for the verts using
  // these faces. There are technically 4, but we only update the 3 that are not
  // going to be deleted.
  std::size_t verts[] = {this->prev(heBase).vert, heBase.vert,
                         this->prev(heTwin).vert};

  // Move the base vertex (arbitrarily) to the middle of the edge. Could leave
  // it where it is, or do something fancier too.
  m_vertices[heBase.vert] =
      0.5 * (m_vertices[heBase.vert] + m_vertices[heTwin.vert]);

  // Adjust all the twin's 1-ring to link to the vertex we are not going to
  // delete.
  std::size_t heIt = this->twin(this->next(heBase)).next;
  std::size_t heEnd = heBase.twin;
  for (; heIt != heEnd; heIt = this->twin(m_heData[heIt]).next)
    {
      assert(m_heData[heIt].vert == heTwin.vert);

      // Associate to the other vertex now, so we can delete this one.
      m_heData[heIt].vert = heBase.vert;
    }

  // Fix the vert associations if required, picking a non-hole face.
  if (m_vertData[verts[0]] == m_heData[heBorder[1]].twin)
    m_vertData[verts[0]] = (m_heData[heBorder[0]].face != HOLE_INDEX)
                               ? heBorder[0]
                               : m_heData[heBorder[1]].next;
  if (m_vertData[verts[1]] == he || m_vertData[verts[1]] == heTwin.next)
    m_vertData[verts[1]] =
        (m_heData[heBorder[1]].face != HOLE_INDEX) ? heBorder[1] : heBorder[2];
  if (m_vertData[verts[2]] == m_heData[heBorder[2]].twin)
    m_vertData[verts[2]] = (m_heData[heBorder[3]].face != HOLE_INDEX)
                               ? heBorder[3]
                               : m_heData[heBorder[2]].next;

  // "Delete" the other vertex
  m_vertData[heTwin.vert] = HOLE_INDEX;

  // Collapse the two triangles bordering our chosen half-edge by connecting the
  // opposite edges together.
  m_heData[heBorder[0]].twin = heBorder[1];
  m_heData[heBorder[1]].twin = heBorder[0];
  m_heData[heBorder[2]].twin = heBorder[3];
  m_heData[heBorder[3]].twin = heBorder[2];

  // Have to delete the faces in the proper order.
  if (fToDelete[0] < fToDelete[1])
    std::swap(fToDelete[0], fToDelete[1]);

  this->delete_half_edges_impl(heToDelete);
  detail::delete_face(m_faceData, m_heData, fToDelete[0]);
  detail::delete_face(m_faceData, m_heData, fToDelete[1]);

  return verts[1];
}

std::size_t EditMesh::add_face(std::size_t (&v)[3])
{
  std::size_t faceIndex = m_faceData.size();
  std::size_t heIndex = m_heData.size();

  // Find the half-edges on the hole face we are filling. We must either:
  //  1. Find no half-edges, if all vertices are unconnected from the mesh.
  //  3. Find one half-edge, if one of the vertices is not connected to the
  //  existing mesh.
  //  4. Find two half-edges, if we are adding a triangle inside of a polygonal
  //  hole.
  //  2. Find three half-edges, if we are filling an existing triangular hole.
  half_edge *he[] = {this->find_edge(v[0], v[1]), this->find_edge(v[1], v[2]),
                     this->find_edge(v[2], v[0])};

  // Find the first half-edge we need to modify. This is an edge
  std::size_t base = HOLE_INDEX;
  for (std::size_t i = 0; i < 3 && base == HOLE_INDEX; ++i)
    {
      if (he[i])
        {
          assert(he[i]->face == HOLE_INDEX && "Non-manifold mesh detected. "
                                              "Cannot connect to an edge which "
                                              "already has two incident faces "
                                              "(ie. One side must be a hole)");
          if (!he[(i + 2) % 3])
            base = i;
        }
    }

  if (base == HOLE_INDEX)
    {
      // This triangle is not connected to any others, or we completely filled a
      // triangular hole.
      if (he[0] /*|| he[1] || he[2]*/)
        {
          assert(he[0] && he[1] && he[2]);
          assert(he[0]->face == HOLE_INDEX && he[1]->face == HOLE_INDEX &&
                 he[2]->face == HOLE_INDEX);
          assert(&m_heData[he[0]->next] == he[1] &&
                 &m_heData[he[1]->next] == he[2] &&
                 &m_heData[he[2]->next] == he[0]);

          // Update the face index of the triangular hole to convert it to a
          // face.
          he[0]->face = he[1]->face = he[2]->face = faceIndex;
          m_faceData.push_back(he[2]->next);
        }
      else
        {
          assert(!he[0] && !he[1] && !he[2]);
          assert(m_vertData[v[0]] == HOLE_INDEX &&
                 m_vertData[v[1]] == HOLE_INDEX &&
                 m_vertData[v[2]] == HOLE_INDEX &&
                 "Non-manifold mesh detected. Cannot have two hole faces "
                 "incident on a vertex.");

          // Make 3 new half-edges for the triangle, and 3 new half-edges for
          // the hole outside of the triangle.
          half_edge newHe[6];
          detail::init(newHe[0], heIndex + 1, heIndex + 5, v[0], faceIndex);
          detail::init(newHe[1], heIndex + 2, heIndex + 4, v[1], faceIndex);
          detail::init(newHe[2], heIndex, heIndex + 3, v[2], faceIndex);
          detail::init(newHe[3], heIndex + 4, heIndex + 2, v[0], HOLE_INDEX);
          detail::init(newHe[4], heIndex + 5, heIndex + 1, v[2], HOLE_INDEX);
          detail::init(newHe[5], heIndex + 3, heIndex, v[1], HOLE_INDEX);
#ifdef USE_PREV
          newHe[0].prev = heIndex + 2;
          newHe[1].prev = heIndex;
          newHe[2].prev = heIndex + 1;

          newHe[3].prev = heIndex + 5;
          newHe[4].prev = heIndex + 3;
          newHe[5].prev = heIndex + 4;
#endif

          m_vertData[v[0]] = heIndex;
          m_vertData[v[1]] = heIndex + 1;
          m_vertData[v[2]] = heIndex + 2;

          m_faceData.push_back(heIndex);
          m_heData.push_back(newHe[0]);
          m_heData.push_back(newHe[1]);
          m_heData.push_back(newHe[2]);
          m_heData.push_back(newHe[3]);
          m_heData.push_back(newHe[4]);
          m_heData.push_back(newHe[5]);
        }
    }
  else
    {
      std::size_t next = (base + 1) % 3, prev = (base + 2) % 3;
      std::size_t baseIndex =
          static_cast<std::size_t>(he[base] - &m_heData.front());

      assert(!he[prev]);

      if (he[next])
        {
          // We have two edges to steal from the hole, and we need to add two
          // new half-edges
          half_edge newHe[2];
          detail::init(newHe[0], baseIndex, heIndex + 1, v[prev], faceIndex);
          detail::init(newHe[1], he[next]->next, heIndex, v[base], HOLE_INDEX);

#ifdef USE_PREV
          newHe[0].prev = he[base]->next;
          newHe[1].prev = he[base]->prev;

          m_heData[he[base]->prev].next = heIndex + 1;
          m_heData[he[next]->next].prev = heIndex + 1;

          he[next]->next = heIndex;
          he[base]->prev = heIndex;
#else
          // Have to find the previous half_edge in the polygonal hole so we can
          // point it to the new half-edge in the hole.
          half_edge *hePrev = &m_heData[he[next]->next];
          while (&m_heData[hePrev->next] != he[base])
            {
              hePrev = &m_heData[hePrev->next];
              assert(hePrev != he[next]); // To catch weirdness.
            }
          assert(&m_heData[hePrev->next] == he[base]);

          hePrev->next = heIndex + 1;
          he[next]->next = heIndex;
#endif

          // Update the face indices of the half-edges to indicate they are in a
          // triangle now.
          he[base]->face = he[next]->face = faceIndex;

          m_faceData.push_back(heIndex);
          m_heData.push_back(newHe[0]);
          m_heData.push_back(newHe[1]);
        }
      else
        {
          assert(m_vertData[v[prev]] == HOLE_INDEX &&
                 "Non-manifold mesh detected. Cannot have two hole faces "
                 "incident on a vertex.");

          // We have one edge to steal from the hole, and we need to add four
          // new half-edges.
          half_edge newHe[4];
          detail::init(newHe[0], baseIndex, heIndex + 2, v[prev], faceIndex);
          detail::init(newHe[1], heIndex, heIndex + 3, v[next], faceIndex);
          detail::init(newHe[2], heIndex + 3, heIndex, v[base], HOLE_INDEX);
          detail::init(newHe[3], he[base]->next, heIndex + 1, v[prev],
                       HOLE_INDEX);

#ifdef USE_PREV
          newHe[0].prev = heIndex + 1;
          newHe[1].prev = baseIndex;
          newHe[2].prev = he[base]->prev;
          newHe[3].prev = heIndex + 2;

          m_heData[he[base]->prev].next = heIndex + 2;
          m_heData[he[base]->next].prev = heIndex + 3;

          he[base]->prev = heIndex;
          he[base]->next = heIndex + 1;
#else
          // Have to find the previous half_edge in the polyognal hole so we can
          // point it to the new half-edge in the hole.
          half_edge *hePrev = &m_heData[he[base]->next];
          while (&m_heData[hePrev->next] != he[base])
            {
              hePrev = &m_heData[hePrev->next];
              assert(hePrev != he[next]); // To catch weirdness.
            }
          assert(&m_heData[hePrev->next] == he[base]);

          hePrev->next = heIndex + 2;
          he[base]->next = heIndex + 1;
#endif

          // Update the face indices of the half-edges to indicate they are in a
          // triangle now.
          he[base]->face = faceIndex;

          m_vertData[v[prev]] = heIndex;
          m_faceData.push_back(heIndex);
          m_heData.push_back(newHe[0]);
          m_heData.push_back(newHe[1]);
          m_heData.push_back(newHe[2]);
          m_heData.push_back(newHe[3]);
        }
    }

  return faceIndex;
}

void EditMesh::delete_face(std::size_t f)
{
  assert(f < m_faceData.size());

  // We can assume that this face has 3 half-edges.
  std::size_t heIndices[3];
  heIndices[0] = m_faceData[f];

  half_edge *he[3];
  he[0] = &m_heData[heIndices[0]];
  he[1] = &m_heData[he[0]->next];
  he[2] = &m_heData[he[1]->next];

  heIndices[1] = he[0]->next;
  heIndices[2] = he[1]->next;

  assert(he[0]->face == f && he[1]->face == f && he[2]->face == f);
  assert(he[2]->next == m_faceData[f]);

  // Search for an edge that has a neighbor, but its prev edge doesn't. This is
  // a canonical place to construct the algorithm from.
  std::size_t base = HOLE_INDEX;
  for (std::size_t i = 0; i < 3 && base == HOLE_INDEX; ++i)
    {
      if (m_heData[he[i]->twin].face != HOLE_INDEX &&
          m_heData[he[(i + 2) % 3]->twin].face == HOLE_INDEX)
        base = i;
    }

  if (base == HOLE_INDEX)
    {
      if (m_heData[he[0]->twin].face == HOLE_INDEX)
        {
          // This is a lone triangle, so delete its half-edges and the exterior
          // hole surrounding it too.

          // TODO: Remove the floating vertices? Currently we are leaving them.
          m_vertData[he[0]->vert] = HOLE_INDEX;
          m_vertData[he[1]->vert] = HOLE_INDEX;
          m_vertData[he[2]->vert] = HOLE_INDEX;

          // Delete all of the edges (both inside & outside half-edges). Must do
          // this last since indices can change arbitrarily when deleting.
          std::size_t toDelete[] = {heIndices[0], heIndices[1], heIndices[2],
                                    he[0]->twin,  he[1]->twin,  he[2]->twin};

          this->delete_half_edges_impl(toDelete);
          detail::delete_face(m_faceData, m_heData, f);
        }
      else
        {
          // This is an interior triangle. Only have to change the face_index to
          // HOLE_INDEX for these edges.

          // Adjust any vertex references to new edges in non-hole faces.
          if (m_vertData[he[0]->vert] == heIndices[0])
            m_vertData[he[0]->vert] = he[2]->twin;
          if (m_vertData[he[1]->vert] == heIndices[1])
            m_vertData[he[1]->vert] = he[0]->twin;
          if (m_vertData[he[2]->vert] == heIndices[2])
            m_vertData[he[2]->vert] = he[1]->twin;

          // Flag all these half-edges as being a hole now.
          he[0]->face = he[1]->face = he[2]->face = HOLE_INDEX;
          detail::delete_face(m_faceData, m_heData, f);
        }
    }
  else
    {
      std::rotate(he, he + base, he + 3);
      std::rotate(heIndices, heIndices + base, heIndices + 3);
      assert(m_heData[he[0]->twin].face != HOLE_INDEX);
      assert(m_heData[he[2]->twin].face == HOLE_INDEX);

      if (m_heData[he[1]->twin].face != HOLE_INDEX)
        {
// We have one edge to remove, and a hole to connect to.
#ifdef USE_PREV
          he[1]->next = m_heData[he[2]->twin].next;
          he[0]->prev = m_heData[he[2]->twin].prev;
          m_heData[he[1]->next].prev = heIndices[1];
          m_heData[he[0]->prev].next = heIndices[0];
#else
          he[1]->next = m_heData[he[2]->twin].next;

          std::size_t hePrev = he[1]->next;
          while (m_heData[hePrev].next != he[2]->twin)
            hePrev = m_heData[hePrev].next;

          assert(m_heData[hePrev].next == he[2]->twin);
          m_heData[hePrev].next = heIndices[0];
#endif

          assert(m_heData[m_vertData[he[0]->vert]].face != HOLE_INDEX);
          assert(m_heData[m_vertData[he[1]->vert]].face != HOLE_INDEX);
          assert(m_heData[m_vertData[he[2]->vert]].face != HOLE_INDEX);

          // We may need to update the vertices if they referenced the edges we
          // are deleting. Choose new half-edges that are inside
          // non-hole triangles.
          if (m_vertData[he[0]->vert] == heIndices[0])
            m_vertData[he[0]->vert] = m_heData[he[0]->twin].next;
          if (m_vertData[he[1]->vert] == heIndices[1])
            m_vertData[he[1]->vert] = he[0]->twin;
          if (m_vertData[he[2]->vert] == heIndices[2])
            m_vertData[he[2]->vert] = he[1]->twin;

          assert(m_heData[m_vertData[he[0]->vert]].face != HOLE_INDEX);
          assert(m_heData[m_vertData[he[1]->vert]].face != HOLE_INDEX);
          assert(m_heData[m_vertData[he[2]->vert]].face != HOLE_INDEX);

          he[0]->face = he[1]->face = HOLE_INDEX;

          std::size_t toDelete[] = {heIndices[2], he[2]->twin};

          // Delete the edges and face. Must do this last since indices can
          // change arbitrarily when deleting.
          this->delete_half_edges_impl(toDelete);
          detail::delete_face(m_faceData, m_heData, f);
        }
      else
        {
// We have two edges to remove, a vertex that will become floating, and a hole
// to connect to.
#ifdef USE_PREV
          he[0]->next = m_heData[he[1]->twin].next;
          he[0]->prev = m_heData[he[2]->twin].prev;
          m_heData[he[0]->next].prev = heIndices[0];
          m_heData[he[0]->prev].next = heIndices[0];
#else
          he[0]->next = m_heData[he[1]->twin].next;

          std::size_t hePrev = he[0]->next;
          while (m_heData[hePrev].next != he[2]->twin)
            hePrev = m_heData[hePrev].next;

          assert(m_heData[hePrev].next == he[2]->twin);
          m_heData[hePrev].next = heIndices[0];
#endif

          // We may need to update the vertices if they referenced the edges we
          // are deleting. Choose new half-edges that are inside
          // non-hole triangles.
          if (m_vertData[he[1]->vert] == heIndices[1])
            m_vertData[he[1]->vert] = he[0]->twin;
          if (m_vertData[he[0]->vert] == he[2]->twin ||
              m_vertData[he[0]->vert] == heIndices[0])
            m_vertData[he[0]->vert] = m_heData[he[0]->twin].next;
          m_vertData[he[2]->vert] = HOLE_INDEX;

          // Update the face association of the one half-edge we are keeping (it
          // joins the hole).
          he[0]->face = HOLE_INDEX;

          // Delete the edges and face. Must do this last since indices can
          // change arbitrarily when deleting.
          std::size_t toDelete[] = {heIndices[1], heIndices[2], he[1]->twin,
                                    he[2]->twin};

          this->delete_half_edges_impl(toDelete);
          detail::delete_face(m_faceData, m_heData, f);
        }
    }
}

std::size_t EditMesh::split_face_center(std::size_t f,
                                        std::size_t (*pOutFaceIndices)[3])
{
  assert(f < m_faceData.size());
  assert(m_faceData[f] < m_heData.size() && m_heData[m_faceData[f]].face == f);

  std::size_t he[3];
  he[0] = m_faceData[f];
  he[1] = m_heData[he[0]].next;
  he[2] = m_heData[he[1]].next;

  assert(m_heData[he[2]].next == he[0]);
  assert(m_heData[he[0]].vert < m_vertices.size() &&
         m_heData[he[1]].vert < m_vertices.size() &&
         m_heData[he[2]].vert < m_vertices.size());

  // New vert at face center
  Eigen::Vector3d newVert =
      (m_vertices[m_heData[he[0]].vert] + m_vertices[m_heData[he[1]].vert] +
       m_vertices[m_heData[he[2]].vert]) /
      3.0;

  std::size_t newVertIndex = m_vertices.size();
  m_vertices.push_back(newVert);

  // Each half-edge gets associated to a new face, and we add 6 half-edges from
  // the old vertices to the new.
  std::size_t newHeIndex = m_heData.size();
  std::size_t newFaceIndex = m_faceData.size();

  if (pOutFaceIndices)
    {
      (*pOutFaceIndices)[0] = f;
      (*pOutFaceIndices)[1] = newFaceIndex;
      (*pOutFaceIndices)[2] = newFaceIndex + 1;
    }

  // Create six new half-edges connecting the center vertex to the old triangle
  // corners.
  half_edge newHe[6];
  detail::init(newHe[0], newHeIndex + 1, newHeIndex + 3, m_heData[he[1]].vert,
               f);
  detail::init(newHe[1], he[0], newHeIndex + 4, newVertIndex, f);
  detail::init(newHe[2], newHeIndex + 3, newHeIndex + 5, m_heData[he[2]].vert,
               newFaceIndex);
  detail::init(newHe[3], he[1], newHeIndex, newVertIndex, newFaceIndex);
  detail::init(newHe[4], newHeIndex + 5, newHeIndex + 1, m_heData[he[0]].vert,
               newFaceIndex + 1);
  detail::init(newHe[5], he[2], newHeIndex + 2, newVertIndex, newFaceIndex + 1);

  // Connect the old half-edges to the new ones, and update their face
  // association.
  // m_heData[he[0]].face = f;
  m_heData[he[0]].next = newHeIndex;
  m_heData[he[1]].face = newFaceIndex;
  m_heData[he[1]].next = newHeIndex + 2;
  m_heData[he[2]].face = newFaceIndex + 1;
  m_heData[he[2]].next = newHeIndex + 4;

#ifdef USE_PREV
  newHe[0].prev = he[0];
  newHe[1].prev = newHeIndex;
  newHe[2].prev = he[1];
  newHe[3].prev = newHeIndex + 2;
  newHe[4].prev = he[2];
  newHe[5].prev = newHeIndex + 4;

  m_heData[he[0]].prev = newHeIndex + 1;
  m_heData[he[1]].prev = newHeIndex + 3;
  m_heData[he[2]].prev = newHeIndex + 5;
#endif

  m_vertData.push_back(newHeIndex + 3); // Arbitrary from 1, 3 & 5
  m_faceData[f] = he[0];
  m_faceData.push_back(he[1]);
  m_faceData.push_back(he[2]);
  m_heData.push_back(newHe[0]);
  m_heData.push_back(newHe[1]);
  m_heData.push_back(newHe[2]);
  m_heData.push_back(newHe[3]);
  m_heData.push_back(newHe[4]);
  m_heData.push_back(newHe[5]);

  return newVertIndex;
}

void EditMesh::split_boundary_edge(std::size_t heToSplit,
                                   std::size_t (*pOutVertIndices)[2],
                                   std::size_t (*pOutFaceIndices)[3])
{
  assert(heToSplit < m_heData.size());
  assert(m_heData[heToSplit].face != HOLE_INDEX &&
         m_heData[m_heData[heToSplit].twin].face == HOLE_INDEX);

  half_edge &heBase = m_heData[heToSplit];
  half_edge &heNext = m_heData[heBase.next];
  half_edge &hePrev = m_heData[heNext.next];
  half_edge &heTwin = m_heData[heBase.twin];

  Eigen::Vector3d newVert1 =
      m_vertices[heBase.vert] +
      (m_vertices[heNext.vert] - m_vertices[heBase.vert]) / 3.0;
  Eigen::Vector3d newVert2 =
      m_vertices[heBase.vert] +
      (m_vertices[heNext.vert] - m_vertices[heBase.vert]) * (2.0 / 3.0);

  // Construct 2 new faces and 8 new half_edges connecting the new verts to the
  // off-edge vert.
  std::size_t newVertIndex = m_vertices.size();
  std::size_t newHeIndex = m_heData.size();
  std::size_t newFaceIndex = m_faceData.size();

  if (pOutVertIndices)
    {
      (*pOutVertIndices)[0] = newVertIndex;
      (*pOutVertIndices)[1] = newVertIndex + 1;
    }

  if (pOutFaceIndices)
    {
      (*pOutFaceIndices)[0] = heBase.face;
      (*pOutFaceIndices)[1] = newFaceIndex;
      (*pOutFaceIndices)[2] = newFaceIndex + 1;
    }

  half_edge newHe[8];
  detail::init(newHe[0], heNext.next, newHeIndex + 1, newVertIndex,
               heBase.face);
  detail::init(newHe[1], newHeIndex + 4, newHeIndex, hePrev.vert, newFaceIndex);
  detail::init(newHe[2], newHeIndex + 1, newHeIndex + 3, newVertIndex + 1,
               newFaceIndex);
  detail::init(newHe[3], newHeIndex + 5, newHeIndex + 2, hePrev.vert,
               newFaceIndex + 1);
  detail::init(newHe[4], newHeIndex + 2, newHeIndex + 6, newVertIndex,
               newFaceIndex);
  detail::init(newHe[5], heBase.next, heBase.twin, newVertIndex + 1,
               newFaceIndex + 1);
  detail::init(newHe[6], newHeIndex + 7, newHeIndex + 4, newVertIndex + 1,
               HOLE_INDEX);
  detail::init(newHe[7], heTwin.next, heToSplit, newVertIndex, HOLE_INDEX);

#ifdef USE_PREV
  newHe[0].prev = heToSplit;
  newHe[1].prev = newHeIndex + 2;
  newHe[2].prev = newHeIndex + 4;
  newHe[3].prev = heBase.next;
  newHe[4].prev = newHeIndex + 1;
  newHe[5].prev = newHeIndex + 3;
  newHe[6].prev = heBase.twin;
  newHe[7].prev = newHeIndex + 6;

  heNext.prev = newHeIndex + 5;
  m_heData[heTwin.next].prev = newHeIndex +
                               7
#endif

                               heBase.next = newHeIndex;
  heBase.twin = newHeIndex + 7;
  heTwin.next = newHeIndex + 6;
  heTwin.twin = newHeIndex + 5;
  heNext.next = newHeIndex + 3;
  heNext.face = newFaceIndex + 1;

  m_vertices.push_back(newVert1);
  m_vertices.push_back(newVert2);
  m_vertData.push_back(newHeIndex + 4);
  m_vertData.push_back(newHeIndex + 5);
  m_faceData[heBase.face] = heToSplit;
  m_faceData.push_back(newHeIndex + 4);
  m_faceData.push_back(newHeIndex + 5);
  m_heData.push_back(newHe[0]);
  m_heData.push_back(newHe[1]);
  m_heData.push_back(newHe[2]);
  m_heData.push_back(newHe[3]);
  m_heData.push_back(newHe[4]);
  m_heData.push_back(newHe[5]);
  m_heData.push_back(newHe[6]);
  m_heData.push_back(newHe[7]);
}

double EditMesh::get_cotan_weight(const vvert_iterator &it) const
{
  const half_edge *itCur = it.m_cur;
  const half_edge *itTwin = &m_heData[itCur->twin];

  double result = 0;

  Eigen::Vector3d a = this->get_vertex(itTwin->vert);
  Eigen::Vector3d b = this->get_vertex(itCur->vert);

  assert((itCur->face != HOLE_INDEX || itTwin->face != HOLE_INDEX) &&
         "Invalid mesh: edge with no face on either side");

  if (itCur->face != HOLE_INDEX)
    {
      Eigen::Vector3d c = this->get_vertex(this->prev(*itCur).vert);
      Eigen::Vector3d e0 = b - c;
      Eigen::Vector3d e1 = a - c;

      // We use the dot product and norm of the cross product to get cos and sin
      // respectively. cotan = cos / sin
      result += static_cast<double>(e0.dot(e1)) / e0.cross(e1).norm();
    }

  if (itTwin->face != HOLE_INDEX)
    {
      Eigen::Vector3d c = this->get_vertex(this->prev(*itTwin).vert);
      Eigen::Vector3d e0 = a - c;
      Eigen::Vector3d e1 = b - c;

      result += static_cast<double>(e0.dot(e1)) / e0.cross(e1).norm();
    }

  return result;
}

double EditMesh::get_mean_value_weight(const vvert_iterator &it) const
{
  const half_edge *itCur = it.m_cur;
  const half_edge *itTwin = &m_heData[itCur->twin];

  double result = 0;

  Eigen::Vector3d a = this->get_vertex(itTwin->vert);
  Eigen::Vector3d b = this->get_vertex(itCur->vert);

  Eigen::Vector3d e0 = (b - a);
  double eLen = e0.norm();

  e0 /= eLen;

  assert((itCur->face != HOLE_INDEX || itTwin->face != HOLE_INDEX) &&
         "Invalid mesh: edge with no face on either side");

  if (itCur->face != HOLE_INDEX)
    {
      Eigen::Vector3d c = this->get_vertex(this->prev(*itCur).vert);
      Eigen::Vector3d e1 = (c - a).normalized();

      result += std::tan(0.5 * std::acos(e0.dot(e1)));
    }

  if (itTwin->face != HOLE_INDEX)
    {
      Eigen::Vector3d c = this->get_vertex(this->prev(*itTwin).vert);
      Eigen::Vector3d e1 = (c - a).normalized();

      result += std::tan(0.5 * std::acos(e0.dot(e1)));
    }

  return result / eLen;
}

Eigen::Vector3d EditMesh::get_normal(const vface_iterator &it) const
{
  Eigen::Vector3d a = this->get_vertex(it.m_next->vert);
  Eigen::Vector3d b = this->get_vertex(it.m_cur->vert);
  Eigen::Vector3d c = this->get_vertex(m_heData[it.m_cur->next].vert);

  return (a - c).cross(b - c).normalized();
}

Eigen::Vector4d EditMesh::get_plane(const vface_iterator &it) const
{
  Eigen::Vector3d a = this->get_vertex(it.m_next->vert);
  Eigen::Vector3d b = this->get_vertex(it.m_cur->vert);
  Eigen::Vector3d c = this->get_vertex(m_heData[it.m_cur->next].vert);
  Eigen::Vector3d n = (a - c).cross(b - c).normalized();

  // Plane equation Ax + By + Cz + D = 0 -> n.dot( [x,y,z] ) - n.dot( c ) = 0
  return Eigen::Vector4d(n.x(), n.y(), n.z(), -n.dot(c));
}

Eigen::Vector3d EditMesh::get_vnormal(std::size_t vertex) const
{
  vvert_iterator vit;
  init_iterator(vit, vertex);

  Eigen::Vector3d normal;
  normal.setZero();

  Eigen::Vector3d center = this->get_vertex(vertex);
  Eigen::Vector3d vec_prev;
  Eigen::Vector3d vec_curr = this->get_vertex(deref_iterator(vit)) - center;

  advance_iterator(vit);
  vit.m_end = vit.m_cur;
  do
    {
      vec_prev = vec_curr;
      vec_curr = this->get_vertex(deref_iterator(vit)) - center;

      if (m_heData[vit.m_cur->twin].face != HOLE_INDEX)
        normal += vec_curr.cross(vec_prev);
    }
  while (advance_iterator(vit));

  return normal.normalized();
}

Eigen::Vector3d EditMesh::get_fnormal(std::size_t face) const
{
  const half_edge *e1 = &m_heData[m_faceData[face]];
  const half_edge *e2 = &this->next(*e1);
  const half_edge *e3 = &this->next(*e2);

  Eigen::Vector3d a = this->get_vertex(e1->vert);
  Eigen::Vector3d b = this->get_vertex(e2->vert);
  Eigen::Vector3d c = this->get_vertex(e3->vert);

  Eigen::Vector3d normal = (b - a).cross(c - a);

  return normal.normalized();
}

void EditMesh::getIndicesForFace(size_t tri_index,
                                 size_t indicesForFace[3]) const
{
  fvert_iterator fvit;
  init_iterator(fvit, tri_index);

  for (size_t i = 0; i < 3; i++)
    {
      indicesForFace[i] = deref_iterator(fvit);
      advance_iterator(fvit);
    }
}

Eigen::Vector3d EditMesh::getFaceMidpoint(size_t tri_index)
{
  fvert_iterator fvit;
  init_iterator(fvit, tri_index);

  Eigen::Vector3d average = Eigen::Vector3d::Zero();
  int count = 0;
  for (size_t i = 0; i < 3; i++)
    {
      size_t index = deref_iterator(fvit);
      if (index != HOLE_INDEX)
        {
          count++;
          average += get_vertex(index);
        }
      advance_iterator(fvit);
    }
  return average / count;
}

void EditMesh::get_draw_data(float *verts, int *indices) const
{

  // for each face
  for (std::size_t f = 0, iEnd = m_faceData.size(); f < iEnd; ++f)
    {
      const half_edge *he = &m_heData[m_faceData[f]];

      // for each vertex of the face
      for (int v = 0; v < 3; ++v)
        {
          Eigen::Vector3d vert = get_vertex(he->vert);
          indices[3 * f + v] = 3 * f + v;

          // for each component of the vertex
          for (int k = 0; k < 3; ++k)
            {
              // this is where we convert from the double-precision of the data
              // structure to float for the graphics card
              verts[3 * (3 * f + v) + k] = static_cast<float>(vert[k]);
            }
          he = &this->next(*he);
        }
    }
}

void EditMesh::get_draw_normals(float *normals) const
{

  // for each face
  for (std::size_t f = 0, iEnd = m_faceData.size(); f < iEnd; ++f)
    {
      Eigen::Vector3d fnormal = this->get_fnormal(f);

      // for each vertex of the face
      for (int v = 0; v < 3; ++v)
        {

          // for each component of the vertex
          for (int k = 0; k < 3; ++k)
            {
              // this is where we convert from the double-precision of the data
              // structure to float for the graphics card
              normals[3 * (3 * f + v) + k] = static_cast<float>(fnormal[k]);
            }
        }
    }
}

void EditMesh::get_draw_selection(int *selection) const
{
  for (std::size_t i = 0, iEnd = m_faceData.size(); i < iEnd; i++)
    {
      size_t verts[3];
      getIndicesForFace(i, verts);

      for (int j = 0; j < 3; j++)
        {
          selection[3 * i + j] = m_selected[verts[j]];
        }
    }
}

void EditMesh::get_draw_colors(float *colors) const
{
  // for each face
  for (std::size_t f = 0, iEnd = m_faceData.size(); f < iEnd; ++f)
    {
      size_t he_idx = m_faceData[f];

      // for each vertex of the face
      for (int v = 0; v < 3; ++v)
        {

          float color[3] = {-1.0f, -1.0f, -1.0f}; // negative colors will be
                                                  // discarded by the shader and
                                                  // normal lighting will be
                                                  // applied
          // if (f % 20 == 0) { // example: paint every 20th face
          //	color[0] = 0.0f;
          //	color[1] = 1.0f;
          //	color[2] = 1.0f;
          //}

          // if (m_heData[he_idx].vert % 20 == 0) { // example: paint every 20th
          // vertex
          //	color[0] = 1.0f;
          //	color[1] = 0.0f;
          //	color[2] = 1.0f;
          //}

          // for each component of the vertex
          for (int k = 0; k < 3; ++k)
            {
              colors[3 * (3 * f + v) + k] = color[k];
            }
          he_idx = m_heData[he_idx].next;
        }
    }
}

// call instead of init to test edge flip
// easiest way is hacking it into mesh constructor
void EditMesh::test_flip()
{
  std::vector<double> xyz;
  std::vector<std::size_t> faces;

  // four verts
  xyz.push_back(-1);
  xyz.push_back(0);
  xyz.push_back(0);
  xyz.push_back(0);
  xyz.push_back(1);
  xyz.push_back(1);
  xyz.push_back(0);
  xyz.push_back(1);
  xyz.push_back(-1);
  xyz.push_back(1);
  xyz.push_back(0);
  xyz.push_back(0);

  // two triangles
  faces.push_back(0);
  faces.push_back(1);
  faces.push_back(2);
  faces.push_back(2);
  faces.push_back(1);
  faces.push_back(3);

  this->init(xyz, faces);

  half_edge *he = &m_heData[0];
  for (size_t i = 0; i < m_heData.size(); ++i)
    {
      he = &m_heData[i];
      if (he->face != HOLE_INDEX && m_heData[he->twin].face != HOLE_INDEX)
        break;
    }
  flip_edge(*he);
  this->edit_count++;
}

void EditMesh::write_to_obj_stream(std::ostream &stream) const
{
  for (auto &v : m_vertices)
    stream << "v " << v.x() << ' ' << v.y() << ' ' << v.z() << std::endl;
  stream << std::endl;
  for (std::size_t i = 0, iEnd = m_faceData.size(); i < iEnd; ++i)
    {
      fvert_iterator it;
      this->init_iterator(it, i);
      stream << "f ";
      bool isFirst = true;
      do
        {
          if (!isFirst)
            stream << ' ';
          isFirst = false;
          stream << this->deref_iterator(it) + 1;
        }
      while (this->advance_iterator(it));
      stream << std::endl;
    }
}

void EditMesh::verify() const
{
  for (std::size_t i = 0, iEnd = m_faceData.size(); i < iEnd; ++i)
    {
      std::size_t c = 0;

      const half_edge *it = &m_heData[m_faceData[i]];
      assert(it->next != m_faceData[i]);
      while (it->next != m_faceData[i])
        {
          assert(it->face == i);
          assert(it->next != HOLE_INDEX && it->twin != HOLE_INDEX &&
                 it->vert < m_vertData.size());
          assert((m_heData[it->twin].face == HOLE_INDEX ||
                  m_heData[m_heData[it->next].twin].face !=
                      m_heData[it->twin].face) &&
                 "Can't have two edges shared between the same faces!");
          it = &m_heData[it->next];
          assert(++c < 1000000); // This isn't strictly a problem, but probably
                                 // no face has a million verts in it.
        }
    }

  for (std::size_t i = 0, iEnd = m_vertData.size(); i < iEnd; ++i)
    {
      assert(m_vertData[i] == HOLE_INDEX || m_vertData[i] < m_heData.size());
      if (m_vertData[i] != HOLE_INDEX)
        {
          const half_edge *it = &m_heData[m_vertData[i]];
          assert(it->vert == i);
          assert(it->face != HOLE_INDEX &&
                 "By convention, vertices should not reference hole faces");
        }
    }

  for (std::size_t i = 0, iEnd = m_heData.size(); i < iEnd; ++i)
    {
      const half_edge *it = &m_heData[i];
      assert(it->vert < m_vertData.size());
      assert(it->face == HOLE_INDEX || it->face < m_faceData.size());

      assert(it->next < m_heData.size());
      assert(it->next != i);
      assert(m_heData[it->next].face == it->face);
      assert(m_heData[it->next].vert != it->vert);

      assert(it->twin < m_heData.size());
      assert(m_heData[it->twin].twin == i);
      assert(m_heData[it->twin].face != it->face);
      assert(m_heData[it->twin].vert == m_heData[it->next].vert);

#ifdef USE_PREV
      assert(it->prev < m_heData.size());
      assert(it->prev != i);
      assert(m_heData[it->next].prev == i);
      assert(m_heData[it->prev].next == i);
      assert(m_heData[it->prev].face == it->face);
      assert(m_heData[it->prev].vert != it->vert);
#endif
    }
}

void EditMesh::test()
{
  EditMesh m1, m2, m3;

  std::vector<double> v;
  std::vector<std::size_t> f;

  v.push_back(0);
  v.push_back(0);
  v.push_back(0);
  v.push_back(1);
  v.push_back(0);
  v.push_back(0);
  v.push_back(0);
  v.push_back(1);
  v.push_back(0);
  v.push_back(1);
  v.push_back(1);
  v.push_back(0);
  v.push_back(2);
  v.push_back(1);
  v.push_back(0);
  v.push_back(1);
  v.push_back(2);
  v.push_back(0);

  f.push_back(0);
  f.push_back(1);
  f.push_back(2);
  f.push_back(2);
  f.push_back(1);
  f.push_back(3);
  f.push_back(3);
  f.push_back(1);
  f.push_back(4);
  f.push_back(4);
  f.push_back(5);
  f.push_back(3);
  f.push_back(3);
  f.push_back(5);
  f.push_back(2);

  m1.init(v, f);

  for (std::size_t i = 0, iEnd = v.size(); i < iEnd; i += 3)
    assert(m2.add_vertex(v[i], v[i + 1], v[i + 2]) == i / 3);
  for (std::size_t i = 0, iEnd = f.size(); i < iEnd; i += 3)
    assert(m2.add_face(*reinterpret_cast<std::size_t(*)[3]>(&f[i])) == i / 3);

  assert(m1.get_face_size() == m2.get_face_size());
  assert(m1.get_vert_size() == m2.get_vert_size());

  m1.verify();
  m2.verify();

  m2.delete_face(0);
  m2.verify();

  m2.delete_face(0); // Was face 4 originally
  m2.verify();

  m2.delete_face(0); // Was face 3 originally
  m2.verify();

  m2.delete_face(0); // Was face 2 originally
  m2.verify();

  m2.delete_face(0); // Was face 1 originally
  m2.verify();

  assert(m2.get_face_size() == 0);

  m2.add_face(0, 1, 2);
  m2.verify();

  m2.split_face_center(0);
  m2.verify();

  assert(m2.get_face_size() == 3);

  m2.clear();
  m2.add_vertex(0, 0, 0);
  m2.add_vertex(1, 0, 0);
  m2.add_vertex(0, 1, 0);
  m2.add_face(0, 1, 2);
  m2.split_boundary_edge(m2.find_twin(0, 1)->twin);
  m2.verify();

  assert(m2.get_face_size() == 3);

  m3.add_vertex(-1, 0, 0);
  m3.add_vertex(1, 0, 0);
  m3.add_vertex(0, 1, 0);
  m3.add_vertex(0, -1, 0);
  m3.add_vertex(-1, 1, 0);
  m3.add_vertex(-1, -1, 0);
  m3.add_vertex(1, 1, 0);
  m3.add_vertex(1, -1, 0);

  m3.add_face(0, 1, 2);
  m3.add_face(0, 2, 4);
  m3.add_face(0, 4, 5);
  m3.add_face(0, 5, 3);
  m3.add_face(0, 3, 1);

  m3.add_face(1, 3, 7);
  m3.add_face(1, 7, 6);
  m3.add_face(1, 6, 2);
  m3.verify();

  m3.flip_edge(*m3.find_edge(1, 0));
  m3.verify();
  m3.flip_edge(*m3.find_edge(2, 3));
  m3.verify();

  std::size_t newVert = m3.collapse_edge(m3.find_edge(1, 0)->twin);
  assert(newVert != HOLE_INDEX);
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
  m1.add_vertex(-1, 0, 0);
  m1.add_vertex(1, 0, 0);
  m1.add_vertex(0, 1, 0);
  m1.add_vertex(-2, 1, 0);
  m1.add_vertex(-2, -1, 0);
  m1.add_vertex(0, -1, 0);
  m1.add_vertex(2, -1, 0);
  m1.add_vertex(2, 1, 0);
  m1.add_vertex(-2, 2, 0);
  m1.add_vertex(-3, 1, 0);
  m1.add_vertex(-3, -1, 0);
  m1.add_vertex(-2, -2, 0);
  m1.add_vertex(2, -2, 0);
  m1.add_vertex(3, -1, 0);
  m1.add_vertex(3, 1, 0);
  m1.add_vertex(2, 2, 0);

  m1.add_face(0, 1, 2);
  m1.add_face(0, 2, 3);
  m1.add_face(0, 3, 4);
  m1.add_face(0, 4, 5);
  m1.add_face(0, 5, 1);
  m1.add_face(1, 5, 6);
  m1.add_face(1, 6, 7);
  m1.add_face(1, 7, 2);
  m1.add_face(3, 2, 8);
  m1.add_face(3, 8, 9);
  m1.add_face(3, 9, 10);
  m1.add_face(3, 10, 4);
  m1.add_face(4, 10, 11);
  m1.add_face(4, 11, 5);
  m1.add_face(5, 11, 12);
  m1.add_face(5, 12, 6);
  m1.add_face(6, 12, 13);
  m1.add_face(6, 13, 14);
  m1.add_face(6, 14, 7);
  m1.add_face(7, 14, 15);
  m1.add_face(7, 15, 2);
  m1.add_face(2, 15, 8);

  m1.verify();

  std::vector<std::size_t> faces;
  vface_iterator it;
  if (m1.init_iterator(it, 0))
    {
      do
        {
          std::vector<std::size_t>::iterator itInsert = std::lower_bound(
              faces.begin(), faces.end(), m1.deref_iterator(it),
              std::greater<std::size_t>());
          if (itInsert == faces.end() || *itInsert != m1.deref_iterator(it))
            faces.insert(itInsert, m1.deref_iterator(it));
        }
      while (m1.advance_iterator(it));
    }
  if (m1.init_iterator(it, 1))
    {
      do
        {
          std::vector<std::size_t>::iterator itInsert = std::lower_bound(
              faces.begin(), faces.end(), m1.deref_iterator(it),
              std::greater<std::size_t>());
          if (itInsert == faces.end() || *itInsert != m1.deref_iterator(it))
            faces.insert(itInsert, m1.deref_iterator(it));
        }
      while (m1.advance_iterator(it));
    }

  std::swap(faces[faces.size() - 1], faces[faces.size() - 2]);

  for (auto face : faces)
    {
      m1.delete_face(face);
      m1.verify();
    }
}

bool EditMesh::is_safe_addface(std::size_t v1, std::size_t v2, std::size_t v3)
{

  std::set<std::size_t> vv;
  vv.insert(m_vertData[v1]);
  vv.insert(m_vertData[v2]);
  vv.insert(m_vertData[v3]);

  // if there's one disconnected vertex its safe
  // more than one is not
  if (vv.size() < 3)
    return false;
  else if (vv.count(HOLE_INDEX) > 0)
    return true;

  vv.clear();
  vv.insert(v1);
  vv.insert(v2);
  vv.insert(v3);

  int count = 0;
  vvert_iterator vit;
  init_iterator(vit, v1);
  do
    {
      if (vit.m_cur->vert == v2 || vit.m_cur->vert == v3)
        {
          count++;
          break;
        }
    }
  while (this->advance_iterator(vit));

  init_iterator(vit, v2);
  do
    {
      if (vit.m_cur->vert == v1 || vit.m_cur->vert == v3)
        {
          count++;
          break;
        }
    }
  while (this->advance_iterator(vit));

  // if two of the vertices have a next within the triplet, its safe
  if (count > 1)
    return true;

  // else unsafe
  return false;
  // return (v1 != HOLE_INDEX && v2 != HOLE_INDEX) ||
  //       (v2 != HOLE_INDEX && v3 != HOLE_INDEX) ||
  //       (v3 != HOLE_INDEX && v1 != HOLE_INDEX);
}

void EditMesh::example()
{
  // iterate over all halfedges of the mesh
  for (size_t i = 0; i < m_heData.size(); ++i)
    {
      const half_edge *he = &m_heData[i];
    }

  // iterate over all faces of the mesh
  for (size_t f = 0; f < m_faceData.size(); ++f)
    {
      // one of its halfedges is m_faceData[i];
    }

  // iterate over all vertices of the mesh
  for (size_t v = 0; v < m_vertData.size(); ++v)
    {
      // one of its halfedges is m_vertData[i];
    }

  // iterate over first ring of a vertex
  size_t vFrom = 42;
  vvert_iterator it;
  if (!this->init_iterator(it, vFrom))
    { /* error, vertex doesn't exist or has no neighbors */
    }
  do
    {
      size_t vTo = deref_iterator(it);
    }
  while (this->advance_iterator(it));

  // similarly for other iterators (vface, fvert, fface)

  // get the 3 halfedges of a face
  {
    size_t face = 42;
    size_t he_idx = m_faceData[face];
    const half_edge *he1 = &m_heData[he_idx];
    const half_edge *he2 = &next(*he1);
    const half_edge *he3 = &next(*he2);
  }
}

void EditMesh::simplify_vertex_removal(int number_operations) {}

void EditMesh::simplify_edge_collapse(int number_operations) {}
