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
#pragma once
#ifndef OBJECT_MESH_H
#define OBJECT_MESH_H

#include "DrawMesh.h"
#include "EditMesh.h"
#include <memory>

typedef std::shared_ptr<EditMesh> EditMesh_ptr;
typedef std::shared_ptr<DrawMesh> DrawMesh_ptr;

class Model
{
public:
  Model();
  Model(EditMesh_ptr em);
  ~Model();

  void init(RenderState &state);

  void drawMesh(int primitive = GL_TRIANGLES);

  /***********************************************************
   * the reflection interface to get mesh information
   ***********************************************************/
  std::size_t info_sizev() { return m_em->get_vert_size(); }
  std::size_t info_sizef() { return m_em->get_face_size(); }
  Eigen::Vector3d info_vertex(std::size_t i);

  void info_bbox(Eigen::Vector3d &bboxMin, Eigen::Vector3d &bboxMax);
  void getIndicesForFace(size_t tri_index, size_t indicesForFace[3]);

  /***********************************************
   * mesh modification algorithms
   ***********************************************/
  // expose your editmesh functions here

  // selection functions
  void select_vert(size_t index) { m_em->select_vert(index); }
  void deselect_vert(size_t index) { m_em->deselect_vert(index); }
  void deselect_verts() { m_em->deselect_allVerts(); }
  bool isSelected(size_t index) { return m_em->isSelected(index); }

  const EditMesh_ptr get_editMesh() const;
  void set_editMesh(EditMesh_ptr em);

protected:
  void updateDrawMesh();

  /***********************************************
   * Variables for this Instance of the EditMesh
   ***********************************************/
  // the bounding box for this instance of the geometry
  Eigen::Vector3f bboxMin;
  Eigen::Vector3f bboxMax;

  // indicates the last drawn edit mesh, if this does not match
  // the value in m_em then the draw mesh needs to be updated
  int last_drawn;

  DrawMesh_ptr m_dm;
  EditMesh_ptr m_em;
};

inline const EditMesh_ptr Model::get_editMesh() const { return m_em; }

inline void Model::set_editMesh(EditMesh_ptr em) { m_em = em; }

#endif // OBJECT_MESH_H