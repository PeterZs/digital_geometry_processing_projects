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
#include "Model.h"

Model::Model() : last_drawn(-1), m_em(NULL), m_dm(NULL)
{
  bboxMin.setZero();
  bboxMax.setZero();
}

Model::Model(EditMesh_ptr em) : last_drawn(-1), m_em(em), m_dm(NULL)
{
  bboxMin.setZero();
  bboxMax.setZero();
}

Model::~Model()
{
  // Shared pointers don't need to be deleted?
}

void Model::init(RenderState &state)
{
  if (!m_em)
    m_em = EditMesh_ptr(new EditMesh());
  m_dm = DrawMesh_ptr(new DrawMesh(state));
  m_dm->init(1);
}

void Model::drawMesh(int primitive)
{
  if (last_drawn != m_em->get_edit_count())
    updateDrawMesh();

  m_dm->drawMesh(primitive);
}

Eigen::Vector3d Model::info_vertex(std::size_t i)
{
  return m_em->get_vertex(i);
}

void Model::info_bbox(Eigen::Vector3d &bboxMin, Eigen::Vector3d &bboxMax)
{
  bboxMin = m_em->bboxMin;
  bboxMax = m_em->bboxMax;
}

void Model::getIndicesForFace(size_t tri_index, size_t indicesForFace[3])
{
  m_em->getIndicesForFace(tri_index, indicesForFace);
}

void Model::updateDrawMesh()
{
  int i_size = 3 * m_em->get_face_size();
  int v_size = 3 * i_size; // m_em->get_vert_size();

  float *v_data = new float[v_size * 2]; // vertex and normal data
  int *i_data = new int[i_size];         // face indexes
  int *s_data = new int[i_size];         // selection data
  float *c_data = new float[v_size];     // vertex colors

  m_em->get_draw_data(v_data, i_data);
  m_em->get_draw_normals(&v_data[v_size]);
  m_em->get_draw_selection(s_data);
  m_em->get_draw_colors(c_data);

// i don't want to allocate on the heap, nor do I want to
// hard code all the values so this is my disgusting solution
#define DRAW_MESH_NUM_CORE_ATTR 2
#define DRAW_MESH_NUM_EXTRA_ATTR 1

  attrib_info attr_info[DRAW_MESH_NUM_CORE_ATTR];
  attr_info[0].attrib_number = 0; // vertices are 0
  attr_info[0].attrib_size = sizeof(float);
  attr_info[0].data_offset = 0; // data starts at beginning of array
  attr_info[0].data_stride = 0; // data is tightly packed
  attr_info[0].num_comp = 3;    // there are 3 components per vertex position

  attr_info[1].attrib_number = 1; // normals are 1
  attr_info[1].attrib_size = sizeof(float);
  attr_info[1].data_offset =
      v_size * sizeof(float);   // data starts after vertices
  attr_info[1].data_stride = 0; // data is tightly packed
  attr_info[1].num_comp = 3;    // there are 3 components per vertex normal

  m_dm->loadVBuffer(0, sizeof(float) * 2 * v_size, (GLubyte *)v_data, 0,
                    DRAW_MESH_NUM_CORE_ATTR, attr_info);

  attr_info[0].attrib_number = 2; // vertex selection is 2
  attr_info[0].attrib_size = sizeof(int);
  attr_info[0].data_offset = 0; // data starts at the beginning of array
  attr_info[0].data_stride = 0; // data is tightly packed
  attr_info[0].num_comp = 1;    // there is 1 components per vertex selection

  m_dm->loadVBuffer(1, sizeof(int) * i_size, (GLubyte *)s_data, 0,
                    DRAW_MESH_NUM_EXTRA_ATTR, attr_info);

  attr_info[0].attrib_number = 3; // color is 3
  attr_info[0].attrib_size = sizeof(float);
  attr_info[0].data_offset = 0; // data starts at the beginning of array
  attr_info[0].data_stride = 0; // data is tightly packed
  attr_info[0].num_comp = 3;    // there are 3 components (r,g,b) per vertex

  m_dm->loadVBuffer(2, sizeof(float) * v_size, (GLubyte *)c_data, 0, 1,
                    attr_info);

  m_dm->loadIBuffer(i_size, sizeof(int), i_data);

#undef DRAW_MESH_NUM_ATTR
  delete[] v_data;
  delete[] i_data;
  delete[] s_data;
  delete[] c_data;

  last_drawn = m_em->get_edit_count();
}
