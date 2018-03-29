/* Copyright (c) Russell Gillette
 * December 2013
 *
 * CPSC Computer Graphics 2014
 * University of British Columbia
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

/* == Utils.h ==
 *
 * General use common/helper functions
 */

#ifndef UTILS_H
#define UTILS_H

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <Eigen/Core>
#include <glm/glm.hpp>
#include <stdio.h>
#include <string>

using Eigen::Vector3d;

template <class T> inline bool approx_eq(T a, T b, double precision)
{
  return abs(a - b) < precision;
}

inline double clamp(double value, double min, double max)
{
  return std::min(std::max(value, min), max);
}

/* Takes the location and file name of the file to be processed
 * and returns the contents stored within a GLchar array
 */
inline const GLchar *loadFileAsString(const char *path, GLint *size)
{
  GLchar *file_contents = NULL;
  int file_size = 0;

  FILE *f = fopen(path, "rb");

  fseek(f, 0, SEEK_END);
  file_size = ftell(f);

  file_contents = new GLchar[file_size + 1];
  fseek(f, 0, SEEK_SET);

  fread(file_contents, file_size, 1, f);

  fclose(f);
  file_contents[file_size] = 0; // null terminate the file

  if (size != NULL)
    *size = file_size + 1;

  return file_contents;
}

inline void realign_triangle(Vector3d &edge1, Vector3d &edge2)
{
  double edge_len[2]; // edge lengths of the two triangles, used to determine
                      // the 2d representation
  double angle;       // angle between to edge lengths to complete the triangle
                      // representation

  // calculate the edge lengths
  edge_len[0] = edge1.norm();
  edge_len[1] = edge2.norm();

  // calculate the angles between edge vectors
  angle =
      edge1.dot(edge2) /
      (edge_len[0] * edge_len[1]); // cos(theta) to avoid costly calculations

  // realign the triangles to be horizontal
  edge1 = Vector3d(edge_len[0], 0, 0);
  edge2 = Vector3d(edge_len[1] * angle, edge_len[1] * (sqrt(1 - angle * angle)),
                   0); // (e_2*cos(90 - theta), e_2*sin(90 - theta))
}

inline void realign_triangle(Vector3d &vert1, Vector3d &vert2, Vector3d &vert3)
{
  // calculate edge vectors of triangles
  vert2 = vert2 - vert1;
  vert3 = vert3 - vert1;

  // realign the triangles using edges to be horizontal
  vert1 = Vector3d::Zero();
  realign_triangle(vert2, vert3);
}

// returns the greatest positive root
inline Eigen::Matrix2d mat2_sqrt(Eigen::Matrix2d M)
{
  double tau = M(0, 0) + M(1, 1);
  double delta = M(0, 0) * M(1, 1) - M(0, 1) * M(1, 0);
  double s = 0, t = 0;

  if (delta < 0)
    {
      printf("ERROR: imaginary root of square matrix\n");
      return Eigen::Matrix2d::Identity(); // return identity
    }
  s = sqrt(delta);
  t = tau + 2 * s;
  if (t < 0)
    {
      printf("ERROR: imaginary root of square matrix\n");
      return Eigen::Matrix2d::Identity(); // return identity
    }
  t = sqrt(t);

  M(0, 0) += s;
  M(1, 1) += s;
  M /= t;

  return M;
}

// NOTE: stride is in sizeof(float) NOT bytes
inline void triangles_to_file(float **vert_array, int *index_array,
                              int num_elem, int num_frames, char *path,
                              unsigned int stride = 0)
{
  char filename[50] = ""; // allocate placeholder memory on the stack
  FILE *f = NULL;

  int num_tri = num_elem / 3;
  int v_stride = (stride == 0) ? 3 : stride;
  int offset = 0;

  int s_1, s_2, s_3;

  for (int i = 0; i < num_tri; ++i)
    {
      sprintf(filename, "%s/triangle_%d.txt", path, i);
      f = fopen(filename, "w+");

      for (int j = 0; j < num_frames; j++)
        {
          if (index_array) // select out vertices from index array
            {
              s_1 = index_array[i * 3] * v_stride;
              s_2 = index_array[i * 3 + 1] * v_stride;
              s_3 = index_array[i * 3 + 2] * v_stride;
            }
          else // assume triangles are listed sequentially
            {
              s_1 = offset;
              s_2 = s_1 + v_stride;
              s_3 = s_2 + v_stride;
            }

          if (j != 0)
            fprintf(f, "\n");

          fprintf(f, "%f,%f,%f,%f,%f,%f,%f,%f,%f", vert_array[j][s_1],
                  vert_array[j][s_1 + 1], vert_array[j][s_1 + 2],
                  vert_array[j][s_2], vert_array[j][s_2 + 1],
                  vert_array[j][s_2 + 2], vert_array[j][s_3],
                  vert_array[j][s_3 + 1], vert_array[j][s_3 + 2]);
        }
      fclose(f);
      offset += 3 * v_stride;
    }
}

inline void tensor_toColor(double &tensor, glm::vec3 &color)
{
  if (tensor > 1)
    color = glm::vec3(1, 1 / tensor, 1 / tensor);
  else
    color = glm::vec3(tensor, tensor, 1);
}

inline void tensor_toColor1Frame(double &tensor, glm::vec3 &color)
{
  // if you wish to change the scaling, the formula should be:
  // a*tensor - (a-1).
  tensor = 4 * tensor - 3;
  tensor_toColor(tensor, color);
}

inline bool ray_intersect_triangle(const Vector3d &ray_origin,
                                   const Vector3d &point_on_ray,
                                   Vector3d triangle_verts[3])
{
  Vector3d d = point_on_ray - ray_origin;
  Vector3d e1 = triangle_verts[1] - triangle_verts[0];
  Vector3d e2 = triangle_verts[2] - triangle_verts[0];
  Vector3d h = d.cross(e2);
  double a = e1.dot(h);

  if (a > -0.00001 && a < 0.00001)
    return false;

  double f = 1 / a;
  Vector3d s = ray_origin - triangle_verts[0];
  double u = f * s.dot(h);
  if (u < 0.0 || u > 1.0)
    return false;

  Vector3d q = s.cross(e1);
  double v = f * d.dot(q);
  if (v < 0.0 || u + v > 1.0)
    return false;

  double t = f * e2.dot(q);
  if (t > 0.00001)
    return true;

  return false;
}

inline bool
vert_inside_select_box(const Eigen::Vector3d &bot_left_origin,
                       const Eigen::Vector3d &point_on_bot_left_ray,
                       const Eigen::Vector3d &bot_right_origin,
                       const Eigen::Vector3d &point_on_bot_right_ray,
                       const Eigen::Vector3d &top_right_origin,
                       const Eigen::Vector3d &point_on_top_right_ray,
                       const Eigen::Vector3d &top_left_origin,
                       const Eigen::Vector3d &point_on_top_left_ray,
                       Eigen::Vector3d vertex)
{
  Eigen::Vector3d ray, tmp, normal_left, normal_right, normal_bot, normal_top;
  ray = point_on_bot_left_ray - bot_left_origin;
  tmp = top_left_origin - bot_left_origin;
  normal_left = ray.cross(tmp);

  ray = point_on_bot_right_ray - bot_right_origin;
  tmp = bot_left_origin - bot_right_origin;
  normal_bot = ray.cross(tmp);

  ray = point_on_top_right_ray - top_right_origin;
  tmp = bot_right_origin - top_right_origin;
  normal_right = ray.cross(tmp);

  ray = point_on_top_left_ray - top_left_origin;
  tmp = top_right_origin - top_left_origin;
  normal_top = ray.cross(tmp);

  Eigen::Vector3d bot_left_check = vertex - bot_left_origin;
  Eigen::Vector3d top_right_check = vertex - top_right_origin;

  return (bot_left_check.dot(normal_left) > 0 &&
          bot_left_check.dot(normal_bot) > 0 &&
          top_right_check.dot(normal_right) > 0 &&
          top_right_check.dot(normal_top) > 0);
}

/**
 * compute the normal of a triangle specified by its three
 * vertices in CCW order
 */
inline Vector3d getTriangleNormal(Eigen::Vector3d *verts)
{
  Vector3d a = verts[2] - verts[0];
  Vector3d b = verts[2] - verts[1];
  return (a.cross(b)).normalized();
}

inline double compute_triangle_area(Eigen::Vector3d *verts)
{
  //|(v3-v1) X (v3-v2)| / 2
  Vector3d a = verts[2] - verts[0];
  Vector3d b = verts[2] - verts[1];
  return (a.cross(b)).norm() / 2.0;
}

inline void compute_eigenvalues(const Eigen::Matrix2d &U, double &min,
                                double &max)
{
  double eigenval[2];
  double tmp = U(0, 0) + U(1, 1);
  double tmp2 = pow(tmp, 2) - 4 * (U(0, 0) * U(1, 1) - U(0, 1) * U(1, 0));
  tmp2 = (tmp2 < 0) ? 0 : sqrt(tmp2); // NOTE: since matrix is PSD, anything
                                      // negative would be floating point error
  eigenval[0] = (tmp + tmp2) / 2;     // factor
  eigenval[1] = (tmp - tmp2) / 2;

  min = (eigenval[0] < eigenval[1]) ? eigenval[0] : eigenval[1];
  max = (eigenval[0] > eigenval[1]) ? eigenval[0] : eigenval[1];
}

/**
 * compute the Eigenvectors of the matrix U
 * \param U the matrix who's eigenvectors are to be computed
 * \param eig_val eigen values of U
 * \param term if there are two matching eigen values for the matrix, 'term'
 * will indicate
 * which of the two eigen vectors to compute
 */
inline Eigen::Vector2d compute_eigenvector(const Eigen::Matrix2d &U,
                                           double eig_value, bool term = 0)
{
  double norm;
  Eigen::Vector2d v;

  // you only end up with identical eigenvalues in the case that
  // the matrix is diagonal, in which case the eigen values can be
  // set as the unit vectors for that row

  if (std::abs(U(0, 1)) < 1e-10)
    { // if diagonal
      if (std::abs(U(0, 0) - U(1, 1)) < 1e-7)
        { // if eigenvalues are the same
          v[0] = (term) ? 0 : 1;
          v[1] = (term) ? 1 : 0;
        }
      else if (std::abs(eig_value - U(0, 0)) < 1e-7)
        {
          v = Eigen::Vector2d(1, 0);
        }
      else
        {
          v = Eigen::Vector2d(0, 1);
        }
    }
  else
    {
      v[0] = U(0, 1);
      v[1] = eig_value - U(0, 0);
      norm = sqrt(v[0] * v[0] + v[1] * v[1]);

      v[0] = v[0] / norm;
      v[1] = v[1] / norm;
    }

  return v;
}

inline void compute_eigenvectors(const Eigen::Matrix2d &U, double min,
                                 double max, Eigen::Vector2d &v1,
                                 Eigen::Vector2d &v2)
{
  v1 = compute_eigenvector(U, min, 0);
  v2 = compute_eigenvector(U, max, 1);

  if (std::abs(v1.dot(v2)) > 0.1)
    {
      printf("Eigenvectors evaluating to not-orthogonal:\n"
             "U = [%lf, %lf\n     %lf, %lf]\n",
             U(0, 0), U(0, 1), U(1, 0), U(1, 1));
    }
}

/**
 * project the vector "v" onto a plane with normal vector "normal"
 */
inline Vector3d projectPlane(const Vector3d &v, const Vector3d &normal)
{
  double len_norm = normal.norm();
  return v - (v.dot(normal) / (len_norm * len_norm)) * normal;
}

/**
 * build the rotation matrix from a plane with the normal 'from' to a plane
 * with normal 'to'
 \param from - the normal of the plane that the vector is currently in
 \param to   - the normal of the plane that the vector is rotating to
 \return zero matrix if from or to is invalid (ie NAN), otherwise rotation
 matrix
 */
inline Eigen::Matrix3d buildRotation(const Vector3d &from, const Vector3d &to)
{
  Vector3d rot_axis = from.cross(to);
  double cosTheta = from.dot(to);
  double sinTheta = rot_axis.norm();

  bool isInvalid = false;
  for (int i = 0; i < 3; i++)
    isInvalid |= (from[i] != from[i] || to[i] != to[i]);

  if (rot_axis.norm() < 1e-13)
    return Eigen::Matrix3d::Identity();
  else if (isInvalid)
    return Eigen::Matrix3d::Zero();

  rot_axis.normalize();

  // construct the cross-product matrix
  Eigen::Matrix3d k;
  k << 0, -rot_axis[2], rot_axis[1], rot_axis[2], 0, -rot_axis[0], -rot_axis[1],
      rot_axis[0], 0;

  // construct the rotation matrix for the eigenvector
  Eigen::Matrix3d rot = Eigen::Matrix3d::Identity();
  rot += sinTheta * k + (1 - cosTheta) * k * k;
  return rot;
}

/**
 * rotate vector from one plane to another
 \param v    - vector to rotate
 \param from - the normal of the plane that the vector is currently in
 \param to   - the normal of the plane that the vector is rotating to
 */
inline Vector3d rotateVector(const Vector3d &v, const Vector3d &from,
                             const Vector3d &to)
{
  return buildRotation(from, to) * v;
}

/**
 * check if a vector "v" already within the plane of a cone defined by vectors
 * 'left'
 * and 'right' falls within the cone.
 */
inline bool withinCone(Vector3d &v_in, const Vector3d &left,
                       const Vector3d &right)
{

  // if the vector is withing the cone, then both edges will be on
  // opposite sides of the vector and the negative of one of the
  // edge vectors will not be between the two other vectors.
  double either_side =
      (right.cross(v_in).normalized()).dot(left.cross(v_in).normalized());

  double negative_edge = ((-right).cross(v_in).normalized())
                             .dot((-right).cross(left).normalized());

  // if within triangle
  if (either_side < 0 && negative_edge > 0)
    return true;

  return false;
}

/**
 * project the vector "v" onto the plane of a triangle defined by vectors 'left'
 * and 'right'. Will update v only if returns true.
 * \return If the vector is within the cone defined by 'left' and 'right'
 * L  v
 * | /
 * |/___R
 */
inline bool projectCone(Vector3d &v_in, const Vector3d &left,
                        const Vector3d &right)
{
  Vector3d l = left.normalized();
  Vector3d r = right.normalized();
  Vector3d v = projectPlane(v_in, r.cross(l)).normalized();

  if (v != v)
    {
      printf("v is orthogonal to projected plane or left and right are "
             "colinear\n");
      return false;
    }

  if (withinCone(v, left, right))
    {
      v_in = v;
      return true;
    }

  return false;
}

/**
 * find the nearest point on a given line segment [p1,p2] from point pt
 */
inline Vector3d lineSeg_nearestPt(Vector3d p1, Vector3d p2, Vector3d pt)
{
  Vector3d e1 = p2 - p1;
  Vector3d e2 = pt - p1;

  double proj1 = e1.dot(e2);
  double proj2 = e2.dot(e2);
  double t = proj1 / proj2;

  if (proj1 == 0 || t <= 0)
    {
      return p1;
    }
  else if (proj2 == 0 || t >= 1)
    {
      return p2;
    }
  else
    {
      return p1 + t * e1;
    }
}

/**
 * find the orthogonal (shortest) distance to a given line segment [p1,p2] from
 * point pt
 */
inline double lineSeg_distance(Vector3d p1, Vector3d p2, Vector3d pt)
{
  return (pt - lineSeg_nearestPt(p1, p2, pt)).norm();
}

/**
 * compute the distance 'v' along ray0 that the two rays collide,
 * returns a negative value if there is no collision (or collision is negative)
 */
inline double compute_rayCollision(const Vector3d &p0, const Vector3d &d0,
                                   const Vector3d &p1, const Vector3d &d1)
{
  // check for same origin
  if (p0 == p1)
    return 0;

  bool colinear = false;
  double denom = 0;
  double t = 0, v = 0;

  if (std::abs(denom = d0[0] * d1[1] - d0[1] * d1[0]) > 1e-7)
    { // using x and y
      t = (d1[1] * (p1[0] - p0[0]) - d1[0] * (p1[1] - p0[1])) / denom;
      v = (t * d0[1] + (p0[1] - p1[1])) / d1[1];
    }
  else if (std::abs(denom = d0[0] * d1[2] - d0[2] * d1[0]) > 1e-7)
    { // using x and z
      t = (d1[2] * (p1[0] - p0[0]) - d1[0] * (p1[2] - p0[2])) / denom;
      v = (t * d0[2] + (p0[2] - p1[2])) / d1[2];
    }
  else if (std::abs(denom = d0[1] * d1[2] - d0[2] * d1[1]) > 1e-7)
    { // using y and z
      t = (d1[2] * (p1[1] - p0[1]) - d1[1] * (p1[2] - p0[2])) / denom;
      v = (t * d0[2] + (p0[2] - p1[2])) / d1[2];
    }
  else
    { // lines are parallel, and we need to check if they are colinear
      // NOTE: in the case of colinear, the point of intersection is everywhere
      // in order for our functionality to work, we assum the first point of
      // collision
      // i.e. at parameter = 0

      Vector3d diff = p1 - p0;
      Vector3d dist = diff.cwiseQuotient(
          d0); // distance 'dist' along ray0 until hitting point p1
      double u = 0;

      // find the parameter value along ray 0 for one component
      for (int i = 0; i < 3; i++)
        {
          if (d0[i] != 0)
            {
              u = dist[i]; // assume colinear, validated in next step
              t = 0;
            }
        }

      // validate that it is the same for all components
      for (int i = 0; i < 3; i++)
        {
          if ((d0[i] == 0 && diff[i] != 0) || std::abs(dist[i] - u) > 1e-7)
            {                                               // not colinear
              t = std::numeric_limits<double>::quiet_NaN(); // no intersection
            }
        }
    }

  return (v < -1e-14) ? -1 : t;
}

inline void test_rayCollision()
{
  Vector3d p0 = Vector3d(0, 0, 0);
  Vector3d p1 = Vector3d(0, 1, 0);
  Vector3d p2 = Vector3d(1, 0, 0);
  Vector3d p3 = Vector3d(1, 1, 0);
  Vector3d p4 = Vector3d(1, -1, 0);
  Vector3d p5 = Vector3d(4, 0, 0);

  Vector3d North = p1;
  Vector3d NorthW = Vector3d(-1, 1, 0);
  Vector3d NorthE = p3;
  Vector3d NorthWW = Vector3d(-1, 0.5, 0);
  Vector3d NorthEE = Vector3d(1, 0.5, 0);
  Vector3d South = -North;
  Vector3d SouthW = -NorthE;
  Vector3d SouthWW = Vector3d(-1, -0.5, 0);
  Vector3d SouthEE = Vector3d(1, -0.5, 0);

  // case 1
  // originating from same point, 0 < angle < 90
  double test1 = compute_rayCollision(p0, North, p0, NorthE);
  assert(test1 == 0.0);

  // case 2
  // originating from different points
  // parallel
  double test2 = compute_rayCollision(p0, NorthE, p1, NorthE);
  assert("should return Nan" && test2 != test2);

  // case 3
  // originating from different points
  // collide on ray0
  double test3 = compute_rayCollision(p0, North, p2, NorthWW);
  assert(test3 == 0.5);

  // case 4
  // originating from different points
  // collide on ray1
  double test4 = compute_rayCollision(p0, NorthEE, p2, North);
  assert("will be 1, because the vectors are not unit legnth" && test4 == 1.0);
  assert("length should be about 1.118 times unit lenght of ray0" &&
         (NorthEE.norm() - 1.118) < 1e-3);

  // case 5
  // originating from different points
  // collide on negative ray0
  // direction vectors 0 < angle < 90
  double test5 = compute_rayCollision(p0, North, p4, NorthWW);
  assert(test5 < 0);

  // case 6
  // originating from different points
  // collide on negative ray1
  // direction vectors 90 < angle < 180
  double test6 = compute_rayCollision(p0, SouthEE, p1, North);
  assert(test6 < 0);

  // case 7
  // originating from different points
  // collinear rays in same direction
  double test7 = compute_rayCollision(p0, NorthE, p3, NorthE);
  assert(test7 == 1.0);

  // case 8
  // originating from same point
  // collinear rays in same direction
  double test8 = compute_rayCollision(p1, NorthE, p1, NorthE);
  assert(test8 == 0.0);

  // case 9
  // originating from different points
  // colinear with rays in opposite directions
  double test9 = compute_rayCollision(p0, NorthE, p3, SouthW);
  assert(test9 == 1.0);

  // case 10
  // rays facing each other, collide after the length of passed in ray (but on
  // ray)
  double test10 = compute_rayCollision(p0, NorthE, p5, NorthW);
  assert(test10 == 2.0);
}

inline void test_projectCone()
{
  Vector3d e1 = Vector3d(1, 0, 1);
  Vector3d e2 = Vector3d(1, 0, -1);

  Vector3d ray1 = Vector3d(0.5, 1, 0);
  Vector3d ray2 = Vector3d(-0.5, 1, -0.5);

  Vector3d tmp;

  // case 1
  // axis diagonal between x and z, on x-z plane
  // ray diagonal on x-y plane
  tmp = ray1;
  assert(projectCone(tmp, e1, e2));
  assert(tmp[0] == 0.5 && tmp[1] == 0);

  // case 2
  // ray outside of the same triangle as case 1
  tmp = ray2;
  assert(!projectCone(tmp, e1, e2));
}

// http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
inline double point_to_line_d_squared(Vector3d line_point_1,
                                      Vector3d line_point_2, Vector3d point)
{
  Vector3d numerator_vec =
      (line_point_2 - line_point_1).cross(line_point_1 - point);
  double numerator = numerator_vec.dot(numerator_vec);
  double denominator =
      (line_point_2 - line_point_1).dot(line_point_2 - line_point_1);
  return numerator / denominator;
}
#endif // UTILS_H