/* Copyright (c) Russell Gillette
 * December 2013
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

#define GLFW_INCLUDE_GLU

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <iostream>

#include <math.h>
#include <stdio.h>

#ifndef GLM_FORCE_RADIANS
#define GLM_FORCE_RADIANS
#endif
#include <glm/ext.hpp>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <imgui.h>
#include <imgui_impl_glfw_gl3.h>

#include "ControlState.h"
#include "DrawMesh.h"
#include "EditMesh.h"
#include "MeshUtils.h"
#include "RenderState.h"
#include "ShaderUtils.h"
#include "TextureTypes.h"
#include "TextureUtils.h"
#include "Utils.h"
#include "WorldState.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "contrib/stb_image_write.h"

// === Globals ===
WorldState *w_state;
RenderState *r_state[2];
Model *g_mesh;
DrawMesh *g_axis; // NOTE: only a single axis

// === Mesh Files ===
int mesh_curr = -1;
int mesh_file_size = 8; // size of the array below
const char *mesh_files[] = {"Mesh/camel.obj",    "Mesh/camel_simple.obj",
                            "Mesh/cow_head.obj", "Mesh/cow1.obj",
                            "Mesh/cow2.obj",     "Mesh/cube.obj",
                            "Mesh/horse.obj",    "Mesh/octopus.obj"};

// main loop, does everything: poll events, update world, render
void mainloop()
{
  ControlState &c_state = ControlState::singleton();
  
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

  glfwPollEvents();

  static bool simplify_window_open = true;
  static int simplify_n_operations = 1;

  ImGui_ImplGlfwGL3_NewFrame();
  if (ImGui::BeginMainMenuBar())
    {
      if (ImGui::BeginMenu("File"))
        {
          if (ImGui::MenuItem("Load Next", "L"))
            {
              c_state.op = EDIT_LOAD_NEXT;
            }
          ImGui::Separator();
          for (int i = 0; i < mesh_file_size; ++i)
            {
              if (ImGui::MenuItem(mesh_files[i], nullptr, (i == mesh_curr)))
                {
                  c_state.op = EDIT_RELOAD;
                  mesh_curr = i;
                }
            }
          ImGui::Separator();
          if (ImGui::MenuItem("Save Image", "S"))
            {
              c_state.op = EDIT_SAVE_IMAGE;
            }
          if (ImGui::MenuItem("Quit", "Q"))
            {
              glfwSetWindowShouldClose(c_state.window, GL_TRUE);
            }
          ImGui::EndMenu();
        }

      if (ImGui::BeginMenu("View"))
        {
          if (ImGui::MenuItem("Next", "N"))
            {
              c_state.view_mode = c_state.view_mode + 1 > VIEW_ALL
                                      ? VIEW_FACES
                                      : c_state.view_mode + 1;
            }
          if (ImGui::MenuItem("Reset View", "R"))
            {
              c_state.clearViewDeltas();
            }
          ImGui::Separator();
          if (ImGui::MenuItem("Faces", "F",
                              (c_state.view_mode & VIEW_FACES) > 0))
            c_state.view_mode ^= VIEW_FACES;
          if (ImGui::MenuItem("Edges", "E",
                              (c_state.view_mode & VIEW_EDGES) > 0))
            c_state.view_mode ^= VIEW_EDGES;
          if (ImGui::MenuItem("Vertices", "V",
                              (c_state.view_mode & VIEW_VERTS) > 0))
            c_state.view_mode ^= VIEW_VERTS;
          if (ImGui::MenuItem("Axis", "A", c_state.view_axis))
            c_state.view_axis = !c_state.view_axis;
          ImGui::EndMenu();
        }

      // CS 524: Feel free to edit and add entries below, remember that you also
      // have to set keyboard events in ControlState.cpp (if you want them)
      if (ImGui::BeginMenu("Edit"))
        {
          if (ImGui::MenuItem("Clear Selection", "C"))
            {
              c_state.op = EDIT_CLEAR_SELECTION;
            }
          if (ImGui::MenuItem("Simplify", "", simplify_window_open))
            {
              simplify_window_open = !simplify_window_open;
            }
          if (ImGui::MenuItem("Item1", "1"))
            {
              printf("Item1 clicked\n"); /* c_state.op = EDIT_something; */
            }
          if (ImGui::MenuItem("Item2", "2"))
            {
              printf("Item2 clicked\n"); /* c_state.op = EDIT_something; */
            }
          if (ImGui::MenuItem("Item3", "3"))
            {
              printf("Item3 clicked\n"); /* c_state.op = EDIT_something; */
            }
          if (ImGui::MenuItem("Item4", "4"))
            {
              printf("Item4 clicked\n"); /* c_state.op = EDIT_something; */
            }
          if (ImGui::MenuItem("item5", "5"))
            {
              printf("Item5 clicked\n"); /* c_state.op = EDIT_something; */
            }
          ImGui::EndMenu();
        }

      // CS 524: You can also have buttons directly in the main bar
      if (ImGui::Button("Item6"))
        {
          printf("Item6 clicked\n");
        }
      ImGui::SameLine(ImGui::GetWindowWidth() - 80);
      ImGui::Text("(%.0f fps)", ImGui::GetIO().Framerate);
      ImGui::EndMainMenuBar();
    }

  // CS 524: Or you can have a separate window
  if (simplify_window_open)
    {
      ImGui::Begin("Mesh Simplification");
      // feel free to add more controls
      ImGui::SliderInt("Number operations", &simplify_n_operations, 1, 100);
      if (ImGui::Button("Vertex Removal Simplification"))
        {
        }
      if (ImGui::Button("Edge Collapse Simplification"))
        {
        }
      ImGui::End();
    }

  // Clear the buffer we will draw into.
  glClearColor(0.549f, 0.47f, 0.937f, 1.0f);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // Setup camera
  w_state->resetProjection(c_state.aspectRatio());
  w_state->updateView(c_state.viewTheta, c_state.viewPhi, c_state.viewDepth,
                      c_state.viewPan);

  /***********************************
   * Apply pending operation
   ***********************************/
  switch (c_state.op)
    {
    case EDIT_LOAD_NEXT:
      mesh_curr = (mesh_curr + 1) % mesh_file_size;
    // fallthrough
    case EDIT_RELOAD:
      delete g_mesh;
      g_mesh = loadModelFromFile(*r_state[0], mesh_files[mesh_curr]);
      c_state.op = EDIT_NONE; // reset the operation
      break;
    case EDIT_SAVE_IMAGE:
      break; // cannot parse it here, since the frame is not rendered yet
    case EDIT_CLEAR_SELECTION:
      g_mesh->get_editMesh()->deselect_allVerts();
      c_state.op = EDIT_NONE;
      break;

    case EDIT_SIMPLIFIY_VERTEX_REMOVAL:
      g_mesh->get_editMesh()->simplify_vertex_removal(simplify_n_operations);
      c_state.op = EDIT_NONE;
      break;

    case EDIT_SIMPLIFY_EDGE_COLLAPSE:
      g_mesh->get_editMesh()->simplify_edge_collapse(simplify_n_operations);
      c_state.op = EDIT_NONE;
      break;

    // CS 524: you can add entry points for your functions here.
    // Remember to set c_state.op either in ControlState (keyboard event
    // handling)
    // or in the GUI event handlers above, or both!

    case EDIT_NONE:
    default:
      break;
    }

  /***********************************
   * XYZ Axis Code
   ***********************************/
  if (c_state.view_axis)
    {
      w_state->useProgram(0);
      w_state->loadLights();

      // Draw X axis in red
      w_state->loadColorMaterial(glm::vec4(1, 0, 0, 1));
      w_state->loadObjectTransforms(glm::rotate(
          glm::mat4(), static_cast<float>(-M_PI_2), glm::vec3(0, 0, 1)));
      g_axis->drawMesh();

      // Draw Y axis in green
      w_state->loadColorMaterial(glm::vec4(0, 1, 0, 1));
      w_state->loadTransforms();
      g_axis->drawMesh();

      // Draw Z axis in blue
      w_state->loadColorMaterial(glm::vec4(0, 0, 1, 1));
      w_state->loadObjectTransforms(glm::rotate(
          glm::mat4(), static_cast<float>(M_PI_2), glm::vec3(1, 0, 0)));
      g_axis->drawMesh();
    }

  /*************************************
   * Update Selection
   *************************************/
  glm::vec3 s_bl, s_tr;
  c_state.getMouseSelection(s_bl, s_tr);
  if (c_state.select_dirty)
    {
      w_state->useProgram(2);

      // determine which vertices are in the selection box
      GLint select_viewport[4];
      glGetIntegerv(GL_VIEWPORT, select_viewport);
      glm::vec3 bl = glm::unProject(
          glm::vec3(s_bl.x, s_bl.y, 0), w_state->view * w_state->model,
          w_state->projection,
          glm::vec4(0.0, 0.0, select_viewport[2], select_viewport[3]));
      glm::vec3 bl_ray = glm::unProject(
          glm::vec3(s_bl.x, s_bl.y, 1), w_state->view * w_state->model,
          w_state->projection,
          glm::vec4(0.0, 0.0, select_viewport[2], select_viewport[3]));
      glm::vec3 br = glm::unProject(
          glm::vec3(s_tr.x, s_bl.y, 0), w_state->view * w_state->model,
          w_state->projection,
          glm::vec4(0.0, 0.0, select_viewport[2], select_viewport[3]));
      glm::vec3 br_ray = glm::unProject(
          glm::vec3(s_tr.x, s_bl.y, 1), w_state->view * w_state->model,
          w_state->projection,
          glm::vec4(0.0, 0.0, select_viewport[2], select_viewport[3]));
      glm::vec3 tr = glm::unProject(
          glm::vec3(s_tr.x, s_tr.y, 0), w_state->view * w_state->model,
          w_state->projection,
          glm::vec4(0.0, 0.0, select_viewport[2], select_viewport[3]));
      glm::vec3 tr_ray = glm::unProject(
          glm::vec3(s_tr.x, s_tr.y, 1), w_state->view * w_state->model,
          w_state->projection,
          glm::vec4(0.0, 0.0, select_viewport[2], select_viewport[3]));
      glm::vec3 tl = glm::unProject(
          glm::vec3(s_bl.x, s_tr.y, 0), w_state->view * w_state->model,
          w_state->projection,
          glm::vec4(0.0, 0.0, select_viewport[2], select_viewport[3]));
      glm::vec3 tl_ray = glm::unProject(
          glm::vec3(s_bl.x, s_tr.y, 1), w_state->view * w_state->model,
          w_state->projection,
          glm::vec4(0.0, 0.0, select_viewport[2], select_viewport[3]));

      int vert_size = g_mesh->info_sizev();
      for (int i = 0; i < vert_size; i++)
        {
          Eigen::Vector3d vert = g_mesh->info_vertex(i);
          if (vert_inside_select_box(
                  Eigen::Vector3d(bl.x, bl.y, bl.z),
                  Eigen::Vector3d(bl_ray.x, bl_ray.y, bl_ray.z),
                  Eigen::Vector3d(br.x, br.y, br.z),
                  Eigen::Vector3d(br_ray.x, br_ray.y, br_ray.z),
                  Eigen::Vector3d(tr.x, tr.y, tr.z),
                  Eigen::Vector3d(tr_ray.x, tr_ray.y, tr_ray.z),
                  Eigen::Vector3d(tl.x, tl.y, tl.z),
                  Eigen::Vector3d(tl_ray.x, tl_ray.y, tl_ray.z), vert))
            {
              if (c_state.modShft)
                {
                  g_mesh->deselect_vert(i);
                }
              else
                {
                  g_mesh->select_vert(i);
                }
            }
        }

      c_state.select_dirty = false;
    }

  /***********************************
   * Mesh Code
   ***********************************/
  w_state->useProgram(1);
  // load values into shader
  glm::vec2 screen(c_state.width, c_state.height);

  glUniform1i(glGetUniformLocation(w_state->getCurrentProgram(), "view_mode"),
              c_state.view_mode);
  glUniform2fv(glGetUniformLocation(w_state->getCurrentProgram(), "scale"), 1,
               glm::value_ptr(screen));
  w_state->loadTransforms();
  w_state->loadMaterials();
  w_state->loadLights();
  w_state->loadTextures();

  if (g_mesh)
    g_mesh->drawMesh();

  /*************************************
   * Draw Selection Box
   *************************************/
  if (c_state.select_active)
    {
      w_state->useProgram(2);

      glUniform3fv(
          glGetUniformLocation(w_state->getCurrentProgram(), "bot_left"), 1,
          glm::value_ptr(s_bl));
      glUniform3fv(
          glGetUniformLocation(w_state->getCurrentProgram(), "top_right"), 1,
          glm::value_ptr(s_tr));

      drawBox(s_bl.x, s_bl.y, s_tr.x, s_tr.y, c_state.width, c_state.height);
    }

  /*************************************
   * Draw GUI
   *************************************/
  ImGui::Render();

  // glfwSwapInterval(0); // uncomment to make the renderer not wait for v-sync,
  // pushing frames as quickly as possible, mostly for benchmarking (and fun).
  glfwSwapBuffers(c_state.window);

  /*************************************
   * Capture screenshot
   *************************************/
  if (c_state.op == EDIT_SAVE_IMAGE)
    {
      const unsigned int bytesPerPixel = 3; // RGB
      unsigned char *pixels =
          new unsigned char[bytesPerPixel * c_state.width * c_state.height];
      glReadPixels(0, 0, c_state.width, c_state.height, GL_RGB,
                   GL_UNSIGNED_BYTE, pixels);
      for (int y = 0; y < c_state.height / 2; ++y)
        {
          const int y2 = c_state.height - y - 1;
          for (int x = 0; x < c_state.width; ++x)
            {
              const int offset1 = bytesPerPixel * (x + y * c_state.width);
              const int offset2 = bytesPerPixel * (x + y2 * c_state.width);
              std::swap(pixels[offset1 + 0], pixels[offset2 + 0]);
              std::swap(pixels[offset1 + 1], pixels[offset2 + 1]);
              std::swap(pixels[offset1 + 2], pixels[offset2 + 2]);
            }
        }
      stbi_write_png("capture.png", c_state.width, c_state.height, 3, pixels,
                     0);
      delete[] pixels;
      c_state.op = EDIT_NONE;
    }
}

// setup
int main(int argc, char *argv[])
{
  EditMesh::test();

  GLenum err = 0;
  /*********************************************
   * GLFW SETUP
   *********************************************/
  err = glfwInit();
  if (!err)
    {
      fputs("Failed to load the GLFW library", stderr);
      exit(EXIT_FAILURE);
    }

// Without these functions glew will load the compatibility api,
// which might not support OpenGL 3.3. I have no idea what windows
// does.
#if defined(__APPLE__) || defined(LINUX)
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#endif

  /*********************************************
   * STATE SETUP (initialize gl context)
   *********************************************/
  // must be setup before glew so that a valid openGL
  // context exists (created with the window)

  w_state = new WorldState();
  ControlState::singleton().init(*w_state);

/*********************************************
 * GLEW SETUP
 *********************************************/
// GLEW is experimental with the core api. Maybe it is better to
// use some other api loader which hopefully is not as shitty.
#if defined(__APPLE__) || defined(LINUX)
  glewExperimental = GL_TRUE;
#endif
  err = glewInit();
  if (err != GLEW_OK)
    {
      fputs("Failed to initialize the GLEW library", stderr);
      exit(EXIT_FAILURE);
    }

  /*********************************************
   * STATE SETUP (construct render states)
   *********************************************/
  // must be setup after glew so that GL array
  // objects exist

  r_state[0] = new RenderState();
  r_state[1] = new RenderState();

  /*********************************************
   * SHADER SETUP
   *********************************************/
  // read default shaders from file
  GLuint shaderProgram[3] = {0};
  GLuint shaders[3] = {0};

  // create axis shader program
  buildShader(GL_VERTEX_SHADER, "shaders/axes.vs.glsl", shaders[0]);
  buildShader(GL_FRAGMENT_SHADER, "shaders/default.fs.glsl", shaders[1]);
  shaderProgram[0] = buildProgram(2, shaders);

  // create the shaders for the mesh
  buildShader(GL_VERTEX_SHADER, "shaders/mesh.vs.glsl", shaders[0]);
  buildShader(GL_GEOMETRY_SHADER, "shaders/wireframe.gs.glsl", shaders[2]);
  buildShader(GL_FRAGMENT_SHADER, "shaders/wireframe.fs.glsl", shaders[1]);
  shaderProgram[1] = buildProgram(3, shaders);

  // load shaders to render selection
  buildShader(GL_VERTEX_SHADER, "shaders/passthrough.vs.glsl", shaders[0]);
  buildShader(GL_FRAGMENT_SHADER, "shaders/select.fs.glsl", shaders[1]);
  shaderProgram[2] = buildProgram(2, shaders);

  // bind shader program
  w_state->setProgram(0, shaderProgram[0]);
  w_state->setProgram(1, shaderProgram[1]);
  w_state->setProgram(2, shaderProgram[2]);
  w_state->useProgram(0);

  // setup the transform matrices and uniform variables
  w_state->loadTransforms();
  w_state->loadLights();
  w_state->loadMaterials();

  /*********************************************
   * SETUP IMGUI
   *********************************************/
  ImGui_ImplGlfwGL3_Init(ControlState::singleton().window, false);

  /*********************************************
   * LOAD MESH
   *********************************************/
  // instruct the mainloop to load the mesh at the first iteration
  mesh_curr =
      2; // you can set the default model to load, or -1 for none at all.
  ControlState::singleton().op = EDIT_RELOAD;

  g_axis = createAxis(*r_state[1], 1);

  /*********************************************
   * SET GL STATE
   *********************************************/
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  /*********************************************
   * MAIN LOOP
   *********************************************/
  ControlState::printHelp();

  while (!glfwWindowShouldClose(ControlState::singleton().window))
    mainloop();

  /*********************************************
   * CLEAN UP
   *********************************************/
  delete g_mesh;
  delete g_axis;

  ImGui_ImplGlfwGL3_Shutdown();

  glfwTerminate();

  exit(EXIT_SUCCESS);
}
