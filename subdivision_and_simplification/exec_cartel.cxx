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

#include "Cartel/contrib/imgui/imgui.h"
#include "Cartel/contrib/imgui/imgui_impl_glfw_gl3.h"
#include "Cartel/ControlState.h"
#include "Cartel/DrawMesh.h"
#include "Cartel/MeshUtils.h"
#include "Cartel/RenderState.h"
#include "Cartel/ShaderUtils.h"
#include "Cartel/TextureTypes.h"
#include "Cartel/TextureUtils.h"
#include "Cartel/Utils.h"
#include "Cartel/WorldState.h"

#include "edit_mesh.hxx"
#include "mesh_io.hxx"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "Cartel/contrib/stb_image_write.h"

// === Globals ===
WorldState *w_state;
RenderState *r_state[2];
hooshi::EditMesh *g_mesh;
int draw_id = 0;
DrawMesh *g_viz_mesh;
DrawMesh *g_axis; // NOTE: only a single axis
const char * mesh_file;

void update_viz_mesh()
{
  if( draw_id != g_mesh->get_edit_count() )
    {
      draw_id = g_mesh->get_edit_count();
      int i_size = 3 * g_mesh->get_face_size();
      int v_size = 3 * i_size; // g_mesh->get_vert_size();

      float *v_data = new float[v_size * 2]; // vertex and normal data
      int *i_data = new int[i_size];         // face indexes
      int *s_data = new int[i_size];         // selection data
      float *c_data = new float[v_size];     // vertex colors

      g_mesh->get_draw_data(v_data, i_data);
      g_mesh->get_draw_normals(&v_data[v_size]);
      g_mesh->get_draw_selection(s_data);
      // g_mesh->get_draw_colors(c_data);
      std::fill(c_data, c_data+v_size, -1.0f); 

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

      g_viz_mesh->loadVBuffer(0, sizeof(float) * 2 * v_size, (GLubyte *)v_data, 0,
			DRAW_MESH_NUM_CORE_ATTR, attr_info);

      attr_info[0].attrib_number = 2; // vertex selection is 2
      attr_info[0].attrib_size = sizeof(int);
      attr_info[0].data_offset = 0; // data starts at the beginning of array
      attr_info[0].data_stride = 0; // data is tightly packed
      attr_info[0].num_comp = 1;    // there is 1 components per vertex selection

      g_viz_mesh->loadVBuffer(1, sizeof(int) * i_size, (GLubyte *)s_data, 0,
			DRAW_MESH_NUM_EXTRA_ATTR, attr_info);

      attr_info[0].attrib_number = 3; // color is 3
      attr_info[0].attrib_size = sizeof(float);
      attr_info[0].data_offset = 0; // data starts at the beginning of array
      attr_info[0].data_stride = 0; // data is tightly packed
      attr_info[0].num_comp = 3;    // there are 3 components (r,g,b) per vertex

      g_viz_mesh->loadVBuffer(2, sizeof(float) * v_size, (GLubyte *)c_data, 0, 1,
			attr_info);

      g_viz_mesh->loadIBuffer(i_size, sizeof(int), i_data);

#undef DRAW_MESH_NUM_ATTR
      delete[] v_data;
      delete[] i_data;
      delete[] s_data;
      delete[] c_data;
    }

}

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
          if (ImGui::MenuItem("Subdivide sqrt(3)", ""))
            {
	      g_mesh->subdivide_sqrt3();
              g_mesh->flag_edited();
            }
          if (ImGui::MenuItem("Simplify", "", simplify_window_open))
            {
              simplify_window_open = !simplify_window_open;
	      if(!simplify_window_open) g_mesh->finalize_simplification();
            }
          ImGui::EndMenu();
        }

      ImGui::SameLine(ImGui::GetWindowWidth() - 80);
      ImGui::Text("(%.0f fps)", ImGui::GetIO().Framerate);
      ImGui::EndMainMenuBar();
    }

  if (simplify_window_open && g_mesh)
    {
      ImGui::Begin("Mesh Simplification");
      ImGui::SliderInt("Number operations", &simplify_n_operations, 1, 5000);

      char b1_text[500], b2_text[500];
      if(g_mesh->is_simplification_in_progress() == 0 )
	{
	  sprintf(b1_text, "Vertex Removal Simplification");
	  sprintf(b2_text, "Edge Collapse Simplification");
	}
      else if (g_mesh->is_simplification_in_progress() == 1 )
	{
	  sprintf(b1_text, "Vertex Removal Simplification --- Active");
	  sprintf(b2_text, "Edge Collapse Simplification");
	}
      else if (g_mesh->is_simplification_in_progress() == 2 )
	{
	  sprintf(b1_text, "Vertex Removal Simplification");
	  sprintf(b2_text, "Edge Collapse Simplification --- Active");
	}	
      if (ImGui::Button(b1_text))
        {
	  if( g_mesh->is_simplification_in_progress() == 0 )
	    {
	      g_mesh->init_simplification(1);
	    }
	  else if( g_mesh->is_simplification_in_progress() == 1 )
	    {
	      g_mesh->finalize_simplification();
	    }
        }
      if (ImGui::Button(b2_text))
        {
	  if( g_mesh->is_simplification_in_progress() == 0 )
	    {
	      g_mesh->init_simplification(2);
	    }
	  else if( g_mesh->is_simplification_in_progress() == 2 )
	    {
	      g_mesh->finalize_simplification();
	    }
        }
      if(  g_mesh->is_simplification_in_progress() )
	{
	  if (ImGui::Button("Simplify"))
	    {
	      for (int i = 0 ; i < simplify_n_operations ; ++i)
		{
		  if(!g_mesh->simplify()) break;;
		}
	      g_mesh->flag_edited();
	    }
	  if (ImGui::Button("Unsimplify"))
	    {
	      for (int i = 0 ; i < simplify_n_operations ; ++i)
		{
		  if(!g_mesh->restore_last_simplification_step()) break;
		}
	      g_mesh->flag_edited();
	    }
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
    case EDIT_RELOAD:
      printf("kdajflajdljf\n");
      if(g_mesh) delete g_mesh;
      if(g_viz_mesh) delete g_viz_mesh;
      //
      g_mesh = new hooshi::EditMesh();
      hooshi::MeshIO(*g_mesh).read_auto(mesh_file);
      g_mesh->flag_edited();
      g_viz_mesh = new DrawMesh(*r_state[0]);
      g_viz_mesh->init(1);
      //
      draw_id = -1;
      c_state.op = EDIT_NONE; // reset the operation
      break;
    case EDIT_SAVE_IMAGE:
      break; // cannot parse it here, since the frame is not rendered yet
    case EDIT_CLEAR_SELECTION:
      // Don't do anything
      // g_mesh->get_editMesh()->deselect_allVerts();
      c_state.op = EDIT_NONE;
      break;

    case EDIT_SIMPLIFIY_VERTEX_REMOVAL:
      //g_mesh->simplify_vertex_removal(simplify_n_operations);
      //c_state.op = EDIT_NONE;
      break;

    case EDIT_SIMPLIFY_EDGE_COLLAPSE:
      //g_mesh->simplify_edge_collapse(simplify_n_operations);
      //c_state.op = EDIT_NONE;
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

  if (g_viz_mesh)
    {
      update_viz_mesh();
      g_viz_mesh->drawMesh();
    }

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
  if( argc < 2 )
    {
      printf("%s [Mesh_file]\n", argv[0]);
      exit(1);
    }
  mesh_file = argv[1];
  g_mesh = NULL;
  g_viz_mesh = NULL;
    
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
  buildShader(GL_VERTEX_SHADER, "Cartel/shaders/axes.vs.glsl", shaders[0]);
  buildShader(GL_FRAGMENT_SHADER, "Cartel/shaders/default.fs.glsl", shaders[1]);
  shaderProgram[0] = buildProgram(2, shaders);

  // create the shaders for the mesh
  buildShader(GL_VERTEX_SHADER, "Cartel/shaders/mesh.vs.glsl", shaders[0]);
  buildShader(GL_GEOMETRY_SHADER, "Cartel/shaders/wireframe.gs.glsl", shaders[2]);
  buildShader(GL_FRAGMENT_SHADER, "Cartel/shaders/wireframe.fs.glsl", shaders[1]);
  shaderProgram[1] = buildProgram(3, shaders);

  // load shaders to render selection
  buildShader(GL_VERTEX_SHADER, "Cartel/shaders/passthrough.vs.glsl", shaders[0]);
  buildShader(GL_FRAGMENT_SHADER, "Cartel/shaders/select.fs.glsl", shaders[1]);
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
  if(g_mesh) delete g_mesh;
  if(g_viz_mesh) delete g_viz_mesh;
  
  delete g_axis;

  ImGui_ImplGlfwGL3_Shutdown();

  glfwTerminate();

  exit(EXIT_SUCCESS);
}
