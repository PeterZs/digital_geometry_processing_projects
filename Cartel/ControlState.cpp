/* Copyright (c) Russell Gillette
 * December 2013
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

#include "ControlState.h"

#include <algorithm>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>

#include <imgui.h>
#include <imgui_impl_glfw_gl3.h>

ControlState ControlState::m_singleton;

ControlState::~ControlState() { glfwDestroyWindow(window); }

int ControlState::init(WorldState &w)
{
  this->w = &w;

  width = 800;
  height = 600;
  /* As of right now we only have one window */
  window = glfwCreateWindow(width, height, "Cartel", NULL, NULL);
  if (!window)
    {
      glfwTerminate();
      fputs("failed to initialize window", stderr);
      return 1; // error
    }
  glfwMakeContextCurrent(window);

  // bind all callbacks
  glfwSetKeyCallback(window, key_callback);
  glfwSetCharCallback(window, char_callback);
  glfwSetFramebufferSizeCallback(window, reshape_callback);
  glfwSetCursorPosCallback(window, mousePos_callback);
  glfwSetCursorEnterCallback(window, mouseEnter_callback);
  glfwSetMouseButtonCallback(window, mouseBtn_callback);
  glfwSetScrollCallback(window, mouseScroll_callback);

  return 0;
}

int ControlState::deltaArrLR() { return arrR - arrL; }

int ControlState::deltaArrUD() { return arrD - arrU; }

void ControlState::clearViewDeltas()
{
  viewTheta = 0;
  viewPhi = 0;
  viewDepth = 0.2f;
  viewPan = glm::vec3(0, 0, 0);
}


/*****************************************************************************
 * Passive Callback functions
 *****************************************************************************/
void ControlState::printHelp()
{
  printf("==== Help ====\n"
         "___Control___\n"
         "m - switch between modes (view and select)\n\n"
         "___Command___\n"
         "q, esc - exit\n"
         "h - help (you have already figured this out)\n\n"
         "___View___\n"
         "left click and drag      - adjusts view\n"
         "shft+left click and drag - pans the view\n"
         "scroll wheel             - zoom\n");
}

// error callback for GLFW
void ControlState::error_callback(int error, const char *desc) { fputs(desc, stderr); }

// callback when window is resized
void ControlState::reshape_callback(GLFWwindow *window, int w, int h)
{
  singleton().height = h;
  singleton().width = w;

  glViewport(0, 0, (GLint)w, (GLint)h);
}

/*****************************************************************************
 * Active Callback functions
 *****************************************************************************/
void ControlState::char_callback(GLFWwindow *window, unsigned int c)
{
  if (ImGui::GetIO().WantTextInput)
    {
      ImGui_ImplGlfwGL3_CharCallback(window, c);
      return;
    }
}

// callback when a key is pressed
void ControlState::key_callback(GLFWwindow *window, int key, int scancode, int action,
				int mode)
{
  if (ImGui::GetIO().WantTextInput)
    {
      ImGui_ImplGlfwGL3_KeyCallback(window, key, scancode, action, mode);
      return;
    }

  switch (key)
    {
    case GLFW_KEY_LEFT:
      singleton().arrL = (action == GLFW_RELEASE) ? 0 : 1;
      break;
    case GLFW_KEY_RIGHT:
      singleton().arrR = (action == GLFW_RELEASE) ? 0 : 1;
      break;
    case GLFW_KEY_UP:
      singleton().arrU = (action == GLFW_RELEASE) ? 0 : 1;
      break;
    case GLFW_KEY_DOWN:
      singleton().arrD = (action == GLFW_RELEASE) ? 0 : 1;
      break;

    case GLFW_KEY_ESCAPE: // see also Q
    case GLFW_KEY_Q:
      glfwSetWindowShouldClose(window, GL_TRUE);
      break;

    case GLFW_KEY_LEFT_SHIFT:
    case GLFW_KEY_RIGHT_SHIFT:
      singleton().modShft = (action == GLFW_RELEASE) ? 0 : 1;
      break;

    case GLFW_KEY_H:
      printHelp();
      break;
    case GLFW_KEY_L:
      if (action == GLFW_RELEASE)
        singleton().op = EDIT_LOAD_NEXT;
      break;
    case GLFW_KEY_M:
      if (action == GLFW_RELEASE)
        singleton().mode = (RENDER_MODE)((singleton().mode + 1) % MODE_MAX);
      break;
    case GLFW_KEY_S:
      if (action == GLFW_RELEASE)
        singleton().op = EDIT_SAVE_IMAGE;
      break;

    case GLFW_KEY_N:
      if (action == GLFW_RELEASE)
        singleton().view_mode = singleton().view_mode + 1 > VIEW_ALL
	  ? VIEW_FACES
	  : singleton().view_mode + 1;
	  break;
    case GLFW_KEY_F:
      if (action == GLFW_RELEASE)
        singleton().view_mode ^= VIEW_FACES;
      break;
    case GLFW_KEY_E:
      if (action == GLFW_RELEASE)
        singleton().view_mode ^= VIEW_EDGES;
      break;
    case GLFW_KEY_V:
      if (action == GLFW_RELEASE)
        singleton().view_mode ^= VIEW_VERTS;
      break;
    case GLFW_KEY_A:
      if (action == GLFW_RELEASE)
        singleton().view_axis = !singleton().view_axis;
      break;
    case GLFW_KEY_R:
      if (action == GLFW_RELEASE)
        singleton().clearViewDeltas();
      break;

    case GLFW_KEY_C:
      if (action == GLFW_RELEASE)
        singleton().op = EDIT_CLEAR_SELECTION;
      break;
    case GLFW_KEY_D:
      if (action == GLFW_RELEASE)
        singleton().op = EDIT_DEBUG;
      break;

    case GLFW_KEY_0:
      if (action == GLFW_RELEASE)
        {
          singleton().op = EDIT_SQRT3_SUBDIV;
        }
      break;
    case GLFW_KEY_1:
      if (action == GLFW_RELEASE)
        {
        }
      break;
    case GLFW_KEY_2:
      if (action == GLFW_RELEASE)
        {
        }
      break;
    case GLFW_KEY_3:
      if (action == GLFW_RELEASE)
        {
        }
      break;
    case GLFW_KEY_4:
      if (action == GLFW_RELEASE)
        {
        }
      break;
    case GLFW_KEY_5:
      if (action == GLFW_RELEASE)
        {
        }
      break;
    }
}

// callback when a mouse button is pressed
void ControlState::mouseBtn_callback(GLFWwindow *win, int button, int action, int mod)
{
  if (ImGui::GetIO().WantCaptureMouse)
    return;

  /* TODO: any controls relative to pressing the mouse buttons goes here */

  if (button == GLFW_MOUSE_BUTTON_LEFT)
    {
      singleton().mouseBtnL = (action == GLFW_PRESS) ? 1 : 0;
    }
  else if (button == GLFW_MOUSE_BUTTON_RIGHT)
    {
      if (action == GLFW_PRESS)
        {
          singleton().mouseBtnR = true;
          singleton().select_active = true;
        }
      else
        {
          singleton().mouseBtnR = false;
          singleton().select_active = false;
          singleton().select_dirty = true;
        }
    }
  else if (button == GLFW_MOUSE_BUTTON_MIDDLE)
    singleton().mouseBtnC = (action == GLFW_PRESS) ? 1 : 0;

  if (action == GLFW_PRESS)
    singleton().mouseEvent = true;
}

// callback when the mouse is moved. This will be called
// ALOT keep it as light as possible!!
void ControlState::mousePos_callback(GLFWwindow *win, double x, double y)
{
  if (ImGui::GetIO().WantCaptureMouse)
    return;

  float fx = static_cast<float>(x);
  float fy = static_cast<float>(y);

  // screen Y coords are inverted.
  fy = singleton().height - fy;

  // currently used to update camera angles if mouse pressed
  if (singleton().mouseBtnL)
    {
      // Calculate change from last known mouse positon.
      float dx = fx - singleton().mouseX;
      float dy = fy - singleton().mouseY;

      if (singleton().modShft) // update screen pan
        {
          float perc = 100 / singleton().viewDepth;
          singleton().viewPan =
	    singleton().viewPan + glm::vec3(dx / perc, dy / perc, 0);
        }
      else // Update viewing angles.
        {
          singleton().viewTheta = std::fmod(
					singleton().viewTheta + glm::radians(360.0f) + dx / 100 / 2.0f,
					glm::radians(360.0f));
          singleton().viewPhi = std::min(
				     glm::radians(90.0f),
				     std::max(glm::radians(-90.0f), singleton().viewPhi - dy / 100));
        }
    }

  singleton().mouseX = fx;
  singleton().mouseY = fy;

  if (singleton().mouseBtnR && singleton().mouseEvent)
    {
      singleton().select_start = glm::vec3(fx, fy, 0);
      singleton().select_end = singleton().select_start;
      singleton().mouseEvent = false;
    }
  else if (singleton().mouseBtnR)
    singleton().select_end = glm::vec3(fx, fy, 0);
}

void ControlState::mouseScroll_callback(GLFWwindow *win, double x_offset,
					double y_offset)
{
  if (ImGui::GetIO().WantCaptureMouse)
    return;

  // since we would read from mouseScroll, set viewDepth
  // and then clear mouseScroll, I shall refrain from even
  // setting it for this usage
  singleton().viewDepth -= static_cast<float>(y_offset) / 30;
}

void ControlState::mouseEnter_callback(GLFWwindow *win, int entered)
{
  singleton().mouseInWindow = entered > 0;
}
