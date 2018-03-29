#Available make targets : `make`, `make clean`, `make run`
#
#To compile on OSX make sure you have the following headers
#- Eigen in / usr / local / include / Eigen
#- glm in / usr / local / include / GLFW
#- glew headers in / usr / local / include / GL
#- glfw headers in / usr / local / include / GLFW
#
#Eigen and glm are header - only libraries,
#you can just download them
#
#For GLEW you need to download the source, make &&sudo make install
#To install GLFW you probably need cmake and then make &&make install

CC = g++
CFLAGS := -std=c++11

# OSX
ifeq ("$(shell uname)", "Darwin")
  LFLAGS := -lglfw -framework Cocoa -framework CoreVideo -framework OpenGL -framework IOKit -lGLEW
  INC_DIR :=  -I/usr/include/eigen3 -I./contrib/imgui
# LINUX
else
  INC_DIR :=  -isystem/usr/include/eigen3 -I./contrib/imgui 
  LFLAGS = $(shell pkg-config --static --libs glew glfw3 gl)
  CFLAGS := ${CFLAGS} -DLINUX
endif

BIN = Cartel
SRC := contrib/imgui/imgui_impl_glfw_gl3.cpp \
       contrib/imgui/imgui.cpp \
       contrib/imgui/imgui_draw.cpp  \
       Texture.cpp \
       ShaderUtils.cpp \
       ControlState.cpp \
       VBuffer.cpp \
       WorldState.cpp \
       EditMesh.cpp \
       Model.cpp \
       DrawMesh.cpp \
       IBuffer.cpp \
       main.cpp 

OBJ = $(SRC:%.cpp=build/%.o)

default: build
build: $(BIN) 

$(BIN): build/main.o $(OBJ)
	$(CC) -o $@ $^  $(LFLAGS) 

build/%.o: %.cpp
	@ mkdir -p $(@D)
	$(CC)  $(CFLAGS)  -c $< -o $@ $(INC_DIR)

clean:
	rm -rf *~ $(OBJ) $(BIN) build

run: build execute

execute:
	./$(BIN)

depend:
	makedepend -Y ${SRC} &> /dev/null

.PHONY: execute run clean 


# DO NOT DELETE

contrib/imgui/imgui_impl_glfw_gl3.o: contrib/imgui/imgui.h
contrib/imgui/imgui_impl_glfw_gl3.o: contrib/imgui/imconfig.h
contrib/imgui/imgui_impl_glfw_gl3.o: contrib/imgui/imgui_impl_glfw_gl3.h
contrib/imgui/imgui.o: contrib/imgui/imgui.h contrib/imgui/imconfig.h
contrib/imgui/imgui.o: contrib/imgui/imgui_internal.h
contrib/imgui/imgui.o: contrib/imgui/stb_textedit.h
contrib/imgui/imgui_draw.o: contrib/imgui/imgui.h contrib/imgui/imconfig.h
contrib/imgui/imgui_draw.o: contrib/imgui/imgui_internal.h
contrib/imgui/imgui_draw.o: contrib/imgui/stb_textedit.h
contrib/imgui/imgui_draw.o: contrib/imgui/stb_rect_pack.h
contrib/imgui/imgui_draw.o: contrib/imgui/stb_truetype.h
Texture.o: Texture.h
ShaderUtils.o: ShaderUtils.h Utils.h
ControlState.o: ControlState.h WorldState.h TextureTypes.h Texture.h
ControlState.o: contrib/imgui/imgui.h contrib/imgui/imconfig.h
ControlState.o: contrib/imgui/imgui_impl_glfw_gl3.h
VBuffer.o: VBuffer.h RenderState.h
WorldState.o: WorldState.h TextureTypes.h Texture.h
EditMesh.o: EditMesh.h
Model.o: Model.h DrawMesh.h IBuffer.h RenderState.h VBuffer.h EditMesh.h
DrawMesh.o: DrawMesh.h IBuffer.h RenderState.h VBuffer.h Utils.h
IBuffer.o: IBuffer.h RenderState.h
main.o: contrib/imgui/imgui.h contrib/imgui/imconfig.h
main.o: contrib/imgui/imgui_impl_glfw_gl3.h ControlState.h WorldState.h
main.o: TextureTypes.h Texture.h DrawMesh.h IBuffer.h RenderState.h VBuffer.h
main.o: EditMesh.h MeshUtils.h Model.h OBJLoader.h ShaderUtils.h
main.o: TextureUtils.h Utils.h contrib/stb_image_write.h
