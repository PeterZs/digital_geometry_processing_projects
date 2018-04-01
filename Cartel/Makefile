# Available make targets : `make`, `make clean`, `make run`
#
# To compile on OSX or LINUX you need
#- Eigen 
#- glm 
#- glew 
#- glfw 

CC = g++
CFLAGS := -std=c++11 -Wall
CFLAGS_DBG := -g -O0 -DDEBUG -Wpedantic
CFLAGS_OPT := -O3 -mtune=native -DNDEBUG -Wpedantic
EXT := _opt


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

ifeq ("${METHOD}", "dbg")
  CFLAGS := ${CFLAGS} ${CFLAGS_DBG}
  EXT := _dbg
else
  CFLAGS := ${CFLAGS} ${CFLAGS_DBG}
  EXT := _opt
endif

BIN = Cartel${EXT}
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
       IBuffer.cpp 

OBJ = $(SRC:%.cpp=build/%.o${EXT})

default: build

build: ${BIN} 

execute:
	env LIBGL_ALWAYS_SOFTWARE=1 ./${BIN}

run: build execute


${BIN}: build/main.o${EXT} ${OBJ}
	${CC} -o $@ $^  ${LFLAGS}

build/%.o${EXT}: %.cpp
	@ mkdir -p $(@D)
	${CC}  ${CFLAGS} -c $< -o $@ ${INC_DIR}

clean:
	rm -rf *~ ${OBJ}  ${BIN} build

Cartel${EXT}.a: ${OBJ}
	ar rcs $@ ${OBJ} 

lib:  Cartel${EXT}.a


.PHONY: execute run clean depend

# Comment the following lines out if they produce errors.  What they
# do is that they allow you to run make dpened This will let make
# recognize the dependencies so that you don't have to run make clean
# everytime.  They need having makedepend.
#
# This won't work on MAC. There is the alternative of g++ -MM on MAC,
# but it is slightly more difficult to use.
depend:
	 makedepend -- ${INC_DIR} -- -o.o_dbg -Y${SRC} -pbuild/ -f- > Makefile.depend
	 makedepend -- ${INC_DIR} -- -o.o_opt -Y${SRC} -pbuild/ -f- >> Makefile.depend


include Makefile.depend

