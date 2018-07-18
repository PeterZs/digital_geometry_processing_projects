# Run ``make depend'' first.
# For consistency use $() for commands and ${} for variables.

## ------------------------------------------------------
##              DIRECOTIRES
## ------------------------------------------------------

SRC=./src/
BIN=./bin/
INC=./include/
BUILD=./build/

## ------------------------------------------------------
##              FLAGS
## ------------------------------------------------------
CC = g++ 
LDFLAGS = -std=c++11
DEPEND= g++

ifeq (${METHOD},opt)
CFLAGS=-I${INC} -O3 -DNDEBUG -march=native
EXT=opt	
else
CFLAGS = -I${INC} -g -Wall -pedantic -Wextra -Wno-unused-parameter
EXT=dbg
endif

## ------------------------------------------------------
##              HEADERS, SOURCES, AND OBJECT FILES
## ------------------------------------------------------
SRCFILES =                \
	${SRC}kdtree.cpp      \
	${SRC}main.cpp     	  \
	${SRC}object.cpp   	  \
	${SRC}parser.cpp   	  \
	${SRC}raytracer.cpp	  
INCFILES =             	  \
	${INC}basic.hpp    	  \
	${INC}image.hpp    	  \
	${INC}object.hpp   	  \
	${INC}raytracer.hpp	  \
	${INC}GetPot       	  \
	${INC}kdtree.hpp   	  \
	${INC}parser.hpp   	  \
	${INC}scene.hpp		  
OBJFILES = $(patsubst ${SRC}%.cpp, ${BUILD}%.o_${EXT}, ${SRCFILES})

## ------------------------------------------------------
##              DEFAULT COMMAND 
## ------------------------------------------------------

all: ${BIN}raytrace.${EXT}

## ------------------------------------------------------
##              RULES FOR CREATING THE EXECUTABLE
## ------------------------------------------------------

${BIN}raytrace.${EXT}: ${OBJFILES}
	${CC} ${LDFLAGS} ${CFLAGS} ${OBJFILES} -o $@

## ------------------------------------------------------
##              RULES FOR CREATING THE OBJECT FILES
## ------------------------------------------------------

${BUILD}%.o_${EXT}:  ${SRC}%.cpp
	${CC} ${LDFLAGS} ${CFLAGS} -c $< -o $@


## ------------------------------------------------------
##                        PHONY COMMANDS
## ------------------------------------------------------
.PHONY: clean test echo depend

clean:
	rm -rf ${BUILD}* ${BIN}* *~ scenes/*.bmp

test:
	make  -j 3
	${BIN}raytrace.dbg scenes/cube.ray
	make -j 3 METHOD=opt
	${BIN}raytrace.opt scenes/basic.ray


echo:
	@echo ${SRCFILES}
	@echo ""
	@echo ${OBJFILES}
	@echo ""
	@echo ${INCFILES}

Makefile: ;

## ------------------------------------------------------
##                         MAKEDEPEND
## ------------------------------------------------------

# Old and garbage
# depend:
# 	 makedepend -- ${INC_DIR} -- -o.o_dbg -Y${SRC} -pbuild/ -f- > Makefile.depend
# 	 makedepend -- ${INC_DIR} -- -o.o_opt -Y${SRC} -pbuild/ -f- >> Makefile.depend

## Working very well.
depend:
	mkdir -p ${BIN}
	mkdir -p ${BUILD}
	rm -f ${BUILD}Makefile.depend_${EXT}
	$(foreach srcfile, ${SRCFILES}, ${DEPEND} ${CFLAGS} -MM ${srcfile} -MT $(patsubst %.cpp, ${BUILD}%.o_${EXT}, $(notdir ${srcfile})) >> ${BUILD}Makefile.depend_${EXT};)


-include ${BUILD}Makefile.depend_${EXT}
