SRC=./src/
BIN=./bin/
INC=./include/
BUILD=./build/

CGAL_FLAGS= -DCGAL_USE_GMP -DCGAL_USE_MPFR -isystem		\
/usr/include/x86_64-linux-gnu -fstack-protector-strong -Wformat	\
-Werror=format-security -D_FORTIFY_SOURCE=2 -frounding-math

CGAL_LIBS= -Wl,-Bsymbolic-functions -Wl,-z,relro -Wl,--as-needed	\
-rdynamic -lmpfr -lgmp							\
/usr/lib/x86_64-linux-gnu/libCGAL_Core.so.11.0.1			\
/usr/lib/x86_64-linux-gnu/libCGAL.so.11.0.1 -lboost_thread		\
-lboost_system -lpthread

CC = g++ 
LDFLAGS = -std=c++11 -DARMA_DONT_USE_WRAPPER $(CGAL_FLAGS)
INCFLAGS =  -I$(INC)  
LIBS= -lopenblas -llapack $(CGAL_LIBS)
ifeq ($(METHOD),opt)
CFLAGS=$(INCFLAGS) -O3 -DNDEBUG  
EXT=opt	
else
CFLAGS =$(INCFLAGS) -g -Wall -pedantic -Wextra -Wno-unused-parameter 
EXT=dbg
endif

SRCFILES_ALL = $(wildcard $(SRC)*.cxx)
SRCFILES_EXEC = $(wildcard $(SRC)*.exe.cxx)
SRCFILES=  $(filter-out $(SRCFILES_EXEC), $(SRCFILES_ALL))
OBJFILES = $(patsubst $(SRC)%.cxx, $(BUILD)%.o_$(EXT), $(SRCFILES))
INCFILES = $(wildcard $(INC)*.hxx)
EXEFILES = $(patsubst $(SRC)%.cxx, $(BIN)%.$(EXT), $(SRCFILES_EXEC))

all: remesh

tests:  $(BIN)tests.exe.$(EXT)
flip: $(BIN)flip.exe.$(EXT)
area: $(BIN)area.exe.$(EXT)
split: $(BIN)split.exe.$(EXT)
collapse: $(BIN)collapse.exe.$(EXT)
make_sphere: $(BIN)make_sphere.exe.$(EXT)
remesh: $(BIN)remesh.exe.$(EXT)
patch: $(BIN)patch.exe.$(EXT)

$(BIN)%.exe.$(EXT): $(BUILD)%.exe.o_$(EXT) $(OBJFILES)
	$(CC) $(LDFLAGS) $(CFLAGS) $(OBJFILES) $<  $(LIBS) -o $@ 

$(BUILD)%.o_$(EXT): $(SRC)%.cxx $(INCFILES)
	$(CC) -c $(LDFLAGS) $(CFLAGS) $< -o $@

## ------------------------------------------------------
.PHONY: clean test echo clean-vtk

.PRECIOUS: $(OBJFILES)

clean:
	rm -rf $(BUILD)* $(BIN)* *vtk *~ \#*\#

clean-vtk:
	rm -rf *vtk

test:

echo:
	@echo $(EXEFILES)
	@echo ""
	@echo $(SRCFILES_EXEC)
	@echo ""
	@echo $(SRCFILES)
	@echo ""
	@echo $(OBJFILES)
	@echo ""
	@echo $(INCFILES)

Makefile: ;
