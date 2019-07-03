#############################################################################
## Building utilities useful to LowerEmbargo.
## To customize this makefile for your use, just modify the path to Trilinos'
## include Makefile.export file below.


include /home/kdeweese/src/trilinos-12.10.1-Source/trilinos_build/Makefile.export.Trilinos_install
#############################################################################

ifndef DCC
DCC=$(Trilinos_CXX_COMPILER)
endif

ifndef DEPS
DEPS=-M
endif

include $(ALLOBJFILES:.o=.d)
CXX_COMPILER_FLAGS = $(Trilinos_CXX_COMPILER_FLAGS)


#############################################################################
BELSRC = belos.cpp
BELOBJ = $(BELSRC:.cpp=.o)

BELEVSRC = belos_kosz_evolve.cpp
BELEVOBJ = $(BELEVSRC:.cpp=.o)

BELEVPSRC = belos_kosz_evolve_print.cpp
BELEVPOBJ = $(BELEVPSRC:.cpp=.o)

BELEVP2SRC = belos_kosz_evolve_print2.cpp
BELEVP2OBJ = $(BELEVP2SRC:.cpp=.o)

BELEVRSRC = belos_kosz_evolve_ratio.cpp
BELEVROBJ = $(BELEVRSRC:.cpp=.o)

BELCHKSRC = belos_kosz_check.cpp
BELCHKOBJ = $(BELCHKSRC:.cpp=.o)

BELCHKOSRC = belos_kosz_checkone.cpp
BELCHKOOBJ = $(BELCHKOSRC:.cpp=.o)

BELMUTSRC = belos_kosz_mutate.cpp
BELMUTOBJ = $(BELMUTSRC:.cpp=.o)


ALLOBJFILES = $(BELSRC) $(BELEVSRS) $(BELEVRSRS) $(BELCHKSRC) $(BELCHKOSRC) $(BELMUTSRC)
#############################################################################


%.o:    %.cpp
	@echo "Compiling $<..."
	$(Trilinos_CXX_COMPILER) $(CXX_COMPILER_FLAGS) \
	-c $<  \
	$(Trilinos_INCLUDE_DIRS) $(Trilinos_TPL_INCLUDE_DIRS) 

%.d:    %.cpp
	@echo "(Re)Building dependency for $<..."
	$(SHELL) -ec '$(DCC) $(DEPS) $(DEFS) $< \
	$(Trilinos_INCLUDE_DIRS) $(Trilinos_TPL_INCLUDE_DIRS) \
	| sed '\''s/$*\.o/& $@/g'\'' > $@'

help:
	@echo "Usage:"
	@echo "  make all"
	@echo "  make belos"
	@echo "  make belos_kosz_evolve"
	@echo "  make belos_kosz_evolve_ratio"
	@echo "  make beloscheck"
	@echo "  make clean"

all:	belos extract_large_cc testConvertMatrix test2d belosdbg beloscheck

belos: $(BELOBJ:.o=.d) $(BELSRC:.cpp=.o) $(BELOBJ)
	@echo "Linking belos.exe..."
	$(Trilinos_CXX_COMPILER) $(CXX_COMPILER_FLAGS) \
	$(BELOBJ) -o belos.exe \
	$(Trilinos_INCLUDE_DIRS) \
	$(Trilinos_LIBRARY_DIRS) $(Trilinos_LIBRARIES) \
	$(Trilinos_TPL_INCLUDE_DIRS) $(Trilinos_TPL_LIBRARIES) \
	$(Trilinos_EXTRA_LD_FLAGS)


belos_kosz_evolve: $(BELEVOBJ:.o=.d) $(BELEVSRC:.cpp=.o) $(BELEVOBJ)
	@echo "Linking belos_kosz_evolve.exe..."
	$(Trilinos_CXX_COMPILER) $(CXX_COMPILER_FLAGS) \
	$(BELEVOBJ) -o belos_kosz_evolve.exe \
	$(Trilinos_INCLUDE_DIRS) \
	$(Trilinos_LIBRARY_DIRS) $(Trilinos_LIBRARIES) \
	$(Trilinos_TPL_INCLUDE_DIRS) $(Trilinos_TPL_LIBRARIES) \
	$(Trilinos_EXTRA_LD_FLAGS) -lstdc++fs

belos_kosz_evolve_print: $(BELEVPOBJ:.o=.d) $(BELEVPSRC:.cpp=.o) $(BELEVPOBJ)
	@echo "Linking belos_kosz_evolve_print.exe..."
	$(Trilinos_CXX_COMPILER) $(CXX_COMPILER_FLAGS) \
	$(BELEVPOBJ) -o belos_kosz_evolve_print.exe \
	$(Trilinos_INCLUDE_DIRS) \
	$(Trilinos_LIBRARY_DIRS) $(Trilinos_LIBRARIES) \
	$(Trilinos_TPL_INCLUDE_DIRS) $(Trilinos_TPL_LIBRARIES) \
	$(Trilinos_EXTRA_LD_FLAGS) -lstdc++fs

belos_kosz_evolve_print2: $(BELEVP2OBJ:.o=.d) $(BELEVP2SRC:.cpp=.o) $(BELEVP2OBJ)
	@echo "Linking belos_kosz_evolve_print2.exe..."
	$(Trilinos_CXX_COMPILER) $(CXX_COMPILER_FLAGS) \
	$(BELEVP2OBJ) -o belos_kosz_evolve_print2.exe \
	$(Trilinos_INCLUDE_DIRS) \
	$(Trilinos_LIBRARY_DIRS) $(Trilinos_LIBRARIES) \
	$(Trilinos_TPL_INCLUDE_DIRS) $(Trilinos_TPL_LIBRARIES) \
	$(Trilinos_EXTRA_LD_FLAGS) -lstdc++fs


belos_kosz_evolve_ratio: $(BELEVROBJ:.o=.d) $(BELEVRSRC:.cpp=.o) $(BELEVROBJ)
	@echo "Linking belos_kosz_evolve_ratio.exe..."
	$(Trilinos_CXX_COMPILER) $(CXX_COMPILER_FLAGS) \
	$(BELEVROBJ) -o belos_kosz_evolve_ratio.exe \
	$(Trilinos_INCLUDE_DIRS) \
	$(Trilinos_LIBRARY_DIRS) $(Trilinos_LIBRARIES) \
	$(Trilinos_TPL_INCLUDE_DIRS) $(Trilinos_TPL_LIBRARIES) \
	$(Trilinos_EXTRA_LD_FLAGS) -lstdc++fs

belos_kosz_mutate: $(BELMUTOBJ:.o=.d) $(BELMUTSRC:.cpp=.o) $(BELMUTOBJ)
	@echo "Linking belos_kosz_mutate.exe..."
	$(Trilinos_CXX_COMPILER) $(CXX_COMPILER_FLAGS) \
	$(BELMUTOBJ) -o belos_kosz_mutate.exe \
	$(Trilinos_INCLUDE_DIRS) \
	$(Trilinos_LIBRARY_DIRS) $(Trilinos_LIBRARIES) \
	$(Trilinos_TPL_INCLUDE_DIRS) $(Trilinos_TPL_LIBRARIES) \
	$(Trilinos_EXTRA_LD_FLAGS) -lstdc++fs


beloscheck: $(BELCHKOBJ:.o=.d) $(BELCHKSRC:.cpp=.o) $(BELCHKOBJ)
	@echo "Linking beloscheck.exe..."
	$(Trilinos_CXX_COMPILER) $(CXX_COMPILER_FLAGS) \
	$(BELCHKOBJ) -o beloscheck.exe \
	$(Trilinos_INCLUDE_DIRS) \
	$(Trilinos_LIBRARY_DIRS) $(Trilinos_LIBRARIES) \
	$(Trilinos_TPL_INCLUDE_DIRS) $(Trilinos_TPL_LIBRARIES) \
	$(Trilinos_EXTRA_LD_FLAGS) -lstdc++fs

beloscheckone: $(BELCHKOOBJ:.o=.d) $(BELCHKOSRC:.cpp=.o) $(BELCHKOOBJ)
	@echo "Linking beloscheckone.exe..."
	$(Trilinos_CXX_COMPILER) $(CXX_COMPILER_FLAGS) \
	$(BELCHKOOBJ) -o beloscheckone.exe \
	$(Trilinos_INCLUDE_DIRS) \
	$(Trilinos_LIBRARY_DIRS) $(Trilinos_LIBRARIES) \
	$(Trilinos_TPL_INCLUDE_DIRS) $(Trilinos_TPL_LIBRARIES) \
	$(Trilinos_EXTRA_LD_FLAGS) -lstdc++fs


clean:	
	/bin/rm *.o *.d *.exe
