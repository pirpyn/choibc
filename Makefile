#!/usr/bin/make -f
SHELL:=/bin/bash
############################################################################
.PHONY: help
help:
	@echo "Compilation rules"
	@echo "   make                  # Basic compilation"
	@echo "   make CXXC=g++         # compiler: g++ accepted"
	@echo "   make MODE=optim       # optim/dev/debug provides different compilation option"
	@echo "   make SHARED=static    # static/shared for for dynamic or static libraries / binaries"
	@echo "   make all CXXCS=\"g++\" MODES=\"optim debug\" SHAREDS=\"static shared\" # Loop over each CXXC/MODE pair and compile with it"
	@echo "   make test             # compile the unit test"
	@echo ""
	@echo "Running rules"
	@echo "   make run CXXC=... MODE=... ARGS=\"arg1 arg2 ...\"    # Compile then execute the program with ARGS as argument"
	@echo "   make run_test"
	@echo ""
	@echo "Help rules"
	@echo "   make info             # Prints compiler flags, link flags, library used etc."
	@echo ""
	@echo "Cleaning rules"
	@echo "   make clean CXXC=g++ MODE=optim                       # Remove the build folder associated with the CXXC/MODE chosen" 
	@echo "   make clean_all CXXCS=\"g++ ...\" MODES=\"optim ...\" # Same for every CXXC/MODE pair possible" 
	@echo ""


############################################################################
# Sources directories to compile the hoibc library
SRCDIRS:=./src/main ./src/hoibc ./src/alglib ./src/bessel
LIBNAMES:=

EXTENSIONS:=.cpp .cc
SUFFIXES:=.d
# Binaries to build located in SRCDIRS
BINS:=main

# C++ compiler
CXXC:=g++

CXXFLAGS:=
DFLAGS:=-Wall
LDFLAGS:=-llapacke
INCFLAGS:=

#############################################################################
# Variables for buildings
.EXTENSIONS: $(EXTENSIONS) .o

# Processor architecture for separate building
release:=$(shell uname -r)
ifeq ($(CXXC),g++)
fcvers:=$(shell $(CXXC) --version 2>&1 | head -n 1 | cut -d ' ' -f 4)
else ifeq ($(CXXC),clang)
fcvers:=$(shell $(CXXC) --version 2>&1 | head -n 1 | cut -d ' ' -f 3 |cut -d '-' -f 1)
endif

blddir:=./build/$(release)/$(CXXC)/$(fcvers)
lklibdir:=./build/lib
lkbindir:=./build/bin
ifeq ($(DEBUG),1)
	MODE=debug
else
	MODE=optim
endif

MODE=optim
ifeq ($(MODE),debug)
	CXXFLAGS+=$(DFLAGS) -g -O0 -D_DEBUG
else ifeq ($(MODE),dev)
	CXXFLAGS+=-g -O2 -D_DEV
else ifeq ($(MODE),optim)
	CXXFLAGS+=-O2 -D_OPTIM
endif
blddir:=$(blddir)/$(MODE)

SHARED=0
ifeq ($(SHARED),1)
	libext=so
	CXXFLAGS+=-fPIC
	blddir:=$(blddir)/shared
	LDFLAGS+=-Wl,-R(libdir)
else
	libext=a
	blddir:=$(blddir)/static
endif

# Building subdirectories
objdir:=$(blddir)/obj
moddir:=$(blddir)/include
libdir:=$(blddir)/lib
bindir:=$(blddir)/bin

LDFLAGS:=-L$(libdir) $(foreach lib,$(LIBNAMES),$(patsubst %,-l%,$(lib))) $(LDFLAGS)
CXXFLAGS+=$(INCFLAGS)

dirs:=$(blddir) $(objdir) $(moddir) $(libdir) $(bindir) $(lklibdir) $(lkbindir)

# Static libraries to build
libs:=$(foreach dir,$(SRCDIRS),$(patsubst %,$(libdir)/lib%.$(libext),$(basename $(notdir $(dir)))))

# Static libraries to use
libs_depends:=$(patsubst %,$(libdir)/lib%.$(libext),$(LIBNAMES))

# Corresponding sources files
sources:=$(foreach srcdir,$(SRCDIRS),$(foreach ext,$(EXTENSIONS),$(wildcard $(srcdir)/*$(ext))))

# Corresponding objects
objects:=$(foreach file,$(sources),$(patsubst %,$(objdir)/%.o,$(basename $(notdir $(file)))))

# Binaries names to build
bins:=$(foreach bin,$(BINS),$(patsubst %,$(bindir)/%,$(bin)))

# Path to look for sources files
VPATH:=$(SRCDIRS) $(objdir) $(libdir)

#############################################################################
# Rules for 'hoibc'
.DEFAULT_GOAL=hoibc
.PHONY: hoibc
hoibc: info
	@echo "Compiling the HOIBC library"
	@$(MAKE) depend
	@$(MAKE) SRCDIRS=./src/bessel lib -j
	@$(MAKE) SRCDIRS=./src/alglib lib -j
	@$(MAKE) SRCDIRS=./src/hoibc LIBNAMES="alglib bessel" lib -j

.PHONY: main
main: hoibc
	@echo "Compiling the program to compute HOIBC coefficient"
	@$(MAKE) SRCDIRS=./src/main depend
	@$(MAKE) SRCDIRS=./src/main LIBNAMES="hoibc alglib bessel" lib -j
	@$(MAKE) SRCDIRS=./src/main LIBNAMES="main hoibc alglib bessel" prog -j

.PHONY: link
link: main | $(lklibdir) $(lkbindir)
	@echo "Simlinking binaries and libraries"
	@for lib in $(libs); do \
	    rm -rf $(lklibdir)/$$(basename $${lib}); \
	    ln -s $$(readlink -m $${lib}) $(lklibdir); \
	done;
	@echo "Libraries linked at $(lklibdir)"
	@for bin in $(bins); do \
	    rm -rf $(lkbindir)/$$(basename $${bin}); \
	    ln -s $$(readlink -m $${bin}) $(lkbindir); \
	done;
	@echo "Binaries linked at $(lkbindir)"

MODES=$(MODE)
CXXCS=$(CXXC)
SHAREDS=$(SHARED)
.PHONY: all
all:
	@for MODE in $(MODES); do \
	  for CXXC in $(CXXCS); do \
	    for SHARED in $(SHAREDS); do \
	      $(MAKE) CXXC=$${CXXC} MODE=$${MODE} SHARED=$${SHARED} main; \
	    done; \
	  done; \
	done

.PHONY:lib
lib: $(libs) | $(libdir)

.PHONY: prog
prog: $(bins) | $(bindir)

.PHONY: info
info:
	@echo "--------------------------------------------------------------------"
	@echo "CXXC:      $(CXXC) $(fcvers)"
	@echo "CXXFLAGS:  $(CXXFLAGS) $(INCFLAGS)"
	@echo "LDFLAGS:   $(LDFLAGS)"
	@echo "--------------------------------------------------------------------"

$(dirs):
	mkdir -p $@;

#############################################################################

$(objdir)/%.o: %.cpp %.hpp %.d | $(objdir)
	@echo "  $<"
	@$(CXXC) $(CXXFLAGS) -o $@ -c $<

$(objdir)/%.o: %.cc %.hh %.h %.d | $(objdir)
	@echo "  $<"
	@$(CXXC) $(CXXFLAGS) -o $@ -c $<

$(libdir)/%.a: $(objects)
	@echo "Creating $@"
	@ar crs $@ $(objects)

$(libdir)/%.so: $(objects)
	@echo "Creating $@"
	@$(CXXC) -shared -o $@ $^

$(bindir)/%: $(objdir)/%.o $(libs) $(libs_depends)
	@echo "Linking $@"
	@$(CXXC) -o $@ $< $(LDFLAGS)

#############################################################################
# Rules for dependencies
# https://stackoverflow.com/a/313787/8506658

# # Corresponding dependencies
dependencies:=$(objects:%.o=%.d)

#We don't need to clean up when we're making these targets
nodeps:=clean clean_all info

.PHONY: depend
depend: $(dependencies)
	@echo "Dependencies done for $(SRCDIRS)"

# Don't create dependencies when we're cleaning, for instance
ifeq (0, $(words $(findstring $(MAKECMDGOALS), $(nodeps))))
    #Chances are, these files don't exist.  GMake will create them and
    #clean up automatically afterwards
    -include $(dependencies)
endif

# This is the rule for creating the dependency files
$(objdir)/%.d: %.cpp | $(objdir)
	@echo "  $@"
	@$(CXXC) $(CXXFLAGS) -MM -MT '$(patsubst %.d,%.o,$@)' $< -MF $@

$(objdir)/%.d: %.cc | $(objdir)
	@echo "  $@"
	@$(CXXC) $(CXXFLAGS) -MM -MT '$(patsubst %.d,%.o,$@)' $< -MF $@


#############################################################################
TEST_SRCDIR:=./src/test
TEST_BINS:=$(foreach ext,$(EXTENSIONS),$(patsubst $(TEST_SRCDIR)/%$(ext),%,$(wildcard $(TEST_SRCDIR)/*$(ext))))

.PHONY: test
test: hoibc
	@echo "Compiling the unit tests"
	@$(MAKE) SRCDIRS="$(SRCDIRS) $(TEST_SRCDIR)" depend
	@$(MAKE) SRCDIRS=$(TEST_SRCDIR) libs="" BINS="$(TEST_BINS)" LIBNAMES="hoibc" prog

.PHONY: run_test
run_test: test
	@TESTS=( $(foreach bin,$(TEST_BINS),$(bindir)/$(bin)) ); \
	echo "Running the $${#TESTS[@]} tests"; \
	for ((i=0;i<$${#TESTS[@]};i++)); do \
		printf  "[%3d / %3d] " $$(($${i}+1)) $${#TESTS[@]}; \
		./$${TESTS[$$i]}; \
	done

#############################################################################
#############################################################################
# Other rules
.PHONY: clean
clean:
	$(RM) -rf $(blddir)

.PHONY: clean_all
clean_all:
	@for MODE in $(MODES); do \
	  for CXXC in $(CXXCS); do \
	    for SHARED in $(SHAREDS); do \
		    $(MAKE) CXXC=$${CXXC} SHARED=$${SHARED} MODE=$${MODE} clean; \
		  done; \
	  done; \
	done

.PHONY: run
run: main
	@for bin in $(bins); do \
		echo ; \
		echo ">> $(PREFIX) $${bin} $(ARGS)"; \
		echo ; \
		$(PREFIX) $${bin} $(ARGS); \
	done

# %.cpp: ;
%.hpp: ;
Makefile: ;