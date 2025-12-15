# -*- Mode: Makefile -*- 

### This makefile produces an executable for the NPZ to EBIS direct converter
ebase = npz_to_ebis_cpp

DIM=3
USE_EB=TRUE

## location of 'lib' directory
CHOMBO_HOME = ../../../Chombo/lib

##
## names of Chombo libraries needed by this program, in order of search.
##
LibNames = Workshop EBAMRElliptic EBAMRTools EBTools AMRElliptic AMRTools BoxTools

# relative paths to source code directories
base_dir = .
src_dirs =

# Explicitly set source files to avoid duplication
vpath %.cpp .

override VERBOSE = # #program doesnt handle -q option

include $(CHOMBO_HOME)/mk/Make.test 