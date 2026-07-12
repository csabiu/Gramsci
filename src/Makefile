# for gfortran
#
BINDIR="../bin/"

FLAGSALWAYS = -DDOUBLE -ffree-line-length-none

#FLAGSOPT=
#-O3 -ffast-math
#FLAGSDEBUG= -g -fbounds-check -O0

#F90=gfortran -DGFORTRAN -fopenmp -Ofast -ffast-math -faggressive-loop-optimizations -ftree-loop-if-convert
F90=gfortran -DGFORTRAN -fopenmp -Ofast -ffast-math -march=native -faggressive-loop-optimizations -ftree-loop-if-convert

# macOS SDK path for linker (fixes -lSystem not found)
SDKPATH=$(shell xcrun --show-sdk-path 2>/dev/null)
LDFLAGS=$(if $(SDKPATH),-L$(SDKPATH)/usr/lib,)

# choose debugging or optimization
FLAGS= ${FLAGSMACHINE} ${FLAGSALWAYS} ${INC} ${LIBS} ${FLAGSOPT}

MY_DIR=`basename ${PWD}`
#------------------------------------

%.o: %.f90
	$(F90) $(FLAGS) -c $< -o $@

%.o: %.F90
	$(F90) $(FLAGS) -c $< -o $@

all:	gramsci

# Module dependency chain
MODULES = extension.o kdtree2.o node_module.o sorting_module.o \
          config_module.o io_module.o analytic_randoms_module.o \
          graph_utils_module.o \
          graph_module.o query_2pcf_module.o query_3pcf_module.o \
          query_4pcf_module.o

gramsci: $(MODULES) gramsci.o
	${F90} ${FLAGS} ${LDFLAGS} -o gramsci $(MODULES) gramsci.o && \
	mv $@ $(BINDIR)/gramsci

# Explicit dependency rules for module compilation order
config_module.o: kdtree2.o node_module.o extension.o
io_module.o: kdtree2.o config_module.o
analytic_randoms_module.o: kdtree2.o config_module.o
graph_utils_module.o: kdtree2.o config_module.o
graph_module.o: kdtree2.o node_module.o sorting_module.o config_module.o graph_utils_module.o
query_2pcf_module.o: kdtree2.o config_module.o analytic_randoms_module.o
query_3pcf_module.o: kdtree2.o config_module.o graph_utils_module.o analytic_randoms_module.o
query_4pcf_module.o: kdtree2.o config_module.o graph_utils_module.o analytic_randoms_module.o
gramsci.o: $(MODULES)

clean:
	/bin/rm -f *.o *.mod
