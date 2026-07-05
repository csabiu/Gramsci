# ===========================================================================
# OpenCL (portable GPU) build of GRAMSCI — gfortran + the OpenCL runtime.
#
# Works on:
#   - macOS / Apple Silicon & Intel  (Apple's OpenCL.framework; fp32 only)
#   - Linux with any OpenCL ICD       (NVIDIA / AMD / Intel / POCL)
#
# Usage:
#   make            # build bin/gramsci_cl
#   make clean
#
# Requirements: gfortran (Homebrew GCC on macOS) and an OpenCL runtime.
# The CPU modules from ../src are recompiled here so the build is self-contained.
# ===========================================================================
BINDIR = ../bin
SRCDIR = ../src

F90 = gfortran

# -march=native is the only CPU-specific flag.  A brand-new Apple core (e.g. a
# newer M-series) paired with a lagging GCC may not recognise it, so probe the
# compiler and drop the flag if it is rejected.  Keeps the build portable across
# M1/M2/M3/M4 and Intel Macs without sacrificing tuning where it is supported.
MARCH := $(shell printf 'program t\nend\n' | $(F90) -march=native -x f95 -fsyntax-only - >/dev/null 2>&1 && echo -march=native)

FFLAGS = -DGFORTRAN -DDOUBLE -fopenmp -Ofast -ffast-math $(MARCH) \
         -faggressive-loop-optimizations -ftree-loop-if-convert -ffree-line-length-none

# ---- Platform-specific OpenCL linkage ----
UNAME := $(shell uname)
ifeq ($(UNAME),Darwin)
  SDKPATH   := $(shell xcrun --show-sdk-path 2>/dev/null)
  LDFLAGS   := $(if $(SDKPATH),-L$(SDKPATH)/usr/lib,)
  # gfortran's linker needs the SDK framework search path (-F) to find OpenCL.
  OPENCL_LIB = $(if $(SDKPATH),-F$(SDKPATH)/System/Library/Frameworks,) -framework OpenCL
else
  LDFLAGS   :=
  OPENCL_LIB = -lOpenCL
endif

# ---- Object lists ----
# CPU modules (compiled from ../src with the same gfortran -> .mod compatible)
CPU_OBJS = extension.o kdtree2.o node_module.o sorting_module.o \
           config_module.o io_module.o analytic_randoms_module.o \
           graph_utils_module.o \
           graph_module.o query_2pcf_module.o \
           query_3pcf_module.o query_4pcf_module.o

# OpenCL backend modules (in src_opencl/)
CL_OBJS  = cl_module.o cl_kernels_module.o cl_env_module.o csr_cl_module.o \
           query_3pcf_cl_module.o query_4pcf_cl_module.o

ALL_OBJS = $(CPU_OBJS) $(CL_OBJS) gramsci_cl.o

all: gramsci_cl

# ---- Kernel embedding: kernels.cl -> cl_kernels_module.F90 ----
cl_kernels_module.F90: kernels.cl embed_kernels.sh
	sh embed_kernels.sh kernels.cl cl_kernels_module.F90

# ---- OpenCL backend source files (src_opencl/) ----
cl_module.o: cl_module.F90
	$(F90) $(FFLAGS) -c $< -o $@
cl_kernels_module.o: cl_kernels_module.F90
	$(F90) $(FFLAGS) -c $< -o $@
cl_env_module.o: cl_env_module.F90
	$(F90) $(FFLAGS) -c $< -o $@
csr_cl_module.o: csr_cl_module.F90
	$(F90) $(FFLAGS) -c $< -o $@
query_3pcf_cl_module.o: query_3pcf_cl_module.F90
	$(F90) $(FFLAGS) -c $< -o $@
query_4pcf_cl_module.o: query_4pcf_cl_module.F90
	$(F90) $(FFLAGS) -c $< -o $@
gramsci_cl.o: gramsci_cl.F90
	$(F90) $(FFLAGS) -c $< -o $@

# ---- CPU source files (from ../src) ----
extension.o: $(SRCDIR)/extension.F90
	$(F90) $(FFLAGS) -c $< -o $@
kdtree2.o: $(SRCDIR)/kdtree2.F90
	$(F90) $(FFLAGS) -c $< -o $@
node_module.o: $(SRCDIR)/node_module.F90
	$(F90) $(FFLAGS) -c $< -o $@
sorting_module.o: $(SRCDIR)/sorting_module.F90
	$(F90) $(FFLAGS) -c $< -o $@
config_module.o: $(SRCDIR)/config_module.F90
	$(F90) $(FFLAGS) -c $< -o $@
io_module.o: $(SRCDIR)/io_module.F90
	$(F90) $(FFLAGS) -c $< -o $@

analytic_randoms_module.o: $(SRCDIR)/analytic_randoms_module.F90
	$(F90) $(FFLAGS) -c $< -o $@
graph_utils_module.o: $(SRCDIR)/graph_utils_module.F90
	$(F90) $(FFLAGS) -c $< -o $@
graph_module.o: $(SRCDIR)/graph_module.F90
	$(F90) $(FFLAGS) -c $< -o $@
query_2pcf_module.o: $(SRCDIR)/query_2pcf_module.F90
	$(F90) $(FFLAGS) -c $< -o $@
query_3pcf_module.o: $(SRCDIR)/query_3pcf_module.F90
	$(F90) $(FFLAGS) -c $< -o $@
query_4pcf_module.o: $(SRCDIR)/query_4pcf_module.F90
	$(F90) $(FFLAGS) -c $< -o $@

gramsci_cl: $(ALL_OBJS)
	$(F90) $(FFLAGS) $(LDFLAGS) -o gramsci_cl $(ALL_OBJS) $(OPENCL_LIB)
	cp gramsci_cl $(BINDIR)/gramsci_cl
	@echo "Built: $(BINDIR)/gramsci_cl"

# ---- Module dependency rules ----
config_module.o: kdtree2.o node_module.o extension.o
io_module.o: kdtree2.o config_module.o
analytic_randoms_module.o: kdtree2.o config_module.o
graph_utils_module.o: kdtree2.o config_module.o
graph_module.o: kdtree2.o node_module.o sorting_module.o config_module.o graph_utils_module.o
query_2pcf_module.o: kdtree2.o config_module.o analytic_randoms_module.o
query_3pcf_module.o: kdtree2.o config_module.o graph_utils_module.o analytic_randoms_module.o
query_4pcf_module.o: kdtree2.o config_module.o graph_utils_module.o analytic_randoms_module.o

cl_env_module.o: cl_module.o cl_kernels_module.o
csr_cl_module.o: config_module.o
query_3pcf_cl_module.o: config_module.o csr_cl_module.o cl_env_module.o query_3pcf_module.o
query_4pcf_cl_module.o: config_module.o csr_cl_module.o cl_env_module.o query_4pcf_module.o
gramsci_cl.o: $(CPU_OBJS) $(CL_OBJS)

clean:
	/bin/rm -f *.o *.mod gramsci_cl cl_kernels_module.F90
