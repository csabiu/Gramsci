# GPU build of GRAMSCI using nvfortran + OpenACC
#
# Requirements: NVIDIA HPC SDK (nvfortran)
#   https://developer.nvidia.com/hpc-sdk
#
# Usage:
#   make all                    # build with default GPU target (cc86 = Ada Lovelace / RTX 3090 Ti)
#   make GPU=cc70               # Volta
#   make GPU=cc80               # Ampere
#   make GPU=cc90               # Hopper
#   make MANAGED=1              # use CUDA managed memory (simpler, slightly slower)
#   make clean
#
BINDIR = ../bin
SRCDIR = ../src

# Compiler
F90 = nvfortran

# GPU compute capability (override on command line: make GPU=cc70)
# Default is cc86 (RTX 3090 Ti, Ada Lovelace sm_86).
GPU ?= cc86

# CUDA managed memory: set MANAGED=1 for automatic host/device migration.
# With MANAGED=0 (default) explicit DATA regions in the source control
# all transfers, giving better performance.
ifeq ($(MANAGED),1)
  GPU_MEM = ,managed
else
  GPU_MEM =
endif

# Compiler flags
ACC_FLAGS  = -acc -gpu=$(GPU)$(GPU_MEM) -Minfo=accel
BASE_FLAGS = -DDOUBLE -Mfreeform
OPT_FLAGS  = -O3 -fast
FLAGS      = $(BASE_FLAGS) $(OPT_FLAGS) $(ACC_FLAGS) -mp

# ---- Object lists ----

# CPU modules (compiled from ../src/ with nvfortran to generate local .mod files;
# gfortran .mod files in ../src/ are NOT compatible with nvfortran)
CPU_OBJS = extension.o kdtree2.o node_module.o sorting_module.o \
           config_module.o io_module.o analytic_randoms_module.o \
           graph_utils_module.o \
           graph_module.o query_2pcf_module.o \
           query_3pcf_module.o query_4pcf_module.o

# GPU-specific modules (compiled from src_gpu/)
GPU_OBJS = csr_module.o query_3pcf_gpu_module.o query_4pcf_gpu_module.o

ALL_OBJS = $(CPU_OBJS) $(GPU_OBJS) gramsci_gpu.o

# ---- Build rules ----

all: gramsci_gpu

# GPU-specific source files (in src_gpu/)
csr_module.o: csr_module.F90
	$(F90) $(FLAGS) -c $< -o $@

query_3pcf_gpu_module.o: query_3pcf_gpu_module.F90
	$(F90) $(FLAGS) -c $< -o $@

query_4pcf_gpu_module.o: query_4pcf_gpu_module.F90
	$(F90) $(FLAGS) -c $< -o $@

gramsci_gpu.o: gramsci_gpu.F90
	$(F90) $(FLAGS) -c $< -o $@

# CPU source files (compiled from ../src/ to generate nvfortran-compatible .mod files)
extension.o: $(SRCDIR)/extension.F90
	$(F90) $(FLAGS) -c $< -o $@

kdtree2.o: $(SRCDIR)/kdtree2.F90
	$(F90) $(FLAGS) -c $< -o $@

node_module.o: $(SRCDIR)/node_module.F90
	$(F90) $(FLAGS) -c $< -o $@

sorting_module.o: $(SRCDIR)/sorting_module.F90
	$(F90) $(FLAGS) -c $< -o $@

config_module.o: $(SRCDIR)/config_module.F90
	$(F90) $(FLAGS) -c $< -o $@

io_module.o: $(SRCDIR)/io_module.F90
	$(F90) $(FLAGS) -c $< -o $@

analytic_randoms_module.o: $(SRCDIR)/analytic_randoms_module.F90
	$(F90) $(FLAGS) -c $< -o $@

graph_utils_module.o: $(SRCDIR)/graph_utils_module.F90
	$(F90) $(FLAGS) -c $< -o $@

graph_module.o: $(SRCDIR)/graph_module.F90
	$(F90) $(FLAGS) -c $< -o $@

query_2pcf_module.o: $(SRCDIR)/query_2pcf_module.F90
	$(F90) $(FLAGS) -c $< -o $@

query_3pcf_module.o: $(SRCDIR)/query_3pcf_module.F90
	$(F90) $(FLAGS) -c $< -o $@

query_4pcf_module.o: $(SRCDIR)/query_4pcf_module.F90
	$(F90) $(FLAGS) -c $< -o $@

gramsci_gpu: $(ALL_OBJS)
	$(F90) $(FLAGS) -o gramsci_gpu $(ALL_OBJS)
	cp gramsci_gpu $(BINDIR)/gramsci_gpu
	@echo "Built: $(BINDIR)/gramsci_gpu"

# ---- Module dependency rules ----
# (nvfortran needs .mod files before a module's users are compiled)

config_module.o: kdtree2.o node_module.o extension.o
io_module.o: kdtree2.o config_module.o
analytic_randoms_module.o: kdtree2.o config_module.o
graph_utils_module.o: kdtree2.o config_module.o
graph_module.o: kdtree2.o node_module.o sorting_module.o config_module.o graph_utils_module.o
query_2pcf_module.o: kdtree2.o config_module.o analytic_randoms_module.o
query_3pcf_module.o: kdtree2.o config_module.o graph_utils_module.o analytic_randoms_module.o
query_4pcf_module.o: kdtree2.o config_module.o graph_utils_module.o analytic_randoms_module.o

csr_module.o: config_module.o
query_3pcf_gpu_module.o: config_module.o csr_module.o query_3pcf_module.o
query_4pcf_gpu_module.o: config_module.o csr_module.o query_4pcf_module.o
gramsci_gpu.o: $(CPU_OBJS) $(GPU_OBJS)

clean:
	/bin/rm -f *.o *.mod gramsci_gpu
