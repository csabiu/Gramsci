# for gfortran
#
BINDIR="../bin/"

FLAGSALWAYS = -DDOUBLE -ffree-line-length-none 
#-Wall -ffree-form -DDOUBLE
FLAGSOPT=
#-O3 -ffast-math 
FLAGSDEBUG= -g -fbounds-check -O0

F90=/usr/local/bin/mpif90 -DMPI -DGFORTRAN -fopenmp -Ofast -ffast-math -march=native -mtune=intel -faggressive-loop-optimizations -ftree-loop-if-convert 
#F90=gfortran-8 -DGFORTRAN -fopenmp -Ofast -ffast-math -march=native -mtune=intel -faggressive-loop-optimizations -ftree-loop-if-convert 
#F90=ifort -O3 -fopenmp -132 -parallel -xSSE4.2 -m64 -fast
#F90=ifort -O0 -g -traceback -debug all -gen-interfaces -warn interfaces -fopenmp -132 
#-O3 -fopenmp -132 -xSSE4.2 -m64 -fast
#

#F90MPI=mpif90 -O3 -DMPI -fopenmp -132 
#F90=mpif90 -O3 -DMPI -fopenmp -132 


#LIBS=-L/usr/local/openmpi-1.6.3-mx-intel12/lib/
INC=-I/usr/local/include/ 
#-I/applic/compilers/gcc/4.9.3/mpi/mvapich2/2.1/include/

# choose debugging or optimization
FLAGS= ${FLAGSMACHINE} ${FLAGSALWAYS} ${INC} ${LIBS} ${FLAGSOPT}  #  change the last to ${FLAGSOPT} or ${FLAGSDEBUG}
MPIFLAGS= ${FLAGSMACHINE} ${FLAGSALWAYS} ${INC} ${LIBS} ${FLAGSOPT}  #  change the last to ${FLAGSOPT} or ${FLAGSDEBUG}

MY_DIR=`basename ${PWD}`
#------------------------------------

%.o: %.f90
	$(F90) $(FLAGS) -c $< -o $@

%.o: %.F90
	$(F90) $(FLAGS) -c $< -o $@

all:	gramsci gramsci2 domain_decomposition

gramsci: gramsci.o kdtree2.o extension.o
	${F90} ${FLAGS} -o gramsci gramsci.o kdtree2.o extension.o ; \
	mv $@ $(BINDIR)/gramsci

gramsci.o: kdtree2.o extension.o

gramsci2: gramsci2.o kdtree2.o extension.o
	${F90} ${FLAGS} -o gramsci2 gramsci2.o kdtree2.o extension.o ; \
	mv $@ $(BINDIR)/gramsci2

gramsci2.o: kdtree2.o extension.o

domain_decomposition: domain_decomposition.o
	${F90} ${FLAGS} -o domain_decomposition domain_decomposition.o ; \
	mv $@ $(BINDIR)/domain_decomposition

clean:
	/bin/rm -f *.o *.mod