PROGRAM = fracsca
EXEFNAME = $(PROGRAM).x

# 1. siingle wavelength
SRCS = opacfractal.f90 geofractal.f90 call.f90

# 2. multiple wavelengths
#SRCS = opacfractal.f90 geofractal.f90 call2.f90

# 3. multiple wavelengths + grain size distribution
#SRCS = opacfractal.f90 geofractal.f90 call3.f90

OBJS = $(SRCS:.f90=.o) 

#debug  = true
#openmp = true

FC = gfortran
#FC = ifort

ifeq ($(FC),gfortran)
   ifeq ($(debug),true)
      DEBUGGING  = -fbounds-check -fbacktrace -O0 -ffpe-trap=invalid,zero,overflow -fcheck=all
   endif
   ifeq ($(openmp),true)
      OMP= -fopenmp
   endif
else 
   ifeq ($(debug),true)
      DEBUGGING = -check all -traceback -check bounds -O0 -g -check -fpe1
   endif
   ifeq ($(openmp),true)
      OMP= -openmp -fp-model strict
   endif
endif
FCFLAGS = $(DEBUGGING) $(OMP)

LD = $(FC)
LDFLAGS = $(FCFLAGS)
INCLUDE = 
LIBS = 
MAN = 
SHELL = /bin/sh
ALLFILES = $(SRCS) $(DATA) Makefile

.SUFFIXES: .f90

#
# Target
#
all: $(PROGRAM)

.f90.o:
	$(FC) $(FCFLAGS) -c $< $(FCFLAGS) -o $@

tarball:
	@tar zcvf $(PROGRAM).tgz $(ALLFILES)

clean:
	@rm -f $(OBJS) $(OBJS:.o=.mod) $(EXEFNAME)

cleanall:
	@rm -f *.o *.out *.sca

#
# Target Command
#
$(PROGRAM): $(OBJS) $(INCLUDE)
	@echo "Linking $(PROGRAM)..."
	@$(LD) $(LDFLAGS) -o $(EXEFNAME) $(OBJS) $(LIBS)
	@echo "Done."

$(PROGRAM).f: $(SRCS) $(INCLUDE)
	@cat $(SRCS) $(INCLUDE) > $(PROGRAM).f
