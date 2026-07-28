# ───────────────────────────────────────────────────────────────────
# This is a Makefile to compile an ECGPACK code. Here are some examples of compilations
# with this Makefile:
#    make release COMPILER=gfortran MACHINE=linux-generic PREC=8 LINALG=netlib EXEFILE=ecg
#    make release COMPILER=gfortran MACHINE=ubuntu-generic PREC=8 EXEFILE=ecg
#    make debug COMPILER=gfortran MACHINE=irgetas PREC=10 LINALG=netlib EXEFILE=ecg
#    make release COMPILER=ifort MACHINE=puma PREC=8 EXEFILE=ecg
#    make debug COMPILER=ifx MACHINE=muon PREC=16 LINALG=netlib EXEFILE=ecg
#    make release COMPILER=nvfortran MACHINE=shabyt PREC=8 EXEFILE=ecg
# ───────────────────────────────────────────────────────────────────

# Note that if using a specific compiler (mpif90 wrapper) requires loading a module on
# your system (this is quite common in HPC environments) then you first need to load
# that module before running the make command.

# Set the default compiler
# Supported options: gfortran, ifort, ifx, nvfortran
# Note that instead of thees compilers their corresponding MPI wrappers are actually 
# invoked (e.g. mpif90 instead of gfortran; mpiifort instead of ifort, etc.)
COMPILER=gfortran

# Set the default machine where the compilation takes place 
# Supported options: linux-generic, ubuntu-generic, irgetas, shabyt, muon, puma, ocelote, elgato
MACHINE=linux-generic

# Set the default precision (real kind: 8, 10, or 16)
# Value PREC=8 corresponds to standard double precision
# Value PREC=10 corresponds to extended precision
# Value PREC=16 corresponds to quadruple precision
# Note that PREC=10 or PREC=16 support is compiler and platform dependent and is not universal. 
PREC=8

# LINALG selects the BLAS/LAPACK implementation used by the energy codes. The
# off-diagonal matrix-element codes do not link an external BLAS/LAPACK library,
# so here LINALG is only a dummy argument; it must be left at its default value
# netlib (which is also the only value allowed for PREC=10 and PREC=16).
LINALG=netlib

# Set the default name for the generated executable file
EXEFILE=ecg

# Directory where all source files are located
SRC_DIR=./src

# Set compiler flags for each machine and compiler type (for both release and debug builds).

ifneq ($(filter $(MACHINE),linux-generic ubuntu-generic irgetas shabyt muon puma ocelote elgato),)
  # gfortran and Open MPI (default on most Linux systems)
  ifeq ($(COMPILER), gfortran)
    FC=mpif90
    export MPICOMPFLAGS=
    export MPILINKFLAGS=
    export FCFLAGS=-J$(CONFIG)
    export LINKFLAGS=-J$(CONFIG)
    ifeq ($(CONFIG), debug)
      export FCFLAGS+=-C -O0 -g -fbounds-check -fbacktrace -finit-real=nan -fcheck=all -fimplicit-none -ffpe-trap=invalid,zero,overflow
      export LINKFLAGS+=-C -O0 -g -fbounds-check -fbacktrace -finit-real=nan -fcheck=all -fimplicit-none -ffpe-trap=invalid,zero,overflow
    endif
    ifeq ($(CONFIG), release)
      export FCFLAGS+=-O3 -march=native
      export LINKFLAGS+=-O3 -march=native
    endif
    ifeq ($(PREC), 8)
      COMPILER_SUPPORTED=yes
    endif
    ifeq ($(PREC), 10)
      ifeq ($(LINALG),netlib)
        COMPILER_SUPPORTED=yes
      endif
    endif
    ifeq ($(PREC), 16)
      ifeq ($(LINALG),netlib)
        COMPILER_SUPPORTED=yes
      endif
    endif
  endif
  # Intel Fortran (ifort) and Intel MPI
  ifeq ($(COMPILER), ifort)
    FC=mpiifort
    export MPICOMPFLAGS=
    export MPILINKFLAGS=
    export FCFLAGS=-module $(CONFIG)
    export LINKFLAGS=-module $(CONFIG)
    ifeq ($(CONFIG), debug)
      export FCFLAGS+=-C -O0 -g -check bounds -check all -check uninit -traceback -ftrapuv -debug full -fpe0
      export LINKFLAGS+=-C -O0 -g -check bounds -check all -check uninit -traceback -ftrapuv -debug full -fpe0
    endif
    ifeq ($(CONFIG), release)
      export FCFLAGS+=-O3 -ip -xhost -fp-model precise
      export LINKFLAGS+=-O3 -ip -xhost -fp-model precise
    endif
    ifeq ($(PREC), 8)
      COMPILER_SUPPORTED=yes
    endif
    ifeq ($(PREC), 16)
      ifeq ($(LINALG),netlib)
        COMPILER_SUPPORTED=yes
      endif
    endif
  endif
  # Intel Fortran (LLVM-based ifx) and Intel MPI
  ifeq ($(COMPILER), ifx)
    FC=mpiifx
    export MPICOMPFLAGS=
    export MPILINKFLAGS=
    export FCFLAGS=-module $(CONFIG)
    export LINKFLAGS=-module $(CONFIG)
    ifeq ($(CONFIG), debug)
      export FCFLAGS+=-O0 -g -check all -check uninit -traceback -ftrapuv -fpe0
      export LINKFLAGS+=-O0 -g -check all -check uninit -traceback -ftrapuv -fpe0
    endif
    ifeq ($(CONFIG), release)
      export FCFLAGS+=-O3 -xhost -fp-model precise
      export LINKFLAGS+=-O3 -xhost -fp-model precise
    endif
    ifeq ($(PREC), 8)
      COMPILER_SUPPORTED=yes
    endif
    ifeq ($(PREC), 16)
      ifeq ($(LINALG),netlib)
        COMPILER_SUPPORTED=yes
      endif
    endif
  endif
  # NVIDIA HPC SDK (nvfortran) and Open MPI
  ifeq ($(COMPILER), nvfortran)
    FC=mpif90
    export MPICOMPFLAGS=
    export MPILINKFLAGS=
    export FCFLAGS=-module $(CONFIG)
    export LINKFLAGS=-module $(CONFIG)
    ifeq ($(CONFIG), debug)
      export FCFLAGS+=-O0 -g -Mbounds -Ktrap=fp -traceback
      export LINKFLAGS+=-O0 -g -Mbounds -Ktrap=fp -traceback
    endif
    ifeq ($(CONFIG), release)
      export FCFLAGS+=-O3 -tp=native
      export LINKFLAGS+=-O3 -tp=native
    endif
    ifeq ($(PREC), 8)
      COMPILER_SUPPORTED=yes
    endif
    ifeq ($(PREC), 10)
      COMPILER_SUPPORTED=no
    endif
    ifeq ($(PREC), 16)
      COMPILER_SUPPORTED=no
    endif
  endif
endif

# Creating lists of objects for convenience

export BARE_OBJS=wp_def_$(PREC).o globvars.o matelem.o workproc.o main.o
export OBJS=$(addprefix $(CONFIG)/,$(BARE_OBJS))

.PHONY: release debug all clean cleaner cleanest cleandebug cleanrelease ptarg

release:
ifndef CONFIG
	-@mkdir -p $(@)
	@$(MAKE) $(MAKEFILE) --no-print-directory $(@) CONFIG=$(@)
else
	@$(MAKE) $(MAKEFILE) --no-print-directory $(CONFIG)/$(EXEFILE)
endif

debug:
ifndef CONFIG
	-@mkdir -p $(@)
	@$(MAKE) $(MAKEFILE) --no-print-directory $(@) CONFIG=$(@)
else
	@$(MAKE) $(MAKEFILE) --no-print-directory $(CONFIG)/$(EXEFILE)
endif

all: release debug

cleanrelease:
	-rm -rf release/*.o release/*.mod	
cleandebug:
	-rm -rf debug/*.o debug/*.mod
clean: cleanrelease cleandebug
cleaner:
	-rm -rf debug/*.o debug/*.mod debug/$(EXEFILE)
	-rm -rf release/*.o release/*.mod release/$(EXEFILE)
cleanest:
	-rm -rf debug release

# This phony target is just used to output some info and error check 

ptarg:
	@echo CONFIG=$(CONFIG)
ifeq ($(COMPILER_SUPPORTED), yes)
	@echo Using COMPILER=$(COMPILER) MACHINE=$(MACHINE) PREC=$(PREC) EXEFILE=$(CONFIG)/$(EXEFILE)
else
	@echo Error - there is no support for COMPILER=$(COMPILER) MACHINE=$(MACHINE) PREC=$(PREC)
	@false
endif

# The real targets

$(CONFIG)/$(EXEFILE): ptarg $(OBJS)
	$(FC) $(OBJS) -o $(CONFIG)/$(EXEFILE) $(LINKFLAGS) $(MPILINKFLAGS)

$(CONFIG)/wp_def_$(PREC).o: $(SRC_DIR)/wp_def_$(PREC).f90
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/wp_def_$(PREC).f90 -o $(CONFIG)/wp_def_$(PREC).o

$(CONFIG)/globvars.o: $(SRC_DIR)/globvars.f90 $(CONFIG)/wp_def_$(PREC).o
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/globvars.f90 -o $(CONFIG)/globvars.o

$(CONFIG)/matelem.o: $(SRC_DIR)/matelem.f90 $(CONFIG)/globvars.o
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/matelem.f90 -o $(CONFIG)/matelem.o

$(CONFIG)/workproc.o: $(SRC_DIR)/workproc.f90 $(CONFIG)/matelem.o
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/workproc.f90 -o $(CONFIG)/workproc.o

$(CONFIG)/main.o: $(SRC_DIR)/main.f90 $(CONFIG)/workproc.o
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/main.f90 -o $(CONFIG)/main.o

