# ───────────────────────────────────────────────────────────────────
# This is a Makefile to compile an ECGPACK code. Here are some examples of compilations
# with this Makefile:
#    make release COMPILER=gfortran MACHINE=linux-generic PREC=8 LINALG=openblas EXEFILE=ecg
#    make release COMPILER=gfortran MACHINE=linux-generic PREC=8 LINALG=netlib OPENMP=yes EXEFILE=ecg
#    make release COMPILER=gfortran MACHINE=ubuntu-generic PREC=8 LINALG=lblas EXEFILE=ecg
#    make debug COMPILER=gfortran MACHINE=irgetas PREC=10 LINALG=netlib EXEFILE=ecg
#    make release COMPILER=ifort MACHINE=puma PREC=8 LINALG=mkl EXEFILE=ecg
#    make debug COMPILER=ifx MACHINE=muon PREC=16 LINALG=netlib EXEFILE=ecg
#    make release COMPILER=nvfortran MACHINE=shabyt PREC=8 LINALG=lblas EXEFILE=ecg
# ───────────────────────────────────────────────────────────────────

# Note that if using a specific BLAS/LAPACK library or compiler (mpif90 wrapper)
# requires loading moduli on your system (this is quite common in HPC environments)
# then you first need to load the BLAS/LAPACK library module and then the compiler/wrapper
# module before running the make command.

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

# Select which BLAS/LAPACK implementation the linear algebra routines are linked
# against via the LINALG argument. Supported values:
#   netlib   - Non-optimized reference BLAS/LAPACK compiled from the bundled modified Netlib
#              source (default). Results in no additions to LINKFLAGS.
#   mkl      - Intel Math Kernel Library.
#   lblas    - Optimized BLAS/LAPACK exposed through the -lblas/-llapack symbolic links.
#   openblas - OpenBLAS library.
#   aocl     - AMD Optimizing CPU Libraries (AOCL-BLAS and AOCL-LAPACK).
# Choosing an optimized library (any value other than netlib) makes effect only for
# PREC=8, because optimized BLAS/LAPACK libraries are not available in extended or
# quadruple precision. For PREC=10 and PREC=16 the value of LINALG must be netlib.
LINALG=netlib

# Enable OpenMP in the bundled qrlinalg BLAS kernels. This is opt-in so the
# existing serial build remains the default. OpenMP and serial objects use
# different directories; changing OPENMP therefore cannot silently reuse
# objects compiled with incompatible flags.
# Supported values: no (serial, default), yes (OpenMP).
OPENMP ?= no

ifneq ($(OPENMP),$(filter $(OPENMP),no yes))
  $(error OPENMP must be no or yes)
endif

OPENMP_SUFFIX := $(if $(filter yes,$(OPENMP)),-omp)
BUILD_DIR := $(CONFIG)$(OPENMP_SUFFIX)

# Set the default name for the generated executable file
EXEFILE=ecg

# Directory where all source files are located
SRC_DIR=./src

# Set compiler flags for each machine and compiler type (for both release and debug builds).

ifneq ($(filter $(MACHINE),linux-generic ubuntu-generic irgetas shabyt muon puma ocelote elgato),)
  # gfortran and Open MPI (default on most Linux systems)
  ifeq ($(COMPILER), gfortran)
    FC=mpif90
    QR_FIXED_FLAGS=-ffixed-line-length-none
    QR_PREPROCESS_FLAGS=-cpp
    export MPICOMPFLAGS=
    export MPILINKFLAGS=
    OPENMP_FLAG=-fopenmp
    export FCFLAGS=-J$(BUILD_DIR)
    export LINKFLAGS=-J$(BUILD_DIR)
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
    QR_FIXED_FLAGS=-extend-source
    QR_PREPROCESS_FLAGS=-fpp
    export MPICOMPFLAGS=
    export MPILINKFLAGS=
    OPENMP_FLAG=-qopenmp
    export FCFLAGS=-module $(BUILD_DIR)
    export LINKFLAGS=-module $(BUILD_DIR)
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
    QR_FIXED_FLAGS=-extend-source
    QR_PREPROCESS_FLAGS=-fpp
    export MPICOMPFLAGS=
    export MPILINKFLAGS=
    OPENMP_FLAG=-qopenmp
    export FCFLAGS=-module $(BUILD_DIR)
    export LINKFLAGS=-module $(BUILD_DIR)
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
    QR_FIXED_FLAGS=-Mextend
    QR_PREPROCESS_FLAGS=-Mpreprocess
    export MPICOMPFLAGS=
    export MPILINKFLAGS=
    OPENMP_FLAG=-mp
    export FCFLAGS=-module $(BUILD_DIR)
    export LINKFLAGS=-module $(BUILD_DIR)
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

  ifeq ($(OPENMP), yes)
    export FCFLAGS+=$(OPENMP_FLAG)
    export LINKFLAGS+=$(OPENMP_FLAG)
  endif

  # Add linker flags for the selected optimized linear algebra library (PREC=8 only).
  # For PREC=10 and PREC=16 only the bundled Netlib source (LINALG=netlib) is used,
  # which requires no additional linker flags.
  ifeq ($(PREC), 8)
    ifeq ($(LINALG), mkl)
      ifeq ($(COMPILER), gfortran)
        export LINKFLAGS+=-L${MKLROOT}/lib -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
      endif
      ifeq ($(COMPILER), nvfortran)
        export LINKFLAGS+=-L${MKLROOT}/lib -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
      endif
      # Intel compilers: let the driver link MKL instead of the manual line below,
      # which avoids the bare "-lm" that fails to resolve on some systems (e.g. puma).
      # ifort uses the older "-mkl" spelling (works across all ifort versions, incl. 2020);
      # ifx (oneAPI only) uses the newer "-qmkl" spelling.
      ifeq ($(COMPILER), ifort)
        # export LINKFLAGS+=-L${MKLROOT}/lib -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
        export LINKFLAGS+=-mkl=sequential
      endif
      ifeq ($(COMPILER), ifx)
        # export LINKFLAGS+=-L${MKLROOT}/lib -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
        export LINKFLAGS+=-qmkl=sequential
      endif
    endif
    ifeq ($(LINALG), lblas)
      export LINKFLAGS+=-llapack -lblas
    endif
    ifeq ($(LINALG), openblas)
      export LINKFLAGS+=-lopenblas
    endif
    ifeq ($(LINALG), aocl)
      export LINKFLAGS+=-L${EBROOTAOCLMINLAPACK}/lib -L${EBROOTAOCLMINBLAS}/lib -lflame -lblis-mt -laoclutils -lstdc++ -lm -lpthread
    endif
  endif
endif

# Creating lists of objects for convenience

ifeq ($(LINALG), netlib)
  export BARE_OBJS_LPKBLS=BLAS.o LAPACK.o
else
  export BARE_OBJS_LPKBLS=
endif
export OBJS_LPKBLS=$(addprefix $(BUILD_DIR)/,$(BARE_OBJS_LPKBLS))

export BARE_OBJS=wp_def_$(PREC).o X1MACH.o dmng.o globvars.o misc.o linalg.o \
              matelem.o matform.o qrupdate_linalg.o qrupdate_error.o \
              qrupdate_real.o qrupdate_complex.o qrupdate.o qrlinalg.o \
              workproc.o main.o spin.o
export OBJS=$(addprefix $(BUILD_DIR)/,$(BARE_OBJS))
export OBJS+=$(OBJS_LPKBLS)

.PHONY: release debug all clean cleaner cleanest cleandebug cleanrelease ptarg

release:
ifndef CONFIG
	-@mkdir -p $(@)$(OPENMP_SUFFIX)
	@$(MAKE) $(MAKEFILE) --no-print-directory $(@) CONFIG=$(@)
else
	@$(MAKE) $(MAKEFILE) --no-print-directory $(BUILD_DIR)/$(EXEFILE)
endif

debug:
ifndef CONFIG
	-@mkdir -p $(@)$(OPENMP_SUFFIX)
	@$(MAKE) $(MAKEFILE) --no-print-directory $(@) CONFIG=$(@)
else
	@$(MAKE) $(MAKEFILE) --no-print-directory $(BUILD_DIR)/$(EXEFILE)
endif

all: release debug

cleanrelease:
	-rm -rf release/*.o release/*.mod release-omp/*.o release-omp/*.mod
cleandebug:
	-rm -rf debug/*.o debug/*.mod debug-omp/*.o debug-omp/*.mod
clean: cleanrelease cleandebug
cleaner:
	-rm -rf debug/*.o debug/*.mod debug/$(EXEFILE)
	-rm -rf debug-omp/*.o debug-omp/*.mod debug-omp/$(EXEFILE)
	-rm -rf release/*.o release/*.mod release/$(EXEFILE)
	-rm -rf release-omp/*.o release-omp/*.mod release-omp/$(EXEFILE)
cleanest:
	-rm -rf debug debug-omp release release-omp

# This phony target is just used to output some info and error check

ptarg:
	@echo CONFIG=$(CONFIG)
ifeq ($(COMPILER_SUPPORTED), yes)
	@echo Using COMPILER=$(COMPILER) MACHINE=$(MACHINE) PREC=$(PREC) LINALG=$(LINALG) OPENMP=$(OPENMP) EXEFILE=$(BUILD_DIR)/$(EXEFILE)
else
	@echo Error - there is no support for COMPILER=$(COMPILER) MACHINE=$(MACHINE) PREC=$(PREC) LINALG=$(LINALG)
	@false
endif

# The real targets

$(BUILD_DIR)/$(EXEFILE): ptarg $(OBJS)
	$(FC) $(OBJS) -o $(BUILD_DIR)/$(EXEFILE) $(LINKFLAGS) $(MPILINKFLAGS)

$(BUILD_DIR)/wp_def_$(PREC).o: $(SRC_DIR)/wp_def_$(PREC).f90
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/wp_def_$(PREC).f90 -o $(BUILD_DIR)/wp_def_$(PREC).o

$(BUILD_DIR)/BLAS.o: $(SRC_DIR)/BLAS.f $(BUILD_DIR)/wp_def_$(PREC).o
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) $(QR_FIXED_FLAGS) $(QR_PREPROCESS_FLAGS) -DQRLINALG_WP=$(PREC) -c $(SRC_DIR)/BLAS.f -o $(BUILD_DIR)/BLAS.o

$(BUILD_DIR)/LAPACK.o: $(SRC_DIR)/LAPACK.f $(BUILD_DIR)/BLAS.o
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) $(QR_FIXED_FLAGS) -c $(SRC_DIR)/LAPACK.f -o $(BUILD_DIR)/LAPACK.o

$(BUILD_DIR)/X1MACH.o: $(SRC_DIR)/X1MACH.f90 $(BUILD_DIR)/wp_def_$(PREC).o
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/X1MACH.f90 -o $(BUILD_DIR)/X1MACH.o

$(BUILD_DIR)/dmng.o: $(SRC_DIR)/dmng.f $(BUILD_DIR)/X1MACH.o
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/dmng.f -o $(BUILD_DIR)/dmng.o

$(BUILD_DIR)/globvars.o: $(SRC_DIR)/globvars.f90 $(BUILD_DIR)/wp_def_$(PREC).o
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/globvars.f90 -o $(BUILD_DIR)/globvars.o

$(BUILD_DIR)/linalg.o: $(SRC_DIR)/linalg.f90 $(BUILD_DIR)/globvars.o $(OBJ_LPKBLS)
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/linalg.f90 -o $(BUILD_DIR)/linalg.o

$(BUILD_DIR)/spin.o: $(SRC_DIR)/spin.f90 $(BUILD_DIR)/wp_def_$(PREC).o
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/spin.f90 -o $(BUILD_DIR)/spin.o

$(BUILD_DIR)/misc.o: $(SRC_DIR)/misc.f90 $(BUILD_DIR)/globvars.o
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/misc.f90 -o $(BUILD_DIR)/misc.o

$(BUILD_DIR)/matelem.o: $(SRC_DIR)/matelem.f90 $(BUILD_DIR)/globvars.o $(BUILD_DIR)/spin.o
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/matelem.f90 -o $(BUILD_DIR)/matelem.o

$(BUILD_DIR)/matform.o: $(SRC_DIR)/matform.f90 $(BUILD_DIR)/matelem.o
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/matform.f90 -o $(BUILD_DIR)/matform.o

$(BUILD_DIR)/qrupdate_linalg.o: $(SRC_DIR)/qrupdate/qrupdate_linalg.f90 $(BUILD_DIR)/wp_def_$(PREC).o
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/qrupdate/qrupdate_linalg.f90 -o $(BUILD_DIR)/qrupdate_linalg.o

$(BUILD_DIR)/qrupdate_error.o: $(SRC_DIR)/qrupdate/qrupdate_error.f90 $(BUILD_DIR)/qrupdate_linalg.o
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/qrupdate/qrupdate_error.f90 -o $(BUILD_DIR)/qrupdate_error.o

$(BUILD_DIR)/qrupdate_real.o: $(SRC_DIR)/qrupdate/qrupdate_real.f90 $(BUILD_DIR)/qrupdate_error.o
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/qrupdate/qrupdate_real.f90 -o $(BUILD_DIR)/qrupdate_real.o

$(BUILD_DIR)/qrupdate_complex.o: $(SRC_DIR)/qrupdate/qrupdate_complex.f90 $(BUILD_DIR)/qrupdate_error.o
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/qrupdate/qrupdate_complex.f90 -o $(BUILD_DIR)/qrupdate_complex.o

$(BUILD_DIR)/qrupdate.o: $(SRC_DIR)/qrupdate/qrupdate.f90 $(BUILD_DIR)/qrupdate_real.o $(BUILD_DIR)/qrupdate_complex.o
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/qrupdate/qrupdate.f90 -o $(BUILD_DIR)/qrupdate.o

$(BUILD_DIR)/qrlinalg.o: $(SRC_DIR)/qrlinalg.f90 $(BUILD_DIR)/qrupdate.o
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) $(QR_PREPROCESS_FLAGS) -c $(SRC_DIR)/qrlinalg.f90 -o $(BUILD_DIR)/qrlinalg.o

$(BUILD_DIR)/workproc.o: $(SRC_DIR)/workproc.f90 $(BUILD_DIR)/misc.o $(BUILD_DIR)/linalg.o $(BUILD_DIR)/matform.o $(BUILD_DIR)/matelem.o $(BUILD_DIR)/dmng.o $(BUILD_DIR)/qrlinalg.o $(OBJS_LPKBLS)
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/workproc.f90 -o $(BUILD_DIR)/workproc.o

$(BUILD_DIR)/main.o: $(SRC_DIR)/main.f90 $(BUILD_DIR)/workproc.o
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/main.f90 -o $(BUILD_DIR)/main.o

