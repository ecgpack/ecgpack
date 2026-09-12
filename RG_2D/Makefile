# ───────────────────────────────────────────────────────────────────
# This is a Makefile to compile an ECGPACK code. Here are some examples of compilations
# with this Makefile:
#    make release COMPILER=gfortran MACHINE=linux-generic PREC=8 LINALG=openblas EXEFILE=ecg
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

# GPU backend. Requires COMPILER=nvfortran (the NVIDIA HPC SDK): the host Fortran
# and the device kernels are both compiled by nvfortran -cuda, so there is a single
# toolchain (nvfortran + HPC-X MPI). Load the nvhpc module first. GPU builds are
# PREC=8 only and use LINALG=netlib so the bundled reference BLAS/LAPACK are
# compiled by nvfortran (avoids linking a system optimized library against nvfortran).
#
#   make release COMPILER=nvfortran MACHINE=shabyt USE_CUDA=yes LINALG=netlib   # arch from MACHINE
#
# USE_CUDA: build the native CUDA Fortran GPU backend (src/gpu_backend.f90,
# nvfortran -cuda) -- orchestration, matrix-element kernels, device lifecycle, and
# the cuSOLVER eigensolver, all Fortran in one module. Ported from RG_0S.
USE_CUDA=no
# CUDA_ARCH: the GPU compute capability (nvfortran -gpu=ccXX). It defaults
# from MACHINE so you don't have to remember it per host; override by passing
# CUDA_ARCH=... explicitly for any other GPU.
#   sm_70 = V100 (shabyt)   sm_80 = A100   sm_86 = RTX 3060 Ti (aurora)   sm_90 = H100 (irgetas)
CUDA_ARCH_shabyt = sm_70
CUDA_ARCH_aurora = sm_86
CUDA_ARCH_irgetas = sm_90
CUDA_ARCH = $(or $(CUDA_ARCH_$(MACHINE)),sm_86)

# CUF_GPU is the nvfortran -gpu=ccXX arch, derived from CUDA_ARCH (sm_70 -> cc70).
CUF_GPU=$(subst sm_,cc,$(CUDA_ARCH))
# CUF_MAXREG: cap registers/thread for the GPU device code. This is a global
# -gpu cap (not per-kernel launch_bounds) because MatrixElementsHS_RG_2D would
# otherwise be a non-inlined device callee that launch_bounds cannot reach.
# 255 (no effective cap) + device LTO is inherited from the RG_0S backend, where
# lower caps were measured slower AND (for the optimized routine) shown to
# MISCOMPILE on NVHPC 25.9 -- see RG_0S/Makefile for the full history. Do not
# lower this for RG_2D without an RG_2D-specific correctness gate.
CUF_MAXREG=255

# Set the default name for the generated executable file
EXEFILE=ecg

# Directory where all source files are located
SRC_DIR=./src

# Set compiler flags for each machine and compiler type (for both release and debug builds).

ifneq ($(filter $(MACHINE),linux-generic ubuntu-generic irgetas shabyt aurora muon puma ocelote elgato),)
  # gfortran and Open MPI (default on most Linux systems)
  ifeq ($(COMPILER), gfortran)
    FC=mpif90
    export MPICOMPFLAGS=
    export MPILINKFLAGS=
    export FCFLAGS=-J$(CONFIG) -cpp
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
    export FCFLAGS=-module $(CONFIG) -fpp
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
    export FCFLAGS=-module $(CONFIG) -fpp
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
  # NVIDIA HPC SDK (nvfortran) and Open MPI. Also the host toolchain for the CUDA
  # backend (USE_CUDA); -Mpreprocess handles the #ifdef directives in matform.f90
  # and the -DUSE_CUDA flag the backend adds. (nvfortran has no -fpp switch.)
  ifeq ($(COMPILER), nvfortran)
    FC=mpif90
    export MPICOMPFLAGS=
    export MPILINKFLAGS=
    export FCFLAGS=-module $(CONFIG) -Mpreprocess -Kieee
    export LINKFLAGS=-module $(CONFIG) -Kieee
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
export OBJS_LPKBLS=$(addprefix $(CONFIG)/,$(BARE_OBJS_LPKBLS))

export BARE_OBJS=wp_def_$(PREC).o X1MACH.o dmng.o globvars.o misc.o linalg.o \
              matelem.o matform.o workproc.o main.o spin.o
export OBJS=$(addprefix $(CONFIG)/,$(BARE_OBJS))
export OBJS+=$(OBJS_LPKBLS)

# CUDA backend (USE_CUDA=yes). nvfortran -cuda compiles the device kernels and pulls
# in the CUDA runtime; -cudalib links cuSOLVER/cuBLAS. Requires COMPILER=nvfortran and PREC=8.
export CUDA_OBJS=
GPU_DEP=
ifeq ($(USE_CUDA), yes)
  ifneq ($(COMPILER), nvfortran)
    $(error USE_CUDA=yes requires COMPILER=nvfortran)
  endif
  ifneq ($(PREC), 8)
    $(error USE_CUDA=yes is only supported with PREC=8)
  endif
  # Native CUDA Fortran GPU build: the whole backend -- run-time orchestration,
  # the matrix-element kernels + host launchers, the device lifecycle, and the
  # cuSOLVER eigensolver -- lives in one module, gpu_backend.f90, compiled by
  # nvfortran -cuda. cuSOLVER/cuBLAS are linked for the eigensolver. -DUSE_CUDA guards BOTH
  # the backend (gpu_backend.f90 + the dispatch touchpoints) and the
  # attributes(host,device) on MatrixElementsHS_RG_2D so the shared physics is device-capable.
  export FCFLAGS+=-DUSE_CUDA -cuda -gpu=$(CUF_GPU),maxregcount:$(CUF_MAXREG),lto
  export LINKFLAGS+=-cuda -gpu=$(CUF_GPU),maxregcount:$(CUF_MAXREG),lto -cudalib=cusolver,cublas
  export BARE_OBJS+=gpu_backend.o
  export OBJS=$(addprefix $(CONFIG)/,$(BARE_OBJS)) $(OBJS_LPKBLS)
  GPU_DEP=$(CONFIG)/gpu_backend.o
  GPUBACK_MODDEP=$(CONFIG)/matelem.o
endif

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
	@echo Using COMPILER=$(COMPILER) MACHINE=$(MACHINE) PREC=$(PREC) LINALG=$(LINALG) EXEFILE=$(CONFIG)/$(EXEFILE)$(if $(filter yes,$(USE_CUDA)), USE_CUDA=yes CUDA_ARCH=$(CUDA_ARCH),)
else
	@echo Error - there is no support for COMPILER=$(COMPILER) MACHINE=$(MACHINE) PREC=$(PREC) LINALG=$(LINALG)
	@false
endif

# The real targets

$(CONFIG)/$(EXEFILE): ptarg $(OBJS) $(CUDA_OBJS)
	$(FC) $(OBJS) $(CUDA_OBJS) -o $(CONFIG)/$(EXEFILE) $(LINKFLAGS) $(MPILINKFLAGS)

$(CONFIG)/wp_def_$(PREC).o: $(SRC_DIR)/wp_def_$(PREC).f90
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/wp_def_$(PREC).f90 -o $(CONFIG)/wp_def_$(PREC).o

$(CONFIG)/BLAS.o: $(SRC_DIR)/BLAS.f $(CONFIG)/wp_def_$(PREC).o
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/BLAS.f -o $(CONFIG)/BLAS.o

$(CONFIG)/LAPACK.o: $(SRC_DIR)/LAPACK.f $(CONFIG)/BLAS.o
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/LAPACK.f -o $(CONFIG)/LAPACK.o

$(CONFIG)/X1MACH.o: $(SRC_DIR)/X1MACH.f90 $(CONFIG)/wp_def_$(PREC).o
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/X1MACH.f90 -o $(CONFIG)/X1MACH.o

$(CONFIG)/dmng.o: $(SRC_DIR)/dmng.f $(CONFIG)/X1MACH.o
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/dmng.f -o $(CONFIG)/dmng.o

$(CONFIG)/globvars.o: $(SRC_DIR)/globvars.f90 $(CONFIG)/wp_def_$(PREC).o
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/globvars.f90 -o $(CONFIG)/globvars.o

$(CONFIG)/linalg.o: $(SRC_DIR)/linalg.f90 $(CONFIG)/globvars.o $(OBJ_LPKBLS)
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/linalg.f90 -o $(CONFIG)/linalg.o

$(CONFIG)/spin.o: $(SRC_DIR)/spin.f90 $(CONFIG)/wp_def_$(PREC).o
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/spin.f90 -o $(CONFIG)/spin.o

$(CONFIG)/misc.o: $(SRC_DIR)/misc.f90 $(CONFIG)/globvars.o
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/misc.f90 -o $(CONFIG)/misc.o

$(CONFIG)/matelem.o: $(SRC_DIR)/matelem.f90 $(CONFIG)/globvars.o $(CONFIG)/spin.o
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/matelem.f90 -o $(CONFIG)/matelem.o

$(CONFIG)/gpu_backend.o: $(SRC_DIR)/gpu_backend.f90 $(CONFIG)/globvars.o $(GPUBACK_MODDEP)
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/gpu_backend.f90 -o $(CONFIG)/gpu_backend.o

$(CONFIG)/matform.o: $(SRC_DIR)/matform.f90 $(CONFIG)/matelem.o $(GPU_DEP)
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/matform.f90 -o $(CONFIG)/matform.o

$(CONFIG)/workproc.o: $(SRC_DIR)/workproc.f90 $(CONFIG)/misc.o $(CONFIG)/linalg.o $(CONFIG)/matform.o $(CONFIG)/matelem.o $(CONFIG)/dmng.o $(GPU_DEP) $(OBJS_LPKBLS)
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/workproc.f90 -o $(CONFIG)/workproc.o

$(CONFIG)/main.o: $(SRC_DIR)/main.f90 $(CONFIG)/workproc.o $(GPU_DEP)
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/main.f90 -o $(CONFIG)/main.o
