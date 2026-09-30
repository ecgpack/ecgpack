MODULE globvars

  ! This module contains declarations of global variables and constants.

  USE wp_def
  ! include 'mpif.h'
  USE mpi

  IMPLICIT NONE


  !=============================================================
  ! Global parameters
  !=============================================================
  INTEGER, PARAMETER :: Verbose = 3
  ! Screen output is gated as (Verbose >= n). Verbose is a PARAMETER, so
  ! the compiler removes the suppressed WRITE statements entirely. Errors
  ! and fatal diagnostics are never gated. Levels (cumulative):
  !   0 : errors, results (energies, expectation values, saved files)
  !   1 : + warnings, routine start/finish lines, summaries and statistics
  !   2 : + progress messages (matrix elements, eigenproblem, swap file,
  !         reallocations), per-candidate lines of the stochastic search,
  !         the echo of the data file
  !   3 : + tracing of the inverse-iteration solver (GSEPIIS) and of the
  !         matrix-element assembly (matform): one block per solve/call
  !   4 : + one line per basis-function pair inside matelem (floods)


  !===============================================================
  !                   Numerical constants
  !===============================================================

  REAL(wp), PARAMETER :: &
            ZERO = 0.E0_wp, &
            ONE = 1.E0_wp, &
            TWO = 2.E0_wp, &
            THREE = 3.E0_wp, &
            FOUR = 4.E0_wp, &
            FIVE = 5.E0_wp, &
            SIX = 6.E0_wp, &
            SEVEN = 7.E0_wp, &
            EIGHT = 8.E0_wp, &
            NINE = 9.E0_wp, &
            TEN = 10.0_wp, &
            ONEHALF = ONE/TWO, &
            ONETHIRD = ONE/THREE, &
            ONEFOURTH = ONE/FOUR, &
            THREEHALF = THREE/TWO, &
            PIm12 = 0.564189583547756286948079451560773E0_wp, &
            PI12 = 1.77245385090551602729816748334115E0_wp, &
            Glob_Pi = 3.1415926535897932384626433832795029E0_wp, &
            PI2 = 9.86960440108935861883449099987615E0_wp, &
            PI52 = 17.4934183276248628462628216798716E0_wp, &
            PI3 = 31.0062766802998201754763150671014E0_wp, &
            PI72 = 54.9571945042393159051835760441529E0_wp, &
            PI4 = 97.4090910340024372364403326887051E0_wp, &
            PI92 = 172.653118516423593924947460372192E0_wp, &
            PI5 = 306.019684785281453262741310043436E0_wp, &
            PI112 = 542.405768750564264410946425441073E0_wp, &
            PI6 = 961.389193575304437030219443652420E0_wp, &
            PI132 = 1704.01797837149693758999361474139E0_wp, &
            PI7 = 3020.29322777679206751420649307204E0_wp, &
            PI152 = 5353.33036243682596560701525479963E0_wp, &
            SQRTPI = 1.7724538509055160272981674833411452E0_wp, &
            ! physical constants
            Glob_EulerConst = 0.57721566490153286060651209008240E0_wp, &
            Glob_FineStructConst = 7.2973525693E-03_wp  ! CODATA 2018


  !================================================================
  ! Messages
  !================================================================

  INTEGER :: Glob_Warning = 0


  !================================================================
  ! Physical system: size and particle count
  !================================================================
  ! Number of particles this build is compiled for: Glob_AllowedNumOfParticles is
  ! a PARAMETER of module wp_def (wp_def_$(PREC).f90). ReadIOFile requires the
  ! data file to specify EXACTLY this many particles, so it is an exact count
  ! rather than an upper bound.

  ! Number of pseudoparticles this build is compiled for (N - 1)
  INTEGER, PARAMETER :: Glob_AllowedNumOfPseudoParticles = &
                            Glob_AllowedNumOfParticles - 1

  ! Number of pseudoparticles (N - 1), fixed at compile time as in the
  ! reference, so that arrays and loops sized by it have constant bounds
  INTEGER, PARAMETER :: Glob_n = Glob_AllowedNumOfPseudoParticles


  !================================================================
  ! Basis type
  !================================================================

  ! Short name for the type of basis this build implements. The data file may
  ! carry an optional BASIS_TYPE line; when it does, ReadIOFile checks it
  ! against this constant and refuses a file written for a different basis.
  CHARACTER(5), PARAMETER :: Glob_BasisType = 'PG_0S'

  ! .true. when the data file that was read actually contained a BASIS_TYPE
  ! line, so that SaveResults writes the line back out only if it was there.
  LOGICAL :: Glob_BasisTypeSupplied = .FALSE.


  !================================================================
  ! ECG basis: parameter counts and normalization prefactors
  !================================================================

  ! np = n(n+1)/2
  !   Number of independent elements in a symmetric (n x n) matrix.
  INTEGER, PARAMETER :: Glob_np = Glob_n*(Glob_n+1)/2

  ! Total number of nonlinear parameters per basis function.
  ! For real P-Gaussians, Glob_npt = Glob_np.
  INTEGER, PARAMETER :: Glob_npt = Glob_np

  ! Normalization prefactors
  REAL(wp), PARAMETER :: Glob_2Raised3n2 = TWO**((3*Glob_n)/TWO)       ! 2^(3n/2)
  REAL(wp), PARAMETER :: Glob_PiRaised3n2 = Glob_Pi**((3*Glob_n)/TWO)  ! pi^(3n/2)

  ! Nonlinear parameters (Cholesky factor L_k) of every basis function.
  ! Shape: (Glob_npt, Glob_CurrBasisSize)
  REAL(wp), ALLOCATABLE, DIMENSION(:, :), SAVE :: Glob_NonlinParam

  ! Basis function numbering (Glob_CurrBasisSize entries)
  INTEGER, ALLOCATABLE, DIMENSION(:), SAVE :: Glob_FuncNum

  ! Powers of the r-premultiplier for each basis function
  ! (0 .. Glob_MaxPowerAllowed; 0 means no premultiplier).
  INTEGER, ALLOCATABLE, DIMENSION(:), SAVE :: Glob_PWR

  ! Glob_IsIndexFixed specifies whether the z-indices (r-premultiplier powers)
  ! of all basis functions are held fixed at one common value instead of being
  ! generated and optimized per function. Set from the optional FIXED_INDEX
  ! line of the data file.
  LOGICAL :: Glob_IsIndexFixed = .FALSE.

  ! When Glob_IsIndexFixed is .true. this holds the common z-index value.
  ! Otherwise it is not used.
  INTEGER :: Glob_IndexFixedValue

  INTEGER :: Glob_VariedParam


  !================================================================
  ! Power limits for the r-premultiplier
  !================================================================
  ! Bounded by the tables in data_gamma.f90, indexed with the power/2 for
  ! gm1 and with (mk+ml)/2 for gm2/gm3/lngamma:
  !     Glob_gm1(0:P/2,0:P/2)   Glob_gm2(0:P)   Glob_gm3(0:P)   Glob_lngamma(0:P)
  ! Glob_MaxPowerPossible must not exceed the P the tables were generated
  ! for (250, gen_data_gamma.py); ReadIOFile rejects larger powers.

  ! Largest r-premultiplier power the code can in principle handle
  INTEGER, PARAMETER :: Glob_MaxPowerPossible = 250

  ! Largest r-premultiplier power actually permitted at run time
  INTEGER, PARAMETER :: Glob_MaxPowerAllowed = 250


  !================================================================
  ! Current basis size and current energy
  !================================================================

  ! Current size of the basis
  INTEGER :: Glob_CurrBasisSize

  ! Current energy value
  REAL(wp) :: Glob_CurrEnergy


  !================================================================
  ! Eigenvalue solver parameters
  !================================================================

  ! Which eigenvalue to target (used when GSEPSolutionMethod = 'G')
  INTEGER :: Glob_WhichEigenvalue

  ! Approximate eigenvalue (used when GSEPSolutionMethod = 'I').
  ! Typically Glob_ApproxEnergy = Glob_CurrEnergy * Glob_InvItParameter.
  REAL(wp) :: Glob_ApproxEnergy

  ! Requested accuracy for the eigenvalue problem (method 'I')
  REAL(wp) :: Glob_EigvalTol

  ! Factor by which the approximate eigenvalue is multiplied
  ! when method 'I' is used
  REAL(wp) :: Glob_InvItParameter

  ! Accuracy trackers (method 'I' only)
  REAL(wp) :: Glob_LastEigvalTol   ! last solver call
  REAL(wp) :: Glob_BestEigvalTol   ! best accuracy seen
  REAL(wp) :: Glob_WorstEigvalTol  ! worst accuracy seen


  !================================================================
  ! Optimization history (one record per basis function)
  !================================================================

  TYPE :: Glob_HistoryStep
    REAL(wp) :: Energy
    INTEGER  :: CyclesDone
    INTEGER  :: InitFuncAtLastStep
    INTEGER  :: NumOfEnergyEvalDuringFullOpt
  END TYPE Glob_HistoryStep

  ! Per-function history of energy and basis-building/optimization
  TYPE(Glob_HistoryStep), ALLOCATABLE, DIMENSION(:), SAVE :: Glob_History


  !================================================================
  ! Masses
  !================================================================

  ! Masses of the particles (NOT pseudoparticles). Size = Glob_n + 1.
  REAL(wp), ALLOCATABLE, DIMENSION(:), SAVE :: Glob_Mass

  ! Total mass of the system
  REAL(wp) :: Glob_MassTotal

  ! .true. if masses of all pseudoparticles (particles 2..n+1) are equal.
  ! Enables simplifications in matrix-element evaluation.
  LOGICAL :: Glob_ArePseudoParticleMassesTheSame


  !================================================================
  ! Young operator (permutational symmetry)
  !================================================================

  ! Length of the Young-operator string
  INTEGER, PARAMETER :: Glob_YOperatorStringLength = 255

  ! Symbolic expression of the Young operator, read from input
  CHARACTER(Glob_YOperatorStringLength) :: Glob_YOperatorString


  !================================================================
  ! Random-generator parameters (used by GenerateTrialParam)
  !================================================================

  REAL(wp) :: Glob_RG_p1  ! distribution control
  REAL(wp) :: Glob_RG_s1  ! scale 1
  REAL(wp) :: Glob_RG_s2  ! scale 2

  INTEGER :: Glob_NPointsForRndTrF
  INTEGER :: Glob_RndTrNumArray(2, 100)

  INTEGER                                  :: Glob_TrialFuncGenMethod = 4       ! Controls distribution used by GenerateTrialParam
  INTEGER                                  :: Glob_MaxFuncEvalForCyclOpt = 100  ! Maximum function evaluations for cyclic optimization (set from field G of OPT_CYCLE)
  INTEGER, ALLOCATABLE, DIMENSION(:), SAVE :: Glob_NTrials                      ! Number of random trials per basis size


  !================================================================
  ! Basis Building and Optimization Program (BBOP)
  !================================================================

  ! Maximal length of file names used anywhere in the code
  INTEGER, PARAMETER :: Glob_FileNameLength = 70

  ! One step of the BBOP script
  TYPE :: Glob_BBOPStep
    CHARACTER(9)                   :: Action
    CHARACTER(1)                   :: GSEPSolutionMethod
    INTEGER                        :: A
    INTEGER                        :: B
    INTEGER                        :: C
    INTEGER                        :: D
    INTEGER                        :: E
    INTEGER                        :: F
    INTEGER                        :: G
    INTEGER                        :: H
    REAL(wp)                       :: Q
    REAL(wp)                       :: R
    CHARACTER(Glob_FileNameLength) :: FileName1
    CHARACTER(Glob_FileNameLength) :: FileName2
    CHARACTER(Glob_FileNameLength) :: FileName3
    CHARACTER(Glob_FileNameLength) :: FileName4
  END TYPE Glob_BBOPStep

  ! Number of BBOP steps
  INTEGER :: Glob_NumOfBBOPSteps

  ! The BBOP script itself
  TYPE(Glob_BBOPStep), ALLOCATABLE, DIMENSION(:), SAVE :: Glob_BBOP

  ! .true. if the BBOP contains any cyclic-optimization step
  LOGICAL :: Glob_IsOptCycleScripted = .FALSE.

  ! .true. while each process is working on its OWN trial function instead
  ! of all processes cooperating on one. Every collective in the energy
  ! evaluation path must be skipped while this is set, because the
  ! processes no longer execute the same sequence of operations - a
  ! collective left enabled here hangs the run. Set and cleared only by
  ! the stochastic-selection loops in BasisEnl_*; the matching linear
  ! algebra switch is linalg_setlocal.
  LOGICAL :: Glob_LocalWorkMode = .FALSE.


  !================================================================
  ! File names
  !================================================================

  ! Name of the input/output data file
  CHARACTER(Glob_FileNameLength) :: Glob_DataFileName = 'inout.txt'
  CHARACTER(Glob_FileNameLength) :: Glob_ErrMsgFileName = 'errormsg.txt'

  CHARACTER(Glob_FileNameLength) :: Glob_CorrFuncGridFileName = 'cfgrid.txt'
  CHARACTER(Glob_FileNameLength) :: Glob_CorrFuncFileName = 'corrfunc.txt'


  !================================================================
  ! MPI
  !================================================================

  ! MPI_WP, the MPI datatype corresponding to real(wp), is a PARAMETER of
  ! module wp_def (wp_def_$(PREC).f90).

  INTEGER :: Glob_ProcID      ! rank of this process
  INTEGER :: Glob_MPIErrCode  ! ierr returned by MPI calls


  !================================================================
  ! Eigenvalue solver: LAPACK tolerance
  !================================================================

  ! Used to tune routine DSYGVX (from LAPACK) accuracy
  ! Set .TRUE. to run a one-off finite-difference check of the overlap
  ! penalty gradient at the start of full optimization (diagnostic only).
  LOGICAL  :: Glob_CheckPenaltyGradient = .FALSE.
  REAL(wp) :: Glob_AbsTolForDSYGVX


  !================================================================
  ! Mass matrix and density-vector coefficients
  !================================================================

  ! Glob_MassMatrix is the mass matrix, M
  REAL(wp), ALLOCATABLE, DIMENSION(:, :), SAVE :: Glob_MassMatrix

  ! Vector Glob_bvc is used for computing particle densities. Its
  ! components depend on the masses of the particles.
  REAL(wp), ALLOCATABLE, DIMENSION(:, :), SAVE :: Glob_bvc


  !================================================================
  ! Pair-permutation (transposition) matrices
  !================================================================

  ! 4-D array holding every pair-permutation matrix:
  !   Glob_Transposit(1:Glob_n, 1:Glob_n, 1, 2)  is P12
  !   Glob_Transposit(1:Glob_n, 1:Glob_n, 5, 5)  is P55
  INTEGER, ALLOCATABLE, DIMENSION(:, :, :, :), SAVE :: Glob_Transposit


  !================================================================
  ! Young operator Y and Y^{+}Y
  !================================================================

  ! Number of independent terms in the Y and Y^{+}Y operators
  INTEGER :: Glob_NumYTerms
  INTEGER :: Glob_NumYHYTerms

  ! Permutation matrices for each independent term of Y and Y^{+}Y.
  !   Glob_YMatr(1:Glob_n, 1:Glob_n, k)    -- k-th term of Y
  !   Glob_YHYMatr(1:Glob_n, 1:Glob_n, k)  -- k-th term of Y^{+}Y
  REAL(wp), ALLOCATABLE, DIMENSION(:, :, :), SAVE :: Glob_YMatr
  REAL(wp), ALLOCATABLE, DIMENSION(:, :, :), SAVE :: Glob_YHYMatr

  ! Coefficients (signs / weights) of each permutation in Y and Y^{+}Y
  REAL(wp), ALLOCATABLE, DIMENSION(:), SAVE :: Glob_YCoeff
  REAL(wp), ALLOCATABLE, DIMENSION(:), SAVE :: Glob_YHYCoeff


  !================================================================
  ! Identical-particle and equivalent-pair bookkeeping
  !================================================================

  ! Number of identical-particle sets in the system
  INTEGER :: Glob_NumOfIdentPartSets

  ! Number of particles in each identical-particle set.
  ! Shape: (Glob_NumOfIdentPartSets)
  INTEGER, ALLOCATABLE, DIMENSION(:), SAVE :: Glob_NumOfPartInIdentPartSet

  ! List of identical particles (their numbers) in each set.
  ! Entries Glob_IdentPartList(1:Glob_NumOfPartInIdentPartSet(j), j)
  ! contain particle numbers that belong to set j.
  INTEGER, ALLOCATABLE, DIMENSION(:, :), SAVE :: Glob_IdentPartList

  ! Number of equivalent pseudoparticle-pair sets (note: j,j is
  ! also a pair even though it involves only pseudoparticle j).
  INTEGER :: Glob_NumOfNoneqvPairSets

  ! Number of pairs in each equivalent-pair set.
  ! Shape: (Glob_NumOfNoneqvPairSets)
  INTEGER, ALLOCATABLE, DIMENSION(:), SAVE :: Glob_NumOfPairsInEqvPairSet

  ! List of equivalent pairs. Entries
  !   Glob_EqvPairList(1:2, 1:Glob_NumOfPairsInEqvPairSet(j), j)
  ! contain the pairs belonging to set j. The first index (1 or 2)
  ! designates the first / second particle of the pair.
  INTEGER, ALLOCATABLE, DIMENSION(:, :, :), SAVE :: Glob_EqvPairList


  !================================================================
  !   External LAPACK function needed by ProgramDataInit
  !------------------------------------------------------------------
  !  Not a Glob_ variable, but ProgramDataInit calls DLAMCH('S')
  !  to set Glob_AbsTolForDSYGVX, so the external declaration
  !  must also be present somewhere in globvars.f90.
  !================================================================

  REAL(wp), EXTERNAL :: DLAMCH
  INTEGER, EXTERNAL  :: ILAENV


  !================================================================
  ! Eigenvalue solver: matrices and workspace
  !================================================================

  ! Method currently used to solve the GSEP.
  ! 'G' = DSYGVX (LAPACK),  'I' = inverse iteration,  'U' = undefined
  CHARACTER(1) :: Glob_GSEPSolutionMethod = 'U'

  ! Leading dimension of Glob_H and Glob_S (= max basis size)
  INTEGER :: Glob_HSLeadDim

  ! Hamiltonian matrix
  REAL(wp), ALLOCATABLE, DIMENSION(:, :), SAVE :: Glob_H

  ! Overlap matrix
  REAL(wp), ALLOCATABLE, DIMENSION(:, :), SAVE :: Glob_S

  ! Diagonal elements of the Hamiltonian matrix
  REAL(wp), ALLOCATABLE, DIMENSION(:), SAVE :: Glob_diagH

  ! Diagonal elements of the overlap matrix
  REAL(wp), ALLOCATABLE, DIMENSION(:), SAVE :: Glob_diagS

  ! Eigenvector
  REAL(wp), ALLOCATABLE, DIMENSION(:), SAVE :: Glob_c

  ! Inverse diagonal elements from Cholesky factorization of
  ! H - Glob_ApproxEnergy*S (used when GSEPSolutionMethod = 'I')
  REAL(wp), ALLOCATABLE, DIMENSION(:), SAVE :: Glob_invD

  ! Number of eigenvalues of the pencil (H,S) lying BELOW the inverse-iteration
  ! shift Glob_ApproxEnergy. Obtained for free from the LDL^T factorization of
  ! H-shift*S that GSEPIIS computes anyway: by Sylvester's law of inertia the
  ! number of negative diagonal entries of D equals the number of eigenvalues
  ! below the shift. Set by GSEPIIS after every solve; -1 if the factorization
  ! failed, in which case the inertia is undefined.
  INTEGER :: Glob_NumEvalsBelowShift = 0

  ! Index (1 = lowest) of the eigenvalue GSEPIIS actually returned, or -1
  ! when unknown. Inverse iteration converges to the eigenvalue nearest
  ! the shift; with m eigenvalues below the shift it is lambda_m or
  ! lambda_{m+1}, so the index is exact. Compared against
  ! Glob_WhichEigenvalue to notice that the INVIT path slid onto another
  ! state.
  INTEGER :: Glob_LastEigIndex = -1

  ! .TRUE. once the shift has been reported as targeting the wrong state, so
  ! the warning is printed once per run rather than once per energy evaluation.
  LOGICAL :: Glob_WrongStateReported = .FALSE.

  ! Previous eigenvector (initial guess for inverse iteration)
  REAL(wp), ALLOCATABLE, DIMENSION(:), SAVE :: Glob_LastEigvector

  ! Maximum number of iterations allowed in GSEPIIS
  INTEGER :: Glob_MaxIterForGSEPIIS = 30

  ! Largest GSEPIIS residual still accepted as a usable eigenpair. This is
  ! a safety net, not an accuracy target: the attainable residual is set
  ! by the conditioning of S and degrades with the basis size, so a value
  ! near EIGVAL_TOLERANCE turns normal ill-conditioning into a fatal
  ! error (1e-6 did). 1e-3 rejects only genuine garbage; set it huge for
  ! the reference behaviour without the test.
  REAL(wp) :: Glob_EigvalTolUsable = 1.0E-3_wp

  ! Maximum number of energy-evaluation failures allowed during
  ! optimization of nonlinear parameters
  INTEGER, PARAMETER :: Glob_MaxEnergyFailsAllowed = 5

  ! Workspace size for DSYGVX (precomputed as (NB+3)*HSLeadDim)
  INTEGER :: Glob_LWorkForDSYGVX

  ! Real workspace for DSYGVX (WORK, size Glob_LWorkForDSYGVX)
  REAL(wp), ALLOCATABLE, DIMENSION(:), SAVE :: Glob_WorkForDSYGVX

  ! Integer workspace for DSYGVX (IWORK, size 5*HSLeadDim)
  INTEGER, ALLOCATABLE, DIMENSION(:), SAVE :: Glob_IWorkForDSYGVX


  !================================================================
  ! Matrix element derivatives
  !================================================================

  ! Derivatives of H and S w.r.t. nonlinear parameters.
  !   Glob_D(1:np, i, j)       = dH_{i+nfru,j} / dvechL_{i+nfru}
  !   Glob_D(np+1:2*np, i, j)  = dS_{i+nfru,j} / dvechL_{i+nfru}
  ! Index i ranges from 1 to Glob_nfo, index j from 1 to Glob_nfa.
  REAL(wp), ALLOCATABLE, DIMENSION(:, :, :), SAVE :: Glob_D

  ! Derivatives contracted with the K-th function:
  !   Glob_dHxKda(1:nfa, 1:npt)  and  Glob_dSxKda(1:nfa, 1:npt)
  REAL(wp), ALLOCATABLE, DIMENSION(:, :), SAVE :: Glob_dHxKda, Glob_dSxKda

  ! Full derivative arrays for all-parameter optimization:
  !   Glob_dHda(1:N, 1:N, 1:npt)  and  Glob_dSda(1:N, 1:N, 1:npt)
  REAL(wp), ALLOCATABLE, DIMENSION(:, :, :), SAVE :: Glob_dHda, Glob_dSda


  !================================================================
  ! Temporary work arrays
  !================================================================

  REAL(wp), ALLOCATABLE, DIMENSION(:), SAVE :: Glob_WorkArrayR1
  REAL(wp), ALLOCATABLE, DIMENSION(:), SAVE :: Glob_WorkArrayR2
  REAL(wp), ALLOCATABLE, DIMENSION(:), SAVE :: Glob_WorkArrayR3
  INTEGER, ALLOCATABLE, DIMENSION(:), SAVE  :: Glob_WorkArrayI1
  INTEGER, ALLOCATABLE, DIMENSION(:), SAVE  :: Glob_WorkArrayI2
  INTEGER, ALLOCATABLE, DIMENSION(:), SAVE  :: Glob_IntWorkArrForSaveResults
  REAL(wp), ALLOCATABLE, DIMENSION(:), SAVE :: Glob_WkGR
  REAL(wp), ALLOCATABLE, DIMENSION(:), SAVE :: Glob_WorkForGSEPIIS
  ! Buffers used to store and send/receive the Hamiltonian and the
  ! overlap matrix elements. Length Glob_HSBuffLen.
  REAL(wp), ALLOCATABLE, DIMENSION(:), SAVE :: Glob_HklBuff1, Glob_HklBuff2
  REAL(wp), ALLOCATABLE, DIMENSION(:), SAVE :: Glob_SklBuff1, Glob_SklBuff2

  ! Buffers used to store and send/receive the DERIVATIVES of the
  ! Hamiltonian and the overlap matrix elements. Shape
  ! (2*Glob_npt, Glob_HSBuffLen): the k-index pair goes into the Dk
  ! buffers and the l-index pair into the Dl buffers, with 1 and 2
  ! being the send and receive halves of the exchange.
  REAL(wp), ALLOCATABLE, DIMENSION(:, :), SAVE :: Glob_DkBuff1, Glob_DkBuff2
  REAL(wp), ALLOCATABLE, DIMENSION(:, :), SAVE :: Glob_DlBuff1, Glob_DlBuff2


  !================================================================
  ! MPI communication buffers
  !================================================================

  ! Number of MPI processes
  INTEGER :: Glob_NumOfProcs

  ! Length of H/S communication buffers
  INTEGER :: Glob_HSBuffLen

  ! Buffers for sending/receiving H and S matrix elements
  REAL(wp), ALLOCATABLE, DIMENSION(:), SAVE :: Glob_HBuffer1, Glob_HBuffer2
  REAL(wp), ALLOCATABLE, DIMENSION(:), SAVE :: Glob_SBuffer1, Glob_SBuffer2

  ! Buffers for sending/receiving derivative matrix elements
  REAL(wp), ALLOCATABLE, DIMENSION(:, :), SAVE :: Glob_dHBuffer1, Glob_dSBuffer1
  REAL(wp), ALLOCATABLE, DIMENSION(:, :), SAVE :: Glob_dHBuffer2, Glob_dSBuffer2
  REAL(wp), ALLOCATABLE, DIMENSION(:, :), SAVE :: Glob_dHiBuffer1, Glob_dSiBuffer1
  REAL(wp), ALLOCATABLE, DIMENSION(:, :), SAVE :: Glob_dHiBuffer2, Glob_dSiBuffer2
  REAL(wp), ALLOCATABLE, DIMENSION(:, :), SAVE :: Glob_dHjBuffer1, Glob_dSjBuffer1
  REAL(wp), ALLOCATABLE, DIMENSION(:, :), SAVE :: Glob_dHjBuffer2, Glob_dSjBuffer2


  !================================================================
  ! Basis enlargement and optimization parameters
  !================================================================

  ! Basis size currently being attempted
  INTEGER :: Glob_nfa

  ! Number of functions being added or optimized simultaneously
  INTEGER :: Glob_nfo

  ! Number of functions that remain unchanged (= Glob_nfa - Glob_nfo)
  INTEGER :: Glob_nfru

  ! Number of the last blacklisted basis function (sorted).
  ! Defines the size of array Glob_Blacklisted.
  INTEGER :: Glob_lbf

  ! Specifies whether a particular function is excluded from cyclic
  ! optimization. Check i <= Glob_lbf before accessing.
  LOGICAL, ALLOCATABLE, DIMENSION(:), SAVE :: Glob_Blacklisted

  ! Scale factor applied to the gradient in optimization routines
  REAL(wp) :: Glob_GradientScaleFactor = ONE

  ! Maximum allowed fraction of trial failures during random
  ! selection (e.g. 0.15 means no more than 15% failures)
  ! NOTE the _wp suffix. Without it the literal is a DEFAULT REAL, i.e.
  ! SINGLE precision: 0.15 is rounded to 24 bits and only then widened to
  ! real(wp), so the constant is 0.150000005960464477539 rather than
  ! 0.15 - a relative error of 4e-8 in a build whose whole point is 1e-19.
  ! Every real literal in this file carries the suffix for that reason.
  REAL(wp), PARAMETER :: Glob_MaxFracOfTrialFailsAllowed = 0.15_wp

  ! How many times basis enlargement / cyclic optimization may
  ! repeat if the generated function ends up linearly dependent
  INTEGER, PARAMETER :: Glob_BadOverlapOrLinCoeffLim = 10

  ! Largest cancellation factor C = sum|c_k S_k| / |<phi|Y+Y|phi>| a NEW basis
  ! function may have. Above it the Young operator has almost annihilated the
  ! function and its normalized matrix elements have lost log10(C) digits
  ! (1e4: four digits). BasisEnlG/I redraw a block until every new function
  ! passes; SaveHSRaw flags such functions in the basis health table.
  REAL(wp), PARAMETER :: Glob_MaxSelfOverlapCancel = 1.0E+04_wp

  ! Maximum 2-norm allowed for D * step in the very first
  ! optimization step attempted by DRMNG
  REAL(wp), PARAMETER :: Glob_MaxScStepAllowedInOpt = 0.001_wp

  ! Overlap penalty control
  LOGICAL  :: Glob_OverlapPenaltyAllowed = .FALSE.
  REAL(wp) :: Glob_OverlapPenaltyThreshold2 = ONE
  REAL(wp) :: Glob_MaxOverlapPenalty = ONE
  REAL(wp) :: Glob_TotalOverlapPenalty = ZERO

  ! Relative threshold for scaling of nonlinear parameters
  REAL(wp), PARAMETER :: Glob_OptScalingThreshold = 1.0E+06_wp

  ! Print nonlinear parameters in CycleOptX routines?
  LOGICAL, PARAMETER :: Glob_AreParamPrintedInCycleOptX = .TRUE.

  ! Detail level for elimination-routine output
  INTEGER, PARAMETER :: Glob_ElimRoutPrintSpec = 2

  ! Maximal allowed values for Glob_np and Glob_npt
  INTEGER, PARAMETER :: Glob_np_MaxAllowed = Glob_n*(Glob_n+1)/2

  INTEGER, PARAMETER :: Glob_npt_MaxAllowed = Glob_np_MaxAllowed

  REAL(wp) :: Glob_OptTol = ZERO


  !================================================================
  ! Routine call counters and timing
  !================================================================

  ! Energy-routine call counters
  INTEGER :: Glob_EnergyGACounter = 0  ! EnergyGA calls
  INTEGER :: Glob_EnergyGBCounter = 0  ! EnergyGB calls
  INTEGER :: Glob_EnergyIACounter = 0  ! EnergyIA and EnergyIAM calls
  INTEGER :: Glob_EnergyIBCounter = 0  ! EnergyIB calls

  ! Temporary counters for computing average inverse-iteration counts
  INTEGER :: Glob_InvItTempCounter1 = 0
  INTEGER :: Glob_InvItTempCounter2 = 0

  ! Accumulated CPU time since program start
  REAL(wp) :: Glob_TimeSinceStart = 0.0_wp

  ! Minimum seconds between intermediate saves during full optimization
  REAL(wp) :: Glob_MinSavingIntervForFullOpt = 60.0_wp

  ! Function evaluations accumulated in the previous full optimization
  ! (carried over when full optimization is restarted)
  INTEGER :: Glob_NumOfFEDuringPrevFO = 0

  ! How many initial cyclic-optimization steps require mandatory saving
  INTEGER, PARAMETER :: Glob_MinMandSavSteps = 5

  ! Time of current execution
  REAL(wp) :: Glob_CurrExecTime


  !================================================================
  ! Additional file names and I/O controls
  !================================================================

  CHARACTER(Glob_FileNameLength) :: Glob_SwapFileName = 'swapfile.dat'
  CHARACTER(Glob_FileNameLength) :: Glob_ReallocFileName = 'REALloc.dat'
  CHARACTER(Glob_FileNameLength) :: Glob_BlackListFileName = 'blacklist.txt'
  CHARACTER(Glob_FileNameLength) :: Glob_ExpValFileName = 'expvals.txt'
  ! What the OVERLAP_D step prints is mirrored into this file, the way
  ! ExpcVals mirrors into expvals.txt. The screen gets a summary plus the
  ! extreme eigenvalues; the file gets the complete spectrum. This is only
  ! the DEFAULT: an OVERLAP_D line may name another file as its 4th field.
  CHARACTER(Glob_FileNameLength) :: Glob_OverlapFileName = 'overlap.txt'
  ! The per-function basis-health table printed by the SAVE_HS_R step is
  ! mirrored into this file, the way OVERLAP_D mirrors into overlap.txt.
  CHARACTER(Glob_FileNameLength)            :: Glob_HealthFileName = 'basis_health.txt'
  CHARACTER(Glob_FileNameLength), PARAMETER :: Glob_FileNameNone = 'none'

  LOGICAL :: Glob_UseSwapFile = .TRUE.
  LOGICAL :: Glob_UseReallocFile = .FALSE.

  ! Controls what data is saved in the Hessian file during full
  ! optimization. If .true., the scaling vector D is also saved.
  LOGICAL, PARAMETER :: Glob_FullOptSaveD = .TRUE.


  !================================================================
  ! Charge matrix and miscellaneous
  !================================================================

  ! Glob_ChargeMatrix is the charge-product matrix in the form matelem.f90
  ! reads: indexed 1..n, the diagonal holds q_i*q_0 (pseudoparticle i with the
  ! reference particle) and the strict lower triangle q_i*q_j. ProgramDataInit
  ! fills it from Glob_ScaledPseudoChargeMatrix; the old workproc.f90 fills it
  ! from Glob_PseudoCharge directly.
  REAL(wp), ALLOCATABLE, DIMENSION(:, :), SAVE :: Glob_ChargeMatrix

  ! Glob_PseudoChargeMatrix is the matrix of pseudocharge products, and
  ! Glob_ScaledPseudoChargeMatrix is its scaled version. The two differ only
  ! when the repulsion or attraction strengths are scaled, i.e. when the data
  ! file carried a REPULSION_SCALING_PARAM or ATTRACTION_SCALING_PARAM line.
  !
  ! Both are indexed 0..n, NOT 1..n: index 0 is the reference particle, so the
  ! reference-pseudoparticle products sit in row/column 0 instead of being
  ! folded into the diagonal the way Glob_ChargeMatrix does it.
  REAL(wp), ALLOCATABLE, DIMENSION(:, :), SAVE :: Glob_PseudoChargeMatrix
  REAL(wp), ALLOCATABLE, DIMENSION(:, :), SAVE :: Glob_ScaledPseudoChargeMatrix

  INTEGER :: Glob_CorrFuncNPoints
  LOGICAL :: Glob_IsCorrFuncNeeded = .FALSE.


  ! Current BBOP step
  INTEGER :: Glob_CurrBBOPStep = 0

  ! Glob_PseudoCharge0 is the charge of the reference particle, q0
  REAL(wp) :: Glob_PseudoCharge0

  ! Glob_PseudoCharge is the charges of pseudoparticles, qi
  REAL(wp), ALLOCATABLE, DIMENSION(:), SAVE :: Glob_PseudoCharge

  ! Glob_RepulsionScalingParam and Glob_AttractionScalingParam (may range from
  ! 0 to inf; default is 1) are parameters that change the repulsion and
  ! attraction strength between particles.
  ! Glob_RepulsionScalingParamPlus and Glob_RepulsionScalingParamMinus are
  ! additional scaling parameters that scale the repulsion between positive and
  ! between negative charges respectively.
  ! Each is set from its own optional line of the data file; the matching
  !*Supplied flag records whether that line was present, so that SaveResults
  ! writes back only the lines that were read.
  REAL(wp) :: Glob_RepulsionScalingParam = 1.0_wp
  REAL(wp) :: Glob_RepulsionScalingParamPlus = 1.0_wp
  REAL(wp) :: Glob_RepulsionScalingParamMinus = 1.0_wp
  REAL(wp) :: Glob_AttractionScalingParam = 1.0_wp
  LOGICAL  :: Glob_RepScalParamSupplied = .FALSE.
  LOGICAL  :: Glob_RepScalParamPlusSupplied = .FALSE.
  LOGICAL  :: Glob_RepScalParamMinusSupplied = .FALSE.
  LOGICAL  :: Glob_AttrScalParamSupplied = .FALSE.


  !================================================================
  ! correlation function
  !================================================================

  ! Grid points for the correlation function
  REAL(wp), ALLOCATABLE, DIMENSION(:), SAVE :: Glob_CorrFuncGrid


END MODULE globvars
