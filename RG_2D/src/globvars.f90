module globvars
!This module contains declarations of global variables
!and constants.
  use wp_def
  implicit none

!=============================================================
!Numerical constants Glob_MaxEnergyFailsAllowed Glob_MaxIterForGSEPIIS
!=============================================================
  real(wp),parameter :: &
    ZERO=0.E0_wp,     &
    ONE=1.E0_wp,      &
    TWO=2.E0_wp,      &
    THREE=3.E0_wp,    &
    FOUR=4.E0_wp,     &
    FIVE=5.E0_wp,     &
    SIX=6.E0_wp,      &
    SEVEN=7.E0_wp,    &
    EIGHT=8.E0_wp,    &
    NINE=9.E0_wp,     &
    TEN=10.0_wp,      &
    ONEHALF=ONE/TWO,     &
    ONETHIRD=ONE/THREE,  &
    ONEFOURTH=ONE/FOUR,  &
    ONEFIFTH=ONE/FIVE,  &
    THREEHALF=THREE/TWO, &
    Glob_Pi=3.1415926535897932384626433832795029E0_wp,     &
    Glob_SqrtPi=1.7724538509055160272981674833411452E0_wp, &
    Glob_FineStructConst=7.2973525643E-03_wp,         &  !CODATA 2022 value 0.0072973525643(11)
    Glob_EulerConst=0.57721566490153286060651209008240E0_wp

!=============================================================
!Global parameters
!=============================================================

!Maximal number of pseudoparticles allowed
  integer,parameter :: Glob_AllowedNumOfPseudoParticles= &
                       Glob_AllowedNumOfParticles-1

!=============================================================
!Global variables
!These should be set when the program starts
!=============================================================

  character(5),parameter  :: Glob_BasisType='RG_2D'  !Short name for the type of basis used
  logical  ::  Glob_BasisTypeSupplied=.false.

  integer       Glob_n          !Number of pseudoparticles

!np=n(n+1)/2 - number of independent parameters in a
!symmetric matrix of size (n x n)
  integer       Glob_np

!Total number of nonlinear parameters per basis function
!In case of real L=1 Gaussians Glob_npt=Glob_np
  integer       Glob_npt

!Glob_np_MaxAllowed and Glob_npt_MaxAllowed determine the
!maximal allowed values for Glob_np and Glob_npt
  integer,parameter :: Glob_np_MaxAllowed= &
                       Glob_AllowedNumOfPseudoParticles*(Glob_AllowedNumOfPseudoParticles+1)/2
  integer,parameter :: Glob_npt_MaxAllowed=Glob_np_MaxAllowed

  real(wp)   Glob_2Raised3n2  !2^(3n/2)
  real(wp)   Glob_PiRaised3n2 !pi^(3n/2)

!Glob_CurrBasisSize is a variable whose value equals the current
!size of the basis
  integer       Glob_CurrBasisSize

!Glob_CurrBasisSize is a variable whose value equals the current
!energy
  real(wp)       Glob_CurrEnergy

!Glob_WhichEigenvalue is a variable that specifies which eigenvalue
!needs to be found. It is used only when the general symmetric
!eigenvalue problem solution method is 'G'
  integer       Glob_WhichEigenvalue

!Glob_ApproxEnergy is a good approximation to the exact eigenvalue.
!It is used only when the general symmetric eigenvalue problem
!solution method is 'I'
  real(wp)       Glob_ApproxEnergy

!Glob_EigvalTol is a variable that specifies the accuracy of solving
!the eigenvalue problem. It is actually used only when the general
!symmetric eigenvalue problem solution method is 'I'
  real(wp)     Glob_EigvalTol

!Glob_InvItParameter is a variable that specifies the factor by which
!an approximate eigenvalue is multiplied when the general symmetric
!eigenvalue problem solution method is 'I'
  real(wp)     Glob_InvItParameter

!Glob_MassMatrix is the mass matrix, M
  real(wp),allocatable,dimension(:,:),save ::  Glob_MassMatrix

!Glob_Mass is the masses of particles (not pseudoparticles!), M_i
  real(wp),allocatable,dimension(:),save ::  Glob_Mass

!Glob_MassTotal is the total mass of the system (all particles)
  real(wp)   Glob_MassTotal

!Glob_dmva2 is a constant depending on the masses of particles,
!which is used in the evaluation of drachmanized mass-velocity correction
  real(wp)   Glob_dmva21, Glob_dmva22

!Glob_dmvB, Glob_dmvM, and Glob_dmvMB are constant diagonal matrices used in
!the evaluation of drachmanized mass-velocity correction. their elements depend
!of the masses of particles
  real(wp)   Glob_dmvM(Glob_AllowedNumOfPseudoParticles,Glob_AllowedNumOfPseudoParticles)
  real(wp)   Glob_dmvB(Glob_AllowedNumOfPseudoParticles,Glob_AllowedNumOfPseudoParticles)
  real(wp)   Glob_dmvMB(Glob_AllowedNumOfPseudoParticles,Glob_AllowedNumOfPseudoParticles)

!Glob_PseudoChargeMatrix is the matrix consisting of pseudocharge products
!Glob_ScaledPseudoChargeMatrix is the scaled version of Glob_PseudoChargeMatrix
!(they dffer only when the repulsion or attraction strengths are scaled)
  real(wp),allocatable,dimension(:,:),save ::  Glob_PseudoChargeMatrix
  real(wp),allocatable,dimension(:,:),save ::  Glob_ScaledPseudoChargeMatrix

!Glob_PseudoCharge is the charges of pseudoparticles, qi
  real(wp),allocatable,dimension(:),save ::  Glob_PseudoCharge

!Glob_PseudoCharge0 is the charge of the reference particle, q0
  real(wp)   Glob_PseudoCharge0

!Glob_RepulsionScalingParam and Glob_AttractionScalingParam
!(may range from 0 to inf; default is 1) are parameters
!that change the repulsion and attraction strength between particles
!Glob_RepulsionScalingParamPlus and Glob_RepulsionScalingParamMinus
!are additional scaling parameters that scale the repulsion between
!positive and negative charges.
  real(wp)  ::  Glob_RepulsionScalingParam=1.0_wp
  real(wp)  ::  Glob_RepulsionScalingParamPlus=1.0_wp
  real(wp)  ::  Glob_RepulsionScalingParamMinus=1.0_wp
  real(wp)  :: Glob_AttractionScalingParam=1.0_wp
  logical  ::  Glob_RepScalParamSupplied=.false.
  logical  ::  Glob_RepScalParamPlusSupplied=.false.
  logical  ::  Glob_RepScalParamMinusSupplied=.false.
  logical  ::  Glob_AttrScalParamSupplied=.false.
  logical  ::  Glob_ArePseudoParticleMassesTheSame = .true.

!Glob_YOperatorStringLength defines the length of string
!Glob_YOperatorString
  integer, parameter  ::  Glob_YOperatorStringLength=255

!Glob_YOperatorString is a string that contains the symbolic
!expression or the Young operator that is read from an
!input/output file
  character(Glob_YOperatorStringLength)  Glob_YOperatorString

!Glob_LastEigvalTol is a variable that stores the accuracy
!reached last time the eigenvalue problem solver was called.
!It is used only when the general symmetric eigenvalue problem
!solution method is 'I'
  real(wp)   Glob_LastEigvalTol

!Glob_BestEigvalTol is a variable that stores the best accuracy
!of solving the eigenvalue problem for a certain period of time.
!It is used only when the general symmetric eigenvalue problem
!solution method is 'I'
  real(wp)   Glob_BestEigvalTol

!Glob_WorstEigvalTol is a variable that stores the best accuracy
!of solving the eigenvalue problem for a certain period of time.
!It is used only when the general symmetric eigenvalue problem
!solution method is 'I'
  real(wp)   Glob_WorstEigvalTol

!Array Glob_NonlinParam contains the nonlinear parameters of basis
!functions (elements of the Cholesky matrix, L_k)
  real(wp),allocatable,dimension(:,:),save :: Glob_NonlinParam

!Array Glob_Index contains the indices of the x,y-premultipliers
!of the basis functions. The indices generally range from 1 to Glob_n
  integer,allocatable,dimension(:,:),save :: Glob_Index

!Logical array Glob_IsIndexFixed of size 2 specifies if each index of 
!basis functions is fixed (true) or can be varied (false).
  logical :: Glob_IsIndexFixed(2) = .false.

!If Glob_IsIndexFixed(i)=.true. then the corresponding fixed index value 
!is stored in Glob_IndexFixedValue(i). Otherwise, Glob_IndexFixedValue is not used.
  integer :: Glob_IndexFixedValue(2)

!Glob_VectorCouplingScheme defines the vector coupling scheme used in the process of basis generation. 
!It can take values 1 (two particles in p-states counple to a D-state), 2 (one particle in a d-state), 
!or 0 (any of the two previous scheme is allowed). The default value is 0.  
  integer :: Glob_VectorCouplingScheme=0

!Variable Glob_IsVectorCouplingSchemeProvided is .true. if the user has provided a specific vector 
!coupling scheme. Otherwise, it is .false.
   logical :: Glob_IsVectorCouplingSchemeProvided=.false.

!Array Glob_FuncNum contains the basis function numbers
  integer,allocatable,dimension(:),save :: Glob_FuncNum

!Definition of a new type, Glob_HistoryStep. The data of this
!is used to store the energy and other information obtained
!when the basis size was equal to certain number of functions.

  type Glob_HistoryStep
    real(wp)   Energy
    integer       CyclesDone
    integer       InitFuncAtLastStep
    integer       NumOfEnergyEvalDuringFullOpt
  endtype Glob_HistoryStep

!Glob_History is an array of type Glob_HistoryStep that contains
!information about the energy and basis building and optimization
  type(Glob_HistoryStep),allocatable,dimension(:),save :: Glob_History

!4D array Glob_Transposit contains all pair permutation matrices
!(transpositions). The structure is as follows:
!  Glob_Transposit(1:Glob_n,1:Glob_n,1,2) -- matrix that corresponds to P12
!  Glob_Transposit(1:Glob_n,1:Glob_n,5,5) -- matrix that corresponds to P55
  integer,allocatable,dimension(:,:,:,:),save  :: Glob_Transposit

!Variables Glob_NumYTerms and Glob_NumYHYTerms are the number of
!independent terms in the Y and Y^{\dagger}Y operator respectively
  integer Glob_NumYTerms
  integer Glob_NumYHYTerms

!3D arrays Glob_YMatr and Glob_YHYMatr contains all matrices for Y and
!Y^{\dagger}Y operators.
!The structure is as follows:
!   Glob_YHYMatr(1:Glob_n,1:Glob_n,5) is the matrix corresponding to the
!   5-th term of Y^{\dagger}Y operator
  real(wp),allocatable,dimension(:,:,:),save  :: Glob_YMatr
  real(wp),allocatable,dimension(:,:,:),save  :: Glob_YHYMatr

!Arrays Glob_YCoeff and Glob_YHYCoeff contain all coefficients
!(coefficients of permutations) in the Y and Y^{\dagger}Y operators
  real(wp),allocatable,dimension(:),save  :: Glob_YCoeff
  real(wp),allocatable,dimension(:),save  :: Glob_YHYCoeff

!Variable Glob_GSEPSolutionMethod defines which method
!is currently used to solve the general symmetric eigenvalue
!problem. Initially wet set it to 'U', which means undefined
  character(1)  :: Glob_GSEPSolutionMethod='U'

!Variable Glob_HSLeadDim is the size (leading dimension)
!of H and S matrices, as well as the leading dimension of
!the array that contains derivatives
  integer   Glob_HSLeadDim

!Array Glob_H is used to store the Hamiltonian matrix
  real(wp),allocatable,dimension(:,:),save  :: Glob_H

!Array Glob_S is used to store the overlap matrix
  real(wp),allocatable,dimension(:,:),save  :: Glob_S

!Array Glob_diagH is used to store the diagonal elements
!of the Hamiltonian matrix
  real(wp),allocatable,dimension(:),save :: Glob_diagH

!Array Glob_diagS is used to store the diagonal elements
!of the overlap matrix
  real(wp),allocatable,dimension(:),save :: Glob_diagS

!Array Glob_D is used to store the derivatives of the Hamiltonian
!and the overlap matrix elements. The first index represents
!nonlinear parameters. It ranges from 1 to 2*Glob_np
!if np=Glob_np then
! Glob_D(1:np,i,j)         is dHijdvechLi
! Glob_D(np+1:2*np,i,j)    is dSijdvechLi
!If only partial optimization is carried out, that is if only
!last several functions are optimized (= Glob_nfo) then index 2 (i)
!ranges from 1 to Glob_nfo although in fact it should have
!been equal i+Glob_nfru.  This is done to avoid allocating
!memory that will not be used. So
!Glob_D(1:Glob_np,i,j) is dH_{i+Glob_nfru,j}/dvechL_{i+Glob_nfru}
! ....
!Index j ranges from 1 to Glob_nfa
  real(wp),allocatable,dimension(:,:,:),save  :: Glob_D

!Array Glob_c is used to store the eigenvector
  real(wp),allocatable,dimension(:),save :: Glob_c

!Array Glob_invD is used only when Glob_GSEPSolutionMethod='I'.
!It contains inverse elements of matrix D when
!H-Glob_ApproxEnergy*S is factorized in the Cholesky form, L*D*LT.
  real(wp),allocatable,dimension(:),save :: Glob_invD

!Array Glob_LastEigvector contains the eigenvector obtained previous
!time when GSEP was solved using the inverse iteration method. It can
!be used as a initial guess for solving GSEP again, in which case
!fewer iterations may be needed.
  real(wp),allocatable,dimension(:),save :: Glob_LastEigvector

!Variable Glob_HSBuffLen defines the length of arrays
!Glob_HBuffer1,Glob_HBuffer2,Glob_SBuffer1, and Glob_SBuffer2
  integer  Glob_HSBuffLen

!These 2 arrays are used as buffers to store and send/receive
!the Hamiltonian and the overlap matrix elements
  real(wp),allocatable,dimension(:),save ::  &
    Glob_HklBuff1,Glob_HklBuff2,Glob_SklBuff1,Glob_SklBuff2

!This array is used as buffer to store and send/receive
!the derivatives of the Hamiltonian and the overlap matrix
!elements
  real(wp),allocatable,dimension(:,:),save :: &
    Glob_DkBuff1,Glob_DkBuff2,Glob_DlBuff1,Glob_DlBuff2

!Basis size currently being attempted
  integer       Glob_nfa

!Number of functions being added or optimized simultaneously
  integer       Glob_nfo

!This variable is introduced for the sake of convenience. It should
!be set to be Glob_nfa-Glob_nfo every time when Glob_nfa or Glob_nfo
!change their values. This is basically the number of functions to
!remain unchanged
  integer       Glob_nfru

!Glob_lbf is a variable whose value is equal to the number of the last
!blacklisted basis function (assuming they are sorted). Thus it defines the size
!of array Glob_Blacklisted
  integer  :: Glob_lbf

!Logical array Glob_Blacklisted contains entries that specify if a particular
!function is supposed to be optimized (this only concerns the cyclic optimization
!routines, not the full optimization ones). Note that the size of this array
!is generally smaller than the current basis size, Glob_CurrBasisSize, so
!before accessing Glob_Blacklisted(i) it should be checked whether i>=Glob_lbf.
  logical,allocatable,dimension(:),save :: Glob_Blacklisted

!Vector Glob_bvc is used for computing particle densities. Its components
!depend on the masses of particles
  real(wp),allocatable,dimension(:,:),save :: Glob_bvc

!This variable is used to tune routine DSYGVX (from LAPACK) accuracy
  real(wp)  Glob_AbsTolForDSYGVX

!Variable Glob_MaxIterForGSEPIIS set the limit of iterations allowed
!in GSEPIIS
  integer  ::  Glob_MaxIterForGSEPIIS=30

!This parameter defines what is the allowed fraction of
!Eigenvalue problem solution failures for random selection
!process. For example, when Glob_MaxFracOfTrialFailsAllowed=0.15
!then it means no more than 15% of failures allowed. If this
!number is exceeded at some step then the calculations are stopped.
  real(wp),parameter  ::  Glob_MaxFracOfTrialFailsAllowed=0.15_wp

!This parameter defines how many failures in energy or gradient
!evaluation are allowed during optimization of nonlinear parameters
  integer,parameter  :: Glob_MaxEnergyFailsAllowed=10 !

!This parameter defines how many times a basis enlargement
!routine or a cyclic optimization routine are allowed to repeat
!random trial and optimization process if the generated
!function/functions end up being linearly dependent with other
!functions in the basis (overlap is close to 1.0) or any linear
!parameters by magnitude are greater then threshold.
!If this limit is exceeded then the routine stops and shows an
!error message.
  integer,parameter  :: Glob_BadOverlapOrLinCoeffLim=10

!This parameter sets the maximum 2-norm allowed for D (Scaling vector)
!dot product the very first step that DRMNG attempts. Tuning this value
!is important as it may affect the optimization process very significantly
  real(wp),parameter :: Glob_MaxScStepAllowedInOpt=0.001_wp

!Glob_OverlapPenaltyAllowed is set to be .true. when overlaps are subject
!of penalty for exceeding some threshold
  logical  ::  Glob_OverlapPenaltyAllowed=.false.

!Glob_OverlapPenaltyThreshold2 defines the absolute square of the overlap threshold.
!In case when  Glob_OverlapPenaltyAllowed=.false. and Glob_OverlapPenaltyThreshold2<ONE
!a penalty may be added to the total energy if any of the pair overlaps
!exceeds Glob_OverlapPenaltyThreshold by magnitude. This global variable
!is changed depending on the user input.
  real(wp)  ::  Glob_OverlapPenaltyThreshold2=ONE

!Glob_MaxOverlapPenalty defines the scale of pair overlap penalties. This global
!variable is changed depending on the user input. In general, using a value that is
!of the same order of magnitude as the total energy is recommended.
  real(wp)  ::  Glob_MaxOverlapPenalty=ONE

!Glob_TotalOverlapPenalty is equal to the sum of all pair overlap penalties
  real(wp)  ::  Glob_TotalOverlapPenalty=ZERO

!This parameter sets a relative threshold for scaling
!of nonlinear parameters
  real(wp),parameter :: Glob_OptScalingThreshold=1.0E+06_wp

!This logical parameter defines if the nonlinear parameters of functions
!being optimized are printed on the screen in CycleOptX routines
  logical, parameter :: Glob_AreParamPrintedInCycleOptX=.true.

!This parameter defines the level of information printed
!in routines that eliminate basis functions
  integer,parameter  :: Glob_ElimRoutPrintSpec=2

!=============================================================
!The following variables are the parameters of the random generator
!that is used to generate new basis function candidates. For the
!description of these see comments in random generator routine -
!GenerateTrialParam
!=============================================================
  real(wp)  Glob_RG_p1
  real(wp)  Glob_RG_s1
  real(wp)  Glob_RG_s2

!=============================================================
!The following variables serve as counters. They are equal the
!number of times certain routines have been called.
!=============================================================
  integer  :: Glob_EnergyGACounter=0
  integer  :: Glob_EnergyGBCounter=0
  integer  :: Glob_EnergyIACounter=0
  integer  :: Glob_EnergyIBCounter=0
!Glob_InvItTempCounter1 and Glob_InvItTempCounter1 are used
!for calculating the average number of inverse iterations
!in certain cycles (when Glob_GSEPSolutionMethod='I')
  integer  :: Glob_InvItTempCounter1=0
  integer  :: Glob_InvItTempCounter2=0

!=============================================================
!These variables are used to measure time or to define certain
!time intervals.
!=============================================================
  real(4)  :: Glob_TimeSinceStart=0.0

!=============================================================
!The following arrays and variables provide temporary work space
!for different subroutines
!=============================================================
  integer,allocatable,dimension(:)  :: Glob_IntWorkArrForSaveResults
  integer                                      Glob_LWorkForDSYGVX
  real(wp),allocatable,dimension(:)  :: Glob_WorkForDSYGVX
  integer,allocatable,dimension(:)  :: Glob_IWorkForDSYGVX
  real(wp),allocatable,dimension(:)  :: Glob_WkGR
  real(wp),allocatable,dimension(:)  :: Glob_WorkForGSEPIIS
!=============================================================
!Data used to store basis building and optimization program
!=============================================================
!Variable Glob_FileNameLength defines the maximal
!length of file names used throughout code
  integer,parameter :: Glob_FileNameLength=70
!The definition of a new type, Glob_BBOPStep
  type Glob_BBOPStep
    character(9)  Action
    character(1)  GSEPSolutionMethod
    integer       A
    integer       B
    integer       C
    integer       D
    integer       E
    integer       F
    integer       G
    integer       H
    real(wp)   Q
    real(wp)   R
    character(Glob_FileNameLength) FileName1
    character(Glob_FileNameLength) FileName2
    character(Glob_FileNameLength) FileName3
    character(Glob_FileNameLength) FileName4
  endtype Glob_BBOPStep

!This is just a constant string that is used to fill the parameter
!list in some routines
  character(Glob_FileNameLength),parameter :: Glob_FileNameNone='none'

!Number of Basis Building and Optimization Program Steps
  integer    Glob_NumOfBBOPSteps

!Array Glob_BBOP contains the steps of Basis Building
!and Optimization Program
  type(Glob_BBOPStep),allocatable,dimension(:),save :: Glob_BBOP

!Current step of BBOP
  integer  :: Glob_CurrBBOPStep=0

!Variable Glob_IsOptCycleScripted is .true. if the BBOP contains cyclic optimization
!routines. Otherwise it is .false.
  logical :: Glob_IsOptCycleScripted=.false.

!=============================================================
!Global variables for working with files
!=============================================================

!Glob_DataFileName is the name of the input/output file
!Glob_Glob_BlackListFileName is the name of the file containing
!the list of functions that are not supposed to be optimized
!(this concerns only cyclic optimization, routines OptCycleG and OptCycleI).
  character(Glob_FileNameLength)  :: Glob_DataFileName='inout.txt'
  character(Glob_FileNameLength)  :: Glob_SwapFileName='swapfile.dat'
  character(Glob_FileNameLength)  :: Glob_ReallocFileName='realloc.dat'
  character(Glob_FileNameLength)  :: Glob_BlackListFileName='blacklist.txt'
  character(Glob_FileNameLength)  :: Glob_ExpValFileName='expvals.txt'

!These are the variables that specify whether swapping
!is allowed and whether a temporary file should be used
!for array reallocation (the latter reduces memory fragmentation)
  logical :: Glob_UseSwapFile=.true.
  logical :: Glob_UseReallocFile=.false.

!Parameter Glob_FullOptSaveD controls what data should be
!saved in Hessian file during Full Optimization.
!If Glob_FullOptSaveD=.false. then only the Hessian is saved
!If Glob_FullOptSaveD=.true. then besides the Hessian, the
!program also saves the scaling vector D (used by the optimization
!routine)
!In the second case, after restart of the full optimization
!the scaling coefficients will be exactly the same as they were
!at the very first start, while when Glob_FullOptSaveD=.false.
!these coefficients are set to some values that depend on the
!current values of nonlinear parameters.
  logical,parameter :: Glob_FullOptSaveD=.true.

!Parameter Glob_MinMandSavSteps defines how many initial steps
!in cyclic optimization routines (OptCycleG, OptCycleG) should
!be followed by mandatory data saving (even if the user specifies
!rare saving frequency in those routines)
  integer,parameter :: Glob_MinMandSavSteps=5

!=============================================================
!Global variables that contain information about
!identical particles in the system. This information might
!be used to automatically symmetrize the expectation values
!of two-particle operators. However, an additional symmetrization
!may be required in case when the system under consideration has
!more complicated symmetry of the Young operator than that
!dictated by permutational symmetry of particles
!=============================================================

!Glob_NumOfIdentPartSets is the number of identical particle sets
!in the system
  integer   Glob_NumOfIdentPartSets

!Array Glob_NumOfPartInIdentPartSet(1:Glob_NumOfIdentPartSets)
!contains the number of particles in each set of identical
!particles
  integer,allocatable,dimension(:),save :: Glob_NumOfPartInIdentPartSet

!Array Glob_IdentPartList contains the list of identical particles
!(their numbers) in each set of identical particles. For example,
!entries Glob_IdentPartList(1:Glob_NumOfPartInIdentPartSet(j),j)
!contain particle numbers that belong to set j.
  integer,allocatable,dimension(:,:),save :: Glob_IdentPartList

!Glob_NumOfEqvPairSets is the number of equivalent pseudoparticle
!pair sets in the system (note that j,j is also considered to be a pair
!even though it involves only pseudoparticle j)
  integer   Glob_NumOfNoneqvPairSets

!Array Glob_NumOfPairsInEqvPairSet(1:Glob_NumOfNoneqvPairSets) contains
!the number of pairs in each set of equivalent pairs
  integer,allocatable,dimension(:),save :: Glob_NumOfPairsInEqvPairSet

!Array Glob_EqvPairList contains the list of
!equivalent pairs in each set of equivalent pairs of particles.
!For example, entries
!Glob_IdentPartList(1:2,1:Glob_NumOfPairsInEqvPairSet(j),j)
!contain pairs that belong to set j. First index changes from 1 to 2
!designating first and second particle in the pair.
  integer,allocatable,dimension(:,:,:),save :: Glob_EqvPairList

!=============================================================
!Global variables used with MPI routines
!=============================================================

  integer     Glob_NumOfProcs !Number of parallel processes
  integer     Glob_ProcID     !The ID of a particular process (ranges
!from 0 to Glob_NumOfProcs-1
  integer     Glob_MPIErrCode !Error code for MPI routines

!=============================================================
!Declaring external functions
!=============================================================

  real(wp), external :: DLAMCH
  integer, external  :: ILAENV

end module globvars
