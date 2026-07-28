module globvars
!This module contains declarations of global variables
!and constants.
  use wp_def
  implicit none

!=============================================================
!Numerical constants
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

  integer       Glob_n          !Number of pseudoparticles

!np=n(n+1)/2 - number of independent parameters in a
!symmetric matrix of size (n x n)
  integer       Glob_np

!Total number of nonlinear parameters per basis function
!In cases of real L=0 or L=1 Gaussians Glob_npt=Glob_np
  integer       Glob_npt

!Glob_np_MaxAllowed and Glob_npt_MaxAllowed determine the
!maximal allowed values for Glob_np and Glob_npt
  integer,parameter :: Glob_np_MaxAllowed= &
                       Glob_AllowedNumOfPseudoParticles*(Glob_AllowedNumOfPseudoParticles+1)/2
  integer,parameter :: Glob_npt_MaxAllowed=Glob_np_MaxAllowed

  real(wp)   Glob_2Raised3n2  !2^(3n/2)
  real(wp)   Glob_PiRaised3n2 !pi^(3n/2)

!Glob_MassMatrix is the mass matrix, M
  real(wp),allocatable,dimension(:,:),save  ::  Glob_MassMatrix

!Glob_Mass is the masses of particles (not pseudoparticles!), M_i
  real(wp),allocatable,dimension(:),save  ::  Glob_Mass

!Glob_MassTotal is the total mass of the system (all particles)
  real(wp)   Glob_MassTotal

!Glob_dmva2 is a constant depending on the masses of particles,
!which is used in the evaluation of drachmanized mass-velocity correction
  real(wp)   Glob_dmva2

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
  real(wp)  ::  Glob_AttractionScalingParam=1.0_wp
  logical  ::  Glob_RepScalParamSupplied=.false.
  logical  ::  Glob_RepScalParamPlusSupplied=.false.
  logical  ::  Glob_RepScalParamMinusSupplied=.false.
  logical  ::  Glob_AttrScalParamSupplied=.false.

!Glob_YOperatorStringLength defines the length of string
!Glob_YOperatorString
  integer, parameter  ::  Glob_YOperatorStringLength=255

!Glob_YOperatorString is a string that contains the symbolic
!expression or the Young operator that is read from an
!input/output file
  character(Glob_YOperatorStringLength)  Glob_YOperatorString

!Array Glob_NonlinParam contains the nonlinear parameters of basis
!functions (elements of the Cholesky matrix, L_k)
  real(wp),allocatable,dimension(:,:),save  :: Glob_NonlinParam

!Array Glob_ZIndex contains the indices of the z-premultiplier
!of the basis functions. The indices generally range from 1 to Glob_n
  integer,allocatable,dimension(:),save  :: Glob_ZIndex
  integer,allocatable,dimension(:,:),save  :: Glob_Index0, Glob_Index1

!Array Glob_FuncNum contains the basis function numbers
  integer,allocatable,dimension(:),save  :: Glob_FuncNum

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

!Array Glob_H is used to store the Hamiltonian matrix
  real(wp),allocatable,dimension(:,:),save  :: Glob_H

!Array Glob_S is used to store the overlap matrix
  real(wp),allocatable,dimension(:,:),save  :: Glob_S

!Array Glob_diagH is used to store the diagonal elements
!of the Hamiltonian matrix
  real(wp),allocatable,dimension(:),save  :: Glob_diagH

!Array Glob_diagS is used to store the diagonal elements
!of the overlap matrix
  real(wp),allocatable,dimension(:),save  :: Glob_diagS,Glob_diagS0,Glob_diagS1

!Array Glob_c is used to store the eigenvector
  real(wp),allocatable,dimension(:),save  :: Glob_c

!Vector Glob_bvc is used for computing particle densities. Its components
!depend on the masses of particles
  real(wp),allocatable,dimension(:,:),save  :: Glob_bvc

!=============================================================
!These variables are used to measure time or to define certain
!time intervals.
!=============================================================
  real(4)  :: Glob_TimeSinceStart=0.0
!=============================================================
!                        New in RGL01
!=============================================================

!Variable Glob_FileNameLength defines the maximal
!length of file names used throughout code
  integer,parameter :: Glob_FileNameLength=70

!The definition of a new type, Glob_DRMCStep
! type Glob_DRMCStep
  ! character(9)  Action
  ! integer       A
  ! integer       B
  ! integer       C
  ! integer       D
  ! integer       E
  ! integer       F
  ! integer       G
  ! integer       H
  ! real(wp)   Q
  ! real(wp)   R
  ! character(Glob_FileNameLength) FileName1
  ! character(Glob_FileNameLength) FileName2
  ! character(Glob_FileNameLength) FileName3
  ! character(Glob_FileNameLength) FileName4
  ! character(Glob_FileNameLength) FileName5
  ! character(Glob_FileNameLength) FileName6
! endtype Glob_DRMCStep

! type(Glob_DRMCStep),allocatable,dimension(:),save :: Glob_DRMC

! !Number of Data Reading and Matrix Calculator Program Steps
! integer    Glob_NumOfDRMCSteps,Glob_CurrDRMCStep

!Glob_CurrBasisSize is a variable whose value equals the current
!size of the basis
integer  :: Glob_CurrBasisSize0
integer  :: Glob_CurrBasisSize1
!!!integer,allocatable,dimension(:),save  :: Glob_CurrBasisSizeInDRMCSteps

  real(wp)  :: Glob_ExpVals1, Glob_ExpVals2

! CURRENT_ENERGYs from wave function files
  real(wp)  :: Glob_CurrEnergy0,Glob_CurrEnergy1

!Array Glob_S is used to store the overlap matrices for L=0 and L=1 cases
  real(wp),allocatable,dimension(:,:),save  :: Glob_S0, Glob_S1

!Array Glob_diag_S is used to store the diagonal elements
!of the overlap matrix
  real(wp),allocatable,dimension(:),save  :: Glob_diag_S0,Glob_diag_S1

  character(Glob_YOperatorStringLength)  :: Glob_YOperatorString0, Glob_YOperatorString1

!Array Glob_c is used to store the eigenvector
  real(wp),allocatable,dimension(:),save  :: Glob_c0, Glob_c1

!Array Glob_FuncNum contains the basis function numbers
  integer,allocatable,dimension(:),save  :: Glob_FuncNum0, Glob_FuncNum1

!Array Glob_NonlinParam contains the nonlinear parameters of basis
!functions (elements of the Cholesky matrix, L_k)
  real(wp),allocatable,dimension(:,:),save  :: Glob_NonlinParam0, Glob_NonlinParam1

!Variables Glob_NumYTerms and Glob_NumYHYTerms are the number of
!independent terms in the Y operator respectively
  integer Glob_NumYTerms0, Glob_NumYTerms1
  integer Glob_NumYHYTerms0, Glob_NumYHYTerms1

!3D arrays Glob_YMatr  contains all matrices for Y operators.
!The structure is as follows:
!   Glob_YMatr(1:Glob_n,1:Glob_n,5) is the matrix corresponding to the
!   5-th term of Y operator.
  real(wp),allocatable,dimension(:,:,:),save  :: Glob_YMatr0, Glob_YMatr1
  real(wp),allocatable,dimension(:,:,:),save  :: Glob_YHYMatr0, Glob_YHYMatr1

!Arrays Glob_YCoeff contain all coefficients
!(coefficients of permutations) in the Y operators
  real(wp),allocatable,dimension(:),save  :: Glob_YCoeff0, Glob_YCoeff1
  real(wp),allocatable,dimension(:),save  :: Glob_YHYCoeff0, Glob_YHYCoeff1

  integer Glob_NumFactY0, Glob_NumFactY1
  integer,allocatable,dimension(:),save  :: Glob_NumTermsInYOpFact0, Glob_NumTermsInYOpFact1
  character(Glob_YOperatorStringLength),allocatable,dimension(:) :: Glob_YOpStr0, Glob_YOpStr1

!=============================================================

!This is just a constant string that is used to fill the parameter
!list in some routines
  character(Glob_FileNameLength),parameter  :: Glob_FileNameNone='none'

!=============================================================
!Global variables for working with files
!=============================================================

!Glob_DataFileName is the name of the input/output file
!Glob_Glob_BlackListFileName is the name of the file containing
!the list of functions that are not supposed to be optimized
!(this concerns only cyclic optimization, routines OptCycleG and OptCycleI).
  character(Glob_FileNameLength)  :: Glob_DataFileName='inout.txt'
  character(Glob_FileNameLength)  :: Glob_WFfile0='wf_state0.txt'
  character(Glob_FileNameLength)  :: Glob_WFfile1='wf_state1.txt'
  character(Glob_FileNameLength)  :: Glob_ExpValFileName='expvals.txt'
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
  integer,allocatable,dimension(:),save  :: Glob_NumOfPartInIdentPartSet

!Array Glob_IdentPartList contains the list of identical particles
!(their numbers) in each set of identical particles. For example,
!entries Glob_IdentPartList(1:Glob_NumOfPartInIdentPartSet(j),j)
!contain particle numbers that belong to set j.
  integer,allocatable,dimension(:,:),save  :: Glob_IdentPartList

!Glob_NumOfEqvPairSets is the number of equivalent pseudoparticle
!pair sets in the system (note that j,j is also considered to be a pair
!even though it involves only pseudoparticle j)
  integer   Glob_NumOfNoneqvPairSets

!Array Glob_NumOfPairsInEqvPairSet(1:Glob_NumOfNoneqvPairSets) contains
!the number of pairs in each set of equivalent pairs
  integer,allocatable,dimension(:),save  :: Glob_NumOfPairsInEqvPairSet

!Array Glob_EqvPairList contains the list of
!equivalent pairs in each set of equivalent pairs of particles.
!For example, entries
!Glob_IdentPartList(1:2,1:Glob_NumOfPairsInEqvPairSet(j),j)
!contain pairs that belong to set j. First index changes from 1 to 2
!designating first and second particle in the pair.
  integer,allocatable,dimension(:,:,:),save  :: Glob_EqvPairList

!=============================================================
!Global variables used with MPI routines
!=============================================================

  integer     Glob_NumOfProcs !Number of parallel processes
  integer     Glob_ProcID     !The ID of a particular process (ranges
!from 0 to Glob_NumOfProcs-1
  integer     Glob_MPIErrCode !Error code for MPI routines

  integer :: Glob_selectTransition !Select the transition for which the transition amplitude is computed
  logical :: Glob_spinFileWasOpened=.false. !Whether the file containing spin functions was opened successfully

end module globvars

