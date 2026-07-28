module matelem
!Module matelem contains subroutines for computing
!matrix elements with real L=0 and L=1 Gaussians.
  use globvars
  implicit none

contains

  subroutine OverlapMatrixElement_RG_2D(m_k, mm_k, vechLk, P, Skk)
    !This subroutine computes symmetry adapted matrix element with
    !two real L=2 correlated Gaussians:
    !
    !fk = (v'_k * r) (w'_k * r) exp[-r'(Lk*Lk')r]
    !
    !m_k and mm_k are integers between 1 and n (n is the number of
    !pseudoparticles). Symmetry adaption is applied to the ket using
    !permutation matrix P
    !
    !Input:
    !   m_k, mm_k :: integers that determine which pseudoparticles carry l=1 momentum
    !   vechLk :: Array of length (n(n+1)/2) of
    !     exponential parameters.
    !   P  :: The symmetry permutation matrix of size n x n
    !Output:
    !   Skk         ::        Overlap matrix element (normalized)

    !Arguments
    integer,intent(in)          :: m_k, mm_k
    real(wp),intent(in)      :: vechLk(Glob_np)
    real(wp),intent(in)      :: P(Glob_n,Glob_n)
    real(wp),intent(out)     :: Skk

    !Parameters (These are needed to declare static arrays. Using static
    !arrays makes the function call a little faster in comparison with
    !the case when arrays are dynamically allocated in stack)
    integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles
    integer,parameter :: nnp=nn*(nn+1)/2

    !Local variables
    integer           n, np
    real(wp)       vl(nn), bl(nn), vkinv_tAkl(nn), bkinv_tAkl(nn)
    real(wp)       Lk(nn,nn),Ll(nn,nn)
    real(wp)       Ak(nn,nn),tAl(nn,nn),tAkl(nn,nn)
    real(wp)       inv_tAkl(nn,nn)
    real(wp)       W1(nn,nn)
    real(wp)       temp1, temp2
    real(wp)       det_tAkl
    real(wp)       tau3, tau33, tau333, tau334, m, m1, m3
    integer           i,j,k,q,t,indx

    n=Glob_n
    np=Glob_np
    !First we build matrices Lk, Ll, Ak, Al from vechLk, vechLl.
    indx=0
    do i=1,n
      do j=i,n
        indx=indx+1
        Lk(i,j)=ZERO
        Lk(j,i)=vechLk(indx)
        Ll(i,j)=ZERO
        Ll(j,i)=vechLk(indx)
      enddo
    enddo

    do i=1,n
      do j=i,n
        temp1=ZERO
        do k=1,i
          temp1=temp1+Lk(i,k)*Lk(j,k)
        enddo
        Ak(i,j)=temp1
        Ak(j,i)=temp1
        temp1=ZERO
        do k=1,i
          temp1=temp1+Ll(i,k)*Ll(j,k)
        enddo
        tAl(i,j)=temp1
        tAl(j,i)=temp1
      enddo
    enddo

    !Then we permute elements of Al to account for
    !the action of the permutation matrix
    !tAl=P'*Al*P
    !We also form matrix tAkl=Ak+tAl
    do i=1,n
      do j=1,n
        temp1=ZERO
        do k=1,n
          temp1=temp1+P(k,j)*tAl(k,i)
        enddo
        W1(j,i)=temp1
      enddo
    enddo
    do i=1,n
      do j=i,n
        temp1=ZERO
        do k=1,n
          temp1=temp1+W1(i,k)*P(k,j)
        enddo
        tAl(i,j)=temp1
        tAl(j,i)=temp1
        tAkl(i,j)=Ak(i,j)+temp1
        tAkl(j,i)=tAkl(i,j)
      enddo
    enddo

    !After this we can do Cholesky factorization of tAkl.
    !The Cholesky factor will be temporarily stored in the
    !lower triangle of W1
    det_tAkl=ONE
    !temp1=ZERO
    do i=1,n
      do j=i,n
        temp1=tAkl(i,j)
        do k=i-1,1,-1
          temp1=temp1-W1(i,k)*W1(j,k)
        enddo
        if (i==j) then
          W1(i,i)=sqrt(temp1)
          det_tAkl=det_tAkl*temp1
        else
          W1(j,i)=temp1/W1(i,i)
          W1(i,j)=ZERO
        endif
      enddo
    enddo

    !Inverting tAkl using its Cholesky factor (stored in W1)
    !and placing the result into inv_tAkl
    do i=1,n
      W1(i,i)=ONE/W1(i,i)
      do j=i+1,n
        temp1=ZERO
        do k=i,j-1
          temp1=temp1-W1(j,k)*W1(k,i)
        enddo
        W1(j,i)=temp1/W1(j,j)
      enddo
    enddo

    do i=1,n
      do j=i,n
        temp1=ZERO
        do k=j,n
          temp1=temp1+W1(k,i)*W1(k,j)
        enddo
        inv_tAkl(i,j)=temp1
        inv_tAkl(j,i)=temp1
      enddo
    enddo

    !Computing vl=P'*vl, bl=P'*bl
    do i=1,n
      vl(i)=P(m_k,i)
      bl(i)=P(mm_k,i)
    enddo
    !Compute vkinv_tAkl=vk'*inv_tAkl, bkinv_tAkl=bk'*inv_tAkl
    do i=1,n
      vkinv_tAkl(i)=inv_tAkl(m_k,i)
      bkinv_tAkl(i)=inv_tAkl(mm_k,i)
    enddo

    !Compute tau3=vkinv_tAkl*vl, tau33=bkinv_tAkl*bl
    tau3=ZERO
    tau33=ZERO
    tau333=ZERO
    tau334=ZERO
    do i=1,n
      tau3=tau3+vkinv_tAkl(i)*vl(i)
      tau33=tau33+bkinv_tAkl(i)*bl(i)
      tau333=tau333+vkinv_tAkl(i)*bl(i)
      tau334=tau334+bkinv_tAkl(i)*vl(i)
    enddo
    m1=tau3*tau33
    m3=tau333*tau334
    m=m1+m3  !look formula 40 in document

    !Evaluating overlap
    !temp1=ZERO
    temp1=FOUR*det_tAkl*sqrt(det_tAkl)
    Skk=Glob_PiRaised3n2*m/temp1

  end subroutine OverlapMatrixElement_RG_2D

  subroutine OverlapMatrixElement_RG_2P(m_k, mm_k, vechLk, P, Skk)
    !This subroutine computes symmetry adapted matrix element with
    !two real L=2 correlated Gaussians:
    !
    !fk = (v'_k * r) (w'_k * r) exp[-r'(Lk*Lk')r]
    !
    !m_k and mm_k are integers between 1 and n (n is the number of
    !pseudoparticles). Symmetry adaption is applied to the ket using
    !permutation matrix P
    !
    !Input:
    !   m_k, mm_k :: integers that determine which pseudoparticles carry l=1 momentum
    !   vechLk :: Array of length (n(n+1)/2) of
    !     exponential parameters.
    !   P  :: The symmetry permutation matrix of size n x n
    !Output:
    !   Skk         ::        Overlap matrix element (normalized)

    !Arguments
    integer,intent(in)          :: m_k, mm_k
    real(wp),intent(in)      :: vechLk(Glob_np)
    real(wp),intent(in)      :: P(Glob_n,Glob_n)
    real(wp),intent(out)     :: Skk

    !Parameters (These are needed to declare static arrays. Using static
    !arrays makes the function call a little faster in comparison with
    !the case when arrays are dynamically allocated in stack)
    integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles
    integer,parameter :: nnp=nn*(nn+1)/2

    !Local variables
    integer           n, np
    real(wp)       vl(nn), bl(nn), vkinv_tAkl(nn), bkinv_tAkl(nn)
    real(wp)       Lk(nn,nn),Ll(nn,nn)
    real(wp)       Ak(nn,nn),tAl(nn,nn),tAkl(nn,nn)
    real(wp)       inv_tAkl(nn,nn)
    real(wp)       W1(nn,nn)
    real(wp)       temp1, temp2
    real(wp)       det_tAkl
    real(wp)       tau3, tau33, tau333, tau334, m, m1, m3
    integer           i,j,k,q,t,indx
    n=Glob_n
    np=Glob_np
    !First we build matrices Lk, Ll, Ak, Al from vechLk, vechLl.
    indx=0
    do i=1,n
      do j=i,n
        indx=indx+1
        Lk(i,j)=ZERO
        Lk(j,i)=vechLk(indx)
        Ll(i,j)=ZERO
        Ll(j,i)=vechLk(indx)
      enddo
    enddo

    do i=1,n
      do j=i,n
        temp1=ZERO
        do k=1,i
          temp1=temp1+Lk(i,k)*Lk(j,k)
        enddo
        Ak(i,j)=temp1
        Ak(j,i)=temp1
        temp1=ZERO
        do k=1,i
          temp1=temp1+Ll(i,k)*Ll(j,k)
        enddo
        tAl(i,j)=temp1
        tAl(j,i)=temp1
      enddo
    enddo

    !Then we permute elements of Al to account for
    !the action of the permutation matrix
    !tAl=P'*Al*P
    !We also form matrix tAkl=Ak+tAl
    do i=1,n
      do j=1,n
        temp1=ZERO
        do k=1,n
          temp1=temp1+P(k,j)*tAl(k,i)
        enddo
        W1(j,i)=temp1
      enddo
    enddo
    do i=1,n
      do j=i,n
        temp1=ZERO
        do k=1,n
          temp1=temp1+W1(i,k)*P(k,j)
        enddo
        tAl(i,j)=temp1
        tAl(j,i)=temp1
        tAkl(i,j)=Ak(i,j)+temp1
        tAkl(j,i)=tAkl(i,j)
      enddo
    enddo

    !After this we can do Cholesky factorization of tAkl.
    !The Cholesky factor will be temporarily stored in the
    !lower triangle of W1
    det_tAkl=ONE
    !temp1=ZERO
    do i=1,n
      do j=i,n
        temp1=tAkl(i,j)
        do k=i-1,1,-1
          temp1=temp1-W1(i,k)*W1(j,k)
        enddo
        if (i==j) then
          W1(i,i)=sqrt(temp1)
          det_tAkl=det_tAkl*temp1
        else
          W1(j,i)=temp1/W1(i,i)
          W1(i,j)=ZERO
        endif
      enddo
    enddo

    !Inverting tAkl using its Cholesky factor (stored in W1)
    !and placing the result into inv_tAkl
    do i=1,n
      W1(i,i)=ONE/W1(i,i)
      do j=i+1,n
        temp1=ZERO
        do k=i,j-1
          temp1=temp1-W1(j,k)*W1(k,i)
        enddo
        W1(j,i)=temp1/W1(j,j)
      enddo
    enddo

    do i=1,n
      do j=i,n
        temp1=ZERO
        do k=j,n
          temp1=temp1+W1(k,i)*W1(k,j)
        enddo
        inv_tAkl(i,j)=temp1
        inv_tAkl(j,i)=temp1
      enddo
    enddo

    !Computing vl=P'*vl, bl=P'*bl
    do i=1,n
      vl(i)=P(m_k,i)
      bl(i)=P(mm_k,i)
    enddo

    !Compute vkinv_tAkl=vk'*inv_tAkl, bkinv_tAkl=bk'*inv_tAkl
    do i=1,n
      vkinv_tAkl(i)=inv_tAkl(m_k,i)
      bkinv_tAkl(i)=inv_tAkl(mm_k,i)
    enddo

    !Compute tau3=vkinv_tAkl*vl, tau33=bkinv_tAkl*bl
    tau3=ZERO
    tau33=ZERO
    tau333=ZERO
    tau334=ZERO
    do i=1,n
      tau3=tau3+vkinv_tAkl(i)*vl(i)
      tau33=tau33+bkinv_tAkl(i)*bl(i)
      tau333=tau333+vkinv_tAkl(i)*bl(i)
      tau334=tau334+bkinv_tAkl(i)*vl(i)
    enddo
    m1=tau3*tau33
    m3=tau333*tau334
    m=m1-m3  !look formula 40 in document

    !Evaluating overlap
    !temp1=ZERO
    temp1=FOUR*det_tAkl*sqrt(det_tAkl)
    Skk=Glob_PiRaised3n2*m/temp1

  end subroutine OverlapMatrixElement_RG_2P

  subroutine spinPreCalc(n, nFactorial, parityFactor, SSNCmassChargeCoefficient, SOmassChargeCoefficient, &
              AMMmassChargeCoefficient, AMMFinmassChargeCoefficient, ketMatrix, &
              spatialYoung0, spatialYoung1, SSNCspinME, SiMinusME, SiPlusME, SziME, spinFreeME)
    use spinStuff
    implicit none

    !input vars:
    character(len = maxLen), intent(in) :: spatialYoung0, spatialYoung1
    integer, intent(in) :: n, nFactorial

    !output vars:
    real(wp), dimension(nFactorial), intent(out) :: parityFactor
    real(wp), dimension(n, n, 4), intent(out) :: SOmassChargeCoefficient, AMMmassChargeCoefficient, AMMFinMassChargeCoefficient
    real(wp), dimension(n, n), intent(out) :: SSNCmassChargeCoefficient
    real(wp), dimension(n, n, nFactorial), intent(out) :: ketMatrix
    real(wp), dimension(n, n, nFactorial), intent(out) :: SSNCspinME
    real(wp), dimension(nFactorial, 2), intent(out) :: spinFreeME
    real(kind = wp), dimension(n, nFactorial), intent(out) :: SiMinusME, SiPlusME, SziME

    ! local variables
    integer :: i, j, k, l, m
    character(len = maxLen) :: mySpatialYoung0, mySpatialYoung1
    integer, dimension(nFactorial) :: parities
    integer, dimension(n, n, nFactorial) :: allPermutations
    real(kind = wp), dimension(:), allocatable :: finalSpinFunction0, finalSpinFunction1
    integer, dimension(:, :), allocatable :: primitives0, primitives1
    integer :: numberOfPrimitives0, numberOfPrimitives1

    SOmassChargeCoefficient = ZERO
    SSNCmassChargeCoefficient = ZERO
    do i = 1, n
      SOmassChargeCoefficient(i, i, 1) = -ONEHALF * Glob_PseudoCharge0 * Glob_PseudoCharge(i) / Glob_Mass(i + 1) * &
                                         (ONE / Glob_Mass(i + 1) + TWO / Glob_Mass(1))
    enddo

    do i = 1, n
      SOmassChargeCoefficient(i, i, 2) = -Glob_PseudoCharge0 * Glob_PseudoCharge(i) / &
                                         (Glob_Mass(i + 1) * Glob_Mass(1))
    enddo

    do i = 1, n
      do j = 1, n
        SOmassChargeCoefficient(i, j, 3) = -Glob_PseudoCharge(i) * Glob_PseudoCharge(j) / &
                                           (Glob_Mass(i + 1) * Glob_Mass(j + 1))
        SSNCmassChargeCoefficient(i, j) = Glob_PseudoCharge(i) * Glob_PseudoCharge(j) / &
                                          (Glob_Mass(i + 1) * Glob_Mass(j + 1))
      enddo
    enddo

    do i = 1, n
      do j = 1, n
        SOmassChargeCoefficient(i, j, 4) = -ONEHALF * Glob_PseudoCharge(i) * Glob_PseudoCharge(j) / &
                                           (Glob_Mass(i + 1)**TWO)
      enddo
    enddo

     AMMmassChargeCoefficient = ZERO
    AMMFinmassChargeCoefficient = ZERO
    do i = 1, n
      AMMmassChargeCoefficient(i, i, 1) = -ONEHALF * Glob_PseudoCharge0 * Glob_PseudoCharge(i) / Glob_Mass(i + 1)**TWO
      AMMFinmassChargeCoefficient(i, i, 1) = -ONEHALF * Glob_PseudoCharge0 * Glob_PseudoCharge(i) / Glob_Mass(i + 1) * &
                                             (ONE / Glob_Mass(1) + ONE / Glob_Mass(i + 1))
    enddo

    do i = 1, n
      AMMFinmassChargeCoefficient(i, i, 2) = -ONEHALF * Glob_PseudoCharge0 * Glob_PseudoCharge(i) / &
                                             (Glob_Mass(i + 1) * Glob_Mass(1))
    enddo

    do i = 1, n
      do j = 1, n
        AMMmassChargeCoefficient(i, j, 3) = -ONEHALF * Glob_PseudoCharge(i) * Glob_PseudoCharge(j) / &
                                            (Glob_Mass(i + 1) * Glob_Mass(j + 1))
        AMMFinmassChargeCoefficient(i, j, 3) = AMMmassChargeCoefficient(i, j, 3)
      enddo
    enddo

    do i = 1, n
      do j = 1, n
        AMMmassChargeCoefficient(i, j, 4) = -ONEHALF * Glob_PseudoCharge(i) * Glob_PseudoCharge(j) / &
                                            (Glob_Mass(i + 1)**TWO)
        AMMFinmassChargeCoefficient(i, j, 4) = AMMmassChargeCoefficient(i, j, 4)
      enddo
    enddo

    ! now we deal with the spin stuff
    ! rename the particles
    mySpatialYoung0 = spatialYoung0
    do i = 1, maxLen
      if (mySpatialYoung0(i:i) == 'P') then
        read(mySpatialYoung0(i + 1:i + 1), *) k
        read(mySpatialYoung0(i + 2:i + 2), *) j
        write(mySpatialYoung0(i + 1:i + 1), '(i1)') k - 1
        write(mySpatialYoung0(i + 2:i + 2), '(i1)') j - 1
      endif
    enddo

    mySpatialYoung1 = spatialYoung1
    do i = 1, maxLen
      if (mySpatialYoung1(i:i) == 'P') then
        read(mySpatialYoung1(i + 1:i + 1), *) k
        read(mySpatialYoung1(i + 2:i + 2), *) j
        write(mySpatialYoung1(i + 1:i + 1), '(i1)') k - 1
        write(mySpatialYoung1(i + 2:i + 2), '(i1)') j - 1
      endif
    enddo

    call generatePermutationMatrices(allPermutations, n, nFactorial, parities)

    !Singlet wf for 1,2 transitions, triplet for 3
    call getSpinFunction(n, nFactorial, mySpatialYoung1, allPermutations, parities, &
                         finalSpinFunction0, primitives0, numberOfPrimitives0)

    !Triplet wf for 1,2 transitions, singlet for 3
    call getSpinFunction(n, nFactorial, mySpatialYoung0, allPermutations, parities, &
                         finalSpinFunction1, primitives1, numberOfPrimitives1)

    call getSpinOpMeanValues(n, nFactorial, allPermutations, finalSpinFunction0, finalSpinFunction1, primitives0, primitives1, &
                             numberOfPrimitives0, numberOfPrimitives1, spinFreeME, SSNCspinME, SiMinusME, SiPlusME, SziME)

    ketMatrix = ZERO
    do i = 1, nFactorial
      do k = 1, n
        do l = 1, n
          ! note the transposition here
          ketMatrix(k, l, i) = real(allPermutations(l, k, i), kind=wp)
        enddo
      enddo

    enddo

    do i = 1, nFactorial
      parityFactor(i) = real(parities(i), kind=wp)
    enddo

  end subroutine spinPreCalc

  subroutine spinDependentMatrixElements(m_k, m_l, mm_k, mm_l, vechLk, vechLl, Pket, &
                                         SOspinME, SSNCspinME, SSNCmassChargeCoefficient, SOmassChargeCoefficient, &
                                         AMMmassChargeCoefficient, AMMFinmassChargeCoefficient, &
                                         SSNCkl, SO1kl, SO2kl, AMM1kl, AMM2kl, AMM1Finkl, AMM2Finkl)
    !This subroutine computes symmetry adapted matrix element
    !with two real L=1 correlated Gaussians. These matrix element
    !is used in calculations of expectation values.

    !Input:
    !   m_k,m_l,mm_k, mm_l :: integers that determine which x or y-components is in the
    !                premultiplier of the Gaussian
    !   vechLk, vechLl :: Arrays of length (n(n+1)/2) of exponential parameters.

    !Output (all matrix elements are computed with normalized functions):

    !   SSNCkl :: Non-contact spin-spin term (without the factor of alpha**2)
    !   SO1kl, SO2kl  :: Spin-Orbit corrections (without the factor of alpha**2)
    !         1 and 2 stay for spin-same orbit and spin-another orbit contributions
    !   AMM1kl, AMM2kl  :: AMM corrections (without the factor of alpha**2)
    !         1 and 2 stay for spin-same orbit and spin-another orbit contributions

    !Arguments
    integer,intent(in)       :: m_k, m_l, mm_k, mm_l
    real(wp),intent(in)   :: vechLk(Glob_np), vechLl(Glob_np)
    real(wp),intent(in)   :: Pket(Glob_n,Glob_n)

    real(wp), intent(out)  :: SO1kl, SO2kl, AMM1kl, AMM2kl, AMM1Finkl, AMM2Finkl, SSNCkl
    !Parameters (These are needed to declare static arrays. Using static
    !arrays makes the function call a little faster in comparison with
    !the case when arrays are dynamically allocated in stack)
    integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles
    integer,parameter :: nnp=nn*(nn+1)/2
    real(wp),intent(in)   :: SSNCspinME(Glob_n, Glob_n), &
                                SOspinME(Glob_n), &
                                SOmassChargeCoefficient(Glob_n, Glob_n, 4), &
                                AMMmassChargeCoefficient(Glob_n, Glob_n, 4), &
                                AMMFinmassChargeCoefficient(Glob_n, Glob_n, 4), &
                                SSNCmassChargeCoefficient(Glob_n, Glob_n)

    !Local variables
    integer           n, np
    integer           tvk(nn),tvl(nn)
    real(wp)       Lk(nn,nn),Ll(nn,nn),inv_Lk(nn,nn),inv_Ll(nn,nn)
    real(wp)       tAk(nn,nn),tAl(nn,nn),tAkl(nn,nn)
    real(wp)       inv_tAkl(nn,nn)

    real(wp)       W1(nn,nn)
    real(wp)       temp1, temp2, temp3, temp1010, temp1001, temp0110, temp0101, t1, t2, det_tAkl
    integer :: i, j, k, indx

    integer :: pm_k, pm_l, pmm_k, pmm_l ! new non-zero components of v_k and v_l
    real(wp) :: commonFactorSO, commonFactorSS, gamma, gamma_diag, jiAlAklinvVl, &
                   jjAlAklinvVl, localEps

    ! V-quantities
    real(wp) :: VkAklinvVl, jiAkAklinvVl, jiAlAklinvVk, jiAklinvVk, jiAklinvVl, &
                   jjAkAklinvVl, jjAlAklinvVk, jjAklinvVk, jjAklinvVl

    ! W-quantities
    real(wp) :: WkAklinvWl, jiAklinvWk, &
                   jiAkAklinvWl, jiAlAklinvWk, jiAklinvWl, &
                   jjAkAklinvWl, jjAlAklinvWk, jjAklinvWk, jjAklinvWl

    !mixed quantities
    real(wp) :: VkAklinvWl, WkAklinvVl

    integer :: indexI, indexJ ! indices enumerating particles from H_SO and AMM operators

    localEps = 1.0e-14_wp ! if the corresponding spin mean value is less then localEps, we don't calculate the spatial part

    ! basically copy-paste from the old ExpecVals subroutine
    n=Glob_n
    np=Glob_np

    !First we build matrices Lk, Ll, Ak, Al from vechLk, vechLl.
    indx=0
    do i=1,n
      do j=i,n
        indx=indx+1
        Lk(i,j)=ZERO
        Lk(j,i)=vechLk(indx)
        Ll(i,j)=ZERO
        Ll(j,i)=vechLl(indx)
      enddo
    enddo

    do i=1,n
      do j=i,n
        temp1=ZERO
        do k=1,i
          temp1=temp1+Lk(i,k)*Lk(j,k)
        enddo
        tAk(i,j)=temp1
        tAk(j,i)=temp1
        temp1=ZERO
        do k=1,i
          temp1=temp1+Ll(i,k)*Ll(j,k)
        enddo
        tAl(i,j)=temp1
        tAl(j,i)=temp1
      enddo
    enddo

    !Then we permute elements of Al to account for
    !the action of the permutation matrix
    !  tAl=Pket'*Al*Pket

    !We also form matrix tAkl=tAk+tAl
    do i=1,n
      do j=1,n
        temp1=ZERO
        temp2=ZERO
        do k=1,n
          temp1=temp1+Pket(k,j)*tAl(k,i)
        enddo
        W1(j,i)=temp1
      enddo
    enddo
    !tAl=W1*Pket

    do i=1,n
      do j=i,n
        temp1=ZERO
        temp2=ZERO
        do k=1,n
          temp1=temp1+W1(j,k)*Pket(k,i)
        enddo
        tAl(j,i)=temp1
        tAl(i,j)=temp1
        tAkl(j,i)=temp1+tAk(j,i)
        tAkl(i,j)=temp1+tAk(i,j)
      enddo
    enddo

    !After this we can do Cholesky factorization of tAkl.
    !The Cholesky factor will be temporarily stored in the
    !lower triangle of W1
    det_tAkl=ONE
    do i=1,n
      do j=i,n
        temp1=tAkl(i,j)
        do k=i-1,1,-1
          temp1=temp1-W1(i,k)*W1(j,k)
        enddo
        if (i==j) then
          W1(i,i)=sqrt(temp1)
          det_tAkl=det_tAkl*temp1
        else
          W1(j,i)=temp1/W1(i,i)
          W1(i,j)=ZERO
        endif
      enddo
    enddo

    !Inverting tAkl using its Cholesky factor (stored in W1)
    !and placing the result into inv_tAkl
    do i=1,n
      W1(i,i)=ONE/W1(i,i)
      do j=i+1,n
        temp1=ZERO
        do k=i,j-1
          temp1=temp1-W1(j,k)*W1(k,i)
        enddo
        W1(j,i)=temp1/W1(j,j)
      enddo
    enddo

    do i=1,n
      do j=i,n
        temp1=ZERO
        do k=j,n
          temp1=temp1+W1(k,i)*W1(k,j)
        enddo
        inv_tAkl(i,j)=temp1
        inv_tAkl(j,i)=temp1
      enddo
    enddo

    ! new code from here

    ! new m_k and m_l
    ! new v_l = (P TRANSPOSED) * v_l
    pm_k = m_k
    pmm_k = mm_k

    pm_l = m_l
    pmm_l = mm_l
    do i = 1, n
      if (abs(Pket(m_l, i) - 1.0_wp) < 1.0e-13_wp) pm_l = i
      if (abs(Pket(mm_l, i) - 1.0_wp) < 1.0e-13_wp) pmm_l = i
    enddo


    ! 3P - > 1D: in commonfactorSO (1/sqrt(2)) - for consistent normalization with Skl (P-state) x 1/sqrt(2) from spin-tensor part)
    ! 3P - > 3D: in common factorSO (1/sqrt(2)) - for consistent normalization with Skl (P-state) x _(-(5/6)) from
    ! reduced matelem
    ! 3P - > 3D: in common factorSS (1/sqrt(2)) - for consistent normalization with Skl (P-state) x 1/sqrt(2) -
    ! (H_SS)_-1 prefactor x (1/sqrt(6)) - (H_SS)_0 prefactor x (3/sqrt(5/2)) - from reduced matelem
    commonFactorSS = 0
    if (Glob_selectTransition == 1) commonFactorSO = ONEHALF * Glob_PiRaised3n2 / (Glob_SqrtPi * det_tAkl * sqrt(det_tAkl))
    if (Glob_selectTransition == 2) then
      commonFactorSO = -sqrt(FIVE / (12._wp)) * Glob_PiRaised3n2 / (Glob_SqrtPi * det_tAkl * sqrt(det_tAkl))
      commonFactorSS = sqrt(15._wp) / TWO * Glob_PiRaised3n2 / (Glob_SqrtPi * det_tAkl * sqrt(det_tAkl))
    endif
    if (Glob_selectTransition == 3) then
      commonFactorSO=-ONEHALF * sqrt(FIVE / THREE) * Glob_PiRaised3n2 / (Glob_SqrtPi * det_tAkl * sqrt(det_tAkl))
    endif
    if (Glob_selectTransition == 4) then
      commonFactorSO = ONEHALF * Glob_PiRaised3n2 / (Glob_SqrtPi * det_tAkl * sqrt(det_tAkl))
      commonFactorSS = ONEHALF * Glob_PiRaised3n2 / (Glob_SqrtPi * det_tAkl * sqrt(det_tAkl))
    endif

    SSNCkl = ZERO

    SO1kl = ZERO
    SO2kl = ZERO

    AMM1kl = ZERO
    AMM2kl = ZERO

    AMM1Finkl = ZERO
    AMM2Finkl = ZERO

    WkAklinvWl = inv_tAkl(pmm_k, pmm_l)
    WkAklinvVl = inv_tAkl(pmm_k, pm_l)
    VkAklinvWl = inv_tAkl(pm_k, pmm_l)
    VkAklinvVl = inv_tAkl(pm_k, pm_l)

    do indexI = 1, n

      !if ((selectTransition == 1 .or. selectTransition == 3) .and. abs(SOspinME(indexI)) < localEps) cycle

      ! gamma diagonal coefficient
      gamma_diag = ONE / sqrt(inv_tAkl(indexI, indexI))

      ! calculating all the traces we need
      ! tr(Axy') is computed as (y, Ax) everywhere
      ! variable names: jiAlAklinvVk = (j^i, A_l A_{kl}^(-1) v_k) (names doesn't account for permutations)

      jiAkAklinvVl = ZERO
      do i = 1, n
        jiAkAklinvVl = jiAkAklinvVl + tAk(indexI, i) * inv_tAkl(i, pm_l)
      enddo
      jiAkAklinvWl = ZERO
      do i = 1, n
        jiAkAklinvWl = jiAkAklinvWl + tAk(indexI, i) * inv_tAkl(i, pmm_l)
      enddo

      jiAlAklinvVk = ZERO
      do i = 1, n
        jiAlAklinvVk = jiAlAklinvVk + tAl(indexI, i) * inv_tAkl(i, pm_k)
      enddo
      jiAlAklinvWk = ZERO
      do i = 1, n
        jiAlAklinvWk = jiAlAklinvWk + tAl(indexI, i) * inv_tAkl(i, pmm_k)
      enddo

      jiAklinvVl = inv_tAkl(indexI, pm_l)
      jiAklinvVk = inv_tAkl(indexI, pm_k)
      jiAklinvWl = inv_tAkl(indexI, pmm_l)
      jiAklinvWk = inv_tAkl(indexI, pmm_k)

      ! I term -> diagonal (spin-same orbit) matrix element (f[ii, ii])
      t1 = jiAklinvVk * jiAkAklinvVl + jiAlAklinvVk * jiAklinvVl
      temp1010 =  gamma_diag**3 / THREE * t1 *  WkAklinvWl - &
                 gamma_diag**5 / FIVE *  t1 * (jiAklinvWk * jiAklinvWl)

      ! Vl <-> Wl
      t1 = jiAklinvVk * jiAkAklinvWl + jiAlAklinvVk * jiAklinvWl
      temp1001 = gamma_diag**3 / THREE * t1 * WkAklinvVl - &
                 gamma_diag**5 / FIVE *t1 * (jiAklinvWk * jiAklinvVl)

      ! Vk <-> Wk
      t1 = jiAklinvWk * jiAkAklinvVl + jiAlAklinvWk * jiAklinvVl
      temp0110 = gamma_diag**3 / THREE * t1 * VkAklinvWl - &
                 gamma_diag**5 / FIVE * t1* (jiAklinvVk * jiAklinvWl)

      ! Vl <-> Wl, Vk <-> Wk
      t1 = jiAklinvWk * jiAkAklinvWl + jiAlAklinvWk * jiAklinvWl
      temp0101 =  gamma_diag**3 / THREE * t1 * VkAklinvVl - &
                 gamma_diag**5 / FIVE * t1 * (jiAklinvVk * jiAklinvVl)

      temp1 =  temp0101 + temp0110 -temp1010 - temp1001

      SO1kl = SO1kl + SOspinME(indexI) * SOmassChargeCoefficient(indexI, indexI, 1) * temp1
      AMM1kl = AMM1kl + SOspinME(indexI) * AMMmassChargeCoefficient(indexI, indexI, 1) * temp1
      AMM1Finkl = AMM1Finkl + SOspinME(indexI) * AMMFinmassChargeCoefficient(indexI, indexI, 1) * temp1

      ! these traces are needed for spin-other orbit contribution and SSNC
      do indexJ = 1, n
        if (indexI == indexJ) cycle

        gamma = ONE / sqrt(inv_tAkl(indexI, indexI) + inv_tAkl(indexJ, indexJ) - &
                           inv_tAkl(indexI, indexJ) - inv_tAkl(indexJ, indexI))

        jjAkAklinvVl = ZERO
        do i = 1, n
          jjAkAklinvVl = jjAkAklinvVl + tAk(indexJ, i) * inv_tAkl(i, pm_l)
        enddo
        jjAkAklinvWl = ZERO
        do i = 1, n
          jjAkAklinvWl = jjAkAklinvWl + tAk(indexJ, i) * inv_tAkl(i, pmm_l)
        enddo

        jjAlAklinvVk = ZERO
        do i = 1, n
          jjAlAklinvVk = jjAlAklinvVk + tAl(indexJ, i) * inv_tAkl(i, pm_k)
        enddo
        jjAlAklinvWk = ZERO
        do i = 1, n
          jjAlAklinvWk = jjAlAklinvWk + tAl(indexJ, i) * inv_tAkl(i, pmm_k)
        enddo

        jjAklinvVk = inv_tAkl(indexJ, pm_k)
        jjAklinvWk = inv_tAkl(indexJ, pmm_k)
        jjAklinvVl = inv_tAkl(indexJ, pm_l)
        jjAklinvWl = inv_tAkl(indexJ, pmm_l)

         !! II term -> f[ii, ij]
        t1 = jiAklinvVk * jjAkAklinvVl + jjAlAklinvVk * jiAklinvVl
        temp1010 = gamma_diag**3 / THREE * t1 * WkAklinvWl - &
                   gamma_diag**5 / FIVE * t1 * (jiAklinvWk * jiAklinvWl)

        !Vl <-> Wl
        t1 = jiAklinvVk * jjAkAklinvWl + jjAlAklinvVk * jiAklinvWl
        temp1001 = gamma_diag**3 / THREE * t1 * WkAklinvVl - &
                   gamma_diag**5 / FIVE * t1 * (jiAklinvWk * jiAklinvVl)

        !Vk <-> Wk
        t1 = jiAklinvWk * jjAkAklinvVl + jjAlAklinvWk * jiAklinvVl
        temp0110 = gamma_diag**3 / THREE * t1 * VkAklinvWl - &
                   gamma_diag**5 / FIVE * t1 * (jiAklinvVk * jiAklinvWl)

        !Vk <-> Wk, Vl <-> Wl
        t1 = jiAklinvWk * jjAkAklinvWl + jjAlAklinvWk * jiAklinvWl
        temp0101 = gamma_diag**3 / THREE * t1 * VkAklinvVl - &
                   gamma_diag**5 / FIVE * t1 * (jiAklinvVk * jiAklinvVl)

        temp1 =  temp0101 + temp0110 -temp1010 - temp1001

        SO2kl = SO2kl + SOspinME(indexI) * SOmassChargeCoefficient(indexI, indexI, 2) * temp1
        AMM2Finkl = AMM2Finkl + SOspinME(indexI) * AMMFinmassChargeCoefficient(indexI, indexI, 2) * temp1

         !! III term -> f[ij, jj]
        t1 = jjAkAklinvVl * (jjAklinvVk - jiAklinvVk) + jjAlAklinvVk * (jjAklinvVl - jiAklinvVl)
        t2 = jjAklinvWk * jjAklinvWl + jiAklinvWk * jiAklinvWl - jiAklinvWk * jjAklinvWl - jjAklinvWk * jiAklinvWl
        temp1010 = gamma**3 / THREE * t1 * WkAklinvWl - &
                   gamma**5 / FIVE * t1 * t2

        !Vl <-> Wl
        t1 = jjAkAklinvWl * (jjAklinvVk - jiAklinvVk) + jjAlAklinvVk * (jjAklinvWl - jiAklinvWl)
        t2 = jjAklinvWk * jjAklinvVl + jiAklinvWk * jiAklinvVl - jiAklinvWk * jjAklinvVl - jjAklinvWk * jiAklinvVl
        temp1001 = gamma**3 / THREE * t1 * WkAklinvVl - &
                   gamma**5 / FIVE * t1 * t2

        !Vk <-> Wk
        t1 = jjAkAklinvVl * (jjAklinvWk - jiAklinvWk) + jjAlAklinvWk * (jjAklinvVl - jiAklinvVl)
        t2 = jjAklinvVk * jjAklinvWl + jiAklinvVk * jiAklinvWl - jiAklinvVk * jjAklinvWl - jjAklinvVk * jiAklinvWl
        temp0110 = gamma**3 / THREE * t1 * VkAklinvWl - &
                   gamma**5 / FIVE * t1 * t2

         !!Vl <-> Wl, Vk <-> Wk
        t1 = jjAkAklinvWl * (jjAklinvWk - jiAklinvWk) + jjAlAklinvWk * (jjAklinvWl - jiAklinvWl)
        t2 = jjAklinvVk * jjAklinvVl + jiAklinvVk * jiAklinvVl - jiAklinvVk * jjAklinvVl - jjAklinvVk * jiAklinvVl
        temp0101 = gamma**3 / THREE * t1 * VkAklinvVl - &
                   gamma**5 / FIVE * t1 * t2

        temp1 =  temp0101 + temp0110 -temp1010 - temp1001

        SO2kl = SO2kl + SOspinME(indexI) * SOmassChargeCoefficient(indexI, indexJ, 3) * temp1
        AMM2kl = AMM2kl + SOspinME(indexI) * AMMmassChargeCoefficient(indexI, indexJ, 3) * temp1
        AMM2Finkl = AMM2Finkl + SOspinME(indexI) * AMMFinmassChargeCoefficient(indexI, indexJ, 3) * temp1

         !! IV term -> f[ij, ii] (i <->j of the III term)
        t1 = jiAkAklinvVl * (jiAklinvVk - jjAklinvVk) + jiAlAklinvVk * (jiAklinvVl - jjAklinvVl)
        t2 = jjAklinvWk * jjAklinvWl + jiAklinvWk * jiAklinvWl - jiAklinvWk * jjAklinvWl - jjAklinvWk * jiAklinvWl
        temp1010 = gamma**3 / THREE * t1 * WkAklinvWl - &
                   gamma**5 / FIVE * t1 * t2

        !Vl <-> Wl
        t1 = jiAkAklinvWl * (jiAklinvVk - jjAklinvVk) + jiAlAklinvVk * (jiAklinvWl - jjAklinvWl)
        t2 = jjAklinvWk * jjAklinvVl + jiAklinvWk * jiAklinvVl - jiAklinvWk * jjAklinvVl - jjAklinvWk * jiAklinvVl
        temp1001 = gamma**3 / THREE * t1 * WkAklinvVl - &
                   gamma**5 / FIVE * t1 * t2

        !Vk <-> Wk
        t1 = jiAkAklinvVl * (jiAklinvWk - jjAklinvWk) + jiAlAklinvWk * (jiAklinvVl - jjAklinvVl)
        t2 = jjAklinvVk * jjAklinvWl + jiAklinvVk * jiAklinvWl - jiAklinvVk * jjAklinvWl - jjAklinvVk * jiAklinvWl
        temp0110 = gamma**3 / THREE * t1 * VkAklinvWl - &
                   gamma**5 / FIVE * t1 * t2

         !!Vl <-> Wl, Vk <-> Wk
        t1 = jiAkAklinvWl * (jiAklinvWk - jjAklinvWk) + jiAlAklinvWk * (jiAklinvWl - jjAklinvWl)
        t2 = jjAklinvVk * jjAklinvVl + jiAklinvVk * jiAklinvVl - jiAklinvVk * jjAklinvVl - jjAklinvVk * jiAklinvVl
        temp0101 = gamma**3 / THREE * t1 * VkAklinvVl - &
                   gamma**5 / FIVE * t1 * t2

        temp1 =  temp0101 + temp0110 - temp1010 - temp1001

        SO2kl = SO2kl + SOspinME(indexI) * SOmassChargeCoefficient(indexI, indexJ, 4) * temp1
        AMM2kl = AMM2kl + SOspinME(indexI) * AMMmassChargeCoefficient(indexI, indexJ, 4) * temp1
        AMM2Finkl = AMM2Finkl + SOspinME(indexI) * AMMFinmassChargeCoefficient(indexI, indexJ, 4) * temp1

         !!SSNC term
        !only 3P -> 3D transition contributes, indexI > indexJ
        if (Glob_selectTransition == 1 .or. Glob_selectTransition == 3 .or. &
            indexJ <= indexI .or. abs(SSNCspinME(indexI, indexJ)) < localEps ) cycle

        !SSNC term
        !g1010
        t1 = jiAklinvVk * jiAklinvVl + jjAklinvVk * jjAklinvVl - jiAklinvVk * jjAklinvVl - jjAklinvVk * jiAklinvVl
        t2 = jiAklinvWk * jiAklinvWl + jjAklinvWk * jjAklinvWl - jiAklinvWk * jjAklinvWl - jjAklinvWk * jiAklinvWl
        temp1010 = gamma**5 / FIVE * (t1 * WkAklinvWl) - &
                   gamma**7 / SEVEN * (t1 * t2)

        !Vl <-> Wl
        t1 = jiAklinvVk * jiAklinvWl + jjAklinvVk * jjAklinvWl - jiAklinvVk * jjAklinvWl - jjAklinvVk * jiAklinvWl
        t2 = jiAklinvWk * jiAklinvVl + jjAklinvWk * jjAklinvVl - jiAklinvWk * jjAklinvVl - jjAklinvWk * jiAklinvVl
        temp1001 = gamma**5 / FIVE * (t1 * WkAklinvVl) - &
                   gamma**7 / SEVEN * (t1 * t2)

        !Vk <-> Wk
        t1 = jiAklinvWk * jiAklinvVl + jjAklinvWk * jjAklinvVl - jiAklinvWk * jjAklinvVl - jjAklinvWk * jiAklinvVl
        t2 = jiAklinvVk * jiAklinvWl + jjAklinvVk * jjAklinvWl - jiAklinvVk * jjAklinvWl - jjAklinvVk * jiAklinvWl
        temp0110 = gamma**5 / FIVE * (t1 * VkAklinvWl) - &
                   gamma**7 / SEVEN * (t1 * t2)

        !Vl <-> Wl, Vk <-> Wk
        t1 = jiAklinvVk * jiAklinvVl + jjAklinvVk * jjAklinvVl - jiAklinvVk * jjAklinvVl - jjAklinvVk * jiAklinvVl
        t2 = jiAklinvWk * jiAklinvWl + jjAklinvWk * jjAklinvWl - jiAklinvWk * jjAklinvWl - jjAklinvWk * jiAklinvWl
        temp0101 = gamma**5 / FIVE * (t2 * VkAklinvVl) - &
                   gamma**7 / SEVEN * (t1 * t2)

        temp1 = -temp1010 - temp1001 + temp0110 + temp0101
        SSNCkl = SSNCkl + SSNCspinME(indexI, indexJ) * SSNCmassChargeCoefficient(indexI, indexJ) * temp1

      enddo ! indexJ cycle
    enddo ! indexI cycle

    SSNCkl = SSNCkl * commonFactorSS
    SO1kl = SO1kl * commonFactorSO
    SO2kl = SO2kl * commonFactorSO
    AMM1kl = AMM1kl * commonFactorSO
    AMM2kl = AMM2kl * commonFactorSO
    AMM1Finkl = AMM1Finkl * commonFactorSO
    AMM2Finkl = AMM2Finkl * commonFactorSO

  end subroutine spinDependentMatrixElements

end module matelem

