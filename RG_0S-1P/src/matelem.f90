module matelem
!Module matelem contains subroutines for computing
!matrix elements with real L=0 and L=1 Gaussians.
  use globvars
  implicit none

contains

  subroutine OverLapElementS0(vechLk, P, Skk)

!
! This subroutine computes symmetry adapted matrix element <P*fk|P*fk>
! with two real L=0 correlated Gaussians

! fk =  exp[-r'(Lk*Lk')r]
!
! Symmetry adaption is applied to the bra using
! permutation matrix Pk and to the ket using permutation matrix Pl.
!
! Input :
!                vechLl                ::        Arrays of length (n(n+1)/2) of exponential parameters.
!                P                  ::        The symmetry permutation matrix of size n x n
!
! Output:
!                Skk                        ::        < P*fk | P*fk >
!
!                      <P*fk|Y*fk>               abs(det_Lk)^3
!   < P*fk | P*fk > = --------------= 2^(3*n/2) ----------------
!                        <fk|fk>                  det_tAkk^3/2

!Arguments
    real(wp),intent(in)   :: vechLk(Glob_np)
    real(wp),intent(in)   :: P(Glob_n,Glob_n)
    real(wp),intent(out)  :: Skk

!Parameters (These are needed to declare static arrays. Using static
!arrays makes the function call a little faster in comparison with
!the case when arrays are dynamically allocated in stack)
    integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles

!Local variables
    real(wp),allocatable,dimension(:)         :: vechLl

    integer           n, np
    real(wp)       Lk(nn,nn), Ll(nn,nn), inv_Lk(nn,nn), inv_Ll(nn,nn)
    real(wp)       tAk(nn,nn), tAl(nn,nn), tAkl(nn,nn), Ak(nn,nn)
    real(wp)       inv_tAkl(nn,nn)
    real(wp)       W1(nn,nn), W2(nn,nn)
    real(wp)       temp1, temp2
    real(wp)       det_Lk, det_Ll, det_tAkl
    integer           i,j,k,indx,q

    n=Glob_n
    np=Glob_np

    allocate(vechLl(n))
    vechLl=vechLk

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

!The determinants of Lk and Ll are just
!the products of their diagonal elements
    det_Lk=ONE
    det_Ll=ONE
    do i=1,n
      det_Lk=det_Lk*Lk(i,i)
      det_Ll=det_Ll*Ll(i,i)
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

!Evaluating overlap
    temp2=abs(det_Ll*det_Lk)
    temp1=temp2/det_tAkl
    Skk=Glob_2Raised3n2*temp1*sqrt(temp1)

  end subroutine OverLapElementS0

  SUBROUTINE NormalizedOverlapMatElem_RG_1P(m_k, vechLk, P, Skk)
!==============================================================================================
! This subroutine calculates the diagonal overlap for L=1 correlated Gaussians:
!
! fk = z{m_k} * exp[-r'(Lk*Lk')r]
!
!                  <P*fk|P*fk>               ||L_kk||^3       v_k' * inv_tAkk * tv_k
! <psi_k|psi_k> = -------------- = 2^(3n/2) ------------- * --------------------------
!                    <fk|fk>                 |tAkk|^{3/2}     v_k' * inv_Akk * v_k
!
! where tAkk = Ak + P'*Ak*P,  tv_k = P'*v_k
!
!----------------------------------------------------------------------------------------------
! Input:
!   m_k      :: integer that determines which z-component is in the premultiplier
!   vechLk   :: Array of length n(n+1)/2 of exponential parameters
!   P        :: The symmetry permutation matrix of size n x n
!
! Output:
!   Skk      :: <psi_k|psi_k>
!==============================================================================================

    IMPLICIT NONE

    !---------------------------------------------------------------------------
    ! ARGUMENTS
    !---------------------------------------------------------------------------
    INTEGER, INTENT(IN)      :: m_k
    REAL(wp), INTENT(IN)  :: vechLk(Glob_np)
    REAL(wp), INTENT(IN)  :: P(Glob_n, Glob_n)
    REAL(wp), INTENT(OUT) :: Skk

    !---------------------------------------------------------------------------
    ! LOCAL VARIABLES
    !---------------------------------------------------------------------------
    INTEGER                  :: n
    INTEGER                  :: i, j, k, indx

    REAL(wp)              :: Lk(Glob_n, Glob_n)
    REAL(wp)              :: tAk(Glob_n, Glob_n), tAkk(Glob_n, Glob_n)
    REAL(wp)              :: inv_tAkk(Glob_n, Glob_n), inv_Akk(Glob_n, Glob_n)
    REAL(wp)              :: W1(Glob_n, Glob_n)

    REAL(wp)              :: tvk(Glob_n)

    REAL(wp)              :: temp1
    REAL(wp)              :: det_Lk, det_tAkk, tau3

    !===========================================================================
    n = Glob_n

    !---------------------------------------------------------------------------
    ! 1. Build lower-triangular Lk from vechLk
    !---------------------------------------------------------------------------
    indx = 0
    DO i = 1, n
      DO j = i, n
        indx = indx + 1
        Lk(i, j) = ZERO
        Lk(j, i) = vechLk(indx)
      END DO
    END DO

    !---------------------------------------------------------------------------
    ! 2. Compute inv_Akk using Cholesky factor Lk
    !    inv_Akk = (1/2) * Lk^{-T} * Lk^{-1} = (2*Ak)^{-1} = Akk^{-1}
    !---------------------------------------------------------------------------
    W1(:,:) = ZERO

    DO i = 1, n
      W1(i,i) = ONE / Lk(i,i)
      DO j = i + 1, n
        temp1 = ZERO
        DO k = i, j - 1
          temp1 = temp1 - Lk(j,k) * W1(k,i)
        END DO
        W1(j,i) = temp1 / Lk(j,j)
      END DO
    END DO

    DO i = 1, n
      DO j = i, n
        temp1 = ZERO
        DO k = j, n
          temp1 = temp1 + W1(k,i) * W1(k,j)
        END DO
        inv_Akk(i,j) = ONEHALF * temp1
        inv_Akk(j,i) = ONEHALF * temp1
      END DO
    END DO

    !---------------------------------------------------------------------------
    ! 3. Compute Ak = Lk * Lk'
    !---------------------------------------------------------------------------
    DO i = 1, n
      DO j = i, n
        temp1 = ZERO
        DO k = 1, i
          temp1 = temp1 + Lk(i, k) * Lk(j,k)
        END DO
        tAk(i,j) = temp1
        tAk(j,i) = temp1
      END DO
    END DO

    !---------------------------------------------------------------------------
    ! 4. Apply symmetry permutation and form tAkk = Ak + P' * Ak * P
    !---------------------------------------------------------------------------
    W1(:,:) = ZERO

    ! W1 = P' * Ak
    DO i = 1, n
      DO j = 1, n
        temp1 = ZERO
        DO k = 1, n
          temp1 = temp1 + P(k,j) * tAk(k,i)
        END DO
        W1(j,i) = temp1
      END DO
    END DO

    ! tAkk = Ak + P' * Ak * P
    DO i = 1, n
      DO j = i, n
        temp1 = ZERO
        DO k = 1, n
          temp1 = temp1 + W1(i, k) * P(k,j)
        END DO
        tAkk(i,j) = tAk(i,j) + temp1
        tAkk(j,i) = tAkk(i,j)
      END DO
    END DO

    !---------------------------------------------------------------------------
    ! 5. Determinant of Lk (product of diagonal elements)
    !---------------------------------------------------------------------------
    det_Lk = ONE
    DO i = 1, n
      det_Lk = det_Lk * Lk(i,i)
    END DO

    !---------------------------------------------------------------------------
    ! 6. Cholesky factorization of tAkk -> W1, and determinant
    !---------------------------------------------------------------------------
    det_tAkk = ONE

    DO i = 1, n
      DO j = i, n
        temp1 = tAkk(i,j)
        DO k = i - 1, 1, -1
          temp1 = temp1 - W1(i, k) * W1(j,k)
        END DO
        IF (i == j) THEN
          W1(i,i) = SQRT(temp1)

          IF (temp1 <= 0.0_wp) THEN
            WRITE(*,*) 'ERROR: tAkk not SPD at i=',i
            call MPI_Abort(MPI_COMM_WORLD, 1, Glob_MPIErrCode) !error stop
          END IF

          det_tAkk = det_tAkk * temp1
        ELSE
          W1(j,i) = temp1 / W1(i,i)
          W1(i,j) = ZERO
        END IF
      END DO
    END DO

    !---------------------------------------------------------------------------
    ! 7. Invert tAkk using its Cholesky factor (stored in W1)
    !---------------------------------------------------------------------------
    DO i = 1, n
      W1(i,i) = ONE / W1(i,i)
      DO j = i + 1, n
        temp1 = ZERO
        DO k = i, j - 1
          temp1 = temp1 - W1(j,k) * W1(k,i)
        END DO
        W1(j,i) = temp1 / W1(j,j)
      END DO
    END DO

    DO i = 1, n
      DO j = i, n
        temp1 = ZERO
        DO k = j, n
          temp1 = temp1 + W1(k,i) * W1(k,j)
        END DO
        inv_tAkk(i,j) = temp1
        inv_tAkk(j,i) = temp1
      END DO
    END DO

    !---------------------------------------------------------------------------
    ! 8. Compute tau3 = v_k' * inv_tAkk * tv_k
    !    where v_k = e_{m_k} (unit vector) and tv_k = P(m_k, :)
    !---------------------------------------------------------------------------
    DO i = 1, n
      tvk(i) = P(m_k, i)
    END DO

    tau3 = ZERO
    DO i = 1, n
      tau3 = tau3 + inv_tAkk(m_k, i) * tvk(i)
    END DO

    !---------------------------------------------------------------------------
    ! 9. Final overlap calculation
    !    Skk = 2^{3n/2} * ||L_kk||^3 / |tAkk|^{3/2} * tau3 / inv_Akk(m_k, m_k)
    !---------------------------------------------------------------------------
    temp1 = (det_Lk * det_Lk) / det_tAkk
    temp1 = ABS(temp1) * SQRT(ABS(temp1))

    Skk = Glob_2Raised3n2 * temp1 * tau3 / inv_Akk(m_k, m_k)

  END SUBROUTINE NormalizedOverlapMatElem_RG_1P

  subroutine MatElemTranDipoleMoment_RG_0S_1P(ml, vechLk, vechLl, Pk, Pl, TranDipolLength_kl,TranDipolVelocity_kl)

!This subroutine computes symmetry adapted matrix element with
!a real L=0 and a real L=1 correlated Gaussians:
!
!fk = exp[-r'(Lk*Lk')r]
!fl = z_{ml} exp[-r'(Ll*Ll')r]
!
!where m_l is some integer between 1 and n (n is the number of
!pseudoparticles). Symmetry adaption is applied to the bra using
!permutation matrix Pk and to the ket using permutation matrix Pl.

!Pk fk = exp[-r' {Pk'*(Lk*Lk')*Pk} r]
!Pl fl = (Pl z_{m_l}) exp[-r' {Pl'*(Ll*Ll')*Pl} r]

! this subroutine computes the following expressions
! which is a part of transition dipole momentum calculation.
!
!
!
!====================================================================
!                                Transition dipole integral in Length form
!====================================================================
!
!
!                                            m_i         < Pk*fk0 | z_i | Pl*fl1 >
! TranDipolLength_kl = SUM{i}(q_i - Q_tot * -----)*-------------------------------------------
!                                            m0      Sqrt( <fk0|fk0> ) * Sqrt( <fl1|fl1> )
!
!
!      < Pk*fk0 | z_i | Pl*fl1 >               2^(3*n/2)     (abs(det_Lk))^1.5 * (abs(det_Ll))^1.5        vi'*inv_tAkl*vl
! ----------------------------------------- = ----------- * ---------------------------------------- * -----------------------
!   Sqrt( <fk0|fk0> ) * Sqrt( <fl1|fl1> )       sqrt(2)                (det_tAkl)^1.5                   Sqrt(vl'*inv_All*vl)
!
!
!
!
!
!====================================================================
!                                Transition dipole integral in Length form
!====================================================================
!
!
!                                    < Pk*fk0 | P(z_i) | Pl*fl1 >
! TranDipolLength_kl = SUM{i}*-------------------------------------------
!                                   Sqrt( <fk0|fk0> ) * Sqrt( <fl1|fl1> )
!
!
!     < Pk*fk0 | P(z_i) | Pl*fl1 >             2^(3*n/2)     (abs(det_Lk))^1.5 * (abs(det_Ll))^1.5       vi'*tAk*inv_tAkl*vl
! ----------------------------------------- = ----------- * ---------------------------------------- * -----------------------
!   Sqrt( <fk0|fk0> ) * Sqrt( <fl1|fl1> )       sqrt(2)                (det_tAkl)^1.5                    Sqrt(vl'*inv_All*vl)
!
!

!Input:
!        m_l                                ::        integer that determine which z-component is in the
!                                                premultiplier of the Gaussian
!        vechLk, vechLl        ::        Arrays of length (n(n+1)/2) of exponential parameters.
!        Pk, Pl                        ::        The symmetry permutation matrices of size n x n
!
!Output:
!        TDkl                          ::        Matrix element (normalized)

!Arguments
    integer, intent(in)          :: ml
    real(wp), intent(in)      :: vechLk(Glob_np), vechLl(Glob_np)
    real(wp), intent(in)      :: Pk(Glob_n,Glob_n), Pl(Glob_n,Glob_n)
    real(wp), intent(out)     :: TranDipolLength_kl, TranDipolVelocity_kl

!Parameters (These are needed to declare static arrays. Using static
!arrays makes the function call a little faster in comparison with
!the case when arrays are dynamically allocated in stack)
    integer,parameter            :: nn=Glob_AllowedNumOfPseudoParticles
    integer,parameter            :: nnp=nn*(nn+1)/2

!Local variables
    integer           :: n, np, Qtotal
    real(wp)       :: Lk(nn,nn),Ll(nn,nn)
    real(wp)       :: inv_Lk(nn,nn),inv_Ll(nn,nn)
    real(wp)       :: det_Lk,det_Ll
    real(wp)       :: tAk(nn,nn),tAl(nn,nn),tAkl(nn,nn)
    real(wp)       :: inv_Akk(nn,nn),inv_All(nn,nn),inv_tAkl(nn,nn),tAk_inv_tAkl(nn,nn)
    real(wp)       :: det_tAkl
    real(wp)       :: tvl(nn),vi(nn),vi_tAk(nn),vi_tAk_inv_tAkl(nn)

    integer           :: i,j,k,indx,ii
    real(wp)       :: temp0, temp1, temp2, temp3, charge_mass0
    real(wp)       :: W1(nn,nn),W2(nn,nn)

    n=Glob_n
    np=Glob_np

!First we build matrices Lk, Ll, Ak, Al from vechLk, vechLl.
    indx=0
    Do i=1,n
      Do j=i,n
        indx=indx+1
        Lk(i,j)=ZERO
        Lk(j,i)=vechLk(indx)
        Ll(i,j)=ZERO
        Ll(j,i)=vechLl(indx)
      EndDo
    EndDo

    Do i=1,n
      Do j=i,n
        temp1=ZERO
        Do k=1,i
          temp1=temp1+Lk(i,k)*Lk(j,k)
        EndDo
        tAk(i,j)=temp1
        tAk(j,i)=temp1
        temp1=ZERO
        Do k=1,i
          temp1=temp1+Ll(i,k)*Ll(j,k)
        EndDo
        tAl(i,j)=temp1
        tAl(j,i)=temp1
      EndDo
    EndDo

!Then we permute elements of Ak and Al to account for
!the action of the permutation matrix
!  tAl=Pl'*Al*Pl
!  tAk=Pk'*Ak*Pk
!We also form matrix tAkl=tAk+tAl
    Do i=1,n
      Do j=1,n
        temp1=ZERO
        temp2=ZERO
        Do k=1,n
          temp1=temp1+Pl(k,j)*tAl(k,i)
          temp2=temp2+tAk(j,k)*Pk(k,i)
        EndDo
        W1(j,i)=temp1
        W2(j,i)=temp2
      EndDo
    EndDo
!tAl=W1*Pl
!tAk=Pk'*W2
    Do i=1,n
      Do j=i,n
        temp1=ZERO
        temp2=ZERO
        Do k=1,n
          temp1=temp1+W1(j,k)*Pl(k,i)
          temp2=temp2+Pk(k,j)*W2(k,i)
        EndDo
        tAl(j,i)=temp1
        tAl(i,j)=temp1
        tAk(j,i)=temp2
        tAk(i,j)=temp2
        tAkl(j,i)=temp1+temp2
        tAkl(i,j)=temp1+temp2
      EndDo
    EndDo

!The determinants of Lk and Ll are just
!the products of their diagonal elements
    det_Lk=ONE
    det_Ll=ONE
    Do i=1,n
      det_Lk=det_Lk*Lk(i,i)
      det_Ll=det_Ll*Ll(i,i)
    EndDo

!After this we can do Cholesky factorization of tAkl.
!The Cholesky factor will be temporarily stored in the
!lower triangle of W1
    det_tAkl=ONE
    Do i=1,n
      Do j=i,n
        temp1=tAkl(i,j)
        Do k=i-1,1,-1
          temp1=temp1-W1(i,k)*W1(j,k)
        EndDo
        IF (i==j) then
          W1(i,i)=sqrt(temp1)
          det_tAkl=det_tAkl*temp1
        Else
          W1(j,i)=temp1/W1(i,i)
          W1(i,j)=ZERO
        EndIF
      EndDo
    EndDo

!Inverting tAkl using its Cholesky factor (stored in W1)
!and placing the result into inv_tAkl
    Do i=1,n
      W1(i,i)=ONE/W1(i,i)
      Do j=i+1,n
        temp1=ZERO
        Do k=i,j-1
          temp1=temp1-W1(j,k)*W1(k,i)
        EndDo
        W1(j,i)=temp1/W1(j,j)
      EndDo
    EndDo

    Do i=1,n
      Do j=i,n
        temp1=ZERO
        Do k=j,n
          temp1=temp1+W1(k,i)*W1(k,j)
        EndDo
        inv_tAkl(i,j)=temp1
        inv_tAkl(j,i)=temp1
      EndDo
    EndDo

!Finding the inverse of Akk and All using their Cholesky factors
!The result is placed in inv_Akk and inv_All
    Do i=1,n
      ! W1(i,i)=ONE/Lk(i,i)
      W2(i,i)=ONE/Ll(i,i)
      Do j=i+1,n
        ! temp1=ZERO
        temp2=ZERO
        Do k=i,j-1
          ! temp1=temp1-Lk(j,k)*W1(k,i)
          temp2=temp2-Ll(j,k)*W2(k,i)
        EndDo
        ! W1(j,i)=temp1/Lk(j,j)
        W2(j,i)=temp2/Ll(j,j)
      EndDo
    EndDo

    Do i=1,n
      Do j=i,n
        ! temp1=ZERO
        temp2=ZERO
        Do k=j,n
          ! temp1=temp1+W1(k,i)*W1(k,j)
          temp2=temp2+W2(k,i)*W2(k,j)
        EndDo
        ! inv_Akk(i,j)=ONEHALF*temp1
        ! inv_Akk(j,i)=ONEHALF*temp1
        inv_All(i,j)=ONEHALF*temp2
        inv_All(j,i)=ONEHALF*temp2
      EndDo
    EndDo

!tvl=Pl'*vl
    Do i=1,n
      tvl(i)=Pl(ml,i)
    EndDo

!====================================================
!                        Evaluating Matrix Elements
!====================================================

!         (abs(det_Lk))^1.5 * (abs(det_Ll))^1.5
! temp1= ----------------------------------------
!                    (det_tAkl)^1.5

    temp1 = abs(det_Ll*det_Lk) / det_tAkl
    temp1 = temp1 * sqrt(abs(temp1))

!                       2^(3*n/2)            temp1
! TranDipolLength_kl = ----------- * -----------------------
!                        sqrt(2)      Sqrt(vl'*inv_All*vl)

    TranDipolLength_kl = Glob_2Raised3n2 * temp1 / sqrt(TWO*inv_All(ml,ml))
    TranDipolVelocity_kl = TranDipolLength_kl * TWO

!
!                                                                m_i
! TranDipolLength_kl = TranDipolLength_kl *SUM{i}(q_i - Q_tot * -----)* vi'*inv_tAkl*vl
!                                                                m0
!
    Qtotal = Glob_PseudoCharge0

    Do i=1,n
      Qtotal = Qtotal + Glob_PseudoCharge(i)
    EndDo

    temp1 = ZERO

    Do i=1,n                                                         !pseudo-particles

      temp2 = ZERO
      Do j=1,n                                                 !trace elements
        temp2 = temp2 + inv_tAkl(i,j)* tvl(j)
      EndDo

      temp0 = Glob_PseudoCharge(i) - Qtotal*Glob_Mass(i+1)/Glob_MassTotal
      temp1 = temp1 + temp0 * temp2

    EndDo

    TranDipolLength_kl = TranDipolLength_kl * temp1

!
!
! TranDipolLength_kl = TranDipolLength_kl * SUM{i} vi' * tAk inv_tAkl * vl
!
!

    temp0 = ZERO
    temp1 = ZERO
    temp2 = ZERO
    temp3 = ZERO

    charge_mass0 = Glob_PseudoCharge0 / Glob_Mass(1)

    vi_tAk = ZERO
    vi_tAk_inv_tAkl = ZERO

    Do i=1,n

      vi = ZERO
      vi(i) = ONE

      temp0 = charge_mass0 - Glob_PseudoCharge(i)/Glob_Mass(i+1)

      Do k=1,n

        temp1 = ZERO
        Do j=1,n
          temp1 = temp1 + vi(j) * tAk(j,k)
        EndDo
        vi_tAk(k) = temp1

      EndDo

      Do k=1,n

        temp2 = ZERO
        Do j=1,n
          temp2 = temp2 + vi_tAk(j) * inv_tAkl(j,k)
        EndDo
        vi_tAk_inv_tAkl(k) = temp2

      EndDo

      Do k=1,n

        temp3 = temp3 + vi_tAk_inv_tAkl(k) * tvl(k) * temp0

      EndDo

    EndDo

    TranDipolVelocity_kl = TranDipolVelocity_kl * temp3

  end subroutine MatElemTranDipoleMoment_RG_0S_1P

end module matelem

