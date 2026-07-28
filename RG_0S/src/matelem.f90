module matelem
!Module matelem contains subroutines for computing
!matrix elements with real L=0 Gaussians without normalization.
  use globvars
  implicit none

contains

  subroutine MatrixElementsHS_RG_0S(vechLk, vechLl, P, &
                            Hkl, Skl, Dk, Dl, grad_k, grad_l)
!This subroutine computes symmetry adapted matrix element with
!two real L=0 correlated Gaussians:
!
!fk =  exp[-r'(Lk*Lk')r]
!
!Symmetry adaption is applied to the ket using
!permutation matrices Glob_YHYMatr(:,:,1:Glob_NumYHYTerms)
!
!Input:
!   vechLk, vechLl :: Arrays of length (n(n+1)/2) of exponential parameters.
!   P   :: The symmetry permutation matrix of size n x n
!   grad_k, grad_l :: Gradient flags
!   grad_k=.true.  means that dHkldvechLk, dSkldvechLk need to be computed.
!   grad_l=.true.  means that dHkldvechLl, dSkldvechLl need to be computed.
!Output:
!   Hkl         ::        Hamiltonian term (normalized)
!   Skl         ::        Overlap matrix element (normalized)
!   Dk,Dl:: derivatives of normalized Hkl and Skl wrt vechLk
!           and vechLl respectively. They are ordered in the
!           following manner:
!           Dk=(dHkldvechLk,dSkldvechLk)
!           Dl=(dHkldvechLl,dSkldvechLl)
!
!This is an optimized rewrite of the original routine (the computed
!values are identical to the original ones up to roundoff). The part
!up to Hkl is essentially the original code; the whole gradient part
!is new. The main algorithmic change is the elimination of the
!original O(n^5) loop over particle pairs in the gradient of the
!potential energy: each pair term tQ=inv_tAkl*Jij*inv_tAkl is a
!congruence of the sparse matrix Jij, so the weighted sum of all of
!them collapses analytically into
!   B = inv_tAkl*C*inv_tAkl,  C = sum_ij c_ij*Jij,
!where C is assembled in O(n^2), which makes the whole gradient part
!O(n^3). Details are explained in the comments in the body.

!Arguments
    real(wp),intent(in)      :: vechLk(Glob_np), vechLl(Glob_np)
    real(wp),intent(in)      :: P(Glob_n,Glob_n)
    real(wp),intent(out)     :: Skl,Hkl
    real(wp),intent(out)     :: Dk(2*Glob_np),Dl(2*Glob_np)
    logical,intent(in)          :: grad_k, grad_l

    integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles

!Local variables
    integer           n, np
    real(wp)       Lk(nn,nn), Ll(nn,nn), PT(nn,nn)
    real(wp)       Ak(nn,nn), tAl(nn,nn), tAkl(nn,nn)
    real(wp)       inv_tAkl(nn,nn), inv_ttAkl(nn,nn)
    real(wp)       inv_tAkltAlM(nn,nn)
    real(wp)       tr_inv_tAklJij32(nn,nn)
    real(wp)       F(nn,nn), G(nn,nn)
    real(wp)       Cmat(nn,nn), Bmat(nn,nn), Z(nn,nn)
    real(wp)       W1(nn,nn), W2(nn,nn)
    real(wp)       temp1, temp2, temp3, temp4, temp5
    real(wp)       det_tAkl
    real(wp)       Tkl, Vkl, cV, HklOverSkl
    integer           i, j, k, indx

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
!det_Lk=ONE
!det_Ll=ONE
!do i=1,n
!  det_Lk=det_Lk*Lk(i,i)
!  det_Ll=det_Ll*Ll(i,i)
!enddo

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

!temp1=abs(det_Ll*det_Lk)/det_tAkl
!Skl=Glob_2Raised3n2*temp1*sqrt(temp1)
    Skl=Glob_PiRaised3n2/(det_tAkl*sqrt(det_tAkl))  !new line

!Doing multiplication W2=inv_tAkl*tAl
    do i=1,n
      do j=1,n
        temp1=ZERO
        do k=1,n
          temp1=temp1+inv_tAkl(j,k)*tAl(k,i)
        enddo
        W2(j,i)=temp1
      enddo
    enddo

!Doing multiplication inv_tAkltAlM=inv_tAkl*tAl*M=W2*M
    do i=1,n
      do j=1,n
        temp1=ZERO
        do k=1,n
          temp1=temp1+W2(j,k)*Glob_MassMatrix(k,i)
        enddo
        inv_tAkltAlM(j,i)=temp1
      enddo
    enddo

!Computing kinetic energy, Tkl=tr[inv_tAkltAlM*Ak]
    Tkl=ZERO
    do i=1,n
      temp1=ZERO
      do k=1,n
        temp1=temp1+inv_tAkltAlM(i,k)*Ak(k,i)
      enddo
      Tkl=Tkl+temp1
    enddo
    Tkl=SIX*Skl*Tkl

!Evaluating potential energy, Vkl, and tr[invCkl*Jij]^(-3/2)
!The lower triangle of array trinvCklJij32
!will contain the corresponding quantities. The latter are needed
!only for the gradients, so in the gradientless case a leaner loop
!(one division per particle pair less) is used.
    temp1=(TWO/Glob_SqrtPi)*Skl
    Vkl=ZERO
    if (grad_k.or.grad_l) then
      do i=1,n
        temp3=inv_tAkl(i,i)
        temp4=sqrt(temp3)
        tr_inv_tAklJij32(i,i)=1/(temp4*temp3)
        temp5=temp1/temp4
        Vkl=Vkl+Glob_ScaledPseudoChargeMatrix(i,0)*temp5
      enddo
      do i=1,n
        do j=i+1,n
          temp3=inv_tAkl(i,i)+inv_tAkl(j,j)-inv_tAkl(j,i)-inv_tAkl(j,i)
          temp4=sqrt(temp3)
          tr_inv_tAklJij32(j,i)=1/(temp4*temp3)
          temp5=temp1/temp4
          Vkl=Vkl+Glob_ScaledPseudoChargeMatrix(i,j)*temp5
        enddo
      enddo
    else
      do i=1,n
        temp5=temp1/sqrt(inv_tAkl(i,i))
        Vkl=Vkl+Glob_ScaledPseudoChargeMatrix(i,0)*temp5
      enddo
      do i=1,n
        do j=i+1,n
          temp3=inv_tAkl(i,i)+inv_tAkl(j,j)-inv_tAkl(j,i)-inv_tAkl(j,i)
          temp5=temp1/sqrt(temp3)
          Vkl=Vkl+Glob_ScaledPseudoChargeMatrix(i,j)*temp5
        enddo
      enddo
    endif

    Hkl=Tkl+Vkl

!In the gradientless case we are done
    if (.not.(grad_k.or.grad_l)) return

!Now we compute the gradients. In what follows
!  W2 = inv_tAkl*tAl  and  K = inv_tAkltAlM = W2*M
!(both are already available from the computation of Tkl).
!
!Gradient of Skl (lower triangles, in vech order):
!  dSkl/dvechLk = -3*Skl*inv_tAkl*Lk
!  dSkl/dvechLl = -3*Skl*(P*inv_tAkl*P')*Ll
!They are stored directly in the tails Dk(np+1:2*np), Dl(np+1:2*np)
!of the output arrays.
!
!Gradient of Tkl. It needs the two matrices
!  F = inv_tAkl*tAl*M*tAl*inv_tAkl = K*W2'
!  G = P*(inv_tAkl*Ak*M*Ak*inv_tAkl)*P'
!Because inv_tAkl*Ak = I-W2 (as Ak+tAl=tAkl), the matrix in the
!parentheses of G equals (I-W2)*M*(I-W2') = M - K - K' + F, i.e. it
!costs only O(n^2) once F is known.
!
!Gradient of Vkl. The original code summed, over all particle pairs,
!the matrices tQ=inv_tAkl*Jij*inv_tAkl (each multiplied by Lk or by
!P...P'*Ll) with coefficients c_ij=q_ij*tr[inv_tAkl*Jij]^(-3/2),
!which is O(n^5) overall. Since sum_ij c_ij*Jij = C is an n x n
!matrix assembled in O(n^2), the whole sum is simply
!  B = inv_tAkl*C*inv_tAkl
!and the gradient contributions are (2/sqrt(pi))*Skl*B*Lk and
!(2/sqrt(pi))*Skl*(P*B*P')*Ll: everything is O(n^3).
!
!The T and V contributions are combined into a single matrix
!  Z = 12*Skl*F + (2/sqrt(pi))*Skl*B          (for the bra)
!  Z = 12*Skl*(M-K-K'+F) + (2/sqrt(pi))*Skl*B (for the ket)
!so that only one triangular product with Lk (resp. one congruence
!with P and one triangular product with Ll) is needed:
!  dHkl/dvechLk = (Hkl/Skl)*dSkl/dvechLk + vech[Z*Lk]
!  dHkl/dvechLl = (Hkl/Skl)*dSkl/dvechLl + vech[(P*Z*P')*Ll]

    cV=(TWO/Glob_SqrtPi)*Skl
    HklOverSkl=Hkl/Skl

    if (grad_k) then
      !dSkldvechLk: lower triangle of -3*Skl*inv_tAkl*Lk. The
      !symmetry of inv_tAkl and the lower triangularity of Lk make
      !the inner loop a dot product of two contiguous columns.
      temp2=-THREE*Skl
      indx=0
      do i=1,n
        do j=i,n
          indx=indx+1
          temp1=ZERO
          do k=i,n
            temp1=temp1+inv_tAkl(k,j)*Lk(k,i)
          enddo
          Dk(np+indx)=temp2*temp1
        enddo
      enddo
    endif

    if (grad_l) then
      !PT=P' is stored explicitly so that all the products with P
      !and P' below can be done with contiguous column access
      do j=1,n
        do i=1,n
          PT(i,j)=P(j,i)
        enddo
      enddo
      !calculating inv_ttAkl=P*inv_tAkl*P'
      !W1=inv_tAkl*PT
      do j=1,n
        do i=1,n
          temp1=ZERO
          do k=1,n
            temp1=temp1+inv_tAkl(k,i)*PT(k,j)
          enddo
          W1(i,j)=temp1
        enddo
      enddo
      !inv_ttAkl=PT'*W1 (only the upper triangle, then mirrored)
      do j=1,n
        do i=1,j
          temp1=ZERO
          do k=1,n
            temp1=temp1+PT(k,i)*W1(k,j)
          enddo
          inv_ttAkl(i,j)=temp1
          inv_ttAkl(j,i)=temp1
        enddo
      enddo
      !dSkldvechLl: lower triangle of -3*Skl*inv_ttAkl*Ll
      temp2=-THREE*Skl
      indx=0
      do i=1,n
        do j=i,n
          indx=indx+1
          temp1=ZERO
          do k=i,n
            temp1=temp1+inv_ttAkl(k,j)*Ll(k,i)
          enddo
          Dl(np+indx)=temp2*temp1
        enddo
      enddo
    endif

!F=K*W2' (only the upper triangle is computed, then mirrored)
    do j=1,n
      do i=1,j
        temp1=ZERO
        do k=1,n
          temp1=temp1+inv_tAkltAlM(i,k)*W2(j,k)
        enddo
        F(i,j)=temp1
        F(j,i)=temp1
      enddo
    enddo

!Assembling C=sum_ij c_ij*Jij, c_ij=q_ij*tr[inv_tAkl*Jij]^(-3/2)
!(Jii is a matrix whose only nonzero element is (i,i)=1; Jij, i/=j,
!has elements (i,i)=(j,j)=1, (i,j)=(j,i)=-1)
    do j=1,n
      do i=1,n
        Cmat(i,j)=ZERO
      enddo
    enddo
    do i=1,n
      Cmat(i,i)=Glob_ScaledPseudoChargeMatrix(0,i)*tr_inv_tAklJij32(i,i)
    enddo
    do i=1,n
      do j=i+1,n
        temp1=Glob_ScaledPseudoChargeMatrix(i,j)*tr_inv_tAklJij32(j,i)
        Cmat(i,i)=Cmat(i,i)+temp1
        Cmat(j,j)=Cmat(j,j)+temp1
        Cmat(j,i)=Cmat(j,i)-temp1
        Cmat(i,j)=Cmat(i,j)-temp1
      enddo
    enddo

    do j=1,n
      do i=1,n
        temp1=ZERO
        do k=1,n
          temp1=temp1+Cmat(k,i)*inv_tAkl(k,j)
        enddo
        W1(i,j)=temp1
      enddo
    enddo
    do j=1,n
      do i=1,j
        temp1=ZERO
        do k=1,n
          temp1=temp1+W1(k,i)*inv_tAkl(k,j)
        enddo
        Bmat(i,j)=temp1
        Bmat(j,i)=temp1
      enddo
    enddo

    if (grad_k) then
      temp2=12*Skl
      do j=1,n
        do i=1,n
          Z(i,j)=temp2*F(i,j)+cV*Bmat(i,j)
        enddo
      enddo
      indx=0
      do i=1,n
        do j=i,n
          indx=indx+1
          temp1=ZERO
          do k=i,n
            temp1=temp1+Z(k,j)*Lk(k,i)
          enddo
          Dk(indx)=HklOverSkl*Dk(np+indx)+temp1
        enddo
      enddo
    endif

    if (grad_l) then
      temp2=12*Skl
      do j=1,n
        do i=1,n
          Z(i,j)=temp2*(Glob_MassMatrix(i,j)-inv_tAkltAlM(i,j) &
                        -inv_tAkltAlM(j,i)+F(i,j))+cV*Bmat(i,j)
        enddo
      enddo
      do j=1,n
        do i=1,n
          temp1=ZERO
          do k=1,n
            temp1=temp1+Z(k,i)*PT(k,j)
          enddo
          W2(i,j)=temp1
        enddo
      enddo
      do j=1,n
        do i=1,j
          temp1=ZERO
          do k=1,n
            temp1=temp1+PT(k,i)*W2(k,j)
          enddo
          G(i,j)=temp1
          G(j,i)=temp1
        enddo
      enddo
      indx=0
      do i=1,n
        do j=i,n
          indx=indx+1
          temp1=ZERO
          do k=i,n
            temp1=temp1+G(k,j)*Ll(k,i)
          enddo
          Dl(indx)=HklOverSkl*Dl(np+indx)+temp1
        enddo
      enddo
    endif

  end subroutine MatrixElementsHS_RG_0S

  subroutine MatrixElementsAll_RG_0S(vechLk, vechLl, Pbra, Pket, &
                                       Hkl, Skl, Tkl, Vkl, rm2kl, rmkl, rkl, r2kl, deltarkl, drach_deltarkl, &
                             MVkl, drach_MVkl1, drach_MVkl2, drach_MVkl3, Darwinkl, drach_Darwinkl, OOkl, rmrmkl, del2kl, prvalkl, &
                                       wf2originkl, NumCFGridPoints, CFGrid, CFkl, NumDensGridPoints, DensGrid, Denskl, &
                                       AreCorrFuncNeeded, ArePartDensNeeded, AreMCorrFuncNeeded, AreMPartDensNeeded)
!This subroutine computes symmetry adapted matrix elements
!with two real L=0 correlated Gaussians. These matrix elements
!are used in calculations of expectation values.
!Symmetry adaption is applied to the bra and ket using permutation matrices Pbra and Pket
!
!Input:
!   vechLk, vechLl :: Arrays of length (n(n+1)/2) of exponential parameters.
!   Pbra    :: The symmetry permutation matrix of size n x n that is applied to bra
!   Pket    :: The symmetry permutation matrix of size n x n that is applied to ket
!Output (all matrix elements are computed with normalized functions):
!   Hkl             ::        Hamiltonian
!   Skl             ::        Overlap
!   Tkl      :: Kinetic energy
!   Vkl      :: Potential energy
!   rm2kl    :: 1/r_i^2, 1/r_{ij}^2
!   rmkl     :: 1/r_i, 1/r_{ij}
!   rkl      :: r_i, r_{ij}
!   r2kl     :: r_i^2, r_{ij}^2
! deltarkl   :: delta(r_i), delta(r_{ij})
! drach_deltarkl:: Drachmanized delta(r_i), delta(r_{ij})
!   MVkl     :: Mass-velocity correction (without the factor of alpha**2)
! drach_MVkl :: Drachmanized mass-velocity correction (without the factor of alpha**2)
! Darwinkl  :: Darwin correction (without the factor of alpha**2)
! drach_Darwinkl:: Drachmanized Darwin correction (without the factor of alpha**2)
!   OOkl    :: Orbit-Orbit correction (without the factor of alpha**2)
! rmrmkl    :: 1/(r_{ij}*r_{pq})
! del2kl    :: delta(r_{ij})delta(r_{pq}) when r_{ij}/=r_{pq}
! prvalkl   :: P(1/r^3_ij) - principal values of matrix element 1/r^3_ij  (appears in the Araki-Sucker term for QED correction)
! wf2originkl:: n-particle density at all-particle coalescence point (absolute square of the wave function at the origin)
!NumCFGridPoints   :: Number of grid points for correlation function calculations
!CFGrid            :: Array containing grid points where matrix elements of
!                     correlation functions should be computed
!CFkl              :: Matrix elements of correlation functions
!NumDensGridPoints :: Number of grid points for particle density calculations
!DensGrid          :: Array containing grid points where matrix elements of
!                     particle densities should be computed
!Denskl            :: Matrix elements of particle densities
!AreCorrFuncNeeded :: flag indicating whether matrix elements of correlation
!                     functions need to be computed
!ArePartDensNeeded :: flag indicating whether matrix elements of particle
!                     densities need to be computed

!Arguments
    real(wp),intent(in)   :: vechLk(Glob_np), vechLl(Glob_np)
    real(wp),intent(in)   :: Pbra(Glob_n,Glob_n),Pket(Glob_n,Glob_n)
    real(wp),intent(out)  :: Hkl,Skl,Tkl,Vkl,MVkl,drach_MVkl1,drach_MVkl2,drach_MVkl3,Darwinkl,drach_Darwinkl,OOkl
    real(wp),intent(out)  :: rm2kl(Glob_n,Glob_n),rmkl(Glob_n,Glob_n)
    real(wp),intent(out)  :: rkl(Glob_n,Glob_n),r2kl(Glob_n,Glob_n)
    real(wp),intent(out)  :: deltarkl(Glob_n,Glob_n)
    real(wp),intent(out)  :: drach_deltarkl(Glob_n,Glob_n)
    real(wp),intent(out)  :: prvalkl(Glob_n,Glob_n)
    real(wp),intent(out)  :: rmrmkl(Glob_n,Glob_n,Glob_n,Glob_n)
    real(wp),intent(out)  :: del2kl(Glob_n,Glob_n,Glob_n,Glob_n)
    real(wp),intent(out)  :: wf2originkl
    integer,intent(in)       :: NumCFGridPoints,NumDensGridPoints
    real(wp),intent(in)   :: CFGrid(NumCFGridPoints),DensGrid(NumDensGridPoints)
    real(wp),intent(out)  :: CFkl(Glob_n*(Glob_n+1)/2,NumCFGridPoints)
    real(wp),intent(out)  :: Denskl(Glob_n+1,NumDensGridPoints)
    logical,intent(in)       :: AreCorrFuncNeeded,ArePartDensNeeded,AreMCorrFuncNeeded, AreMPartDensNeeded

!Parameters (These are needed to declare static arrays. Using static
!arrays makes the function call a little faster in comparison with
!the case when arrays are dynamically allocated in stack)
    integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles

!Local variables
    integer           n, np
    real(wp)       Lk(nn,nn), Ll(nn,nn)
    real(wp)       tAk(nn,nn), tAl(nn,nn), tAkl(nn,nn)
    real(wp)       inv_tAk(nn,nn), inv_tAl(nn,nn), inv_tAkl(nn,nn)
    real(wp)       inv_tAkltAlM(nn,nn), inv_invtAkinvtAl(nn,nn)
    real(wp)       AtAl(nn,nn), GAl(nn,nn), Cmat(nn,nn)
    real(wp)       AXs(nn,nn), YAk(nn,nn), AYsk(nn,nn)
    real(wp)       W1(nn,nn), W2(nn,nn), W3(nn,nn), W4(nn,nn), W5(nn,nn)
    real(wp)       sk(nn), sl(nn), Ask(nn), Asl(nn)
    real(wp)       temp1, temp2, temp3, temp4, temp5, temp6, temp7
    real(wp)       tr1, tr2, tr3, tr4, tr5
    real(wp)       alpha, beta, gamm, delt, zeta
    real(wp)       trMtAl, trUMtAl, trXAl, trAXs, trYAk, trAYs, trAYsAXs
    real(wp)       det_tAkl, det_tAk, det_tAl, det_invtAkinvtAl
    integer           i,j,k,indx,p,q
    real(wp)       TrAJ(nn,nn),sqrtTrAJ(nn,nn),TrAJAJ(nn,nn,nn,nn),MTrAJ(nn,nn),sqrtMTrAJ(nn,nn)
    real(wp)       Mass_For_Darwin(0:nn)
    real(wp)       V2kl, MSkl

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

!Then we permute elements of Ak and Al to account for
!the action of the permutation matrix
!  tAl=Pket'*Al*Pket
!  tAk=Pbra'*Ak*Pbra
!We also form matrix tAkl=tAk+tAl
    do i=1,n
      do j=1,n
        temp1=ZERO
        temp2=ZERO
        do k=1,n
          temp1=temp1+Pket(k,j)*tAl(k,i)
          temp2=temp2+tAk(j,k)*Pbra(k,i)
        enddo
        W1(j,i)=temp1
        W2(j,i)=temp2
      enddo
    enddo
!tAl=W1*Pket
!tAk=Pbra'*W2
    do i=1,n
      do j=i,n
        temp1=ZERO
        temp2=ZERO
        do k=1,n
          temp1=temp1+W1(j,k)*Pket(k,i)
          temp2=temp2+Pbra(k,j)*W2(k,i)
        enddo
        tAl(j,i)=temp1
        tAl(i,j)=temp1
        tAk(j,i)=temp2
        tAk(i,j)=temp2
        tAkl(j,i)=temp1+temp2
        tAkl(i,j)=temp1+temp2
      enddo
    enddo

!I will delete not-necessary parts of the code
!The determinants of Lk and Ll are just
!the products of their diagonal elements
!det_Lk=ONE
!det_Ll=ONE
!do i=1,n
!  det_Lk=det_Lk*Lk(i,i)
!  det_Ll=det_Ll*Ll(i,i)
!enddo

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

!Do Cholesky factorization of tAk and tAk and then invert
    if (AreMCorrFuncNeeded.or.AreMPartDensNeeded) then
      det_tAk=ONE
      det_tAl=ONE
      do i=1,n
        do j=i,n
          temp1=tAk(i,j)
          temp2=tAl(i,j)
          do k=i-1,1,-1
            temp1=temp1-W1(i,k)*W1(j,k)
            temp2=temp2-W2(i,k)*W2(j,k)
          enddo
          if (i==j) then
            W1(i,i)=sqrt(temp1)
            det_tAk=det_tAk*temp1
            W2(i,i)=sqrt(temp2)
            det_tAl=det_tAl*temp2
          else
            W1(j,i)=temp1/W1(i,i)
            W1(i,j)=ZERO
            W2(j,i)=temp2/W2(i,i)
            W2(i,j)=ZERO
          endif
        enddo
      enddo

!Inverting tAk and tAl using its Cholesky factors (stored in W1, W2)
!and placing the result into inv_tAk, inv_tAl
      do i=1,n
        W1(i,i)=ONE/W1(i,i)
        W2(i,i)=ONE/W2(i,i)
        do j=i+1,n
          temp1=ZERO
          temp2=ZERO
          do k=i,j-1
            temp1=temp1-W1(j,k)*W1(k,i)
            temp2=temp2-W2(j,k)*W2(k,i)
          enddo
          W1(j,i)=temp1/W1(j,j)
          W2(j,i)=temp2/W2(j,j)
        enddo
      enddo

      do i=1,n
        do j=i,n
          temp1=ZERO
          temp2=ZERO
          do k=j,n
            temp1=temp1+W1(k,i)*W1(k,j)
            temp2=temp2+W2(k,i)*W2(k,j)
          enddo
          inv_tAk(i,j)=temp1
          inv_tAl(i,j)=temp2
          inv_tAk(j,i)=temp1
          inv_tAl(j,i)=temp2
        enddo
      enddo

!Now calculate inv_invtAkinvtAl
      det_invtAkinvtAl=ONE
      do i=1,n
        do j=i,n
          temp1=inv_tAk(i,j)+inv_tAl(i,j)
          do k=i-1,1,-1
            temp1=temp1-W1(i,k)*W1(j,k)
          enddo
          if (i==j) then
            W1(i,i)=sqrt(temp1)
            det_invtAkinvtAl=det_invtAkinvtAl*temp1
          else
            W1(j,i)=temp1/W1(i,i)
            W1(i,j)=ZERO
          endif
        enddo
      enddo

!Inverting invtAk+invtAl
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
          inv_invtAkinvtAl(i,j)=temp1
          inv_invtAkinvtAl(j,i)=temp1
        enddo
      enddo
    endif

!Evaluating overlap
!temp2=abs(det_Ll*det_Lk)
!temp1=temp2/det_tAkl
!Skl=Glob_2Raised3n2*temp1*sqrt(temp1)
!wf2originkl=Glob_2Raised3n2*(temp2*sqrt(temp2))/(PI**(THREE*n/TWO))
    wf2originkl=ONE
    Skl=Glob_PiRaised3n2/(det_tAkl*sqrt(det_tAkl))  !new line

    if(AreMCorrFuncNeeded.or.AreMPartDensNeeded) then
      temp1=1/det_tAk/det_tAl/det_invtAkinvtAl
      MSkl=Glob_PiRaised3n2*temp1*sqrt(temp1)
    endif

!Doing multiplication AtAl=inv_tAkl*tAl (kept: several sections below
!rely on it and on the identity inv_tAkl*tAk = I - AtAl)
    do i=1,n
      do j=1,n
        temp1=ZERO
        do k=1,n
          temp1=temp1+inv_tAkl(j,k)*tAl(k,i)
        enddo
        AtAl(j,i)=temp1
      enddo
    enddo

!Doing multiplication inv_tAkltAlM=inv_tAkl*tAl*M=AtAl*M
    do i=1,n
      do j=1,n
        temp1=ZERO
        do k=1,n
          temp1=temp1+AtAl(j,k)*Glob_MassMatrix(k,i)
        enddo
        inv_tAkltAlM(j,i)=temp1
      enddo
    enddo

!GAl = tAl*inv_tAkl*tAl = tAl*AtAl (symmetric). Its elements are the
!bilinear forms tAl(:,i)'*inv_tAkl*tAl(:,j) needed by the orbit-orbit
!and mass-velocity sections below.
    do i=1,n
      do j=i,n
        temp1=ZERO
        do k=1,n
          temp1=temp1+tAl(j,k)*AtAl(k,i)
        enddo
        GAl(j,i)=temp1
        GAl(i,j)=temp1
      enddo
    enddo

!Computing kinetic energy, Tkl=tr[inv_tAkltAlM*Ak]
    Tkl=ZERO
    do i=1,n
      temp1=ZERO
      do k=1,n
        temp1=temp1+inv_tAkltAlM(i,k)*tAk(k,i)
      enddo
      Tkl=Tkl+temp1
    enddo
    Tkl=SIX*Skl*Tkl

!Evaluating potential energy, Vkl,
!(1/r_{ij})_kl, (r_{ij})_kl, (r_{ij}^2)_kl
!and delta(r_{ij})_kl
    temp1=(TWO/Glob_SqrtPi)*Skl
    temp2=THREEHALF*Skl
    temp3=Skl/(Glob_Pi*Glob_SqrtPi)
    Vkl=ZERO
    do i=1,n
      TrAJ(i,i)=inv_tAkl(i,i)
      sqrtTrAJ(i,i)=sqrt(TrAJ(i,i))
      temp4=sqrtTrAJ(i,i)
      temp5=TrAJ(i,i)
      rmkl(i,i)=temp1/temp4
      rkl(i,i)=temp1*temp4
      r2kl(i,i)=temp2*temp5
      deltarkl(i,i)=temp3/(temp4*temp5)
      prvalkl(i,i)=(temp1/(temp4*temp5))*(Glob_EulerConst+log(temp5))
      Vkl=Vkl+Glob_ScaledPseudoChargeMatrix(i,0)*rmkl(i,i)
    enddo
    do i=1,n
      do j=i+1,n
        TrAJ(i,j)=inv_tAkl(i,i)+inv_tAkl(j,j)-inv_tAkl(j,i)-inv_tAkl(j,i)
        TrAJ(j,i)=TrAJ(i,j)
        sqrtTrAJ(j,i)=sqrt(TrAJ(j,i))
        sqrtTrAJ(i,j)=sqrtTrAJ(j,i)
        temp4=sqrtTrAJ(j,i)
        temp5=TrAJ(j,i)
        rmkl(i,j)=temp1/temp4
        rmkl(j,i)=rmkl(i,j)
        rkl(i,j)=temp1*temp4
        rkl(j,i)=rkl(i,j)
        r2kl(i,j)=temp2*temp5
        r2kl(j,i)=r2kl(i,j)
        deltarkl(i,j)=temp3/(temp4*temp5)
        deltarkl(j,i)=deltarkl(i,j)
        prvalkl(i,j)=(temp1/(temp4*temp5))*(Glob_EulerConst+log(temp5))
        prvalkl(j,i)=prvalkl(i,j)
        Vkl=Vkl+Glob_ScaledPseudoChargeMatrix(i,j)*rmkl(i,j)
      enddo
    enddo
    Hkl=Tkl+Vkl

    if (AreMCorrFuncNeeded) then
      do i=1,n
        MTrAJ(i,i)=inv_invtAkinvtAl(i,i)*4
        sqrtMTrAJ(i,i)=sqrt(MTrAJ(i,i))
      enddo
      do i=1,n
        do j=i+1,n
          MTrAJ(i,j)=(inv_invtAkinvtAl(i,i)+inv_invtAkinvtAl(j,j)-inv_invtAkinvtAl(j,i)-inv_invtAkinvtAl(j,i))*4
          MTrAJ(j,i)=MTrAJ(i,j)
          sqrtMTrAJ(j,i)=sqrt(MTrAJ(j,i))
          sqrtMTrAJ(i,j)=sqrtMTrAJ(j,i)
        enddo
      enddo
    end if

!Evaluating tr[inv_tAkl Jij inv_tAkl Jpq]
    do i=1,n
      temp1=inv_tAkl(i,i)*inv_tAkl(i,i)
      TrAJAJ(i,i,i,i)=temp1
      do p=i+1,n
        do q=p+1,n
          temp2=inv_tAkl(p,i)-inv_tAkl(q,i)
          temp1=temp2*temp2
          TrAJAJ(i,i,p,q)=temp1
          TrAJAJ(i,i,q,p)=temp1
          TrAJAJ(p,q,i,i)=temp1
          TrAJAJ(q,p,i,i)=temp1
        enddo
      enddo
      do j=i+1,n
        temp1=inv_tAkl(j,i)*inv_tAkl(j,i)
        TrAJAJ(j,j,i,i)=temp1
        TrAJAJ(i,i,j,j)=temp1
        do p=i,n
          temp2=inv_tAkl(p,i)-inv_tAkl(p,j)
          temp1=temp2*temp2
          TrAJAJ(i,j,p,p)=temp1
          TrAJAJ(j,i,p,p)=temp1
          TrAJAJ(p,p,i,j)=temp1
          TrAJAJ(p,p,j,i)=temp1
          do q=p+1,n
            temp2=inv_tAkl(q,i)-inv_tAkl(p,i)-inv_tAkl(q,j)+inv_tAkl(p,j)
            temp1=temp2*temp2
            TrAJAJ(i,j,p,q)=temp1
            TrAJAJ(j,i,p,q)=temp1
            TrAJAJ(i,j,q,p)=temp1
            TrAJAJ(j,i,q,p)=temp1
            TrAJAJ(p,q,i,j)=temp1
            TrAJAJ(p,q,j,i)=temp1
            TrAJAJ(q,p,i,j)=temp1
            TrAJAJ(q,p,j,i)=temp1
          enddo
        enddo
      enddo
    enddo

!This is a slow old version of the previous loop that computes TrAJAJ
!do i=1,n
!  do j=1,n
!    do p=1,n
!      do q=1,n
!        W1(1:n,1:n)=ZERO; W1(i,j)=-ONE; W1(j,i)=-ONE; W1(i,i)=ONE; W1(j,j)=ONE
!        W2(1:n,1:n)=ZERO; W2(p,q)=-ONE; W2(q,p)=-ONE; W2(p,p)=ONE; W2(q,q)=ONE;
!        W3(1:n,1:n)=matmul(inv_tAkl(1:n,1:n),W1(1:n,1:n))
!        W4(1:n,1:n)=matmul(inv_tAkl(1:n,1:n),W2(1:n,1:n))
!        W5(1:n,1:n)=matmul(W3(1:n,1:n),W4(1:n,1:n))
!        TrAJAJ(i,j,p,q)=trace(n,W5)
!      enddo
!    enddo
!  enddo
!enddo

!Evaluating (1/r_{ij}*1/r_{pq}))_kl and <delta(r_{ij})delta(r_{pq})>_kl
!in one fused loop (they run over the same index combinations and share
!TrAJ/TrAJAJ). Within the loop ranges (i<=j, p>=i, q>=p) the "same pair"
!condition (p==i.and.q==j).or.(p==j.and.q==i) reduces to p==i.and.q==j.
!Those are exactly the elements of del2kl that must be zero, so writing
!them explicitly makes the full-array zeroing of del2kl unnecessary.
    temp1=4*Skl/Glob_Pi
    temp6=2*Skl
    temp7=Skl/(Glob_Pi**3)
    do i=1,n
      do j=i,n
        do p=i,n
          do q=p,n
            if ((p==i).and.(q==j)) then
              temp2=temp6/TrAJ(i,j)
              rmrmkl(i,j,p,q)=temp2
              rmrmkl(j,i,p,q)=temp2
              rmrmkl(i,j,q,p)=temp2
              rmrmkl(j,i,q,p)=temp2
              rmrmkl(p,q,i,j)=temp2
              rmrmkl(p,q,j,i)=temp2
              rmrmkl(q,p,i,j)=temp2
              rmrmkl(q,p,j,i)=temp2
              del2kl(i,j,p,q)=ZERO
              del2kl(i,j,q,p)=ZERO
              del2kl(j,i,p,q)=ZERO
              del2kl(j,i,q,p)=ZERO
              del2kl(p,q,i,j)=ZERO
              del2kl(p,q,j,i)=ZERO
              del2kl(q,p,i,j)=ZERO
              del2kl(q,p,j,i)=ZERO
            else
              temp3=TrAJAJ(i,j,p,q)
              if (temp3>ZERO) then
                temp4=sqrt(temp3)
                temp5=temp4/(sqrtTrAJ(i,j)*sqrtTrAJ(p,q))
                temp2=temp1*asin(temp5)/temp4
              else
                temp2=temp1/(sqrtTrAJ(i,j)*sqrtTrAJ(p,q))
              endif
              rmrmkl(i,j,p,q)=temp2
              rmrmkl(j,i,p,q)=temp2
              rmrmkl(i,j,q,p)=temp2
              rmrmkl(j,i,q,p)=temp2
              rmrmkl(p,q,i,j)=temp2
              rmrmkl(p,q,j,i)=temp2
              rmrmkl(q,p,i,j)=temp2
              rmrmkl(q,p,j,i)=temp2
              temp2=trAJ(i,j)*trAJ(p,q)-temp3
              temp3=temp7/(temp2*sqrt(temp2))
              del2kl(i,j,p,q)=temp3
              del2kl(i,j,q,p)=temp3
              del2kl(j,i,p,q)=temp3
              del2kl(j,i,q,p)=temp3
              del2kl(p,q,i,j)=temp3
              del2kl(p,q,j,i)=temp3
              del2kl(q,p,i,j)=temp3
              del2kl(q,p,j,i)=temp3
            endif
          enddo
        enddo
      enddo
    enddo

!Extracting rm2kl from rmrmkl
    do i=1,n
      do j=1,n
        rm2kl(j,i)=rmrmkl(j,i,j,i)
      enddo
    enddo

!Evaluating drachmanized delta(r_{ij})_kl
!W3(i,j) will contain all matrix elements r'(tAk*M*tAl)r / rij.
!With X=tAk*M*tAl one has
!  inv_tAkl*X*inv_tAkl = (I-AtAl)*M*AtAl'   (using inv_tAkl*tAkl = I)
!  tr[inv_tAkl*X]      = tr[M*tAl] - tr[(AtAl*M)*tAl]
!so instead of building X and calling ME_rXr_over_rij_all (which spends
!O(n^2) per pair recomputing rows of inv_tAkl*X) everything follows from
!two matrix products, O(n^2) trace contractions, and O(1) per pair.
    do i=1,n
      do j=1,n
        temp1=ZERO
        do k=1,n
          temp1=temp1+Glob_MassMatrix(j,k)*AtAl(i,k)
        enddo
        W1(j,i)=temp1
      enddo
    enddo
    do i=1,n
      do j=1,n
        temp1=ZERO
        do k=1,n
          temp1=temp1+AtAl(j,k)*W1(k,i)
        enddo
        Cmat(j,i)=W1(j,i)-temp1
      enddo
    enddo
    trMtAl=ZERO
    trUMtAl=ZERO
    do i=1,n
      do k=1,n
        trMtAl=trMtAl+Glob_MassMatrix(i,k)*tAl(k,i)
        trUMtAl=trUMtAl+inv_tAkltAlM(i,k)*tAl(k,i)
      enddo
    enddo
    temp2=3*(trMtAl-trUMtAl)
    do i=1,n
      do j=i,n
        if (i==j) then
          temp3=Cmat(i,i)
        else
          temp3=Cmat(i,i)+Cmat(j,j)-Cmat(i,j)-Cmat(j,i)
        endif
        W3(j,i)=rmkl(j,i)*(temp2-temp3/TrAJ(j,i))*ONEHALF
        W3(i,j)=W3(j,i)
      enddo
    enddo
!Loop that computes all drachmanized delta(r_{ij})_kl  as well as V^2_kl
    V2kl=ZERO
    do p=1,n
      do q=p,n
        temp1=ZERO
        do i=1,n
          temp1=temp1+Glob_ScaledPseudoChargeMatrix(0,i)*rmrmkl(p,q,i,i)
          do j=i+1,n
            temp1=temp1+Glob_ScaledPseudoChargeMatrix(i,j)*rmrmkl(p,q,i,j)
          enddo
        enddo
        if (p==q) then
          temp4=2*Glob_Pi*Glob_MassMatrix(p,p)
          temp5=Glob_ScaledPseudoChargeMatrix(0,p)
        else
          temp4=2*Glob_Pi*(Glob_MassMatrix(p,p)+Glob_MassMatrix(q,q) &
                      -Glob_MassMatrix(q,p)-Glob_MassMatrix(q,p))
          temp5=Glob_ScaledPseudoChargeMatrix(p,q)
        endif
        drach_deltarkl(q,p)=(Glob_CurrEnergy*rmkl(q,p)-temp1-4*W3(q,p))/temp4
        drach_deltarkl(p,q)=drach_deltarkl(q,p)
        V2kl=V2kl+temp5*temp1
      enddo
    enddo

!Evaluating Orbit-Orbit (OO) matrix element (without the factor of alpha**2)
!The matrices the original code assembled for every pair,
!  W1 = c*b' + b*c' + tr1*a*d' + d*(tr1*b' + tr3*c')
!  W2 = a*b'
!  W3 = d*c'
!with a=tAkl(:,j), b=tAl(:,j), c=tAl(:,i), and d=e_j (first sum) or
!d=e_j-e_i (second sum), are sums of rank-one terms u*v'. Every trace
!formed from them inside ME_rXr_over_rij and ME_rXr_rYr_over_rij
!therefore reduces to O(1) scalars:
!  tr[inv_tAkl*u*v'] = v'*inv_tAkl*u
!  d'*inv_tAkl*u*v'*inv_tAkl*d = (d'*inv_tAkl*u)*(v'*inv_tAkl*d)
!and (with J=d*d', which equals Jjj resp. Jij up to the irrelevant sign
!of d) all the needed building blocks are
!  inv_tAkl*a = e_j (exactly, since inv_tAkl*tAkl = I)
!  alpha = b'*inv_tAkl*c = GAl(j,i)
!  beta  = b'*inv_tAkl*d,  gamm = c'*inv_tAkl*d  (elements of AtAl)
!  delt  = d'*inv_tAkl*d = TrAJ,  a'*inv_tAkl*b = tAl(j,j),
!  a'*inv_tAkl*c = tAl(j,i) = tr1,  a'*inv_tAkl*d = 1.
!In terms of these the two integrals are
!  temp1 = rm*(3*(2*alpha+s) - (2*beta*gamm+delt*s)/delt)/2,
!          s = tr1 + tr1*beta + tr3*gamm
!  temp2 = rm*(9*gamm*zeta + 3*(beta*tr1+alpha)
!          - (3*(gamm*beta+zeta*delt*gamm)
!             + (2*beta*gamm+alpha*delt+beta*tr1*delt))/delt
!          + 3*beta*delt*gamm/delt**2)/4,  zeta = tAl(j,j)
!which makes the whole OO evaluation O(n^2) instead of O(n^5).
    OOkl=ZERO
!First double loop for OO (integral pair index (j,j), d=e_j)
    do i=1,n
      do j=1,n
        tr1=tAl(j,i)
        tr3=3*tAl(j,j)
        alpha=GAl(j,i)
        beta=AtAl(j,j)
        gamm=AtAl(j,i)
        delt=TrAJ(j,j)
        zeta=tAl(j,j)
        temp4=rmkl(j,j)
        temp5=tr1+tr1*beta+tr3*gamm
        temp1=temp4*(3*(2*alpha+temp5)-(2*beta*gamm+delt*temp5)/delt)*ONEHALF
        tr4=beta
        tr5=delt*gamm
        temp2=temp4*( 9*gamm*zeta + 3*(beta*tr1+alpha) &
             -(3*(gamm*tr4+zeta*tr5)+(2*beta*gamm+alpha*delt+beta*tr1*delt))/delt &
             +3*tr4*tr5/(delt*delt) )/4
        temp3=-12*tr1*temp4+4*temp1-8*temp2
        OOkl=OOkl-temp3*Glob_ScaledPseudoChargeMatrix(j,0)/Glob_Mass(j+1)
      enddo
    enddo
    OOkl=OOkl/Glob_Mass(1)

!Second double loop for OO (integral pair (i,j), d=e_j-e_i)
    do i=1,n
      do j=i+1,n
        tr1=tAl(j,i)
        tr3=3*tAl(j,j)
        alpha=GAl(j,i)
        beta=AtAl(j,j)-AtAl(i,j)
        gamm=AtAl(j,i)-AtAl(i,i)
        delt=TrAJ(i,j)
        zeta=tAl(j,j)
        temp4=rmkl(i,j)
        temp5=tr1+tr1*beta+tr3*gamm
        temp1=temp4*(3*(2*alpha+temp5)-(2*beta*gamm+delt*temp5)/delt)*ONEHALF
        tr4=beta
        tr5=delt*gamm
        temp2=temp4*( 9*gamm*zeta + 3*(beta*tr1+alpha) &
             -(3*(gamm*tr4+zeta*tr5)+(2*beta*gamm+alpha*delt+beta*tr1*delt))/delt &
             +3*tr4*tr5/(delt*delt) )/4
        temp3=-12*tr1*temp4+4*temp1-8*temp2
        OOkl=OOkl+ &
              temp3*Glob_ScaledPseudoChargeMatrix(i,j)/(Glob_Mass(i+1)*Glob_Mass(j+1))
      enddo
    enddo
    OOkl=OOkl/2

!Evaluating mass-velocity matrix element.
!All matrices involved are rank one: with J the all-ones matrix,
!tAk*J*tAk = sk*sk' where sk = tAk*(1,..,1)' (and likewise for tAl),
!and tAk*Jii*tAk = tAk(:,i)*tAk(i,:). Hence
!  tr[inv_tAkl*u*u']            = u'*inv_tAkl*u
!  tr[inv_tAkl*u*u'*inv_tAkl*v*v'] = (u'*inv_tAkl*v)**2
!For the per-particle terms the dot products follow from
!inv_tAkl*tAk = I - AtAl and GAl = tAl*inv_tAkl*tAl:
!  tAk(:,i)'*inv_tAkl*tAk(:,i) = tAk(i,i) - tAl(i,i) + GAl(i,i)
!  tAl(:,i)'*inv_tAkl*tAl(:,i) = GAl(i,i)
!  tAk(:,i)'*inv_tAkl*tAl(:,i) = tAl(i,i) - GAl(i,i)
!This makes the whole section O(n^2) instead of O(n^4).
    do p=1,n
      temp1=ZERO
      temp2=ZERO
      do q=1,n
        temp1=temp1+tAk(q,p)
        temp2=temp2+tAl(q,p)
      enddo
      sk(p)=temp1
      sl(p)=temp2
    enddo
    tr1=ZERO
    tr2=ZERO
    do p=1,n
      tr1=tr1+sk(p)
      tr2=tr2+sl(p)
    enddo
    do p=1,n
      temp1=ZERO
      temp2=ZERO
      do q=1,n
        temp1=temp1+inv_tAkl(p,q)*sk(q)
        temp2=temp2+inv_tAkl(p,q)*sl(q)
      enddo
      Ask(p)=temp1
      Asl(p)=temp2
    enddo
    tr4=ZERO
    tr5=ZERO
    temp3=ZERO
    do p=1,n
      tr4=tr4+sk(p)*Ask(p)
      tr5=tr5+sl(p)*Asl(p)
      temp3=temp3+sk(p)*Asl(p)
    enddo
    temp3=temp3*temp3
    temp2=tr1*tr5+tr2*tr4

    MVkl=36*Skl*(tr4*tr5 + (TWO/THREE)*temp3 - temp2 + tr1*tr2) &
          /(Glob_Mass(1)*Glob_Mass(1)*Glob_Mass(1))

!sum for mass-velocity from 1 to n
    do i=1,n
      tr1=tAk(i,i)
      tr2=tAl(i,i)
      tr5=GAl(i,i)
      tr4=tAk(i,i)-tAl(i,i)+tr5
      temp3=(tAl(i,i)-tr5)*(tAl(i,i)-tr5)
      temp2=tr1*tr5+tr2*tr4
      MVkl=MVkl+36*Skl*(tr4*tr5 + (TWO/THREE)*temp3 - temp2 + tr1*tr2) &
            /(Glob_Mass(i+1)*Glob_Mass(i+1)*Glob_Mass(i+1))
    enddo

    MVkl=-MVkl/8

!Evaluation of the drachmanized mass-velocity.
!The original helpers myME_dXd_dYd/myME_dXd/myME_over_rij_dXd all
!sandwich a fixed matrix between tAl (X side) or tAk (Y side) and then
!form traces with inv_tAkl. Using inv_tAkl*tAl = AtAl and
!inv_tAkl*tAk = I - AtAl those traces are formed here directly, the
!X = Glob_dmvM pieces are shared between the two myME_dXd_dYd calls,
!and the per-pair myME_over_rij_dXd calls (which redid four O(n^3)
!products for the same X = Glob_dmvB every time) collapse into two
!products plus O(1) work per pair. The commonFactor of the helpers,
!Glob_PiRaised3n2/(det_tAkl*sqrt(det_tAkl)), is Skl.
!X side (X = Glob_dmvM): W1 = X*tAl, AXs = inv_tAkl*(tAl*X*tAl) = AtAl*W1
    do i=1,n
      do j=1,n
        temp1=ZERO
        do k=1,n
          temp1=temp1+Glob_dmvM(j,k)*tAl(k,i)
        enddo
        W1(j,i)=temp1
      enddo
    enddo
    trXAl=ZERO
    do i=1,n
      trXAl=trXAl+W1(i,i)
    enddo
    do i=1,n
      do j=1,n
        temp1=ZERO
        do k=1,n
          temp1=temp1+AtAl(j,k)*W1(k,i)
        enddo
        AXs(j,i)=temp1
      enddo
    enddo
    trAXs=ZERO
    do i=1,n
      trAXs=trAXs+AXs(i,i)
    enddo
!Y side (Y = Glob_dmvM): YAk = Y*tAk, AYsk = inv_tAkl*(tAk*Y*tAk) = (I-AtAl)*YAk
    do i=1,n
      do j=1,n
        temp1=ZERO
        do k=1,n
          temp1=temp1+Glob_dmvM(j,k)*tAk(k,i)
        enddo
        YAk(j,i)=temp1
      enddo
    enddo
    trYAk=ZERO
    do i=1,n
      trYAk=trYAk+YAk(i,i)
    enddo
    do i=1,n
      do j=1,n
        temp1=ZERO
        do k=1,n
          temp1=temp1+AtAl(j,k)*YAk(k,i)
        enddo
        AYsk(j,i)=YAk(j,i)-temp1
      enddo
    enddo
    trAYs=ZERO
    do i=1,n
      trAYs=trAYs+AYsk(i,i)
    enddo
    trAYsAXs=ZERO
    do i=1,n
      do j=1,n
        trAYsAXs=trAYsAXs+AYsk(i,j)*AXs(j,i)
      enddo
    enddo
    temp1 = 16._wp*Skl*(NINE/FOUR*trAXs*trAYs+THREE/TWO*trAYsAXs) &
            - 36._wp*Skl*trXAl*trAYs - 36._wp*Skl*trYAk*trAXs &
            + 36._wp*Skl*trYAk*trXAl &
            - V2kl - Glob_CurrEnergy*Glob_CurrEnergy*Skl + 2*Glob_CurrEnergy*Vkl
    if (.not. Glob_ArePseudoParticleMassesTheSame) then
!Second myME_dXd_dYd call: only the Y side changes (Y = Glob_dmvMB)
      do i=1,n
        do j=1,n
          temp2=ZERO
          do k=1,n
            temp2=temp2+Glob_dmvMB(j,k)*tAk(k,i)
          enddo
          YAk(j,i)=temp2
        enddo
      enddo
      trYAk=ZERO
      do i=1,n
        trYAk=trYAk+YAk(i,i)
      enddo
      do i=1,n
        do j=1,n
          temp2=ZERO
          do k=1,n
            temp2=temp2+AtAl(j,k)*YAk(k,i)
          enddo
          AYsk(j,i)=YAk(j,i)-temp2
        enddo
      enddo
      trAYs=ZERO
      do i=1,n
        trAYs=trAYs+AYsk(i,i)
      enddo
      trAYsAXs=ZERO
      do i=1,n
        do j=1,n
          trAYsAXs=trAYsAXs+AYsk(i,j)*AXs(j,i)
        enddo
      enddo
      temp2 = 16._wp*Skl*(NINE/FOUR*trAXs*trAYs+THREE/TWO*trAYsAXs) &
              - 36._wp*Skl*trXAl*trAYs - 36._wp*Skl*trYAk*trAXs &
              + 36._wp*Skl*trYAk*trXAl &
              - V2kl - Glob_CurrEnergy*Glob_CurrEnergy*Skl + 2*Glob_CurrEnergy*Vkl
!myME_dXd(Glob_dmvB): 6*Skl*(tr[inv_tAkl*tAl*B*tAl] - tr[B*tAl]).
!W1 = B*tAl; the first trace is the O(n^2) contraction of AtAl with W1.
      do i=1,n
        do j=1,n
          temp3=ZERO
          do k=1,n
            temp3=temp3+Glob_dmvB(j,k)*tAl(k,i)
          enddo
          W1(j,i)=temp3
        enddo
      enddo
      temp5=ZERO
      temp6=ZERO
      do i=1,n
        temp5=temp5+W1(i,i)
        do j=1,n
          temp6=temp6+AtAl(i,j)*W1(j,i)
        enddo
      enddo
      temp2 = temp2 + Glob_CurrEnergy*SIX*Skl*(temp6-temp5)
!All-pairs part of myME_over_rij_dXd(Glob_dmvB,p,q,...):
!inv_tAkl*(tAl*B*tAl)*inv_tAkl = AtAl*B*AtAl' = AtAl*W4 with
!W4 = B*AtAl'; per pair only gamma=1/sqrtTrAJ and an O(1) contraction
!of W5 = AtAl*W4 remain. trAXs-trXAl of the helper is temp6-temp5.
      do i=1,n
        do j=1,n
          temp3=ZERO
          do k=1,n
            temp3=temp3+Glob_dmvB(j,k)*AtAl(i,k)
          enddo
          W4(j,i)=temp3
        enddo
      enddo
      do i=1,n
        do j=1,n
          temp3=ZERO
          do k=1,n
            temp3=temp3+AtAl(j,k)*W4(k,i)
          enddo
          W5(j,i)=temp3
        enddo
      enddo
      temp4=SIX*TWO*Skl/Glob_SqrtPi
      do i=1,n
        do j=i,n
          if (i==j) then
            temp3=W5(i,i)
            temp7=Glob_ScaledPseudoChargeMatrix(0,i)
          else
            temp3=W5(i,i)+W5(j,j)-W5(i,j)-W5(j,i)
            temp7=Glob_ScaledPseudoChargeMatrix(i,j)
          endif
          temp2=temp2-temp7*temp4*( (temp6-temp5)/sqrtTrAJ(j,i) &
                -(ONE/THREE)*temp3/(sqrtTrAJ(j,i)*TrAJ(j,i)) )
        enddo
      enddo
    endif
    drach_MVkl1 = temp1*Glob_dmva21 + MVkl
    drach_MVkl2 = temp1*Glob_dmva22 + MVkl
    drach_MVkl3 = temp2*Glob_dmva22 + MVkl

!Evaluation of the Darwin correction
    Mass_For_Darwin(0)=Glob_Mass(1)
    Mass_For_Darwin(1:n)=Glob_Mass(2:n+1)
    Darwinkl=ZERO
    do i=1,n
      Darwinkl=Darwinkl+(   &
                ONE/(Mass_For_Darwin(0)*Mass_For_Darwin(0)) &
                +ONE/(Mass_For_Darwin(i)*Mass_For_Darwin(i)) &
                )*Glob_ScaledPseudoChargeMatrix(0,i)*deltarkl(i,i)
    enddo
    do i=1,n
      do j=1,n
        if(j/=i) then
          Darwinkl=Darwinkl+   &
                    ONE/(Mass_For_Darwin(i)*Mass_For_Darwin(i)) &
                    *Glob_ScaledPseudoChargeMatrix(i,j)*deltarkl(i,j)
        endif
      enddo
    enddo
    Darwinkl=-Darwinkl*Glob_Pi/2
!Evaluation of the drachmanized Darwin correction
    drach_Darwinkl=ZERO
    do i=1,n
      drach_Darwinkl=drach_Darwinkl+(   &
                      ONE/(Mass_For_Darwin(0)*Mass_For_Darwin(0)) &
                      +ONE/(Mass_For_Darwin(i)*Mass_For_Darwin(i)) &
                      )*Glob_ScaledPseudoChargeMatrix(0,i)*drach_deltarkl(i,i)
    enddo
    do i=1,n
      do j=1,n
        if(j/=i) then
          drach_Darwinkl=drach_Darwinkl+   &
                          ONE/(Mass_For_Darwin(i)*Mass_For_Darwin(i)) &
                          *Glob_ScaledPseudoChargeMatrix(i,j)*drach_deltarkl(i,j)
        endif
      enddo
    enddo
    drach_Darwinkl=-drach_Darwinkl*Glob_Pi/2

!Evaluation of correlation functions
    if (AreCorrFuncNeeded) then
      temp1=Skl/(Glob_Pi*Glob_SqrtPi)
      p=0
      do i=1,n
        do j=i,n
          p=p+1
          temp3=temp1/(sqrtTrAJ(j,i)*TrAJ(j,i))
          temp5=ONE/TrAJ(j,i)
          do k=1,NumCFGridPoints
            temp2=CFGrid(k)*CFGrid(k)
            CFkl(p,k)=temp3*exp(-temp2*temp5)
            !CFkl(p,k)=temp2*temp3*exp(-temp2*temp5)  !Multiplied by \xi^2
          enddo
        enddo
      enddo
    endif

    if (ArePartDensNeeded) then
      temp1=Skl/(Glob_Pi*Glob_SqrtPi)
      do i=1,n+1
        temp3=ZERO
        do p=1,n
          temp3=temp3+Glob_bvc(p,i)*Glob_bvc(p,i)*inv_tAkl(p,p)
          do q=p+1,n
            temp3=temp3+2*Glob_bvc(q,i)*Glob_bvc(p,i)*inv_tAkl(q,p)
          enddo
        enddo
        temp4=temp1/(sqrt(temp3)*temp3)
        temp5=ONE/temp3
        do k=1,NumDensGridPoints
          temp2=DensGrid(k)*DensGrid(k)
          Denskl(i,k)=temp4*exp(-temp2*temp5)
          !Denskl(i,k)=temp2*temp4*exp(-temp2*temp5) !Multiplied by \xi^2
        enddo
      enddo
    endif

    if (AreMCorrFuncNeeded) then
      temp1=MSkl/(Glob_Pi*Glob_SqrtPi)
      p=0
      do i=1,n
        do j=i,n
          p=p+1
          temp3=temp1/(sqrtMTrAJ(j,i)*MTrAJ(j,i))
          do k=1,NumCFGridPoints
            temp2=CFGrid(k)*CFGrid(k)
            CFkl(p,k)=temp3*exp(-temp2/MTrAJ(j,i))
            !CFkl(p,k)=temp2*temp3*exp(-temp2/TrAJ(j,i))  !Multiplied by \xi^2
          enddo
        enddo
      enddo
    end if

    if (AreMPartDensNeeded) then
      temp1=MSkl/(Glob_Pi*Glob_SqrtPi)
      do i=1,n+1
        temp3=ZERO
        do p=1,n
          temp3=temp3+Glob_bvc(p,i)*Glob_bvc(p,i)*inv_invtAkinvtAl(p,p)
          do q=p+1,n
            temp3=temp3+2*Glob_bvc(q,i)*Glob_bvc(p,i)*inv_invtAkinvtAl(q,p)
          enddo
        enddo
        temp3=temp3*4
        temp4=temp1/(sqrt(temp3)*temp3)
        do k=1,NumDensGridPoints
          temp2=DensGrid(k)*DensGrid(k)
          Denskl(i,k)=temp4*exp(-temp2/temp3)
          !Denskl(i,k)=temp2*temp4*exp(-temp2/temp3) !Multiplied by \xi^2
        enddo
      enddo
    end if
  end subroutine MatrixElementsAll_RG_0S

  function trace(k,M)
    real(wp) trace
    integer k
    integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles
    real(wp) M(nn,nn)
    integer i
    trace=ZERO
    do i=1,k
      trace=trace+M(i,i)
    enddo
  end function trace

  function ME_rXr_over_rij(X,i,j,inv_tAkl,ME_1_over_rij,TrAJ)
!Function ME_rXr_over_rij computes the following
!matrix element with real L=0 Gaussians phi_k and phi_l:
!<phi_k| r'Xr/r_{ij} |phi_l>
!Here X is an arbitrary (i.e. nonsymmetric) real matrix.
!Index i can be equal to j. In the latter case
!<phi_k| r'Xr/r_{i} |phi_l> is computed
!Input:
!            X  :: n x n real matrix
!           i,j :: indices denoting i and j.
!        inv_tAkl :: n x n real matrix where the inverse of Ak+tAl is stored
! ME_1_over_rij :: the value of <phi_k| 1/r_{ij} |phi_l> matrix element
!          TrAJ :: the value of Tr[inv_tAkl Jij]
!Note that n=Glob_n and nn=Glob_AllowedNumOfPseudoParticles. Although
!all arrays (both arguments and local ones) are static and have dimension
!nn x nn, only n x n subarrays are referenced.
    real(wp) ME_rXr_over_rij
    integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles

!Arguments:
    real(wp)  X(nn,nn),inv_tAkl(nn,nn),ME_1_over_rij,TrAJ
    integer      i,j
!Local variables
    real(wp)  AXi(nn),AXj(nn)
    real(wp)  TrAX,TrAXAJ
    integer      n,m,p,q

    n=Glob_n

    TrAX=ZERO
    do m=1,n
      do p=1,n
        TrAX=TrAX+inv_tAkl(m,p)*X(p,m)
      enddo
    enddo

    do m=1,n
      AXi(m)=ZERO
      do p=1,n
        AXi(m)=AXi(m)+inv_tAkl(i,p)*X(p,m)
      enddo
    enddo
    if (j/=i) then
      do m=1,n
        AXj(m)=ZERO
        do p=1,n
          AXj(m)=AXj(m)+inv_tAkl(j,p)*X(p,m)
        enddo
      enddo
    endif

    if (i==j) then
      TrAXAJ=ZERO
      do m=1,n
        TrAXAJ=TrAXAJ+AXi(m)*inv_tAkl(m,i)
      enddo
    else
      TrAXAJ=ZERO
      do m=1,n
        TrAXAJ=TrAXAJ+(AXi(m)-AXj(m))*(inv_tAkl(m,i)-inv_tAkl(m,j))
      enddo
    endif

    ME_rXr_over_rij=ME_1_over_rij*(3*TrAX-TrAXAJ/TrAJ)/2

  end function ME_rXr_over_rij

  subroutine ME_rXr_over_rij_all(X,inv_tAkl,rmkl,TrAJ,ME)
!Subroutine ME_rXr_over_rij_all computes the following
!matrix elements with real L=0 Gaussians phi_k and phi_l:
!<phi_k| r'Xr/r_{ij} |phi_l>    for all combinations i,j=1..n at the same time
!Index i can be equal to j. In the latter case
!<phi_k| r'Xr/r_{i} |phi_l> is computed
!Input:
!            X  :: n x n real matrix
!      inv_tAkl :: n x n real matrix where the inverse of Ak+tAl is stored
!          rmkl :: the values of <phi_k| 1/r_{ij} |phi_l> matrix element
!          TrAJ :: the values of Tr[inv_tAkl Jij]
!Output:
!            ME :: n x n real matrix where all computed matrix elements are returned
!Note that n=Glob_n and nn=Glob_AllowedNumOfPseudoParticles. Although
!all arrays (both arguments and local ones) are static and have dimension
!nn x nn, only n x n subarrays are referenced.
    integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles

!Arguments:
    real(wp)  X(nn,nn),inv_tAkl(nn,nn),rmkl(Glob_n,Glob_n),TrAJ(nn,nn),ME(nn,nn)
!Local variables
    integer      i,j
    real(wp)  AXi(nn),AXj(nn)
    real(wp)  TrAX,TrAXAJ
    integer      n,m,p,q

    n=Glob_n

    TrAX=ZERO
    do m=1,n
      do p=1,n
        TrAX=TrAX+inv_tAkl(m,p)*X(p,m)
      enddo
    enddo

    do i=1,n
      do j=i,n
        do m=1,n
          AXi(m)=ZERO
          do p=1,n
            AXi(m)=AXi(m)+inv_tAkl(i,p)*X(p,m)
          enddo
        enddo
        if (j/=i) then
          do m=1,n
            AXj(m)=ZERO
            do p=1,n
              AXj(m)=AXj(m)+inv_tAkl(j,p)*X(p,m)
            enddo
          enddo
        endif

        if (i==j) then
          TrAXAJ=ZERO
          do m=1,n
            TrAXAJ=TrAXAJ+AXi(m)*inv_tAkl(m,i)
          enddo
        else
          TrAXAJ=ZERO
          do m=1,n
            TrAXAJ=TrAXAJ+(AXi(m)-AXj(m))*(inv_tAkl(m,i)-inv_tAkl(m,j))
          enddo
        endif

        ME(i,j)=rmkl(i,j)*(3*TrAX-TrAXAJ/TrAJ(i,j))/2
        ME(j,i)=ME(i,j)

      enddo
    enddo

  end subroutine ME_rXr_over_rij_all

  function ME_rXr_rYr_over_rij(X,Y,i,j,inv_tAkl,ME_1_over_rij,TrAJ)
!Function ME_rXr_rYr_over_rij computes the following
!matrix element with real L=0 Gaussians phi_k and phi_l:
!<phi_k| (r'Xr)(r'Yr)/r_{ij} |phi_l>
!Here X and Y are arbitrary (i.e. nonsymmetric) real matrices.
!Index i can be equal to j. In the latter case
!<phi_k| (r'Xr)(r'Yr)/r_{i} |phi_l> is computed
!Input:
!            X  :: n x n real matrix
!            Y  :: n x n real matrix
!           i,j :: indices denoting i and j.
!        inv_tAkl :: n x n real matrix where the inverse of Ak+tAl is stored
! ME_1_over_rij :: the value of <phi_k| 1/r_{ij} |phi_l> matrix element
!          TrAJ :: the value of Tr[inv_tAkl Jij]
!Note that n=Glob_n and nn=Glob_AllowedNumOfPseudoParticles. Although
!all arrays (both arguments and local ones) are static and have dimension
!nn x nn, only n x n subarrays are referenced.
    real(wp) ME_rXr_rYr_over_rij
    integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles

!Arguments:
    real(wp)  X(nn,nn),Y(nn,nn),inv_tAkl(nn,nn),ME_1_over_rij,TrAJ
    integer      i,j
!Local variables
    real(wp)  Ys(nn,nn),Xs(nn,nn),AX(nn,nn),AY(nn,nn),AXAYi(nn),AXAYj(nn)
    real(wp)  TrAX,TrAY,TrAXAY,TrAXAJ,TrAYAJ,TrAXAYAJ
    integer      n,m,p,q

    n=Glob_n

    do p=1,n
      Ys(p,p)=Y(p,p)
      do q=p+1,n
        Ys(p,q)=(Y(p,q)+Y(q,p))/2
        Ys(q,p)=Ys(p,q)
      enddo
    enddo
    do p=1,n
      Xs(p,p)=X(p,p)
      do q=p+1,n
        Xs(p,q)=(X(p,q)+X(q,p))/2
        Xs(q,p)=Xs(p,q)
      enddo
    enddo

    do p=1,n
      do q=1,n
        AY(p,q)=ZERO
        do m=1,n
          AY(p,q)=AY(p,q)+inv_tAkl(p,m)*Ys(m,q)
        enddo
      enddo
    enddo
    do p=1,n
      do q=1,n
        AX(p,q)=ZERO
        do m=1,n
          AX(p,q)=AX(p,q)+inv_tAkl(p,m)*Xs(m,q)
        enddo
      enddo
    enddo

    do m=1,n
      AXAYi(m)=ZERO
      do p=1,n
        AXAYi(m)=AXAYi(m)+AX(i,p)*AY(p,m)
      enddo
    enddo
    if (j/=i) then
      do m=1,n
        AXAYj(m)=ZERO
        do p=1,n
          AXAYj(m)=AXAYj(m)+AX(j,p)*AY(p,m)
        enddo
      enddo
    endif

    TrAY=ZERO
    TrAX=ZERO
    do m=1,n
      TrAY=TrAY+AY(m,m)
      TrAX=TrAX+AX(m,m)
    enddo
    TrAXAY=ZERO
    do m=1,n
      do p=1,n
        TrAXAY=TrAXAY+AX(m,p)*AY(p,m)
      enddo
    enddo

    if (i==j) then
      TrAXAJ=ZERO
      do m=1,n
        TrAXAJ=TrAXAJ+AX(i,m)*inv_tAkl(m,i)
      enddo
      TrAYAJ=ZERO
      do m=1,n
        TrAYAJ=TrAYAJ+AY(i,m)*inv_tAkl(m,i)
      enddo
      TrAXAYAJ=ZERO
      do m=1,n
        TrAXAYAJ=TrAXAYAJ+AXAYi(m)*inv_tAkl(m,i)
      enddo
    else
      TrAXAJ=ZERO
      do m=1,n
        TrAXAJ=TrAXAJ+(AX(i,m)-AX(j,m))*(inv_tAkl(m,i)-inv_tAkl(m,j))
      enddo
      TrAYAJ=ZERO
      do m=1,n
        TrAYAJ=TrAYAJ+(AY(i,m)-AY(j,m))*(inv_tAkl(m,i)-inv_tAkl(m,j))
      enddo
      TrAXAYAJ=ZERO
      do m=1,n
        TrAXAYAJ=TrAXAYAJ+(AXAYi(m)-AXAYj(m))*(inv_tAkl(m,i)-inv_tAkl(m,j))
      enddo
    endif

    ME_rXr_rYr_over_rij=ME_1_over_rij*(   &
                         9*TrAY*TrAX + 6*TrAXAY           &
                         - ( 3*(TrAY*TrAXAJ + TrAX*TrAYAJ) + 4*TrAXAYAJ )/TrAJ   &
                         + 3*TrAXAJ*TrAYAJ/(TrAJ*TrAJ)    &
                         )/4

  end function ME_rXr_rYr_over_rij

  function myME_dXd_dYd(X,Y,inv_tAkl,tAk,tAl,det_tAkl)

    real(wp) myME_dXd_dYd, det_tAkl
    integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles
    real(wp) X(nn,nn),Y(nn,nn),inv_tAkl(nn,nn),tAk(nn,nn),tAl(nn,nn)

!Local variables:
    integer :: i,j,k,n
    real(wp) :: temp, trAXs, trAYs, trXAl, trYAk, trAYsAXs, commonFactor
    real(wp) :: XAl(nn,nn), Xs(nn,nn), YAk(nn,nn), Ys(nn,nn),&
                   AXs(nn,nn), YsAXs(nn,nn)
    real(wp) :: term1, term2, term3, term4

    n=Glob_n
!!! Q-part  !!!
! Build Xs matrix
    XAl = ZERO
    do i=1,n
      do j=1,n
        temp = ZERO
        do k=1,n
          temp = temp + X(i,k)*tAl(k,j)
        enddo
        XAl(i,j) = temp
      enddo
    enddo

    Xs = ZERO
    do i=1,n
      do j=1,n
        temp = ZERO
        do k=1,n
          temp = temp + tAl(i,k) * XAl(k,j)
        enddo
        Xs(i,j) = temp
      enddo
    enddo

    YAk = ZERO
    do i=1,n
      do j=1,n
        temp = ZERO
        do k=1,n
          temp = temp + Y(i,k)*tAk(k,j)
        enddo
        YAk(i,j) = temp
      enddo
    enddo

    Ys = ZERO
    do i=1,n
      do j=1,n
        temp = ZERO
        do k=1,n
          temp = temp + tAk(i,k) * YAk(k,j)
        enddo
        Ys(i,j) = temp
      enddo
    enddo

!Symmetrize
    do i = 1,n
      do j = i+1,n
        temp=ONEHALF*(Xs(j,i)+Xs(i,j))
        Xs(j,i) = temp
        Xs(i,j) = temp
      enddo
    enddo

    do i = 1,n
      do j = i+1,n
        temp=ONEHALF*(Ys(j,i)+Ys(i,j))
        Ys(j,i) = temp
        Ys(i,j) = temp
      enddo
    enddo
!End Symmetrize

    AXs = ZERO
    do i=1,n
      do j=1,n
        temp = ZERO
        do k=1,n
          temp = temp + inv_tAkl(i,k)*Xs(k,j)
        enddo
        AXs(i,j) = temp
      enddo
    enddo

    YsAXs = ZERO
    do i=1,n
      do j=1,n
        temp = ZERO
        do k=1,n
          temp = temp + Ys(i,k)*AXs(k,j)
        enddo
        YsAXs(i,j) = temp
      enddo
    enddo

    trAXs = ZERO
    do i=1,n
      do j=1,n
        trAXs = trAXs + inv_tAkl(i,j)*Xs(j,i)
      enddo
    enddo

    trAYs = ZERO
    do i=1,n
      do j=1,n
        trAYs = trAYs + inv_tAkl(i,j)*Ys(j,i)
      enddo
    enddo

    trXAl = ZERO
    do i=1,n
      trXAl = trXAl + XAl(i,i)
    enddo

    trYAk = ZERO
    do i=1,n
      trYAk = trYAk + YAk(i,i)
    enddo

    trAYsAXs = ZERO
    do i=1,n
      do j=1,n
        trAYsAXs = trAYsAXs + inv_tAkl(i,j)*YsAXs(j,i)
      enddo
    enddo

    commonFactor = Glob_PiRaised3n2/(det_tAkl*sqrt(det_tAkl))

    term1 = 16._wp*commonFactor*(&
            NINE/FOUR*trAXs*trAYs+ THREE/TWO*trAYsAXs)
    term2 = -36._wp*commonFactor*trXAl*trAYs
    term3 = -36._wp*commonFactor*trYAk*trAXs
    term4 = 36._wp*commonFactor*trYAk*trXAl

    myME_dXd_dYd = term1 + term2 + term3 + term4

  end function myME_dXd_dYd

  function myME_dXd(X,inv_tAkl,tAl,det_tAkl)

    real(wp) myME_dXd, det_tAkl
    integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles
    real(wp) X(nn,nn),inv_tAkl(nn,nn),tAl(nn,nn)

!Local variables:
    integer :: i,j,k,n
    real(wp) :: temp, trXAl, trAXs, commonFactor
    real(wp) :: XAl(nn,nn), Xs(nn,nn)

    n=Glob_n
!!! Q-part  !!!
    !Build Xs matrix
    XAl = ZERO
    do i=1,n
      do j=1,n
        temp = ZERO
        do k=1,n
          temp = temp + X(i,k)*tAl(k,j)
        enddo
        XAl(i,j) = temp
      enddo
    enddo

    Xs = ZERO
    do i=1,n
      do j=1,n
        temp = ZERO
        do k=1,n
          temp = temp + tAl(i,k) * XAl(k,j)
        enddo
        Xs(i,j) = temp
      enddo
    enddo

!Symmetrize
    do i = 1,n
      do j = i+1,n
        temp=ONEHALF*(Xs(j,i)+Xs(i,j))
        Xs(j,i) = temp
        Xs(i,j) = temp
      enddo
    enddo

    trXAl = ZERO
    do i=1,n
      trXAl = trXAl + XAl(i,i)
    enddo

    trAXs = ZERO
    do i=1,n
      do j=1,n
        trAXs = trAXs + inv_tAkl(i,j)*Xs(j,i)
      enddo
    enddo

    commonFactor = SIX*Glob_PiRaised3n2/(det_tAkl*sqrt(det_tAkl))

    myME_dXd = commonFactor*(trAXs - trXAl)

  end function myME_dXd

  function myME_over_rij_dXd(X,p,q,inv_tAkl,tAl,det_tAkl)
!Function ME_1_over_rij_dXd computes the following
!matrix element with real L=0 Gaussians phi_k and phi_l:
!<phi_k| (1/r_{ij})(nabla_r'*X*nabla_r) |phi_l>
!Here X is arbitrary (i.e. nonsymmetric) real matrix.
!Index i can be equal to j. In the latter case
!<phi_k| (1/r_i)(nabla_r'*X*nabla_r) |phi_l>is computed
!Input:
!            X  :: n x n real matrix
!           i,j :: indices denoting i and j.
!      inv_tAkl :: n x n real matrix where the inverse of Ak+tAl is stored
! ME_1_over_rij :: the value of <phi_k| 1/r_{ij} |phi_l> matrix element
!          TrAJ :: the value of Tr[inv_tAkl Jij]
!Note that n=Glob_n and nn=Glob_AllowedNumOfPseudoParticles. Although
!all arrays (both arguments and local ones) are static and have dimension
!nn x nn, only n x n subarrays are referenced.

!Input parameters:
    real(wp) myME_over_rij_dXd, det_tAkl
    integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles
    real(wp) :: X(nn,nn),inv_tAkl(nn,nn),tAl(nn,nn)
    integer p,q

!Local variables:
    integer :: i,j,k,n
    real(wp) :: temp, gamma, trAXs, jijAXsAjij, trXAl, commonFactor
    real(wp) :: XAl(nn,nn), Xs(nn,nn), XsA(nn,nn), AXsA(nn,nn)

    n=Glob_n
!!! Q-part  !!!
    !Build Xs matrix
    XAl = ZERO
    do i=1,n
      do j=1,n
        temp = ZERO
        do k=1,n
          temp = temp + X(i,k)*tAl(k,j)
        enddo
        XAl(i,j) = temp
      enddo
    enddo

    Xs = ZERO
    do i=1,n
      do j=1,n
        temp = ZERO
        do k=1,n
          temp = temp + tAl(i,k) * XAl(k,j)
        enddo
        Xs(i,j) = temp
      enddo
    enddo

!Symmetrize
    do i = 1,n
      do j = i+1,n
        temp=ONEHALF*(Xs(j,i)+Xs(i,j))
        Xs(j,i) = temp
        Xs(i,j) = temp
      enddo
    enddo

    XsA = ZERO
    do i=1,n
      do j=1,n
        temp = ZERO
        do k=1,n
          temp = temp + Xs(i,k)*inv_tAkl(k,j)
        enddo
        XsA(i,j) = temp
      enddo
    enddo

    AXsA = ZERO
    do i=1,n
      do j=1,n
        temp = ZERO
        do k=1,n
          temp = temp + inv_tAkl(i,k)*XsA(k,j)
        enddo
        AXsA(i,j) = temp
      enddo
    enddo

    jijAXsAjij = ZERO

    if (p == q) then
      jijAXsAjij = AXsA(p,p)
      gamma = inv_tAkl(p,p)
      gamma = ONE/sqrt(gamma)
    else
      jijAXsAjij = AXsA(p,p) + AXsA(q,q) - AXsA(p,q) - AXsA(q,p)
      gamma = inv_tAkl(p,p) + inv_tAkl(q,q) - inv_tAkl(p,q) - inv_tAkl(q,p)
      gamma = ONE/sqrt(gamma)
    endif

    trAXs = ZERO
    do i=1,n
      do j=1,n
        trAXs = trAXs + inv_tAkl(i,j)*Xs(j,i)
      enddo
    enddo

    trXAl = ZERO
    do i=1,n
      trXAl = trXAl + XAl(i,i)
    enddo

    commonFactor = SIX*TWO*Glob_PiRaised3n2/(Glob_SqrtPi*det_tAkl*sqrt(det_tAkl))

    myME_over_rij_dXd = commonFactor*(&
                        gamma*(trAXs - trXAl) - ONE/THREE*gamma**3*jijAXsAjij)

  end function myME_over_rij_dXd

  function ME_1_over_rij_dXd(X,i,j,inv_tAkl,tAl,ME_1_over_rij,TrAJ)
!Function ME_1_over_rij_dXd computes the following
!matrix element with real L=0 Gaussians phi_k and phi_l:
!<phi_k| (1/r_{ij})(nabla_r'*X*nabla_r) |phi_l>
!Here X and Y are arbitrary (i.e. nonsymmetric) real matrices.
!Index i can be equal to j. In the latter case
!<phi_k| (1/r_i)(nabla_r'*X*nabla_r) |phi_l>is computed
!Input:
!            X  :: n x n real matrix
!            Y  :: n x n real matrix
!           i,j :: indices denoting i and j.
!      inv_tAkl :: n x n real matrix where the inverse of Ak+tAl is stored
! ME_1_over_rij :: the value of <phi_k| 1/r_{ij} |phi_l> matrix element
!          TrAJ :: the value of Tr[inv_tAkl Jij]
!Note that n=Glob_n and nn=Glob_AllowedNumOfPseudoParticles. Although
!all arrays (both arguments and local ones) are static and have dimension
!nn x nn, only n x n subarrays are referenced.
!Input parameters:
    real(wp) ME_1_over_rij_dXd
    integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles
    real(wp) X(nn,nn),inv_tAkl(nn,nn),tAl(nn,nn),ME_1_over_rij,TrAJ
!Local variables:
    integer i,j,p,q,k,n
    real(wp) M(nn,nn),Z(nn,nn)
    real(wp) tr1,t

    n=Glob_n
!Compute Z=tAl*X, tr1=trace[Z], and M=tAl*X*tAl
    tr1=ZERO
    do p=1,n
      do q=1,n
        t=ZERO
        do k=1,n
          t=t+tAl(q,k)*X(k,p)
        enddo
        Z(q,p)=t
      enddo
      tr1=tr1+Z(p,p)
    enddo
    do p=1,n
      do q=1,n
        t=ZERO
        do k=1,n
          t=t+Z(q,k)*tAl(k,p)
        enddo
        M(q,p)=t
      enddo
    enddo
!Compute the matrix element
    ME_1_over_rij_dXd=4*ME_rXr_over_rij(M,i,j,inv_tAkl,ME_1_over_rij,TrAJ)-6*tr1*ME_1_over_rij
  end function ME_1_over_rij_dXd

  subroutine ME_1_over_rij_dXd_all(X,inv_tAkl,tAl,rmkl,TrAJ,ME)
!Subroutine ME_1_over_rij_dXd computes
!matrix elements with real L=0 Gaussians phi_k and phi_l:
!<phi_k| (1/r_{ij})(nabla_r'*X*nabla_r) |phi_l>
!for all combinations of i and j indexes.
!Here X and Y are arbitrary (i.e. nonsymmetric) real matrices.
!When i is equal to j it stands for the matrix elements
!<phi_k| (1/r_i)(nabla_r'*X*nabla_r) |phi_l>
!Input:
!            X  :: n x n real matrix
!            Y  :: n x n real matrix
!           i,j :: indices denoting i and j.
!      inv_tAkl :: n x n real matrix where the inverse of Ak+tAl is stored
!          rmkl :: the values of <phi_k| 1/r_{ij} |phi_l> matrix element
!          TrAJ :: the values of Tr[inv_tAkl Jij]
!Output:
!            ME :: n x n real matrix where all computed matrix elements are returned
!Note that n=Glob_n and nn=Glob_AllowedNumOfPseudoParticles. Although
!all arrays (both arguments and local ones) are static and have dimension
!nn x nn, only n x n subarrays are referenced.
!Input parameters:
    integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles
    real(wp) X(nn,nn),inv_tAkl(nn,nn),tAl(nn,nn),rmkl(Glob_n,Glob_n),TrAJ(nn,nn)
    real(wp) ME(nn,nn)
!Local variables:
    integer i,j,p,q,k,n
    real(wp) M(nn,nn),Z(nn,nn)
    real(wp) tr1,t

    n=Glob_n
!Compute Z=tAl*X, tr1=trace[Z], and M=tAl*X*tAl
    tr1=ZERO
    do p=1,n
      do q=1,n
        t=ZERO
        do k=1,n
          t=t+tAl(q,k)*X(k,p)
        enddo
        Z(q,p)=t
      enddo
      tr1=tr1+Z(p,p)
    enddo
    do p=1,n
      do q=1,n
        t=ZERO
        do k=1,n
          t=t+Z(q,k)*tAl(k,p)
        enddo
        M(q,p)=t
      enddo
    enddo
!Compute all matrix elements <phi_k| r'Mr/r_{ij} |phi_l>
    call ME_rXr_over_rij_all(M,inv_tAkl,rmkl,TrAJ,Z)
    t=6*tr1
    do i=1,n
      do j=i,n
        ME(j,i)=4*Z(j,i)-t*rmkl(j,i)
        ME(i,j)=ME(j,i)
      enddo
    enddo

  end subroutine ME_1_over_rij_dXd_all

!!Old, slow, and simple version
!function ME_1_over_rij_dXd(X,i,j,inv_tAkl,tAl,ME_1_over_rij,TrCJ)
!real(wp) ME_1_over_rij_dXd
!integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles
!real(wp) X(nn,nn),inv_tAkl(nn,nn),tAl(nn,nn),ME_1_over_rij,TrCJ
!integer i,j,n
!real(wp) M(nn,nn),Z(nn,nn)
!real(wp) tr1
!n=Glob_n
!M(1:n,1:n)=matmul(tAl(1:n,1:n),matmul(X(1:n,1:n),tAl(1:n,1:n)))
!Z(1:n,1:n)=matmul(tAl(1:n,1:n),X(1:n,1:n))
!tr1=trace(n,Z)
!ME_1_over_rij_dXd=4*ME_rXr_over_rij(M,i,j,inv_tAkl,ME_1_over_rij,TrCJ)-6*tr1*ME_1_over_rij
!end function ME_1_over_rij_dXd

  subroutine spinPreCalc(n, nFactorial, parityFactor, SSFmassChargeCoefficient, &
                         AnihMassChargeCoefficient, ketMatrix, spatialYoung, &
                         positronPosition, numberOfSpinFunctions, spinFreeME, SiSjME)
    use spinStuff
    implicit none

    character(len = maxLen), intent(in) :: spatialYoung
    integer, intent(in) :: n, nFactorial

    real(wp), dimension(nFactorial), intent(out) :: parityFactor
    real(wp), dimension(n, n, nFactorial), intent(out) :: ketMatrix
    real(wp), dimension(n, n), intent(out) :: SSFmassChargeCoefficient, AnihMassChargeCoefficient
    integer, intent(out) :: positronPosition, numberOfSpinFunctions

    real(wp), dimension(nFactorial), intent(out) :: spinFreeME
    real(kind = wp), dimension(n, n, 2, nFactorial), intent(out) :: SiSjME

    ! local variables
    integer :: i, j, k, l, m
    character(len = maxLen) :: mySpatialYoung
    integer, dimension(nFactorial) :: parities
    integer, dimension(n, n, nFactorial) :: allPermutations

    SSFmassChargeCoefficient = ZERO
    do i = 1, n
      do j = 1, n
        SSFmassChargeCoefficient(i, j) = -Glob_PseudoCharge(i) * Glob_PseudoCharge(j) / &
                                         (Glob_Mass(i + 1) * Glob_Mass(j + 1)) * EIGHT * Glob_Pi / THREE
      enddo
    enddo

    AnihMassChargeCoefficient = ZERO
    do i = 1, n
      do j = 1, n
        AnihMassChargeCoefficient(i, j) = -Glob_PseudoCharge(i) * Glob_PseudoCharge(j) / &
                                          (Glob_Mass(i + 1) * Glob_Mass(j + 1)) * TWO * Glob_Pi
      enddo
    enddo

    ! now we deal with the spin stuff

    ! rename the particles
    mySpatialYoung = spatialYoung
    do i = 1, maxLen
      if (mySpatialYoung(i:i) == 'P') then
        read(mySpatialYoung(i + 1:i + 1), *) k
        read(mySpatialYoung(i + 2:i + 2), *) j
        write(mySpatialYoung(i + 1:i + 1), '(i1)') k - 1
        write(mySpatialYoung(i + 2:i + 2), '(i1)') j - 1
      endif
    enddo

    call getSpinOperatorsMeanValues(n, nFactorial, mySpatialYoung, positronPosition, numberOfSpinFunctions, &
                                    allPermutations, parities, spinFreeME, SiSjMe)

    ketMatrix = ZERO
    do i = 1, nFactorial

      do k = 1, n
        do l = 1, n

          ketMatrix(k, l, i) = real(allPermutations(l, k, i))
          ! note the transposition here

        enddo
      enddo

    enddo

    do i = 1, nFactorial
      parityFactor(i) = real(parities(i))
    enddo

  end subroutine

  subroutine OverlapMatrixElement_RG_0S(vechLk, P, Skk)
!This subroutine computes symmetry adapted matrix element with
!two real L=0 correlated Gaussians:
!
!fk =  exp[-r'(Lk*Lk')r]
!
!Symmetry adaption is applied to the ket using
!permutation matrices Glob_YHYMatr(:,:,1:Glob_NumYHYTerms)
!
!Input:
!   vechLk :: Array of length (n(n+1)/2) of exponential parameters.
!   P   :: The symmetry permutation matrix of size n x n
!Output:
!   Skk         ::        Overlap matrix element (normalized)

!Arguments
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
    real(wp)       Lk(nn,nn), Ll(nn,nn)
    real(wp)       Ak(nn,nn), tAl(nn,nn), tAkl(nn,nn)
    real(wp)       inv_tAkl(nn,nn)
    real(wp)       W1(nn,nn)
    real(wp)       temp1
    real(wp)       det_tAkl
    integer           i, j, k, indx

    n=Glob_n
    np=Glob_np
!First we build matrices Lk, Ll, Ak, Al from vechLk, vechLl.
    indx=0
    do i=1,n
    do j=i,n
      indx=indx+1
      Lk(i,j)=ZERO
      Lk(j,i)=vechLk(indx)
      Ll(i,j)=Lk(i,j)
      Ll(j,i)=Lk(j,i)
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

!Evaluating overlap

!temp1=abs(det_Ll*det_Lk)/det_tAkl
!Skl=Glob_2Raised3n2*temp1*sqrt(temp1)
    Skk=Glob_PiRaised3n2/(det_tAkl*sqrt(det_tAkl))  !new line

  end subroutine OverlapMatrixElement_RG_0S

end module matelem
