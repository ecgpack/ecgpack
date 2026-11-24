module matelem
!Module matelem contains subroutines for computing
!matrix elements with complex L=0 Gaussians without normalization.
  use globvars
  use misc
  implicit none

contains

#ifdef USE_CUDA
  attributes(host,device) &
#endif
  subroutine MatrixElementsHS_CG_0S(Lk, Ll, Ck, Cl, P, &
                            Hkl, Skl, Dk, Dl, grad_k, grad_l)
!This subroutine computes symmetry adapted matrix element with
!two complex L=0 correlated Gaussians:
!
!fk =  exp[-r'(Lk*Lk'+i*Bk)r] = exp[-r'Ck r]
!
!Function fk is the bra (it enters the integrals complex conjugated)
!and fl is the ket. Symmetry adaption is applied to the ket using
!permutation matrices Glob_YHYMatr(:,:,1:Glob_NumYHYTerms)
!
!Input:
!   Lk, Ll :: The lower-triangular Cholesky-style parameter matrices (the
!             unpacked vechLk/vechLl). They depend only on one basis
!             function, so the caller precomputes them once per function per
!             sweep (PrecomputeMatrices_L_C below) instead of this routine
!             unpacking them for every (pair x term) call.
!   Ck, Cl :: Ck=Lk*Lk'+i*Bk, Cl=Ll*Ll'+i*Bl (complex symmetric matrices)
!             -- precomputed for the same reason.
!   P   :: The symmetry permutation matrix of size n x n
!   grad_k, grad_l :: Gradient flags
!   grad_k=.true.  means that dHkldvechLk, dHkldvechBk, dSkldvechLk,
!                  dSkldvechBk need to be computed.
!   grad_l=.true.  means that dHkldvechLl, dHkldvechBl, dSkldvechLl,
!                  dSkldvechBl need to be computed.
!Output:
!   Hkl         ::        Hamiltonian term (not normalized)
!   Skl         ::        Overlap matrix element (not normalized)
!   Dk,Dl:: derivatives of Hkl and Skl wrt the nonlinear parameters
!           of the bra (function k) and of the ket (function l)
!           respectively. They are ordered in the following manner:
!           Dk=(dHkldvechLk,dHkldvechBk,dSkldvechLk,dSkldvechBk)
!           Dl=(dHkldvechLl,dHkldvechBl,dSkldvechLl,dSkldvechBl)
!
!The formulas are those of S. Bubin and L. Adamowicz, J. Chem. Phys.
!124, 224317 (2006), written for Gaussians that are not normalized.
!With the notations
!   tCl  = P'*Cl*P,   tCkl = conjg(Ck)+tCl   (complex symmetric)
!they read
!   Skl = pi^(3n/2)/det(tCkl)^(3/2)
!   Tkl = 6*Skl*tr[M*tCl*inv_tCkl*conjg(Ck)]
!   Vkl = (2/sqrt(pi))*Skl*sum_ij q_ij/sqrt(tr[inv_tCkl*Jij])
!
!The structure of the routine and the scheme used for the gradient
!are the same as in the routine for real L=0 Gaussians,
!MatrixElementsHS_RG_0S (code RG_0S). In particular, the sum over
!particle pairs in the gradient of the potential energy is collapsed
!analytically into
!   B = inv_tCkl*C*inv_tCkl,  C = sum_ij c_ij*Jij,
!where C is assembled in O(n^2), which makes the whole gradient part
!O(n^3). Details are explained in the comments in the body.
!All complex matrices that appear here are symmetric (not hermitian),
!so the algebra is identical to that of the real case.

!Arguments
    integer,parameter :: nn=Glob_n
    real(wp),intent(in)      :: Lk(nn,nn), Ll(nn,nn)
    complex(wp),intent(in)   :: Ck(nn,nn), Cl(nn,nn)
    real(wp),intent(in)      :: P(nn,nn)
    complex(wp),intent(out)  :: Skl,Hkl
    complex(wp),intent(out)  :: Dk(4*Glob_np),Dl(4*Glob_np)
    logical,intent(in)          :: grad_k, grad_l

!Local variables
    integer           n, np, np2, np3
    real(wp)       PT(nn,nn)
    complex(wp)    tCl(nn,nn), tCkl(nn,nn)
    complex(wp)    inv_tCkl(nn,nn), inv_ttCkl(nn,nn)
    complex(wp)    inv_tCkltClM(nn,nn)
    complex(wp)    tr_inv_tCklJij32(nn,nn)
    complex(wp)    F(nn,nn), G(nn,nn)
    complex(wp)    Cmat(nn,nn), Bmat(nn,nn), Z(nn,nn)
    complex(wp)    W1(nn,nn), W2(nn,nn)
    complex(wp)    temp1, temp2, temp3, temp4, temp5
    complex(wp)    sqrtdet_tCkl
    complex(wp)    Tkl, Vkl, cV, HklOverSkl
    integer           i, j, k, indx

    n=Glob_n
    np=Glob_np
    np2=2*np
    np3=3*np
!Lk, Ll, Ck, Cl arrive precomputed (hoisted to once per function per sweep).

!Then we permute elements of Cl to account for
!the action of the permutation matrix
!tCl=P'*Cl*P
!We also form matrix tCkl=conjg(Ck)+tCl
    do i=1,nn
      do j=1,nn
        temp1=ZERO
        do k=1,nn
          temp1=temp1+P(k,j)*Cl(k,i)
        enddo
        W1(j,i)=temp1
      enddo
    enddo
    do i=1,nn
      do j=i,nn
        temp1=ZERO
        do k=1,nn
          temp1=temp1+W1(i,k)*P(k,j)
        enddo
        tCl(i,j)=temp1
        tCl(j,i)=temp1
        tCkl(i,j)=conjg(Ck(i,j))+temp1
        tCkl(j,i)=tCkl(i,j)
      enddo
    enddo

!After this we can do pseudo-Cholesky factorization of tCkl
!(pseudo because tCkl is complex symmetric, so it is factorized as
!a product of a lower triangular matrix and its transpose, not as a
!product of a lower triangular matrix and its hermitian conjugate).
!The pseudo-Cholesky factor will be temporarily stored in the
!lower triangle of W1.
!
!The real part of tCkl is positive definite. For such a matrix the
!real part of each pivot (the quantity whose square root is taken
!below) is positive, so that the principal values of the square
!roots are the proper ones. The product of the diagonal elements of
!the factor is then the properly defined square root of det(tCkl)
!(the one obtained by the analytic continuation from the real case,
!as required by the Gaussian integral formula). Evaluating the
!square root of the product of the pivots instead could yield a
!wrong sign when the phase of det(tCkl) exceeds pi by magnitude,
!which becomes possible when n>2.
    sqrtdet_tCkl=ONE
    do i=1,nn
      do j=i,nn
        temp1=tCkl(i,j)
        do k=i-1,1,-1
          temp1=temp1-W1(i,k)*W1(j,k)
        enddo
        if (i==j) then
          W1(i,i)=sqrt(temp1)
          sqrtdet_tCkl=sqrtdet_tCkl*W1(i,i)
        else
          W1(j,i)=temp1/W1(i,i)
          W1(i,j)=ZERO
        endif
      enddo
    enddo

!Inverting tCkl using its pseudo-Cholesky factor (stored in W1)
!and placing the result into inv_tCkl
    do i=1,nn
      W1(i,i)=ONE/W1(i,i)
      do j=i+1,nn
        temp1=ZERO
        do k=i,j-1
          temp1=temp1-W1(j,k)*W1(k,i)
        enddo
        W1(j,i)=temp1/W1(j,j)
      enddo
    enddo

    do i=1,nn
      do j=i,nn
        temp1=ZERO
        do k=j,nn
          temp1=temp1+W1(k,i)*W1(k,j)
        enddo
        inv_tCkl(i,j)=temp1
        inv_tCkl(j,i)=temp1
      enddo
    enddo

!Evaluating overlap
    Skl=Glob_PiRaised3n2/(sqrtdet_tCkl*sqrtdet_tCkl*sqrtdet_tCkl)

!Doing multiplication W2=inv_tCkl*tCl
    do i=1,nn
      do j=1,nn
        temp1=ZERO
        do k=1,nn
          temp1=temp1+inv_tCkl(j,k)*tCl(k,i)
        enddo
        W2(j,i)=temp1
      enddo
    enddo

!Doing multiplication inv_tCkltClM=inv_tCkl*tCl*M=W2*M
    do i=1,nn
      do j=1,nn
        temp1=ZERO
        do k=1,nn
          temp1=temp1+W2(j,k)*Glob_MassMatrix(k,i)
        enddo
        inv_tCkltClM(j,i)=temp1
      enddo
    enddo

!Computing kinetic energy, Tkl=tr[inv_tCkltClM*conjg(Ck)]
    Tkl=ZERO
    do i=1,nn
      temp1=ZERO
      do k=1,nn
        temp1=temp1+inv_tCkltClM(i,k)*conjg(Ck(k,i))
      enddo
      Tkl=Tkl+temp1
    enddo
    Tkl=SIX*Skl*Tkl

!Evaluating potential energy, Vkl, and tr[inv_tCkl*Jij]^(-3/2)
!The lower triangle of array tr_inv_tCklJij32
!will contain the corresponding quantities. The latter are needed
!only for the gradients, so in the gradientless case a leaner loop
!(one division per particle pair less) is used.
!The real part of tr[inv_tCkl*Jij] is always positive, so the
!principal value of the square root is the proper one.
    temp1=(TWO/Glob_SqrtPi)*Skl
    Vkl=ZERO
    if (grad_k.or.grad_l) then
      do i=1,nn
        temp3=inv_tCkl(i,i)
        temp4=sqrt(temp3)
        tr_inv_tCklJij32(i,i)=1/(temp4*temp3)
        temp5=temp1/temp4
        Vkl=Vkl+Glob_ScaledPseudoChargeMatrix(i,0)*temp5
      enddo
      do i=1,nn
        do j=i+1,nn
          temp3=inv_tCkl(i,i)+inv_tCkl(j,j)-inv_tCkl(j,i)-inv_tCkl(j,i)
          temp4=sqrt(temp3)
          tr_inv_tCklJij32(j,i)=1/(temp4*temp3)
          temp5=temp1/temp4
          Vkl=Vkl+Glob_ScaledPseudoChargeMatrix(i,j)*temp5
        enddo
      enddo
    else
      do i=1,nn
        temp5=temp1/sqrt(inv_tCkl(i,i))
        Vkl=Vkl+Glob_ScaledPseudoChargeMatrix(i,0)*temp5
      enddo
      do i=1,nn
        do j=i+1,nn
          temp3=inv_tCkl(i,i)+inv_tCkl(j,j)-inv_tCkl(j,i)-inv_tCkl(j,i)
          temp5=temp1/sqrt(temp3)
          Vkl=Vkl+Glob_ScaledPseudoChargeMatrix(i,j)*temp5
        enddo
      enddo
    endif

    Hkl=Tkl+Vkl

!In the gradientless case we are done
    if (.not.(grad_k.or.grad_l)) return

!Now we compute the gradients. In what follows
!  W2 = inv_tCkl*tCl  and  K = inv_tCkltClM = W2*M
!(both are already available from the computation of Tkl).
!
!Gradient of Skl (lower triangles, in vech order):
!  dSkl/dvechLk = -3*Skl*inv_tCkl*Lk
!  dSkl/dvechLl = -3*Skl*(P*inv_tCkl*P')*Ll
!  dSkl/dvechBk =  (3/2)*i*Skl*(2*inv_tCkl-diag(inv_tCkl))
!  dSkl/dvechBl = -(3/2)*i*Skl*(2*P*inv_tCkl*P'-diag(P*inv_tCkl*P'))
!They are stored directly in the tails Dk(2*np+1:4*np), Dl(2*np+1:4*np)
!of the output arrays.
!
!Gradient of Tkl. It needs the two matrices
!  F = inv_tCkl*tCl*M*tCl*inv_tCkl = K*W2'
!  G = P*(inv_tCkl*conjg(Ck)*M*conjg(Ck)*inv_tCkl)*P'
!Because inv_tCkl*conjg(Ck) = I-W2 (as conjg(Ck)+tCl=tCkl), the matrix
!in the parentheses of G equals (I-W2)*M*(I-W2') = M - K - K' + F, i.e.
!it costs only O(n^2) once F is known.
!
!Gradient of Vkl. A straightforward implementation sums, over all
!particle pairs, the matrices tQ=inv_tCkl*Jij*inv_tCkl (each multiplied
!by Lk or by P...P'*Ll) with coefficients
!c_ij=q_ij*tr[inv_tCkl*Jij]^(-3/2), which is O(n^5) overall. Since
!sum_ij c_ij*Jij = C is an n x n matrix assembled in O(n^2), the whole
!sum is simply
!  B = inv_tCkl*C*inv_tCkl
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
!The derivatives with respect to the elements of Bk and Bl are
!obtained from the same matrices at no additional cost:
!  dHkl/dvechBk = (Hkl/Skl)*dSkl/dvechBk - (i/2)*vech[2*Z-diag(Z)]
!  dHkl/dvechBl = (Hkl/Skl)*dSkl/dvechBl
!                 + (i/2)*vech[2*P*Z*P'-diag(P*Z*P')]

    cV=(TWO/Glob_SqrtPi)*Skl
    HklOverSkl=Hkl/Skl

    if (grad_k) then
      !dSkldvechLk: lower triangle of -3*Skl*inv_tCkl*Lk. The
      !symmetry of inv_tCkl and the lower triangularity of Lk make
      !the inner loop a dot product of two contiguous columns.
      temp2=-THREE*Skl
      indx=0
      do i=1,nn
        do j=i,nn
          indx=indx+1
          temp1=ZERO
          do k=i,nn
            temp1=temp1+inv_tCkl(k,j)*Lk(k,i)
          enddo
          Dk(np2+indx)=temp2*temp1
        enddo
      enddo
      !dSkldvechBk: (3/2)*i*Skl*vech[2*inv_tCkl-diag(inv_tCkl)]
      temp2=cmplx(ZERO,THREEHALF,wp)*Skl
      temp3=temp2+temp2
      indx=0
      do i=1,nn
        indx=indx+1
        Dk(np3+indx)=temp2*inv_tCkl(i,i)
        do j=i+1,nn
          indx=indx+1
          Dk(np3+indx)=temp3*inv_tCkl(j,i)
        enddo
      enddo
    endif

    if (grad_l) then
      !PT=P' is stored explicitly so that all the products with P
      !and P' below can be done with contiguous column access
      do j=1,nn
        do i=1,nn
          PT(i,j)=P(j,i)
        enddo
      enddo
      !calculating inv_ttCkl=P*inv_tCkl*P'
      !W1=inv_tCkl*PT
      do j=1,nn
        do i=1,nn
          temp1=ZERO
          do k=1,nn
            temp1=temp1+inv_tCkl(k,i)*PT(k,j)
          enddo
          W1(i,j)=temp1
        enddo
      enddo
      !inv_ttCkl=PT'*W1 (only the upper triangle, then mirrored)
      do j=1,nn
        do i=1,j
          temp1=ZERO
          do k=1,nn
            temp1=temp1+PT(k,i)*W1(k,j)
          enddo
          inv_ttCkl(i,j)=temp1
          inv_ttCkl(j,i)=temp1
        enddo
      enddo
      !dSkldvechLl: lower triangle of -3*Skl*inv_ttCkl*Ll
      temp2=-THREE*Skl
      indx=0
      do i=1,nn
        do j=i,nn
          indx=indx+1
          temp1=ZERO
          do k=i,nn
            temp1=temp1+inv_ttCkl(k,j)*Ll(k,i)
          enddo
          Dl(np2+indx)=temp2*temp1
        enddo
      enddo
      !dSkldvechBl: -(3/2)*i*Skl*vech[2*inv_ttCkl-diag(inv_ttCkl)]
      temp2=cmplx(ZERO,-THREEHALF,wp)*Skl
      temp3=temp2+temp2
      indx=0
      do i=1,nn
        indx=indx+1
        Dl(np3+indx)=temp2*inv_ttCkl(i,i)
        do j=i+1,nn
          indx=indx+1
          Dl(np3+indx)=temp3*inv_ttCkl(j,i)
        enddo
      enddo
    endif

!F=K*W2' (only the upper triangle is computed, then mirrored)
    do j=1,nn
      do i=1,j
        temp1=ZERO
        do k=1,nn
          temp1=temp1+inv_tCkltClM(i,k)*W2(j,k)
        enddo
        F(i,j)=temp1
        F(j,i)=temp1
      enddo
    enddo

!Assembling C=sum_ij c_ij*Jij, c_ij=q_ij*tr[inv_tCkl*Jij]^(-3/2)
!(Jii is a matrix whose only nonzero element is (i,i)=1; Jij, i/=j,
!has elements (i,i)=(j,j)=1, (i,j)=(j,i)=-1)
    do j=1,nn
      do i=1,nn
        Cmat(i,j)=ZERO
      enddo
    enddo
    do i=1,nn
      Cmat(i,i)=Glob_ScaledPseudoChargeMatrix(0,i)*tr_inv_tCklJij32(i,i)
    enddo
    do i=1,nn
      do j=i+1,nn
        temp1=Glob_ScaledPseudoChargeMatrix(i,j)*tr_inv_tCklJij32(j,i)
        Cmat(i,i)=Cmat(i,i)+temp1
        Cmat(j,j)=Cmat(j,j)+temp1
        Cmat(j,i)=Cmat(j,i)-temp1
        Cmat(i,j)=Cmat(i,j)-temp1
      enddo
    enddo

    do j=1,nn
      do i=1,nn
        temp1=ZERO
        do k=1,nn
          temp1=temp1+Cmat(k,i)*inv_tCkl(k,j)
        enddo
        W1(i,j)=temp1
      enddo
    enddo
    do j=1,nn
      do i=1,j
        temp1=ZERO
        do k=1,nn
          temp1=temp1+W1(k,i)*inv_tCkl(k,j)
        enddo
        Bmat(i,j)=temp1
        Bmat(j,i)=temp1
      enddo
    enddo

    if (grad_k) then
      temp2=12*Skl
      do j=1,nn
        do i=1,nn
          Z(i,j)=temp2*F(i,j)+cV*Bmat(i,j)
        enddo
      enddo
      !dHkldvechLk
      indx=0
      do i=1,nn
        do j=i,nn
          indx=indx+1
          temp1=ZERO
          do k=i,nn
            temp1=temp1+Z(k,j)*Lk(k,i)
          enddo
          Dk(indx)=HklOverSkl*Dk(np2+indx)+temp1
        enddo
      enddo
      !dHkldvechBk
      temp2=cmplx(ZERO,-ONEHALF,wp)
      temp3=cmplx(ZERO,-ONE,wp)
      indx=0
      do i=1,nn
        indx=indx+1
        Dk(np+indx)=HklOverSkl*Dk(np3+indx)+temp2*Z(i,i)
        do j=i+1,nn
          indx=indx+1
          Dk(np+indx)=HklOverSkl*Dk(np3+indx)+temp3*Z(j,i)
        enddo
      enddo
    endif

    if (grad_l) then
      temp2=12*Skl
      do j=1,nn
        do i=1,nn
          Z(i,j)=temp2*(Glob_MassMatrix(i,j)-inv_tCkltClM(i,j) &
                        -inv_tCkltClM(j,i)+F(i,j))+cV*Bmat(i,j)
        enddo
      enddo
      do j=1,nn
        do i=1,nn
          temp1=ZERO
          do k=1,nn
            temp1=temp1+Z(k,i)*PT(k,j)
          enddo
          W2(i,j)=temp1
        enddo
      enddo
      do j=1,nn
        do i=1,j
          temp1=ZERO
          do k=1,nn
            temp1=temp1+PT(k,i)*W2(k,j)
          enddo
          G(i,j)=temp1
          G(j,i)=temp1
        enddo
      enddo
      !dHkldvechLl
      indx=0
      do i=1,nn
        do j=i,nn
          indx=indx+1
          temp1=ZERO
          do k=i,nn
            temp1=temp1+G(k,j)*Ll(k,i)
          enddo
          Dl(indx)=HklOverSkl*Dl(np2+indx)+temp1
        enddo
      enddo
      !dHkldvechBl
      temp2=cmplx(ZERO,ONEHALF,wp)
      temp3=cmplx(ZERO,ONE,wp)
      indx=0
      do i=1,nn
        indx=indx+1
        Dl(np+indx)=HklOverSkl*Dl(np3+indx)+temp2*G(i,i)
        do j=i+1,nn
          indx=indx+1
          Dl(np+indx)=HklOverSkl*Dl(np3+indx)+temp3*G(j,i)
        enddo
      enddo
    endif

  end subroutine MatrixElementsHS_CG_0S

  subroutine PrecomputeMatrices_L_C(npt, Nmax, NonlinParam, Lk, Ck)
!This subroutine precomputes the Lk and Ck=Lk*Lk'+i*Bk matrices for all basis functions, so that
!the matrix elements can be computed faster without having to recompute these matrices each time.
!Input:
!  npt: number of nonlinear parameters per basis function (npt=2*n(n+1)/2)
!  Nmax: number of basis functions
!  NonlinParam: nonlinear parameters for all basis functions (npt x Nmax).
!               For each function they are ordered as (vechLk,vechBk)
!Output:
!  Lk: Lk matrices for all basis functions (nn x nn x Nmax)
!  Ck: Ck matrices for all basis functions (nn x nn x Nmax)
    integer,parameter     :: nn=Glob_n
    integer, intent(in)   :: npt, Nmax
    real(wp),intent(in)   :: NonlinParam(npt,Nmax)
    real(wp),intent(out)  :: Lk(nn,nn,Nmax)
    complex(wp),intent(out) :: Ck(nn,nn,Nmax)
    integer  :: f,i,j,k,indx,np
    real(wp) :: temp1
    np=npt/2
    do f=1,Nmax
      indx=0
      do i=1,nn
        do j=i,nn
          indx=indx+1
          Lk(i,j,f)=ZERO
          Lk(j,i,f)=NonlinParam(indx,f)
        enddo
      enddo
      indx=0
      do i=1,nn
        do j=i,nn
          indx=indx+1
          temp1=ZERO
          do k=1,i
            temp1=temp1+Lk(i,k,f)*Lk(j,k,f)
          enddo
          Ck(i,j,f)=cmplx(temp1,NonlinParam(np+indx,f),wp)
          Ck(j,i,f)=Ck(i,j,f)
        enddo
      enddo
    enddo
  end subroutine PrecomputeMatrices_L_C

  subroutine MatrixElementsAll_CG_0S(Lk, Ll, Ck, Cl, Pbra, Pket, &
                                       Hkl, Skl, Tkl, Vkl, rm2kl, rmkl, rkl, r2kl, deltarkl, drach_deltarkl, &
                             MVkl, drach_MVkl1, drach_MVkl2, drach_MVkl3, Darwinkl, drach_Darwinkl, OOkl, rmrmkl, del2kl, prvalkl, &
                                       wf2originkl, NumCFGridPoints, CFGrid, CFkl, NumDensGridPoints, DensGrid, Denskl, &
                                       AreCorrFuncNeeded, ArePartDensNeeded, AreMCorrFuncNeeded, AreMomDensNeeded)
!This subroutine computes symmetry adapted matrix elements
!with two complex L=0 correlated Gaussians. These matrix elements
!are used in calculations of expectation values.
!Symmetry adaption is applied to the bra and ket using permutation matrices Pbra and Pket
!
!Input:
!   Lk, Ll  :: The lower-triangular matrices Lk and Ll (see PrecomputeMatrices_L_C)
!   Ck, Cl  :: Complex symmetric matrices Ck=Lk*Lk'+i*Bk and Cl=Ll*Ll'+i*Bl
!              (see PrecomputeMatrices_L_C)
!   Pbra    :: The symmetry permutation matrix of size n x n that is applied to bra
!   Pket    :: The symmetry permutation matrix of size n x n that is applied to ket
!Output (all matrix elements are computed with functions that are not normalized;
!function k is the bra, i.e. the function that is complex conjugated):
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
!CFkl              :: Matrix elements of coordinate or momentum-space correlation functions
!NumDensGridPoints :: Number of grid points for particle density calculations
!DensGrid          :: Array containing grid points where matrix elements of
!                     particle densities should be computed
!Denskl            :: Matrix elements of coordinate or momentum-space particle densities
!AreCorrFuncNeeded :: flag indicating whether matrix elements of correlation
!                     functions need to be computed
!ArePartDensNeeded :: flag indicating whether matrix elements of particle
!                     densities need to be computed
!AreMCorrFuncNeeded :: flag indicating whether matrix elements of momentum correlation
!                     functions need to be computed
!AreMomDensNeeded :: flag indicating whether matrix elements of momentum
!                     densities need to be computed
!
!The semi-relativistic kinetic energy, which is available for real Gaussians
!(code RG_0S), is not computed here: with complex Gaussians it requires
!the modified Bessel function of the second kind of a complex argument.

!Arguments
    integer,parameter :: nn=Glob_n
    real(wp),intent(in)      :: Lk(nn,nn),Ll(nn,nn)
    complex(wp),intent(in)   :: Ck(nn,nn),Cl(nn,nn)
    real(wp),intent(in)      :: Pbra(Glob_n,Glob_n),Pket(Glob_n,Glob_n)
    complex(wp),intent(out)  :: Hkl,Skl,Tkl,Vkl,MVkl,drach_MVkl1,drach_MVkl2,drach_MVkl3,Darwinkl,drach_Darwinkl,OOkl
    complex(wp),intent(out)  :: rm2kl(Glob_n,Glob_n),rmkl(Glob_n,Glob_n)
    complex(wp),intent(out)  :: rkl(Glob_n,Glob_n),r2kl(Glob_n,Glob_n)
    complex(wp),intent(out)  :: deltarkl(Glob_n,Glob_n)
    complex(wp),intent(out)  :: drach_deltarkl(Glob_n,Glob_n)
    complex(wp),intent(out)  :: prvalkl(Glob_n,Glob_n)
    complex(wp),intent(out)  :: rmrmkl(Glob_n,Glob_n,Glob_n,Glob_n)
    complex(wp),intent(out)  :: del2kl(Glob_n,Glob_n,Glob_n,Glob_n)
    complex(wp),intent(out)  :: wf2originkl
    integer,intent(in)       :: NumCFGridPoints,NumDensGridPoints
    real(wp),intent(in)      :: CFGrid(NumCFGridPoints),DensGrid(NumDensGridPoints)
    complex(wp),intent(out)  :: CFkl(Glob_n*(Glob_n+1)/2,NumCFGridPoints)
    complex(wp),intent(out)  :: Denskl(Glob_n+1,NumDensGridPoints)
    logical,intent(in)       :: AreCorrFuncNeeded,ArePartDensNeeded
    logical,intent(in)       :: AreMCorrFuncNeeded,AreMomDensNeeded

!Local variables
    integer           n, np
    complex(wp)    tCk(nn,nn), tCl(nn,nn), tCkl(nn,nn)
    complex(wp)    inv_tCkl(nn,nn)
    complex(wp)    inv_tCkltClM(nn,nn), inv_invtCkinvtCl(nn,nn)
    complex(wp)    CtCl(nn,nn), GCl(nn,nn), Cmat(nn,nn)
    complex(wp)    CXs(nn,nn), YCk(nn,nn), CYsk(nn,nn)
    complex(wp)    W1(nn,nn), W2(nn,nn), W3(nn,nn), W4(nn,nn), W5(nn,nn)
    complex(wp)    sk(nn), sl(nn), Csk(nn), Csl(nn)
    complex(wp)    temp1, temp2, temp3, temp4, temp5, temp6, temp7
    complex(wp)    tr1, tr2, tr3, tr4, tr5
    complex(wp)    alpha, beta, gamm, delt, zeta
    complex(wp)    trMtCl, trUMtCl, trXCl, trCXs, trYCk, trCYs, trCYsCXs
    complex(wp)    sqrtdet_tCkl
    integer           i,j,k,indx,p,q
    complex(wp)    TrCJ(nn,nn),sqrtTrCJ(nn,nn),TrCJCJ(nn,nn,nn,nn)
    complex(wp)    Eta(nn+1,nn+1)
    real(wp)       Mass_For_Darwin(0:nn)
    complex(wp)    V2kl

    n=Glob_n
    np=Glob_np
!Then we permute elements of Ck and Cl to account for
!the action of the permutation matrices
!  tCl=Pket'*Cl*Pket
!  tCk=Pbra'*conjg(Ck)*Pbra
!Notice that array tCk contains the complex conjugated matrix: function
!k is the bra and enters all integrals complex conjugated. With this
!definition of tCk all the expressions below look exactly the same as
!the corresponding expressions for real Gaussians with matrices tAk and
!tAl (they are the analytic continuation of the latter).
!We also form matrix tCkl=tCk+tCl
    do i=1,nn
      do j=1,nn
        temp1=ZERO
        temp2=ZERO
        do k=1,nn
          temp1=temp1+Pket(k,j)*Cl(k,i)
          temp2=temp2+conjg(Ck(j,k))*Pbra(k,i)
        enddo
        W1(j,i)=temp1
        W2(j,i)=temp2
      enddo
    enddo
!tCl=W1*Pket
!tCk=Pbra'*W2
    do i=1,nn
      do j=i,nn
        temp1=ZERO
        temp2=ZERO
        do k=1,nn
          temp1=temp1+W1(j,k)*Pket(k,i)
          temp2=temp2+Pbra(k,j)*W2(k,i)
        enddo
        tCl(j,i)=temp1
        tCl(i,j)=temp1
        tCk(j,i)=temp2
        tCk(i,j)=temp2
        tCkl(j,i)=temp1+temp2
        tCkl(i,j)=temp1+temp2
      enddo
    enddo

!After this we can do pseudo-Cholesky factorization of tCkl
!(tCkl is complex symmetric, so it is factorized as a product of a
!lower triangular matrix and its transpose). The pseudo-Cholesky
!factor will be temporarily stored in the lower triangle of W1.
!The product of the diagonal elements of the factor is the properly
!defined square root of det(tCkl) (see the comments in subroutine
!MatrixElementsHS_CG_0S).
    sqrtdet_tCkl=ONE
    do i=1,nn
      do j=i,nn
        temp1=tCkl(i,j)
        do k=i-1,1,-1
          temp1=temp1-W1(i,k)*W1(j,k)
        enddo
        if (i==j) then
          W1(i,i)=sqrt(temp1)
          sqrtdet_tCkl=sqrtdet_tCkl*W1(i,i)
        else
          W1(j,i)=temp1/W1(i,i)
          W1(i,j)=ZERO
        endif
      enddo
    enddo

!Inverting tCkl using its pseudo-Cholesky factor (stored in W1)
!and placing the result into inv_tCkl
    do i=1,nn
      W1(i,i)=ONE/W1(i,i)
      do j=i+1,nn
        temp1=ZERO
        do k=i,j-1
          temp1=temp1-W1(j,k)*W1(k,i)
        enddo
        W1(j,i)=temp1/W1(j,j)
      enddo
    enddo

    do i=1,nn
      do j=i,nn
        temp1=ZERO
        do k=j,nn
          temp1=temp1+W1(k,i)*W1(k,j)
        enddo
        inv_tCkl(i,j)=temp1
        inv_tCkl(j,i)=temp1
      enddo
    enddo

!Evaluating overlap
    wf2originkl=ONE
    Skl=Glob_PiRaised3n2/(sqrtdet_tCkl*sqrtdet_tCkl*sqrtdet_tCkl)

!Doing multiplication CtCl=inv_tCkl*tCl (kept: several sections below
!rely on it and on the identity inv_tCkl*tCk = I - CtCl)
    do i=1,nn
      do j=1,nn
        temp1=ZERO
        do k=1,nn
          temp1=temp1+inv_tCkl(j,k)*tCl(k,i)
        enddo
        CtCl(j,i)=temp1
      enddo
    enddo

!Doing multiplication inv_tCkltClM=inv_tCkl*tCl*M=CtCl*M
    do i=1,nn
      do j=1,nn
        temp1=ZERO
        do k=1,nn
          temp1=temp1+CtCl(j,k)*Glob_MassMatrix(k,i)
        enddo
        inv_tCkltClM(j,i)=temp1
      enddo
    enddo

!GCl = tCl*inv_tCkl*tCl = tCl*CtCl (symmetric). Its elements are the
!bilinear forms tCl(:,i)'*inv_tCkl*tCl(:,j) needed by the orbit-orbit
!and mass-velocity sections below.
    do i=1,nn
      do j=i,nn
        temp1=ZERO
        do k=1,nn
          temp1=temp1+tCl(j,k)*CtCl(k,i)
        enddo
        GCl(j,i)=temp1
        GCl(i,j)=temp1
      enddo
    enddo

!Computing kinetic energy, Tkl=tr[inv_tCkltClM*tCk]
    Tkl=ZERO
    do i=1,nn
      temp1=ZERO
      do k=1,nn
        temp1=temp1+inv_tCkltClM(i,k)*tCk(k,i)
      enddo
      Tkl=Tkl+temp1
    enddo
    Tkl=SIX*Skl*Tkl

!Evaluating potential energy, Vkl,
!(1/r_{ij})_kl, (r_{ij})_kl, (r_{ij}^2)_kl
!and delta(r_{ij})_kl
!The real part of tr[inv_tCkl*Jij] is always positive, so the principal
!values of the square root and of the logarithm are the proper ones.
    temp1=(TWO/Glob_SqrtPi)*Skl
    temp2=THREEHALF*Skl
    temp3=Skl/(Glob_Pi*Glob_SqrtPi)
    Vkl=ZERO
    do i=1,nn
      TrCJ(i,i)=inv_tCkl(i,i)
      sqrtTrCJ(i,i)=sqrt(TrCJ(i,i))
      temp4=sqrtTrCJ(i,i)
      temp5=TrCJ(i,i)
      rmkl(i,i)=temp1/temp4
      rkl(i,i)=temp1*temp4
      r2kl(i,i)=temp2*temp5
      deltarkl(i,i)=temp3/(temp4*temp5)
      prvalkl(i,i)=(temp1/(temp4*temp5))*(Glob_EulerConst+log(temp5))
      Vkl=Vkl+Glob_ScaledPseudoChargeMatrix(i,0)*rmkl(i,i)
    enddo
    do i=1,nn
      do j=i+1,nn
        TrCJ(i,j)=inv_tCkl(i,i)+inv_tCkl(j,j)-inv_tCkl(j,i)-inv_tCkl(j,i)
        TrCJ(j,i)=TrCJ(i,j)
        sqrtTrCJ(j,i)=sqrt(TrCJ(j,i))
        sqrtTrCJ(i,j)=sqrtTrCJ(j,i)
        temp4=sqrtTrCJ(j,i)
        temp5=TrCJ(j,i)
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

!Evaluating tr[inv_tCkl Jij inv_tCkl Jpq]
    do i=1,nn
      temp1=inv_tCkl(i,i)*inv_tCkl(i,i)
      TrCJCJ(i,i,i,i)=temp1
      do p=i+1,nn
        do q=p+1,nn
          temp2=inv_tCkl(p,i)-inv_tCkl(q,i)
          temp1=temp2*temp2
          TrCJCJ(i,i,p,q)=temp1
          TrCJCJ(i,i,q,p)=temp1
          TrCJCJ(p,q,i,i)=temp1
          TrCJCJ(q,p,i,i)=temp1
        enddo
      enddo
      do j=i+1,nn
        temp1=inv_tCkl(j,i)*inv_tCkl(j,i)
        TrCJCJ(j,j,i,i)=temp1
        TrCJCJ(i,i,j,j)=temp1
        do p=i,nn
          temp2=inv_tCkl(p,i)-inv_tCkl(p,j)
          temp1=temp2*temp2
          TrCJCJ(i,j,p,p)=temp1
          TrCJCJ(j,i,p,p)=temp1
          TrCJCJ(p,p,i,j)=temp1
          TrCJCJ(p,p,j,i)=temp1
          do q=p+1,nn
            temp2=inv_tCkl(q,i)-inv_tCkl(p,i)-inv_tCkl(q,j)+inv_tCkl(p,j)
            temp1=temp2*temp2
            TrCJCJ(i,j,p,q)=temp1
            TrCJCJ(j,i,p,q)=temp1
            TrCJCJ(i,j,q,p)=temp1
            TrCJCJ(j,i,q,p)=temp1
            TrCJCJ(p,q,i,j)=temp1
            TrCJCJ(p,q,j,i)=temp1
            TrCJCJ(q,p,i,j)=temp1
            TrCJCJ(q,p,j,i)=temp1
          enddo
        enddo
      enddo
    enddo

!This is a slow old version of the previous loop that computes TrCJCJ
!do i=1,nn
!  do j=1,nn
!    do p=1,nn
!      do q=1,nn
!        W1(1:n,1:n)=ZERO; W1(i,j)=-ONE; W1(j,i)=-ONE; W1(i,i)=ONE; W1(j,j)=ONE
!        W2(1:n,1:n)=ZERO; W2(p,q)=-ONE; W2(q,p)=-ONE; W2(p,p)=ONE; W2(q,q)=ONE;
!        W3(1:n,1:n)=matmul(inv_tCkl(1:n,1:n),W1(1:n,1:n))
!        W4(1:n,1:n)=matmul(inv_tCkl(1:n,1:n),W2(1:n,1:n))
!        W5(1:n,1:n)=matmul(W3(1:n,1:n),W4(1:n,1:n))
!        TrCJCJ(i,j,p,q)=trace(n,W5)
!      enddo
!    enddo
!  enddo
!enddo

!Evaluating (1/r_{ij}*1/r_{pq}))_kl and <delta(r_{ij})delta(r_{pq})>_kl
!in one fused loop (they run over the same index combinations and share
!TrCJ/TrCJCJ). Within the loop ranges (i<=j, p>=i, q>=p) the "same pair"
!condition (p==i.and.q==j).or.(p==j.and.q==i) reduces to p==i.and.q==j.
!Those are exactly the elements of del2kl that must be zero, so writing
!them explicitly makes the full-array zeroing of del2kl unnecessary.
    temp1=4*Skl/Glob_Pi
    temp6=2*Skl
    temp7=Skl/(Glob_Pi**3)
    do i=1,nn
      do j=i,nn
        do p=i,nn
          do q=p,nn
            if ((p==i).and.(q==j)) then
              temp2=temp6/TrCJ(i,j)
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
              temp3=TrCJCJ(i,j,p,q)
              if (temp3/=COMPLEXZERO) then
                !asin(x)/x is an even function of x, so it does not
                !matter which of the two values of the square root of
                !temp3 is used here
                temp4=sqrt(temp3)
                temp5=temp4/(sqrtTrCJ(i,j)*sqrtTrCJ(p,q))
                temp2=temp1*asin(temp5)/temp4
              else
                temp2=temp1/(sqrtTrCJ(i,j)*sqrtTrCJ(p,q))
              endif
              rmrmkl(i,j,p,q)=temp2
              rmrmkl(j,i,p,q)=temp2
              rmrmkl(i,j,q,p)=temp2
              rmrmkl(j,i,q,p)=temp2
              rmrmkl(p,q,i,j)=temp2
              rmrmkl(p,q,j,i)=temp2
              rmrmkl(q,p,i,j)=temp2
              rmrmkl(q,p,j,i)=temp2
              temp2=TrCJ(i,j)*TrCJ(p,q)-temp3
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
    do i=1,nn
      do j=1,nn
        rm2kl(j,i)=rmrmkl(j,i,j,i)
      enddo
    enddo

!Evaluating drachmanized delta(r_{ij})_kl
!W3(i,j) will contain all matrix elements r'(tCk*M*tCl)r / rij.
!With X=tCk*M*tCl one has
!  inv_tCkl*X*inv_tCkl = (I-CtCl)*M*CtCl'   (using inv_tCkl*tCkl = I)
!  tr[inv_tCkl*X]      = tr[M*tCl] - tr[(CtCl*M)*tCl]
!so instead of building X and calling ME_rXr_over_rij_all (which spends
!O(n^2) per pair recomputing rows of inv_tCkl*X) everything follows from
!two matrix products, O(n^2) trace contractions, and O(1) per pair.
    do i=1,nn
      do j=1,nn
        temp1=ZERO
        do k=1,nn
          temp1=temp1+Glob_MassMatrix(j,k)*CtCl(i,k)
        enddo
        W1(j,i)=temp1
      enddo
    enddo
    do i=1,nn
      do j=1,nn
        temp1=ZERO
        do k=1,nn
          temp1=temp1+CtCl(j,k)*W1(k,i)
        enddo
        Cmat(j,i)=W1(j,i)-temp1
      enddo
    enddo
    trMtCl=ZERO
    trUMtCl=ZERO
    do i=1,nn
      do k=1,nn
        trMtCl=trMtCl+Glob_MassMatrix(i,k)*tCl(k,i)
        trUMtCl=trUMtCl+inv_tCkltClM(i,k)*tCl(k,i)
      enddo
    enddo
    temp2=3*(trMtCl-trUMtCl)
    do i=1,nn
      do j=i,nn
        if (i==j) then
          temp3=Cmat(i,i)
        else
          temp3=Cmat(i,i)+Cmat(j,j)-Cmat(i,j)-Cmat(j,i)
        endif
        W3(j,i)=rmkl(j,i)*(temp2-temp3/TrCJ(j,i))*ONEHALF
        W3(i,j)=W3(j,i)
      enddo
    enddo
!Loop that computes all drachmanized delta(r_{ij})_kl  as well as V^2_kl
    V2kl=ZERO
    do p=1,nn
      do q=p,nn
        temp1=ZERO
        do i=1,nn
          temp1=temp1+Glob_ScaledPseudoChargeMatrix(0,i)*rmrmkl(p,q,i,i)
          do j=i+1,nn
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
!with a=tCkl(:,j), b=tCl(:,j), c=tCl(:,i), and d=e_j (first sum) or
!d=e_j-e_i (second sum), are sums of rank-one terms u*v'. Every trace
!formed from them inside ME_rXr_over_rij and ME_rXr_rYr_over_rij
!therefore reduces to O(1) scalars:
!  tr[inv_tCkl*u*v'] = v'*inv_tCkl*u
!  d'*inv_tCkl*u*v'*inv_tCkl*d = (d'*inv_tCkl*u)*(v'*inv_tCkl*d)
!and (with J=d*d', which equals Jjj resp. Jij up to the irrelevant sign
!of d) all the needed building blocks are
!  inv_tCkl*a = e_j (exactly, since inv_tCkl*tCkl = I)
!  alpha = b'*inv_tCkl*c = GCl(j,i)
!  beta  = b'*inv_tCkl*d,  gamm = c'*inv_tCkl*d  (elements of CtCl)
!  delt  = d'*inv_tCkl*d = TrCJ,  a'*inv_tCkl*b = tCl(j,j),
!  a'*inv_tCkl*c = tCl(j,i) = tr1,  a'*inv_tCkl*d = 1.
!In terms of these the two integrals are
!  temp1 = rm*(3*(2*alpha+s) - (2*beta*gamm+delt*s)/delt)/2,
!          s = tr1 + tr1*beta + tr3*gamm
!  temp2 = rm*(9*gamm*zeta + 3*(beta*tr1+alpha)
!          - (3*(gamm*beta+zeta*delt*gamm)
!             + (2*beta*gamm+alpha*delt+beta*tr1*delt))/delt
!          + 3*beta*delt*gamm/delt**2)/4,  zeta = tCl(j,j)
!which makes the whole OO evaluation O(n^2) instead of O(n^5).
    OOkl=ZERO
!First double loop for OO (integral pair index (j,j), d=e_j)
    do i=1,nn
      do j=1,nn
        tr1=tCl(j,i)
        tr3=3*tCl(j,j)
        alpha=GCl(j,i)
        beta=CtCl(j,j)
        gamm=CtCl(j,i)
        delt=TrCJ(j,j)
        zeta=tCl(j,j)
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
    do i=1,nn
      do j=i+1,nn
        tr1=tCl(j,i)
        tr3=3*tCl(j,j)
        alpha=GCl(j,i)
        beta=CtCl(j,j)-CtCl(i,j)
        gamm=CtCl(j,i)-CtCl(i,i)
        delt=TrCJ(i,j)
        zeta=tCl(j,j)
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
!tCk*J*tCk = sk*sk' where sk = tCk*(1,..,1)' (and likewise for tCl),
!and tCk*Jii*tCk = tCk(:,i)*tCk(i,:). Hence
!  tr[inv_tCkl*u*u']            = u'*inv_tCkl*u
!  tr[inv_tCkl*u*u'*inv_tCkl*v*v'] = (u'*inv_tCkl*v)**2
!For the per-particle terms the dot products follow from
!inv_tCkl*tCk = I - CtCl and GCl = tCl*inv_tCkl*tCl:
!  tCk(:,i)'*inv_tCkl*tCk(:,i) = tCk(i,i) - tCl(i,i) + GCl(i,i)
!  tCl(:,i)'*inv_tCkl*tCl(:,i) = GCl(i,i)
!  tCk(:,i)'*inv_tCkl*tCl(:,i) = tCl(i,i) - GCl(i,i)
!This makes the whole section O(n^2) instead of O(n^4).
    do p=1,nn
      temp1=ZERO
      temp2=ZERO
      do q=1,nn
        temp1=temp1+tCk(q,p)
        temp2=temp2+tCl(q,p)
      enddo
      sk(p)=temp1
      sl(p)=temp2
    enddo
    tr1=ZERO
    tr2=ZERO
    do p=1,nn
      tr1=tr1+sk(p)
      tr2=tr2+sl(p)
    enddo
    do p=1,nn
      temp1=ZERO
      temp2=ZERO
      do q=1,nn
        temp1=temp1+inv_tCkl(p,q)*sk(q)
        temp2=temp2+inv_tCkl(p,q)*sl(q)
      enddo
      Csk(p)=temp1
      Csl(p)=temp2
    enddo
    tr4=ZERO
    tr5=ZERO
    temp3=ZERO
    do p=1,nn
      tr4=tr4+sk(p)*Csk(p)
      tr5=tr5+sl(p)*Csl(p)
      temp3=temp3+sk(p)*Csl(p)
    enddo
    temp3=temp3*temp3
    temp2=tr1*tr5+tr2*tr4

    MVkl=36*Skl*(tr4*tr5 + (TWO/THREE)*temp3 - temp2 + tr1*tr2) &
          /(Glob_Mass(1)*Glob_Mass(1)*Glob_Mass(1))

!sum for mass-velocity from 1 to n
    do i=1,nn
      tr1=tCk(i,i)
      tr2=tCl(i,i)
      tr5=GCl(i,i)
      tr4=tCk(i,i)-tCl(i,i)+tr5
      temp3=(tCl(i,i)-tr5)*(tCl(i,i)-tr5)
      temp2=tr1*tr5+tr2*tr4
      MVkl=MVkl+36*Skl*(tr4*tr5 + (TWO/THREE)*temp3 - temp2 + tr1*tr2) &
            /(Glob_Mass(i+1)*Glob_Mass(i+1)*Glob_Mass(i+1))
    enddo

    MVkl=-MVkl/8

!Evaluation of the drachmanized mass-velocity.
!The original helpers myME_dXd_dYd/myME_dXd/myME_over_rij_dXd all
!sandwich a fixed matrix between tCl (X side) or tCk (Y side) and then
!form traces with inv_tCkl. Using inv_tCkl*tCl = CtCl and
!inv_tCkl*tCk = I - CtCl those traces are formed here directly, the
!X = Glob_dmvM pieces are shared between the two myME_dXd_dYd calls,
!and the per-pair myME_over_rij_dXd calls (which redid four O(n^3)
!products for the same X = Glob_dmvB every time) collapse into two
!products plus O(1) work per pair. The commonFactor of the helpers
!is Skl.
!X side (X = Glob_dmvM): W1 = X*tCl, CXs = inv_tCkl*(tCl*X*tCl) = CtCl*W1
    do i=1,nn
      do j=1,nn
        temp1=ZERO
        do k=1,nn
          temp1=temp1+Glob_dmvM(j,k)*tCl(k,i)
        enddo
        W1(j,i)=temp1
      enddo
    enddo
    trXCl=ZERO
    do i=1,nn
      trXCl=trXCl+W1(i,i)
    enddo
    do i=1,nn
      do j=1,nn
        temp1=ZERO
        do k=1,nn
          temp1=temp1+CtCl(j,k)*W1(k,i)
        enddo
        CXs(j,i)=temp1
      enddo
    enddo
    trCXs=ZERO
    do i=1,nn
      trCXs=trCXs+CXs(i,i)
    enddo
!Y side (Y = Glob_dmvM): YCk = Y*tCk, CYsk = inv_tCkl*(tCk*Y*tCk) = (I-CtCl)*YCk
    do i=1,nn
      do j=1,nn
        temp1=ZERO
        do k=1,nn
          temp1=temp1+Glob_dmvM(j,k)*tCk(k,i)
        enddo
        YCk(j,i)=temp1
      enddo
    enddo
    trYCk=ZERO
    do i=1,nn
      trYCk=trYCk+YCk(i,i)
    enddo
    do i=1,nn
      do j=1,nn
        temp1=ZERO
        do k=1,nn
          temp1=temp1+CtCl(j,k)*YCk(k,i)
        enddo
        CYsk(j,i)=YCk(j,i)-temp1
      enddo
    enddo
    trCYs=ZERO
    do i=1,nn
      trCYs=trCYs+CYsk(i,i)
    enddo
    trCYsCXs=ZERO
    do i=1,nn
      do j=1,nn
        trCYsCXs=trCYsCXs+CYsk(i,j)*CXs(j,i)
      enddo
    enddo
    temp1 = 16._wp*Skl*(NINE/FOUR*trCXs*trCYs+THREE/TWO*trCYsCXs) &
            - 36._wp*Skl*trXCl*trCYs - 36._wp*Skl*trYCk*trCXs &
            + 36._wp*Skl*trYCk*trXCl &
            - V2kl - Glob_CurrEnergy*Glob_CurrEnergy*Skl + 2*Glob_CurrEnergy*Vkl
    if (.not. Glob_ArePseudoParticleMassesTheSame) then
!Second myME_dXd_dYd call: only the Y side changes (Y = Glob_dmvMB)
      do i=1,nn
        do j=1,nn
          temp2=ZERO
          do k=1,nn
            temp2=temp2+Glob_dmvMB(j,k)*tCk(k,i)
          enddo
          YCk(j,i)=temp2
        enddo
      enddo
      trYCk=ZERO
      do i=1,nn
        trYCk=trYCk+YCk(i,i)
      enddo
      do i=1,nn
        do j=1,nn
          temp2=ZERO
          do k=1,nn
            temp2=temp2+CtCl(j,k)*YCk(k,i)
          enddo
          CYsk(j,i)=YCk(j,i)-temp2
        enddo
      enddo
      trCYs=ZERO
      do i=1,nn
        trCYs=trCYs+CYsk(i,i)
      enddo
      trCYsCXs=ZERO
      do i=1,nn
        do j=1,nn
          trCYsCXs=trCYsCXs+CYsk(i,j)*CXs(j,i)
        enddo
      enddo
      temp2 = 16._wp*Skl*(NINE/FOUR*trCXs*trCYs+THREE/TWO*trCYsCXs) &
              - 36._wp*Skl*trXCl*trCYs - 36._wp*Skl*trYCk*trCXs &
              + 36._wp*Skl*trYCk*trXCl &
              - V2kl - Glob_CurrEnergy*Glob_CurrEnergy*Skl + 2*Glob_CurrEnergy*Vkl
!myME_dXd(Glob_dmvB): 6*Skl*(tr[inv_tCkl*tCl*B*tCl] - tr[B*tCl]).
!W1 = B*tCl; the first trace is the O(n^2) contraction of CtCl with W1.
      do i=1,nn
        do j=1,nn
          temp3=ZERO
          do k=1,nn
            temp3=temp3+Glob_dmvB(j,k)*tCl(k,i)
          enddo
          W1(j,i)=temp3
        enddo
      enddo
      temp5=ZERO
      temp6=ZERO
      do i=1,nn
        temp5=temp5+W1(i,i)
        do j=1,nn
          temp6=temp6+CtCl(i,j)*W1(j,i)
        enddo
      enddo
      temp2 = temp2 + Glob_CurrEnergy*SIX*Skl*(temp6-temp5)
!All-pairs part of myME_over_rij_dXd(Glob_dmvB,p,q,...):
!inv_tCkl*(tCl*B*tCl)*inv_tCkl = CtCl*B*CtCl' = CtCl*W4 with
!W4 = B*CtCl'; per pair only gamma=1/sqrtTrCJ and an O(1) contraction
!of W5 = CtCl*W4 remain. trCXs-trXCl of the helper is temp6-temp5.
      do i=1,nn
        do j=1,nn
          temp3=ZERO
          do k=1,nn
            temp3=temp3+Glob_dmvB(j,k)*CtCl(i,k)
          enddo
          W4(j,i)=temp3
        enddo
      enddo
      do i=1,nn
        do j=1,nn
          temp3=ZERO
          do k=1,nn
            temp3=temp3+CtCl(j,k)*W4(k,i)
          enddo
          W5(j,i)=temp3
        enddo
      enddo
      temp4=SIX*TWO*Skl/Glob_SqrtPi
      do i=1,nn
        do j=i,nn
          if (i==j) then
            temp3=W5(i,i)
            temp7=Glob_ScaledPseudoChargeMatrix(0,i)
          else
            temp3=W5(i,i)+W5(j,j)-W5(i,j)-W5(j,i)
            temp7=Glob_ScaledPseudoChargeMatrix(i,j)
          endif
          temp2=temp2-temp7*temp4*( (temp6-temp5)/sqrtTrCJ(j,i) &
                -(ONE/THREE)*temp3/(sqrtTrCJ(j,i)*TrCJ(j,i)) )
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
    do i=1,nn
      Darwinkl=Darwinkl+(   &
                ONE/(Mass_For_Darwin(0)*Mass_For_Darwin(0)) &
                +ONE/(Mass_For_Darwin(i)*Mass_For_Darwin(i)) &
                )*Glob_ScaledPseudoChargeMatrix(0,i)*deltarkl(i,i)
    enddo
    do i=1,nn
      do j=1,nn
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
    do i=1,nn
      drach_Darwinkl=drach_Darwinkl+(   &
                      ONE/(Mass_For_Darwin(0)*Mass_For_Darwin(0)) &
                      +ONE/(Mass_For_Darwin(i)*Mass_For_Darwin(i)) &
                      )*Glob_ScaledPseudoChargeMatrix(0,i)*drach_deltarkl(i,i)
    enddo
    do i=1,nn
      do j=1,nn
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
      do i=1,nn
        do j=i,nn
          p=p+1
          temp3=temp1/(sqrtTrCJ(j,i)*TrCJ(j,i))
          temp5=ONE/TrCJ(j,i)
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
      do i=1,nn+1
        temp3=ZERO
        do p=1,nn
          temp3=temp3+Glob_bvc(p,i)*Glob_bvc(p,i)*inv_tCkl(p,p)
          do q=p+1,nn
            temp3=temp3+2*Glob_bvc(q,i)*Glob_bvc(p,i)*inv_tCkl(q,p)
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

    !inv_invtCkinvtCl=(inv(tCk)+inv(tCl))^{-1}
    !                =tCk*inv(tCkl)*tCl=tCk*CtCl, symmetrized.
    !The equivalent form tCl-tCl*inv(tCkl)*tCl is not used: when tCl has
    !much larger eigenvalues than tCk it subtracts nearly equal numbers,
    !and roundoff can make the momentum widths wrong or even non-positive.
    !The real part of this matrix is positive definite (as the matrix is
    !the inverse of a sum of two matrices with positive definite real
    !parts), so the real parts of the momentum widths below are positive.
    do i=1,nn
      do j=i,nn
        temp1=ZERO
        temp2=ZERO
        do k=1,nn
          temp1=temp1+tCk(i,k)*CtCl(k,j)
          temp2=temp2+tCk(j,k)*CtCl(k,i)
        enddo
        inv_invtCkinvtCl(i,j)=ONEHALF*(temp1+temp2)
        inv_invtCkinvtCl(j,i)=inv_invtCkinvtCl(i,j)
      enddo
    enddo
    ! Construct
    !   Eta(I,J) = 4*a_I^T R_KL a_J
    ! where
    !   a_1  = -(1,1,...,1)^T,
    !   a_I  = e_(I-1),  I = 2,...,n+1.
    ! Eta has dimensions (n+1,n+1).
    Eta(1,1) = ZERO
    do i = 1, nn
      temp1 = ZERO
      do j = 1, nn
        temp2 = inv_invtCkinvtCl(j,i)
        ! Eta(j+1,i+1) = R_KL(j,i)
        Eta(j+1,i+1) = temp2
        ! Required for Eta(1,i+1)
        temp1 = temp1 + temp2
      enddo
      ! Eta(1,i+1) = -1^T R_KL e_i
      Eta(1,i+1) = -temp1
      ! R_KL is symmetric, so Eta(i+1,1) is identical.
      Eta(i+1,1) = -temp1
      ! Eta(1,1) = sum over all elements of R_KL.
      Eta(1,1) = Eta(1,1) + temp1
    enddo


    if (AreMCorrFuncNeeded) then
      temp1= ONEFOURTH * ONEHALF * Skl/(Glob_Pi*Glob_SqrtPi)
      p = 0
      do i = 1, nn
          do j = i, nn
              p = p + 1
              if (i == j) then
                  temp2 = Eta(1,1) + Eta(i+1,i+1) &
                        - TWO*Eta(1,i+1)
              else
                  temp2 = Eta(i+1,i+1) + Eta(j+1,j+1) &
                        - TWO*Eta(i+1,j+1)
              endif

              if (real(temp2,wp) <= ZERO) then
                write(*,'(1x,a,1x,i0)') 'Non-positive momentum correlation width on MPI rank',Glob_ProcID
                call MPI_Abort(MPI_COMM_WORLD,1,Glob_MPIErrCode)
              endif

              temp3 = temp1/(temp2*sqrt(temp2))
              do k = 1, NumCFGridPoints
                  temp4 = CFGrid(k)**2
                  CFkl(p,k) = temp3*exp(- ONEFOURTH * temp4/temp2)
              enddo
          enddo
      enddo
    end if

    if (AreMomDensNeeded) then
      temp1= ONEFOURTH * ONEHALF * Skl/(Glob_Pi*Glob_SqrtPi)
      !cycle through all particles
      do i = 1, nn+1
        temp2 = Eta(i,i)

        if (real(temp2,wp) <= ZERO) then
          write(*,'(1x,a,1x,i0)') 'Non-positive momentum density width on MPI rank',Glob_ProcID
          call MPI_Abort(MPI_COMM_WORLD,1,Glob_MPIErrCode)
        endif

        !prefactor: Skl * (1/(4*pi*eta))^(3/2)
        temp3 = temp1 / (temp2 * sqrt(temp2))
        !compute the matrix element value at all density grid points
        do k = 1, NumDensGridPoints
          !value of xi^2 at a specific point
          temp4 = DensGrid(k)**2
          Denskl(i,k) = temp3 * exp(- ONEFOURTH * temp4 / temp2)
        enddo
      enddo
    endif

  end subroutine MatrixElementsAll_CG_0S

  function trace(k,M)
    complex(wp) trace
    integer k
    integer,parameter :: nn=Glob_n
    complex(wp) M(nn,nn)
    integer i
    trace=ZERO
    do i=1,k
      trace=trace+M(i,i)
    enddo
  end function trace

  function ME_rXr_over_rij(X,i,j,inv_tCkl,ME_1_over_rij,TrCJ)
!Function ME_rXr_over_rij computes the following
!matrix element with complex L=0 Gaussians phi_k and phi_l:
!<phi_k| r'Xr/r_{ij} |phi_l>
!Here X is an arbitrary (i.e. nonsymmetric) complex matrix.
!Index i can be equal to j. In the latter case
!<phi_k| r'Xr/r_{i} |phi_l> is computed
!Input:
!            X  :: n x n complex matrix
!           i,j :: indices denoting i and j.
!        inv_tCkl :: n x n complex matrix where the inverse of tCkl is stored
! ME_1_over_rij :: the value of <phi_k| 1/r_{ij} |phi_l> matrix element
!          TrCJ :: the value of Tr[inv_tCkl Jij]
!All arrays (both arguments and local ones) are static and have dimension
!nn x nn, where nn=Glob_n.
    complex(wp) ME_rXr_over_rij
    integer,parameter :: nn=Glob_n

!Arguments:
    complex(wp)  X(nn,nn),inv_tCkl(nn,nn),ME_1_over_rij,TrCJ
    integer      i,j
!Local variables
    complex(wp)  CXi(nn),CXj(nn)
    complex(wp)  TrCX,TrCXCJ
    integer      n,m,p,q

    n=Glob_n

    TrCX=ZERO
    do m=1,nn
      do p=1,nn
        TrCX=TrCX+inv_tCkl(m,p)*X(p,m)
      enddo
    enddo

    do m=1,nn
      CXi(m)=ZERO
      do p=1,nn
        CXi(m)=CXi(m)+inv_tCkl(i,p)*X(p,m)
      enddo
    enddo
    if (j/=i) then
      do m=1,nn
        CXj(m)=ZERO
        do p=1,nn
          CXj(m)=CXj(m)+inv_tCkl(j,p)*X(p,m)
        enddo
      enddo
    endif

    if (i==j) then
      TrCXCJ=ZERO
      do m=1,nn
        TrCXCJ=TrCXCJ+CXi(m)*inv_tCkl(m,i)
      enddo
    else
      TrCXCJ=ZERO
      do m=1,nn
        TrCXCJ=TrCXCJ+(CXi(m)-CXj(m))*(inv_tCkl(m,i)-inv_tCkl(m,j))
      enddo
    endif

    ME_rXr_over_rij=ME_1_over_rij*(3*TrCX-TrCXCJ/TrCJ)/2

  end function ME_rXr_over_rij

  subroutine ME_rXr_over_rij_all(X,inv_tCkl,rmkl,TrCJ,ME)
!Subroutine ME_rXr_over_rij_all computes the following
!matrix elements with complex L=0 Gaussians phi_k and phi_l:
!<phi_k| r'Xr/r_{ij} |phi_l>    for all combinations i,j=1..n at the same time
!Index i can be equal to j. In the latter case
!<phi_k| r'Xr/r_{i} |phi_l> is computed
!Input:
!            X  :: n x n complex matrix
!      inv_tCkl :: n x n complex matrix where the inverse of tCkl is stored
!          rmkl :: the values of <phi_k| 1/r_{ij} |phi_l> matrix element
!          TrCJ :: the values of Tr[inv_tCkl Jij]
!Output:
!            ME :: n x n complex matrix where all computed matrix elements are returned
!All arrays (both arguments and local ones) are static and have dimension
!nn x nn, where nn=Glob_n.
    integer,parameter :: nn=Glob_n

!Arguments:
    complex(wp)  X(nn,nn),inv_tCkl(nn,nn),rmkl(Glob_n,Glob_n),TrCJ(nn,nn),ME(nn,nn)
!Local variables
    integer      i,j
    complex(wp)  CXi(nn),CXj(nn)
    complex(wp)  TrCX,TrCXCJ
    integer      n,m,p,q

    n=Glob_n

    TrCX=ZERO
    do m=1,nn
      do p=1,nn
        TrCX=TrCX+inv_tCkl(m,p)*X(p,m)
      enddo
    enddo

    do i=1,nn
      do j=i,nn
        do m=1,nn
          CXi(m)=ZERO
          do p=1,nn
            CXi(m)=CXi(m)+inv_tCkl(i,p)*X(p,m)
          enddo
        enddo
        if (j/=i) then
          do m=1,nn
            CXj(m)=ZERO
            do p=1,nn
              CXj(m)=CXj(m)+inv_tCkl(j,p)*X(p,m)
            enddo
          enddo
        endif

        if (i==j) then
          TrCXCJ=ZERO
          do m=1,nn
            TrCXCJ=TrCXCJ+CXi(m)*inv_tCkl(m,i)
          enddo
        else
          TrCXCJ=ZERO
          do m=1,nn
            TrCXCJ=TrCXCJ+(CXi(m)-CXj(m))*(inv_tCkl(m,i)-inv_tCkl(m,j))
          enddo
        endif

        ME(i,j)=rmkl(i,j)*(3*TrCX-TrCXCJ/TrCJ(i,j))/2
        ME(j,i)=ME(i,j)

      enddo
    enddo

  end subroutine ME_rXr_over_rij_all

  function ME_rXr_rYr_over_rij(X,Y,i,j,inv_tCkl,ME_1_over_rij,TrCJ)
!Function ME_rXr_rYr_over_rij computes the following
!matrix element with complex L=0 Gaussians phi_k and phi_l:
!<phi_k| (r'Xr)(r'Yr)/r_{ij} |phi_l>
!Here X and Y are arbitrary (i.e. nonsymmetric) complex matrices.
!Index i can be equal to j. In the latter case
!<phi_k| (r'Xr)(r'Yr)/r_{i} |phi_l> is computed
!Input:
!            X  :: n x n complex matrix
!            Y  :: n x n complex matrix
!           i,j :: indices denoting i and j.
!        inv_tCkl :: n x n complex matrix where the inverse of tCkl is stored
! ME_1_over_rij :: the value of <phi_k| 1/r_{ij} |phi_l> matrix element
!          TrCJ :: the value of Tr[inv_tCkl Jij]
!All arrays (both arguments and local ones) are static and have dimension
!nn x nn, where nn=Glob_n.
    complex(wp) ME_rXr_rYr_over_rij
    integer,parameter :: nn=Glob_n

!Arguments:
    complex(wp)  X(nn,nn),Y(nn,nn),inv_tCkl(nn,nn),ME_1_over_rij,TrCJ
    integer      i,j
!Local variables
    complex(wp)  Ys(nn,nn),Xs(nn,nn),CX(nn,nn),CY(nn,nn),CXCYi(nn),CXCYj(nn)
    complex(wp)  TrCX,TrCY,TrCXCY,TrCXCJ,TrCYCJ,TrCXCYCJ
    integer      n,m,p,q

    n=Glob_n

    do p=1,nn
      Ys(p,p)=Y(p,p)
      do q=p+1,nn
        Ys(p,q)=(Y(p,q)+Y(q,p))/2
        Ys(q,p)=Ys(p,q)
      enddo
    enddo
    do p=1,nn
      Xs(p,p)=X(p,p)
      do q=p+1,nn
        Xs(p,q)=(X(p,q)+X(q,p))/2
        Xs(q,p)=Xs(p,q)
      enddo
    enddo

    do p=1,nn
      do q=1,nn
        CY(p,q)=ZERO
        do m=1,nn
          CY(p,q)=CY(p,q)+inv_tCkl(p,m)*Ys(m,q)
        enddo
      enddo
    enddo
    do p=1,nn
      do q=1,nn
        CX(p,q)=ZERO
        do m=1,nn
          CX(p,q)=CX(p,q)+inv_tCkl(p,m)*Xs(m,q)
        enddo
      enddo
    enddo

    do m=1,nn
      CXCYi(m)=ZERO
      do p=1,nn
        CXCYi(m)=CXCYi(m)+CX(i,p)*CY(p,m)
      enddo
    enddo
    if (j/=i) then
      do m=1,nn
        CXCYj(m)=ZERO
        do p=1,nn
          CXCYj(m)=CXCYj(m)+CX(j,p)*CY(p,m)
        enddo
      enddo
    endif

    TrCY=ZERO
    TrCX=ZERO
    do m=1,nn
      TrCY=TrCY+CY(m,m)
      TrCX=TrCX+CX(m,m)
    enddo
    TrCXCY=ZERO
    do m=1,nn
      do p=1,nn
        TrCXCY=TrCXCY+CX(m,p)*CY(p,m)
      enddo
    enddo

    if (i==j) then
      TrCXCJ=ZERO
      do m=1,nn
        TrCXCJ=TrCXCJ+CX(i,m)*inv_tCkl(m,i)
      enddo
      TrCYCJ=ZERO
      do m=1,nn
        TrCYCJ=TrCYCJ+CY(i,m)*inv_tCkl(m,i)
      enddo
      TrCXCYCJ=ZERO
      do m=1,nn
        TrCXCYCJ=TrCXCYCJ+CXCYi(m)*inv_tCkl(m,i)
      enddo
    else
      TrCXCJ=ZERO
      do m=1,nn
        TrCXCJ=TrCXCJ+(CX(i,m)-CX(j,m))*(inv_tCkl(m,i)-inv_tCkl(m,j))
      enddo
      TrCYCJ=ZERO
      do m=1,nn
        TrCYCJ=TrCYCJ+(CY(i,m)-CY(j,m))*(inv_tCkl(m,i)-inv_tCkl(m,j))
      enddo
      TrCXCYCJ=ZERO
      do m=1,nn
        TrCXCYCJ=TrCXCYCJ+(CXCYi(m)-CXCYj(m))*(inv_tCkl(m,i)-inv_tCkl(m,j))
      enddo
    endif

    ME_rXr_rYr_over_rij=ME_1_over_rij*(   &
                         9*TrCY*TrCX + 6*TrCXCY           &
                         - ( 3*(TrCY*TrCXCJ + TrCX*TrCYCJ) + 4*TrCXCYCJ )/TrCJ   &
                         + 3*TrCXCJ*TrCYCJ/(TrCJ*TrCJ)    &
                         )/4

  end function ME_rXr_rYr_over_rij

  function myME_dXd_dYd(X,Y,inv_tCkl,tCk,tCl,Skl)
!Here tCk must contain the complex conjugated (and permuted) matrix
!of the bra, tCkl=tCk+tCl, and Skl is the overlap matrix element.

    complex(wp) myME_dXd_dYd, Skl
    integer,parameter :: nn=Glob_n
    complex(wp) X(nn,nn),Y(nn,nn),inv_tCkl(nn,nn),tCk(nn,nn),tCl(nn,nn)

!Local variables:
    integer :: i,j,k,n
    complex(wp) :: temp, trCXs, trCYs, trXCl, trYCk, trCYsCXs, commonFactor
    complex(wp) :: XCl(nn,nn), Xs(nn,nn), YCk(nn,nn), Ys(nn,nn),&
                   CXs(nn,nn), YsCXs(nn,nn)
    complex(wp) :: term1, term2, term3, term4

    n=Glob_n
!!! Q-part  !!!
! Build Xs matrix
    XCl = ZERO
    do i=1,nn
      do j=1,nn
        temp = ZERO
        do k=1,nn
          temp = temp + X(i,k)*tCl(k,j)
        enddo
        XCl(i,j) = temp
      enddo
    enddo

    Xs = ZERO
    do i=1,nn
      do j=1,nn
        temp = ZERO
        do k=1,nn
          temp = temp + tCl(i,k) * XCl(k,j)
        enddo
        Xs(i,j) = temp
      enddo
    enddo

    YCk = ZERO
    do i=1,nn
      do j=1,nn
        temp = ZERO
        do k=1,nn
          temp = temp + Y(i,k)*tCk(k,j)
        enddo
        YCk(i,j) = temp
      enddo
    enddo

    Ys = ZERO
    do i=1,nn
      do j=1,nn
        temp = ZERO
        do k=1,nn
          temp = temp + tCk(i,k) * YCk(k,j)
        enddo
        Ys(i,j) = temp
      enddo
    enddo

!Symmetrize
    do i = 1,nn
      do j = i+1,nn
        temp=ONEHALF*(Xs(j,i)+Xs(i,j))
        Xs(j,i) = temp
        Xs(i,j) = temp
      enddo
    enddo

    do i = 1,nn
      do j = i+1,nn
        temp=ONEHALF*(Ys(j,i)+Ys(i,j))
        Ys(j,i) = temp
        Ys(i,j) = temp
      enddo
    enddo
!End Symmetrize

    CXs = ZERO
    do i=1,nn
      do j=1,nn
        temp = ZERO
        do k=1,nn
          temp = temp + inv_tCkl(i,k)*Xs(k,j)
        enddo
        CXs(i,j) = temp
      enddo
    enddo

    YsCXs = ZERO
    do i=1,nn
      do j=1,nn
        temp = ZERO
        do k=1,nn
          temp = temp + Ys(i,k)*CXs(k,j)
        enddo
        YsCXs(i,j) = temp
      enddo
    enddo

    trCXs = ZERO
    do i=1,nn
      do j=1,nn
        trCXs = trCXs + inv_tCkl(i,j)*Xs(j,i)
      enddo
    enddo

    trCYs = ZERO
    do i=1,nn
      do j=1,nn
        trCYs = trCYs + inv_tCkl(i,j)*Ys(j,i)
      enddo
    enddo

    trXCl = ZERO
    do i=1,nn
      trXCl = trXCl + XCl(i,i)
    enddo

    trYCk = ZERO
    do i=1,nn
      trYCk = trYCk + YCk(i,i)
    enddo

    trCYsCXs = ZERO
    do i=1,nn
      do j=1,nn
        trCYsCXs = trCYsCXs + inv_tCkl(i,j)*YsCXs(j,i)
      enddo
    enddo

    commonFactor = Skl

    term1 = 16._wp*commonFactor*(&
            NINE/FOUR*trCXs*trCYs+ THREE/TWO*trCYsCXs)
    term2 = -36._wp*commonFactor*trXCl*trCYs
    term3 = -36._wp*commonFactor*trYCk*trCXs
    term4 = 36._wp*commonFactor*trYCk*trXCl

    myME_dXd_dYd = term1 + term2 + term3 + term4

  end function myME_dXd_dYd

  function myME_dXd(X,inv_tCkl,tCl,Skl)
!Here Skl is the overlap matrix element.

    complex(wp) myME_dXd, Skl
    integer,parameter :: nn=Glob_n
    complex(wp) X(nn,nn),inv_tCkl(nn,nn),tCl(nn,nn)

!Local variables:
    integer :: i,j,k,n
    complex(wp) :: temp, trXCl, trCXs, commonFactor
    complex(wp) :: XCl(nn,nn), Xs(nn,nn)

    n=Glob_n
!!! Q-part  !!!
    !Build Xs matrix
    XCl = ZERO
    do i=1,nn
      do j=1,nn
        temp = ZERO
        do k=1,nn
          temp = temp + X(i,k)*tCl(k,j)
        enddo
        XCl(i,j) = temp
      enddo
    enddo

    Xs = ZERO
    do i=1,nn
      do j=1,nn
        temp = ZERO
        do k=1,nn
          temp = temp + tCl(i,k) * XCl(k,j)
        enddo
        Xs(i,j) = temp
      enddo
    enddo

!Symmetrize
    do i = 1,nn
      do j = i+1,nn
        temp=ONEHALF*(Xs(j,i)+Xs(i,j))
        Xs(j,i) = temp
        Xs(i,j) = temp
      enddo
    enddo

    trXCl = ZERO
    do i=1,nn
      trXCl = trXCl + XCl(i,i)
    enddo

    trCXs = ZERO
    do i=1,nn
      do j=1,nn
        trCXs = trCXs + inv_tCkl(i,j)*Xs(j,i)
      enddo
    enddo

    commonFactor = SIX*Skl

    myME_dXd = commonFactor*(trCXs - trXCl)

  end function myME_dXd

  function myME_over_rij_dXd(X,p,q,inv_tCkl,tCl,Skl)
!Function ME_1_over_rij_dXd computes the following
!matrix element with complex L=0 Gaussians phi_k and phi_l:
!<phi_k| (1/r_{ij})(nabla_r'*X*nabla_r) |phi_l>
!Here X is arbitrary (i.e. nonsymmetric) complex matrix.
!Index i can be equal to j. In the latter case
!<phi_k| (1/r_i)(nabla_r'*X*nabla_r) |phi_l>is computed
!Input:
!            X  :: n x n complex matrix
!           i,j :: indices denoting i and j.
!      inv_tCkl :: n x n complex matrix where the inverse of tCkl is stored
!           Skl :: the value of the overlap matrix element <phi_k|phi_l>
!All arrays (both arguments and local ones) are static and have dimension
!nn x nn, where nn=Glob_n.

!Input parameters:
    complex(wp) myME_over_rij_dXd, Skl
    integer,parameter :: nn=Glob_n
    complex(wp) :: X(nn,nn),inv_tCkl(nn,nn),tCl(nn,nn)
    integer p,q

!Local variables:
    integer :: i,j,k,n
    complex(wp) :: temp, gamma, trCXs, jijCXsCjij, trXCl, commonFactor
    complex(wp) :: XCl(nn,nn), Xs(nn,nn), XsC(nn,nn), CXsC(nn,nn)

    n=Glob_n
!!! Q-part  !!!
    !Build Xs matrix
    XCl = ZERO
    do i=1,nn
      do j=1,nn
        temp = ZERO
        do k=1,nn
          temp = temp + X(i,k)*tCl(k,j)
        enddo
        XCl(i,j) = temp
      enddo
    enddo

    Xs = ZERO
    do i=1,nn
      do j=1,nn
        temp = ZERO
        do k=1,nn
          temp = temp + tCl(i,k) * XCl(k,j)
        enddo
        Xs(i,j) = temp
      enddo
    enddo

!Symmetrize
    do i = 1,nn
      do j = i+1,nn
        temp=ONEHALF*(Xs(j,i)+Xs(i,j))
        Xs(j,i) = temp
        Xs(i,j) = temp
      enddo
    enddo

    XsC = ZERO
    do i=1,nn
      do j=1,nn
        temp = ZERO
        do k=1,nn
          temp = temp + Xs(i,k)*inv_tCkl(k,j)
        enddo
        XsC(i,j) = temp
      enddo
    enddo

    CXsC = ZERO
    do i=1,nn
      do j=1,nn
        temp = ZERO
        do k=1,nn
          temp = temp + inv_tCkl(i,k)*XsC(k,j)
        enddo
        CXsC(i,j) = temp
      enddo
    enddo

    jijCXsCjij = ZERO

    if (p == q) then
      jijCXsCjij = CXsC(p,p)
      gamma = inv_tCkl(p,p)
      gamma = ONE/sqrt(gamma)
    else
      jijCXsCjij = CXsC(p,p) + CXsC(q,q) - CXsC(p,q) - CXsC(q,p)
      gamma = inv_tCkl(p,p) + inv_tCkl(q,q) - inv_tCkl(p,q) - inv_tCkl(q,p)
      gamma = ONE/sqrt(gamma)
    endif

    trCXs = ZERO
    do i=1,nn
      do j=1,nn
        trCXs = trCXs + inv_tCkl(i,j)*Xs(j,i)
      enddo
    enddo

    trXCl = ZERO
    do i=1,nn
      trXCl = trXCl + XCl(i,i)
    enddo

    commonFactor = SIX*TWO*Skl/Glob_SqrtPi

    myME_over_rij_dXd = commonFactor*(&
                        gamma*(trCXs - trXCl) - ONE/THREE*gamma**3*jijCXsCjij)

  end function myME_over_rij_dXd

  function ME_1_over_rij_dXd(X,i,j,inv_tCkl,tCl,ME_1_over_rij,TrCJ)
!Function ME_1_over_rij_dXd computes the following
!matrix element with complex L=0 Gaussians phi_k and phi_l:
!<phi_k| (1/r_{ij})(nabla_r'*X*nabla_r) |phi_l>
!Here X and Y are arbitrary (i.e. nonsymmetric) complex matrices.
!Index i can be equal to j. In the latter case
!<phi_k| (1/r_i)(nabla_r'*X*nabla_r) |phi_l>is computed
!Input:
!            X  :: n x n complex matrix
!            Y  :: n x n complex matrix
!           i,j :: indices denoting i and j.
!      inv_tCkl :: n x n complex matrix where the inverse of tCkl is stored
! ME_1_over_rij :: the value of <phi_k| 1/r_{ij} |phi_l> matrix element
!          TrCJ :: the value of Tr[inv_tCkl Jij]
!All arrays (both arguments and local ones) are static and have dimension
!nn x nn, where nn=Glob_n.
!Input parameters:
    complex(wp) ME_1_over_rij_dXd
    integer,parameter :: nn=Glob_n
    complex(wp) X(nn,nn),inv_tCkl(nn,nn),tCl(nn,nn),ME_1_over_rij,TrCJ
!Local variables:
    integer i,j,p,q,k,n
    complex(wp) M(nn,nn),Z(nn,nn)
    complex(wp) tr1,t

    n=Glob_n
!Compute Z=tCl*X, tr1=trace[Z], and M=tCl*X*tCl
    tr1=ZERO
    do p=1,nn
      do q=1,nn
        t=ZERO
        do k=1,nn
          t=t+tCl(q,k)*X(k,p)
        enddo
        Z(q,p)=t
      enddo
      tr1=tr1+Z(p,p)
    enddo
    do p=1,nn
      do q=1,nn
        t=ZERO
        do k=1,nn
          t=t+Z(q,k)*tCl(k,p)
        enddo
        M(q,p)=t
      enddo
    enddo
!Compute the matrix element
    ME_1_over_rij_dXd=4*ME_rXr_over_rij(M,i,j,inv_tCkl,ME_1_over_rij,TrCJ)-6*tr1*ME_1_over_rij
  end function ME_1_over_rij_dXd

  subroutine ME_1_over_rij_dXd_all(X,inv_tCkl,tCl,rmkl,TrCJ,ME)
!Subroutine ME_1_over_rij_dXd computes
!matrix elements with complex L=0 Gaussians phi_k and phi_l:
!<phi_k| (1/r_{ij})(nabla_r'*X*nabla_r) |phi_l>
!for all combinations of i and j indexes.
!Here X and Y are arbitrary (i.e. nonsymmetric) complex matrices.
!When i is equal to j it stands for the matrix elements
!<phi_k| (1/r_i)(nabla_r'*X*nabla_r) |phi_l>
!Input:
!            X  :: n x n complex matrix
!            Y  :: n x n complex matrix
!           i,j :: indices denoting i and j.
!      inv_tCkl :: n x n complex matrix where the inverse of tCkl is stored
!          rmkl :: the values of <phi_k| 1/r_{ij} |phi_l> matrix element
!          TrCJ :: the values of Tr[inv_tCkl Jij]
!Output:
!            ME :: n x n complex matrix where all computed matrix elements are returned
!All arrays (both arguments and local ones) are static and have dimension
!nn x nn, where nn=Glob_n.
!Input parameters:
    integer,parameter :: nn=Glob_n
    complex(wp) X(nn,nn),inv_tCkl(nn,nn),tCl(nn,nn),rmkl(Glob_n,Glob_n),TrCJ(nn,nn)
    complex(wp) ME(nn,nn)
!Local variables:
    integer i,j,p,q,k,n
    complex(wp) M(nn,nn),Z(nn,nn)
    complex(wp) tr1,t

    n=Glob_n
!Compute Z=tCl*X, tr1=trace[Z], and M=tCl*X*tCl
    tr1=ZERO
    do p=1,nn
      do q=1,nn
        t=ZERO
        do k=1,nn
          t=t+tCl(q,k)*X(k,p)
        enddo
        Z(q,p)=t
      enddo
      tr1=tr1+Z(p,p)
    enddo
    do p=1,nn
      do q=1,nn
        t=ZERO
        do k=1,nn
          t=t+Z(q,k)*tCl(k,p)
        enddo
        M(q,p)=t
      enddo
    enddo
!Compute all matrix elements <phi_k| r'Mr/r_{ij} |phi_l>
    call ME_rXr_over_rij_all(M,inv_tCkl,rmkl,TrCJ,Z)
    t=6*tr1
    do i=1,nn
      do j=i,nn
        ME(j,i)=4*Z(j,i)-t*rmkl(j,i)
        ME(i,j)=ME(j,i)
      enddo
    enddo

  end subroutine ME_1_over_rij_dXd_all

  subroutine spinPreCalc(n, nFactorial, parityFactor, SSFmassChargeCoefficient, &
                         AnihMassChargeCoefficient, ketMatrix, spatialYoung, &
                         positronPosition, numberOfSpinFunctions, spinFreeME, SiSjME)
    use spinStuff
    implicit none
    integer,parameter :: nn=Glob_n

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
    do i = 1, nn
      do j = 1, nn
        SSFmassChargeCoefficient(i, j) = -Glob_PseudoCharge(i) * Glob_PseudoCharge(j) / &
                                         (Glob_Mass(i + 1) * Glob_Mass(j + 1)) * EIGHT * Glob_Pi / THREE
      enddo
    enddo

    AnihMassChargeCoefficient = ZERO
    do i = 1, nn
      do j = 1, nn
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

      do k = 1, nn
        do l = 1, nn

          ketMatrix(k, l, i) = real(allPermutations(l, k, i))
          ! note the transposition here

        enddo
      enddo

    enddo

    do i = 1, nFactorial
      parityFactor(i) = real(parities(i))
    enddo

  end subroutine


  subroutine OverlapMatrixElement_CG_0S(Paramk, P, Skk)
!This subroutine computes symmetry adapted overlap matrix element
!of a complex L=0 correlated Gaussian
!
!fk =  exp[-r'(Lk*Lk'+i*Bk)r]
!
!with itself. Symmetry adaption is applied to the ket using
!a permutation matrix P.
!
!Input:
!   Paramk :: Array of length 2*(n(n+1)/2) of exponential parameters.
!             They are ordered as Paramk=(vechLk,vechBk)
!   P   :: The symmetry permutation matrix of size n x n
!Output:
!   Skk         ::        Overlap matrix element (not normalized)

!Arguments
    real(wp),intent(in)      :: Paramk(2*Glob_np)
    real(wp),intent(in)      :: P(Glob_n,Glob_n)
    complex(wp),intent(out)  :: Skk

!Parameters (These are needed to declare static arrays. Using static
!arrays makes the function call a little faster in comparison with
!the case when arrays are dynamically allocated in stack)
    integer,parameter :: nn=Glob_n

!Local variables
    real(wp)       Lk(nn,nn,1)
    complex(wp)    Ck(nn,nn,1), tCkl(nn,nn)
    complex(wp)    W1(nn,nn)
    complex(wp)    temp1
    complex(wp)    sqrtdet_tCkl
    integer           i, j, k

!First we build matrices Lk and Ck from Paramk
    call PrecomputeMatrices_L_C(2*Glob_np,1,Paramk,Lk,Ck)

!Then we permute elements of Ck to account for
!the action of the permutation matrix
!tCl=P'*Ck*P
!We also form matrix tCkl=conjg(Ck)+tCl
    do i=1,nn
    do j=1,nn
      temp1=ZERO
      do k=1,nn
        temp1=temp1+P(k,j)*Ck(k,i,1)
      enddo
      W1(j,i)=temp1
    enddo
    enddo
    do i=1,nn
    do j=i,nn
      temp1=ZERO
      do k=1,nn
        temp1=temp1+W1(i,k)*P(k,j)
      enddo
      tCkl(i,j)=conjg(Ck(i,j,1))+temp1
      tCkl(j,i)=tCkl(i,j)
    enddo
    enddo

!After this we can do pseudo-Cholesky factorization of tCkl.
!The pseudo-Cholesky factor will be temporarily stored in the
!lower triangle of W1. The product of its diagonal elements is
!the properly defined square root of det(tCkl) (see the comments
!in subroutine MatrixElementsHS_CG_0S).
    sqrtdet_tCkl=ONE
    do i=1,nn
    do j=i,nn
      temp1=tCkl(i,j)
      do k=i-1,1,-1
        temp1=temp1-W1(i,k)*W1(j,k)
      enddo
      if (i==j) then
        W1(i,i)=sqrt(temp1)
        sqrtdet_tCkl=sqrtdet_tCkl*W1(i,i)
      else
        W1(j,i)=temp1/W1(i,i)
        W1(i,j)=ZERO
      endif
    enddo
    enddo

!Evaluating overlap
    Skk=Glob_PiRaised3n2/(sqrtdet_tCkl*sqrtdet_tCkl*sqrtdet_tCkl)

  end subroutine OverlapMatrixElement_CG_0S

end module matelem
