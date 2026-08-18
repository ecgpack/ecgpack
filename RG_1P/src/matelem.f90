module matelem
!Module matelem contains subroutines for computing
!matrix elements with real L=1 Gaussians without normalization
  use globvars
  implicit none

contains

  subroutine MatrixElementsHS_RG_1P(m_k, m_l, vechLk, vechLl, P, &
                              Hkl, Skl, Dk, Dl, grad_k, grad_l)
!This subroutine computes symmetry adapted matrix element with
!two real L=1 correlated Gaussians:
!
!fk = z_{m_k} exp[-r'(Lk*Lk')r]
!
!where m_k is some integer between 1 and n (n is the number of
!pseudoparticles). Symmetry adaption is applied to the ket using
!permutation matrix P
!
!Input:
!   m_k, m_l :: integers that determine which z-component is in the
!                premultiplier of the Gaussian
!   vechLk, vechLl :: Arrays of length (n(n+1)/2) of
!     exponential parameters.
!   P  :: The symmetry permutation matrix of size n x n
!   grad_k, grad_l :: Gradient flags
!   grad_k=.true.  means that dHkldvechLk, dSkldvechLk need to be computed.
!   grad_l=.true.  means that dHkldvechLl, dSkldvechLl need to be computed.
!Output:
!   Hkl         ::        Hamiltonian term (normalized)
!   Skl         ::        Overlap matrix element (normalized)
!   Dk,Dl:: derivatives of normalized Hkl and Skl wrt Paramk
!           and Paraml respectively. They are ordered in the
!           following manner:
!           Dk=(dHkldvechLk,dSkldvechLk)
!           Dl=(dHkldvechLl,dSkldvechLl)

!Arguments
    integer,intent(in)          :: m_k,m_l
    real(wp),intent(in)      :: vechLk(Glob_np), vechLl(Glob_np)
    real(wp),intent(in)      :: P(Glob_n,Glob_n)
    real(wp),intent(out)     :: Skl,Hkl
    real(wp),intent(out)     :: Dk(2*Glob_np),Dl(2*Glob_np)
    logical,intent(in)          :: grad_k, grad_l

!Parameters (These are needed to declare static arrays. Using static
!arrays makes the function call a little faster in comparison with
!the case when arrays are dynamically allocated in stack)
    integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles
    integer        perm(nn)
    logical        Pisperm

!Local variables
    integer           n, np
    integer           tvl(nn)
    real(wp)       Lk(nn,nn),Ll(nn,nn)
    real(wp)       Ak(nn,nn),tAl(nn,nn),tAkl(nn,nn)
    real(wp)       inv_tAkl(nn,nn)
    real(wp)       inv_tAkltAl(nn,nn),inv_tAkltAlM(nn,nn)
    real(wp)       inv_tAklAk(nn,nn),inv_tAklAkM(nn,nn)
    real(wp)       eta1(nn,nn),sqrt_eta1(nn,nn),eta2(nn,nn),Rkl(nn,nn)
    real(wp)       W1(nn,nn),W3(nn,nn)
    real(wp)       twosym_tFkl(nn,nn),twosym_tGkl(nn,nn)
    real(wp)       tKkl(nn,nn),Ssym(nn,nn),Fkl(nn,nn)
    real(wp)       Cmat(nn,nn),Bmat(nn,nn),Zsym(nn,nn)
    real(wp)       inv_tAkltvl(nn),vkinv_tAkl(nn),vkinv_tAkltAlM(nn)
    real(wp)       u1(nn),u2(nn),u3(nn),dv1(nn),dv2(nn),pv(nn),rv(nn)
    real(wp)       temp1, temp2, temp3, temp4, temp5, temp6
    real(wp)       det_tAkl
    real(wp)       tau1,tau2,tau3,sigma,HklOverSkl
    real(wp)       Tkl, Vkl
    integer           i,j,k,indx

    n=Glob_n
    np=Glob_np
!First we build matrices Lk, Ll, Ak, Al from vechLk, vechLl.
    indx=0
    do i=1,nn
      do j=i,nn
        indx=indx+1
        Lk(i,j)=ZERO
        Lk(j,i)=vechLk(indx)
        Ll(i,j)=ZERO
        Ll(j,i)=vechLl(indx)
      enddo
    enddo

    do i=1,nn
      do j=i,nn
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
!Optimization: for every atomic symmetry projection P is a permutation
!matrix -- one +1 per column, zeros elsewhere -- and multiplying by one does
!not compute anything, it shuffles rows and columns:
!  (P'*A*P)(i,j) = A(perm(i),perm(j))
!so the two O(n^3) congruences below collapse into an O(n^2) gather. The
!scan that builds perm costs n^2 comparisons. If ANY column is not a single
!+1 (the -1 columns that occur in molecular/positronic systems) the original
!dense congruence is used instead, so the routine stays fully general.
    Pisperm=.true.
    do j=1,nn
      k=0
      do i=1,nn
        if (P(i,j)==ONE) then
          if (k/=0) Pisperm=.false.
          k=i
        elseif (P(i,j)/=ZERO) then
          Pisperm=.false.
        endif
      enddo
      if (k==0) then
        Pisperm=.false.
        k=j
      endif
      perm(j)=k
    enddo
    if (Pisperm) then
      !W1 is scratch: tAl cannot be permuted in place
      do i=1,nn
        do j=i,nn
          W1(i,j)=tAl(perm(i),perm(j))
        enddo
      enddo
      do i=1,nn
        do j=i,nn
          temp1=W1(i,j)
          tAl(i,j)=temp1
          tAl(j,i)=temp1
          tAkl(i,j)=Ak(i,j)+temp1
          tAkl(j,i)=tAkl(i,j)
        enddo
      enddo
    else
    do i=1,nn
      do j=1,nn
        temp1=ZERO
        do k=1,nn
          temp1=temp1+P(k,j)*tAl(k,i)
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
        tAl(i,j)=temp1
        tAl(j,i)=temp1
        tAkl(i,j)=Ak(i,j)+temp1
        tAkl(j,i)=tAkl(i,j)
      enddo
    enddo
    endif

!I delete non-necessary parts of the code

!The determinants of Lk and Ll are just
!the products of their diagonal elements
!det_Lk=ONE
!det_Ll=ONE
!do i=1,nn
!  det_Lk=det_Lk*Lk(i,i)
!  det_Ll=det_Ll*Ll(i,i)
!enddo

!After this we can do Cholesky factorization of tAkl.
!The Cholesky factor will be temporarily stored in the
!lower triangle of W1
    det_tAkl=ONE
    do i=1,nn
      do j=i,nn
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
        inv_tAkl(i,j)=temp1
        inv_tAkl(j,i)=temp1
      enddo
    enddo

!I delete non-necessary parts of the code

!Finding the inverse of Akk and All using their Cholesky factors
!The result is placed in inv_Akk and inv_All
!do i=1,nn
!  W1(i,i)=ONE/Lk(i,i)
!  W2(i,i)=ONE/Ll(i,i)
!  do j=i+1,nn
!    temp1=ZERO
!    temp2=ZERO
!    do k=i,j-1
!      temp1=temp1-Lk(j,k)*W1(k,i)
!      temp2=temp2-Ll(j,k)*W2(k,i)
!    enddo
!    W1(j,i)=temp1/Lk(j,j)
!    W2(j,i)=temp2/Ll(j,j)
!  enddo
!enddo

!do i=1,nn
!  do j=i,nn
!     temp1=ZERO
!     temp2=ZERO
!    do k=j,nn
!       temp1=temp1+W1(k,i)*W1(k,j)
!       temp2=temp2+W2(k,i)*W2(k,j)
!     enddo
!     inv_Akk(i,j)=ONEHALF*temp1
!         inv_Akk(j,i)=ONEHALF*temp1
!     inv_All(i,j)=ONEHALF*temp2
!         inv_All(j,i)=ONEHALF*temp2
!   enddo
!enddo

!Computing tvl=P'*vl
    do i=1,nn
      tvl(i)=P(m_l,i)
    enddo

!Compute inv_tAkltvl = inv_tAkl * tvl
    do i=1,nn
      temp1=ZERO
      do j=1,nn
        temp1=temp1+inv_tAkl(j,i)*tvl(j)
      enddo
      inv_tAkltvl(i)=temp1
    enddo

!Compute vkinv_tAkl=vk'*inv_tAkl
    do i=1,nn
      vkinv_tAkl(i)=inv_tAkl(m_k,i)
    enddo

!Compute tau3=vkinv_tAkl*tvl
    tau3=ZERO
    do i=1,nn
      tau3=tau3+vkinv_tAkl(i)*tvl(i)
    enddo

!Evaluating overlap
!temp1=abs(det_Ll*det_Lk)/det_tAkl
    temp1=det_tAkl*sqrt(det_tAkl)
!Skl=Glob_2Raised3n2*tau3*temp1*sqrt(temp1/(inv_Akk(m_k,m_k)*inv_All(m_l,m_l)))
    Skl=Glob_PiRaised3n2*tau3/(TWO*temp1)

!Doing multiplication inv_tAkltAl=inv_tAkl*tAl
    do i=1,nn
      do j=1,nn
        temp1=ZERO
        do k=1,nn
          temp1=temp1+inv_tAkl(j,k)*tAl(k,i)
        enddo
        inv_tAkltAl(j,i)=temp1
      enddo
    enddo

!Doing multiplication inv_tAkltAlM=inv_tAkltAl*M
    do i=1,nn
      do j=1,nn
        temp1=ZERO
        do k=1,nn
          temp1=temp1+inv_tAkltAl(j,k)*Glob_MassMatrix(k,i)
        enddo
        inv_tAkltAlM(j,i)=temp1
      enddo
    enddo

!Computing tau1=tr[inv_tAkltAlM*Ak]
    tau1=ZERO
    do i=1,nn
      temp1=ZERO
      do k=1,nn
        temp1=temp1+inv_tAkltAlM(i,k)*Ak(k,i)
      enddo
      tau1=tau1+temp1
    enddo
!Computing tau2 = vk'*inv_tAkltAlM*Ak*inv_tAkltvl
!We do it by multiplying twice the row-vector on the left
!by a matrix on the right and computing a dot product in the end.

!vkinv_tAkltAlM'=vk'*inv_tAkltAlM
    do i=1,nn
      vkinv_tAkltAlM(i)=inv_tAkltAlM(m_k,i)
    enddo
!u1=vkinv_tAkltAlM'*Ak
!tau2=u1'*inv_tAkltvl (storage for u1 as such is not needed, we use temp1=u1(i))
    tau2=ZERO
    do i=1,nn
      temp1=ZERO
      !nvfortran 25.9/26.3 miscompile this nest when BOTH bounds are the
      !compile-time nn: at nparticles=5, -O2/-O3, Hkl comes back ~50% wrong
      !with no warning. Either bound as the runtime n avoids it, and n == nn
      !always. Measured cost ~0.5% (gfortran, median over 42 cells).
      do j=1,n
        temp1=temp1+vkinv_tAkltAlM(j)*Ak(j,i)
      enddo
      tau2=tau2+temp1*inv_tAkltvl(i)
    enddo

!Evaluating the kinetic energy
    Tkl=Skl*(SIX*tau1+FOUR*tau2/tau3)

!Evaluating eta1(i,j), sqrt_eta1(i,j), eta2(i,j), Rkl(i,j),
!and the potential energy. Notice that only the lower triangles
!of eta1, sqrt_eta1, eta2, and Rkl are filled.
    Vkl=ZERO
    temp1=Skl*(TWO/Glob_SqrtPi)
    do i=1,nn
      temp2=inv_tAkl(i,i)
      temp3=sqrt(temp2)
      eta1(i,i)=temp2
      sqrt_eta1(i,i)=temp3
      !eta2(i,i)=vk'*inv_tAkl*Jii*inv_tAkl*tvl. Since Jii=e_i*e_i',
      !the matrix inv_tAkl*Jii*inv_tAkl is rank one and eta2 is just
      !the product of two already available numbers.
      temp4=vkinv_tAkl(i)*inv_tAkltvl(i)
      eta2(i,i)=temp4
      Rkl(i,i)=temp1*(ONE-temp4/(THREE*temp2*tau3))/temp3
      Vkl=Vkl+Glob_ScaledPseudoChargeMatrix(i,0)*Rkl(i,i)
    enddo
    do i=1,nn
      do j=i+1,nn
        temp2=inv_tAkl(i,i)+inv_tAkl(j,j)-inv_tAkl(j,i)-inv_tAkl(j,i)
        temp3=sqrt(temp2)
        eta1(j,i)=temp2
        sqrt_eta1(j,i)=temp3
        !eta2(j,i)=vk'*inv_tAkl*Jij*inv_tAkl*tvl. Since
        !Jij=(e_i-e_j)*(e_i-e_j)', the matrix inv_tAkl*Jij*inv_tAkl
        !is rank one and eta2 is just the product of two already
        !available vector-element differences.
        temp4=(vkinv_tAkl(i)-vkinv_tAkl(j))*(inv_tAkltvl(i)-inv_tAkltvl(j))
        eta2(j,i)=temp4
        Rkl(j,i)=temp1*(ONE-temp4/(THREE*temp2*tau3))/temp3
        Vkl=Vkl+Glob_ScaledPseudoChargeMatrix(i,j)*Rkl(j,i)
      enddo
    enddo

    Hkl=Tkl+Vkl

!Now we compute the gradients. This part is an optimized rewrite of
!the original code (results identical up to roundoff). Notation:
!  u = inv_tAkltvl = inv_tAkl*tvl,  v = vkinv_tAkl' = inv_tAkl*vk,
!  W2 = inv_tAkltAl = inv_tAkl*tAl, K = inv_tAkltAlM = W2*M,
!  Ssym = tKkl+tKkl' = u*v'+v*u'.
!
!Gradient of Vkl. The original code looped over all particle pairs
!and for each pair built the n x n matrix
!  twosym_tQkl = c1*( c2*a*a' + (1/(3*tau3))*(a*s'+s*a')
!                     -(eta2/(3*tau3^2))*Ssym ),
!  a = inv_tAkl*(e_i-e_j),  s = Ssym*(e_i-e_j),
!then multiplied it by Lk (and, for the ket, sandwiched it between P
!and P' and multiplied by Ll) - O(n^5) overall. The weighted sum of
!all pair terms can instead be assembled analytically in O(n^3):
!  Bmat = inv_tAkl*Cmat*inv_tAkl               (a*a' terms; Cmat
!         accumulates the pair weights like a charge matrix)
!       + pv*u'+u*pv'+rv*v'+v*rv'              (a*s'+s*a' terms,
!         because s = beta*u+gamma*v with per-pair scalars beta,gamma,
!         so pv=inv_tAkl*dv1 and rv=inv_tAkl*dv2 with the vectors
!         dv1,dv2 accumulating the scalar weights)
!       - sigma*Ssym                           (last terms)
!
!Gradient of Tkl. It needs the matrices
!  6*K*W2'                          = 6*Fkl              (bra)
!  6*(inv_tAkl*Ak*M)*(inv_tAkl*Ak)' = 6*(M-K-K'+Fkl)     (ket)
!where the second equality holds because inv_tAkl*Ak = I-W2
!(as Ak+tAl=tAkl), so the ket matrix costs O(n^2) once Fkl is known.
!
!The T and V contributions are combined into a single symmetric
!matrix Zsym so that only one triangular product with Lk (resp. one
!congruence with P and one triangular product with Ll) is needed:
!  dHkl/dvechLk = (Hkl/Skl)*dSkl/dvechLk + Skl*vech[Zsym*Lk]
!  dHkl/dvechLl = (Hkl/Skl)*dSkl/dvechLl + Skl*vech[(P*Zsym*P')*Ll]

    if (grad_k.or.grad_l) then
      !Evaluating matrix tKkl = inv_tAkltvl * vkinv_tAkl'
      !which will be used a lot below, and its symmetrization Ssym
      do i=1,nn
        do j=1,nn
          tKkl(i,j)=inv_tAkltvl(i)*vkinv_tAkl(j)
        enddo
      enddo
      do i=1,nn
        do j=1,i
          Ssym(i,j)=tKkl(i,j)+tKkl(j,i)
          Ssym(j,i)=Ssym(i,j)
        enddo
      enddo
      !Evaluating twosym_tFkl = tFkl + tFkl' , where tFkl = (3/2)*inv_tAkl + tKkl/tau3
      do i=1,nn
        do j=1,i
          twosym_tFkl(i,j)=THREE*inv_tAkl(j,i)+Ssym(i,j)/tau3
          twosym_tFkl(j,i)=twosym_tFkl(i,j)
        enddo
      enddo
      !Evaluating Fkl=inv_tAkltAlM*inv_tAkltAl'
      !(only the upper triangle, then mirrored)
      do j=1,nn
        do i=1,j
          temp1=ZERO
          do k=1,nn
            temp1=temp1+inv_tAkltAlM(i,k)*inv_tAkltAl(j,k)
          enddo
          Fkl(i,j)=temp1
          Fkl(j,i)=temp1
        enddo
      enddo
      HklOverSkl=Hkl/Skl
    endif

!Gradient of Skl

    if (grad_k) then
      !Evaluating -Skl*vech((twosym_tFkl)*Lk)'
      indx=0
      do i=1,nn
        do j=i,nn
          temp1=ZERO
          do k=i,nn
            temp1=temp1-twosym_tFkl(k,j)*Lk(k,i)
          enddo
          indx=indx+1
          Dk(Glob_np+indx)=Skl*temp1
        enddo
      enddo
    endif

    if (grad_l) then
      !Evaluating twosym_tGkl = P * twosym_tFkl *P'
      do i=1,nn
        do j=1,nn
          temp1=ZERO
          do k=1,nn
            temp1=temp1+P(i,k)*twosym_tFkl(k,j)
          enddo
          W1(i,j)=temp1
        enddo
      enddo
      do i=1,nn
        do j=1,i
          temp1=ZERO
          do k=1,nn
            temp1=temp1+W1(i,k)*P(j,k)
          enddo
          twosym_tGkl(i,j)=temp1
          twosym_tGkl(j,i)=temp1
        enddo
      enddo
      !Evaluating -Skl*vech((twosym_tGkl)*Ll)'
      indx=0
      do i=1,nn
        do j=i,nn
          temp1=ZERO
          do k=i,nn
            temp1=temp1-twosym_tGkl(k,j)*Ll(k,i)
          enddo
          indx=indx+1
          Dl(Glob_np+indx)=Skl*temp1
        enddo
      enddo
    endif

!Assembling the pair sum Bmat for the gradient of Vkl
!(see the comment above)

    if (grad_k.or.grad_l) then
      do j=1,nn
        do i=1,nn
          Cmat(i,j)=ZERO
        enddo
        dv1(j)=ZERO
        dv2(j)=ZERO
      enddo
      sigma=ZERO
      temp3=ONETHIRD/tau3
      !terms with Jii (interaction with the reference particle)
      do i=1,nn
        temp1=(TWO/Glob_SqrtPi)/(eta1(i,i)*sqrt_eta1(i,i))
        temp5=Glob_ScaledPseudoChargeMatrix(i,0)*temp1
        temp2=ONE-eta2(i,i)/(eta1(i,i)*tau3)
        Cmat(i,i)=Cmat(i,i)+temp5*temp2
        temp6=temp5*temp3
        dv1(i)=dv1(i)+temp6*vkinv_tAkl(i)
        dv2(i)=dv2(i)+temp6*inv_tAkltvl(i)
        sigma=sigma+temp6*eta2(i,i)/tau3
      enddo
      !terms with Jij (interparticle interactions)
      do i=1,nn
        do j=i+1,nn
          temp1=(TWO/Glob_SqrtPi)/(eta1(j,i)*sqrt_eta1(j,i))
          temp5=Glob_ScaledPseudoChargeMatrix(i,j)*temp1
          temp2=ONE-eta2(j,i)/(eta1(j,i)*tau3)
          temp4=temp5*temp2
          Cmat(i,i)=Cmat(i,i)+temp4
          Cmat(j,j)=Cmat(j,j)+temp4
          Cmat(j,i)=Cmat(j,i)-temp4
          Cmat(i,j)=Cmat(i,j)-temp4
          temp6=temp5*temp3
          temp4=temp6*(vkinv_tAkl(i)-vkinv_tAkl(j))
          dv1(i)=dv1(i)+temp4
          dv1(j)=dv1(j)-temp4
          temp4=temp6*(inv_tAkltvl(i)-inv_tAkltvl(j))
          dv2(i)=dv2(i)+temp4
          dv2(j)=dv2(j)-temp4
          sigma=sigma+temp6*eta2(j,i)/tau3
        enddo
      enddo
      !W1=Cmat*inv_tAkl (both factors are symmetric)
      do j=1,nn
        do i=1,nn
          temp1=ZERO
          do k=1,nn
            temp1=temp1+Cmat(k,i)*inv_tAkl(k,j)
          enddo
          W1(i,j)=temp1
        enddo
      enddo
      !upper triangle of inv_tAkl*Cmat*inv_tAkl=W1'*inv_tAkl, mirrored
      do j=1,nn
        do i=1,j
          temp1=ZERO
          do k=1,nn
            temp1=temp1+W1(k,i)*inv_tAkl(k,j)
          enddo
          Bmat(i,j)=temp1
          Bmat(j,i)=temp1
        enddo
      enddo
      !pv=inv_tAkl*dv1 and rv=inv_tAkl*dv2
      do i=1,nn
        temp1=ZERO
        temp2=ZERO
        do k=1,nn
          temp1=temp1+inv_tAkl(k,i)*dv1(k)
          temp2=temp2+inv_tAkl(k,i)*dv2(k)
        enddo
        pv(i)=temp1
        rv(i)=temp2
      enddo
      !adding the rank-one and Ssym parts
      do j=1,nn
        do i=1,nn
          Bmat(i,j)=Bmat(i,j)+pv(i)*inv_tAkltvl(j)+inv_tAkltvl(i)*pv(j) &
                    +rv(i)*vkinv_tAkl(j)+vkinv_tAkl(i)*rv(j)-sigma*Ssym(i,j)
        enddo
      enddo
    endif

!Gradient of Hkl with respect to vechLk

    if (grad_k) then
      !Computing u1'=vkinv_tAkltAlM'*inv_tAkltAl'
      do i=1,nn
        temp1=ZERO
        do j=1,nn
          temp1=temp1+vkinv_tAkltAlM(j)*inv_tAkltAl(i,j)
        enddo
        u1(i)=temp1
      enddo
      !Computing u2=inv_tAkltAlM*Ak*inv_tAkltvl
      do i=1,nn
        temp1=ZERO
        do j=1,nn
          temp1=temp1+Ak(i,j)*inv_tAkltvl(j)
        enddo
        u3(i)=temp1
      enddo
      !nvfortran 25.9 miscompiles RG_1P's gradient path when this and three
      !neighbouring nests are all constant-bound (nparticles>=7, -cuda builds
      !only). Any one of the four as the runtime n avoids it.
      do i=1,n
        temp1=ZERO
        do j=1,nn
          temp1=temp1+inv_tAkltAlM(i,j)*u3(j)
        enddo
        u2(i)=temp1
      enddo
      !Zsym = tUkl+tUkl'+Bmat, where
      !tUkl = 6*Fkl + (4/tau3)*(inv_tAkltvl*u1'-u2*vkinv_tAkl') + (4*tau2/tau3^2)*tKkl
      temp1=FOUR/tau3
      temp2=temp1*tau2/tau3
      do j=1,nn
        do i=1,nn
          Zsym(i,j)=12*Fkl(i,j) &
                    +temp1*(inv_tAkltvl(i)*u1(j)+u1(i)*inv_tAkltvl(j) &
                            -u2(i)*vkinv_tAkl(j)-vkinv_tAkl(i)*u2(j)) &
                    +temp2*Ssym(i,j)+Bmat(i,j)
        enddo
      enddo
      !Evaluating (Hkl/Skl)*dSkldvechLk' + Skl*vech(Zsym*Lk)'
      indx=0
      do i=1,nn
        do j=i,nn
          temp1=ZERO
          do k=i,nn
            temp1=temp1+Zsym(k,j)*Lk(k,i)
          enddo
          indx=indx+1
          Dk(indx)=Skl*temp1+HklOverSkl*Dk(Glob_np+indx)
        enddo
      enddo
    endif

!Gradient of Hkl with respect to vechLl

    if (grad_l) then
      !Computing inv_tAklAk=inv_tAkl*Ak and inv_tAklAkM=inv_tAkl*Ak*M.
      !Both come for free: inv_tAkl*Ak=I-inv_tAkltAl (as Ak+tAl=tAkl),
      !hence inv_tAkl*Ak*M=M-inv_tAkltAlM.
      do i=1,nn
        do j=1,nn
          inv_tAklAk(j,i)=-inv_tAkltAl(j,i)
          inv_tAklAkM(j,i)=Glob_MassMatrix(j,i)-inv_tAkltAlM(j,i)
        enddo
        inv_tAklAk(i,i)=inv_tAklAk(i,i)+ONE
      enddo
      !Computing u1=inv_tAklAkM*inv_tAklAk'*tvl
      do i=1,nn
        temp1=ZERO
        do j=1,nn
          temp1=temp1+inv_tAklAk(j,i)*tvl(j)
        enddo
        u3(i)=temp1
      enddo
      do i=1,nn
        temp1=ZERO
        do j=1,nn
          temp1=temp1+inv_tAklAkM(i,j)*u3(j)
        enddo
        u1(i)=temp1
      enddo
      !Computing u2'=vkinv_tAkltAlM'*inv_tAklAk'
      do i=1,nn
        temp1=ZERO
        do j=1,nn
          temp1=temp1+vkinv_tAkltAlM(j)*inv_tAklAk(i,j)
        enddo
        u2(i)=temp1
      enddo
      !Zsym = tWkl+tWkl'+Bmat (the inner matrix, before the P congruence), where
      !tWkl = 6*inv_tAklAkM*inv_tAklAk' + (4/tau3)*(u1*vkinv_tAkl'-inv_tAkltvl*u2')
      !       + (4*tau2/tau3^2)*tKkl
      !and 6*inv_tAklAkM*inv_tAklAk' = 6*(M-K-K'+Fkl) (see the comment above)
      temp1=FOUR/tau3
      temp2=temp1*tau2/tau3
      do j=1,nn
        do i=1,nn
          Zsym(i,j)=12*(Glob_MassMatrix(i,j)-inv_tAkltAlM(i,j) &
                        -inv_tAkltAlM(j,i)+Fkl(i,j)) &
                    +temp1*(u1(i)*vkinv_tAkl(j)+vkinv_tAkl(i)*u1(j) &
                            -inv_tAkltvl(i)*u2(j)-u2(i)*inv_tAkltvl(j)) &
                    +temp2*Ssym(i,j)+Bmat(i,j)
        enddo
      enddo
      !Congruence W3=P*Zsym*P' (only the lower triangle, then mirrored)
      do i=1,nn
        do j=1,nn
          temp1=ZERO
          do k=1,nn
            temp1=temp1+P(i,k)*Zsym(k,j)
          enddo
          W1(i,j)=temp1
        enddo
      enddo
      do i=1,nn
        do j=1,i
          temp1=ZERO
          do k=1,nn
            temp1=temp1+W1(i,k)*P(j,k)
          enddo
          W3(i,j)=temp1
          W3(j,i)=temp1
        enddo
      enddo
      !Evaluating (Hkl/Skl)*dSkldvechLl' + Skl*vech(W3*Ll)'
      indx=0
      do i=1,nn
        do j=i,nn
          temp1=ZERO
          do k=i,nn
            temp1=temp1+W3(k,j)*Ll(k,i)
          enddo
          indx=indx+1
          Dl(indx)=Skl*temp1+HklOverSkl*Dl(Glob_np+indx)
        enddo
      enddo
    endif

  end subroutine MatrixElementsHS_RG_1P

  subroutine MatrixElementsAll_RG_1P(m_k, m_l, vechLk, vechLl, Pbra, Pket, &
                                         Hkl, Skl, Tkl, Vkl, rm2kl, rmkl, rkl, r2kl, deltarkl, drach_deltarkl, &
                                         MVkl, drach_MVkl1, drach_MVkl2, Darwinkl, drach_Darwinkl, OOkl, rmrmkl, prvalkl, &
                                         NumCFGridPoints, CFGrid, CFkl, NumDensGridPoints, DensGrid, Denskl, &
                                         AreCorrFuncNeeded, ArePartDensNeeded, AreMCorrFuncNeeded, AreMPartDensNeeded)
!This subroutine computes symmetry adapted matrix elements
!with two real L=1 correlated Gaussians. These matrix elements
!are used in calculations of expectation values.
!Symmetry adaptation is applied to the bra and ket using permutation matrices Pbra and Pket
!
!Input:
!   m_k, m_l :: integers that determine which z-component is in the
!       premultiplier of the Gaussian
!   vechLk, vechLl :: Arrays of length (n(n+1)/2) of exponential parameters.
!   Pbra :: The symmetry permutation matrix of size n x n that is applied to bra
!   Pket :: The symmetry permutation matrix of size n x n that is applied to ket
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
! prvalkl   :: P(1/r^3_ij) - principal values of matrix element 1/r^3_ij  (appears in the Araki-Sucker term for QED correction)
!NumCFGridPoints    :: Number of grid points for correlation function calculations
!CFGrid             :: Array containing grid points where matrix elements of
!                      correlation functions should be computed
!CFkl               :: Matrix elements of correlation functions
!NumDensGridPoints  :: Number of grid points for particle density calculations
!DensGrid           :: Array containing grid points where matrix elements of
!                      particle densities should be computed
!Denskl             :: Matrix elements of particle densities
!AreCorrFuncNeeded  :: flag indicating whether matrix elements of correlation
!                      functions need to be computed
!ArePartDensNeeded  :: flag indicating whether matrix elements of particle
!                      densities need to be computed
!AreMCorrFuncNeeded :: flag indicating whether matrix elements of momentum correlation
!                      functions need to be computed
!AreMPartDensNeeded :: flag indicating whether matrix elements of particle
!                      momentum densities need to be computed

!Arguments
    integer,intent(in)       :: m_k,m_l
    real(wp),intent(in)   :: vechLk(Glob_np), vechLl(Glob_np)
    real(wp),intent(in)   :: Pbra(Glob_n,Glob_n),Pket(Glob_n,Glob_n)
    real(wp),intent(out)  :: Hkl,Skl,Tkl,Vkl,MVkl,drach_MVkl1,drach_MVkl2,Darwinkl,drach_Darwinkl,OOkl
    real(wp),intent(out)  :: rm2kl(Glob_n,Glob_n),rmkl(Glob_n,Glob_n)
    real(wp),intent(out)  :: rkl(Glob_n,Glob_n),r2kl(Glob_n,Glob_n)
    real(wp),intent(out)  :: deltarkl(Glob_n,Glob_n)
    real(wp),intent(out)  :: drach_deltarkl(Glob_n,Glob_n)
    real(wp),intent(out)  :: prvalkl(Glob_n,Glob_n)
    real(wp),intent(out)  :: rmrmkl(Glob_n,Glob_n,Glob_n,Glob_n)
    integer,intent(in)       :: NumCFGridPoints,NumDensGridPoints
    real(wp),intent(in)   :: CFGrid(2,NumCFGridPoints),DensGrid(2,NumDensGridPoints)
    real(wp),intent(out)  :: CFkl(Glob_n*(Glob_n+1)/2,NumCFGridPoints)
    real(wp),intent(out)  :: Denskl(Glob_n+1,NumDensGridPoints)
    logical,intent(in)       :: AreCorrFuncNeeded,ArePartDensNeeded,AreMCorrFuncNeeded,AreMPartDensNeeded

!Parameters (These are needed to declare static arrays. Using static
!arrays makes the function call a little faster in comparison with
!the case when arrays are dynamically allocated in stack)
    integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles
    integer,parameter :: nnp=nn*(nn+1)/2

!Local variables
    integer           n,np
    integer           tvk(nn),tvl(nn)
    real(wp)       Lk(nn,nn),Ll(nn,nn)
    real(wp)       inv_tAk(nn,nn),inv_tAl(nn,nn),tAk(nn,nn),tAl(nn,nn),tAkl(nn,nn)
    real(wp)       inv_tAkl(nn,nn), inv_tAkltAl(nn,nn)
    real(wp)       inv_invtAkinvtAl(nn,nn),tvkinv_tAk(nn),inv_tAltvl(nn)
    real(wp)       eta2(nn,nn),meta2(nn,nn),inv_tAkltAlM(nn,nn)
    real(wp)       W1(nn,nn),W2(nn,nn)
    real(wp)       inv_tAkltvl(nn),tvkinv_tAkl(nn),tvkinv_tAkltAlM(nn),u1(nn)
    real(wp)       tvkinv_tAkinv_invtAkinvtAl(nn),inv_invtAkinvtAlinv_tAltvl(nn)
    real(wp)       temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,temp9
    real(wp)       temp10,temp11,temp12,temp13,temp14,threshold,tr1, tr2, tr3, tr4
    real(wp)       det_tAkl, det_tAk, det_tAl, det_invtAkinvtAl
    real(wp)       tau1,tau2,tau3,inv_tau3 ,V2kl, tau4, MSkl
    integer           i,j,k,t,indx,p,q
    real(wp)       TrAJ(nn,nn),sqrtTrAJ(nn,nn),TrAJAJ(nn,nn,nn,nn),MTrAJ(nn,nn),sqrtMTrAJ(nn,nn)
    real(wp)       jAj(nn,nn,nn,nn),jAtvl(nn,nn),tvkAj(nn,nn),Mass_For_Darwin(0:nn)
!Precomputed matrices/vectors for the fast (rank-one based) evaluation
!of the orbit-orbit, mass-velocity, and drachmanized delta sections
    real(wp)       GAl(nn,nn),tAlAvk(nn),tAlAvl(nn)
    real(wp)       UXd(nn,nn),Qsd(nn,nn),gvo(nn),gvc(nn),gv3(nn),gv4(nn),wvd(nn)
    real(wp)       dtheta,dlam,dkap,dom,dch,dhh,dmm,dt1,dt2
    real(wp)       sk(nn),sl(nn),Ask(nn),Asl(nn)
    real(wp)       qkk,qll,qkl,vkqk,vkql,qkvl,qlvl,wvk,wvl,wqk,wql
    real(wp)       trAvltuk,trAtulvk,trAultuk,trKx,trLx,merr,dwd2
    real(wp)       s_ab,s_ad,s_ak,s_al,s_bc,s_bd,s_bk,s_bl,s_cd,s_ck,s_cl,s_dd,s_dk,s_dl
    real(wp)       tX1,tXJ1,tXV1,kXd1,dXl1,tX2,tXJ2,tXV2,kXd2,dXl2
    real(wp)       tY3,tYJ3,tYV3,kYd3,dYl3,tX4,tXJ4,tY5,tYJ5,tY7,tYJ7
    real(wp)       tXY23,tXYJ23,tXYV23,tYXV23,kXYd23,kYXd23,dXYl23,dYXl23
    real(wp)       tXY35,tXYJ35,tXY27,tXYJ27,tJVs,oc,odd,vkj,vli,vlj

    n=Glob_n
    np=Glob_np
!First we build matrices Lk, Ll, Ak, Al from vechLk, vechLl.
    indx=0
    do i=1,nn
      do j=i,nn
        indx=indx+1
        Lk(i,j)=ZERO
        Lk(j,i)=vechLk(indx)
        Ll(i,j)=ZERO
        Ll(j,i)=vechLl(indx)
      enddo
    enddo

    do i=1,nn
      do j=i,nn
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
    do i=1,nn
      do j=1,nn
        temp1=ZERO
        temp2=ZERO
        do k=1,nn
          temp1=temp1+Pket(k,j)*tAl(k,i)
          temp2=temp2+tAk(j,k)*Pbra(k,i)
        enddo
        W1(j,i)=temp1
        W2(j,i)=temp2
      enddo
    enddo
!tAl=W1*Pket
!tAk=Pbra'*W2
    do i=1,nn
      do j=i,nn
        temp1=ZERO
        temp2=ZERO
        do k=1,nn
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

!I delete non-necessary parts of the code
!The determinants of Lk and Ll are just
!the products of their diagonal elements
!det_Lk=ONE
!det_Ll=ONE
!do i=1,nn
!  det_Lk=det_Lk*Lk(i,i)
!  det_Ll=det_Ll*Ll(i,i)
!enddo

!After this we can do Cholesky factorization of tAkl.
!The Cholesky factor will be temporarily stored in the
!lower triangle of W1
    det_tAkl=ONE
    do i=1,nn
      do j=i,nn
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
        inv_tAkl(i,j)=temp1
        inv_tAkl(j,i)=temp1
      enddo
    enddo

!Do Cholesky factorization of tAk and tAk and then invert
    if (AreMCorrFuncNeeded.or.AreMPartDensNeeded) then
      det_tAk=ONE
      det_tAl=ONE
      do i=1,nn
        do j=i,nn
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
      do i=1,nn
        W1(i,i)=ONE/W1(i,i)
        W2(i,i)=ONE/W2(i,i)
        do j=i+1,nn
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

      do i=1,nn
        do j=i,nn
          temp1=ZERO
          temp2=ZERO
          do k=j,nn
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
      do i=1,nn
        do j=i,nn
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
          inv_invtAkinvtAl(i,j)=temp1
          inv_invtAkinvtAl(j,i)=temp1
        enddo
      enddo
    endif

!I delete non-necessary parts of the code
!Finding the inverse of Akk and All using their Cholesky factors
!The result is placed in inv_Akk and inv_All
!do i=1,nn
!  W1(i,i)=ONE/Lk(i,i)
!  W2(i,i)=ONE/Ll(i,i)
!  do j=i+1,nn
!    temp1=ZERO
!    temp2=ZERO
!    do k=i,j-1
!      temp1=temp1-Lk(j,k)*W1(k,i)
!      temp2=temp2-Ll(j,k)*W2(k,i)
!    enddo
!    W1(j,i)=temp1/Lk(j,j)
!    W2(j,i)=temp2/Ll(j,j)
!  enddo
!enddo

!do i=1,nn
!  do j=i,nn
!     temp1=ZERO
!     temp2=ZERO
!     do k=j,nn
!       temp1=temp1+W1(k,i)*W1(k,j)
!       temp2=temp2+W2(k,i)*W2(k,j)
!     enddo
!     inv_Akk(i,j)=ONEHALF*temp1
!         inv_Akk(j,i)=ONEHALF*temp1
!     inv_All(i,j)=ONEHALF*temp2
!         inv_All(j,i)=ONEHALF*temp2
!   enddo
!enddo

!Computing tvk=Pbra'*vk and tvl=Pket'*vl
    do i=1,nn
      tvk(i)=Pbra(m_k,i)
      tvl(i)=Pket(m_l,i)
    enddo

!Compute inv_tAkltvl = inv_tAkl * tvl
    do i=1,nn
      temp1=ZERO
      do j=1,nn
        temp1=temp1+inv_tAkl(j,i)*tvl(j)
      enddo
      inv_tAkltvl(i)=temp1
    enddo

!Compute tvkinv_tAkl=tvk'*inv_tAkl
    do i=1,nn
      temp1=ZERO
      do j=1,nn
        temp1=temp1+tvk(j)*inv_tAkl(j,i)
      enddo
      tvkinv_tAkl(i)=temp1
    enddo

!Compute tau3=tvkinv_tAkl*tvl
    tau3=ZERO
    do i=1,nn
      tau3=tau3+tvkinv_tAkl(i)*tvl(i)
    enddo

    if (AreMCorrFuncNeeded.or.AreMPartDensNeeded)  then
      !Compute tvkinv_tAk
      do i=1,nn
        temp1=ZERO
        do j=1,nn
          temp1=temp1+tvk(j)*inv_tAk(j,i)
        enddo
        tvkinv_tAk(i)=temp1
      enddo

      !Compute inv_tAltvl
      do i=1,nn
        temp1=ZERO
        do j=1,nn
          temp1=temp1+inv_tAl(j,i)*tvl(j)
        enddo
        inv_tAltvl(i)=temp1
      enddo

      !Compute inv_invtAkinvtAlinv_tAltvl
      do i=1,nn
        temp1=ZERO
        do j=1,nn
          temp1=temp1+inv_invtAkinvtAl(j,i)*inv_tAltvl(j)
        enddo
        inv_invtAkinvtAlinv_tAltvl(i)=temp1
      enddo

      !Compute tvkinv_tAkinv_invtAkinvtAl
      do i=1,nn
        temp1=ZERO
        do j=1,nn
          temp1=temp1+inv_invtAkinvtAl(j,i)*tvkinv_tAk(j)
        enddo
        tvkinv_tAkinv_invtAkinvtAl(i)=temp1
      enddo

      !Compute tau4 = tvkinv_tAk*inv_invtAkinvtAl*inv_tAltvl
      tau4=ZERO
      do i=1,nn
        do j=1,nn
          tau4=tau4+tvkinv_tAk(i)*inv_invtAkinvtAl(i,j)*inv_tAltvl(j)
        enddo
      enddo
    end if

!Evaluating overlap
!temp1=abs(det_Ll*det_Lk)/det_tAkl
!Skl=Glob_2Raised3n2*tau3*temp1*sqrt(temp1/(inv_Akk(m_k,m_k)*inv_All(m_l,m_l)))
    temp1=det_tAkl*sqrt(det_tAkl)
    Skl=Glob_PiRaised3n2*tau3/(TWO*temp1)

    if(AreMCorrFuncNeeded.or.AreMPartDensNeeded) then
      temp1=det_tAk*det_tAl*det_invtAkinvtAl
      temp2=temp1*sqrt(temp1)
      MSkl=Glob_PiRaised3n2*tau4/(TWO*temp2)
    endif

!Doing multiplication inv_tAkltAl=inv_tAkl*tAl, inv_tAkltAk=inv_tAkl*tAk
    do i=1,nn
      do j=1,nn
        temp1=ZERO
        do k=1,nn
          temp1=temp1+inv_tAkl(j,k)*tAl(k,i)
        enddo
        inv_tAkltAl(j,i)=temp1
      enddo
    enddo

!Doing multiplication inv_tAkltAlM=inv_tAkltAl*M
    do i=1,nn
      do j=1,nn
        temp1=ZERO
        do k=1,nn
          temp1=temp1+inv_tAkltAl(j,k)*Glob_MassMatrix(k,i)
        enddo
        inv_tAkltAlM(j,i)=temp1
      enddo
    enddo

!GAl = tAl*inv_tAkl*tAl = tAl*inv_tAkltAl (symmetric), and the vectors
!tAlAvk = tAl*inv_tAkl*tvk, tAlAvl = tAl*inv_tAkl*tvl. Their elements
!supply the O(1) building blocks of the rank-one based orbit-orbit and
!mass-velocity sections below.
    do i=1,nn
      do j=i,nn
        temp1=ZERO
        do k=1,nn
          temp1=temp1+tAl(j,k)*inv_tAkltAl(k,i)
        enddo
        GAl(j,i)=temp1
        GAl(i,j)=temp1
      enddo
    enddo
    do i=1,nn
      temp1=ZERO
      temp2=ZERO
      do k=1,nn
        temp1=temp1+tAl(i,k)*tvkinv_tAkl(k)
        temp2=temp2+tAl(i,k)*inv_tAkltvl(k)
      enddo
      tAlAvk(i)=temp1
      tAlAvl(i)=temp2
    enddo

!Evaluating vector-matrix-vector products
!j^{ij}' inv_tAkl tvl
!tvk' inv_tAkl j^{ij}
!(moved up: the eta2 evaluation below uses them)
    do j=1,nn
      do i=1,nn
        if (i==j) then
          jAtvl(i,i)=inv_tAkltvl(i)
          tvkAj(i,i)=tvkinv_tAkl(i)
        else
          jAtvl(i,j)=inv_tAkltvl(i)-inv_tAkltvl(j)
          tvkAj(i,j)=tvkinv_tAkl(i)-tvkinv_tAkl(j)
        endif
      enddo
    enddo

!Computing tau1=tr[inv_tAkltAlM*tAk]
    tau1=ZERO
    do i=1,nn
      temp1=ZERO
      do k=1,nn
        temp1=temp1+inv_tAkltAlM(i,k)*tAk(k,i)
      enddo
      tau1=tau1+temp1
    enddo
!Computing tau2 = tvk'*inv_tAkltAlM*tAk*inv_tAkltvl
!We do it by multiplying twice the row-vector on the left
!by a matrix on the right and computing a dot product in the end.

!tvkinv_tAkltAlM'=tvk'*inv_tAkltAlM
    do i=1,nn
      temp1=ZERO
      do j=1,nn
        temp1=temp1+tvk(j)*inv_tAkltAlM(j,i)
      enddo
      tvkinv_tAkltAlM(i)=temp1
    enddo

!u1=tvkinv_tAkltAlM'*tAk
!tau2=u1'*inv_tAkltvl (storage for u1 as such is not needed, we use temp1=u1(i))
    tau2=ZERO
    do i=1,nn
      temp1=ZERO
      do j=1,nn
        temp1=temp1+tvkinv_tAkltAlM(j)*tAk(j,i)
      enddo
      tau2=tau2+temp1*inv_tAkltvl(i)
    enddo

!Evaluating the kinetic energy
    Tkl=Skl*(SIX*tau1+FOUR*tau2/tau3)

!Evaluating TrAJ(i,j)=Tr[inv_tAkl*Jij], sqrtTrAJ(i,j)=sqrt(TrAJ(i,j)),
!           eta2(i,j)=tvk'*inv_tAkl*Jij*inv_tAkl*tvl,
!Vkl, (1/r_{ij}^2)_kl, (1/r_{ij})_kl, (r_{ij})_kl, (r_{ij}^2)_kl
!and delta(r_{ij})_kl
    Vkl=ZERO
    temp5=Skl*TWO
    temp1=temp5/Glob_SqrtPi
    temp8=Skl/(Glob_Pi*Glob_SqrtPi)
    do i=1,nn
      temp2=inv_tAkl(i,i)
      TrAJ(i,i)=temp2
      temp3=sqrt(temp2)
      sqrtTrAJ(i,i)=temp3
      !eta2 = tvk'*inv_tAkl*Jii*inv_tAkl*tvl: Jii=e_i*e_i' is rank one,
      !so eta2 = (tvk'*inv_tAkl*e_i)*(e_i'*inv_tAkl*tvl)
      temp4=tvkAj(i,i)*jAtvl(i,i)
      eta2(i,i)=temp4
      temp7=temp4/(temp2*tau3)
      temp6=temp7/THREE
      temp9=ONE-temp7
      rm2kl(i,i)=temp5*(ONE-2*temp6)/temp2
      rmkl(i,i)=temp1*(ONE-temp6)/temp3
      Vkl=Vkl+Glob_ScaledPseudoChargeMatrix(i,0)*rmkl(i,i)
      rkl(i,i)=temp1*(ONE+temp6)*temp3
      r2kl(i,i)=Skl*THREEHALF*(ONE+2*temp6)*temp2
      temp10=temp8/(temp2*temp3)
      deltarkl(i,i)=temp10*temp9
      prvalkl(i,i)=Glob_Pi*temp10*( TWO*temp9*(Glob_EulerConst+log(temp2)) + temp7*FOUR/THREE )
      !prvalkl(i,i)=(temp1/(temp2*temp3))*temp9*(Glob_EulerConst+log(temp2))+(FOUR*Skl/3)*temp7/(temp2*temp3)
      !prvalkl(i,i)=(temp9*(Glob_EulerConst+log(temp2))/Glob_SqrtPi+temp7*TWO/THREE)*temp5/(temp2*temp3)
      !prvalkl(i,i)=(temp1/(temp2*temp3))*(1-temp7)*(Glob_EulerConst+log(temp2))
      !prvalkl(i,i)=TWO*PI*(Glob_EulerConst+log(temp2))*deltarkl(i,i)
    enddo
    do i=1,nn
      do j=i+1,nn
        temp2=inv_tAkl(i,i)+inv_tAkl(j,j)-inv_tAkl(j,i)-inv_tAkl(j,i)
        TrAJ(i,j)=temp2
        TrAJ(j,i)=temp2
        temp3=sqrt(temp2)
        sqrtTrAJ(j,i)=temp3
        sqrtTrAJ(i,j)=temp3
        !eta2 = tvk'*inv_tAkl*Jij*inv_tAkl*tvl: Jij is rank one, so
        !eta2 = (tvk'*inv_tAkl*(e_i-e_j))*((e_i-e_j)'*inv_tAkl*tvl)
        temp4=tvkAj(i,j)*jAtvl(i,j)
        eta2(j,i)=temp4
        eta2(i,j)=temp4
        temp7=temp4/(temp2*tau3)
        temp6=temp7/THREE
        temp9=ONE-temp7
        rm2kl(j,i)=temp5*(ONE-2*temp6)/temp2
        rm2kl(i,j)=rm2kl(j,i)
        rmkl(j,i)=temp1*(ONE-temp6)/temp3
        rmkl(i,j)=rmkl(j,i)
        Vkl=Vkl+Glob_ScaledPseudoChargeMatrix(i,j)*rmkl(j,i)
        rkl(j,i)=temp1*(ONE+temp6)*temp3
        rkl(i,j)=rkl(j,i)
        r2kl(j,i)=Skl*THREEHALF*(ONE+2*temp6)*temp2
        r2kl(i,j)=r2kl(j,i)
        temp10=temp8/(temp2*temp3)
        deltarkl(j,i)=temp10*temp9
        deltarkl(i,j)=deltarkl(j,i)
        prvalkl(j,i)=Glob_Pi*temp10*( TWO*temp9*(Glob_EulerConst+log(temp2)) + temp7*FOUR/THREE )
        !prvalkl(j,i)=(temp1/(temp2*temp3))*temp9*(Glob_EulerConst+log(temp2))+(FOUR*Skl/3)*temp7/(temp2*temp3)
        !prvalkl(j,i)=(temp9*(Glob_EulerConst+log(temp2))/Glob_SqrtPi+temp7*TWO/THREE)*temp5/(temp2*temp3)
        !prvalkl(j,i)=(temp1/(temp2*temp3))*(1-temp7)*(Glob_EulerConst+log(temp2))
        !prvalkl(j,i)=TWO*PI*(Glob_EulerConst+log(temp2))*deltarkl(j,i)
        prvalkl(i,j)=prvalkl(j,i)
      enddo
    enddo

    if (AreMCorrFuncNeeded) then
      do i=1,nn
        MTrAJ(i,i)=inv_invtAkinvtAl(i,i)
        sqrtMTrAJ(i,i)=sqrt(MTrAJ(i,i))
        do q=1,nn
          temp4=ZERO
          do k=1,nn
            temp4=temp4+tvkinv_tAk(k)*inv_invtAkinvtAl(k,i)*inv_invtAkinvtAl(q,i)
          enddo
          u1(q)=temp4
        enddo
        !eta2=u1'*tvl
        temp4=ZERO
        do k=1,nn
          temp4=temp4+u1(k)*inv_tAltvl(k)
        enddo
        meta2(i,i)=temp4
      enddo
      do i=1,nn
        do j=i+1,nn
          MTrAJ(i,j)=inv_invtAkinvtAl(i,i)+inv_invtAkinvtAl(j,j)-inv_invtAkinvtAl(j,i)-inv_invtAkinvtAl(j,i)
          MTrAJ(j,i)=MTrAJ(i,j)
          sqrtMTrAJ(j,i)=sqrt(MTrAJ(j,i))
          sqrtMTrAJ(i,j)=sqrtMTrAJ(j,i)
        enddo
      enddo
    end if

    Hkl=Tkl+Vkl

!Evaluating tr[inv_tAkl Jij inv_tAkl Jpq] and
!j^{ij}' inv_tAkl j^{p,q}  where j^{ij}=e^i-e^j and
!the only nonzero element of e^i is the i-th element
    do i=1,nn
      temp2=inv_tAkl(i,i)
      temp1=temp2*temp2
      TrAJAJ(i,i,i,i)=temp1
      jAj(i,i,i,i)=temp2
      do p=i+1,nn
        do q=p+1,nn
          temp2=inv_tAkl(p,i)-inv_tAkl(q,i)
          temp1=temp2*temp2
          TrAJAJ(i,i,p,q)=temp1
          TrAJAJ(i,i,q,p)=temp1
          TrAJAJ(p,q,i,i)=temp1
          TrAJAJ(q,p,i,i)=temp1
          jAj(i,i,p,q)=temp2
          jAj(i,i,q,p)=-temp2
          jAj(p,q,i,i)=temp2
          jAj(q,p,i,i)=-temp2
        enddo
      enddo
      do j=i+1,nn
        temp2=inv_tAkl(j,i)
        temp1=temp2*temp2
        TrAJAJ(j,j,i,i)=temp1
        TrAJAJ(i,i,j,j)=temp1
        jAj(j,j,i,i)=temp2
        jAj(i,i,j,j)=temp2
        do p=i,nn
          temp2=inv_tAkl(p,i)-inv_tAkl(p,j)
          temp1=temp2*temp2
          TrAJAJ(i,j,p,p)=temp1
          TrAJAJ(j,i,p,p)=temp1
          TrAJAJ(p,p,i,j)=temp1
          TrAJAJ(p,p,j,i)=temp1
          jAj(i,j,p,p)=temp2
          jAj(j,i,p,p)=-temp2
          jAj(p,p,i,j)=temp2
          jAj(p,p,j,i)=-temp2
          do q=p+1,nn
            temp2=inv_tAkl(p,i)-inv_tAkl(q,i)-inv_tAkl(p,j)+inv_tAkl(q,j)
            temp1=temp2*temp2
            TrAJAJ(i,j,p,q)=temp1
            TrAJAJ(j,i,p,q)=temp1
            TrAJAJ(i,j,q,p)=temp1
            TrAJAJ(j,i,q,p)=temp1
            TrAJAJ(p,q,i,j)=temp1
            TrAJAJ(p,q,j,i)=temp1
            TrAJAJ(q,p,i,j)=temp1
            TrAJAJ(q,p,j,i)=temp1
            jAj(i,j,p,q)=temp2
            jAj(j,i,p,q)=-temp2
            jAj(i,j,q,p)=-temp2
            jAj(j,i,q,p)=temp2
            jAj(p,q,i,j)=temp2
            jAj(p,q,j,i)=-temp2
            jAj(q,p,i,j)=-temp2
            jAj(q,p,j,i)=temp2
          enddo
        enddo
      enddo
    enddo

!!Evaluating tr [inv_tAkl J_{ij} inv_tAkl J_{pq} inv_tAkl V] where V = tvk tvl'
!!and tr [inv_tAkl J_{ij} inv_tAkl J_{pq} inv_tAkl J_{ij} inv_tAkl V]
!do i=1,nn
! do j=i,nn
!   do p=1,nn
!     do q=p,nn
!       temp1=tvkAj(i,j)*jAj(i,j,p,q)*jAtvl(p,q)
!       temp2=tvkAj(i,j)*jAj(i,j,p,q)*jAj(p,q,i,j)*jAtvl(i,j)
!       TrAJAJAV(i,j,p,q)=temp1
!       TrAJAJAV(j,i,p,q)=temp1
!       TrAJAJAV(i,j,q,p)=temp1
!       TrAJAJAV(j,i,q,p)=temp1
!       omega(i,j,p,q)=temp2
!       omega(j,i,p,q)=temp2
!       omega(i,j,q,p)=temp2
!       omega(j,i,q,p)=temp2
!     enddo
!   enddo
! enddo
!enddo

!Evaluation of (1/r_{ij}*1/r_{pq}))_kl is not implemented yet
    temp1=4*Skl/(3*Glob_Pi)
    do i=1,nn
      do j=i,nn
        do p=i,nn
          do q=p,nn   !try q=max(p,j),n - it may speed things up a little
            if (((p==i).and.(q==j)).or.((p==j).and.(q==i))) then
              temp2=rm2kl(i,j)
              rmrmkl(i,j,p,q)=temp2
              rmrmkl(j,i,p,q)=temp2
              rmrmkl(i,j,q,p)=temp2
              rmrmkl(j,i,q,p)=temp2
              rmrmkl(p,q,i,j)=temp2
              rmrmkl(p,q,j,i)=temp2
              rmrmkl(q,p,i,j)=temp2
              rmrmkl(q,p,j,i)=temp2
            else
              !jAj(nn,nn,nn,nn),jAtvl(nn,nn),tvkAj(nn,nn)
              temp2=TrAJ(i,j)          !sigma    a
              temp3=TrAJ(p,q)          !chi      b
              temp4=jAj(p,q,i,j)       !         c
              !tau3                    !zeta     w
              temp5=tvkAj(i,j)         !         f
              temp6=tvkAj(p,q)         !         g
              temp7=jAtvl(i,j)         !         u
              temp8=jAtvl(p,q)         !         v
              temp9 = sqrt(temp2*temp3)
              temp10 = temp4/temp9
              threshold = ONE - temp10*temp10
              if ( threshold  < TEN * EPSILON(threshold) )then
                ! employed expression
                ! expansion_term_2 = 2*Skl* (Sqrt(a*b)*g*u + Sqrt(a*b)*f*v - 3*a*b*tau3) / (3*a*b*Sqrt(a*b)*tau3)
                temp11 = temp9 * temp6 * temp7                                                              ! Sqrt(a*b)*g*u
                temp12 = temp9 * temp5 * temp8                                                              ! Sqrt(a*b)*f*v
                temp13 = THREE * temp2 * temp3 * tau3                                                       ! 3*a*b*tau3
                temp14 = (TWO * Skl) * (temp11 + temp12 - temp13) / (temp13 * temp9)
              else
                temp11=sqrt(ONE-temp10*temp10)                                                                              ! sqrt(1-x^2)
                temp12=(temp6*temp7+temp5*temp8)/(temp2*temp3*tau3)                                                         ! h
                temp13=(temp2*temp6*temp8+temp3*temp5*temp7)/(temp2*temp3*temp9*tau3)                                       ! t
                temp14=(temp1/temp11)*(THREE/temp9-temp13+(temp12-temp10*THREE/temp9)*ftransaux(temp10))
              endif
              rmrmkl(i,j,p,q)=temp14
              rmrmkl(j,i,p,q)=temp14
              rmrmkl(i,j,q,p)=temp14
              rmrmkl(j,i,q,p)=temp14
              rmrmkl(p,q,i,j)=temp14
              rmrmkl(p,q,j,i)=temp14
              rmrmkl(q,p,i,j)=temp14
              rmrmkl(q,p,j,i)=temp14
            endif
          enddo
        enddo
      enddo
    enddo

!write(*,*) '===================================='
!write(*,'(1x,a,e23.16)') 'Skl=',Skl
!write(*,*) 'tvk='
!do i=1,nn
!    write(*,'(2x,i3)') tvk(i)
!enddo
!write(*,*) 'tvl='
!do i=1,nn
!    write(*,'(2x,i3)') tvl(i)
!enddo
!write(*,*) 'tAk='
!do i=1,nn
!  do j=1,nn
!    write(*,'(2x,f23.20)',advance='no') tAk(i,j)
!  enddo
!  write(*,*)
!enddo
!write(*,*) 'tAl='
!do i=1,nn
!  do j=1,nn
!    write(*,'(2x,f23.20)',advance='no') tAl(i,j)
!  enddo
!  write(*,*)
!enddo
!write(*,*) 'inv_tAkl='
!do i=1,nn
!  do j=1,nn
!    write(*,'(2x,f23.20)',advance='no') inv_tAkl(i,j)
!  enddo
!  write(*,*)
!enddo
!write(*,*) 'M='
!do i=1,nn
!  do j=1,nn
!    write(*,'(2x,f23.20)',advance='no') Glob_MassMatrix(i,j)
!  enddo
!  write(*,*)
!enddo
!write(*,*) '1/(rijrpq)='
!do i=1,nn
!  do j=1,nn
!    do p=1,nn
!      do q=1,nn
!        write(*,'(2x,i3,i3,i3,i3,2x,e23.16)') i,j,p,q,rmrmkl(i,j,p,q)
!      enddo
!    enddo
!  enddo
!enddo
!write(*,*) 'ME_d_X_over_rij_d='

!Loop that computes all drachmanized delta(r_{ij})_kl  as well as V^2_kl
!The original code called ME_d_X_over_rij_d(p,q,Glob_dmvM,...) for every
!pair, redoing four O(n^3) matrix products with the same X=Glob_dmvM
!each time. All pair-independent pieces are precomputed here instead,
!using inv_tAkl*tAl = inv_tAkltAl (=:U) and inv_tAkl*tAk = I-U:
!  Qsd    : symmetrized inv_tAkl*(tAl*X*tAk)*inv_tAkl = sym(U*X*(I-U'))
!           (the symmetrization accounts for the tAk*X*tAl term)
!  dtheta : tr[inv_tAkl*tAl*X*tAk] = tr[(U*X)*tAk]
!  dlam   : tvkinv_tAkl'*(tAl*X*tAk)*inv_tAkltvl = tAlAvk'*X*(tvl-tAlAvl)
!  gvo    : (inv_tAkl*(tAl*X*tAk)'*tvk_A + inv_tAkl*(tAk*X*tAl)'*tvk_A)/2
!  gvc    : (inv_tAkl*(tAl*X*tAk)*Avl + inv_tAkl*(tAk*X*tAl)*Avl)/2
!  gv3    : inv_tAkl*tAk*X*tvl = (I-U)*(X*tvl)
!  gv4    : inv_tAkl*tAl*X*tvk = U*(X*tvk)
!after which each pair costs O(1) (contractions with inv_tAkl*(e_p-e_q)
!become differences of vector elements).
    do i=1,nn
      do j=1,nn
        temp1=ZERO
        do k=1,nn
          temp1=temp1+inv_tAkltAl(j,k)*Glob_dmvM(k,i)
        enddo
        UXd(j,i)=temp1
      enddo
    enddo
    dtheta=ZERO
    do i=1,nn
      do k=1,nn
        dtheta=dtheta+UXd(i,k)*tAk(k,i)
      enddo
    enddo
    do i=1,nn
      do j=1,nn
        temp1=ZERO
        do k=1,nn
          temp1=temp1+UXd(j,k)*inv_tAkltAl(i,k)
        enddo
        W1(j,i)=UXd(j,i)-temp1
      enddo
    enddo
    do i=1,nn
      do j=i,nn
        temp1=ONEHALF*(W1(j,i)+W1(i,j))
        Qsd(j,i)=temp1
        Qsd(i,j)=temp1
      enddo
    enddo
    do i=1,nn
      temp1=ZERO
      do k=1,nn
        temp1=temp1+Glob_dmvM(i,k)*tAlAvk(k)
      enddo
      wvd(i)=temp1
    enddo
    do i=1,nn
      temp1=ZERO
      do k=1,nn
        temp1=temp1+inv_tAkltAl(i,k)*wvd(k)
      enddo
      gvo(i)=wvd(i)-temp1
    enddo
    do i=1,nn
      temp1=ZERO
      do k=1,nn
        temp1=temp1+Glob_dmvM(i,k)*(tvk(k)-tAlAvk(k))
      enddo
      wvd(i)=temp1
    enddo
    do i=1,nn
      temp1=ZERO
      do k=1,nn
        temp1=temp1+inv_tAkltAl(i,k)*wvd(k)
      enddo
      gvo(i)=ONEHALF*(gvo(i)+temp1)
    enddo
    dlam=ZERO
    do i=1,nn
      temp1=ZERO
      do k=1,nn
        temp1=temp1+Glob_dmvM(i,k)*(tvl(k)-tAlAvl(k))
      enddo
      wvd(i)=temp1
      dlam=dlam+tAlAvk(i)*temp1
    enddo
    do i=1,nn
      temp1=ZERO
      do k=1,nn
        temp1=temp1+inv_tAkltAl(i,k)*wvd(k)
      enddo
      gvc(i)=temp1
    enddo
    do i=1,nn
      temp1=ZERO
      do k=1,nn
        temp1=temp1+Glob_dmvM(i,k)*tAlAvl(k)
      enddo
      wvd(i)=temp1
    enddo
    do i=1,nn
      temp1=ZERO
      do k=1,nn
        temp1=temp1+inv_tAkltAl(i,k)*wvd(k)
      enddo
      gvc(i)=ONEHALF*(gvc(i)+wvd(i)-temp1)
    enddo
    do i=1,nn
      temp1=ZERO
      do k=1,nn
        temp1=temp1+Glob_dmvM(i,k)*tvl(k)
      enddo
      wvd(i)=temp1
    enddo
    do i=1,nn
      temp1=ZERO
      do k=1,nn
        temp1=temp1+inv_tAkltAl(i,k)*wvd(k)
      enddo
      gv3(i)=wvd(i)-temp1
    enddo
    do i=1,nn
      temp1=ZERO
      do k=1,nn
        temp1=temp1+Glob_dmvM(i,k)*tvk(k)
      enddo
      wvd(i)=temp1
    enddo
    do i=1,nn
      temp1=ZERO
      do k=1,nn
        temp1=temp1+inv_tAkltAl(i,k)*wvd(k)
      enddo
      gv4(i)=temp1
    enddo

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
        temp4=ZERO
        temp5=ZERO
        if (p==q) then
          temp4=2*Glob_Pi*Glob_MassMatrix(p,p)
          temp5=Glob_ScaledPseudoChargeMatrix(0,p)
          dkap=Qsd(p,p)
          dt1=inv_tAkltvl(p)
          dt2=tvkinv_tAkl(p)
          dom=gvo(p)*dt1
          dch=dt2*gvc(p)
          dhh=dt2*gv3(p)
          dmm=gv4(p)*dt1
        else
          temp4=2*Glob_Pi*(Glob_MassMatrix(p,p)+Glob_MassMatrix(q,q) &
                      -Glob_MassMatrix(p,q)-Glob_MassMatrix(p,q))
          temp5=Glob_ScaledPseudoChargeMatrix(p,q)
          dkap=Qsd(p,p)+Qsd(q,q)-2*Qsd(q,p)
          dt1=inv_tAkltvl(p)-inv_tAkltvl(q)
          dt2=tvkinv_tAkl(p)-tvkinv_tAkl(q)
          dom=(gvo(p)-gvo(q))*dt1
          dch=dt2*(gvc(p)-gvc(q))
          dhh=dt2*(gv3(p)-gv3(q))
          dmm=(gv4(p)-gv4(q))*dt1
        endif

        !temp2 = ME_d_X_over_rij_d(p,q,Glob_dmvM,...) assembled from the
        !precomputed pieces (same final formula as in the routine)
        temp3=(4*Skl)/(15*Glob_SqrtPi*tau3*trAJ(p,q)*trAJ(p,q)*sqrt(trAJ(p,q)))
        temp2=temp3*( 15*dkap*eta2(p,q) + &
               5*trAJ(p,q)*(6*dlam*trAJ(p,q)+9*tau3*dtheta*trAJ(p,q) &
              -3*dkap*tau3-2*dom-2*dch+dhh+dmm-3*dtheta*eta2(p,q)) )
        drach_deltarkl(p,q)=(Glob_CurrEnergy*rmkl(p,q)-temp1-temp2)/temp4
        drach_deltarkl(q,p)=drach_deltarkl(p,q)
        V2kl=V2kl+temp5*temp1
      enddo
    enddo

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

!Mass-velocity correction.
!ME_dWd2 is called by the original code with rank-one matrices W: the
!all-ones matrix (= w*w' with w=(1,..,1)') and Eii (= e_i*e_i'). For
!W = w*w' every matrix inside ME_dWd2 is rank one (tAk*W*tAk = qk*qk'
!with qk = tAk*w, etc.), so all its traces reduce to dot products of
!qk = tAk*w and ql = tAl*w with inv_tAkl. For w = e_i even those are
!O(1) lookups (using inv_tAkl*tAk = I - inv_tAkltAl):
!  qk'*inv_tAkl*qk = tAk(i,i)-tAl(i,i)+GAl(i,i)
!  ql'*inv_tAkl*ql = GAl(i,i),  qk'*inv_tAkl*ql = tAl(i,i)-GAl(i,i)
!  tvk'*inv_tAkl*qk = tvk(i)-tAlAvk(i),  tvk'*inv_tAkl*ql = tAlAvk(i)
!  qk'*inv_tAkl*tvl = tvl(i)-tAlAvl(i),  ql'*inv_tAkl*tvl = tAlAvl(i)
!Only the general-matrix call with W = Glob_dmvM (drachmanized MV)
!still goes through ME_dWd2 itself.
    inv_tau3=1/tau3
!W = all ones: qk=sk, ql=sl are the row sum vectors of tAk, tAl
    do p=1,nn
      temp1=ZERO
      temp2=ZERO
      do q=1,nn
        temp1=temp1+tAk(q,p)
        temp2=temp2+tAl(q,p)
      enddo
      sk(p)=temp1
      sl(p)=temp2
    enddo
    do p=1,nn
      temp1=ZERO
      temp2=ZERO
      do q=1,nn
        temp1=temp1+inv_tAkl(p,q)*sk(q)
        temp2=temp2+inv_tAkl(p,q)*sl(q)
      enddo
      Ask(p)=temp1
      Asl(p)=temp2
    enddo
    qkk=ZERO
    qll=ZERO
    qkl=ZERO
    vkqk=ZERO
    vkql=ZERO
    qkvl=ZERO
    qlvl=ZERO
    wvk=ZERO
    wvl=ZERO
    wqk=ZERO
    wql=ZERO
    do p=1,nn
      qkk=qkk+sk(p)*Ask(p)
      qll=qll+sl(p)*Asl(p)
      qkl=qkl+sk(p)*Asl(p)
      vkqk=vkqk+tvk(p)*Ask(p)
      vkql=vkql+tvk(p)*Asl(p)
      qkvl=qkvl+sk(p)*inv_tAkltvl(p)
      qlvl=qlvl+sl(p)*inv_tAkltvl(p)
      wvk=wvk+tvk(p)
      wvl=wvl+tvl(p)
      wqk=wqk+sk(p)
      wql=wql+sl(p)
    enddo
    trAvltuk=4*wvk*qkvl+6*wqk*tau3
    trAtulvk=4*wvl*vkql+6*wql*tau3
    trAultuk=16*wvk*wvl*qkl+24*wvk*wql*qkvl+24*wqk*wvl*vkql+36*wqk*wql*tau3
    trKx=vkqk*(4*wvl*qkl+6*wql*qkvl)
    trLx=(4*wvk*qkl+6*wqk*vkql)*qlvl
    merr=Skl*( THREEHALF*THREEHALF*qkk*qll + THREEHALF*qkl*qkl &
         +(THREEHALF*(qkk*vkql*qlvl+qll*vkqk*qkvl) &
         +vkqk*qkl*qlvl+vkql*qkl*qkvl)*inv_tau3 )
    dwd2=16*merr+Skl*inv_tau3*( trAultuk-SIX*(qkk*trAtulvk+qll*trAvltuk) &
         -FOUR*(trKx+trLx) )
    MVkl=dwd2/(Glob_Mass(1)*Glob_Mass(1)*Glob_Mass(1))
!W = Eii, i=1..n: everything is an O(1) lookup
    do i=1,nn
      qkk=tAk(i,i)-tAl(i,i)+GAl(i,i)
      qll=GAl(i,i)
      qkl=tAl(i,i)-GAl(i,i)
      vkqk=tvk(i)-tAlAvk(i)
      vkql=tAlAvk(i)
      qkvl=tvl(i)-tAlAvl(i)
      qlvl=tAlAvl(i)
      wvk=tvk(i)
      wvl=tvl(i)
      wqk=tAk(i,i)
      wql=tAl(i,i)
      trAvltuk=4*wvk*qkvl+6*wqk*tau3
      trAtulvk=4*wvl*vkql+6*wql*tau3
      trAultuk=16*wvk*wvl*qkl+24*wvk*wql*qkvl+24*wqk*wvl*vkql+36*wqk*wql*tau3
      trKx=vkqk*(4*wvl*qkl+6*wql*qkvl)
      trLx=(4*wvk*qkl+6*wqk*vkql)*qlvl
      merr=Skl*( THREEHALF*THREEHALF*qkk*qll + THREEHALF*qkl*qkl &
           +(THREEHALF*(qkk*vkql*qlvl+qll*vkqk*qkvl) &
           +vkqk*qkl*qlvl+vkql*qkl*qkvl)*inv_tau3 )
      dwd2=16*merr+Skl*inv_tau3*( trAultuk-SIX*(qkk*trAtulvk+qll*trAvltuk) &
           -FOUR*(trKx+trLx) )
      MVkl=MVkl+dwd2/(Glob_Mass(i+1)*Glob_Mass(i+1)*Glob_Mass(i+1))
    enddo
    MVkl=-MVkl/8
    drach_MVkl1=ZERO
    drach_MVkl2=ZERO

    temp1=ME_dWd2(Glob_dmvM,tAk,tAl,inv_tAkl,tvk,tvl,inv_tAkltvl,tvkinv_tAkl,inv_tau3,Skl) &
           -V2kl-Glob_CurrEnergy*Glob_CurrEnergy*Skl+2*Glob_CurrEnergy*Vkl
    drach_MVkl1 = temp1*Glob_dmva21 + MVkl
    drach_MVkl2 = temp1*Glob_dmva22 + MVkl

!Evaluating Orbit-Orbit (OO) matrix element (without the factor of alpha**2)
!All seven matrices W1..W7 that the original code assembled per pair are
!sums of rank-one terms u*v' built from a=tAkl(:,j), b=tAl(:,j),
!c=tAl(:,i), d (=e_j in the first sum, e_j-e_i in the second), tvk and
!tvl:
!  W1 = c*b' + b*c' + tr1*a*d' + d*(tr1*b'+tr3*c')
!  W2 = W6 = a*b'
!  W3 = d*c'
!  W4 = 2*tvl(i)*b*tvk' + 2*tvl(j)*c*tvk' + tr1*tvk(j)*tvl*d'
!       + (2*tr1*tvl(j)+tr3*tvl(i))*d*tvk' + tvl(i)*tvk(j)*d*b'
!  W5 = tvl(j)*a*tvk' + tvl(j)*b*tvk' + tvk(j)*tvl*b'
!  W7 = tvl(i)*d*tvk'
!Every trace formed from their symmetrized versions inside the integral
!routines (ME_rXr_over_rij, ME_rXr_rYr_over_rij, SG_ME_rXr_over_rij,
!SG_ME_rXr_rYr_over_rij) reduces to O(1) combinations of the scalars
!x'*inv_tAkl*y, x,y in {a,b,c,d,tvk,tvl}, denoted s_xy below. All of
!them come from GAl, inv_tAkltAl, tAlAvk, tAlAvl, tvkinv_tAkl,
!inv_tAkltvl, TrAJ and tAl itself (note inv_tAkl*a = e_j exactly, and
!the J matrix of the integral pair equals d*d'). The trace expressions
!were generated symbolically from the rank-one expansions; the final
!formulas are verbatim those of the integral routines. This makes the
!whole OO evaluation O(n^2) instead of O(n^5).
    OOkl=ZERO
!First double loop for OO (integral pair (j,j), d=e_j)
    do i=1,nn
      do j=1,nn
        tr1=tAl(j,i)
        tr3=3*tAl(j,j)
        tr4=tvl(j)*tvk(j)
        vkj=tvk(j)
        vli=tvl(i)
        vlj=tvl(j)
        s_ab=tAl(j,j)
        s_ad=ONE
        s_ak=tvk(j)
        s_al=tvl(j)
        s_bc=GAl(j,i)
        s_bd=inv_tAkltAl(j,j)
        s_bk=tAlAvk(j)
        s_bl=tAlAvl(j)
        s_cd=inv_tAkltAl(j,i)
        s_ck=tAlAvk(i)
        s_cl=tAlAvl(i)
        s_dd=TrAJ(j,j)
        s_dk=tvkinv_tAkl(j)
        s_dl=inv_tAkltvl(j)
        !building-block traces of the symmetrized rank-one expansions
        tX1=s_ad*tr1+2*s_bc+s_bd*tr1+s_cd*tr3
        tXJ1=s_ad*s_dd*tr1+2*s_bd*s_cd+s_bd*s_dd*tr1+s_cd*s_dd*tr3
        tXV1=ONEHALF*s_ak*s_dl*tr1+ONEHALF*s_al*s_dk*tr1+s_bk*s_cl &
            +ONEHALF*s_bk*s_dl*tr1+s_bl*s_ck+ONEHALF*s_bl*s_dk*tr1 &
            +ONEHALF*s_ck*s_dl*tr3+ONEHALF*s_cl*s_dk*tr3
        kXd1=ONEHALF*s_ad*s_dk*tr1+ONEHALF*s_ak*s_dd*tr1+s_bd*s_ck &
            +ONEHALF*s_bd*s_dk*tr1+s_bk*s_cd+ONEHALF*s_bk*s_dd*tr1 &
            +ONEHALF*s_cd*s_dk*tr3+ONEHALF*s_ck*s_dd*tr3
        dXl1=ONEHALF*s_ad*s_dl*tr1+ONEHALF*s_al*s_dd*tr1+s_bd*s_cl &
            +ONEHALF*s_bd*s_dl*tr1+s_bl*s_cd+ONEHALF*s_bl*s_dd*tr1 &
            +ONEHALF*s_cd*s_dl*tr3+ONEHALF*s_cl*s_dd*tr3
        tX2=s_ab
        tXJ2=s_ad*s_bd
        tXV2=ONEHALF*(s_ak*s_bl+s_al*s_bk)
        kXd2=ONEHALF*(s_ad*s_bk+s_ak*s_bd)
        dXl2=ONEHALF*(s_ad*s_bl+s_al*s_bd)
        tY3=s_cd
        tYJ3=s_cd*s_dd
        tYV3=ONEHALF*(s_ck*s_dl+s_cl*s_dk)
        kYd3=ONEHALF*(s_cd*s_dk+s_ck*s_dd)
        dYl3=ONEHALF*(s_cd*s_dl+s_cl*s_dd)
        tX4=s_bd*vkj*vli+2*s_bk*vli+2*s_ck*vlj+2*s_dk*tr1*vlj+s_dk*tr3*vli+s_dl*tr1*vkj
        tXJ4=s_bd*s_dd*vkj*vli+2*s_bd*s_dk*vli+2*s_cd*s_dk*vlj &
            +2*s_dd*s_dk*tr1*vlj+s_dd*s_dk*tr3*vli+s_dd*s_dl*tr1*vkj
        tY5=s_ak*vlj+s_bk*vlj+s_bl*vkj
        tYJ5=s_ad*s_dk*vlj+s_bd*s_dk*vlj+s_bd*s_dl*vkj
        tY7=s_dk*vli
        tYJ7=s_dd*s_dk*vli
        tXY23=ONEHALF*(tr1*s_bd+s_ad*s_bc)
        tXYJ23=ONEFOURTH*tr1*s_bd*s_dd+ONEFOURTH*s_ad*s_bc*s_dd+ONEHALF*s_ad*s_bd*s_cd
        tXYV23=ONEFOURTH*(tr1*s_bk*s_dl+s_ad*s_bk*s_cl+s_ak*s_bc*s_dl+s_ak*s_bd*s_cl)
        tYXV23=ONEFOURTH*(tr1*s_bl*s_dk+s_ad*s_bl*s_ck+s_al*s_bc*s_dk+s_al*s_bd*s_ck)
        kXYd23=ONEFOURTH*(tr1*s_bk*s_dd+s_ad*s_bk*s_cd+s_ak*s_bc*s_dd+s_ak*s_bd*s_cd)
        kYXd23=ONEFOURTH*tr1*s_bd*s_dk+ONEFOURTH*s_ad*s_bc*s_dk+ONEHALF*s_ad*s_bd*s_ck
        dXYl23=ONEFOURTH*tr1*s_bd*s_dl+ONEFOURTH*s_ad*s_bc*s_dl+ONEHALF*s_ad*s_bd*s_cl
        dYXl23=ONEFOURTH*(tr1*s_bl*s_dd+s_ad*s_bl*s_cd+s_al*s_bc*s_dd+s_al*s_bd*s_cd)
        tXY35=ONEHALF*(tr1*s_dk*vlj+s_ad*s_ck*vlj+s_bc*s_dk*vlj+s_bc*s_dl*vkj &
             +s_bd*s_ck*vlj+s_bd*s_cl*vkj)
        tXYJ35=ONEFOURTH*tr1*s_dd*s_dk*vlj+ONEHALF*s_ad*s_cd*s_dk*vlj &
              +ONEFOURTH*s_ad*s_ck*s_dd*vlj+ONEFOURTH*s_bc*s_dd*s_dk*vlj &
              +ONEFOURTH*s_bc*s_dd*s_dl*vkj+ONEHALF*s_bd*s_cd*s_dk*vlj &
              +ONEHALF*s_bd*s_cd*s_dl*vkj+ONEFOURTH*s_bd*s_ck*s_dd*vlj &
              +ONEFOURTH*s_bd*s_cl*s_dd*vkj
        tXY27=ONEHALF*vli*(s_ad*s_bk+s_ak*s_bd)
        tXYJ27=ONEFOURTH*vli*(2*s_ad*s_bd*s_dk+s_ad*s_bk*s_dd+s_ak*s_bd*s_dd)
        !the integrals (same formulas as in the integral routines)
        odd=ONE/s_dd
        oc=(TWO/Glob_SqrtPi)*Skl*inv_tau3*odd*sqrt(odd)
        tJVs=s_dk*s_dl
        temp1=oc*( THREEHALF*s_dd*tau3*tX1-ONEHALF*(tX1*tJVs+tau3*tXJ1) &
              +s_dd*tXV1+ONEHALF*odd*tJVs*tXJ1-ONETHIRD*(s_dk*dXl1+s_dl*kXd1) )
        temp2=oc*( NINE*ONEFOURTH*s_dd*tau3*tX2*tY3 &
              -THREE*ONEFOURTH*(tX2*tY3*tJVs+tau3*tY3*tXJ2+tau3*tX2*tYJ3) &
              +THREEHALF*(s_dd*tY3*tXV2+s_dd*tX2*tYV3+s_dd*tau3*tXY23) &
              -ONEHALF*(tJVs*tXY23+tXV2*tYJ3+tXJ2*tYV3) &
              +THREE*ONEFOURTH*odd*(tY3*tJVs*tXJ2+tX2*tJVs*tYJ3+tau3*tXJ2*tYJ3) &
              -FIVE*ONEFOURTH*odd*odd*tJVs*tXJ2*tYJ3 &
              -ONEHALF*(tX2*(s_dk*dYl3)+tau3*tXYJ23+tY3*(s_dk*dXl2) &
                       +tY3*(s_dl*kXd2)+tX2*(s_dl*kYd3)+tau3*tXYJ23) &
              +s_dd*tXYV23+s_dd*tYXV23 &
              +ONEHALF*odd*(tXJ2*(s_dk*dYl3)+tJVs*tXYJ23+tYJ3*(s_dk*dXl2) &
                           +tYJ3*(s_dl*kXd2)+tXJ2*(s_dl*kYd3)+tJVs*tXYJ23) &
              -ONETHIRD*(kYd3*dXl2+kXYd23*s_dl+kYXd23*s_dl+s_dk*dXYl23 &
                        +s_dk*dYXl23+kXd2*dYl3) )
        temp3=oc*(THREE*s_dd*tY3-tYJ3)
        temp4=ONETHIRD*oc*(THREE*s_dd*tX4-tXJ4)
        temp5=ONETHIRD*THREE*oc*( THREEHALF*s_dd*tY3*tY5-ONEHALF*(tY5*tYJ3+tY3*tYJ5) &
              +s_dd*tXY35-ONETHIRD*(tXYJ35+tXYJ35)+ONEHALF*odd*tYJ3*tYJ5 )
        temp6=ONETHIRD*THREE*oc*( THREEHALF*s_dd*tX2*tY7-ONEHALF*(tY7*tXJ2+tX2*tYJ7) &
              +s_dd*tXY27-ONETHIRD*(tXYJ27+tXYJ27)+ONEHALF*odd*tXJ2*tYJ7 )
        temp7=-12*tr1*rmkl(j,j)+4*temp1-8*temp2-2*tr4*temp3-2*temp4+4*temp5+4*temp6
        OOkl=OOkl-temp7*Glob_ScaledPseudoChargeMatrix(j,0)/Glob_Mass(j+1)
      enddo
    enddo
    OOkl=OOkl/Glob_Mass(1)

!Second double loop for OO (integral pair (i,j), d=e_j-e_i)
    do i=1,nn
      do j=i+1,nn
        tr1=tAl(j,i)
        tr3=3*tAl(j,j)
        tr4=tvl(j)*tvk(j)
        vkj=tvk(j)
        vli=tvl(i)
        vlj=tvl(j)
        s_ab=tAl(j,j)
        s_ad=ONE
        s_ak=tvk(j)
        s_al=tvl(j)
        s_bc=GAl(j,i)
        s_bd=inv_tAkltAl(j,j)-inv_tAkltAl(i,j)
        s_bk=tAlAvk(j)
        s_bl=tAlAvl(j)
        s_cd=inv_tAkltAl(j,i)-inv_tAkltAl(i,i)
        s_ck=tAlAvk(i)
        s_cl=tAlAvl(i)
        s_dd=TrAJ(i,j)
        s_dk=tvkinv_tAkl(j)-tvkinv_tAkl(i)
        s_dl=inv_tAkltvl(j)-inv_tAkltvl(i)
        !building-block traces of the symmetrized rank-one expansions
        tX1=s_ad*tr1+2*s_bc+s_bd*tr1+s_cd*tr3
        tXJ1=s_ad*s_dd*tr1+2*s_bd*s_cd+s_bd*s_dd*tr1+s_cd*s_dd*tr3
        tXV1=ONEHALF*s_ak*s_dl*tr1+ONEHALF*s_al*s_dk*tr1+s_bk*s_cl &
            +ONEHALF*s_bk*s_dl*tr1+s_bl*s_ck+ONEHALF*s_bl*s_dk*tr1 &
            +ONEHALF*s_ck*s_dl*tr3+ONEHALF*s_cl*s_dk*tr3
        kXd1=ONEHALF*s_ad*s_dk*tr1+ONEHALF*s_ak*s_dd*tr1+s_bd*s_ck &
            +ONEHALF*s_bd*s_dk*tr1+s_bk*s_cd+ONEHALF*s_bk*s_dd*tr1 &
            +ONEHALF*s_cd*s_dk*tr3+ONEHALF*s_ck*s_dd*tr3
        dXl1=ONEHALF*s_ad*s_dl*tr1+ONEHALF*s_al*s_dd*tr1+s_bd*s_cl &
            +ONEHALF*s_bd*s_dl*tr1+s_bl*s_cd+ONEHALF*s_bl*s_dd*tr1 &
            +ONEHALF*s_cd*s_dl*tr3+ONEHALF*s_cl*s_dd*tr3
        tX2=s_ab
        tXJ2=s_ad*s_bd
        tXV2=ONEHALF*(s_ak*s_bl+s_al*s_bk)
        kXd2=ONEHALF*(s_ad*s_bk+s_ak*s_bd)
        dXl2=ONEHALF*(s_ad*s_bl+s_al*s_bd)
        tY3=s_cd
        tYJ3=s_cd*s_dd
        tYV3=ONEHALF*(s_ck*s_dl+s_cl*s_dk)
        kYd3=ONEHALF*(s_cd*s_dk+s_ck*s_dd)
        dYl3=ONEHALF*(s_cd*s_dl+s_cl*s_dd)
        tX4=s_bd*vkj*vli+2*s_bk*vli+2*s_ck*vlj+2*s_dk*tr1*vlj+s_dk*tr3*vli+s_dl*tr1*vkj
        tXJ4=s_bd*s_dd*vkj*vli+2*s_bd*s_dk*vli+2*s_cd*s_dk*vlj &
            +2*s_dd*s_dk*tr1*vlj+s_dd*s_dk*tr3*vli+s_dd*s_dl*tr1*vkj
        tY5=s_ak*vlj+s_bk*vlj+s_bl*vkj
        tYJ5=s_ad*s_dk*vlj+s_bd*s_dk*vlj+s_bd*s_dl*vkj
        tY7=s_dk*vli
        tYJ7=s_dd*s_dk*vli
        tXY23=ONEHALF*(tr1*s_bd+s_ad*s_bc)
        tXYJ23=ONEFOURTH*tr1*s_bd*s_dd+ONEFOURTH*s_ad*s_bc*s_dd+ONEHALF*s_ad*s_bd*s_cd
        tXYV23=ONEFOURTH*(tr1*s_bk*s_dl+s_ad*s_bk*s_cl+s_ak*s_bc*s_dl+s_ak*s_bd*s_cl)
        tYXV23=ONEFOURTH*(tr1*s_bl*s_dk+s_ad*s_bl*s_ck+s_al*s_bc*s_dk+s_al*s_bd*s_ck)
        kXYd23=ONEFOURTH*(tr1*s_bk*s_dd+s_ad*s_bk*s_cd+s_ak*s_bc*s_dd+s_ak*s_bd*s_cd)
        kYXd23=ONEFOURTH*tr1*s_bd*s_dk+ONEFOURTH*s_ad*s_bc*s_dk+ONEHALF*s_ad*s_bd*s_ck
        dXYl23=ONEFOURTH*tr1*s_bd*s_dl+ONEFOURTH*s_ad*s_bc*s_dl+ONEHALF*s_ad*s_bd*s_cl
        dYXl23=ONEFOURTH*(tr1*s_bl*s_dd+s_ad*s_bl*s_cd+s_al*s_bc*s_dd+s_al*s_bd*s_cd)
        tXY35=ONEHALF*(tr1*s_dk*vlj+s_ad*s_ck*vlj+s_bc*s_dk*vlj+s_bc*s_dl*vkj &
             +s_bd*s_ck*vlj+s_bd*s_cl*vkj)
        tXYJ35=ONEFOURTH*tr1*s_dd*s_dk*vlj+ONEHALF*s_ad*s_cd*s_dk*vlj &
              +ONEFOURTH*s_ad*s_ck*s_dd*vlj+ONEFOURTH*s_bc*s_dd*s_dk*vlj &
              +ONEFOURTH*s_bc*s_dd*s_dl*vkj+ONEHALF*s_bd*s_cd*s_dk*vlj &
              +ONEHALF*s_bd*s_cd*s_dl*vkj+ONEFOURTH*s_bd*s_ck*s_dd*vlj &
              +ONEFOURTH*s_bd*s_cl*s_dd*vkj
        tXY27=ONEHALF*vli*(s_ad*s_bk+s_ak*s_bd)
        tXYJ27=ONEFOURTH*vli*(2*s_ad*s_bd*s_dk+s_ad*s_bk*s_dd+s_ak*s_bd*s_dd)
        !the integrals (same formulas as in the integral routines)
        odd=ONE/s_dd
        oc=(TWO/Glob_SqrtPi)*Skl*inv_tau3*odd*sqrt(odd)
        tJVs=s_dk*s_dl
        temp1=oc*( THREEHALF*s_dd*tau3*tX1-ONEHALF*(tX1*tJVs+tau3*tXJ1) &
              +s_dd*tXV1+ONEHALF*odd*tJVs*tXJ1-ONETHIRD*(s_dk*dXl1+s_dl*kXd1) )
        temp2=oc*( NINE*ONEFOURTH*s_dd*tau3*tX2*tY3 &
              -THREE*ONEFOURTH*(tX2*tY3*tJVs+tau3*tY3*tXJ2+tau3*tX2*tYJ3) &
              +THREEHALF*(s_dd*tY3*tXV2+s_dd*tX2*tYV3+s_dd*tau3*tXY23) &
              -ONEHALF*(tJVs*tXY23+tXV2*tYJ3+tXJ2*tYV3) &
              +THREE*ONEFOURTH*odd*(tY3*tJVs*tXJ2+tX2*tJVs*tYJ3+tau3*tXJ2*tYJ3) &
              -FIVE*ONEFOURTH*odd*odd*tJVs*tXJ2*tYJ3 &
              -ONEHALF*(tX2*(s_dk*dYl3)+tau3*tXYJ23+tY3*(s_dk*dXl2) &
                       +tY3*(s_dl*kXd2)+tX2*(s_dl*kYd3)+tau3*tXYJ23) &
              +s_dd*tXYV23+s_dd*tYXV23 &
              +ONEHALF*odd*(tXJ2*(s_dk*dYl3)+tJVs*tXYJ23+tYJ3*(s_dk*dXl2) &
                           +tYJ3*(s_dl*kXd2)+tXJ2*(s_dl*kYd3)+tJVs*tXYJ23) &
              -ONETHIRD*(kYd3*dXl2+kXYd23*s_dl+kYXd23*s_dl+s_dk*dXYl23 &
                        +s_dk*dYXl23+kXd2*dYl3) )
        temp3=oc*(THREE*s_dd*tY3-tYJ3)
        temp4=ONETHIRD*oc*(THREE*s_dd*tX4-tXJ4)
        temp5=ONETHIRD*THREE*oc*( THREEHALF*s_dd*tY3*tY5-ONEHALF*(tY5*tYJ3+tY3*tYJ5) &
              +s_dd*tXY35-ONETHIRD*(tXYJ35+tXYJ35)+ONEHALF*odd*tYJ3*tYJ5 )
        temp6=ONETHIRD*THREE*oc*( THREEHALF*s_dd*tX2*tY7-ONEHALF*(tY7*tXJ2+tX2*tYJ7) &
              +s_dd*tXY27-ONETHIRD*(tXYJ27+tXYJ27)+ONEHALF*odd*tXJ2*tYJ7 )
        temp7=-12*tr1*rmkl(i,j)+4*temp1-8*temp2-2*tr4*temp3-2*temp4+4*temp5+4*temp6
        OOkl=OOkl+&
              temp7*Glob_ScaledPseudoChargeMatrix(i,j)/(Glob_Mass(i+1)*Glob_Mass(j+1))
      enddo
    enddo
    OOkl=OOkl/2

!Evaluation of correlation functions
    if (AreCorrFuncNeeded) then
      temp1=Skl/(Glob_Pi*Glob_SqrtPi)
      p=0
      do i=1,nn
        do j=i,nn
          p=p+1
          temp2=temp1/(sqrtTrAJ(j,i)*TrAJ(j,i))
          temp3=-1/TrAJ(j,i)
          temp4=2/TrAJ(j,i)
          temp5=eta2(j,i)/(TrAJ(j,i)*tau3)
          do k=1,NumCFGridPoints
            temp6=CFGrid(2,k)*CFGrid(2,k)         !this is \xi_z^2
            temp7=temp6+CFGrid(1,k)*CFGrid(1,k)   !this is  \xi^2
            temp8=ONE+(temp4*temp6-ONE)*temp5
            CFkl(p,k)=temp2*temp8*exp(temp7*temp3)
          enddo
        enddo
      enddo
    endif

    if (ArePartDensNeeded) then
      temp1=Skl/(Glob_Pi*Glob_SqrtPi)
      do i=1,nn+1
        temp2=ZERO
        do p=1,nn
          temp2=temp2+Glob_bvc(p,i)*Glob_bvc(p,i)*inv_tAkl(p,p)
          do q=p+1,nn
            temp2=temp2+2*Glob_bvc(q,i)*Glob_bvc(p,i)*inv_tAkl(q,p)
          enddo
        enddo
        temp3=ZERO
        temp4=ZERO
        do p=1,nn
          temp3=temp3+tvkinv_tAkl(p)*Glob_bvc(p,i)
          temp4=temp4+Glob_bvc(p,i)*inv_tAkltvl(p)
        enddo
        temp5=temp3*temp4/(temp2*tau3)
        temp6=-1/temp2
        temp7=2/temp2
        temp8=temp1/(sqrt(temp2)*temp2)
        do k=1,NumDensGridPoints
          temp9=DensGrid(2,k)*DensGrid(2,k)          !this is  \xi_z^2
          temp10=temp9+DensGrid(1,k)*DensGrid(1,k)   !this is -\xi^2
          temp11=ONE+(temp7*temp9-ONE)*temp5
          Denskl(i,k)=temp8*temp11*exp(temp10*temp6)
        enddo
      enddo
    endif

    if (AreMCorrFuncNeeded) then
      temp1=MSkl/(Glob_Pi*Glob_SqrtPi)
      p=0
      do i=1,nn
        do j=i,nn
          p=p+1
          temp2=temp1/(sqrtMTrAJ(j,i)*MTrAJ(j,i))/8
          temp3=-1/MTrAJ(j,i)/4
          temp4=1/MTrAJ(j,i)/2
          temp5=meta2(j,i)/(MTrAJ(j,i)*tau4)
          do k=1,NumCFGridPoints
            temp6=CFGrid(2,k)*CFGrid(2,k)         !this is \xi_z^2
            temp7=temp6+CFGrid(1,k)*CFGrid(1,k)   !this is  \xi^2
            temp8=ONE+(temp4*temp6-ONE)*temp5
            CFkl(p,k)=temp2*temp8*exp(temp7*temp3)
          enddo
        enddo
      enddo
    endif

    if (AreMPartDensNeeded) then
      temp1=MSkl/(Glob_Pi*Glob_SqrtPi)
      do i=1,nn+1
        temp2=ZERO
        do p=1,nn
          temp2=temp2+Glob_bvc(p,i)*Glob_bvc(p,i)*inv_invtAkinvtAl(p,p)
          do q=p+1,nn
            temp2=temp2+2*Glob_bvc(q,i)*Glob_bvc(p,i)*inv_invtAkinvtAl(q,p)
          enddo
        enddo
        temp3=ZERO
        temp4=ZERO
        do p=1,nn
          temp3=temp3+tvkinv_tAkinv_invtAkinvtAl(p)*Glob_bvc(p,i)
          temp4=temp4+Glob_bvc(p,i)*inv_invtAkinvtAlinv_tAltvl(p)
        enddo
        temp5=temp3*temp4/(temp2*tau4)
        temp6=-1/temp2
        temp7=2/temp2
        temp8=temp1/(sqrt(temp2)*temp2)
        do k=1,NumDensGridPoints
          temp9=DensGrid(2,k)*DensGrid(2,k)          !this is  \xi_z^2
          temp10=temp9+DensGrid(1,k)*DensGrid(1,k)   !this is -\xi^2
          temp11=ONE+(temp7*temp9/4-ONE)*temp5
          Denskl(i,k)=temp8*temp11*exp(temp10*temp6/4)/8
        enddo
      enddo
    endif

  end subroutine MatrixElementsAll_RG_1P

  subroutine symmetrize_matrix(W)
!subroutine symmetrize_matrix makes an arbitrary square matrix W
!symmetric by the following procedure:
!W = (1/2)*(W + W')
!Input:
!   W :: n x n real matrix

    integer, parameter :: nn = Glob_AllowedNumOfPseudoParticles
    real(wp)           W(nn, nn), t
    integer               i,j,n

    n = Glob_n

    do i = 1,nn
      do j = i+1,nn
        t=ONEHALF*(W(j,i)+W(i,j))
        W(j,i) = t
        W(i,j) = t
      end do
    end do

  end subroutine symmetrize_matrix

  function ME_rXr(X,inv_tAkl,inv_tAkltvl,tvkinv_tAkl,inv_tau3,Skl)
!function ME_rXr computes the following matrix element:
!<\tilde phi_k| r'Xr |\tilde phi_l>
!Here X is a symmetric real matrix. If matrix X is not symmetric
!then user need to symmetrize it before calling this function.
!Input:
!            X  :: n x n real matrix
!      inv_tAkl :: n x n real matrix where the inverse of tAk+tAl is stored
!   tvkinv_tAkl :: n-component vector where tvk*inv_tAkl is stored
!   inv_tAkltvl :: n-component vector where inv_tAkl*tvl is stored
!      inv_tau3 :: scalar, inv_tau3 = 1/tau3 = 1/tr[inv_tAkl*tvl*tvk']
!           Skl :: scalar, overlap Skl=<\tilde phi_k|\tilde phi_l>
    real(wp)   ME_rXr
    integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles
!Arguments:
    real(wp)   X(nn,nn),inv_tAkl(nn,nn),inv_tAkltvl(nn),tvkinv_tAkl(nn),inv_tau3,Skl
!Local variables:
    integer       i,j,n
    real(wp)   trAX,trAXAtvltvk,temp1,workvec1(nn)

    n=Glob_n
    trAX=ZERO
    do i=1,nn
      do j=1,nn
        trAX=trAX+inv_tAkl(i,j)*X(j,i)
      enddo
    enddo
    do i=1,nn
      temp1=ZERO
      do j=1,nn
        temp1=temp1+tvkinv_tAkl(j)*X(j,i)
      enddo
      workvec1(i)=temp1
    enddo
    trAXAtvltvk=ZERO
    do i=1,nn
      trAXAtvltvk=trAXAtvltvk+workvec1(i)*inv_tAkltvl(i)
    enddo

    ME_rXr=Skl*(THREEHALF*trAX+trAXAtvltvk*inv_tau3)

  end function ME_rXr

  function ME_rXr_rYr(X,Y,inv_tAkl,inv_tAkltvl,tvkinv_tAkl,inv_tau3,Skl)
!function ME_rXr_rYr computes the following matrix element:
!<\tilde phi_k| r'Xr r'Yr |\tilde phi_l>
!Here X and Y are symmetric real matrices. If matrices are not symmetric
!then user needs to symmetrize them before calling this function.
!Input:
!            X  :: n x n real matrix
!            Y  :: n x n real matrix
!      inv_tAkl :: n x n real matrix where the inverse of tAk+tAl is stored
!   tvkinv_tAkl :: n-component vector where tvk*inv_tAkl is stored
!   inv_tAkltvl :: n-component vector where inv_tAkl*tvl is stored
!      inv_tau3 :: scalar, inv_tau3 = 1/tau3 = 1/tr[inv_tAkl*tvl*tvk']
!           Skl :: scalar, overlap Skl=<\tilde phi_k|\tilde phi_l>
    real(wp)   ME_rXr_rYr
    integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles
!Arguments:
    real(wp)   X(nn,nn),Y(nn,nn),inv_tAkl(nn,nn),inv_tAkltvl(nn),tvkinv_tAkl(nn),inv_tau3,Skl
!Local variables:
    integer       i,j,k,n
    real(wp)   trAX,trAY,trAXAY,temp1,temp2
    real(wp)   trAXAtvltvk,trAYAtvltvk,trAXAYAtvltvk,trAYAXAtvltvk
    real(wp)   workvec1(nn),workvec2(nn),workvec3(nn),workvec4(nn)
    real(wp)   AX(nn,nn),AY(nn,nn)

    n=Glob_n
    trAX=ZERO
    trAY=ZERO
    do i=1,nn
      do j=1,nn
        temp1=ZERO
        temp2=ZERO
        do k=1,nn
          temp1=temp1+inv_tAkl(j,k)*X(k,i)
          temp2=temp2+inv_tAkl(j,k)*Y(k,i)
        enddo
        AX(j,i)=temp1
        AY(j,i)=temp2
      enddo
      trAX=trAX+AX(i,i)
      trAY=trAY+AY(i,i)
    enddo
    trAXAY=ZERO
    do i=1,nn
      do k=1,nn
        trAXAY=trAXAY+AX(i,k)*AY(k,i)
      enddo
    enddo

    do i=1,nn
      temp1=ZERO
      temp2=ZERO
      do j=1,nn
        temp1=temp1+tvkinv_tAkl(j)*X(j,i)
        temp2=temp2+tvkinv_tAkl(j)*Y(j,i)
      enddo
      workvec1(i)=temp1
      workvec2(i)=temp2
    enddo
    trAXAtvltvk=ZERO
    trAYAtvltvk=ZERO
    do i=1,nn
      trAXAtvltvk=trAXAtvltvk+workvec1(i)*inv_tAkltvl(i)
      trAYAtvltvk=trAYAtvltvk+workvec2(i)*inv_tAkltvl(i)
    enddo

    do i=1,nn
      temp1=ZERO
      temp2=ZERO
      do j=1,nn
        temp1=temp1+workvec1(j)*AY(j,i)
        temp2=temp2+workvec2(j)*AX(j,i)
      enddo
      workvec3(i)=temp1
      workvec4(i)=temp2
    enddo
    trAXAYAtvltvk=ZERO
    trAYAXAtvltvk=ZERO
    do i=1,nn
      trAXAYAtvltvk=trAXAYAtvltvk+workvec3(i)*inv_tAkltvl(i)
      trAYAXAtvltvk=trAYAXAtvltvk+workvec4(i)*inv_tAkltvl(i)
    enddo

    ME_rXr_rYr=Skl*( THREEHALF*THREEHALF*trAX*trAY + THREEHALF*trAXAY+    &
                    (THREEHALF*(trAX*trAYAtvltvk+trAY*trAXAtvltvk)+trAXAYAtvltvk+trAYAXAtvltvk)*inv_tau3 )

  end function ME_rXr_rYr

  function ME_dWd2(W,tAk,tAl,inv_tAkl,tvk,tvl,inv_tAkltvl,tvkinv_tAkl,inv_tau3,Skl)
!function ME_dWd2 computes the following matrix element:
!<\tilde phi_k| (\nabla_r' W \nabla_r)^2 |\tilde phi_l>
!Here W is a real symmetric matrix. If matrix W is not symmetric
!then user needs to symmetrize it before calling this function.
!Input:
!            W  :: n x n real matrix
!           tAk :: n x n real matrix where \tilde Ak is stored
!           tAl :: n x n real matrix where \tilde Al is stored
!      inv_tAkl :: n x n real matrix where the inverse of tAk+tAl is stored
!   tvkinv_tAkl :: n-component vector where tvk*inv_tAkl is stored
!   inv_tAkltvl :: n-component vector where inv_tAkl*tvl is stored
!      inv_tau3 :: scalar, inv_tau3 = 1/tau3 = 1/tr[inv_tAkl*tvl*tvk']
!           Skl :: scalar, overlap Skl=<\tilde phi_k|\tilde phi_l>
    real(wp)   ME_dWd2
    integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles
!Arguments:
    real(wp)   W(nn,nn),tAk(nn,nn),tAl(nn,nn),inv_tAkl(nn,nn)
    integer       tvk(nn),tvl(nn)
    real(wp)   inv_tAkltvl(nn),tvkinv_tAkl(nn),inv_tau3,Skl
!Local variables:
    integer       i,j,k,n
    real(wp)   temp1,temp2,tuk(nn),tul(nn),tAkWtAk(nn,nn),tAlWtAl(nn,nn)
    real(wp)   WorkMat1(nn,nn),WorkMat2(nn,nn)
    real(wp)   workvec1(nn),workvec2(nn),workvec3(nn),workvec4(nn)
    real(wp)   trWorkMat1,trWorkMat2,trAtAkWtAk,trAtAlWtAl,trAtvltuk,trAtultvk
    real(wp)   trAtAkWtAkAtultvk,trAtAlWtAlAtvltuk,trAtultuk,temp3

    n=Glob_n

!Compute  WorkMat1=W*tAk  WorkMat2=W*tAl  and their traces
    trWorkMat1=ZERO
    trWorkMat2=ZERO
    do i=1,nn
      do j=1,nn
        temp1=ZERO
        temp2=ZERO
        do k=1,nn
          temp1=temp1+W(j,k)*tAk(k,i)            !  tAk(k,i)*W(j,k)
          temp2=temp2+W(j,k)*tAl(k,i)            !  tAl(k,i)*W(j,k)
        enddo
        WorkMat1(j,i)=temp1
        WorkMat2(j,i)=temp2
      enddo
      trWorkMat1=trWorkMat1+WorkMat1(i,i)
      trWorkMat2=trWorkMat2+WorkMat2(i,i)
    enddo

!Compute tuk' = 4*tvk'*W*tAk + 6*tr[W*tAk]*tvk'
!        tul' = 4*tvl'*W*tAl + 6*tr[W*tAl]*tvl'
    do i=1,nn
      temp1=ZERO
      temp2=ZERO
      do j=1,nn
        temp1=temp1+tvk(j)*WorkMat1(j,i)
        temp2=temp2+tvl(j)*WorkMat2(j,i)
      enddo
      tuk(i)=FOUR*temp1+trWorkMat1*tvk(i)*6
      tul(i)=FOUR*temp2+trWorkMat2*tvl(i)*6
    enddo

!Compute tAkWtAk=tAk*W*tAk  tAlWtAl=tAl*W*tAl
    do i=1,nn
      do j=1,i
        temp1=ZERO
        temp2=ZERO
        do k=1,nn
          temp1=temp1+tAk(i,k)*WorkMat1(k,j)
          temp2=temp2+tAl(i,k)*WorkMat2(k,j)
        enddo
        tAkWtAk(i,j)=temp1
        tAkWtAk(j,i)=temp1
        tAlWtAl(i,j)=temp2
        tAlWtAl(j,i)=temp2
      enddo
    enddo

!Compute trAtAkWtAk=tr[inv_tAkl*tAk*W*tAk]
!        trAtAlWtAl=tr[inv_tAkl*tAl*W*tAl]
    trAtAkWtAk=ZERO
    trAtAlWtAl=ZERO
    do i=1,nn
      do j=1,nn
        trAtAkWtAk=trAtAkWtAk+inv_tAkl(i,j)*tAkWtAk(j,i)
        trAtAlWtAl=trAtAlWtAl+inv_tAkl(i,j)*tAlWtAl(j,i)
      enddo
    enddo

!Compute trAtvltuk=tr[inv_tAkl*tvl*tuk], trAtultvk=tr[inv_tAkl*tul*tvk]
    trAtvltuk=ZERO
    trAtultvk=ZERO
    do i=1,nn
      trAtvltuk=trAtvltuk+tuk(i)*inv_tAkltvl(i)
      trAtultvk=trAtultvk+tvkinv_tAkl(i)*tul(i)
    enddo

!Compute
!trAtAkWtAkAtultvk = tr[inv_tAkl*tAkWtAk*inv_tAkl*tul*tvk'] = tvkinv_tAkl'*tAkWtAk*inv_tAkl*tul
!trAtAlWtAlAtvltuk = tr[inv_tAkl*tAlWtAl*inv_tAkl*tvl*tuk'] = tuk'*inv_tAkl*tAlWtAl*inv_tAkltvl
    do i=1,nn
      temp1=ZERO
      temp2=ZERO
      do j=1,nn
        temp1=temp1+tvkinv_tAkl(j)*tAkWtAk(j,i)
        temp2=temp2+tuk(j)*inv_tAkl(j,i)
      enddo
      workvec1(i)=temp1
      workvec2(i)=temp2
    enddo
    do i=1,nn
      temp1=ZERO
      temp2=ZERO
      do j=1,nn
        temp1=temp1+workvec1(j)*inv_tAkl(j,i)
        temp2=temp2+workvec2(j)*tAlWtAl(j,i)
      enddo
      workvec3(i)=temp1
      workvec4(i)=temp2
    enddo
    trAtAkWtAkAtultvk=ZERO
    trAtAlWtAlAtvltuk=ZERO
    do i=1,nn
      trAtAkWtAkAtultvk=trAtAkWtAkAtultvk+workvec3(i)*tul(i)
      trAtAlWtAlAtvltuk=trAtAlWtAlAtvltuk+workvec4(i)*inv_tAkltvl(i)
    enddo

!Compute trAtultuk=tr[inv_tAkl*tul*tuk']
    trAtultuk=ZERO
    do i=1,nn
      trAtultuk=trAtultuk+workvec2(i)*tul(i)
    enddo

    ME_dWd2= 16*ME_rXr_rYr(tAkWtAk,tAlWtAl,inv_tAkl,inv_tAkltvl,tvkinv_tAkl,inv_tau3,Skl)  &
             +Skl*inv_tau3*(trAtultuk-SIX*(trAtAkWtAk*trAtultvk +trAtAlWtAl*trAtvltuk) &
                            -FOUR*(trAtAkWtAkAtultvk+trAtAlWtAlAtvltuk))

  end function ME_dWd2

  function ME_dWd21(X,Glob_B,tAk,tAl,inv_tAkl,tvk,tvl,inv_tAkltvl,tvkinv_tAkl,inv_tau3,Skl)
!function ME_dWd2 computes the following matrix element:
!<\tilde phi_k| (\nabla_r' X \nabla_r)^2 |\tilde phi_l>
!Here X is a real symmetric matrix. If matrix X is not symmetric
!then user needs to symmetrize it before calling this function.
!Input:
!             X :: n x n real matrix (Glob_MassMatrix)
!        Glob_B :: n x n real matrix
!           tAk :: n x n real matrix where \tilde Ak is stored
!           tAl :: n x n real matrix where \tilde Al is stored
!      inv_tAkl :: n x n real matrix where the inverse of tAk+tAl is stored
!   tvkinv_tAkl :: n-component vector where tvk*inv_tAkl is stored
!   inv_tAkltvl :: n-component vector where inv_tAkl*tvl is stored
!      inv_tau3 :: scalar, inv_tau3 = 1/tau3 = 1/tr[inv_tAkl*tvl*tvk']
!           Skl :: scalar, overlap Skl=<\tilde phi_k|\tilde phi_l>
    real(wp)   ME_dWd21
    integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles
!Arguments:
    real(wp)   X(nn,nn),Glob_B(nn,nn),tAk(nn,nn),tAl(nn,nn),inv_tAkl(nn,nn)
    integer       tvk(nn),tvl(nn)
    real(wp)   inv_tAkltvl(nn),tvkinv_tAkl(nn),inv_tau3,Skl
!Local variables:
    integer       i,j,k,n
    real(wp)   temp1,temp2,tuk(nn),tul(nn),tAkXtAk(nn,nn),tAlXtAl(nn,nn)
    real(wp)   WorkMat1(nn,nn),WorkMat2(nn,nn)
    real(wp)   workvec1(nn),workvec2(nn),workvec3(nn),workvec4(nn)
    real(wp)   trWorkMat1,trWorkMat2,trAtAkXtAk,trAtAlXtAl,trAtvltuk,trAtultvk
    real(wp)   trAtAkXtAkAtultvk,trAtAlXtAlAtvltuk,trAtultuk,temp3

    n=Glob_n

!Compute  WorkMat1=X*tAk  WorkMat2=X*tAl  and their traces
    trWorkMat1=ZERO
    trWorkMat2=ZERO
    do i=1,nn
      do j=1,nn
        temp1=ZERO
        temp2=ZERO
        do k=1,nn
          temp1=temp1+X(j,k)*tAk(k,i)
          temp2=temp2+Glob_B(j,k)*tAl(k,i)
        enddo
        WorkMat1(j,i)=temp1
        WorkMat2(j,i)=temp2
      enddo
      trWorkMat1=trWorkMat1+WorkMat1(i,i)
      trWorkMat2=trWorkMat2+WorkMat2(i,i)
    enddo

!Compute tuk' = 4*tvk'*X*tAk + 6*tr[X*tAk]*tvk'
!        tul' = 4*tvl'*X*tAl + 6*tr[X*tAl]*tvl'

    do i=1,nn
      temp1=ZERO
      temp2=ZERO
      do j=1,nn
        temp1=temp1+tvk(j)*WorkMat1(j,i)
        temp2=temp2+tvl(j)*WorkMat2(j,i)
      enddo
      tuk(i)=FOUR*temp1+trWorkMat1*tvk(i)*6
      tul(i)=FOUR*temp2+trWorkMat2*tvl(i)*6
    enddo

!Compute tAkXtAk=tAk*X*tAk  tAlWtAl=tAl*X*tAl
    do i=1,nn
      do j=1,i
        temp1=ZERO
        temp2=ZERO
        do k=1,nn
          temp1=temp1+tAk(i,k)*WorkMat1(k,j)
          temp2=temp2+tAl(i,k)*WorkMat2(k,j)
        enddo
        tAkXtAk(i,j)=temp1
        tAkXtAk(j,i)=temp1
        tAlXtAl(i,j)=temp2
        tAlXtAl(j,i)=temp2
      enddo
    enddo

!Compute trAtAkXtAk=tr[inv_tAkl*tAk*X*tAk]
!        trAtAlXtAl=tr[inv_tAkl*tAl*X*tAl]
    trAtAkXtAk=ZERO
    trAtAlXtAl=ZERO
    do i=1,nn
      do j=1,nn
        trAtAkXtAk=trAtAkXtAk+inv_tAkl(i,j)*tAkXtAk(j,i)
        trAtAlXtAl=trAtAlXtAl+inv_tAkl(i,j)*tAlXtAl(j,i)
      enddo
    enddo

!Compute trAtvltuk=tr[inv_tAkl*tvl*tuk], trAtultvk=tr[inv_tAkl*tul*tvk]
    trAtvltuk=ZERO
    trAtultvk=ZERO
    do i=1,nn
      trAtvltuk=trAtvltuk+tuk(i)*inv_tAkltvl(i)
      trAtultvk=trAtultvk+tvkinv_tAkl(i)*tul(i)
    enddo

!Compute
!trAtAkXtAkAtultvk = tr[inv_tAkl*tAkXtAk*inv_tAkl*tul*tvk'] = tvkinv_tAkl'*tAkXtAk*inv_tAkl*tul
!trAtAlXtAlAtvltuk = tr[inv_tAkl*tAlXtAl*inv_tAkl*tvl*tuk'] = tuk'*inv_tAkl*tAlXtAl*inv_tAkltvl
    do i=1,nn
      temp1=ZERO
      temp2=ZERO
      do j=1,nn
        temp1=temp1+tvkinv_tAkl(j)*tAkXtAk(j,i)
        temp2=temp2+tuk(j)*inv_tAkl(j,i)
      enddo
      workvec1(i)=temp1
      workvec2(i)=temp2
    enddo
    do i=1,nn
      temp1=ZERO
      temp2=ZERO
      do j=1,nn
        temp1=temp1+workvec1(j)*inv_tAkl(j,i)
        temp2=temp2+workvec2(j)*tAlXtAl(j,i)
      enddo
      workvec3(i)=temp1
      workvec4(i)=temp2
    enddo
    trAtAkXtAkAtultvk=ZERO
    trAtAlXtAlAtvltuk=ZERO
    do i=1,nn
      trAtAkXtAkAtultvk=trAtAkXtAkAtultvk+workvec3(i)*tul(i)
      trAtAlXtAlAtvltuk=trAtAlXtAlAtvltuk+workvec4(i)*inv_tAkltvl(i)
    enddo

!Compute trAtultuk=tr[inv_tAkl*tul*tuk']
    trAtultuk=ZERO
    do i=1,nn
      trAtultuk=trAtultuk+workvec2(i)*tul(i)
    enddo

    ME_dWd21= 16*ME_rXr_rYr(tAkXtAk,tAlXtAl,inv_tAkl,inv_tAkltvl,tvkinv_tAkl,inv_tau3,Skl)  &
              +Skl*inv_tau3*(trAtultuk-SIX*(trAtAkXtAk*trAtultvk +trAtAlXtAl*trAtvltuk) &
                             -FOUR*(trAtAkXtAkAtultvk+trAtAlXtAlAtvltuk))

  end function ME_dWd21

  function ME_dXd(X,tvk,tvl,inv_tAkltvl,inv_tAkl,tAk,tAl,inv_tAkltAl,Skl,tau3)
    real(wp)   ME_dXd
    integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles
    real(wp)   X(nn,nn),tAk(nn,nn),tAl(nn,nn),inv_tAkl(nn,nn),inv_tAkltAl(nn,nn)
    integer       i,j,n,k,tvk(nn),tvl(nn)
    real(wp)   inv_tAkltAlX(nn,nn),inv_tAkltAlXtAk(nn,nn),tvkinv_tAkltAlX(nn),inv_tAkltvl(nn)
    real(wp)   temp1, Skl,tau,tau1,tau2,tau3
    n=Glob_n
    do i=1,nn
      do j=1,nn
        temp1=ZERO
        do k=1,nn
          temp1=temp1+inv_tAkltAl(j,k)*X(k,i)
        enddo
        inv_tAkltAlX(j,i)=temp1
      enddo
    enddo
    tau1=ZERO
    do i=1,nn
      temp1=ZERO
      do k=1,nn
        temp1=temp1+inv_tAkltAlX(i,k)*tAk(k,i)
      enddo
      tau1=tau1+temp1
    enddo
    do i=1,nn
      temp1=ZERO
      do j=1,nn
        temp1=temp1+tvk(j)*inv_tAkltAlX(j,i)
      enddo
      tvkinv_tAkltAlX(i)=temp1
    enddo
    tau2=ZERO
    do i=1,nn
      temp1=ZERO
      do j=1,nn
        temp1=temp1+tvkinv_tAkltAlX(j)*tAk(j,i)
      enddo
      tau2=tau2+temp1*inv_tAkltvl(i)
    enddo
    ME_dXd=Skl*(SIX*tau1+FOUR*tau2/tau3)
  end function ME_dXd

  function SG_ME_rXr_over_rij(i,j,X,inv_tAkl,t_V,Skl)
!function SG_ME_rXr_over_rij computes the following matrix element:
!<\tilde psi_k| (r' X r)/r_ij |\tilde psi_l>
!where psi_k = exp[-(r' Ak r)] is a Simple Gaussian wavefunction
!Here X is a some real symmetric matrix. If matrix X is not symmetric
!then user needs to symmetrize it before calling this function.
!Input:
!            X  :: n x n real matrix
!      inv_tAkl :: n x n real matrix where the inverse of tAk+tAl is stored
!           t_V :: scalar, t_V = tau3 = tr[inv_tAkl*tvl*tvk']
!           Skl :: scalar, overlap Skl=<\tilde phi_k|\tilde phi_l>
    real(wp)   SG_ME_rXr_over_rij
    integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles
!Arguments:
    real(wp)   X(nn,nn),inv_tAkl(nn,nn)
    integer       i,j
    real(wp)   t_V,Skl
!Local variables:
    integer       p,q,n
    real(wp)   temp1,temp2,temp3
    real(wp)   Aj(nn),AjX(nn)
    real(wp)   t_J,t_X,t_XJ

    n=Glob_n

!Form Aj=inv_tAkl*ji        j/=i
!     Aj=inv_tAkl*(ji-jj)   j/=i
!Remember that Jii=ji*ji' and Jij=(ji-jj)*(ji-jj)'
    if (i==j) then
      do p=1,nn
        Aj(p)=inv_tAkl(p,i)
      enddo
    else
      do p=1,nn
        Aj(p)=inv_tAkl(p,i)-inv_tAkl(p,j)
      enddo
    endif

!Compute AjX'=Aj'*X
    do p=1,nn
      temp1=ZERO
      do q=1,nn
        temp1=temp1+Aj(q)*X(q,p)
      enddo
      AjX(p)=temp1
    enddo

!Compute t_J=tr[inv_tAkl*Jij]
    if (i==j) then
      t_J=inv_tAkl(i,i)
    else
      t_J=inv_tAkl(i,i)+inv_tAkl(j,j)-inv_tAkl(j,i)-inv_tAkl(j,i)
    endif

!Compute t_XJ=tr[inv_tAkl*X*inv_tAkl*Jij]=AjX'*Aj
    t_XJ=ZERO
    do p=1,nn
      t_XJ=t_XJ+AjX(p)*Aj(p)
    enddo

!Compute t_X=tr[inv_tAkl*X]
    t_X=ZERO
    do p=1,nn
      do q=1,nn
        t_X=t_X+inv_tAkl(q,p)*X(q,p)
      enddo
    enddo

    temp1=(TWO/Glob_SqrtPi)*Skl
    temp2=1/t_V
    temp3=1/t_J
    SG_ME_rXr_over_rij=temp1*temp2*temp3*sqrt(temp3)*(THREE*t_J*t_X - t_XJ)

  end function SG_ME_rXr_over_rij

  function SG_ME_rXr_rYr_over_rij(i,j,X,Y,inv_tAkl,t_V,Skl)
!function ME_rXr_rYr_over_rij computes the following matrix element:
!<\tilde phi_k| (r' X r)(r' Y r)/r_ij |\tilde phi_l>
!Here X and Y are some real symmetric matrices. If matrices X or Y are not symmetric
!then user needs to symmetrize them before calling this function.
!Input:
!            X  :: n x n real matrix
!            Y  :: n x n real matrix
!      inv_tAkl :: n x n real matrix where the inverse of tAk+tAl is stored
!           t_V :: scalar, t_V = tau3 = tr[inv_tAkl*tvl*tvk']
!           Skl :: scalar, overlap Skl=<\tilde phi_k|\tilde phi_l>
    real(wp)   SG_ME_rXr_rYr_over_rij
    integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles
!Arguments:
    real(wp)   X(nn,nn),Y(nn,nn),inv_tAkl(nn,nn)
    integer       i,j
    real(wp)   t_V,Skl
!Local variables:
    integer       p,q,s,n
    real(wp)   temp1,temp2,temp3,temp4,temp5,temp6
    real(wp)   AX(nn,nn),AY(nn,nn)
    real(wp)   Aj(nn),AjX(nn),AjY(nn),AXAj(nn),AYAj(nn)
    real(wp)   t_J,t_X,t_Y
    real(wp)   t_XJ,t_YJ,t_XY
    real(wp)   t_XYJ,t_YXJ

    n=Glob_n

!Form Aj=inv_tAkl*ji        j/=i
!     Aj=inv_tAkl*(ji-jj)   j/=i
!Remember that Jii=ji*ji' and Jij=(ji-jj)*(ji-jj)'
    if (i==j) then
      do p=1,nn
        Aj(p)=inv_tAkl(p,i)
      enddo
    else
      do p=1,nn
        Aj(p)=inv_tAkl(p,i)-inv_tAkl(p,j)
      enddo
    endif

!Compute AjX'=Aj'*X
!    and AjY'=Aj'*Y
    do p=1,nn
      temp1=ZERO
      temp2=ZERO
      do q=1,nn
        temp1=temp1+Aj(q)*X(q,p)
        temp2=temp2+Aj(q)*Y(q,p)
      enddo
      AjX(p)=temp1
      AjY(p)=temp2
    enddo

!Compute AX=inv_tAkl*X  t_X=tr[inv_tAkl*X]
!        AY=inv_tAkl*Y  t_Y=tr[inv_tAkl*Y]
    t_X=ZERO
    t_Y=ZERO
    do p=1,nn
      do q=1,nn
        temp1=ZERO
        temp2=ZERO
        do s=1,nn
          temp1=temp1+inv_tAkl(s,q)*X(p,s)
          temp2=temp2+inv_tAkl(s,q)*Y(p,s)
        enddo
        AX(q,p)=temp1
        AY(q,p)=temp2
      enddo
      t_X=t_X+AX(p,p)
      t_Y=t_Y+AY(p,p)
    enddo

!Compute t_J=tr[inv_tAkl*Jij]
    if (i==j) then
      t_J=inv_tAkl(i,i)
    else
      t_J=inv_tAkl(i,i)+inv_tAkl(j,j)-inv_tAkl(j,i)-inv_tAkl(j,i)
    endif

!Compute t_XY=tr[inv_tAkl*X*inv_tAkl*Y]=tr[AX*AY]
!        AXAj=AX*Aj
!        AYAj=AY*Aj
    t_XY=ZERO
    do p=1,nn
      temp3=ZERO
      temp4=ZERO
      do q=1,nn
        t_XY=t_XY+AX(p,q)*AY(q,p)
        temp3=temp3+AX(p,q)*Aj(q)
        temp4=temp4+AY(p,q)*Aj(q)
      enddo
      AXAj(p)=temp3
      AYAj(p)=temp4
    enddo

!Compute
!t_XJ=tr[inv_tAkl*X*inv_tAkl*Jij]=AjX'*Aj
!t_YJ=tr[inv_tAkl*Y*inv_tAkl*Jij]=AjY'*Aj
!t_XYJ=tr[inv_tAkl*X*inv_tAkl*Y*inv_tAkl*Jij]=AjX'*AYAj

    t_XJ=ZERO
    t_YJ=ZERO
    t_XYJ=ZERO

    do p=1,nn
      t_XJ=t_XJ+AjX(p)*Aj(p)
      t_YJ=t_YJ+AjY(p)*Aj(p)
      t_XYJ=t_XYJ+AjX(p)*AYAj(p)
    enddo

!Compute t_YXJ=tr[inv_tAkl*Y*inv_tAkl*X*inv_tAkl*Jij]
    t_YXJ=t_XYJ

    temp1=(TWO/Glob_SqrtPi)*Skl
    temp2=1/t_V
    temp3=1/t_J
    SG_ME_rXr_rYr_over_rij=THREE*temp1*temp2*temp3*sqrt(temp3)*(  &
                            THREEHALF*t_J*t_X*t_Y - ONEHALF*(t_Y*t_XJ + t_X*t_YJ) + &
                            t_J*t_XY  - ONETHIRD*(t_XYJ + t_YXJ) + ONEHALF*temp3*t_XJ*t_YJ&
                            )

  end function SG_ME_rXr_rYr_over_rij

  function ME_rXr_over_rij(i,j,X,inv_tAkl,tvk,tvl,inv_tAkltvl,tvkinv_tAkl,t_V,Skl)
!function ME_rXr_over_rij computes the following matrix element:
!<\tilde phi_k| (r' X r)/r_ij |\tilde phi_l>
!Here X is a some real symmetric matrix. If matrix X is not symmetric
!then user needs to symmetrize it before calling this function.
!Input:
!            X  :: n x n real matrix
!      inv_tAkl :: n x n real matrix where the inverse of tAk+tAl is stored
!   tvkinv_tAkl :: n-component vector where tvk*inv_tAkl is stored
!   inv_tAkltvl :: n-component vector where inv_tAkl*tvl is stored
!           t_V :: scalar, t_V = tau3 = tr[inv_tAkl*tvl*tvk']
!           Skl :: scalar, overlap Skl=<\tilde phi_k|\tilde phi_l>
    real(wp)   ME_rXr_over_rij
    integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles
!Arguments:
    real(wp)   X(nn,nn),inv_tAkl(nn,nn)
    integer       i,j,tvk(nn),tvl(nn)
    real(wp)   inv_tAkltvl(nn),tvkinv_tAkl(nn),t_V,Skl
!Local variables:
    integer       p,q,n
    real(wp)   Ajtvl,temp1,temp2,temp3
    real(wp)   Aj(nn),AjX(nn)
    real(wp)   t_J,t_X,t_JV,t_XJ,t_XV,t_JXV,t_XJV

    n=Glob_n

!Form Aj=inv_tAkl*ji        j/=i
!     Aj=inv_tAkl*(ji-jj)   j/=i
!Remember that Jii=ji*ji' and Jij=(ji-jj)*(ji-jj)'
    if (i==j) then
      do p=1,nn
        Aj(p)=inv_tAkl(p,i)
      enddo
    else
      do p=1,nn
        Aj(p)=inv_tAkl(p,i)-inv_tAkl(p,j)
      enddo
    endif

!Compute AjX'=Aj'*X
    do p=1,nn
      temp1=ZERO
      do q=1,nn
        temp1=temp1+Aj(q)*X(q,p)
      enddo
      AjX(p)=temp1
    enddo

!Compute t_J=tr[inv_tAkl*Jij]
    if (i==j) then
      t_J=inv_tAkl(i,i)
    else
      t_J=inv_tAkl(i,i)+inv_tAkl(j,j)-inv_tAkl(j,i)-inv_tAkl(j,i)
    endif

!Compute Ajtvl=Aj'*tvl
!        t_XJ=tr[inv_tAkl*X*inv_tAkl*Jij]=AjX'*Aj
!        t_JV=tr[inv_tAkl*Jij*inv_tAkl*tvl*tvk']=tvk'*Aj*Ajtvl
!        t_JXV=tr[inv_tAkl*Jij*inv_tAkl*X*inv_tAkl*tvl*tvk']=tvk'*Aj*AjX'*inv_tAkltvl
    Ajtvl=ZERO
    t_XJ=ZERO
    temp1=ZERO
    temp2=ZERO
    do p=1,nn
      Ajtvl=Ajtvl+Aj(p)*tvl(p)
      t_XJ=t_XJ+AjX(p)*Aj(p)
      temp1=temp1+tvk(p)*Aj(p)
      temp2=temp2+AjX(p)*inv_tAkltvl(p)
    enddo
    t_JV=temp1*Ajtvl
    t_JXV=temp1*temp2

!Compute t_X=tr[inv_tAkl*X]
!        t_XV=tr[inv_tAkl*X*inv_tAkl*tvl*tvk']=tvkinv_tAkl'*X*inv_tAkltvl
!        t_XJV=tr[inv_tAkl*X*inv_tAkl*Jij*inv_tAkl*tvl*tvk']=tvkinv_tAkl'*X*Aj*Ajtvl
    t_X=ZERO
    t_XV=ZERO
    temp2=ZERO
    do p=1,nn
      temp1=ZERO
      do q=1,nn
        t_X=t_X+inv_tAkl(q,p)*X(q,p)
        temp1=temp1+tvkinv_tAkl(q)*X(q,p)
      enddo
      t_XV=t_XV+temp1*inv_tAkltvl(p)
      temp2=temp2+temp1*Aj(p)
    enddo
    t_XJV=temp2*Ajtvl

    temp1=(TWO/Glob_SqrtPi)*Skl
    temp2=1/t_V
    temp3=1/t_J
    ME_rXr_over_rij=temp1*temp2*temp3*sqrt(temp3)*( &
                     THREEHALF*t_J*t_V*t_X - ONEHALF*(t_X*t_JV + t_V*t_XJ) + t_J*t_XV + ONEHALF*temp3*t_JV*t_XJ  &
                     - ONETHIRD*(t_JXV + t_XJV) &
                     )
  end function ME_rXr_over_rij

  function ME_rXr_rYr_over_rij(i,j,X,Y,inv_tAkl,tvk,tvl,inv_tAkltvl,tvkinv_tAkl,t_V,Skl)
!function ME_rXr_rYr_over_rij computes the following matrix element:
!<\tilde phi_k| (r' X r)(r' Y r)/r_ij |\tilde phi_l>
!Here X and Y are some real symmetric matrices. If matrices X or Y are not symmetric
!then user needs to symmetrize them before calling this function.
!Input:
!            X  :: n x n real matrix
!            Y  :: n x n real matrix
!      inv_tAkl :: n x n real matrix where the inverse of tAk+tAl is stored
!   tvkinv_tAkl :: n-component vector where tvk*inv_tAkl is stored
!   inv_tAkltvl :: n-component vector where inv_tAkl*tvl is stored
!           t_V :: scalar, t_V = tau3 = tr[inv_tAkl*tvl*tvk']
!           Skl :: scalar, overlap Skl=<\tilde phi_k|\tilde phi_l>
    real(wp)   ME_rXr_rYr_over_rij
    integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles
!Arguments:
    real(wp)   X(nn,nn),Y(nn,nn),inv_tAkl(nn,nn)
    integer       i,j,tvk(nn),tvl(nn)
    real(wp)   inv_tAkltvl(nn),tvkinv_tAkl(nn),t_V,Skl
!Local variables:
    integer       p,q,s,n
    real(wp)   temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,temp9
    real(wp)   AX(nn,nn),AY(nn,nn)
    real(wp)   Aj(nn),AjX(nn),AjY(nn),AXAj(nn),AYAj(nn)
    real(wp)   AXinv_tAkltvl(nn),AYinv_tAkltvl(nn),tvkinv_tAklX(nn),tvkinv_tAklY(nn)
    real(wp)   Ajtvl,tvkAj
    real(wp)   t_J,t_X,t_Y
    real(wp)   t_JV,t_XJ,t_YJ,t_XV,t_YV,t_XY
    real(wp)   t_JXV,t_JYV,t_XJV,t_YJV,t_XYV,t_YXV,t_XYJ,t_YXJ
    real(wp)   t_YJXV,t_XJYV,t_XYJV,t_YXJV,t_JXYV,t_JYXV

    n=Glob_n

!Form Aj=inv_tAkl*ji        j/=i
!     Aj=inv_tAkl*(ji-jj)   j/=i
!Remember that Jii=ji*ji' and Jij=(ji-jj)*(ji-jj)'
    if (i==j) then
      do p=1,nn
        Aj(p)=inv_tAkl(p,i)
      enddo
    else
      do p=1,nn
        Aj(p)=inv_tAkl(p,i)-inv_tAkl(p,j)
      enddo
    endif

!Compute AjX'=Aj'*X
!    and AjY'=Aj'*Y
    do p=1,nn
      temp1=ZERO
      temp2=ZERO
      do q=1,nn
        temp1=temp1+Aj(q)*X(q,p)
        temp2=temp2+Aj(q)*Y(q,p)
      enddo
      AjX(p)=temp1
      AjY(p)=temp2
    enddo

!Compute AX=inv_tAkl*X  t_X=tr[inv_tAkl*X]
!        AY=inv_tAkl*Y  t_Y=tr[inv_tAkl*Y]
    t_X=ZERO
    t_Y=ZERO
    do p=1,nn
      do q=1,nn
        temp1=ZERO
        temp2=ZERO
        do s=1,nn
          temp1=temp1+inv_tAkl(s,q)*X(p,s)
          temp2=temp2+inv_tAkl(s,q)*Y(p,s)
        enddo
        AX(q,p)=temp1
        AY(q,p)=temp2
      enddo
      t_X=t_X+AX(p,p)
      t_Y=t_Y+AY(p,p)
    enddo

!Compute t_J=tr[inv_tAkl*Jij]
    if (i==j) then
      t_J=inv_tAkl(i,i)
    else
      t_J=inv_tAkl(i,i)+inv_tAkl(j,j)-inv_tAkl(j,i)-inv_tAkl(j,i)
    endif

!Compute t_XY=tr[inv_tAkl*X*inv_tAkl*Y]=tr[AX*AY]
!        AXinv_tAkltvl=AX*inv_tAkltvl
!        AYinv_tAkltvl=AY*inv_tAkltvl
!        AXAj=AX*Aj
!        AYAj=AY*Aj
!        tvkinv_tAklX'=tvkinv_tAkl'*X
!        tvkinv_tAklY'=tvkinv_tAkl'*Y
    t_XY=ZERO
    do p=1,nn
      temp1=ZERO
      temp2=ZERO
      temp3=ZERO
      temp4=ZERO
      temp5=ZERO
      temp6=ZERO
      do q=1,nn
        t_XY=t_XY+AX(p,q)*AY(q,p)
        temp1=temp1+AX(p,q)*inv_tAkltvl(q)
        temp2=temp2+AY(p,q)*inv_tAkltvl(q)
        temp3=temp3+AX(p,q)*Aj(q)
        temp4=temp4+AY(p,q)*Aj(q)
        temp5=temp5+tvkinv_tAkl(q)*X(q,p)
        temp6=temp6+tvkinv_tAkl(q)*Y(q,p)
      enddo
      AXinv_tAkltvl(p)=temp1
      AYinv_tAkltvl(p)=temp2
      AXAj(p)=temp3
      AYAj(p)=temp4
      tvkinv_tAklX(p)=temp5
      tvkinv_tAklY(p)=temp6
    enddo

!Compute
!Ajtvl=Aj'*tvl
!tvkAj=tvk'*Aj
!t_JV=tr[inv_tAkl*Jij*inv_tAkl*tvl*tvk']=tvk'*Aj*Ajtvl
!t_JXV=tr[inv_tAkl*Jij*inv_tAkl*X*inv_tAkl*tvl*tvk']=tvk'*Aj*AjX'*inv_tAkltvl
!t_JYV=tr[inv_tAkl*Jij*inv_tAkl*Y*inv_tAkl*tvl*tvk']=tvk'*Aj*AjY'*inv_tAkltvl
!t_XJ=tr[inv_tAkl*X*inv_tAkl*Jij]=AjX'*Aj
!t_YJ=tr[inv_tAkl*Y*inv_tAkl*Jij]=AjY'*Aj
!t_XYJ=tr[inv_tAkl*X*inv_tAkl*Y*inv_tAkl*Jij]=AjX'*AYAj
!t_XV=tr[inv_tAkl*X*inv_tAkl*tvl*tvk']=tvkinv_tAklX'*inv_tAkltvl
!t_YV=tr[inv_tAkl*Y*inv_tAkl*tvl*tvk']=tvkinv_tAklY'*inv_tAkltvl
!t_XYV=tr[inv_tAkl*X*inv_tAkl*Y*inv_tAkl*tvl*tvk']=tvkinv_tAklX'*AYinv_tAkltvl
!t_YXV=tr[inv_tAkl*Y*inv_tAkl*X*inv_tAkl*tvl*tvk']=tvkinv_tAklY'*AXinv_tAkltvl
!t_XJV=tr[inv_tAkl*X*inv_tAkl*Jij*inv_tAkl*tvl*tvk']=tvkinv_tAklX'*Aj*Ajtvl
!t_YJV=tr[inv_tAkl*Y*inv_tAkl*Jij*inv_tAkl*tvl*tvk']=tvkinv_tAklY'*Aj*Ajtvl
!t_XYJV=tr[inv_tAkl*X*inv_tAkl*Y*inv_tAkl*Jij*inv_tAkl*tvl*tvk']=tvkinv_tAklX'*AYAj*Ajtvl
!t_YXJV=tr[inv_tAkl*Y*inv_tAkl*X*inv_tAkl*Jij*inv_tAkl*tvl*tvk']=tvkinv_tAklY'*AXAj*Ajtvl
!t_XJYV=tr[inv_tAkl*X*inv_tAkl*Jij*inv_tAkl*Y*inv_tAkl*tvl*tvk']=tvkinv_tAklX'*Aj*AjY'*inv_tAkltvl
!t_YJXV=tr[inv_tAkl*Y*inv_tAkl*Jij*inv_tAkl*X*inv_tAkl*tvl*tvk']=tvkinv_tAklY'*Aj*AjX'*inv_tAkltvl
!t_JXYV=tr[inv_tAkl*Jij*inv_tAkl*X*inv_tAkl*Y*inv_tAkl*tvl*tvk']=tvkAj*AjX'*AYinv_tAkltvl
!t_JYXV=tr[inv_tAkl*Jij*inv_tAkl*Y*inv_tAkl*X*inv_tAkl*tvl*tvk']=tvkAj*AjY'*AXinv_tAkltvl
    Ajtvl=ZERO
    tvkAj=ZERO
    temp1=ZERO
    temp2=ZERO
    temp3=ZERO
    t_XJ=ZERO
    t_YJ=ZERO
    t_XYJ=ZERO
    t_XV=ZERO
    t_YV=ZERO
    t_XYV=ZERO
    t_YXV=ZERO
    temp4=ZERO
    temp5=ZERO
    temp6=ZERO
    temp7=ZERO
    temp8=ZERO
    temp9=ZERO
    do p=1,nn
      Ajtvl=Ajtvl+Aj(p)*tvl(p)
      tvkAj=tvkAj+tvk(p)*Aj(p)
      temp1=temp1+tvk(p)*Aj(p)
      temp2=temp2+AjX(p)*inv_tAkltvl(p)
      temp3=temp3+AjY(p)*inv_tAkltvl(p)
      t_XJ=t_XJ+AjX(p)*Aj(p)
      t_YJ=t_YJ+AjY(p)*Aj(p)
      t_XYJ=t_XYJ+AjX(p)*AYAj(p)
      t_XV=t_XV+tvkinv_tAklX(p)*inv_tAkltvl(p)
      t_YV=t_YV+tvkinv_tAklY(p)*inv_tAkltvl(p)
      t_XYV=t_XYV+tvkinv_tAklX(p)*AYinv_tAkltvl(p)
      t_YXV=t_YXV+tvkinv_tAklY(p)*AXinv_tAkltvl(p)
      temp4=temp4+tvkinv_tAklX(p)*Aj(p)
      temp5=temp5+tvkinv_tAklY(p)*Aj(p)
      temp6=temp6+tvkinv_tAklX(p)*AYAj(p)
      temp7=temp7+tvkinv_tAklY(p)*AXAj(p)
      temp8=temp8+AjX(p)*AYinv_tAkltvl(p)
      temp9=temp9+AjY(p)*AXinv_tAkltvl(p)
    enddo
    t_JV=temp1*Ajtvl
    t_JXV=temp1*temp2
    t_JYV=temp1*temp3
    t_XJV=temp4*Ajtvl
    t_YJV=temp5*Ajtvl
    t_XYJV=temp6*Ajtvl
    t_YXJV=temp7*Ajtvl
    t_XJYV=temp4*temp3
    t_YJXV=temp5*temp2
    t_JXYV=tvkAj*temp8
    t_JYXV=tvkAj*temp9

!Compute t_YXJ=tr[inv_tAkl*Y*inv_tAkl*X*inv_tAkl*Jij]
    t_YXJ=t_XYJ

    temp1=(TWO/Glob_SqrtPi)*Skl
    temp2=1/t_V
    temp3=1/t_J
    ME_rXr_rYr_over_rij=temp1*temp2*temp3*sqrt(temp3)*(  &
                         NINE*ONEFOURTH*t_J*t_V*t_X*t_Y   &
                         -THREE*ONEFOURTH*(t_X*t_Y*t_JV + t_V*t_Y*t_XJ + t_V*t_X*t_YJ)  &
                         +THREEHALF*(t_J*t_Y*t_XV + t_J*t_X*t_YV + t_J*t_V*t_XY)  &
                         -ONEHALF*(t_JV*t_XY + t_XV*t_YJ + t_XJ*t_YV)  &
                         +THREE*ONEFOURTH*temp3*(t_Y*t_JV*t_XJ+t_X*t_JV*t_YJ+t_V*t_XJ*t_YJ)  &
                         -FIVE*ONEFOURTH*temp3*temp3*t_JV*t_XJ*t_YJ   &
                         -ONEHALF*(t_X*t_JYV + t_V*t_XYJ + t_Y*t_JXV + t_Y*t_XJV + t_X*t_YJV + t_V*t_YXJ) + t_J*t_XYV + t_J*t_YXV  &
                         +ONEHALF*temp3*(t_XJ*t_JYV + t_JV*t_XYJ + t_YJ*t_JXV + t_YJ*t_XJV + t_XJ*t_YJV + t_JV*t_YXJ) &
                         -ONETHIRD*(t_YJXV + t_XYJV + t_YXJV + t_JXYV + t_JYXV + t_XJYV)  &
                         )
!    t_J    t_X    t_Y  (t_V is passed)
!    t_XJ   t_YJ   t_JV   t_XV   t_YV   t_XY
!    t_JXV  t_JYV  t_XJV  t_YJV  t_XYJ  t_YXJ=t_XYJ  t_XYV  t_YXV

  end function ME_rXr_rYr_over_rij

  function ME_d_X_over_rij_d(i,j,X,tAk,tAl,inv_tAkl,tvk,tvl,inv_tAkltvl, &
                             tvkinv_tAkl,tr_AJ,tr_AV,tr_AJAV,Skl)
!function ME_d_X_over_rij_d computes the following matrix element:
!<\nabla_r \tilde phi_k| (1/r_ij) X | \nabla_r \tilde phi_l>
!Here X is a some real symmetric matrix. If matrix X is not symmetric
!then user needs to symmetrize it before calling this function.
!Input:
!           i,j :: indices of r_ij
!            X  :: n x n real matrix
!           tAk :: n x n real matrix of tAk
!           tAl :: n x n real matrix of tAl
!      inv_tAkl :: n x n real matrix where the inverse of tAk+tAl is stored
!   tvkinv_tAkl :: n-component vector where tvk*inv_tAkl is stored
!   inv_tAkltvl :: n-component vector where inv_tAkl*tvl is stored
!!        tr_AJ :: scalar, tr_AJ = TrAJ(i,j) = tr[inv_tAkl*Jij]
!         tr_AV :: scalar, tr_AV = tau3 = tr[inv_tAkl*tvl*tvk']
!       tr_AJAV :: scalar, tr_AJAV = eta2(i,j) = tr[inv_tAkl*Jij*inv_tAkl*tvl*tvk']
!           Skl :: scalar, overlap Skl=<\tilde phi_k|\tilde phi_l>
    real(wp)   ME_d_X_over_rij_d
    integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles
!Arguments:
    real(wp)   X(nn,nn),tAk(nn,nn),tAl(nn,nn),inv_tAkl(nn,nn)
    integer       i,j,tvk(nn),tvl(nn)
    real(wp)   inv_tAkltvl(nn),tvkinv_tAkl(nn),tr_AJ,tr_AV,tr_AJAV,Skl
!Local variables:
    integer       k,n,s,c
    real(wp)   temp1,temp2
    real(wp)   tAlX(nn,nn),tAlXtAk(nn,nn),Aj(nn),tAkX(nn,nn),tAkXtAl(nn,nn)
    real(wp)   theta,kappa,lambda,omega,chi,theta1,kappa1,chi1,omega1
    real(wp)   h,m,Xtvl(nn),tvkX(nn)

!Doing multiplication Xtvl=X*tvl,tvkX=tvk*X
    n=Glob_n
    do s=1,nn
      temp2=ZERO
      do c=1,nn
        temp2=temp2+X(c,s)*tvl(c)
      enddo
      Xtvl(s)=temp2
    enddo
    do s=1,nn
      temp2=ZERO
      do c=1,nn
        temp2=temp2+tvk(c)*X(c,s)
      enddo
      tvkX(s)=temp2
    enddo
!Doing multiplication tAlX=tAl*X,tAkX=tAk*X
    do s=1,nn
      do c=1,nn
        temp1=ZERO
        temp2=ZERO
        do k=1,nn
          temp1=temp1+tAl(c,k)*X(k,s)
          temp2=temp2+tAk(c,k)*X(k,s)
        enddo
        tAlX(c,s)=temp1
        tAkX(c,s)=temp2
      enddo
    enddo

!Doing multiplication tAlXtAk=tAlX*tAk,tAkXtAl=tAkX*tAl
    do s=1,nn
      do c=1,nn
        temp1=ZERO
        temp2=ZERO
        do k=1,nn
          temp1=temp1+tAlX(c,k)*tAk(k,s)
          temp2=temp2+tAkX(c,k)*tAl(k,s)
        enddo
        tAlXtAk(c,s)=temp1
        tAkXtAl(c,s)=temp2
      enddo
    enddo

!Compute theta=tr[inv_tAkl*(tAl*X*tAk+tAk*X*tAl)]/2
    theta=ZERO
    theta1=ZERO
    do s=1,nn
      temp1=ZERO
      temp2=ZERO
      do k=1,nn
        temp1=temp1+inv_tAkl(s,k)*tAlXtAk(k,s)
        temp2=temp2+inv_tAkl(s,k)*tAkXtAl(k,s)
      enddo
      theta=theta+temp1
      theta1=theta1+temp2
    enddo
    theta=(theta+theta1)/2

!Form Aj=inv_tAkl*ji        j/=i
!     Aj=inv_tAkl*(ji-jj)   j/=i
!Remember that Jii=ji*ji' and Jij=(ji-jj)*(ji-jj)'

    if (i==j) then
      do s=1,nn
        Aj(s)=inv_tAkl(s,i)
      enddo
    else
      do s=1,nn
        Aj(s)=inv_tAkl(s,i)-inv_tAkl(s,j)
      enddo
    endif

!Compute kappa=tr[inv_tAkl*Jij*inv_tAkl*(tAl*X*tAk+tAk*X*tAl)]/2 = Aj' (tAlXtAk +tAlXtAk)Aj/2
    kappa=ZERO
    kappa1=ZERO
    do s=1,nn
      temp1=ZERO
      temp2=ZERO
      do c=1,nn
        temp1=temp1+Aj(c)*tAlXtAk(c,s)*Aj(s)
        temp2=temp2+Aj(c)*tAkXtAl(c,s)*Aj(s)
      enddo
      kappa=kappa+temp1
      kappa1=kappa1+temp2
    enddo
    kappa=(kappa+kappa1)/2

!Compute temp1 = j^{ij}' inv_tAkl tvl = j^{ij} ' inv_tAkltvl
!        temp2 = tvk' inv_tAkl j^{ij} = tvkinv_tAkl' j^{ij}
    if (i==j) then
      temp1=inv_tAkltvl(i)
      temp2=tvkinv_tAkl(i)
    else
      temp1=inv_tAkltvl(i)-inv_tAkltvl(j)
      temp2=tvkinv_tAkl(i)-tvkinv_tAkl(j)
    endif

!Compute the following quantities
!lambda = tr[inv_tAkl*tAl*X*tAk*inv_tAkl*V] = tvkinv_tAkl' tAlXtAk inv_tAkltvl
! omega = tr[inv_tAkl*(tAl*X*tAk+tAk*X*tAl)*inv_tAkl*Jij*inv_tAkl*V]/2 = (tvkinv_tAkl' (tAlXtAk+tAkXtAl) Aj) * temp1/2
!   chi = tr[inv_tAkl*Jij*inv_tAkl*(tAl*X*tAk+tAk*X*tAl)*inv_tAkl*V]/2 = temp2 * (Aj' (tAlXtAk+tAkXtAl) inv_tAkl)/2
!h =  tr[inv_tAkl*Jij*inv_tAkl*tAk*X*tvl*tvk]
!m=   tr[X*tAl*inv_tAkl*Jij*inv_tAkl*tvl*tvk]
    lambda=ZERO
    omega=ZERO
    omega1=ZERO
    chi=ZERO
    chi1=ZERO
    h=ZERO
    m=ZERO
    do s=1,nn
      do c=1,nn
        lambda=lambda+tvkinv_tAkl(c)*tAlXtAk(c,s)*inv_tAkltvl(s)
        omega=omega+tvkinv_tAkl(c)*tAlXtAk(c,s)*Aj(s)
        omega1=omega1+tvkinv_tAkl(c)*tAkXtAl(c,s)*Aj(s)
        m=m+tvkX(c)*tAl(c,s)*Aj(s)
        chi=chi+Aj(c)*tAlXtAk(c,s)*inv_tAkltvl(s)
        chi1=chi1+Aj(c)*tAkXtAl(c,s)*inv_tAkltvl(s)
        h=h+Aj(c)*tAk(c,s)*Xtvl(s)
      enddo
    enddo
    omega=omega*temp1
    omega1=omega1*temp1
    omega=(omega+omega1)/2
    m=m*temp1
    chi=temp2*chi
    chi1=temp2*chi1
    chi=(chi+chi1)/2
    h=temp2*h

    temp1=(4*Skl)/(15*Glob_SqrtPi*tr_AV*tr_AJ*tr_AJ*sqrt(tr_AJ))
    temp2=15*kappa*tr_AJAV + &
           5*tr_AJ*(6*lambda*tr_AJ+9*tr_AV*theta*tr_AJ-3*kappa*tr_AV-2*omega-2*chi+h+m-3*theta*tr_AJAV)
    ME_d_X_over_rij_d=temp1*temp2

  end function ME_d_X_over_rij_d

  function ftransaux(x)
!This function evaluates
!f(x) = ( |x| - sqrt(1-x^2)arccos(sqrt(1-x^2)) ) / (x |x|)
!given a real valued -1<x<1 argument. |x| stands for absolute value.
!A series representation is employed for |x|<xmin
!Depending on the kind parameter (wp=4,8,10,16 - double, extended, or quadruple
!precision) for real numbers, a different xmin is used.
!In all cases the accuracy is close to machine precision corresponding to
!that kind parameter (1-2 last significant figures might be inaccurate in
!the worst case).
    real(wp) ftransaux,x
    real(wp),parameter :: xmin_4=0.30_wp !for single precision
    real(wp),parameter :: xmin_8=0.27_wp !for double precision
    real(wp),parameter :: xmin_10=0.2_wp !for extended precision
    real(wp),parameter :: xmin_16=0.065_wp !for quadruple precision
!Local variables
    real(wp) x2,ax,t,xmin

    ax=abs(x)
    selectcase (wp)
    case(0:4)
      xmin=xmin_4
    case(5:8)
      xmin=xmin_8
    case(9:10)
      xmin=xmin_10
    case(11:16)
      xmin=xmin_16
    endselect
    if (ax<xmin) then
      x2=x*x
      t=(524288.0_wp/50702925.0_wp)
      t=(262144.0_wp/22309287.0_wp)+x2*t
      t=(65536.0_wp/4849845.0_wp)+x2*t
      t=(32768.0_wp/2078505.0_wp)+x2*t
      t=(2048.0_wp/109395.0_wp)+x2*t
      t=(1024.0_wp/45045.0_wp)+x2*t
      t=(256.0_wp/9009.0_wp)+x2*t
      t=(128.0_wp/3465.0_wp)+x2*t
      t=(16.0_wp/315.0_wp)+x2*t
      t=(8.0_wp/105.0_wp)+x2*t
      t=(2.0_wp/15.0_wp)+x2*t
      t=(1.0_wp/3.0_wp)+x2*t
      ftransaux=x*t
    else
      t=sqrt(1.0_wp-x*x)
      ftransaux=(ax-t*acos(t))/(ax*x)
    endif

  end function ftransaux

  subroutine spinPreCalc(n, nFactorial, SziME, parityFactor, SSFmassChargeCoefficient, SSNCmassChargeCoefficient, &
                         SOmassChargeCoefficient, AMMmassChargeCoefficient, &
                         AMMFinmassChargeCoefficient, AnihMassChargeCoefficient, ketMatrix, spatialYoung, &
                         positronPosition, numberOfSpinFunctions, spinFreeME, SiSjME, SSNCspinME)
    use spinStuff
    implicit none

    character(len = maxLen), intent(in) :: spatialYoung
    integer, intent(in) :: n, nFactorial

    real(wp), dimension(nFactorial), intent(out) :: parityFactor
    real(wp), dimension(n, n, 4), intent(out) :: SOmassChargeCoefficient, AMMmassChargeCoefficient, AMMFinmassChargeCoefficient
    real(wp), dimension(n, n, nFactorial), intent(out) :: ketMatrix
    real(wp), dimension(n, n), intent(out) :: SSFmassChargeCoefficient, AnihMassChargeCoefficient, SSNCmassChargeCoefficient

    integer, intent(out) :: positronPosition, numberOfSpinFunctions

    real(wp), dimension(nFactorial), intent(out) :: spinFreeME
    real(wp), dimension(n, 2, nFactorial), intent(out) :: SziME
    real(kind = wp), dimension(n, n, 2, nFactorial), intent(out) :: SiSjME, SSNCspinME

    ! local variables
    integer :: i, j, k, l, m
    character(len = maxLen) :: mySpatialYoung
    integer, dimension(nFactorial) :: parities
    integer, dimension(n, n, nFactorial) :: allPermutations

    SSFmassChargeCoefficient = ZERO
    SSNCmassChargeCoefficient = ZERO
    do i = 1, n
      do j = 1, n
        SSFmassChargeCoefficient(i, j) = -Glob_PseudoCharge(i) * Glob_PseudoCharge(j) / &
                                         (Glob_Mass(i + 1) * Glob_Mass(j + 1)) * EIGHT * Glob_Pi / THREE
        SSNCmassChargeCoefficient(i, j) = Glob_PseudoCharge(i) * Glob_PseudoCharge(j) / &
                                          (Glob_Mass(i + 1) * Glob_Mass(j + 1))

      enddo
    enddo

    SOmassChargeCoefficient = ZERO
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
                                    allPermutations, parities, spinFreeME, SziME, SiSjMe, SSNCspinME)

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

  subroutine spinDependentMatrixElements(m_k, m_l, vechLk, vechLl, Pket, &
                                         SziME, SSNCspinME, SSNCmassChargeCoefficient, SOmassChargeCoefficient, &
                                         AMMmassChargeCoefficient, AMMFinmassChargeCoefficient, SSNCkl, SO1kl, SO2kl, &
                                         AMM1kl, AMM2kl, AMM1Finkl, AMM2Finkl, numberOfSpinFunctions)
    !This subroutine computes symmetry adapted matrix element
    !with two real L=1 correlated Gaussians. These matrix element
    !is used in calculations of expectation values.

    !Input:
    !   m_k, m_l :: integers that determine which z-component is in the
    !       premultiplier of the Gaussian
    !   vechLk, vechLl :: Arrays of length (n(n+1)/2) of exponential parameters.

    !Output (all matrix elements are computed with normalized functions):

    !   SSNCkl :: Non-contact spin-spin term (without the factor of alpha**2)
    !   SO1kl, SO2kl  :: Spin-Orbit corrections (without the factor of alpha**2)
    !         1 and 2 stay for spin-same orbit and spin-another orbit contributions
    !   AMM1kl, AMM2kl  :: AMM corrections (without the factor of alpha**2)
    !         1 and 2 stay for spin-same orbit and spin-another orbit contributions

    !Arguments
    integer,intent(in)       :: m_k, m_l, numberOfSpinFunctions
    real(wp),intent(in)   :: vechLk(Glob_np), vechLl(Glob_np)
    real(wp),intent(in)   :: Pket(Glob_n,Glob_n)

    real(wp), dimension(numberOfSpinFunctions), intent(out)  :: SO1kl, SO2kl, AMM1kl, AMM2kl, &
                                                                   AMM1Finkl, AMM2Finkl, SSNCkl
    !Parameters (These are needed to declare static arrays. Using static
    !arrays makes the function call a little faster in comparison with
    !the case when arrays are dynamically allocated in stack)
    integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles
    integer,parameter :: nnp=nn*(nn+1)/2
    real(wp),intent(in)   :: SSNCspinME(Glob_n, Glob_n, numberOfSpinFunctions), &
                                SziME(Glob_n, numberOfSpinFunctions), &
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
    real(wp)       temp1, temp2, det_tAkl
    integer :: i, j, k, indx

    integer :: pm_k, pm_l ! new non-zero components of v_k and v_l
    real(wp) :: commonFactor, gamma, gamma_diag, jiVl, jiAlAklinvVk, jiAlAklinvVl, jiAklinvVk, jiAklinvVl, &
                   jjAlAklinvVl, jjAlAklinvVk, jjAklinvVk, jjAklinvVl, jjVl, localEps

    integer :: indexI, indexJ ! indices enumerating particles from H_SO and AMM operators

    localEps = 1.0e-14_wp ! if the corresponding spin mean value is less then localEps, we don't calculate the spatial part

    ! basically copy-paste from the old ExpecVals subroutine
    n=Glob_n
    np=Glob_np

    !First we build matrices Lk, Ll, Ak, Al from vechLk, vechLl.
    indx=0
    do i=1,nn
      do j=i,nn
        indx=indx+1
        Lk(i,j)=ZERO
        Lk(j,i)=vechLk(indx)
        Ll(i,j)=ZERO
        Ll(j,i)=vechLl(indx)
      enddo
    enddo

    do i=1,nn
      do j=i,nn
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
    do i=1,nn
      do j=1,nn
        temp1=ZERO
        temp2=ZERO
        do k=1,nn
          temp1=temp1+Pket(k,j)*tAl(k,i)
        enddo
        W1(j,i)=temp1
      enddo
    enddo
    !tAl=W1*Pket

    do i=1,nn
      do j=i,nn
        temp1=ZERO
        temp2=ZERO
        do k=1,nn
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
    do i=1,nn
      do j=i,nn
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
        inv_tAkl(i,j)=temp1
        inv_tAkl(j,i)=temp1
      enddo
    enddo

    ! new code from here

    ! new m_k and m_l
    ! new v_l = (P TRANSPOSED) * v_l
    pm_k = m_k

    pm_l = m_l
    do i = 1, nn
      if (abs(Pket(m_l, i) - 1.0_wp) < 1.0e-13_wp) then ! for integers it would be == 1
        pm_l = i
        exit
      endif
    enddo

    !common factor
    commonFactor = TWO * Glob_PiRaised3n2 / (Glob_SqrtPi * det_tAkl * sqrt(det_tAkl))

    SO1kl = ZERO
    SO2kl = ZERO

    AMM1kl = ZERO
    AMM2kl = ZERO
    AMM1Finkl = ZERO
    AMM2Finkl = ZERO

    do indexI = 1, nn

      i = 0
      do k = 1, numberOfSpinFunctions
        if (abs(SziME(indexI, k)) < localEps) i = i + 1
      enddo
      if (i == numberOfSpinFunctions) cycle ! corresponding <Szi> = 0

      ! gamma diagonal coefficient
      gamma_diag = ONE / sqrt(inv_tAkl(indexI, indexI))

      gamma = gamma_diag ! for spin-same-orbit

      ! kronecker deltas
      jiVl = ZERO
      if (pm_l == indexI) then
        jiVl = ONE
      endif

      ! calculating all the traces we need
      ! tr(Axy') is computed as (y, Ax) everywhere
      ! variable names: jiAlAklinvVk = (j^i, A_l A_{kl}^(-1) v_k) (names doesn't account for permutations)

      jiAlAklinvVl = ZERO
      do i = 1, nn
        jiAlAklinvVl = jiAlAklinvVl + tAl(indexI, i) * inv_tAkl(i, pm_l)
      enddo

      jiAlAklinvVk = ZERO
      do i = 1, nn
        jiAlAklinvVk = jiAlAklinvVk + tAl(indexI, i) * inv_tAkl(i, pm_k)
      enddo

      jiAklinvVl = inv_tAkl(indexI, pm_l)

      jiAklinvVk = inv_tAkl(indexI, pm_k)

      ! diagonal (spin-same orbit) matrix element
      temp1 =  &
        gamma**3 / THREE * (jiVl * jiAklinvVk + &
                            jiAlAklinvVk * jiAklinvVl - &
                            jiAklinvVk * jiAlAklinvVl)

      do k = 1, numberOfSpinFunctions
        SO1kl(k) = SO1kl(k) + SziME(indexI, k) * SOmassChargeCoefficient(indexI, indexI, 1) * temp1
        AMM1kl(k) = AMM1kl(k) + SziME(indexI, k) * AMMmassChargeCoefficient(indexI, indexI, 1) * temp1
        AMM1Finkl(k) = AMM1Finkl(k) + SziME(indexI, k) * AMMFinmassChargeCoefficient(indexI, indexI, 1) * temp1
      enddo

      ! these traces are needed for spin-other orbit contribution

      do indexJ = 1, nn
        if (indexI == indexJ) cycle

        gamma = ONE / sqrt(inv_tAkl(indexI, indexI) + inv_tAkl(indexJ, indexJ) - &
                           inv_tAkl(indexI, indexJ) - inv_tAkl(indexJ, indexI))

        ! we need more traces

        jjAlAklinvVl = ZERO
        do i = 1, nn
          jjAlAklinvVl = jjAlAklinvVl + tAl(indexJ, i) * inv_tAkl(i, pm_l)
        enddo

        jjAlAklinvVk = ZERO
        do i = 1, nn
          jjAlAklinvVk = jjAlAklinvVk + tAl(indexJ, i) * inv_tAkl(i, pm_k)
        enddo

        jjAklinvVk = inv_tAkl(indexJ, pm_k)

        jjAklinvVl = inv_tAkl(indexJ, pm_l)

        ! kronecker deltas
        jjVl = ZERO
        if (pm_l == indexJ) then
          jjVl = ONE
        endif

        ! NOTE for this term we need gamma_diag (tr[Cklinv J_ii]), not gamma (tr[Cklinv J_ij])
        temp1 = &
          gamma_diag**3 / THREE *  (jiAklinvVk * jjVl + &
                                    jjAlAklinvVk * jiAklinvVl - &
                                    jiAklinvVk * jjAlAklinvVl)

        do k = 1, numberOfSpinFunctions
          SO2kl(k) = SO2kl(k) + SziME(indexI, k) * SOmassChargeCoefficient(indexI, indexI, 2) * temp1
          AMM2Finkl(k) = AMM2Finkl(k) + SziME(indexI, k) * AMMFinmassChargeCoefficient(indexI, indexI, 2) * temp1
        enddo

        temp1 = &
          gamma**3 / THREE * (jjVl * (jjAklinvVk - jiAklinvVk) + &
                              jjAlAklinvVk * (jjAklinvVl - jiAklinvVl) + &
                              jjAlAklinvVl * (jiAklinvVk - jjAklinvVk))

        do k = 1, numberOfSpinFunctions
          SO2kl(k) = SO2kl(k) + SziME(indexI, k) * SOmassChargeCoefficient(indexI, indexJ, 3) * temp1
          AMM2kl(k) = AMM2kl(k) + SziME(indexI, k) * AMMmassChargeCoefficient(indexI, indexJ, 3) * temp1
          AMM2Finkl(k) = AMM2Finkl(k) + SziME(indexI, k) * AMMFinmassChargeCoefficient(indexI, indexJ, 3) * temp1
        enddo

        temp1 = &
          gamma**3 / THREE * (jiVl * (jiAklinvVk - jjAklinvVk) + &
                              jiAlAklinvVk * (jiAklinvVl - jjAklinvVl) + &
                              jiAlAklinvVl * (jjAklinvVk - jiAklinvVk))

        do k = 1, numberOfSpinFunctions
          SO2kl(k) = SO2kl(k) + SziME(indexI, k) * SOmassChargeCoefficient(indexI, indexJ, 4) * temp1
          AMM2kl(k) = AMM2kl(k) + SziME(indexI, k) * AMMmassChargeCoefficient(indexI, indexJ, 4) * temp1
          AMM2Finkl(k) = AMM2Finkl(k) + SziME(indexI, k) * AMMFinmassChargeCoefficient(indexI, indexJ, 4) * temp1
        enddo

      enddo ! indexJ cycle

    enddo ! indexI cycle

    ! SS term (separate loop for not to interfere with the case Siz == 0 at previous loop)
    SSNCkl = ZERO
    do indexI = 1,nn
      do indexJ = indexI + 1,nn !(we need only indexJ > indexI)

        gamma = ONE / sqrt(inv_tAkl(indexI, indexI) + inv_tAkl(indexJ, indexJ) - &
                           inv_tAkl(indexI, indexJ) - inv_tAkl(indexJ, indexI))

        jjAklinvVk = inv_tAkl(indexJ, pm_k)
        jjAklinvVl = inv_tAkl(indexJ, pm_l)
        jiAklinvVl = inv_tAkl(indexI, pm_l)
        jiAklinvVk = inv_tAkl(indexI, pm_k)

        temp1 = &
          (gamma**5 / 15._wp) * ( jiAklinvVk * (jiAklinvVl  - jjAklinvVl) + &
                                    jjAklinvVk * (jjAklinvVl - jiAklinvVl) ) !additional factor of 1/sqrt(6) is taken from spin part

        do k = 1, NumberOfSpinFunctions
          if (abs(SSNCspinME(indexI, indexJ, k)) < localEps) cycle
          SSNCkl(k) = SSNCkl(k) + SSNCspinME(indexI, indexJ, k) * SSNCmassChargeCoefficient(indexI, indexJ) * temp1
        enddo

      enddo !indexJ loop
    enddo !indexI loop

    SSNCkl = SSNCkl * commonFactor
    SO1kl = SO1kl * commonFactor
    SO2kl = SO2kl * commonFactor
    AMM1kl = AMM1kl * commonFactor
    AMM2kl = AMM2kl * commonFactor
    AMM1Finkl = AMM1Finkl * commonFactor
    AMM2Finkl = AMM2Finkl * commonFactor

  end subroutine spinDependentMatrixElements

  subroutine OverlapMatrixElement_RG_1P(m_k, vechLk, P, Skk)
!This subroutine computes symmetry adapted matrix element with
!two real L=1 correlated Gaussians:
!
!fk = z_{m_k} exp[-r'(Lk*Lk')r]
!
!where m_k is some integer between 1 and n (n is the number of
!pseudoparticles). Symmetry adaption is applied to the ket using
!permutation matrix P
!
!Input:
!   m_k :: integer that determine which z-component is in the
!                premultiplier of the Gaussian
!   vechLk :: Array of length (n(n+1)/2) of
!     exponential parameters.
!   P  :: The symmetry permutation matrix of size n x n
!Output:
!   Skk         ::        Overlap matrix element (normalized)

!Arguments
    integer,intent(in)          :: m_k
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
    integer           tvl(nn)
    real(wp)       Lk(nn,nn),Ll(nn,nn)
    real(wp)       Ak(nn,nn),tAl(nn,nn),tAkl(nn,nn)
    real(wp)       inv_tAkl(nn,nn)
    real(wp)       W1(nn,nn)
    real(wp)       inv_tAkltvl(nn),vkinv_tAkl(nn)
    real(wp)       temp1
    real(wp)       det_tAkl
    real(wp)       tau3
    integer           i,j,k,q,t,indx

    n=Glob_n
    np=Glob_np
!First we build matrices Lk, Ll, Ak, Al from vechLk, vechLl.
    indx=0
    do i=1,nn
      do j=i,nn
        indx=indx+1
        Lk(i,j)=ZERO
        Lk(j,i)=vechLk(indx)
        Ll(i,j)=Lk(i,j)
        Ll(j,i)=Lk(j,i)
      enddo
    enddo

    do i=1,nn
      do j=i,nn
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
    do i=1,nn
      do j=1,nn
        temp1=ZERO
        do k=1,nn
          temp1=temp1+P(k,j)*tAl(k,i)
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
    do i=1,nn
      do j=i,nn
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
        inv_tAkl(i,j)=temp1
        inv_tAkl(j,i)=temp1
      enddo
    enddo

!Computing tvl=P'*vl
    do i=1,nn
      tvl(i)=P(m_k,i)
    enddo

!Compute inv_tAkltvl = inv_tAkl * tvl
    do i=1,nn
      temp1=ZERO
      do j=1,nn
        temp1=temp1+inv_tAkl(j,i)*tvl(j)
      enddo
      inv_tAkltvl(i)=temp1
    enddo

!Compute vkinv_tAkl=vk'*inv_tAkl
    do i=1,nn
      vkinv_tAkl(i)=inv_tAkl(m_k,i)
    enddo

!Compute tau3=vkinv_tAkl*tvl
    tau3=ZERO
    do i=1,nn
      tau3=tau3+vkinv_tAkl(i)*tvl(i)
    enddo

!Evaluating overlap
!temp1=abs(det_Ll*det_Lk)/det_tAkl
    temp1=det_tAkl*sqrt(det_tAkl)
!Skl=Glob_2Raised3n2*tau3*temp1*sqrt(temp1/(inv_Akk(m_k,m_k)*inv_All(m_l,m_l)))
    Skk=Glob_PiRaised3n2*tau3/(TWO*temp1)

  end subroutine OverlapMatrixElement_RG_1P

end module matelem
