module matelem
!Module matelem contains subroutines for computing
!matrix elements with complex L=0 Gaussians.
  use globvars
  implicit none

contains

  subroutine MatrixElementsHS_CG_0S(Paramk, Paraml, P, Hkl, Skl, Dk, Dl, grad_k, grad_l)
!This subroutine computes symmetry adapted matrix element
!with two complex L=0 correlated Gaussians
!
!fk = exp[-r'(Lk*Lk'+iBk)r]
!
!symmetry adaption is applied to the ket using a permutation matrix P
!(these matrices P are stored in Glob_YHYMatr(:,:,1:Glob_NumYHYTerms))
!
!Input:
!   Paramk,Paraml :: Arrays of length 2*(n(n+1)/2) of
!     exponential parameters. They are ordered in such a way
!     that Paramk=(vechLk,vechBk).
!   P  :: The symmetry permutation matrix of size n x n
!   grad_k, grad_l :: Gradient flags
!   grad_k=.true.  means that dHkldvechLk, dHkldvechBk, dSkldvechLk, dSkldvechBk need to be computed.
!   grad_l=.true.  means that dHkldvechLl, dHkldvechBl, dSkldvechLl, dSkldvechBl need to be computed.
!Output:
!   Hkl         ::        Hamiltonian term (normalized)
!   Skl         ::        Overlap matrix element (normalized)
!   Dk,Dl:: derivatives of normalized Hkl and Skl wrt Paramk
!           and Paraml respectively. They are ordered in the
!           following manner:
!           Dk=(dHkldvechLk,dHkldvechBk,dSkldvechLk,dSkldvechBk)
!           Dl=(dHkldvechLl,dHkldvechBl,dSkldvechLl,dSkldvechBl)

!Arguments
    real(wp),intent(in)      :: Paramk(2*Glob_np), Paraml(2*Glob_np)
    integer,intent(in)          :: P(Glob_n,Glob_n)
    complex(wp),intent(out)  :: Skl,Hkl
    complex(wp),intent(out)  :: Dk(4*Glob_np),Dl(4*Glob_np)
    logical,intent(in)          :: grad_k, grad_l

!Parameters (These are needed to declare static arrays. Using static
!arrays makes the function call a little faster in comparison with
!the case when arrays are dynamically allocated in stack)
    integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles
    integer,parameter :: nnp=nn*(nn+1)/2

!Local variables
    integer           n, np
    real(wp)       vechLk(nnp), vechBk(nnp)
    real(wp)       vechLl(nnp), vechBl(nnp)
    complex(wp)    dHkldvechLk(nnp), dHkldvechBk(nnp)
    complex(wp)    dHkldvechLl(nnp), dHkldvechBl(nnp)
    complex(wp)    dSkldvechLk(nnp), dSkldvechBk(nnp)
    complex(wp)    dSkldvechLl(nnp), dSkldvechBl(nnp)
    real(wp)       Lk(nn,nn), Ll(nn,nn), inv_Lk(nn,nn), inv_Ll(nn,nn)
    real(wp)       Ak(nn,nn), tAl(nn,nn), Bk(nn,nn), tBl(nn,nn)
    complex(wp)    Ck(nn,nn), tCl(nn,nn), tCkl(nn,nn), ttCkl(nn,nn)
    complex(wp)    inv_tCkl(nn,nn), inv_ttCkl(nn,nn)
    complex(wp)    inv_tCkltClM(nn,nn)
    complex(wp)    tr_inv_tCklJij32(nn,nn)
    complex(wp)    F(nn,nn), G(nn,nn)
    complex(wp)    tQ(nn,nn), ttQ(nn,nn)
    complex(wp)    Wc1(nn,nn), Wc2(nn,nn), Wc3(nn,nn)
    complex(wp)    WVcLk(nnp), WVcLl(nnp), WVcBk(nnp), WVcBl(nnp)
    real(wp)       Wr1(nn,nn), Wr2(nn,nn), Wr3(nn,nn), Wr4(nn,nn)
    real(wp)       rtemp, rsum, rnumb
    complex(wp)    ctemp, csum, cnumb, cmpnum,cRklij
    real(wp)       det_Lk, det_Ll
    complex(wp)    det_tCkl
    complex(wp)    Tkl, Vkl
    integer           i, j, k, kk, kkk, q, t, indx
    logical           grad_k, grad_l, grad_kl

    n=Glob_n
    np=Glob_np

!extracting nonlinear parameters
    vechLk(1:np)=Paramk(1:np)
    vechBk(1:np)=Paramk(np+1:np+np)
    vechLl(1:np)=Paraml(1:np)
    vechBl(1:np)=Paraml(np+1:np+np)

!First we build matrices Lk, Ll, Ak, Al, Bk, Bl
!from vechLk, vechLl, vechBk, vechBl.
    indx=0
    do i=1,n
      do j=i,n
        indx=indx+1
        Lk(i,j)=ZERO
        Lk(j,i)=vechLk(indx)
        Ll(i,j)=ZERO
        Ll(j,i)=vechLl(indx)
        Bk(j,i)=vechBk(indx)
        Bk(i,j)=vechBk(indx)
        tBl(j,i)=vechBl(indx)
        tBl(i,j)=vechBl(indx)
      enddo
    enddo

    do i=1,n
      do j=i,n
        rtemp=ZERO
        do k=1,j
          rtemp=rtemp+Lk(i,k)*Lk(j,k)
        enddo
        Ak(i,j)=rtemp
        Ak(j,i)=rtemp
        rtemp=ZERO
        do k=1,j
          rtemp=rtemp+Ll(i,k)*Ll(j,k)
        enddo
        tAl(i,j)=rtemp
        tAl(j,i)=rtemp
      enddo
    enddo

!Then we permute elements of Al and Bl to account for
!the action of the permutation matrix
!  tBl=P'*Bl*P  tAl=P'*Al*P
    do i=1,n
      do j=1,n
        rnumb=ZERO
        rtemp=ZERO
        do k=1,n
          rnumb=rnumb+P(k,i)*tAl(k,j)
          rtemp=rtemp+P(k,i)*tBl(k,j)
        enddo
        Wr1(i,j)=rnumb
        Wr2(i,j)=rtemp
      enddo
    enddo
    do i=1,n
      do j=i,n
        rnumb=ZERO
        rtemp=ZERO
        do k=1,n
          rnumb=rnumb+Wr1(i,k)*P(k,j)
          rtemp=rtemp+Wr2(i,k)*P(k,j)
        enddo
        tAl(i,j)=rnumb
        tAl(j,i)=rnumb
        tBl(i,j)=rtemp
        tBl(j,i)=rtemp
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

!Next, we form Ck, tCl, and tCkl=tCl+Ck
!Notice that array Ck will actually contain Ck^*
    do i=1,n
      do j=i,n
        cnumb=cmplx(Ak(j,i),-Bk(j,i),wp)
        Ck(j,i)=cnumb
        Ck(i,j)=cnumb
        cnumb=cmplx(tAl(j,i),tBl(j,i),wp)
        tCl(j,i)=cnumb
        tCl(i,j)=cnumb
        cnumb=Ck(j,i)+tCl(j,i)
        tCkl(j,i)=cnumb
        tCkl(i,j)=cnumb
      enddo
    enddo

!After this we can do pseudo-Cholesky factorization of tCkl
!(pseudo because it is a product of a lower triangular
!matrix and its transpose, not a product of a lower triangular
!matrix and its hermitian conjugate as the original Cholesky
!factorization)
!Matrix tCkl will be destroyed. Instead we will have
!its pseudo-Cholesky factor stored in the lower triangle
    det_tCkl=COMPLEXONE
    do i=1,n
      do j=i,n
        csum=tCkl(i,j)
        do k=i-1,1,-1
          csum=csum-tCkl(i,k)*tCkl(j,k)
        enddo
        if (i==j) then
          tCkl(i,i)=sqrt(csum)
          det_tCkl=det_tCkl*csum
        else
          tCkl(j,i)=csum/tCkl(i,i)
        endif
      enddo
    enddo

!Inverting tCkl using pseudo-Cholesky factor
!and placing the result into inv_tCkl
    do i=1,n
      tCkl(i,i)=ONE/tCkl(i,i)
      do j=i+1,n
        csum=COMPLEXZERO
        do k=i,j-1
          csum=csum-tCkl(j,k)*tCkl(k,i)
        enddo
        tCkl(j,i)=csum/tCkl(j,j)
      enddo
    enddo

    do i=1,n
      do j=i,n
        ctemp=COMPLEXZERO
        do k=j,n
          ctemp=ctemp+tCkl(k,i)*tCkl(k,j)
        enddo
        inv_tCkl(i,j)=ctemp
        inv_tCkl(j,i)=ctemp
      enddo
    enddo

!Evaluating overlap
    cnumb=abs(det_Ll*det_Lk)/det_tCkl
    Skl=Glob_2Raised3n2*cnumb*sqrt(cnumb)

!Doing multiplication Cl*invCkl*Ck^*
    do i=1,n
      do j=1,n
        ctemp=COMPLEXZERO
        do k=1,n
          ctemp=ctemp+inv_tCkl(i,k)*Ck(k,j)
        enddo
        Wc1(i,j)=ctemp
      enddo
    enddo
    do i=1,n
      do j=1,n
        ctemp=COMPLEXZERO
        do k=1,n
          ctemp=ctemp+tCl(i,k)*Wc1(k,j)
        enddo
        Wc2(i,j)=ctemp
      enddo
    enddo

!Evaluating kinetic energy, Tkl
    ctemp=COMPLEXZERO
    do i=1,n
      do k=1,n
        ctemp=ctemp+Glob_MassMatrix(i,k)*Wc2(k,i)
      enddo
    enddo
    Tkl=6*Skl*ctemp

!Evaluating potential energy, Vkl, and tr[invCkl*Jij]^(-3/2)
!The lower triangle of array trinvCklJij32
!will contain the corresponding quantities.
    cnumb=(TWO/Glob_SqrtPi)*Skl
    Vkl=COMPLEXZERO
    do i=1,n
      cmpnum=inv_tCkl(i,i)
      csum=sqrt(cmpnum)
      tr_inv_tCklJij32(i,i)=1/(csum*cmpnum)
      cRklij=cnumb/csum
      Vkl=Vkl+Glob_ScaledPseudoChargeMatrix(i,0)*cRklij
    enddo
    do i=1,n
      do j=i+1,n
        cmpnum=inv_tCkl(i,i)+inv_tCkl(j,j)-inv_tCkl(j,i)-inv_tCkl(j,i)
        csum=sqrt(cmpnum)
        tr_inv_tCklJij32(j,i)=1/(csum*cmpnum)
        cRklij=cnumb/csum
        Vkl=Vkl+Glob_ScaledPseudoChargeMatrix(i,j)*cRklij
      enddo
    enddo

    Hkl=Tkl+Vkl

!Now we start evaluating gradient
!Gradient of Skl

    if (grad_k) then
      !Multiplying inv_tCkl and Lk and storing it in Wc1
      do i=1,n
        do j=1,n
          csum=COMPLEXZERO
          do k=j,n
            csum=csum+inv_tCkl(i,k)*Lk(k,j)
          enddo
          Wc1(i,j)=csum
        enddo
      enddo
      !Inverting Lk. The inverse of a lower triangular matrix is a lower
      !triangular matrix. The inverse of Lk will be stored in inv_Lk. Upper
      !triangle of inv_Lk is set to be zero.
      do i=1,n
        do j=i,n
          inv_Lk(j,i)=Lk(j,i)
        enddo
      enddo
      do i=1,n
        inv_Lk(i,i)=1/inv_Lk(i,i)
        do j=i+1,n
          inv_Lk(i,j)=ZERO
          rsum=ZERO
          do k=i,j-1
            rsum=rsum-inv_Lk(j,k)*inv_Lk(k,i)
          enddo
          inv_Lk(j,i)=rsum/inv_Lk(j,j)
        enddo
      enddo
      !storing dSkldvechLk
      ctemp=Skl*THREEHALF
      indx=0
      do i=1,n
        do j=i,n
          indx=indx+1
          dSkldvechLk(indx)=(inv_Lk(i,j)-TWO*Wc1(j,i))*ctemp
        enddo
      enddo
      !storing dSkldvechBk
      ctemp=Skl*cmplx(ZERO,THREEHALF,wp)
      cnumb=TWO*ctemp
      indx=0
      do i=1,n
        indx=indx+1
        dSkldvechBk(indx)=inv_tCkl(i,i)*ctemp
        do j=i+1,n
          indx=indx+1
          dSkldvechBk(indx)=inv_tCkl(j,i)*cnumb
        enddo
      enddo
    endif !end if (grad_k)

    if (grad_l) then
      !calculating inv_ttCkl: inv_ttCkl=P*inv_tCkl*P'
      do i=1,n
        do j=1,n
          csum=COMPLEXZERO
          do k=1,n
            csum=csum+P(i,k)*inv_tCkl(k,j)
          enddo
          Wc1(i,j)=csum
        enddo
      enddo
      do i=1,n
        do j=i,n
          csum=COMPLEXZERO
          do k=1,n
            csum=csum+Wc1(i,k)*P(j,k)
          enddo
          inv_ttCkl(i,j)=csum
          inv_ttCkl(j,i)=csum
        enddo
      enddo
      !Multiplying inv_ttCkl and Ll and storing it Wc1
      do i=1,n
        do j=1,n
          csum=COMPLEXZERO
          do k=1,n
            csum=csum+inv_ttCkl(i,k)*Ll(k,j)
          enddo
          Wc1(i,j)=csum
        enddo
      enddo
      !Inverting Ll. The inverse of a lower triangular matrix
      !is a lower triangular matrix. The inverse of Ll will be
      !stored in inv_Ll. Upper triangle of inv_Ll is not set to
      !be zero.
      do i=1,n
        do j=i,n
          inv_Ll(j,i)=Ll(j,i)
        enddo
      enddo
      do i=1,n
        inv_Ll(i,i)=1/inv_Ll(i,i)
        do j=i+1,n
          inv_Ll(i,j)=ZERO
          rsum=ZERO
          do k=i,j-1
            rsum=rsum-inv_Ll(j,k)*inv_Ll(k,i)
          enddo
          inv_Ll(j,i)=rsum/inv_Ll(j,j)
        enddo
      enddo
      !storing dSkldvechLl
      ctemp=Skl*THREEHALF
      indx=0
      do i=1,n
        do j=i,n
          indx=indx+1
          dSkldvechLl(indx)=(inv_Ll(i,j)-TWO*Wc1(j,i))*ctemp
        enddo
      enddo
      !storing dSkldvechBl
      ctemp=Skl*cmplx(ZERO,-THREEHALF,wp)
      cnumb=TWO*ctemp
      indx=0
      do i=1,n
        indx=indx+1
        dSkldvechBl(indx)=inv_ttCkl(i,i)*ctemp
        do j=i+1,n
          indx=indx+1
          dSkldvechBl(indx)=inv_ttCkl(j,i)*cnumb
        enddo
      enddo
    endif !end if (grad_l)

!Gradient of Tkl

    if (grad_l) then
      !Compute matrix G=P*inv_tCkl*Ck*M*Ck*inv_tCkl*P'
      do i=1,n
        do j=1,n
          csum=COMPLEXZERO
          do k=1,n
            csum=csum+Ck(i,k)*Glob_MassMatrix(k,j)
          enddo
          Wc1(i,j)=csum
        enddo
      enddo
      do i=1,n
        do j=1,n
          csum=COMPLEXZERO
          do k=1,n
            csum=csum+Wc1(i,k)*Ck(k,j)
          enddo
          Wc2(i,j)=csum
        enddo
      enddo
      do i=1,n
        do j=1,n
          csum=COMPLEXZERO
          do k=1,n
            csum=csum+inv_tCkl(i,k)*Wc2(k,j)
          enddo
          Wc1(i,j)=csum
        enddo
      enddo
      do i=1,n
        do j=i,n
          csum=COMPLEXZERO
          do k=1,n
            csum=csum+Wc1(i,k)*inv_tCkl(k,j)
          enddo
          Wc2(i,j)=csum
          Wc2(j,i)=csum
        enddo
      enddo
      do i=1,n
        do j=1,n
          csum=COMPLEXZERO
          do k=1,n
            csum=csum+P(i,k)*Wc2(k,j)
          enddo
          Wc1(i,j)=csum
        enddo
      enddo
      do i=1,n
        do j=i,n
          csum=COMPLEXZERO
          do k=1,n
            csum=csum+Wc1(i,k)*P(j,k)
          enddo
          G(i,j)=csum
          G(j,i)=csum
        enddo
      enddo
      !Finding the product G*Ll and
      !putting it in array Wc1.
      do i=1,n
        do j=1,n
          csum=COMPLEXZERO
          do k=j,n
            csum=csum+G(i,k)*Ll(k,j)
          enddo
          Wc1(i,j)=csum
        enddo
      enddo
      !storing dTkldvechLl
      ctemp=Tkl/Skl
      csum=12*Skl
      indx=0
      do i=1,n
        do j=i,n
          indx=indx+1
          dHkldvechLl(indx)=ctemp*dSkldvechLl(indx)+csum*Wc1(j,i)
        enddo
      enddo
      !storing dTkldvechBl
      csum=Skl*cmplx(ZERO,SIX,wp)
      indx=0
      do i=1,n
        indx=indx+1
        dHkldvechBl(indx)=ctemp*dSkldvechBl(indx)+csum*G(i,i)
        do j=i+1,n
          indx=indx+1
          dHkldvechBl(indx)=ctemp*dSkldvechBl(indx)+TWO*csum*G(j,i)
        enddo
      enddo
    endif !end if (grad_l)

    if (grad_k) then
      !Compute matrix F=inv_tCkl*Cl*M*Cl*inv_tCkl
      do i=1,n
        do j=1,n
          csum=COMPLEXZERO
          do k=1,n
            csum=csum+tCl(i,k)*Glob_MassMatrix(k,j)
          enddo
          Wc1(i,j)=csum
        enddo
      enddo
      do i=1,n
        do j=1,n
          csum=COMPLEXZERO
          do k=1,n
            csum=csum+Wc1(i,k)*tCl(k,j)
          enddo
          Wc2(i,j)=csum
        enddo
      enddo
      do i=1,n
        do j=1,n
          csum=COMPLEXZERO
          do k=1,n
            csum=csum+inv_tCkl(i,k)*Wc2(k,j)
          enddo
          Wc1(i,j)=csum
        enddo
      enddo
      do i=1,n
        do j=i,n
          csum=COMPLEXZERO
          do k=1,n
            csum=csum+Wc1(i,k)*inv_tCkl(k,j)
          enddo
          F(i,j)=csum
          F(j,i)=csum
        enddo
      enddo
      !Finding the product F*Lk and
      !putting it in array Wc1.
      do i=1,n
        do j=1,n
          csum=COMPLEXZERO
          do k=j,n
            csum=csum+F(i,k)*Lk(k,j)
          enddo
          Wc1(i,j)=csum
        enddo
      enddo
      !storing dTkldvechLk
      ctemp=Tkl/Skl
      csum=12*Skl
      indx=0
      do i=1,n
        do j=i,n
          indx=indx+1
          dHkldvechLk(indx)=ctemp*dSkldvechLk(indx)+csum*Wc1(j,i)
        enddo
      enddo
      !storing dTkldvechBl
      csum=Skl*cmplx(ZERO,-SIX,wp)
      indx=0
      do i=1,n
        indx=indx+1
        dHkldvechBk(indx)=ctemp*dSkldvechBk(indx)+csum*F(i,i)
        do j=i+1,n
          indx=indx+1
          dHkldvechBk(indx)=ctemp*dSkldvechBk(indx)+TWO*csum*F(j,i)
        enddo
      enddo
    endif !end if (grad_k)

!Gradient of Vkl

    if (grad_l.OR.grad_k) then
      !First we set to zero some work arrays
      if (grad_k) then
        WVcLk(1:np)=COMPLEXZERO
        WVcBk(1:np)=COMPLEXZERO
      endif
      if (grad_l) then
        WVcLl(1:np)=COMPLEXZERO
        WVcBl(1:np)=COMPLEXZERO
      endif
      !Now we can proceed
      do i=1,n
        do j=1,i
          !Evaluating tQ=inv_tCkl*Jij*inv_tCkl
          if (i==j) then
            do k=1,n
              tQ(k,k)=inv_tCkl(k,i)*inv_tCkl(i,k)
              do kk=k+1,n
                cnumb=inv_tCkl(k,i)*inv_tCkl(i,kk)
                tQ(k,kk)=cnumb
                tQ(kk,k)=cnumb
              enddo
            enddo
          else
            do k=1,n
              tQ(k,k)=inv_tCkl(k,i)*inv_tCkl(i,k)+inv_tCkl(k,j)*inv_tCkl(j,k) &
                       -inv_tCkl(k,i)*inv_tCkl(j,k)-inv_tCkl(k,j)*inv_tCkl(i,k)
              do kk=k+1,n
                cnumb=inv_tCkl(k,i)*inv_tCkl(i,kk)+inv_tCkl(k,j)*inv_tCkl(j,kk) &
                       -inv_tCkl(k,i)*inv_tCkl(j,kk)-inv_tCkl(k,j)*inv_tCkl(i,kk)
                tQ(k,kk)=cnumb
                tQ(kk,k)=cnumb
              enddo
            enddo
          endif
          if (grad_k) then
            !Multiplying tQ by Lk
            !and storing the results in the lower
            !triangle of Wc1.
            do k=1,n
              do kk=1,k
                csum=COMPLEXZERO
                do kkk=kk,n
                  csum=csum+tQ(k,kkk)*Lk(kkk,kk)
                enddo
                Wc1(k,kk)=csum
              enddo
            enddo
          endif
          if (grad_l) then
            !Evaluating ttQ=P*tQ*P'
            do k=1,n
              do kk=1,n
                csum=COMPLEXZERO
                do kkk=1,n
                  csum=csum+P(k,kkk)*tQ(kkk,kk)
                enddo
                Wc2(k,kk)=csum
              enddo
            enddo
            do k=1,n
              do kk=k,n
                csum=COMPLEXZERO
                do kkk=1,n
                  csum=csum+Wc2(k,kkk)*P(kk,kkk)
                enddo
                ttQ(k,kk)=csum
                ttQ(kk,k)=csum
              enddo
            enddo
            !Multiplying ttQ by Ll
            !and storing the results in the lower
            !triangle of Wc2.
            do k=1,n
              do kk=1,k
                csum=COMPLEXZERO
                do kkk=kk,n
                  csum=csum+ttQ(k,kkk)*Ll(kkk,kk)
                enddo
                Wc2(k,kk)=csum
              enddo
            enddo
          endif

          !Calculating ij-terms of the sums in the expressions for
          !the dVkldvechXk
          if (i==j) then
            cnumb=Glob_ScaledPseudoChargeMatrix(0,i)*tr_inv_tCklJij32(i,i)
          else
            cnumb=Glob_ScaledPseudoChargeMatrix(i,j)*tr_inv_tCklJij32(i,j)
          endif
          if (grad_k) then
            indx=0
            do k=1,n
              indx=indx+1
              WVcBk(indx)=WVcBk(indx)+cnumb*tQ(k,k)
              WVcLk(indx)=WVcLk(indx)+cnumb*Wc1(k,k)
              do kk=k+1,n
                indx=indx+1
                WVcBk(indx)=WVcBk(indx)+2*cnumb*tQ(kk,k)
                WVcLk(indx)=WVcLk(indx)+cnumb*Wc1(kk,k)
              enddo
            enddo
          endif
          if (grad_l) then
            indx=0
            do k=1,n
              indx=indx+1
              WVcBl(indx)=WVcBl(indx)+cnumb*ttQ(k,k)
              WVcLl(indx)=WVcLl(indx)+cnumb*Wc2(k,k)
              do kk=k+1,n
                indx=indx+1
                WVcBl(indx)=WVcBl(indx)+2*cnumb*ttQ(kk,k)
                WVcLl(indx)=WVcLl(indx)+cnumb*Wc2(kk,k)
              enddo
            enddo
          endif
        enddo
      enddo
      !Multiplying by common numbers and getting the final
      !result for the gradient of Vkl
      cnumb=(TWO/Glob_SqrtPi)*Skl
      ctemp=(IMAGINARYUNIT/Glob_SqrtPi)*Skl
      if (grad_k) then
        dHkldvechLk(1:np)=dHkldvechLk(1:np)+(Vkl/Skl)*dSkldvechLk(1:np)+cnumb*WVcLk(1:np)
        dHkldvechBk(1:np)=dHkldvechBk(1:np)+(Vkl/Skl)*dSkldvechBk(1:np)-ctemp*WVcBk(1:np)
      endif
      if (grad_l) then
        dHkldvechLl(1:np)=dHkldvechLl(1:np)+(Vkl/Skl)*dSkldvechLl(1:np)+cnumb*WVcLl(1:np)
        dHkldvechBl(1:np)=dHkldvechBl(1:np)+(Vkl/Skl)*dSkldvechBl(1:np)+ctemp*WVcBl(1:np)
      endif
    endif

!Packing derivatives into the output array
    if (grad_k) then
      Dk(1:np)=dHkldvechLk(1:np)
      Dk(np+1:2*np)=dHkldvechBk(1:np)
      Dk(2*np+1:3*np)=dSkldvechLk(1:np)
      Dk(3*np+1:4*np)=dSkldvechBk(1:np)
    endif
    if (grad_l) then
      Dl(1:np)=dHkldvechLl(1:np)
      Dl(np+1:2*np)=dHkldvechBl(1:np)
      Dl(2*np+1:3*np)=dSkldvechLl(1:np)
      Dl(3*np+1:4*np)=dSkldvechBl(1:np)
    endif

  end subroutine MatrixElementsHS_CG_0S

  subroutine MatrixElementsAll_CG_0S(Paramk, Paraml, Pbra, Pket, Hkl, Skl, Tkl, Vkl, &
                                       rm2kl, rmkl, rkl, r2kl, deltarkl, drach_deltarkl, MVkl, drach_MVkl, Darwinkl, &
                                       drach_Darwinkl, OOkl, NumCFGridPoints, CFGrid, CFkl, NumDensGridPoints, &
                                       DensGrid, Denskl, AreCorrFuncNeeded, ArePartDensNeeded)
!This subroutine computes symmetry adapted matrix elements
!with two complex correlated Gaussians. These matrix elements
!are used in calculations of expectation values.
!Symmetry adaption is applied to the bra and ket using permutation matrices Pbra and Pket
!
!Input:
!   Paramk,Paraml :: Arrays of length 2*(n(n+1)/2) of exponential
!           parameters. They are ordered in such a way
!           that Paramk=(vechLk,vechBk).
!   Pbra :: The symmetry permutation matrix of size n x n that is applied to bra
!   Pket :: The symmetry permutation matrix of size n x n that is applied to ket
!Output (all matrix elements are computed with normalized functions):
!   Hkl         ::        Hamiltonian
!   Skl         ::        Overlap
!   Tkl  :: Kinetic energy
!   Vkl  :: Potential energy
!   rm2kl :: 1/r_i^2, 1/r_{ij}^2
!   rmkl :: 1/r_i, 1/r_{ij}
!   rkl  :: r_i, r_{ij}
!   r2kl :: r_i^2, r_{ij}^2
!deltarkl:: delta(r_i), delta(r_{ij})
!drach_deltarkl:: Drachmanized delta(r_i), delta(r_{ij})
!   MVkl :: Mass-velocity correction (without the factor of alpha**2)
!drach_MVkl :: Drachmanized mass-velocity correction (without the factor of alpha**2)
!Darwinkl:: Darwin correction (without the factor of alpha**2)
!drach_Darwinkl:: Drachmanized Darwin correction (without the factor of alpha**2)
!   OOkl :: Orbit-Orbit correction (without the factor of alpha**2)
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
    real(wp),intent(in)      :: Paramk(2*Glob_np), Paraml(2*Glob_np)
    integer,intent(in)          :: Pbra(Glob_n,Glob_n),Pket(Glob_n,Glob_n)
    complex(wp),intent(out)  :: Hkl,Skl,Tkl,Vkl,MVkl,drach_MVkl,Darwinkl,drach_Darwinkl,OOkl
    complex(wp),intent(out)  :: rm2kl(Glob_n,Glob_n),rmkl(Glob_n,Glob_n)
    complex(wp),intent(out)  :: rkl(Glob_n,Glob_n),r2kl(Glob_n,Glob_n)
    complex(wp),intent(out)  :: deltarkl(Glob_n,Glob_n)
    complex(wp),intent(out)  :: drach_deltarkl(Glob_n,Glob_n)
    integer,intent(in)          :: NumCFGridPoints,NumDensGridPoints
    real(wp),intent(in)      :: CFGrid(NumCFGridPoints),DensGrid(NumDensGridPoints)
    complex(wp),intent(out)  :: CFkl(Glob_n*(Glob_n+1)/2,NumCFGridPoints)
    complex(wp),intent(out)  :: Denskl(Glob_n+1,NumDensGridPoints)
    logical,intent(in)          :: AreCorrFuncNeeded,ArePartDensNeeded

!Parameters (These are needed to declare static arrays. Using static
!arrays makes the function call a little faster in comparison with
!the case when arrays are dynamically allocated in stack)
    integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles
    integer,parameter :: nnp=nn*(nn+1)/2

!Local variables
    integer           n,np
    real(wp)       vechLk(nnp),vechBk(nnp)
    real(wp)       vechLl(nnp),vechBl(nnp)
    real(wp)       Lk(nn,nn),Ll(nn,nn),inv_Lk(nn,nn),inv_Ll(nn,nn)
    real(wp)       tAk(nn,nn),tAl(nn,nn),tBk(nn,nn),tBl(nn,nn)
    complex(wp)    tCk(nn,nn),tCl(nn,nn),tCkl(nn,nn),tCkl_copy(nn,nn)
    complex(wp)    inv_tCkl(nn,nn)
    complex(wp)    Wc1(nn,nn),Wc2(nn,nn),Wc3(nn,nn),Wc4(nn,nn),Wc5(nn,nn)
    real(wp)       Wr1(nn,nn),Wr2(nn,nn),Wr3(nn,nn),Wr4(nn,nn)
    real(wp)       rtemp,rsum,rnumb
    complex(wp)    ctemp,csum,cnumb,cmpnum,caux,caux1,tr1,tr2,tr3,tr4,tr5
    real(wp)       det_Lk,det_Ll
    complex(wp)    det_tCkl
    integer           i,j,k,kk,kkk,indx,p,q
    complex(wp)    TrCJ(nn,nn),TrCJCJ(nn,nn,nn,nn)
    complex(wp)    rmrmkl(nn,nn,nn,nn)
    real(wp)       Mass_For_Darwin(0:nn)

    n=Glob_n
    np=Glob_np

!Extracting nonlinear parameters
    vechLk(1:np)=Paramk(1:np)
    vechBk(1:np)=Paramk(np+1:np+np)
    vechLl(1:np)=Paraml(1:np)
    vechBl(1:np)=Paraml(np+1:np+np)

!First we build matrices Lk, Ll, Ak, Al, Bk, Bl
!from vechLk, vechLl, vechBk, vechBl.
    indx=0
    do i=1,n
      do j=i,n
        indx=indx+1
        Lk(i,j)=ZERO
        Lk(j,i)=vechLk(indx)
        Ll(i,j)=ZERO
        Ll(j,i)=vechLl(indx)
        tBk(j,i)=vechBk(indx)
        tBk(i,j)=vechBk(indx)
        tBl(j,i)=vechBl(indx)
        tBl(i,j)=vechBl(indx)
      enddo
    enddo

    do i=1,n
      do j=i,n
        rtemp=ZERO
        do k=1,j
          rtemp=rtemp+Lk(i,k)*Lk(j,k)
        enddo
        tAk(i,j)=rtemp
        tAk(j,i)=rtemp
        rtemp=ZERO
        do k=1,j
          rtemp=rtemp+Ll(i,k)*Ll(j,k)
        enddo
        tAl(i,j)=rtemp
        tAl(j,i)=rtemp
      enddo
    enddo

!Then we permute elements of Ak, Bk, Al and Bl to account for
!the action of the permutation matrices
!  tBk=Pbra'*Bk*Pbra  tAk=Pbra'*Ak*Pbra
!  tBl=Pket'*Bl*Pket  tAl=Pket'*Al*Pket
    do i=1,n
      do j=1,n
        rnumb=ZERO
        rtemp=ZERO
        do k=1,n
          rnumb=rnumb+Pbra(k,i)*tAk(k,j)
          rtemp=rtemp+Pbra(k,i)*tBk(k,j)
        enddo
        Wr1(i,j)=rnumb
        Wr2(i,j)=rtemp
      enddo
    enddo
    do i=1,n
      do j=i,n
        rnumb=ZERO
        rtemp=ZERO
        do k=1,n
          rnumb=rnumb+Wr1(i,k)*Pbra(k,j)
          rtemp=rtemp+Wr2(i,k)*Pbra(k,j)
        enddo
        tAk(i,j)=rnumb
        tAk(j,i)=rnumb
        tBk(i,j)=rtemp
        tBk(j,i)=rtemp
      enddo
    enddo
    do i=1,n
      do j=1,n
        rnumb=ZERO
        rtemp=ZERO
        do k=1,n
          rnumb=rnumb+Pket(k,i)*tAl(k,j)
          rtemp=rtemp+Pket(k,i)*tBl(k,j)
        enddo
        Wr1(i,j)=rnumb
        Wr2(i,j)=rtemp
      enddo
    enddo
    do i=1,n
      do j=i,n
        rnumb=ZERO
        rtemp=ZERO
        do k=1,n
          rnumb=rnumb+Wr1(i,k)*Pket(k,j)
          rtemp=rtemp+Wr2(i,k)*Pket(k,j)
        enddo
        tAl(i,j)=rnumb
        tAl(j,i)=rnumb
        tBl(i,j)=rtemp
        tBl(j,i)=rtemp
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

!Next, we form tCk, tCl, and tCkl=tCk+tCl
!Notice that array tCk will actually contain tCk^*
    do i=1,n
      do j=i,n
        cnumb=cmplx(tAk(j,i),-tBk(j,i),wp)
        tCk(j,i)=cnumb
        tCk(i,j)=cnumb
        cnumb=cmplx(tAl(j,i),tBl(j,i),wp)
        tCl(j,i)=cnumb
        tCl(i,j)=cnumb
        cnumb=tCk(j,i)+tCl(j,i)
        tCkl(j,i)=cnumb
        tCkl(i,j)=cnumb
      enddo
    enddo

!After this we can do pseudo-Cholesky factorization of tCkl
!(pseudo because it is a product of a lower triangular
!matrix and its transpose, not a product of a lower triangular
!matrix and its hermitian conjugate as the original Cholesky
!factorization)
!We will be doing all computations in tCkl_copy -- a copy of tCkl
!Matrix tCkl_copy will be destroyed. Instead we will have
!its pseudo-Cholesky factor stored in the lower triangle
    tCkl_copy(1:n,1:n)=tCkl(1:n,1:n)
    det_tCkl=COMPLEXONE
    do i=1,n
      do j=i,n
        csum=tCkl_copy(i,j)
        do k=i-1,1,-1
          csum=csum-tCkl_copy(i,k)*tCkl_copy(j,k)
        enddo
        if (i==j) then
          tCkl_copy(i,i)=sqrt(csum)
          det_tCkl=det_tCkl*csum
        else
          tCkl_copy(j,i)=csum/tCkl_copy(i,i)
        endif
      enddo
    enddo

!Inverting tCkl_copy using pseudo-Cholesky factor
!and placing the result into inv_tCkl
    do i=1,n
      tCkl_copy(i,i)=ONE/tCkl_copy(i,i)
      do j=i+1,n
        csum=COMPLEXZERO
        do k=i,j-1
          csum=csum-tCkl_copy(j,k)*tCkl_copy(k,i)
        enddo
        tCkl_copy(j,i)=csum/tCkl_copy(j,j)
      enddo
    enddo

    do i=1,n
      do j=i,n
        ctemp=COMPLEXZERO
        do k=j,n
          ctemp=ctemp+tCkl_copy(k,i)*tCkl_copy(k,j)
        enddo
        inv_tCkl(i,j)=ctemp
        inv_tCkl(j,i)=ctemp
      enddo
    enddo

!Evaluating overlap
    cnumb=abs(det_Ll*det_Lk)/det_tCkl
    Skl=Glob_2Raised3n2*cnumb*sqrt(cnumb)

!Doing multiplication tCl*invtCkl*tCk^*
    do i=1,n
      do j=1,n
        ctemp=COMPLEXZERO
        do k=1,n
          ctemp=ctemp+inv_tCkl(i,k)*tCk(k,j)
        enddo
        Wc1(i,j)=ctemp
      enddo
    enddo
    do i=1,n
      do j=1,n
        ctemp=COMPLEXZERO
        do k=1,n
          ctemp=ctemp+tCl(i,k)*Wc1(k,j)
        enddo
        Wc2(i,j)=ctemp
      enddo
    enddo

!Evaluating kinetic energy, Tkl
    ctemp=COMPLEXZERO
    do i=1,n
      do k=1,n
        ctemp=ctemp+Glob_MassMatrix(i,k)*Wc2(k,i)
      enddo
    enddo
    Tkl=6*Skl*ctemp

!Evaluating potential energy, Vkl,
!(1/r_{ij})_kl, (r_{ij})_kl, (r_{ij}^2)_kl
!and delta(r_{ij})_kl
    cnumb=(TWO/Glob_SqrtPi)*Skl
    caux=THREEHALF*Skl
    Vkl=COMPLEXZERO
    caux1=Skl/(Glob_Pi*Glob_SqrtPi)
    do i=1,n
      TrCJ(i,i)=inv_tCkl(i,i)
      csum=sqrt(TrCJ(i,i))
      rmkl(i,i)=cnumb/csum
      rkl(i,i)=cnumb*csum
      r2kl(i,i)=caux*TrCJ(i,i)
      deltarkl(i,i)=caux1/(csum*TrCJ(i,i))
      Vkl=Vkl+Glob_ScaledPseudoChargeMatrix(i,0)*rmkl(i,i)
    enddo
    do i=1,n
      do j=i+1,n
        TrCJ(i,j)=inv_tCkl(i,i)+inv_tCkl(j,j)-inv_tCkl(j,i)-inv_tCkl(j,i)
        TrCJ(j,i)=TrCJ(i,j)
        csum=sqrt(TrCJ(i,j))
        rmkl(i,j)=cnumb/csum
        rmkl(j,i)=rmkl(i,j)
        rkl(i,j)=cnumb*csum
        rkl(j,i)=rkl(i,j)
        r2kl(i,j)=caux*TrCJ(i,j)
        r2kl(j,i)=r2kl(i,j)
        deltarkl(i,j)=caux1/(csum*TrCJ(i,j))
        deltarkl(j,i)=deltarkl(i,j)
        Vkl=Vkl+Glob_ScaledPseudoChargeMatrix(i,j)*rmkl(i,j)
      enddo
    enddo
    Hkl=Tkl+Vkl

!Evaluating tr[invtCkl Jij invtCkl Jpq]
    do i=1,n
      caux=inv_tCkl(i,i)*inv_tCkl(i,i)
      TrCJCJ(i,i,i,i)=caux
      do j=i+1,n
        caux=inv_tCkl(j,i)*inv_tCkl(j,i)
        TrCJCJ(j,j,i,i)=caux
        TrCJCJ(i,i,j,j)=caux
        do p=i,n
          ctemp=inv_tCkl(p,i)-inv_tCkl(p,j)
          caux=ctemp*ctemp
          TrCJCJ(i,j,p,p)=caux
          TrCJCJ(j,i,p,p)=caux
          TrCJCJ(p,p,i,j)=caux
          TrCJCJ(p,p,j,i)=caux
          do q=p+1,n
            ctemp=inv_tCkl(q,i)-inv_tCkl(p,i)-inv_tCkl(q,j)+inv_tCkl(p,j)
            caux=ctemp*ctemp
            TrCJCJ(i,j,p,q)=caux
            TrCJCJ(j,i,p,q)=caux
            TrCJCJ(i,j,q,p)=caux
            TrCJCJ(j,i,q,p)=caux
            TrCJCJ(p,q,i,j)=caux
            TrCJCJ(p,q,j,i)=caux
            TrCJCJ(q,p,i,j)=caux
            TrCJCJ(q,p,j,i)=caux
          enddo
        enddo
      enddo
    enddo

!Evaluating (1/r_{ij}*1/r_{pq})_kl
    cnumb=4*Skl/Glob_Pi
    cmpnum=2*Skl
    do i=1,n
      do j=i,n
        do p=i,n
          do q=p,n
            if (((p==i).and.(q==j)).or.((p==j).and.(q==i))) then
              ctemp=cmpnum/TrCJ(i,j)
              rmrmkl(i,j,p,q)=ctemp
              rmrmkl(j,i,p,q)=ctemp
              rmrmkl(i,j,q,p)=ctemp
              rmrmkl(j,i,q,p)=ctemp
              rmrmkl(p,q,i,j)=ctemp
              rmrmkl(p,q,j,i)=ctemp
              rmrmkl(q,p,i,j)=ctemp
              rmrmkl(q,p,j,i)=ctemp
            else
              caux=sqrt(TrCJCJ(i,j,p,q))
              if (caux/=COMPLEXZERO) then
                caux1=caux/sqrt(TrCJ(i,j)*TrCJ(p,q))
                csum=-IMAGINARYUNIT*log(IMAGINARYUNIT*caux1+sqrt(ONE-caux1*caux1))
                ctemp=cnumb*csum/caux
              else
                ctemp=cnumb/sqrt(TrCJ(i,j)*TrCJ(p,q))
              endif
              rmrmkl(i,j,p,q)=ctemp
              rmrmkl(j,i,p,q)=ctemp
              rmrmkl(i,j,q,p)=ctemp
              rmrmkl(j,i,q,p)=ctemp
              rmrmkl(p,q,i,j)=ctemp
              rmrmkl(p,q,j,i)=ctemp
              rmrmkl(q,p,i,j)=ctemp
              rmrmkl(q,p,j,i)=ctemp
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
!Computing tCk^{\dagger}*M*tCl and placing it in Wc2
    do i=1,n
      do j=1,n
        ctemp=COMPLEXZERO
        do k=1,n
          ctemp=ctemp+Glob_MassMatrix(i,k)*tCl(k,j)
        enddo
        Wc1(i,j)=ctemp
      enddo
    enddo
    do i=1,n
      do j=1,n
        ctemp=COMPLEXZERO
        do k=1,n
          ctemp=ctemp+tCk(i,k)*Wc1(k,j)
        enddo
        Wc2(i,j)=ctemp
      enddo
    enddo
!Loop that computes all drachmanized delta(r_{ij})_kl
    do p=1,n
      do q=p,n
        ctemp=COMPLEXZERO
        do i=1,n
          ctemp=ctemp+Glob_ScaledPseudoChargeMatrix(0,i)*rmrmkl(p,q,i,i)
          do j=i+1,n
            ctemp=ctemp+Glob_ScaledPseudoChargeMatrix(i,j)*rmrmkl(p,q,i,j)
          enddo
        enddo
        if (p==q) then
          caux=2*Glob_Pi*Glob_MassMatrix(p,p)
        else
          caux=2*Glob_Pi*(Glob_MassMatrix(p,p)+Glob_MassMatrix(q,q)-Glob_MassMatrix(p,q) &
                     -Glob_MassMatrix(p,q))
        endif
        cnumb=ME_rXr_over_rij(Wc2,p,q,inv_tCkl,rmkl(p,q),TrCJ(p,q))
        drach_deltarkl(p,q)=(Glob_CurrEnergy*rmkl(p,q)-ctemp-4*cnumb)/caux
        drach_deltarkl(q,p)=drach_deltarkl(p,q)
      enddo
    enddo

!Evaluating Orbit-Orbit (OO) matrix element (without the factor of alpha**2)
    OOkl=COMPLEXZERO
!First double loop for OO
    do i=1,n
      do j=1,n
        tr1=tCl(j,i)
        tr2=tr1
        tr3=3*tCl(j,j)
        !Wc1 = Cl Eij Cl + Cl Ejj Eji Cl
        do p=1,n
          do q=1,n
            Wc1(q,p)=tCl(q,i)*tCl(p,j)+tCl(q,j)*tCl(p,i)
          enddo
        enddo
        !Wc1 = Wc1 + Ckl Ejj Cl Eij
        do p=1,n
          Wc1(p,j)=Wc1(p,j)+tCkl(p,j)*tCl(j,i)
        enddo
        !Wc1 = Wc1 + Eji Cl Ejj Cl + tr3 Eji Cl
        do p=1,n
          Wc1(j,p)=Wc1(j,p)+tCl(i,j)*tCl(p,j)+tr3*tCl(p,i)
        enddo
        !Wc2 = Ckl Ejj Cl
        do p=1,n
          do q=1,n
            Wc2(p,q)=tCkl(p,j)*tCl(q,j)
          enddo
        enddo
        !Wc3 = Eji Cl
        Wc3(1:n,1:n)=COMPLEXZERO
        do p=1,n
          Wc3(j,p)=tCl(p,i)
        enddo
        !compute integrals
        cnumb=ME_rXr_over_rij(Wc1,j,j,inv_tCkl,rmkl(j,j),TrCJ(j,j))
        ctemp=ME_rXr_rYr_over_rij(Wc2,Wc3,j,j,inv_tCkl,rmkl(j,j),TrCJ(j,j))
        caux=-6*(tr1+tr2)*rmkl(j,j)+4*cnumb-8*ctemp
        OOkl=OOkl-caux*Glob_ScaledPseudoChargeMatrix(j,0)/Glob_Mass(j+1)
      enddo
    enddo
    OOkl=OOkl/Glob_Mass(1)

!Second double loop for OO
    do i=1,n
      do j=i+1,n
        tr1=tCl(j,i)
        tr2=tr1
        tr3=3*tCl(j,j)
        !Wc1 = Cl Eij Cl + Cl Ejj Eji Cl
        do p=1,n
          do q=1,n
            Wc1(q,p)=tCl(q,i)*tCl(p,j)+tCl(q,j)*tCl(p,i)
          enddo
        enddo
        !Wc1 = Wc1 + Ckl Ejj Cl (Eij - Eii)
        do p=1,n
          Wc1(p,j)=Wc1(p,j)+tCkl(p,j)*tCl(j,i)
          Wc1(p,i)=Wc1(p,i)-tCkl(p,j)*tCl(j,i)
        enddo
        !Wc1 = Wc1 + (Eji - Eii) Cl Ejj Cl + tr3 (Eji - Eii) Cl
        do p=1,n
          cnumb=tCl(i,j)*tCl(p,j)+tr3*tCl(p,i)
          Wc1(j,p)=Wc1(j,p)+cnumb
          Wc1(i,p)=Wc1(i,p)-cnumb
        enddo
        !Wc2 = Ckl Ejj Cl
        do p=1,n
          do q=1,n
            Wc2(p,q)=tCkl(p,j)*tCl(q,j)
          enddo
        enddo
        !Wc3 = (Eji - Eii) Cl
        Wc3(1:n,1:n)=COMPLEXZERO
        do p=1,n
          Wc3(j,p)=tCl(p,i)
          Wc3(i,p)=-tCl(p,i)
        enddo
        !compute integrals
        cnumb=ME_rXr_over_rij(Wc1,i,j,inv_tCkl,rmkl(i,j),TrCJ(i,j))
        ctemp=ME_rXr_rYr_over_rij(Wc2,Wc3,i,j,inv_tCkl,rmkl(i,j),TrCJ(i,j))
        caux=-6*(tr1+tr2)*rmkl(i,j)+4*cnumb-8*ctemp
        OOkl=OOkl+caux*Glob_ScaledPseudoChargeMatrix(i,j)/(Glob_Mass(i+1)*Glob_Mass(j+1))
      enddo
    enddo
    OOkl=OOkl/2

!Evaluating mass-velocity matrix element
!Wc3=J*tCk
    do p=1,n
      ctemp=COMPLEXZERO
      do q=1,n
        ctemp=ctemp+tCk(q,p)
      enddo
      do q=1,n
        Wc3(q,p)=ctemp
      enddo
    enddo
!tr1=tr[J*tCk]
    tr1=COMPLEXZERO
    do p=1,n
      tr1=tr1+Wc3(p,p)
    enddo
!Wc1=tCk*J*tCk
    do p=1,n
      ctemp=COMPLEXZERO
      do i=1,n
        ctemp=ctemp+tCk(p,i)*Wc3(i,p)
      enddo
      Wc1(p,p)=ctemp
      do q=p+1,n
        ctemp=COMPLEXZERO
        do j=1,n
          ctemp=ctemp+tCk(q,j)*Wc3(j,p)
        enddo
        Wc1(q,p)=ctemp
        Wc1(p,q)=ctemp
      enddo
    enddo
!Wc3=J*tCl
    do p=1,n
      ctemp=COMPLEXZERO
      do q=1,n
        ctemp=ctemp+tCl(q,p)
      enddo
      do q=1,n
        Wc3(q,p)=ctemp
      enddo
    enddo
!tr2=tr[J*tCk]
    tr2=COMPLEXZERO
    do p=1,n
      tr2=tr2+Wc3(p,p)
    enddo
!Wc2=tCl*J*tCl
    do p=1,n
      ctemp=COMPLEXZERO
      do i=1,n
        ctemp=ctemp+tCl(p,i)*Wc3(i,p)
      enddo
      Wc2(p,p)=ctemp
      do q=p+1,n
        ctemp=COMPLEXZERO
        do j=1,n
          ctemp=ctemp+tCl(q,j)*Wc3(j,p)
        enddo
        Wc2(q,p)=ctemp
        Wc2(p,q)=ctemp
      enddo
    enddo
!Wc3 = tr1*Wc2+tr2*Wc1 = tr[J*tCk]*tCl*J*tCl + tr[J*tCl]*tCk*J*tCk
    do p=1,n
      Wc3(p,p)=tr1*Wc2(p,p)+tr2*Wc1(p,p)
      do q=p+1,n
        Wc3(q,p)=tr1*Wc2(q,p)+tr2*Wc1(q,p)
        Wc3(p,q)=Wc3(q,p)
      enddo
    enddo
!Wc4=inv_tCkl*Wc1
    do p=1,n
      do q=1,n
        ctemp=COMPLEXZERO
        do j=1,n
          ctemp=ctemp+inv_tCkl(q,j)*Wc1(j,p)
        enddo
        Wc4(q,p)=ctemp
      enddo
    enddo
!tr4=tr[Wc4]
    tr4=COMPLEXZERO
    do p=1,n
      tr4=tr4+Wc4(p,p)
    enddo
!Wc5=inv_tCkl*Wc2
    do p=1,n
      do q=1,n
        ctemp=COMPLEXZERO
        do j=1,n
          ctemp=ctemp+inv_tCkl(q,j)*Wc2(j,p)
        enddo
        Wc5(q,p)=ctemp
      enddo
    enddo
!tr5=tr[Wc5]
    tr5=COMPLEXZERO
    do p=1,n
      tr5=tr5+Wc5(p,p)
    enddo
!cnumb=tr[Wc4*Wc5]
    cnumb=COMPLEXZERO
    do p=1,n
      do q=1,n
        cnumb=cnumb+Wc4(p,q)*Wc5(q,p)
      enddo
    enddo
!caux=tr[inv_tCkl*Wc3]
    caux=COMPLEXZERO
    do p=1,n
      do q=1,n
        caux=caux+inv_tCkl(p,q)*Wc3(q,p)
      enddo
    enddo

    MVkl=36*Skl*(tr4*tr5+(TWO/THREE)*cnumb-caux+tr1*tr2) &
          /(Glob_Mass(1)*Glob_Mass(1)*Glob_Mass(1))

!sum for mass-velocity from 1 to n
    do i=1,n
      !Wc1=tCk*Jii*tCk
      do p=1,n
        Wc1(p,p)=tCk(p,i)*tCk(i,p)
        do q=p+1,n
          Wc1(q,p)=tCk(q,i)*tCk(i,p)
          Wc1(p,q)=Wc1(q,p)
        enddo
      enddo
      !tr1=tr[tCk*Jii]
      tr1=tCk(i,i)
      !Wc2=tCl*Jii*tCl
      do p=1,n
        Wc2(p,p)=tCl(p,i)*tCl(i,p)
        do q=p+1,n
          Wc2(q,p)=tCl(q,i)*tCl(i,p)
          Wc2(p,q)=Wc2(q,p)
        enddo
      enddo
      !tr2=tr[tCl*Jii]
      tr2=tCl(i,i)
      !caux=tr[inv_tCkl*(tr1*Wc2+tr2*Wc1)]
      caux=COMPLEXZERO
      do p=1,n
        do q=1,n
          caux=caux+inv_tCkl(p,q)*(tr1*Wc2(q,p)+tr2*Wc1(q,p))
        enddo
      enddo
      !Wc4=inv_tCkl*Wc1
      do p=1,n
        do q=1,n
          ctemp=COMPLEXZERO
          do j=1,n
            ctemp=ctemp+inv_tCkl(q,j)*Wc1(j,p)
          enddo
          Wc4(q,p)=ctemp
        enddo
      enddo
      !tr4=tr[Wc4]
      tr4=COMPLEXZERO
      do p=1,n
        tr4=tr4+Wc4(p,p)
      enddo
      !Wc5=inv_tCkl*Wc2
      do p=1,n
        do q=1,n
          ctemp=COMPLEXZERO
          do j=1,n
            ctemp=ctemp+inv_tCkl(q,j)*Wc2(j,p)
          enddo
          Wc5(q,p)=ctemp
        enddo
      enddo
      !tr5=tr[Wc5]
      tr5=COMPLEXZERO
      do p=1,n
        tr5=tr5+Wc5(p,p)
      enddo
      !cnumb=tr[Wc4*Wc5]
      cnumb=COMPLEXZERO
      do p=1,n
        do q=1,n
          cnumb=cnumb+Wc4(p,q)*Wc5(q,p)
        enddo
      enddo
      MVkl=MVkl+36*Skl*(tr4*tr5 + (TWO/THREE)*cnumb - caux + tr1*tr2) &
            /(Glob_Mass(i+1)*Glob_Mass(i+1)*Glob_Mass(i+1))
    enddo

    MVkl=-MVkl/8

!Evaluation of the drachmanized mass-velocity correction

!First term of drach_MVkl
!Form Wc1 = J/(m_0^2) + \sum Jii/(m_i^2)
    cnumb=1/(Glob_Mass(1)*Glob_Mass(1))
    do p=1,n
      Wc1(p,p)=cnumb+1/(Glob_Mass(p+1)*Glob_Mass(p+1))
      do q=p+1,n
        Wc1(p,q)=cnumb
        Wc1(q,p)=cnumb
      enddo
    enddo
!Wc2=tCl*Wc1'; tr1=tr[tCl*Wc1'];
    tr1=COMPLEXZERO
    do p=1,n
      do q=1,n
        cnumb=COMPLEXZERO
        do j=1,n
          cnumb=cnumb+tCl(j,q)*Wc1(j,p)
        enddo
        Wc2(q,p)=cnumb
      enddo
      tr1=tr1+Wc2(p,p)
    enddo
!Wc3=Wc2*tCl; remember that Wc3 is symmetric
    do p=1,n
      do q=p,n
        cnumb=COMPLEXZERO
        do j=1,n
          cnumb=cnumb+Wc2(q,j)*tCl(j,p)
        enddo
        Wc3(q,p)=cnumb
        Wc3(p,q)=cnumb
      enddo
    enddo
!tr2=tr[inv_tCkl*Wc3]
    tr2=COMPLEXZERO
    do p=1,n
      do j=1,n
        tr2=tr2+inv_tCkl(j,p)*Wc3(j,p)
      enddo
    enddo
    drach_MVkl=Skl*Glob_CurrEnergy*SIX*(tr2-tr1)/4
    caux=drach_MVkl   ! <-- delete ======================================================
!Second term of drach_MVkl; Here we will reuse Wc3 and tr1 obtained previously
    cnumb=COMPLEXZERO
    do i=1,n
      cnumb=cnumb+Glob_ScaledPseudoChargeMatrix(0,i)*  &
             (4*ME_rXr_over_rij(Wc3,i,i,inv_tCkl,rmkl(i,i),TrCJ(i,i))-6*tr1*rmkl(i,i))
    enddo
    do i=1,n
      do j=i+1,n
        cnumb=cnumb+Glob_ScaledPseudoChargeMatrix(i,j)*   &
               (4*ME_rXr_over_rij(Wc3,i,j,inv_tCkl,rmkl(i,j),TrCJ(i,j))-6*tr1*rmkl(i,j))
      enddo
    enddo
    drach_MVkl=drach_MVkl-cnumb/4
    caux1=drach_MVkl   ! <-- delete ======================================================
!Third term of drach_MVkl;
!Wc1=J*tCk; tr1=tr[Wc1]
    tr1=COMPLEXZERO
    do p=1,n
      cnumb=COMPLEXZERO
      do q=1,n
        cnumb=cnumb+tCk(q,p)
      enddo
      do q=1,n
        Wc1(q,p)=cnumb
      enddo
      tr1=tr1+cnumb
    enddo
!Wc2= tCl * \sum J_{jj} (m_0+m_j)/(m_0^2*m_j^2); tr2=tr[Wc2]
    tr2=COMPLEXZERO
    do p=1,n
      cnumb=(Glob_Mass(1)+Glob_Mass(p+1))/(Glob_Mass(1)*Glob_Mass(1)*Glob_Mass(p+1)*Glob_Mass(p+1))
      do q=1,n
        Wc2(q,p)=cnumb*tCl(q,p)
      enddo
      tr2=tr2+Wc2(p,p)
    enddo
!Wc3=inv_tCkl*tCk*Wc1; tr3=tr[Wc3]
    do p=1,n
      do q=1,n
        cnumb=COMPLEXZERO
        do j=1,n
          cnumb=cnumb+inv_tCkl(q,j)*tCk(j,p)
        enddo
        Wc5(q,p)=cnumb
      enddo
    enddo
    tr3=COMPLEXZERO
    do p=1,n
      do q=1,n
        cnumb=COMPLEXZERO
        do j=1,n
          cnumb=cnumb+Wc5(q,j)*Wc1(j,p)
        enddo
        Wc3(q,p)=cnumb
      enddo
      tr3=tr3+Wc3(p,p)
    enddo
!Wc4=inv_tCkl*Wc2*tCl; tr4=tr[Wc4]
    do p=1,n
      do q=1,n
        cnumb=COMPLEXZERO
        do j=1,n
          cnumb=cnumb+inv_tCkl(q,j)*Wc2(j,p)
        enddo
        Wc5(q,p)=cnumb
      enddo
    enddo
    tr4=COMPLEXZERO
    do p=1,n
      do q=1,n
        cnumb=COMPLEXZERO
        do j=1,n
          cnumb=cnumb+Wc5(q,j)*tCl(j,p)
        enddo
        Wc4(q,p)=cnumb
      enddo
      tr4=tr4+Wc4(p,p)
    enddo
!tr5=tr[Wc3*Wc4]
    tr5=COMPLEXZERO
    do p=1,n
      do j=1,n
        tr5=tr5+Wc3(p,j)*Wc4(j,p)
      enddo
    enddo
    drach_MVkl=drach_MVkl+Skl*(36*tr3*tr4+24*tr5-36*tr1*tr4-36*tr2*tr3+36*tr1*tr2)/8
    ctemp=drach_MVkl  ! <-- delete ======================================================
!Fourth term of drach_MVkl;
    do k=1,n-1
      !tr1=tr[Jkk*tCk]
      tr1=tCk(k,k)
      !Wc2= tCl * \sum_{j>i} J_{jj} (m_i+m_j)/(m_i^2*m_j^2); tr2=tr[Wc2]
      do p=1,k
        do q=1,n
          Wc2(q,p)=COMPLEXZERO
        enddo
      enddo
      tr2=COMPLEXZERO
      do p=k+1,n
        cnumb=(Glob_Mass(k+1)+Glob_Mass(p+1)) &
               /(Glob_Mass(k+1)*Glob_Mass(k+1)*Glob_Mass(p+1)*Glob_Mass(p+1))
        do q=1,n
          Wc2(q,p)=cnumb*tCl(q,p)
        enddo
        tr2=tr2+Wc2(p,p)
      enddo
      !Wc3=inv_tCkl*tCk*Jkk*tCk; tr3=tr[Wc3]
      do p=1,n
        do q=1,n
          Wc5(q,p)=tCk(q,k)*tCk(p,k)
        enddo
      enddo
      tr3=COMPLEXZERO
      do p=1,n
        do q=1,n
          cnumb=COMPLEXZERO
          do j=1,n
            cnumb=cnumb+inv_tCkl(q,j)*Wc5(j,p)
          enddo
          Wc3(q,p)=cnumb
        enddo
        tr3=tr3+Wc3(p,p)
      enddo
      !Wc4=inv_tCkl*Wc2*tCl; tr4=tr[Wc4]
      do p=1,n
        do q=1,n
          cnumb=COMPLEXZERO
          do j=1,n
            cnumb=cnumb+inv_tCkl(q,j)*Wc2(j,p)
          enddo
          Wc5(q,p)=cnumb
        enddo
      enddo
      tr4=COMPLEXZERO
      do p=1,n
        do q=1,n
          cnumb=COMPLEXZERO
          do j=1,n
            cnumb=cnumb+Wc5(q,j)*tCl(j,p)
          enddo
          Wc4(q,p)=cnumb
        enddo
        tr4=tr4+Wc4(p,p)
      enddo
      !tr5=tr[Wc3*Wc4]
      tr5=COMPLEXZERO
      do p=1,n
        do j=1,n
          tr5=tr5+Wc3(p,j)*Wc4(j,p)
        enddo
      enddo
      drach_MVkl=drach_MVkl+Skl*(36*tr3*tr4+24*tr5-36*tr1*tr4-36*tr2*tr3+36*tr1*tr2)/8
    enddo

!Evaluation of the Darwin correction
    Mass_For_Darwin(0)=Glob_Mass(1)
    Mass_For_Darwin(1:n)=Glob_Mass(2:n+1)
!Mass_For_Darwin(0)=10.0D20
!Mass_For_Darwin(1)=10.0D20
!Mass_For_Darwin(2)=10.0D20
    Darwinkl=COMPLEXZERO
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
!
!ME_rXr_over_rij(Wc3,i,j,inv_tCkl,rmkl(i,j),TrCJ(i,j))
!function ME_dXd_over_rij(X,i,j,invCkl,tCl,ME_1_over_rij,TrCJ)
!function ME_dXd_dYd(X,Y,invCkl,tCk,tCl,Skl)
!function ME_dXd(X,invCkl,tCl,Skl)

!First term
    Wc1(1:n,1:n)=ONE/(Glob_Mass(1)*Glob_Mass(1))
    do i=1,n
      Wc1(i,i)=Wc1(i,i)+ONE/(Glob_Mass(i+1)*Glob_Mass(i+1))
    enddo
    drach_Darwinkl=Glob_CurrEnergy*ME_dXd(Wc1,inv_tCkl,tCl,Skl)/4
!Second term
    do i=1,n
      drach_Darwinkl=drach_Darwinkl-Glob_ScaledPseudoChargeMatrix(0,i)*  &
                      ME_dXd_over_rij(Wc1,i,i,inv_tCkl,tCl,rmkl(i,i),TrCJ(i,i))/4
    enddo
    do i=1,n
      do j=i+1,n
        drach_Darwinkl=drach_Darwinkl-Glob_ScaledPseudoChargeMatrix(i,j)*  &
                        ME_dXd_over_rij(Wc1,i,j,inv_tCkl,tCl,rmkl(i,j),TrCJ(i,j))/4
      enddo
    enddo
!Third term
    Wc1(1:n,1:n)=COMPLEXONE
    Wc2(1:n,1:n)=COMPLEXZERO
    do j=1,n
      Wc2(j,j)=Wc2(j,j)+(Glob_Mass(1)+Glob_Mass(j+1)) &
                /(Glob_Mass(1)*Glob_Mass(1)*Glob_Mass(j+1)*Glob_Mass(j+1))
    enddo
    drach_Darwinkl=drach_Darwinkl+ME_dXd_dYd(Wc1,Wc2,inv_tCkl,tCk,tCl,Skl)/8
!Fourth term
    do i=1,n
      do j=i+1,n
        Wc1(1:n,1:n)=COMPLEXZERO
        Wc1(i,i)=COMPLEXONE
        Wc2(1:n,1:n)=COMPLEXZERO
        Wc2(j,j)=COMPLEXONE
        cnumb=(Glob_Mass(i+1)+Glob_Mass(j+1)) &
               /(Glob_Mass(i+1)*Glob_Mass(i+1)*Glob_Mass(j+1)*Glob_Mass(j+1))
        drach_Darwinkl=drach_Darwinkl+cnumb*ME_dXd_dYd(Wc1,Wc2,inv_tCkl,tCk,tCl,Skl)/8
      enddo
    enddo

!Evaluation of correlation functions
    if (AreCorrFuncNeeded) then
      ctemp=Skl/(Glob_Pi*Glob_SqrtPi)
      do k=1,NumCFGridPoints
        p=0
        rnumb=-CFGrid(k)*CFGrid(k)
        do i=1,n
          do j=i,n
            p=p+1
            CFkl(p,k)=ctemp*exp(rnumb/TrCJ(j,i))/(sqrt(TrCJ(j,i))*TrCJ(j,i))
          enddo
        enddo
      enddo
    endif

    if (ArePartDensNeeded) then
      ctemp=Skl/(Glob_Pi*Glob_SqrtPi)
      do k=1,NumDensGridPoints
        rnumb=-DensGrid(k)*DensGrid(k)
        do i=1,n+1
          caux=COMPLEXZERO
          do p=1,n
            caux=caux+Glob_bvc(p,i)*Glob_bvc(p,i)*inv_tCkl(p,p)
            do q=p+1,n
              caux=caux+2*Glob_bvc(q,i)*Glob_bvc(p,i)*inv_tCkl(q,p)
            enddo
          enddo
          Denskl(i,k)=ctemp*exp(rnumb/caux)/(sqrt(caux)*caux)
        enddo
      enddo
    endif

!if (ArePartDensNeeded) then
!  ctemp=Skl/(PI*Glob_SqrtPi)
!    do i=1,n+1
!      caux=COMPLEXZERO
!      do p=1,n
!        caux=caux+Glob_bvc(p,i)*Glob_bvc(p,i)*inv_tCkl(p,p)
!        do q=p+1,n
!          caux=caux+2*Glob_bvc(q,i)*Glob_bvc(p,i)*inv_tCkl(q,p)
!        enddo
!      enddo
!      caux1=ctemp/(sqrt(caux)*caux)
!      do k=1,NumDensGridPoints
!        rnumb=-DensGrid(k)*DensGrid(k)
!        Denskl(i,k)=caux1*exp(rnumb/caux)
!      enddo
!    enddo
!endif

  end subroutine MatrixElementsAll_CG_0S

  function trace(k,M)
    complex(wp) trace
    integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles
    complex(wp) M(nn,nn)
    integer i
    trace=COMPLEXZERO
    do i=1,k
      trace=trace+M(i,i)
    enddo
  end function trace

  function ME_rXr_over_rij(X,i,j,invCkl,ME_1_over_rij,TrCJ)
!Function ME_rXr_over_rij computes the following
!matrix element with complex Gaussians phi_k and phi_l:
!<phi_k| r'Xr/r_{ij} |phi_l>
!Here X is an arbitrary (i.e. nonsymmetric) complex matrix.
!Index i can be equal to j. In the latter case
!<phi_k| r'Xr/r_{i} |phi_l> is computed
!Input:
!            X  :: n x n complex matrix
!           i,j :: indices denoting i and j.
!        invCkl :: n x n complex matrix where the inverse of Ck+tCl is stored
! ME_1_over_rij :: the value of <phi_k| 1/r_{ij} |phi_l> matrix element
!          TrCJ :: the value of Tr[invCkl Jij]
!Note that n=Glob_n and nn=Glob_AllowedNumOfPseudoParticles. Although
!all arrays (both arguments and local ones) are static and have dimension
!nn x nn, only n x n subarrays are referenced.
    complex(wp) ME_rXr_over_rij
    integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles

!Arguments:
    complex(wp) X(nn,nn),invCkl(nn,nn),ME_1_over_rij,TrCJ
    integer i,j
!Local variables
    complex(wp) CXi(nn),CXj(nn)
    complex(wp) TrCX,TrCXCJ
    integer        n,m,p,q

    n=Glob_n

    TrCX=COMPLEXZERO
    do m=1,n
      do p=1,n
        TrCX=TrCX+invCkl(m,p)*X(p,m)
      enddo
    enddo

    do m=1,n
      CXi(m)=COMPLEXZERO
      do p=1,n
        CXi(m)=CXi(m)+invCkl(i,p)*X(p,m)
      enddo
    enddo
    if (j/=i) then
      do m=1,n
        CXj(m)=COMPLEXZERO
        do p=1,n
          CXj(m)=CXj(m)+invCkl(j,p)*X(p,m)
        enddo
      enddo
    endif

    if (i==j) then
      TrCXCJ=COMPLEXZERO
      do m=1,n
        TrCXCJ=TrCXCJ+CXi(m)*invCkl(m,i)
      enddo
    else
      TrCXCJ=COMPLEXZERO
      do m=1,n
        TrCXCJ=TrCXCJ+(CXi(m)-CXj(m))*(invCkl(m,i)-invCkl(m,j))
      enddo
    endif

    ME_rXr_over_rij=ME_1_over_rij*(3*TrCX-TrCXCJ/TrCJ)/2

  end function ME_rXr_over_rij

  subroutine ME_rXr_over_rij_all(X,invCkl,rmkl,TrCJ,ME)
!Subroutine ME_rXr_over_rij_all computes the following
!matrix elements with complex L=0 Gaussians phi_k and phi_l:
!<phi_k| r'Xr/r_{ij} |phi_l>    for all combinations i,j=1..n at the same time
!Index i can be equal to j. In the latter case
!<phi_k| r'Xr/r_{i} |phi_l> is computed
!Input:
!            X  :: n x n complex matrix
!        invCkl :: n x n complex matrix where the inverse of Ck+tCl is stored
!          rmkl :: the values of <phi_k| 1/r_{ij} |phi_l> matrix element
!          TrCJ :: the value of Tr[invCkl Jij]
!Output:
!            ME :: n x n complex matrix where all computed matrix elements are returned
!Note that n=Glob_n and nn=Glob_AllowedNumOfPseudoParticles. Although
!all arrays (both arguments and local ones) are static and have dimension
!nn x nn, only n x n subarrays are referenced.
    integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles

!Arguments:
    complex(wp)  X(nn,nn),invCkl(nn,nn),rmkl(Glob_n,Glob_n),TrCJ(nn,nn),ME(nn,nn)
!Local variables
    integer      i,j
    complex(wp)  CXi(nn),CXj(nn)
    complex(wp)  TrCX,TrCXCJ
    integer      n,m,p,q

    n=Glob_n

    TrCX=COMPLEXZERO
    do m=1,n
      do p=1,n
        TrCX=TrCX+invCkl(m,p)*X(p,m)
      enddo
    enddo

    do i=1,n
      do j=i,n
        do m=1,n
          CXi(m)=COMPLEXZERO
          do p=1,n
            CXi(m)=CXi(m)+invCkl(i,p)*X(p,m)
          enddo
        enddo
        if (j/=i) then
          do m=1,n
            CXj(m)=COMPLEXZERO
            do p=1,n
              CXj(m)=CXj(m)+invCkl(j,p)*X(p,m)
            enddo
          enddo
        endif

        if (i==j) then
          TrCXCJ=COMPLEXZERO
          do m=1,n
            TrCXCJ=TrCXCJ+CXi(m)*invCkl(m,i)
          enddo
        else
          TrCXCJ=COMPLEXZERO
          do m=1,n
            TrCXCJ=TrCXCJ+(CXi(m)-CXj(m))*(invCkl(m,i)-invCkl(m,j))
          enddo
        endif

        ME(i,j)=ONEHALF*rmkl(i,j)*(3*TrCX-TrCXCJ/TrCJ(i,j))
        ME(j,i)=ME(i,j)

      enddo
    enddo

  end subroutine ME_rXr_over_rij_all

  function ME_rXr_rYr_over_rij(X,Y,i,j,invCkl,ME_1_over_rij,TrCJ)
!Function ME_rXr_rYr_over_rij computes the following
!matrix element with complex Gaussians phi_k and phi_l:
!<phi_k| (r'Xr)(r'Yr)/r_{ij} |phi_l>
!Here X and Y are arbitrary (i.e. nonsymmetric) complex matrices.
!Index i can be equal to j. In the latter case
!<phi_k| (r'Xr)(r'Yr)/r_{i} |phi_l> is computed
!Input:
!            X  :: n x n complex matrix
!            Y  :: n x n complex matrix
!           i,j :: indices denoting i and j.
!        invCkl :: n x n complex matrix where the inverse of Ck+tCl is stored
! ME_1_over_rij :: the value of <phi_k| 1/r_{ij} |phi_l> matrix element
!          TrCJ :: the value of Tr[invCkl Jij]
!Note that n=Glob_n and nn=Glob_AllowedNumOfPseudoParticles. Although
!all arrays (both arguments and local ones) are static and have dimension
!nn x nn, only n x n subarrays are referenced.
    complex(wp) ME_rXr_rYr_over_rij
    integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles

!Arguments:
    complex(wp) X(nn,nn),Y(nn,nn),invCkl(nn,nn),ME_1_over_rij,TrCJ
    integer i,j
!Local variables
    complex(wp) Ys(nn,nn),Xs(nn,nn),CX(nn,nn),CY(nn,nn),CXCYi(nn),CXCYj(nn)
    complex(wp) TrCX,TrCY,TrCXCY,TrCXCJ,TrCYCJ,TrCXCYCJ
    integer        n,m,p,q

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
        CY(p,q)=COMPLEXZERO
        do m=1,n
          CY(p,q)=CY(p,q)+invCkl(p,m)*Ys(m,q)
        enddo
      enddo
    enddo
    do p=1,n
      do q=1,n
        CX(p,q)=COMPLEXZERO
        do m=1,n
          CX(p,q)=CX(p,q)+invCkl(p,m)*Xs(m,q)
        enddo
      enddo
    enddo

    do m=1,n
      CXCYi(m)=COMPLEXZERO
      do p=1,n
        CXCYi(m)=CXCYi(m)+CX(i,p)*CY(p,m)
      enddo
    enddo
    if (j/=i) then
      do m=1,n
        CXCYj(m)=COMPLEXZERO
        do p=1,n
          CXCYj(m)=CXCYj(m)+CX(j,p)*CY(p,m)
        enddo
      enddo
    endif

    TrCY=COMPLEXZERO
    TrCX=COMPLEXZERO
    do m=1,n
      TrCY=TrCY+CY(m,m)
      TrCX=TrCX+CX(m,m)
    enddo
    TrCXCY=COMPLEXZERO
    do m=1,n
      do p=1,n
        TrCXCY=TrCXCY+CX(m,p)*CY(p,m)
      enddo
    enddo

    if (i==j) then
      TrCXCJ=COMPLEXZERO
      do m=1,n
        TrCXCJ=TrCXCJ+CX(i,m)*invCkl(m,i)
      enddo
      TrCYCJ=COMPLEXZERO
      do m=1,n
        TrCYCJ=TrCYCJ+CY(i,m)*invCkl(m,i)
      enddo
      TrCXCYCJ=COMPLEXZERO
      do m=1,n
        TrCXCYCJ=TrCXCYCJ+CXCYi(m)*invCkl(m,i)
      enddo
    else
      TrCXCJ=COMPLEXZERO
      do m=1,n
        TrCXCJ=TrCXCJ+(CX(i,m)-CX(j,m))*(invCkl(m,i)-invCkl(m,j))
      enddo
      TrCYCJ=COMPLEXZERO
      do m=1,n
        TrCYCJ=TrCYCJ+(CY(i,m)-CY(j,m))*(invCkl(m,i)-invCkl(m,j))
      enddo
      TrCXCYCJ=COMPLEXZERO
      do m=1,n
        TrCXCYCJ=TrCXCYCJ+(CXCYi(m)-CXCYj(m))*(invCkl(m,i)-invCkl(m,j))
      enddo
    endif

    ME_rXr_rYr_over_rij=ME_1_over_rij*(   &
                         9*TrCY*TrCX + 6*TrCXCY           &
                         - ( 3*(TrCY*TrCXCJ + TrCX*TrCYCJ) + 4*TrCXCYCJ )/TrCJ   &
                         + 3*TrCXCJ*TrCYCJ/(TrCJ*TrCJ)    &
                         )/4

  end function ME_rXr_rYr_over_rij

  function ME_dXd_dYd(X,Y,invCkl,tCk,tCl,Skl)
    complex(wp) ME_dXd_dYd
    integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles
    complex(wp) X(nn,nn),Y(nn,nn),invCkl(nn,nn),tCk(nn,nn),tCl(nn,nn)
    complex(wp) Skl,M(nn,nn),MM(nn,nn),Z(nn,nn)
    integer i,j,n
    complex(wp) tr1,tr2,tr3,tr4,tr5
    n=Glob_n
    M(1:n,1:n)=matmul(X(1:n,1:n),tCk(1:n,1:n))
    tr1=trace(n,M)
    M(1:n,1:n)=matmul(tCl(1:n,1:n),Y(1:n,1:n))
    tr2=trace(n,M)
    M(1:n,1:n)=matmul(invCkl(1:n,1:n),matmul(tCk(1:n,1:n),matmul(X(1:n,1:n),tCk(1:n,1:n))))
    tr3=trace(n,M)
    MM(1:n,1:n)=matmul(invCkl(1:n,1:n),matmul(tCl(1:n,1:n),matmul(Y(1:n,1:n),tCl(1:n,1:n))))
    tr4=trace(n,MM)
    Z(1:n,1:n)=matmul(M(1:n,1:n),MM(1:n,1:n))
    tr5=trace(n,Z)
    ME_dXd_dYd=(36*tr3*tr4+24*tr5-36*tr1*tr4-36*tr2*tr3+36*tr1*tr2)*Skl
  end function ME_dXd_dYd

  function ME_dXd(X,invCkl,tCl,Skl)
    complex(wp) ME_dXd
    integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles
    complex(wp) X(nn,nn),invCkl(nn,nn),tCl(nn,nn),Skl,M(nn,nn)
    integer i,j,n
    complex(wp) tr1,tr2
    n=Glob_n
    M(1:n,1:n)=matmul(tCl(1:n,1:n),X(1:n,1:n))
    tr1=trace(n,M)
    M(1:n,1:n)=matmul(invCkl(1:n,1:n),matmul(tCl(1:n,1:n),matmul(X(1:n,1:n),tCl(1:n,1:n))))
    tr2=trace(n,M)
    ME_dXd=6*(tr2-tr1)*Skl
  end function ME_dXd

  function ME_dXd_over_rij(X,i,j,invCkl,tCl,ME_1_over_rij,TrCJ)
    complex(wp) ME_dXd_over_rij
    integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles
    complex(wp) X(nn,nn),invCkl(nn,nn),tCl(nn,nn),ME_1_over_rij,TrCJ
    integer i,j,n
    complex(wp) M(nn,nn),Z(nn,nn)
    complex(wp) tr1
    n=Glob_n
    M(1:n,1:n)=matmul(tCl(1:n,1:n),matmul(X(1:n,1:n),tCl(1:n,1:n)))
    Z(1:n,1:n)=matmul(tCl(1:n,1:n),X(1:n,1:n))
    tr1=trace(n,Z)
    ME_dXd_over_rij=4*ME_rXr_over_rij(M,i,j,invCkl,ME_1_over_rij,TrCJ)-6*tr1*ME_1_over_rij
  end function ME_dXd_over_rij

  function ME_rXr(X,Overlap,inv_tCkl)
    complex(wp) ME_rXr
    integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles
!Arguments:
    complex(wp) X(nn,nn),Overlap,inv_tCkl(nn,nn)
!Local variables
    complex(wp) TrCX
    integer        n,i,k

    n=Glob_n
    TrCX=COMPLEXZERO

    do i=1,n
      do k=1,n
        TrCX=TrCX+inv_tCkl(i,k)*X(k,i)
      enddo
    enddo

    ME_rXr=THREEHALF*Overlap*TrCX

  end function ME_rXr

  function ME_rXr_rYr(X,Y,Overlap,inv_tCkl)
    complex(wp) ME_rXr_rYr
    integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles
!Arguments:
    complex(wp) X(nn,nn),Y(nn,nn),Overlap,inv_tCkl(nn,nn)
!Local variables
    complex(wp) CX(nn,nn),CY(nn,nn)
    complex(wp) TrCX,TrCY,TrCXCY,ctemp1,ctemp2
    integer        n,i,j,k

    n=Glob_n
    TrCX=COMPLEXZERO
    TrCY=COMPLEXZERO

    do i=1,n
      do j=1,n
        ctemp1=COMPLEXZERO
        ctemp2=COMPLEXZERO
        do k=1,n
          ctemp1=ctemp1+inv_tCkl(j,k)*X(k,i)
          ctemp2=ctemp2+inv_tCkl(j,k)*Y(k,i)
        enddo
        CX(j,i)=ctemp1
        CY(j,i)=ctemp2
      enddo
      TrCX=TrCX+CX(i,i)
      TrCY=TrCY+CY(i,i)
    enddo

    TrCXCY=COMPLEXZERO
    do i=1,n
      do k=1,n
        TrCXCY=TrCXCY+CX(i,k)*CY(k,i)
      enddo
    enddo

    ME_rXr_rYr=Overlap*((NINE/FOUR)*TrCX*TrCY+THREEHALF*TrCXCY)

  end function ME_rXr_rYr

  function ME_rXr_rYr_rZr(X,Y,Z,Overlap,inv_tCkl)
    complex(wp) ME_rXr_rYr_rZr
    integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles
!Arguments:
    complex(wp) X(nn,nn),Y(nn,nn),Z(nn,nn),Overlap,inv_tCkl(nn,nn)
!Local variables
    complex(wp) CX(nn,nn),CY(nn,nn),CZ(nn,nn)
    complex(wp) CXCY(nn,nn),CXCZ(nn,nn)
    complex(wp) TrCX,TrCY,TrCZ,TrCXCY,TrCXCZ,TrCYCZ,TrCXCYCZ,TRCXCZCY
    complex(wp) ctemp1,ctemp2,ctemp3
    integer        n,i,j,k

    n=Glob_n
    TrCX=COMPLEXZERO
    TrCY=COMPLEXZERO
    TrCZ=COMPLEXZERO

    do i=1,n
      do j=1,n
        ctemp1=COMPLEXZERO
        ctemp2=COMPLEXZERO
        ctemp3=COMPLEXZERO
        do k=1,n
          ctemp1=ctemp1+inv_tCkl(j,k)*X(k,i)
          ctemp2=ctemp2+inv_tCkl(j,k)*Y(k,i)
          ctemp3=ctemp3+inv_tCkl(j,k)*Z(k,i)
        enddo
        CX(j,i)=ctemp1
        CY(j,i)=ctemp2
        CZ(j,i)=ctemp3
      enddo
      TrCX=TrCX+CX(i,i)
      TrCY=TrCY+CY(i,i)
      TrCZ=TrCZ+CZ(i,i)
    enddo

    TrCXCY=COMPLEXZERO
    TrCXCZ=COMPLEXZERO

    do i=1,n
      do j=1,n
        ctemp1=COMPLEXZERO
        ctemp2=COMPLEXZERO
        do k=1,n
          ctemp1=ctemp1+CX(j,k)*CY(k,i)
          ctemp2=ctemp2+CX(j,k)*CZ(k,i)
        enddo
        CXCY(j,i)=ctemp1
        CXCZ(j,i)=ctemp2
      enddo
      TrCXCY=TrCXCY+CXCY(i,i)
      TrCXCZ=TrCXCZ+CXCZ(i,i)
    enddo

    TrCYCZ=COMPLEXZERO
    TrCXCYCZ=COMPLEXZERO
    TrCXCZCY=COMPLEXZERO

    do i=1,n
      do k=1,n
        TrCYCZ=TrCYCZ+CY(i,k)*CZ(k,i)
        TrCXCYCZ=TrCXCYCZ+CXCY(i,k)*CZ(k,i)
        TrCXCZCY=TrCXCZCY+CXCZ(i,k)*CY(k,i)
      enddo
    enddo

    ME_rXr_rYr_rZr=Overlap*(  &
                    (27/EIGHT)*TrCX*TrCY*TrCZ + &
                    (NINE/FOUR)*(TrCXCY*TrCZ+TrCXCZ*TrCY+TrCYCZ*TrCX) + &
                    THREEHALF*(TrCXCYCZ+TrCXCZCY) &
                    )

  end function ME_rXr_rYr_rZr

end module matelem
