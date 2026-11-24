module matform
!Module matform contains procedures that form Hamiltonian
!and overlap matrices and related routines
  use matelem
#ifdef USE_CUDA
  use gpu_backend, only: gpu_active, gpu_build_HS, gpu_build_HS_deriv
#endif
  implicit none

contains

  subroutine StoreHS(i,j,Hij,Sij)
!Routine StoreHS stores calculated matrix elements of the
!Hamiltonian and the overlap in proper places of global
!arrays. Upon doing this, the routine normalizes
!matrix elements.
!Important comment:  i must be greater or equal to j.
!Hij and Sij are the matrix elements <phi_i|H|phi_j> and <phi_i|phi_j>,
!i.e. function i is the bra (the function that is complex conjugated).
!The matrices of the Hamiltonian and the overlap are hermitian. Their
!diagonal elements are real.
    integer,intent(in)      ::  i,j
    complex(wp),intent(in)  ::  Hij,Sij
    real(wp)             ::  f

    select case (Glob_GSEPSolutionMethod)
    case('I')
      !Only the lower triangle of array Glob_H (including the diagonal) is
      !used to store H-Glob_ApproxEnergy*S.
      !The entire array Glob_S is used to store S.
      if (i==j) then
        Glob_diagS(i)=real(Sij,wp)
        Glob_S(i,i)=ONE
        Glob_H(i,i)=real(Hij,wp)/Glob_diagS(i)-Glob_ApproxEnergy
      else
        f=1/sqrt(Glob_diagS(i)*Glob_diagS(j))
        Glob_S(i,j)=Sij*f
        Glob_S(j,i)=conjg(Glob_S(i,j))
        Glob_H(i,j)=(Hij-Glob_ApproxEnergy*Sij)*f
      endif
    case('G')
      !In the case when Glob_GSEPSolutionMethod='G' the diagonal
      !of the Hamiltonian matrix is stored in Glob_diagH
      !The diagonal of the overlap is stored in Glob_diagS
      !The lower triangles of arrays Glob_H and
      !Glob_S are used to store H and S
      if (i==j) then
        Glob_diagS(i)=real(Sij,wp)
        Glob_diagH(i)=real(Hij,wp)/Glob_diagS(i)
      else
        f=1/sqrt(Glob_diagS(i)*Glob_diagS(j))
        Glob_S(i,j)=Sij*f
        Glob_H(i,j)=Hij*f
      endif
    case('Q')
      !Q keeps the normalized, unshifted physical H and S matrices in their
      !canonical basis order. The complete lower triangles, including both
      !diagonals, are authoritative. The upper triangles are deliberately not
      !written because qrlinalg reads only the lower triangles and workproc
      !uses the upper storage only when packing a swap file.
      !
      !Glob_diagS retains the raw self-overlap. It is required to normalize
      !every off-diagonal element involving this function and later to form
      !derivatives of normalized matrix elements. Unlike G, Q also stores the
      !physical normalized Hamiltonian diagonal directly in Glob_H.
      if (i==j) then
        Glob_diagS(i)=real(Sij,wp)
        Glob_S(i,i)=ONE
        Glob_H(i,i)=real(Hij,wp)/Glob_diagS(i)
      else
        f=1/sqrt(Glob_diagS(i)*Glob_diagS(j))
        Glob_S(i,j)=Sij*f
        Glob_H(i,j)=Hij*f
      endif
    endselect
  end subroutine StoreHS

  subroutine StoreHSD(i,j,Hij,Sij,Di,Dj)
!Routine StoreHSdHdS stores calculated matrix elements of the
!Hamiltonian and the overlap as well as their derivatives in
!proper places of global arrays. Upon doing this, the routine
!normalizes matrix elements. There must be i>=j
!when routine is called. Di and Dj are the
!derivatives of Hij and Sij with respect to nonlinear parameters
!of i-th function and j-th functions:
!
!Di(1:Glob_np)              is dHijdvechLi
!Di(Glob_np+1:2*Glob_np)    is dHijdvechBi
!Di(2*Glob_np+1:3*Glob_np)  is dSijdvechLi
!Di(3*Glob_np+1:4*Glob_np)  is dSijdvechBi
!
!Dj(1:Glob_np)              is dHijdvechLj
!Dj(Glob_np+1:2*Glob_np)    is dHijdvechBj
!Dj(2*Glob_np+1:3*Glob_np)  is dSijdvechLj
!Dj(3*Glob_np+1:4*Glob_np)  is dSijdvechBj
!
!Here Hij=<phi_i|H|phi_j> and Sij=<phi_i|phi_j>, i.e. function i is
!the bra and function j is the ket. In global array Glob_D the
!function, whose parameters the derivatives are taken with respect
!to, is always the bra. For this reason the complex conjugated
!derivatives with respect to the parameters of the ket (function j)
!are stored: they are the derivatives of Hji=conjg(Hij) and
!Sji=conjg(Sij). When i=j the derivatives of the diagonal matrix
!elements (which are real) with respect to the parameters of function
!i (it is both the bra and the ket then) are stored.

    integer,intent(in)     ::  i,j
    complex(wp),intent(in) ::  Hij,Sij
    complex(wp),intent(in) ::  Di(2*Glob_npt),Dj(2*Glob_npt)
    real(wp)            ::  f

    select case (Glob_GSEPSolutionMethod)
    case('I')
      !Only the lower triangle of array Glob_H (including the diagonal) is
      !used to store H-Glob_ApproxEnergy*S.
      !The entire array Glob_S is used to store S.
      if (i==j) then
        Glob_diagS(i)=real(Sij,wp)
        Glob_S(i,i)=ONE
        Glob_H(i,i)=real(Hij,wp)/Glob_diagS(i)-Glob_ApproxEnergy
        Glob_D(1:2*Glob_npt,i-Glob_nfru,i)=TWO*real(Di(1:2*Glob_npt),wp)/Glob_diagS(i)
      else
        f=1/sqrt(Glob_diagS(i)*Glob_diagS(j))
        Glob_S(i,j)=Sij*f
        Glob_S(j,i)=conjg(Glob_S(i,j))
        Glob_H(i,j)=(Hij-Glob_ApproxEnergy*Sij)*f
        Glob_D(1:2*Glob_npt,i-Glob_nfru,j)=Di(1:2*Glob_npt)*f
        if (j>Glob_nfru) Glob_D(1:2*Glob_npt,j-Glob_nfru,i)=conjg(Dj(1:2*Glob_npt))*f
      endif
    case('G')
      !In the case when Glob_GSEPSolutionMethod='G' the diagonal
      !of the Hamiltonian matrix is stored in Glob_diagH
      !The diagonal of the overlap is stored in Glob_diagS
      !Lower triangles of arrays Glob_H and
      !Glob_S are used to store H and S
      if (i==j) then
        Glob_diagS(i)=real(Sij,wp)
        Glob_diagH(i)=real(Hij,wp)/Glob_diagS(i)
        Glob_D(1:2*Glob_npt,i-Glob_nfru,i)=TWO*real(Di(1:2*Glob_npt),wp)/Glob_diagS(i)
      else
        f=1/sqrt(Glob_diagS(i)*Glob_diagS(j))
        Glob_S(i,j)=Sij*f
        Glob_H(i,j)=Hij*f
        Glob_D(1:2*Glob_npt,i-Glob_nfru,j)=Di(1:2*Glob_npt)*f
        if (j>Glob_nfru) Glob_D(1:2*Glob_npt,j-Glob_nfru,i)=conjg(Dj(1:2*Glob_npt))*f
      endif
    endselect
  end subroutine StoreHSD

  subroutine ComputeMatElem(Nmin,Nmax)
!Subroutine ComputeMatElem computes matrix elements of the
!Hamiltonian and the overlap with basis functions whose number
!ranges from Nmin to Nmax. The derivatives of H and S are NOT
!calculated at all. Routine StoreHS is called to store
!calculated matrix elements in proper global arrays. It is
!assumed that matrix elements of the first Nmin-1 functions are
!already calculated and stored properly when the routine
!is called. Thus, only those matrix elements are computed that
!are not known yet. If all matrix elements are needed then
!one should set Nmin=1.
!  Input parameters :
!   Nmin-1 :: The number of functions whose matrix elements
!             are already known.
!     Nmax :: The number of functions whose matrix elements
!             need to be calculated.
!Since the arrays that hold the matrix elements are complex, they
!are passed to MPI routines as arrays of real numbers of twice the
!length.
!Arguments :
    integer         Nmin,Nmax
!Local variables :
    integer     k,l,i,kk,ll,ii,j,q
    integer,parameter :: nn=Glob_n
    integer     kstart,lstart,kstop,lstop,n,np,npt,nb
    complex(wp) Skl,Hkl
    complex(wp) Ssum,Hsum
!Hoisted per-function inputs of MatrixElementsHS_CG_0S (L and C=L*L'+i*B),
!precomputed once per sweep instead of rebuilt in every (pair x term) call.
    real(wp),allocatable :: Lh(:,:,:)
    complex(wp),allocatable :: Ch(:,:,:)
!These arrays are not actually used but needed for proper calling
!of subroutine MatrixElementsHS_CG_0S.
    complex(wp)  Dk(2*Glob_npt),Dl(2*Glob_npt)

    n=Glob_n
    np=Glob_np
    npt=Glob_npt
    nb=Glob_HSBuffLen
#ifdef USE_CUDA
    if (gpu_active()) then
      call gpu_build_HS(Nmin,Nmax,StoreHS)
      return
    endif
#endif

    allocate(Lh(nn,nn,Nmax),Ch(nn,nn,Nmax))
    call PrecomputeMatrices_L_C(npt,Nmax,Glob_NonlinParam(1:npt,1:Nmax),Lh,Ch)

    Glob_HklBuff1(1:nb)=ZERO
    Glob_SklBuff1(1:nb)=ZERO
    i=0

    do k=Nmin,Nmax
      do l=k,1,-1
        i=i+1
        if (i==1) then
          kstart=k
          lstart=l
        endif
        Hsum=ZERO; Ssum=ZERO
        q=(i-1)*Glob_NumYHYTerms-1
        do j=1,Glob_NumYHYTerms
          if (mod(q+j,Glob_NumOfProcs)==Glob_ProcID) then
            call MatrixElementsHS_CG_0S(Lh(1,1,k),Lh(1,1,l),Ch(1,1,k),Ch(1,1,l), &
#ifdef USE_CUDA
                                        Glob_YHYMatr(1,1,j), &
#else
                                        Glob_YHYMatr(1:n,1:n,j), &
#endif
                                        Hkl,Skl,Dk,Dl,.false.,.false.)
            Hsum=Hsum+Glob_YHYCoeff(j)*Hkl
            Ssum=Ssum+Glob_YHYCoeff(j)*Skl
          endif
        enddo
        Glob_HklBuff1(i)=Hsum
        Glob_SklBuff1(i)=Ssum
        if (i==Glob_HSBuffLen) then
          call MPI_ALLREDUCE(Glob_HklBuff1,Glob_HklBuff2,2*i, &
                             MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
          call MPI_ALLREDUCE(Glob_SklBuff1,Glob_SklBuff2,2*i, &
                             MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
          ii=0
          do kk=kstart,k
            if (kk==kstart) then
              ll=lstart
            else
              ll=kk
            endif
            if (kk==k) then
              lstop=l
            else
              lstop=1
            endif
            do while (ll>=lstop)
              ii=ii+1
              call StoreHS(kk,ll,Glob_HklBuff2(ii),Glob_SklBuff2(ii))
              ll=ll-1
            enddo
          enddo
          i=0
          Glob_HklBuff1(1:nb)=ZERO
          Glob_SklBuff1(1:nb)=ZERO
        endif
      enddo
    enddo
    if (i>0) then
      call MPI_ALLREDUCE(Glob_HklBuff1,Glob_HklBuff2,2*i, &
                         MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
      call MPI_ALLREDUCE(Glob_SklBuff1,Glob_SklBuff2,2*i, &
                         MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
      ii=0
      do kk=kstart,Nmax
        if (kk==kstart) then
          ll=lstart
        else
          ll=kk
        endif
        lstop=1
        do while (ll>=lstop)
          ii=ii+1
          call StoreHS(kk,ll,Glob_HklBuff2(ii),Glob_SklBuff2(ii))
          ll=ll-1
        enddo
      enddo
    endif

  end subroutine ComputeMatElem

  subroutine ComputeMatElemAndDeriv(Nmin,Nmax)
!Subroutine ComputeMatElem computes matrix elements of the
!Hamiltonian and the overlap as well as their derivatives with
!basis functions whose number ranges from Nmin to Nmax. Routine
!StoreHSD is called to store the calculated values in proper
!global arrays. It is assumed that matrix elements of the first
!Nmin-1 functions are already calculated and placed where
!needed when the routine is called. Thus, only those matrix
!elements are computed that are not known yet. If all matrix
!elements are needed then one should set Nmin=1. Also, one needs
!to make sure that the value of the global variables Glob_nfo and
!Glob_nfru are equal to what they should be.
!  Input parameters :
!   Nmin-1 :: The number of functions whose matrix elements
!             are already known.
!     Nmax :: The number of functions whose matrix elements
!             need to be calculated.
!Arguments :
    integer         Nmin,Nmax
!Local variables :
    integer     k,l,i,kk,ll,ii,j,q
    integer,parameter :: nn=Glob_n
    integer     kstart,lstart,kstop,lstop,n,np,npt,npt2,nb
    complex(wp) Skl,Hkl
    complex(wp) Ssum,Hsum
!Hoisted per-function inputs of MatrixElementsHS_CG_0S (see ComputeMatElem).
    real(wp),allocatable :: Lh(:,:,:)
    complex(wp),allocatable :: Ch(:,:,:)
    complex(wp) Dk(2*Glob_npt)
    complex(wp) Dl(2*Glob_npt)
    complex(wp) Dksum(2*Glob_npt)
    complex(wp) Dlsum(2*Glob_npt)
    logical     grad_l

    n=Glob_n
    np=Glob_np
    npt=Glob_npt
    npt2=npt*2
    nb=Glob_HSBuffLen
#ifdef USE_CUDA
    if (gpu_active()) then
      call gpu_build_HS_deriv(Nmin,Nmax,StoreHSD)
      return
    endif
#endif

    allocate(Lh(nn,nn,Nmax),Ch(nn,nn,Nmax))
    call PrecomputeMatrices_L_C(npt,Nmax,Glob_NonlinParam(1:npt,1:Nmax),Lh,Ch)

    Glob_HklBuff1(1:nb)=ZERO
    Glob_SklBuff1(1:nb)=ZERO
    Glob_DkBuff1(1:npt2,1:nb)=ZERO
    if (Glob_nfo>1) Glob_DlBuff1(1:npt2,1:nb)=ZERO
    i=0

    do k=Nmin,Nmax
      do l=k,1,-1
        i=i+1
        if (i==1) then
          kstart=k
          lstart=l
        endif
        Hsum=ZERO
        Ssum=ZERO
        Dksum(1:npt2)=ZERO
        if ((l>Glob_nfru).and.(l/=k)) then
          grad_l=.true.
          Dlsum(1:npt2)=ZERO
        else
          grad_l=.false.
        endif
        q=(i-1)*Glob_NumYHYTerms-1
        do j=1,Glob_NumYHYTerms
          if (mod(q+j,Glob_NumOfProcs)==Glob_ProcID) then
            call MatrixElementsHS_CG_0S(Lh(1,1,k),Lh(1,1,l),Ch(1,1,k),Ch(1,1,l), &
#ifdef USE_CUDA
                                        Glob_YHYMatr(1,1,j), &
#else
                                        Glob_YHYMatr(1:n,1:n,j), &
#endif
                                        Hkl,Skl,Dk,Dl,.true.,grad_l)
            Hsum=Hsum+Glob_YHYCoeff(j)*Hkl
            Ssum=Ssum+Glob_YHYCoeff(j)*Skl
            Dksum(1:npt2)=Dksum(1:npt2)+Glob_YHYCoeff(j)*Dk(1:npt2)
            if ((l>Glob_nfru).and.(l/=k)) Dlsum(1:npt2)=Dlsum(1:npt2)+Glob_YHYCoeff(j)*Dl(1:npt2)
          endif
        enddo
        Glob_HklBuff1(i)=Hsum
        Glob_SklBuff1(i)=Ssum
        Glob_DkBuff1(1:npt2,i)=Dksum(1:npt2)
        if ((l>Glob_nfru).and.(l/=k)) Glob_DlBuff1(1:npt2,i)=Dlsum(1:npt2)
        if (i==Glob_HSBuffLen) then
          call MPI_ALLREDUCE(Glob_HklBuff1,Glob_HklBuff2,2*i, &
                             MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
          call MPI_ALLREDUCE(Glob_SklBuff1,Glob_SklBuff2,2*i, &
                             MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
          call MPI_ALLREDUCE(Glob_DkBuff1,Glob_DkBuff2,2*i*npt2, &
                             MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
          if (Glob_nfo>1) call MPI_ALLREDUCE(Glob_DlBuff1,Glob_DlBuff2,2*i*npt2, &
                                             MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
          ii=0
          do kk=kstart,k
            if (kk==kstart) then
              ll=lstart
            else
              ll=kk
            endif
            if (kk==k) then
              lstop=l
            else
              lstop=1
            endif
            do while (ll>=lstop)
              ii=ii+1
              call StoreHSD(kk,ll,Glob_HklBuff2(ii),Glob_SklBuff2(ii), &
                            Glob_DkBuff2(1:npt2,ii),Glob_DlBuff2(1:npt2,ii))
              ll=ll-1
            enddo
          enddo
          i=0
          Glob_HklBuff1(1:nb)=ZERO
          Glob_SklBuff1(1:nb)=ZERO
          Glob_DkBuff1(1:npt2,1:nb)=ZERO
          if (Glob_nfo>1) Glob_DlBuff1(1:npt2,1:nb)=ZERO
        endif
      enddo
    enddo
    if (i>0) then
      call MPI_ALLREDUCE(Glob_HklBuff1,Glob_HklBuff2,2*i, &
                         MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
      call MPI_ALLREDUCE(Glob_SklBuff1,Glob_SklBuff2,2*i, &
                         MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
      call MPI_ALLREDUCE(Glob_DkBuff1,Glob_DkBuff2,2*i*npt2, &
                         MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
      if (Glob_nfo>1) call MPI_ALLREDUCE(Glob_DlBuff1,Glob_DlBuff2,2*i*npt2, &
                                         MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
      ii=0
      do kk=kstart,Nmax
        if (kk==kstart) then
          ll=lstart
        else
          ll=kk
        endif
        lstop=1
        do while (ll>=lstop)
          ii=ii+1
          call StoreHSD(kk,ll,Glob_HklBuff2(ii),Glob_SklBuff2(ii), &
                        Glob_DkBuff2(1:npt2,ii),Glob_DlBuff2(1:npt2,ii))
          ll=ll-1
        enddo
      enddo
    endif

  end subroutine ComputeMatElemAndDeriv

end module matform
