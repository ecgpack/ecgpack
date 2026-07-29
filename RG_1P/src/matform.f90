module matform
!Module matform contains procedures that form Hamiltonian
!and overlap matrices and related routines
  use matelem
  implicit none

contains

  subroutine StoreHS(i,j,Hij,Sij)
!Routine StoreHS stores calculated matrix elements of the
!Hamiltonian and the overlap in proper places of global
!arrays. Upon doing this, the routine normalizes
!matrix elements.
!Important comment:  i must be greater or equal to j.
    integer,intent(in)      ::  i,j
    real(wp),intent(in)  ::  Hij,Sij
    real(wp)             ::  f

    select case (Glob_GSEPSolutionMethod)
    case('I')
      !Only the lower triangle of array Glob_H (including the diagonal) is
      !used to store H-Glob_ApproxEnergy*S.
      !The entire array Glob_S is used to store S.
      if (i==j) then
        Glob_diagS(i)=Sij
        Glob_S(i,i)=ONE
        Glob_H(i,i)=Hij/Glob_diagS(i)-Glob_ApproxEnergy
      else
        f=1/sqrt(Glob_diagS(i)*Glob_diagS(j))
        Glob_S(i,j)=Sij*f
        Glob_S(j,i)=Glob_S(i,j)
        Glob_H(i,j)=(Hij-Glob_ApproxEnergy*Sij)*f
      endif
    case('G')
      !In the case when Glob_GSEPSolutionMethod='G' the diagonal
      !of the Hamiltonian matrix is stored in Glob_diagH
      !The diagonal of the overlap is stored in Glob_diagS
      !The lower triangles of arrays Glob_H and
      !Glob_S are used to store H and S
      if (i==j) then
        Glob_diagS(i)=Sij
        Glob_diagH(i)=Hij/Glob_diagS(i)
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
!Di(Glob_np+1:2*Glob_np)    is dSijdvechLi
!
!Dj(1:Glob_np)              is dHijdvechLj
!Dj(Glob_np+1:2*Glob_np)    is dSijdvechLj

    integer,intent(in)     ::  i,j
    real(wp),intent(in) ::  Hij,Sij
    real(wp),intent(in) ::  Di(2*Glob_npt),Dj(2*Glob_npt)
    real(wp)            ::  f

    select case (Glob_GSEPSolutionMethod)
    case('I')
      !Only the lower triangle of array Glob_H (including the diagonal) is
      !used to store H-Glob_ApproxEnergy*S.
      !The entire array Glob_S is used to store S.
      if (i==j) then
        Glob_diagS(i)=Sij
        Glob_S(i,i)=ONE
        Glob_H(i,i)=Hij/Glob_diagS(i)-Glob_ApproxEnergy
        Glob_D(1:2*Glob_npt,i-Glob_nfru,i)=TWO*Di(1:2*Glob_npt)/Glob_diagS(i)
      else
        f=1/sqrt(Glob_diagS(i)*Glob_diagS(j))
        Glob_S(i,j)=Sij*f
        Glob_S(j,i)=Glob_S(i,j)
        Glob_H(i,j)=(Hij-Glob_ApproxEnergy*Sij)*f
        Glob_D(1:2*Glob_npt,i-Glob_nfru,j)=Di(1:2*Glob_npt)*f
        if (j>Glob_nfru) Glob_D(1:2*Glob_npt,j-Glob_nfru,i)=Dj(1:2*Glob_npt)*f
      endif
    case('G')
      !In the case when Glob_GSEPSolutionMethod='G' the diagonal
      !of the Hamiltonian matrix is stored in Glob_diagH
      !The diagonal of the overlap is stored in Glob_diagS
      !Lower triangles of arrays Glob_H and
      !Glob_S are used to store H and S
      if (i==j) then
        Glob_diagS(i)=Sij
        Glob_diagH(i)=Hij/Glob_diagS(i)
        Glob_D(1:2*Glob_npt,i-Glob_nfru,i)=TWO*Di(1:2*Glob_npt)/Glob_diagS(i)
      else
        f=1/sqrt(Glob_diagS(i)*Glob_diagS(j))
        Glob_S(i,j)=Sij*f
        Glob_H(i,j)=Hij*f
        Glob_D(1:2*Glob_npt,i-Glob_nfru,j)=Di(1:2*Glob_npt)*f
        if (j>Glob_nfru) Glob_D(1:2*Glob_npt,j-Glob_nfru,i)=Dj(1:2*Glob_npt)*f
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
!Arguments :
    integer     Nmin,Nmax
!Local variables :
    integer     k,l,i,kk,ll,ii,j,q
    integer     kstart,lstart,kstop,lstop,n,np,np1,npt,nb
    real(wp) Paramk(Glob_AllowedNumOfPseudoParticles*(Glob_AllowedNumOfPseudoParticles+1)/2)
    real(wp) Paraml(Glob_AllowedNumOfPseudoParticles*(Glob_AllowedNumOfPseudoParticles+1)/2)
    integer     mk,ml
    real(wp) Skl,Hkl
    real(wp) Ssum,Hsum
!These arrays are not actually used but needed for proper calling
!of subroutine MatrixElementsHS_RG_0S. Thus, one can set some small size
!for them
    real(wp)  Dk(2),Dl(2)

    n=Glob_n
    np=Glob_np
    np1=np+1
    npt=Glob_npt
    nb=Glob_HSBuffLen
    Glob_HklBuff1(1:nb)=ZERO
    Glob_SklBuff1(1:nb)=ZERO
    i=0

    do k=Nmin,Nmax
      Paramk(1:npt)=Glob_NonlinParam(1:npt,k)
      mk=Glob_ZIndex(k)
      do l=k,1,-1
        i=i+1
        if (i==1) then
          kstart=k
          lstart=l
        endif
        Paraml(1:npt)=Glob_NonlinParam(1:npt,l)
        ml=Glob_ZIndex(l)
        Hsum=ZERO; Ssum=ZERO
        q=(i-1)*Glob_NumYHYTerms-1
        do j=1,Glob_NumYHYTerms
          if (mod(q+j,Glob_NumOfProcs)==Glob_ProcID) then
            call MatrixElementsHS_RG_1P(mk,ml,Paramk,Paraml,Glob_YHYMatr(1:n,1:n,j), &
                                  Hkl,Skl,Dk,Dl,.false.,.false.)
            Hsum=Hsum+Glob_YHYCoeff(j)*Hkl
            Ssum=Ssum+Glob_YHYCoeff(j)*Skl
          endif
        enddo
        Glob_HklBuff1(i)=Hsum
        Glob_SklBuff1(i)=Ssum
        if (i==Glob_HSBuffLen) then
          call MPI_ALLREDUCE(Glob_HklBuff1,Glob_HklBuff2,i, &
                             MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
          call MPI_ALLREDUCE(Glob_SklBuff1,Glob_SklBuff2,i, &
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
      call MPI_ALLREDUCE(Glob_HklBuff1,Glob_HklBuff2,i, &
                         MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
      call MPI_ALLREDUCE(Glob_SklBuff1,Glob_SklBuff2,i, &
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
    integer     Nmin,Nmax
!Local variables :
    integer     k,l,i,kk,ll,ii,j,q
    integer     kstart,lstart,kstop,lstop,n,np,npt,npt2,nb
    real(wp) Paramk(Glob_AllowedNumOfPseudoParticles*(Glob_AllowedNumOfPseudoParticles+1)/2)
    real(wp) Paraml(Glob_AllowedNumOfPseudoParticles*(Glob_AllowedNumOfPseudoParticles+1)/2)
    integer     mk,ml
    real(wp) Skl,Hkl
    real(wp) Ssum,Hsum
    real(wp) Dk(Glob_AllowedNumOfPseudoParticles*(Glob_AllowedNumOfPseudoParticles+1))
    real(wp) Dl(Glob_AllowedNumOfPseudoParticles*(Glob_AllowedNumOfPseudoParticles+1))
    real(wp) Dksum(Glob_AllowedNumOfPseudoParticles*(Glob_AllowedNumOfPseudoParticles+1))
    real(wp) Dlsum(Glob_AllowedNumOfPseudoParticles*(Glob_AllowedNumOfPseudoParticles+1))
    logical     grad_l

    n=Glob_n
    np=Glob_np
    npt=Glob_npt
    npt2=np*2
    nb=Glob_HSBuffLen

    Glob_HklBuff1(1:nb)=ZERO
    Glob_SklBuff1(1:nb)=ZERO
    Glob_DkBuff1(1:npt2,1:nb)=ZERO
    if (Glob_nfo>1) Glob_DlBuff1(1:npt2,1:nb)=ZERO
    i=0

    do k=Nmin,Nmax
      Paramk(1:npt)=Glob_NonlinParam(1:npt,k)
      mk=Glob_ZIndex(k)
      do l=k,1,-1
        i=i+1
        if (i==1) then
          kstart=k
          lstart=l
        endif
        Paraml(1:npt)=Glob_NonlinParam(1:npt,l)
        ml=Glob_ZIndex(l)
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
            call MatrixElementsHS_RG_1P(mk,ml,Paramk,Paraml,Glob_YHYMatr(1:n,1:n,j), &
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
          call MPI_ALLREDUCE(Glob_HklBuff1,Glob_HklBuff2,i, &
                             MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
          call MPI_ALLREDUCE(Glob_SklBuff1,Glob_SklBuff2,i, &
                             MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
          call MPI_ALLREDUCE(Glob_DkBuff1,Glob_DkBuff2,i*npt2, &
                             MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
          if (Glob_nfo>1) call MPI_ALLREDUCE(Glob_DlBuff1,Glob_DlBuff2,i*npt2, &
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
      call MPI_ALLREDUCE(Glob_HklBuff1,Glob_HklBuff2,i, &
                         MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
      call MPI_ALLREDUCE(Glob_SklBuff1,Glob_SklBuff2,i, &
                         MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
      call MPI_ALLREDUCE(Glob_DkBuff1,Glob_DkBuff2,i*npt2, &
                         MPI_WP,MPI_SUM,MPI_COMM_WORLD,Glob_MPIErrCode)
      if (Glob_nfo>1) call MPI_ALLREDUCE(Glob_DlBuff1,Glob_DlBuff2,i*npt2, &
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