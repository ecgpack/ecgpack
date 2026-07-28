module matelem
!Module matelem contains subroutines for computing
!matrix elements with real L=1 Gaussians.
  use globvars
  implicit none

contains

  subroutine MatrixElementsAll_RG_2D(m_k, mm_k, m_l, mm_l, vechLk, vechLl, Pbra, Pket, &
                                          Hkl, Skl, Tkl, Vkl, rm2kl, rmkl, rkl, r2kl, deltarkl, drach_deltarkl, &
                                          MVkl, drach_MVkl, Darwinkl, drach_Darwinkl, OOkl, rmrmkl, prvalkl, &
                                          NumCFGridPoints, CFGrid, &
                                          CFkl, NumDensGridPoints, DensGrid, Denskl, AreCorrFuncNeeded, ArePartDensNeeded)

!Arguments
    integer,intent(in)       :: m_k,m_l,mm_k,mm_l
    real(wp),intent(in)   :: vechLk(Glob_np), vechLl(Glob_np)
    real(wp),intent(in)   :: Pbra(Glob_n,Glob_n),Pket(Glob_n,Glob_n)
    real(wp),intent(out)  :: Hkl,Skl,Tkl,Vkl,MVkl,drach_MVkl,Darwinkl,drach_Darwinkl,OOkl
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
    logical,intent(in)       :: AreCorrFuncNeeded,ArePartDensNeeded

    integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles
    integer,parameter :: nnp=nn*(nn+1)/2

!Local variables
    integer           n,np
    integer           tvk(nn),tvl(nn),tbk(nn),tbl(nn)
    real(wp)       Lk(nn,nn),Ll(nn,nn),inv_Lk(nn,nn),inv_Ll(nn,nn),tvl8(nn)
    real(wp)       tAk(nn,nn),tAl(nn,nn),tAkl(nn,nn),rmkl1(nn,nn),Vkl1,Hkl1
    real(wp)       inv_Akk(nn,nn),inv_All(nn,nn),inv_tAkl(nn,nn), inv_tAkltAl(nn,nn)
    real(wp)       eta2(nn,nn),inv_tAkltAlM(nn,nn),eta22(nn,nn)
    real(wp)       eta(nn,nn)
    real(wp)       W1(nn,nn),W2(nn,nn),W3(nn,nn),W4(nn,nn),W5(nn,nn),W6(nn,nn),W7(nn,nn)
    real(wp)       W44(nn,nn),W55(nn,nn),W77(nn,nn)
    real(wp)       W4b(nn,nn),W5b(nn,nn),W7b(nn,nn)
    real(wp)       W44b(nn,nn),W55b(nn,nn),W77b(nn,nn),temp444,temp4444,temp444b,temp4444b
    real(wp)       inv_tAkltvl(nn),tvkinv_tAkl(nn),tvkinv_tAkltAlM(nn),u1(nn)
    real(wp)       inv_tAkltbl(nn),tbkinv_tAkl(nn),tbkinv_tAkltAlM(nn),u11(nn)
    real(wp)       temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,temp9,temp44
    real(wp)       temp55,temp66,temp55b,temp66b,temp44b,temp4b,temp5b,temp6b,temp1b,temp11b
    real(wp)       temp10,temp11,temp12,temp13,temp14,threshold,tr1, tr2, tr3, tr4,tr4vkvl,tr4bkbl,tr4bkvl,tr4vkbl
    real(wp)       det_Lk, det_Ll, det_tAkl,tau1,tau2,tau3,inv_tau3 ,V2kl,tau22,tau33,m
    integer           i,j,k,t,indx,p,q
    real(wp)       TrAJ(nn,nn),sqrtTrAJ(nn,nn),TrAJAJ(nn,nn,nn,nn)
    real(wp)       jAj(nn,nn,nn,nn),jAtvl(nn,nn),tvkAj(nn,nn),Mass_For_Darwin(0:nn)

    real(wp)   m1,m2,m3,tau331,tau332,tau333,tau334,tau221,tau222,tau223,tau224,temp38,temp338,temp36,temp336,temp37,temp337
    real(wp)   temp31,temp331,temp32,temp332,temp33,temp333,temp34,temp334,temp35,temp335,templast
    real(wp)   eta221(nn,nn),eta222(nn,nn),eta223(nn,nn),eta224(nn,nn),u111(nn)
    real(wp)   temp441,temp442,temp443,temp4440,temp4441,temp4442,h,term1,term2
    real(wp)   inv_tAkltAk(nn,nn),inv_tAkltAkM(nn,nn)
    real(wp)   inv_tAkltbk(nn),tvlinv_tAkl(nn),tvlinv_tAkltAkM(nn),tvlinv_tAkltAlM(nn)
    real(wp)   tbltbk(nn,nn),tvltvk(nn,nn),tbltvk(nn,nn),tvltbk(nn,nn),tvktbk(nn,nn),tvltbl(nn,nn)

    !Vars to calculate orbit-orbit interaction 
    real(wp)   ::  tvk_r(nn),tvl_r(nn),tbk_r(nn),tbl_r(nn)
    real(wp)   ::  EMatr(nn,nn), KMatr(nn,nn), DMatr(nn,nn), FMatr(nn,nn), GMatr(nn,nn)


    !Vars for calculating <1/rij 1/pq>
    real(wp) :: a, b, d, fij, fpq, tfij, tfpq, uij, upq, tuij, tupq, phi,  phi_sq, phi_cube, dsqab, &
    acosphi, tau, ttau, myeta, myteta, &
    commonFactor, arccosCommon, a1, a2, a3, a4, aone, atwo, arccosAns
    real(wp) :: R11, R12, R1, R21, R22, R23, R2, R31, R32, R33, R3, R4, R51, R52, R53, R5, R6, &
    ROne, RTwo, radicalCommon, radicalAns, totalAns, commonArccosRadical, xx, &
    RDZeroOne, RDZeroTwo, RDOneOne, RDOneTwo, RDTwoOne, RDTwoTwo, RDThreeOne, RDThreeTwo, &
    RDFourOne, RDFourTwo
    real(wp)   local_eps_for_xx
    !Vars to calculate delta-fucntions directly
    real(wp) :: myalpha, jijAVk, jijAVl, jijAWk, jijAWl
    !Var to set if orbit-orbit correction is needed

    local_eps_for_xx = 1.0e-6_wp





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

!Doing multiplication inv_tAkltAl=inv_tAkl*tAl, inv_tAkltAk=inv_tAkl*tAk
    do i=1,n
      do j=1,n
        temp1=ZERO
        temp2=ZERO
        do k=1,n
          temp1=temp1+inv_tAkl(j,k)*tAl(k,i)
          temp2=temp2+inv_tAkl(j,k)*tAk(k,i)
        enddo
        inv_tAkltAl(j,i)=temp1
        inv_tAkltAk(j,i)=temp2
      enddo
    enddo

!Doing multiplication inv_tAkltAlM=inv_tAkltAl*M
    do i=1,n
      do j=1,n
        temp1=ZERO
        temp2=ZERO
        do k=1,n
          temp1=temp1+inv_tAkltAl(j,k)*Glob_MassMatrix(k,i)
          temp2=temp2+inv_tAkltAk(j,k)*Glob_MassMatrix(k,i)
        enddo
        inv_tAkltAlM(j,i)=temp1
        inv_tAkltAkM(j,i)=temp2
      enddo
    enddo

!Computing tau1=tr[inv_tAkltAlM*tAk]
    tau1=ZERO
    do i=1,n
      temp1=ZERO
      do k=1,n
        temp1=temp1+inv_tAkltAlM(i,k)*tAk(k,i)
      enddo
      tau1=tau1+temp1
    enddo

!Computing tvk=Pbra'*vk and tvl=Pket'*vl
    do i=1,n
      tvk(i)=Pbra(m_k,i)
      tvl(i)=Pket(m_l,i)
      tbk(i)=Pbra(mm_k,i)
      tbl(i)=Pket(mm_l,i)
    enddo

    tvk_r = real(tvk, kind=wp)
    tvl_r = real(tvl, kind=wp)
    tbk_r = real(tbk, kind=wp)
    tbl_r = real(tbl, kind=wp)

!Compute inv_tAkltvl = inv_tAkl * tvl

    do i=1,n
      temp1=ZERO
      temp2=ZERO
      temp3=ZERO
      do j=1,n
        temp1=temp1+inv_tAkl(j,i)*tvl(j)
        temp2=temp2+inv_tAkl(j,i)*tbl(j)
        temp3=temp3+inv_tAkl(j,i)*tbk(j)
      enddo
      inv_tAkltvl(i)=temp1
      inv_tAkltbl(i)=temp2
      inv_tAkltbk(i)=temp3
    enddo

!Compute tvkinv_tAkl=tvk'*inv_tAkl
    do i=1,n
      temp1=ZERO
      temp2=ZERO
      temp3=ZERO
      do j=1,n
        temp1=temp1+tvk(j)*inv_tAkl(j,i)
        temp2=temp2+tbk(j)*inv_tAkl(j,i)
        temp3=temp3+tvl(j)*inv_tAkl(j,i)
      enddo
      tvkinv_tAkl(i)=temp1
      tbkinv_tAkl(i)=temp2
      tvlinv_tAkl(i)=temp3
    enddo

!Compute tau3=tvkinv_tAkl*tvl
    tau3=ZERO
    tau33=ZERO
    tau333=ZERO
    tau334=ZERO
    do i=1,n
      tau3=tau3+tvkinv_tAkl(i)*tvl(i)
      tau33=tau33+tbkinv_tAkl(i)*tbl(i)
      tau333=tau333+tvkinv_tAkl(i)*tbl(i)
      tau334=tau334+tbkinv_tAkl(i)*tvl(i)
    enddo
    m1=tau3*tau33
    m3=tau333*tau334
    m=m1+m3  !look formula 40 in document
!Evaluating overlap
!temp1=abs(det_Ll*det_Lk)/det_tAkl
!Skl=Glob_2Raised3n2*tau3*temp1*sqrt(temp1/(inv_Akk(m_k,m_k)*inv_All(m_l,m_l)))
    temp1=FOUR*det_tAkl*sqrt(det_tAkl)
    Skl=Glob_PiRaised3n2*m/temp1

    do i=1,n
      temp1=ZERO
      temp2=ZERO
      temp3=ZERO
      temp4=ZERO
      do j=1,n
        temp1=temp1+tvk(j)*inv_tAkltAlM(j,i)
        temp2=temp2+tbk(j)*inv_tAkltAlM(j,i)
        temp3=temp3+tvl(j)*inv_tAkltAlM(j,i)
        temp4=temp4+tvl(j)*inv_tAkltAkM(j,i)
      enddo
      tvkinv_tAkltAlM(i)=temp1
      tbkinv_tAkltAlM(i)=temp2
      tvlinv_tAkltAlM(i)=temp3
      tvlinv_tAkltAkM(i)=temp4
    enddo

!u1=tvkinv_tAkltAlM'*tAk
!tau2=u1'*inv_tAkltvl (storage for u1 as such is not needed, we use temp1=u1(i))
    tau2=ZERO
    tau22=ZERO
    tau223=ZERO
    tau224=ZERO
    do i=1,n
      temp1=ZERO
      temp2=ZERO
      do j=1,n
        temp1=temp1+tvkinv_tAkltAlM(j)*tAk(j,i)
        temp2=temp2+tbkinv_tAkltAlM(j)*tAk(j,i)
      enddo
      tau2=tau2+temp1*inv_tAkltvl(i)
      tau22=tau22+temp2*inv_tAkltbl(i)
      tau223=tau223+temp1*inv_tAkltbl(i)
      tau224=tau224+temp2*inv_tAkltvl(i)
    enddo
    h=tau3*tau22+tau33*tau2+tau333*tau224+tau334*tau223
!Evaluating the kinetic energy
    Tkl=Skl*(SIX*tau1+FOUR*h/m)
    Vkl=ZERO
    Vkl1=ZERO
    temp5=Skl*TWO
    temp1=temp5/Glob_SqrtPi
    temp8=Skl/(Glob_Pi*Glob_SqrtPi)
    do i=1,n
      temp2=inv_tAkl(i,i)
      TrAJ(i,i)=temp2
      temp3=sqrt(temp2)
      sqrtTrAJ(i,i)=temp3
      !u1'=tvk'*inv_tAkl*Jii*inv_tAkl
      do q=1,n
        temp4=ZERO
        temp44=ZERO
        temp444=ZERO
        do k=1,n
          temp4=temp4+tvk(k)*inv_tAkl(k,i)*inv_tAkl(q,i)
          temp44=temp44+tbk(k)*inv_tAkl(k,i)*inv_tAkl(q,i)
          temp444=temp444+tvl(k)*inv_tAkl(k,i)*inv_tAkl(q,i)
        enddo
        u1(q)=temp4
        u11(q)=temp44
        u111(q)=temp444
      enddo
      !eta2=u1'*tvl
      temp4=ZERO
      temp44=ZERO
      temp443=ZERO
      temp444=ZERO
      temp4440=ZERO
      temp4441=ZERO
      temp4442=ZERO
      do k=1,n
        temp4=temp4+u1(k)*tvl(k)
        temp44=temp44+u11(k)*tbl(k)
        temp443=temp443+u1(k)*tbl(k)
        temp444=temp444+u11(k)*tvl(k)
      enddo
      eta2(i,i)=temp44
      eta22(i,i)=temp4
      eta223(i,i)= temp444
      eta224(i,i)= temp443
      !temp444=temp4*temp44
      eta(i,i)=temp4*temp44+temp443*temp444
      term1=tau3*temp44+tau33*temp4+tau333*temp444+tau334*temp443
      term2=temp4*temp44+temp443*temp444
      rm2kl(i,i)=temp5*(ONE-TWO*ONETHIRD*term1/(m*temp2) + EIGHT*ONEFIFTH*term2/(THREE*m*temp2*temp2))/temp2
      rmkl(i,i)=temp1*(ONE-ONETHIRD*term1/(m*temp2) + ONEFIFTH*term2/(m*temp2*temp2))/temp3
      !rmkl(i,i)=ME_over_rij(i,i,inv_tAkl,det_tAkl,tvk,tvl,tbk,tbl)
      Vkl=Vkl+Glob_ScaledPseudoChargeMatrix(i,0)*rmkl(i,i)
      !Vkl1=Vkl1+Glob_ScaledPseudoChargeMatrix(i,0)*rmkl1(i,i)
      rkl(i,i)= temp1*temp3*(ONE+ONETHIRD*term1/(m*temp2) - ONEFIFTH*term2/(THREE*m*temp2*temp2))
      !r2kl(i,i)=Skl*THREEHALF*temp2*(ONE+TWO*ONETHIRD*term1/(m*temp2))
      temp10=temp8/(temp2*temp3)
      !deltarkl(i,i)=temp10*(ONE-term1/(m*temp2)+term2/(THREE*m*temp2*temp2))
      !prvalkl(i,i)=PI*temp10*( TWO*(Glob_EulerConst+log(temp2))*(ONE-term1/(m*temp2)+term2/(m*temp2*temp2)) &
      !+ FOUR*(term1-TWO*term2/temp2)/(THREE*m*temp2)+EIGHT*term2/(15*m*temp2*temp2) )
    enddo
    do i=1,n
      do j=i+1,n
        temp2=inv_tAkl(i,i)+inv_tAkl(j,j)-inv_tAkl(j,i)-inv_tAkl(j,i)
        TrAJ(i,j)=temp2
        TrAJ(j,i)=temp2
        temp3=sqrt(temp2)
        sqrtTrAJ(j,i)=temp3
        sqrtTrAJ(i,j)=temp3
        !u1'=tvk'*inv_tAkl*Jij*inv_tAkl
        do q=1,n
          temp4=ZERO
          temp44=ZERO
          temp444=ZERO
          temp10=ZERO
          do k=1,n
            temp4=temp4+tvk(k)*(inv_tAkl(k,i)-inv_tAkl(k,j))*(inv_tAkl(q,i)-inv_tAkl(q,j))
            temp44=temp44+tbk(k)*(inv_tAkl(k,i)-inv_tAkl(k,j))*(inv_tAkl(q,i)-inv_tAkl(q,j))
            temp444=temp444+tvl(k)*(inv_tAkl(k,i)-inv_tAkl(k,j))*(inv_tAkl(q,i)-inv_tAkl(q,j))
          enddo
          u1(q)=temp4
          u11(q)=temp44
          u111(q)=temp444
        enddo
        temp4=ZERO
        temp44=ZERO
        temp443=ZERO
        temp444=ZERO
        temp4440=ZERO
        temp4441=ZERO
        temp4442=ZERO
        do k=1,n
          temp4=temp4+u1(k)*tvl(k)
          temp44=temp44+u11(k)*tbl(k)
          temp443=temp443+u1(k)*tbl(k)
          temp444=temp444+u11(k)*tvl(k)
        enddo
        !temp444=temp4*temp44
        eta2(j,i)=temp44
        eta2(i,j)=temp44
        eta22(j,i)=temp4
        eta22(i,j)=temp4
        eta223(j,i)= temp444
        eta224(j,i)= temp443
        eta(j,i)=temp4*temp44+temp443*temp444
        eta(i,j)=temp4*temp44+temp443*temp444
        term1=tau3*temp44+tau33*temp4+tau333*temp444+tau334*temp443
        term2=temp4*temp44+temp443*temp444
        rm2kl(j,i)=temp5*(ONE-TWO*ONETHIRD*term1/(m*temp2) + EIGHT*ONEFIFTH*term2/(THREE*m*temp2*temp2))/temp2
        rm2kl(i,j)=rm2kl(j,i)
        rmkl(j,i)=temp1*(ONE-ONETHIRD*term1/(m*temp2)+ONEFIFTH*term2/(m*temp2*temp2))/temp3
        !rmkl(j,i)=ME_over_rij(i,j,inv_tAkl,det_tAkl,tvk,tvl,tbk,tbl)
        rmkl(i,j)=rmkl(j,i)
        !rmkl1(i,j)=rmkl1(j,i)
        Vkl=Vkl+Glob_ScaledPseudoChargeMatrix(i,j)*rmkl(j,i)
        !Vkl1=Vkl1+Glob_ScaledPseudoChargeMatrix(i,j)*rmkl1(j,i)
        rkl(j,i)=temp1*temp3*(ONE+ONETHIRD*term1/(m*temp2) - ONEFIFTH*term2/(THREE*m*temp2*temp2))
        rkl(i,j)=rkl(j,i)
        !r2kl(j,i)=Skl*THREEHALF*temp2*(ONE+TWO*ONETHIRD*term1/(m*temp2))
        !r2kl(i,j)=r2kl(j,i)
        temp10=temp8/(temp2*temp3)
        !deltarkl(j,i)=temp10*(ONE-term1/(m*temp2)+term2/(THREE*m*temp2*temp2))
        !deltarkl(i,j)=deltarkl(j,i)
        !prvalkl(j,i)=PI*temp10*( TWO*(Glob_EulerConst+log(temp2))*(ONE-term1/(m*temp2)+term2/(m*temp2*temp2)) &
        !  + FOUR*(term1-TWO*term2/temp2)/(THREE*m*temp2)+EIGHT*term2/(15*m*temp2*temp2) )
        !prvalkl(i,j)=prvalkl(j,i)
      enddo
    enddo
    Hkl=Tkl+Vkl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1/(r_ij*r_pq) is not implemented yet
tau = tau3                                                                                     !tau
ttau = tau33   
myeta =  tau333
myteta = tau334


  temp1=4*Skl/(3*Glob_Pi)
    do i=1,n
      do j=i,n
        do p=i,n
          do q=p,n   !try q=max(p,j),n - it may speed things up a little
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
              if (i==j .and. p/=q) then
                a = inv_tAkl(i,i)
                b = inv_tAkl(p,p) + inv_tAkl(q,q) -  inv_tAkl(p,q) - inv_tAkl(q,p)
                d = inv_tAkl(i, p) - inv_tAkl(i, q)
                fij = tvkinv_tAkl(i)
                tfij = tbkinv_tAkl(i)
                uij = inv_tAkltvl(i)
                tuij = inv_tAkltbl(i)
                fpq =  tvkinv_tAkl(p) - tvkinv_tAkl(q)                                                            !fij
                tfpq = tbkinv_tAkl(p) - tbkinv_tAkl(q)
                upq = inv_tAkltvl(p) - inv_tAkltvl(q)                                                  !uij
                tupq = inv_tAkltbl(p) - inv_tAkltbl(q)
              else if (i/=j .and. p==q) then
                a = inv_tAkl(i,i) + inv_tAkl(j,j) -  inv_tAkl(i,j) - inv_tAkl(j,i)
                b = inv_tAkl(p,p)
                d = inv_tAkl(i, p) - inv_tAkl(j, p)
                fij =  tvkinv_tAkl(i) -   tvkinv_tAkl(j)                                                     !fij
                tfij = tbkinv_tAkl(i) - tbkinv_tAkl(j)
                uij = inv_tAkltvl(i) - inv_tAkltvl(j)                                                             !uij
                tuij = inv_tAkltbl(i) - inv_tAkltbl(j)
                fpq = tvkinv_tAkl(p)
                tfpq = tbkinv_tAkl(p)
                upq = inv_tAkltvl(p)
                tupq =inv_tAkltbl(p)
              else if (i==j .and. p==q) then
                a = inv_tAkl(i,i)
                b = inv_tAkl(p,p)
                d = inv_tAkl(i,p)
                fij = tvkinv_tAkl(i)                                                          !fij
                tfij = tbkinv_tAkl(i)
                uij = inv_tAkltvl(i)                                                          !uij
                tuij = inv_tAkltbl(i)
                fpq = tvkinv_tAkl(p)
                tfpq = tbkinv_tAkl(p)
                upq = inv_tAkltvl(p)
                tupq = inv_tAkltbl(p)
              else
                a = inv_tAkl(i,i) + inv_tAkl(j,j) -  inv_tAkl(i,j) - inv_tAkl(j,i)
                b = inv_tAkl(p,p) + inv_tAkl(q,q) -  inv_tAkl(p,q) - inv_tAkl(q,p)     !a
                d = inv_tAkl(i, p) + inv_tAkl(j, q) - inv_tAkl(i, q) - inv_tAkl(j, p)   !d
                fij = tvkinv_tAkl(i) - tvkinv_tAkl(j)                                                            !fij
                tfij = tbkinv_tAkl(i) - tbkinv_tAkl(j)
                fpq = tvkinv_tAkl(p) - tvkinv_tAkl(q)                                                             !fij
                tfpq = tbkinv_tAkl(p) - tbkinv_tAkl(q)
                uij = inv_tAkltvl(i) - inv_tAkltvl(j)                                                             !uij
                tuij = inv_tAkltbl(i) - inv_tAkltbl(j)
                upq = inv_tAkltvl(p) - inv_tAkltvl(q)                                                             !uij
                tupq = inv_tAkltbl(p) - inv_tAkltbl(q)
              endif
              dsqab = d/(sqrt(a*b))
              xx = dsqab*dsqab
              phi = sqrt(ONE-xx)
              phi_sq = ONE - xx
              phi_cube = phi_sq*phi
              commonFactor = (Glob_PiRaised3n2/Glob_Pi)/(abs(det_tAkl)*sqrt(abs(det_tAkl)))

              if (xx < local_eps_for_xx) then

                !Zero-order correction
                radicalCommon = 1._wp/(45._wp*(a*b)**(2)*sqrt(a*b))
                R1 = 15._wp*a*b*(3._wp*a*b*tau - b*fij*uij - a*fpq*upq)*ttau
                R2 = b*(-15._wp*a*b*tau + 9._wp*b*fij*uij + 5._wp*a*fpq*upq)*tfij*tuij
                R3 = a*(-15._wp*a*b*tau + 5._wp*b*fij*uij + 9._wp*a*fpq*upq)*tfpq*tupq
                RDZeroOne = R1 + R2 + R3

                R1 = 15._wp*a*b*(3._wp*a*b*myeta - b*fij*tuij - a*fpq*tupq)*myteta
                R2 = b*(-15._wp*a*b*myeta + 9._wp*b*fij*tuij + 5._wp*a*fpq*tupq)*tfij*uij
                R3 = a*(-15._wp*a*b*myeta + 5._wp*b*fij*tuij + 9._wp*a*fpq*tupq)*tfpq*upq
                RDZeroTwo = R1 + R2 + R3

                !First-order correction
                R1 = (FIVE*a*b*ttau-THREE*b*tfij*tuij-THREE*a*tfpq*tupq)*(fij*upq+fpq*uij)
                R2 = (FIVE*a*b*tau-THREE*b*fij*uij-THREE*a*fpq*upq)*(tfij*tupq+tfpq*tuij)
                RDOneOne = (R1 + R2)*d

                R1 = (FIVE*a*b*myteta-THREE*b*tfij*uij-THREE*a*tfpq*upq)*(fij*tupq+fpq*tuij)
                R2 = (FIVE*a*b*myeta-THREE*b*fij*tuij-THREE*a*fpq*tupq)*(tfij*upq+tfpq*uij)
                RDOneTwo = (R1 + R2)*d

                !Second-order correction
                R1 = 25._wp*a*b*ttau*(a*b*tau - b*fij*uij - a*fpq*upq)
                R2 = b*tfij*tuij*(-25._wp*a*b*tau+25._wp*b*fij*uij+21._wp*a*fpq*upq)
                R3 = 6._wp*a*b*tfij*tupq*(fpq*uij + fij*upq)
                R4 = -25._wp*(a**2)*b*tau*tfpq*tupq
                R5 = 3._wp*a*b*tfpq*fij*(2._wp*upq*tuij + 7._wp*uij*tupq)
                R6 = a*tfpq*fpq*(6._wp*b*uij*tuij + 25._wp*a*upq*tupq)
                RDTwoOne = (R1 + R2 + R3 + R4 + R5 + R6)*d**2/(150._wp*(a*b)**3*sqrt(a*b))

                R1 = 25._wp*a*b*myteta*(a*b*myeta - b*fij*tuij - a*fpq*tupq)
                R2 = b*tfij*uij*(-25._wp*a*b*myeta+25._wp*b*fij*tuij+21._wp*a*fpq*tupq)
                R3 = 6._wp*a*b*tfij*upq*(fpq*tuij + fij*tupq)
                R4 = -25._wp*(a**2)*b*myeta*tfpq*upq
                R5 = 3._wp*a*b*tfpq*fij*(2._wp*tupq*uij + 7._wp*tuij*upq)
                R6 = a*tfpq*fpq*(6._wp*b*tuij*uij + 25._wp*a*tupq*upq)
                RDTwoTwo = (R1 + R2 + R3 + R4 + R5 + R6)*d**2/(150._wp*(a*b)**3*sqrt(a*b))

                !Third-order correction
                R1 = (a*b*ttau-b*tfij*tuij-a*tfpq*tupq)*(fij*upq+fpq*uij)
                R2 = (a*b*tau-b*fij*uij-a*fpq*upq)*(tfij*tupq+tfpq*tuij)
                RDThreeOne = (R1 + R2)*d**3/(10._wp*(a*b)**3*sqrt(a*b))

                R1 = (a*b*myteta-b*tfij*uij-a*tfpq*upq)*(fij*tupq+fpq*tuij)
                R2 = (a*b*myeta-b*fij*tuij-fpq*tupq)*(tfij*upq+tfpq*uij)
                RDThreeTwo = (R1 + R2)*d**3/(10._wp*(a*b)**3*sqrt(a*b))

                !Fourth-order correction
                R1 = 7._wp*a*b*ttau*(3._wp*a*b*tau - 5._wp*b*fij*uij - 5._wp*a*fpq*upq)
                R2 = b*tfij*tuij*(-35._wp*a*b*tau + 49._wp*b*fij*uij + 45._wp*a*fpq*upq)
                R3 = 20._wp*a*b*tfij*tupq*(fpq*uij+fij*upq)
                R4 = a*tfpq*fpq*(20._wp*b*uij*tuij + 49._wp*a*upq*tupq)
                R5 = -35._wp*(a**2)*b*tau*tupq*tfpq
                R6 = 5._wp*a*b*tfpq*fij*(4._wp*upq*tuij + 9._wp*uij*tupq)
                RDFourOne = (R1 + R2 + R3 + R4 + R5 + R6)*(d**4)/(280._wp*((a*b)**4)*sqrt(a*b))

                R1 = 7._wp*a*b*myteta*(3._wp*a*b*myeta - 5._wp*b*fij*tuij - 5._wp*a*fpq*tupq)
                R2 = b*tfij*uij*(-35._wp*a*b*myeta + 49._wp*b*fij*tuij + 45._wp*a*fpq*tupq)
                R3 = 20._wp*a*b*tfij*upq*(fpq*tuij+fij*tupq)
                R4 = a*tfpq*fpq*(20._wp*b*tuij*uij + 49._wp*a*tupq*upq)
                R5 = -35._wp*(a**2)*b*myeta*upq*tfpq
                R6 = 5._wp*a*b*tfpq*fij*(4._wp*tupq*uij + 9._wp*tuij*upq)
                RDFourTwo = (R1 + R2 + R3 + R4 + R5 + R6)*(d**4)/(280._wp*((a*b)**4)*sqrt(a*b))

                totalAns = commonFactor*radicalCommon*(RDZeroOne + RDZeroTwo + RDOneOne + RDOneTwo) + &
                           commonFactor*(RDTwoOne + RDTwoTwo + RDThreeOne + RDThreeTwo + RDFourOne + RDFourTwo)

              else
                acosphi=asin(abs(dsqab))

                commonArccosRadical = 1._wp/(15._wp*abs(d)**3*sqrt(a*b)*(a*b)**3*phi_cube)
            !!! Calculation of arccos part  !!!
                arccosCommon =-(a*b)**3*sqrt(a*b)*phi_cube

                a1 = 5._wp*d*ttau*(fpq*uij + fij*upq - 3._wp*d*tau)
                a2 = (5._wp*d*tau - 3._wp*fij*upq - 3._wp*fpq*uij)*(tfij*tupq + tfpq*tuij)
                a3 = 2._wp*fij*uij*tfpq*tupq
                a4 = 2._wp*fpq*upq*tfij*tuij
                aone = a1 + a2 + a3 + a4

                !Second term
                a1 = 5._wp*d*myteta*(fpq*tuij + fij*tupq - 3._wp*d*myeta)
                a2 = (5._wp*d*myeta - 3._wp*fij*tupq - 3._wp*fpq*tuij)*(tfij*upq + tfpq*uij)
                a3 = 2._wp*fij*tuij*tfpq*upq
                a4 = 2._wp*fpq*tupq*tfij*uij
                atwo = a1 + a2 + a3 + a4

                arccosAns = arccosCommon * (aone + atwo) * acosphi

            !!! Calculation of radical part  !!!
                radicalCommon=abs(d)

                !First term
                R11 = (a*b*upq - b*d*uij)*fij
                R12 = (a*b*uij - a*d*upq)*fpq
                R1 = 5._wp*(a**2)*(b**2)*d*(phi_sq)*ttau*(R11+R12)
                R21 = (d**2)*((2._wp*(d**2) - 3._wp*a*b)*uij + a*d*upq)*fij
                R22 = 5._wp*(a**2)*b*(d**2)*(phi_sq)*tau
                R23 = a*((d**3)*uij + a*((d**2)-2._wp*a*b)*upq)*fpq
                R2 = -(b**2)*tfij*tuij*(R21 + R22 + R23)
                R31 = -b*((d**3)*uij + a*(3._wp*a*b-4._wp*(d**2))*upq)*fij
                R32 = 5._wp*(a**2)*(b**2)*d*(phi_sq)*tau
                R33 = -a*(b*(3._wp*a*b-4._wp*(d**2))*uij + (d**3)*upq)*fpq
                R3 = a*b*tfij*tupq*(R31 + R32 + R33)
                R4 = a*b*tfpq*tuij*(R31 + R32 + R33)
                R51 = (d**2)*((2._wp*(d**2) - 3._wp*a*b)*upq + b*d*uij)*fpq
                R52 = 5._wp*a*(b**2)*(d**2)*(phi_sq)*tau
                R53 = b*((d**3) * upq + b*((d**2) - 2._wp*a*b)*uij)*fij
                R5 = -(a**2)*tfpq*tupq*(R51 + R52 + R53)
                ROne = R1 + R2 + R3 + R4 + R5

                !Second term
                R11 = (a*b*tupq - b*d*tuij)*fij
                R12 = (a*b*tuij - a*d*tupq)*fpq
                R1 = 5._wp*(a**2)*(b**2)*d*(phi_sq)*myteta*(R11 + R12)
                R21 = (d**2)*((2._wp*(d**2) - 3._wp*a*b)*tuij + a*d*tupq)*fij
                R22 = 5._wp*(a**2)*b*(d**2)*(phi_sq)*myeta
                R23 = a*((d**3)*tuij + a*((d**2)-2._wp*a*b)*tupq)*fpq
                R2 = -(b**2)*tfij*uij*(R21 + R22 + R23)
                R31 = -b*((d**3)*tuij + a*(3._wp*a*b-4._wp*(d**2))*tupq)*fij
                R32 = 5._wp*(a**2)*(b**2)*d*(phi_sq)*myeta
                R33 = -a*(b*(3._wp*a*b-4._wp*(d**2))*tuij + (d**3)*tupq)*fpq
                R3 = a*b*tfij*upq*(R31 + R32 + R33)
                R4 = a*b*tfpq*uij*(R31 + R32 + R33)
                R51 = (d**2)*((2._wp*(d**2) - 3._wp*a*b)*tupq + b*d*tuij)*fpq
                R52 = 5._wp*a*(b**2)*(d**2)*(phi_sq)*myeta
                R53 = b*((d**3) * tupq + b*((d**2) - 2._wp*a*b)*tuij)*fij
                R5 = -(a**2)*tfpq*upq*(R51 + R52 + R53)
                RTwo = R1 + R2 + R3 + R4 + R5

                radicalAns = radicalCommon * (ROne + RTwo)

                totalAns = commonFactor * commonArccosRadical * (arccosAns + radicalAns)

              endif
              rmrmkl(i,j,p,q)=totalAns
              rmrmkl(j,i,p,q)=totalAns
              rmrmkl(i,j,q,p)=totalAns
              rmrmkl(j,i,q,p)=totalAns
              rmrmkl(p,q,i,j)=totalAns
              rmrmkl(p,q,j,i)=totalAns
              rmrmkl(q,p,i,j)=totalAns
              rmrmkl(q,p,j,i)=totalAns
            endif
          enddo
        enddo
      enddo
    enddo



!evaluate delta-functions directly
!V ---- tau
!tV --- myeta
!W ---- ttau
!tW --- myteta
deltarkl = ZERO
do i = 1,n 
  do j = i,n
    if (i==j) then 
      myalpha = inv_tAkl(i,i)
      jijAvk = tvkinv_tAkl(i)
      jijAvl = inv_tAkltvl(i)
      jijAwk = tbkinv_tAkl(i)
      jijAwl = inv_tAkltbl(i)
    else
      myalpha = inv_tAkl(i,i) + inv_tAkl(j,j) - inv_tAkl(i,j) - inv_tAkl(j,i)
      jijAvk = tvkinv_tAkl(i) - tvkinv_tAkl(j)
      jijAvl = inv_tAkltvl(i) - inv_tAkltvl(j)
      jijAwk = tbkinv_tAkl(i) - tbkinv_tAkl(j)
      jijAwl = inv_tAkltbl(i) - inv_tAkltbl(j)
    endif
    R1 = (myalpha**2)*(tau*ttau + myeta*myteta)
    R2 = -myalpha*(jijAvk*jijAvl*ttau + tau*jijAwk*jijAwl+ &
    jijAvk*jijAwl*myteta + myeta*jijAwk*jijAvl)
    R3 = TWO*(jijAvk*jijAvl*jijAwk*jijAwl)

    deltarkl(i,j) = (R1 + R2 + R3)*&
    Glob_Piraised3n2/(FOUR*Glob_Pi*Glob_SqrtPi*(myalpha**3)*sqrt(myalpha)*det_tAkl*sqrt(det_tAkl))

    deltarkl(j,i) = deltarkl(i,j)
  enddo
enddo

!evaluate drachmanized delta-function and V2kl operator
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
        temp4=ZERO
        temp5=ZERO
        if (p==q) then
          temp4=2*Glob_Pi*Glob_MassMatrix(p,p)
          temp5=Glob_ScaledPseudoChargeMatrix(0,p)
        else
          temp4=2*Glob_Pi*(Glob_MassMatrix(p,p)+Glob_MassMatrix(q,q) &
                      -Glob_MassMatrix(p,q)-Glob_MassMatrix(p,q))
          temp5=Glob_ScaledPseudoChargeMatrix(p,q)
        endif

        !temp2=ME_rXr_over_rij(W2,p,q,inv_tAkl,rmkl(p,q),TrAJ(p,q))
        !temp2=ZERO
        temp2 = ME_d_X_over_rij_d(p,q,Glob_dmvM,tAk,tAl,inv_tAkl,det_tAkl,tvk,tvl,tbk,tbl, &
                                  tvkinv_tAkl, tbkinv_tAkl, inv_tAkltvl, inv_tAkltbl)
        drach_deltarkl(p,q)=(0.5*(Glob_CurrEnergy0+Glob_CurrEnergy1)*rmkl(p,q)-temp1-temp2)/temp4
        !drach_deltarkl(p,q)=temp2
        drach_deltarkl(q,p)=drach_deltarkl(p,q)

        V2kl=V2kl+temp5*temp1
      enddo
    enddo


!Evaluating vector-matrix-vector products
!j^{ij}' inv_tAkl tvl
!tvk' inv_tAkl j^{ij}
    do j=1,n
      do i=1,n
        if (i==j) then
          jAtvl(i,i)=inv_tAkltvl(i)
          tvkAj(i,i)=tvkinv_tAkl(i)
        else
          jAtvl(i,j)=inv_tAkltvl(i)-inv_tAkltvl(j)
          tvkAj(i,j)=tvkinv_tAkl(i)-tvkinv_tAkl(j)
        endif
      enddo
    enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

!Mass-velocity correction
    inv_tau3=1/tau3
    W1(1:n,1:n)=ONE
    temp1=Glob_Mass(1)*Glob_Mass(1)*Glob_Mass(1)
!MVkl=ME_dWd2(W1,tAk,tAl,inv_tAkl,tvk,tvl,inv_tAkltvl,tvkinv_tAkl,inv_tau3,Skl)/temp1
    MVkl=dXddYd(W1,W1,tvk,tbk,tvl,tbl,tAl,tAk,inv_tAkl,det_tAkl,tau3,tau33,tau333,tau334,inv_tAkltAl,inv_tAkltAk)/temp1
    W1(1:n,1:n)=ZERO
    do i=1,n
      W1(i,i)=ONE
      temp1=Glob_Mass(i+1)*Glob_Mass(i+1)*Glob_Mass(i+1)
      !MVkl=MVkl+ME_dWd2(W1,tAk,tAl,inv_tAkl,tvk,tvl,inv_tAkltvl,tvkinv_tAkl,inv_tau3,Skl)/temp1
      MVkl=MVkl+dXddYd(W1,W1,tvk,tbk,tvl,tbl,tAl,tAk,inv_tAkl,det_tAkl,tau3,tau33,tau333,tau334,inv_tAkltAl,inv_tAkltAk)/temp1
      W1(i,i)=ZERO
    enddo
    MVkl=-MVkl/8


    temp1=dXddYd(Glob_dmvM,Glob_dmvM,tvk,tbk,tvl,tbl,tAl,tAk,inv_tAkl,det_tAkl,tau3,tau33,tau333,tau334,inv_tAkltAl,inv_tAkltAk)&
           -V2kl-Glob_CurrEnergy1*Glob_CurrEnergy0*Skl+(Glob_CurrEnergy1+Glob_CurrEnergy0)*Vkl

    drach_MVkl = temp1*Glob_dmva2 + MVkl

    !Evaluation of the orbit-orbit interaction
    OOkl = ZERO
    do i=1,n
      !init auxiliary matrices
      EMatr = ZERO 
      EMatr(i,i) = ONE
      call symmetrize_matrix(EMatr)
      KMatr = ZERO 
      KMatr(i,i) = ONE
      DMatr = KMatr
      FMatr = KMatr
      GMatr = KMatr

      OOkl = OOkl - ONEHALF*Glob_ScaledPseudoChargeMatrix(0,i)/&
      (Glob_Mass(1)*Glob_Mass(i+1))*&
      (ME_over_rij_dXd(i,i,EMatr,tAl,inv_tAkl,det_tAkl,tvk_r,tvl_r,tbk_r,tbl_r) - &
      ME_KDFG(i,i,KMatr,DMatr,FMatr,GMatr,tAk,tAl,inv_tAkl,det_tAkl,tvk_r,tvl_r,tbk_r,tbl_r))
    enddo
    do i=1,n 
      do j=1,n
        if (i==j) cycle
        !init auxiliary matrices
        EMatr = ZERO 
        EMatr(i,j) = ONE
        call symmetrize_matrix(EMatr)
        KMatr = ZERO 
        KMatr(i,i) = ONE 
        DMatr = KMatr 
        FMatr = KMatr 
        GMatr = ZERO
        GMatr(i,j) = ONE

        OOkl = OOkl - ONEHALF*Glob_ScaledPseudoChargeMatrix(0,i)/&
        (Glob_Mass(1)*Glob_Mass(i+1))*&
        (ME_over_rij_dXd(i,i,EMatr,tAl,inv_tAkl,det_tAkl,tvk_r,tvl_r,tbk_r,tbl_r) - &
        ME_KDFG(i,i,KMatr,DMatr,FMatr,GMatr,tAk,tAl,inv_tAkl,det_tAkl,tvk_r,tvl_r,tbk_r,tbl_r)) 
      enddo
    enddo
    do i=1,n
      do j=i+1,n
        !init auxiliary matrices
        EMatr = ZERO 
        EMatr(i,j) = ONE
        call symmetrize_matrix(EMatr)
        KMatr = ZERO 
        KMatr(i,j) = ONE
        KMatr(j,j) = -ONE 
        DMatr = ZERO 
        DMatr(j,i) = ONE 
        FMatr = ZERO 
        FMatr(i,i) = ONE 
        GMatr = ZERO 
        GMatr(j,j) = ONE

        OOkl = OOkl + ONEHALF*Glob_ScaledPseudoChargeMatrix(i,j)/&
        (Glob_Mass(i+1)*Glob_Mass(j+1))*&
        (ME_over_rij_dXd(i,j,EMatr,tAl,inv_tAkl,det_tAkl,tvk_r,tvl_r,tbk_r,tbl_r) + &
        ME_KDFG(i,j,KMatr,DMatr,FMatr,GMatr,tAk,tAl,inv_tAkl,det_tAkl,tvk_r,tvl_r,tbk_r,tbl_r))
      enddo
    enddo




!Evaluation of correlation functions
    if (AreCorrFuncNeeded) then
      temp1=Skl/(Glob_Pi*Glob_SqrtPi)
      p=0
      do i=1,n
        do j=i,n
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
      do i=1,n+1
        temp2=ZERO
        do p=1,n
          temp2=temp2+Glob_bvc(p,i)*Glob_bvc(p,i)*inv_tAkl(p,p)
          do q=p+1,n
            temp2=temp2+2*Glob_bvc(q,i)*Glob_bvc(p,i)*inv_tAkl(q,p)
          enddo
        enddo
        temp3=ZERO
        temp4=ZERO
        do p=1,n
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



  end subroutine MatrixElementsAll_RG_2D



  function ME_KDFG(rindexI,rindexJ,KK,DD,FF,GG,tAk,tAl,inv_tAkl,det_tAkl,tvk,tvl,twk,twl)
    
  real(wp)   ME_KDFG
  integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles
  !Arguments:
  real(wp)   KK(nn,nn), DD(nn,nn), FF(nn,nn), GG(nn,nn), &
  tAl(nn,nn), tAk(nn,nn), inv_tAkl(nn,nn), det_tAkl
  real(wp)       tvk(nn),tvl(nn),twk(nn),twl(nn)
  integer       rindexI,rindexJ
  
  !Local variables: 
  integer     ::  i,j,k,n
  real(wp) ::  KG(nn,nn), KGCl(nn,nn), CkD(nn,nn),CkDF(nn,nn), CkDFCl(nn,nn), &
  KGCltFtDCk(nn,nn), DF(nn,nn), DFCl(nn,nn), KGCltFtD(nn,nn), &
  tGtKDF(nn,nn), KGClDFCl(nn,nn), KGClDF(nn,nn), KGCltFtDCl(nn,nn), &
  ClDFCl(nn,nn), ClDF(nn,nn)
  real(wp) ::  KGCl_s(nn,nn), CkDFCl_s(nn,nn), KGCltFtDCk_s(nn,nn), &
  tGtKDF_s(nn,nn), KGClDFCl_s(nn,nn), &
  KGCltFtDCl_s(nn,nn), ClDFCl_s(nn,nn)
  real(wp) ::  CkDFvl(nn), CkDFwl(nn), KGvl(nn), KGwl(nn), &
  KGCltFtDvk(nn), KGCltFtDwk(nn), CltFtDvk(nn),  CltFtDwk(nn), &
  DFvl(nn), DFwl(nn), KGClDFvl(nn), KGClDFwl(nn), ClDFvl(nn), ClDFwl(nn), &
  KGCltFtDvl(nn),  KGCltFtDwl(nn), CltFtDvl(nn), CltFtDwl(nn)
  real(wp) ::  commonFactor, gamma, temp, temp1, temp2, temp3
  real(wp) ::  term1, term2, term3, term4, term5, term6, term7, term8, &
  term9, term10, term11, term12, term13, term14, term15, term16, term17, &
  term18, term19, term20, term21, term22, term23, term24, &
  expr1, expr2, expr3, expr4, expr5, expr6, expr7, expr8, &
  expr9, expr10, expr11, expr12, expr13, expr14, expr15, expr16, expr17, &
  expr18, expr19, expr20, expr21, expr22, &
  bigTerm1, bigTerm2, bigTerm3
  real(wp) ::  vkDFvl, vkDFwl, wkDFvl, wkDFwl, trDFCl

  
  !Term1
  n = Glob_n
  ME_KDFG = ZERO

  KGCltFtDCk = ZERO
  KGCltFtDCk_s = ZERO
 
  !Build KG, KGCl matrices
  KG = ZERO
  do i=1,n 
    do j=1,n 
      temp = ZERO
      do k=1,n
        temp = temp + KK(i,k)*GG(k,j)
      enddo
      KG(i,j) = temp
    enddo
  enddo
  KGCl = ZERO
  do i=1,n 
    do j=1,n 
      temp = ZERO
      do k=1,n
        temp = temp + KG(i,k)*tAl(k,j)
      enddo
      KGCl(i,j) = temp
    enddo
  enddo

  !Build CkD, CkDF, CkDFCl matrices
  !tFtD, CltFtD
  CkD = ZERO
  do i=1,n 
    do j=1,n 
      temp = ZERO
      do k=1,n
        temp = temp + tAk(i,k)*DD(k,j)
      enddo
      CkD(i,j) = temp
    enddo
  enddo
  CkDF = ZERO
  do i=1,n 
    do j=1,n 
      temp = ZERO
      do k=1,n
        temp = temp + CkD(i,k)*FF(k,j)
      enddo
      CkDF(i,j) = temp
    enddo
  enddo
  CkDFCl = ZERO
  do i=1,n 
    do j=1,n 
      temp = ZERO
      do k=1,n
        temp = temp + CkDF(i,k)*tAl(k,j)
      enddo
      CkDFCl(i,j) = temp
    enddo
  enddo
  DF = ZERO 
  do i=1,n 
    do j=1,n 
      temp = ZERO
      do k=1,n
        temp = temp + DD(i,k)*FF(k,j)
      enddo
      DF(i,j) = temp
    enddo
  enddo
  DFCl = ZERO
  do i=1,n 
    do j=1,n 
      temp = ZERO
      do k=1,n
        temp = temp + DF(i,k)*tAl(k,j)
      enddo
      DFCl(i,j) = temp
    enddo
  enddo
  KGCltFtD = ZERO
  do i=1,n 
    do j=1,n 
      temp = ZERO
      do k=1,n
        temp = temp + KG(i,k)*DFCl(j,k) !transposed
      enddo
      KGCltFtD(i,j) = temp
    enddo
  enddo
  KGCltFtDCk = ZERO
  do i=1,n 
    do j=1,n 
      temp = ZERO
      do k=1,n
        temp = temp + KGCltFtD(i,k)*tAk(k,j)
      enddo
      KGCltFtDCk(i,j) = temp
    enddo
  enddo
  tGtKDF = ZERO 
  do i=1,n 
    do j=1,n 
      temp = ZERO
      do k=1,n
        temp = temp + KG(k,i)*DF(k,j) !transpose
      enddo
      tGtKDF(i,j) = temp
    enddo
  enddo
  KGClDFCl = ZERO 
  do i=1,n 
    do j=1,n 
      temp = ZERO
      do k=1,n
        temp = temp + KGCl(i,k)*DFCl(k,j)
      enddo
      KGClDFCl(i,j) = temp
    enddo
  enddo
  KGClDF = ZERO 
  do i=1,n 
    do j=1,n 
      temp = ZERO
      do k=1,n
        temp = temp + KGCl(i,k)*DF(k,j)
      enddo
      KGClDF(i,j) = temp
    enddo
  enddo
  KGCltFtDCl = ZERO 
  do i=1,n 
    do j=1,n 
      temp = ZERO
      do k=1,n
        temp = temp + KGCltFtD(i,k)*tAl(k,j)
      enddo
       KGCltFtDCl(i,j) = temp
    enddo
  enddo
  ClDFCl = ZERO 
  do i=1,n 
    do j=1,n 
      temp = ZERO
      do k=1,n
        temp = temp + tAl(i,k)*DFCl(k,j)
      enddo
      ClDFCl(i,j) = temp
    enddo
  enddo
  ClDF = ZERO 
  do i=1,n 
    do j=1,n 
      temp = ZERO
      do k=1,n
        temp = temp + tAl(i,k)*DF(k,j)
      enddo
      ClDF(i,j) = temp
    enddo
  enddo




  !Symmetrize KGCltFtDCk and CkDFCl
  do i = 1,n
    do j = 1,n
      KGCltFtDCk_s(i,j)=ONEHALF*(KGCltFtDCk(i,j)+KGCltFtDCk(j,i))
      KGCl_s(i,j) = ONEHALF*(KGCl(i,j)+KGCl(j,i))
      CkDFCl_s(i,j) = ONEHALF*(CkDFCl(i,j)+CkDFCl(j,i))
      tGtKDF_s(i,j) = ONEHALF*(tGtKDF(i,j) + tGtKDF(j,i))
      KGClDFCl_s(i,j) = ONEHALF*(KGClDFCl(i,j) + KGClDFCl(j,i))
      KGCltFtDCl_s(i,j) = ONEHALF*(KGCltFtDCl(i,j) + KGCltFtDCl(j,i))
      ClDFCl_s(i,j) =  ONEHALF*(ClDFCl(i,j)+ClDFCl(j,i))
    enddo
  enddo

  !Build modified vectors vl -> CkDFvl, vl -> KGvl
  !wl -> CkDFwl, wl -> KGwl
  !vk -> KGClF'D'vk, wk -> KGClF'D'wk
  CkDFvl = ZERO 
  CkDFwl = ZERO 
  do i=1,n 
    temp = ZERO
    temp1 = ZERO
    do j=1,n 
      temp = temp + CkDF(i,j)*tvl(j)
      temp1 = temp1 + CkDF(i,j)*twl(j)
    enddo
    CkDFvl(i) = temp
    CkDFwl(i) = temp1
  enddo
  KGvl = ZERO
  KGwl = ZERO
  do i=1,n 
    temp = ZERO
    temp1 = ZERO
    do j=1,n 
      temp = temp + KG(i,j)*tvl(j)
      temp1 = temp1 + KG(i,j)*twl(j)
    enddo
    KGvl(i) = temp
    KGwl(i) = temp1
  enddo
  KGCltFtDvk = ZERO
  KGCltFtDwk = ZERO
  KGCltFtDvl = ZERO
  KGCltFtDwl = ZERO
  do i=1,n 
    temp = ZERO
    temp1 = ZERO
    temp2 = ZERO 
    temp3 = ZERO
    do j=1,n 
      temp = temp + KGCltFtD(i,j)*tvk(j)
      temp1 = temp1 + KGCltFtD(i,j)*twk(j)
      temp2 = temp2 + KGCltFtD(i,j)*tvl(j)
      temp3 = temp3 + KGCltFtD(i,j)*twl(j)
    enddo
    KGCltFtDvk(i) = temp
    KGCltFtDwk(i) = temp1
    KGCltFtDvl(i) = temp2 
    KGCltFtDwl(i) = temp3
  enddo
  CltFtDvk = ZERO
  CltFtDwk = ZERO
  CltFtDvl = ZERO
  CltFtDwl = ZERO
  do i=1,n 
    temp = ZERO
    temp1 = ZERO
    temp2 = ZERO 
    temp3 = ZERO 
    do j=1,n 
      temp = temp + DFCl(j,i)*tvk(j) !transpose
      temp1 = temp1 + DFCl(j,i)*twk(j)
      temp2 = temp2 + DFCl(j,i)*tvl(j)
      temp3 = temp3 + DFCl(j,i)*twl(j)
    enddo
    CltFtDvk(i) = temp
    CltFtDwk(i) = temp1
    CltFtDvl(i) = temp2
    CltFtDwl(i) = temp3
  enddo
  DFvl = ZERO 
  DFwl = ZERO
  do i=1,n 
    temp = ZERO 
    temp1 = ZERO
    do j=1,n 
      temp = temp + DF(i,j)*tvl(j)
      temp1 = temp1 + DF(i,j)*twl(j)
    enddo
    DFvl(i) = temp 
    DFwl(i) = temp1 
  enddo
  KGClDFvl = ZERO 
  KGClDFwl = ZERO 
  do i=1,n 
    temp = ZERO 
    temp1 = ZERO
    do j=1,n 
      temp = temp + KGClDF(i,j)*tvl(j)
      temp1 = temp1 + KGClDF(i,j)*twl(j)
    enddo
    KGClDFvl(i) = temp 
    KGClDFwl(i) = temp1 
  enddo
  ClDFvl = ZERO 
  ClDFwl = ZERO 
  do i=1,n 
    temp = ZERO
    temp1 = ZERO
    do j=1,n 
      temp = temp + ClDF(i,j)*tvl(j)
      temp1 = temp1 + ClDF(i,j)*twl(j)
    enddo
    ClDFvl(i) = temp
    ClDFwl(i) = temp1
  enddo

  !Scalars
  vkDFvl = ZERO 
  vkDFwl = ZERO
  wkDFvl = ZERO
  wkDFwl = ZERO
  do i=1,n 
    vkDFvl = vkDFvl + tvk(i)*DFvl(i)
    vkDFwl = vkDFwl + tvk(i)*DFwl(i)
    wkDFvl = wkDFvl + twk(i)*DFvl(i)
    wkDFwl = wkDFwl + twk(i)*DFwl(i)
  enddo
  trDFCl = ZERO 
  do i=1,n 
    trDFCl = trDFCl + DFCl(i,i)
  enddo


  term1 = FOUR*&
    ME_rXr_over_rij_real(rindexI,rindexJ,KGCltFtDCk_s,inv_tAkl,det_tAkl,tvk,tvl,twk,twl)
  term2 = -EIGHT*&
    ME_rXr_rYr_over_rij_real(rindexI,rindexJ,KGCl_s,CkDFCl_s,inv_tAkl,det_tAkl,tvk,tvl,twk,twl)
  term3 = FOUR*&
    ME_rXr_over_rij_real(rindexI,rindexJ,KGCl_s,inv_tAkl,det_tAkl,tvk,CkDFvl,twk,twl)
  term4 = FOUR*&
    ME_rXr_over_rij_real(rindexI,rindexJ,CkDFCl_s,inv_tAkl,det_tAkl,tvk,KGvl,twk,twl)
  term5 = FOUR*&
    ME_rXr_over_rij_real(rindexI,rindexJ,KGCl_s,inv_tAkl,det_tAkl,tvk,tvl,twk,CkDFwl)
  term6 = FOUR*&
    ME_rXr_over_rij_real(rindexI,rindexJ,CkDFCl_s,inv_tAkl,det_tAkl,tvk,tvl,twk,KGwl)
  term7 = -TWO*&
    ME_over_rij_real(rindexI,rindexJ,inv_tAkl,det_tAkl,tvk,CkDFvl,twk,KGwl)
  term8 = -TWO*&
    ME_over_rij_real(rindexI,rindexJ,inv_tAkl,det_tAkl,tvk,KGvl,twk,CkDFwl)
  
  term9 = -TWO*&
    ME_over_rij_real(rindexI,rindexJ,inv_tAkl,det_tAkl,KGCltFtDvk,tvl,twk,twl)
  term10 = FOUR*&
    ME_rXr_over_rij_real(rindexI,rindexJ,KGCl_s,inv_tAkl,det_tAkl,CltFtDvk,tvl,twk,twl)
  term11 = -TWO*vkDFvl*&
    ME_rXr_over_rij_WkWl_real(rindexI,rindexJ,KGCl_s,inv_tAkl,det_tAkl,twk,twl)
  term12 = -TWO*&
    ME_over_rij_real(rindexI,rindexJ,inv_tAkl,det_tAkl,CltFtDvk,KGvl,twk,twl)
  term13 = -TWO*vkDFwl*&
     ME_rXr_over_rij_WkWl_real(rindexI,rindexJ,KGCl_s,inv_tAkl,det_tAkl,twk,tvl)
  term14 = -TWO*&
     ME_over_rij_real(rindexI,rindexJ,inv_tAkl,det_tAkl,CltFtDvk,tvl,twk,KGwl)
  term15 = vkDFvl*&
     ME_over_rij_WkWl_real(rindexI,rindexJ,inv_tAkl,det_tAkl,twk,KGwl)
  term16 = vkDFwl*&
     ME_over_rij_WkWl_real(rindexI,rindexJ,inv_tAkl,det_tAkl,twk,KGvl)

  term17 = -TWO*&
     ME_over_rij_real(rindexI,rindexJ,inv_tAkl,det_tAkl,tvk,tvl,KGCltFtDwk,twl)
  term18 = FOUR*&
     ME_rXr_over_rij_real(rindexI,rindexJ,KGCl_s,inv_tAkl,det_tAkl,tvk,tvl,CltFtDwk,twl)
  term19 = -TWO*wkDFvl*&
     ME_rXr_over_rij_WkWl_real(rindexI,rindexJ,KGCl_s,inv_tAkl,det_tAkl,tvk,twl)
  term20 = -TWO*&
     ME_over_rij_real(rindexI,rindexJ,inv_tAkl,det_tAkl,tvk,KGvl,CltFtDwk,twl)
  term21 = -TWO*WkDFwl*&
     ME_rXr_over_rij_WkWl_real(rindexI,rindexJ,KGCl_s,inv_tAkl,det_tAkl,tvk,tvl)
  term22 = -TWO*&
     ME_over_rij_real(rindexI,rindexJ,inv_tAkl,det_tAkl,tvk,tvl,CltFtDwk,KGwl)
  term23 = wkDFvl*&
     ME_over_rij_WkWl_real(rindexI,rindexJ,inv_tAkl,det_tAkl,tvk,KGwl)
  term24 = WkDFwl*&
     ME_over_rij_WkWl_real(rindexI,rindexJ,inv_tAkl,det_tAkl,tvk,KGvl)

  bigTerm1 = -(term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8 + &
  term9 + term10 + term11 + term12 + term13 + term14 + term15 + term16 + &
  term17 + term18 + term19 + term20 + term21 + term22 + term23 + term24)


  bigTerm2 = -ME_over_rij_dXd(rindexI,rindexJ,tGtKDF_s,tAl,inv_tAkl,det_tAkl,tvk,tvl,twk,twl)

  
  expr1 = 12._wp*trDFCl*&
  ME_rXr_over_rij_real(rindexI,rindexJ,KGCl_s,inv_tAkl,det_tAkl,tvk,tvl,twk,twl)
  expr2 = FOUR*&
  ME_rXr_over_rij_real(rindexI,rindexJ,KGClDFCl_s,inv_tAkl,det_tAkl,tvk,tvl,twk,twl)
  expr3 = -TWO*&
  ME_over_rij_real(rindexI,rindexJ,inv_tAkl,det_tAkl,tvk,KGClDFvl,twk,twl)
  expr4 = -SIX*trDFCl*&
  ME_over_rij_real(rindexI,rindexJ,inv_tAkl,det_tAkl,tvk,KGvl,twk,twl)
  expr5 = -TWO*&
  ME_over_rij_real(rindexI,rindexJ,inv_tAkl,det_tAkl,tvk,tvl,twk,KGClDFwl)
  expr6 = -SIX*trDFCl*&
  ME_over_rij_real(rindexI,rindexJ,inv_tAkl,det_tAkl,tvk,tvl,twk,KGwl)

  expr7 = FOUR*&
  ME_rXr_over_rij_real(rindexI,rindexJ,KGCltFtDCl_s,inv_tAkl,det_tAkl,tvk,tvl,twk,twl)
  expr8 = -EIGHT*&
  ME_rXr_rYr_over_rij_real(rindexI,rindexJ,KGCl_s,ClDFCl_s,inv_tAkl,det_tAkl,tvk,tvl,twk,twl)
  expr9 = FOUR*&
  ME_rXr_over_rij_real(rindexI,rindexJ,KGCl_s,inv_tAkl,det_tAkl,tvk,ClDFvl,twk,twl)
  expr10 = FOUR*&
  ME_rXr_over_rij_real(rindexI,rindexJ,ClDFCl_s,inv_tAkl,det_tAkl,tvk,KGvl,twk,twl)
  expr11 = FOUR*&
  ME_rXr_over_rij_real(rindexI,rindexJ,KGCl_s,inv_tAkl,det_tAkl,tvk,tvl,twk,ClDFwl)
  expr12 = FOUR*&
  ME_rXr_over_rij_real(rindexI,rindexJ,ClDFCl_s,inv_tAkl,det_tAkl,tvk,tvl,twk,KGwl)
  expr13 = -TWO*&
  ME_over_rij_real(rindexI,rindexJ,inv_tAkl,det_tAkl,tvk,ClDFvl,twk,KGwl)
  expr14 = -TWO*&
  ME_over_rij_real(rindexI,rindexJ,inv_tAkl,det_tAkl,tvk,KGvl,twk,ClDFwl)

  expr15 = -TWO*&
  ME_over_rij_real(rindexI,rindexJ,inv_tAkl,det_tAkl,tvk,KGCltFtDvl,twk,twl)
  expr16 = FOUR*&
  ME_rXr_over_rij_real(rindexI,rindexJ,KGCl_s,inv_tAkl,det_tAkl,tvk,CltFtDvl,twk,twl)
  expr17 = ZERO 
  expr18 = -TWO*&
  ME_over_rij_real(rindexI,rindexJ,inv_tAkl,det_tAkl,tvk,CltFtDvl,twk,KGwl)

  expr19 = -TWO*&
  ME_over_rij_real(rindexI,rindexJ,inv_tAkl,det_tAkl,tvk,tvl,twk,KGCltFtDwl)
  expr20 = FOUR*&
  ME_rXr_over_rij_real(rindexI,rindexJ,KGCl_s,inv_tAkl,det_tAkl,tvk,tvl,twk,CltFtDwl)
  expr21 = ZERO 
  expr22 = -TWO*&
  ME_over_rij_real(rindexI,rindexJ,inv_tAkl,det_tAkl,tvk,KGvl,twk,CltFtDwl)

  bigTerm3 = -(expr1 + expr2 + expr3 +expr4 + expr5 + expr6 + expr7 + expr8 + &
  expr9 + expr10 + expr11 + expr12 + expr13 + expr14 + expr15 + expr16 + &
  expr17 + expr18 + expr19 + expr20 + expr21 + expr22)


  ME_KDFG = bigTerm1 + bigTerm2 + bigTerm3


 
end function ME_KDFG

function ME_over_rij_dXd(p,q,X,tAl,inv_tAkl,det_tAkl,tvk,tvl,twk,twl)
    
  real(wp)   ME_over_rij_dXd
  integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles
  !Arguments:
  real(wp)   X(nn,nn),tAl(nn,nn),inv_tAkl(nn,nn),det_tAkl
  integer       p,q
  real(wp)   tvk(nn),tvl(nn),twk(nn),twl(nn)


  !Local variables:
  integer     :: c,s,n,k,i,j
  real(wp) :: ClX(nn,nn), ClXCl(nn,nn), ClXCl_s(nn,nn)
  real(wp) :: commonFactor, gamma, temp, temp1, temp2, temp3
  real(wp) :: ClXvl(nn), ClXwl(nn)
  real(wp) :: trXCl
  real(wp) :: RQ, RVl, RWl, Rtr 

  n = Glob_n
  !Matrices
  ClX = ZERO
  do i=1,n 
    do j=1,n 
      temp = ZERO
      do k=1,n 
        temp = temp + tAl(i,k)*X(k,j)
      enddo
      ClX(i,j) = temp
    enddo
  enddo
  ClXCl = ZERO
  do i=1,n 
    do j=1,n 
      temp = ZERO
      do k=1,n 
        temp = temp + ClX(i,k)*tAl(k,j)
      enddo
      ClXCl(i,j) = temp
    enddo
  enddo
  ClXCl_s = ZERO 
  do i=1,n 
    do j=1,n 
      ClXCl_s(i,j) = ONEHALF*(ClXCl(i,j) + ClXCl(j,i))
    enddo 
  enddo
  
  !Vectors
  ClXvl = ZERO 
  ClXwl = ZERO 
  do i=1,n 
    temp = ZERO 
    temp1 = ZERO 
    do j=1,n 
      temp = temp + ClX(i,j)*tvl(j)
      temp1 = temp1 + ClX(i,j)*twl(j)
    enddo
    ClXvl(i) = temp
    ClXwl(i) = temp1
  enddo

  !Scalars
  trXCl = ZERO 
  do i=1,n 
    do j=1,n 
      trXCl = trXCl + X(i,j)*tAl(j,i)
    enddo
  enddo

  RQ = FOUR*&
  ME_rXr_over_rij_real(p,q,ClXCl_s,inv_tAkl,det_tAkl,tvk,tvl,twk,twl)

  RVl = -FOUR*&
  ME_over_rij_real(p,q,inv_tAkl,det_tAkl,tvk,ClXvl,twk,twl)

  RWl = -FOUR*&
  ME_over_rij_real(p,q,inv_tAkl,det_tAkl,tvk,tvl,twk,ClXwl)

  Rtr = -SIX*trXCl*&
  ME_over_rij_real(p,q,inv_tAkl,det_tAkl,tvk,tvl,twk,twl)

  ME_over_rij_dXd = RQ + RVl + RWl + Rtr

end function ME_over_rij_dXd



function ME_over_rij_WkWl_real(p,q,inv_tAkl,det_tAkl,twk,twl)
  real(wp)  ME_over_rij_WkWl_real

  integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles
  !Arguments:
  integer       p,q
  real(wp)   inv_tAkl(nn,nn),det_tAkl
  real(wp)   twk(nn),twl(nn)

  
  !Local variables:
  integer       n,k,i,j
  real(wp)   XC(nn,nn), CXC(nn,nn)
  real(wp)   W, Wij, jijCwk, jijCwl 
  real(wp)   twkinv_tAkl(nn), inv_tAkltwl(nn)
  real(wp)   commonFactor, gamma, temp, temp1
  real(wp)   I1, I2
  
  
  n = Glob_n
  !Find trA
  twkinv_tAkl = ZERO 
  inv_tAkltwl = ZERO 
  do i=1,n 
    temp = ZERO
    temp1 = ZERO 
    do j=1,n 
      temp = temp + twk(j)*inv_tAkl(j,i)
      temp1 = temp1 + inv_tAkl(i,j)*twl(j)
    enddo
    twkinv_tAkl(i) = temp
    inv_tAkltwl(i) = temp1
  enddo

  W=ZERO
  do i=1,n
    W=W+twkinv_tAkl(i)*twl(i)
  enddo

  if (p == q) then 
    gamma = inv_tAkl(p, p)
    gamma = ONE/sqrt(gamma)
    jijCwk = twkinv_tAkl(p)
    jijCwl = inv_tAkltwl(p)
  else
    gamma = inv_tAkl(p, p) + inv_tAkl(q, q) - inv_tAkl(p, q) -  inv_tAkl(q, p)
    gamma = ONE/sqrt(gamma)
    jijCwk = twkinv_tAkl(p) - twkinv_tAkl(q) 
    jijCwl = inv_tAkltwl(p) - inv_tAkltwl(q) 
  endif
  Wij = jijCwk*jijCwl

  I1 = W*gamma 
  I2 = -Wij*gamma**3/THREE 

  commonFactor = Glob_Piraised3n2/(Glob_SqrtPi*det_tAkl*sqrt(det_tAkl))
  ME_over_rij_WkWl_real = (I1 + I2)*commonFactor

end function ME_over_rij_WkWl_real

function ME_rXr_over_rij_WkWl_real(p,q,X,inv_tAkl,det_tAkl,twk,twl)
  real(wp)  ME_rXr_over_rij_WkWl_real

  integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles
  !Arguments:
  integer       p,q
  real(wp)   X(nn,nn),inv_tAkl(nn,nn),det_tAkl
  real(wp)   twk(nn),twl(nn)

  
  !Local variables:
  integer       n,k,i,j
  real(wp)   XC(nn,nn), CXC(nn,nn)
  real(wp)   trX, trijX, W, Wij, WX, WijX, WXij, WijXij, &
  jijCwk, jijCwl, jijCXCwl, jijCXCwk, jijCXCjij 
  real(wp)   twkinv_tAkl(nn), inv_tAkltwl(nn), CXCwl(nn), CXCwk(nn)
  real(wp)   commonFactor, gamma, temp, temp1
  real(wp)   I1, I2, I3
  
  
  n = Glob_n
  !Find trA
  trX = ZERO
  do i=1,n
    do j=1,n
      trX = trX + inv_tAkl(i,j)*X(j,i) 
    enddo
  enddo


  XC = ZERO
  do i=1,n
    do j=1,n
      temp = ZERO
      do k=1,n
        temp = temp + X(i,k)*inv_tAkl(k,j)
      enddo
      XC(i,j) = temp
    enddo
  enddo

  CXC = ZERO
  do i=1,n
    do j=1,n
      temp = ZERO
      do k=1,n
        temp = temp + inv_tAkl(i,k)*XC(k,j)
      enddo
      CXC(i,j) = temp
    enddo
  enddo


  twkinv_tAkl = ZERO 
  inv_tAkltwl = ZERO 
  do i=1,n 
    temp = ZERO
    temp1 = ZERO 
    do j=1,n 
      temp = temp + twk(j)*inv_tAkl(j,i)
      temp1 = temp1 + inv_tAkl(i,j)*twl(j)
    enddo
    twkinv_tAkl(i) = temp
    inv_tAkltwl(i) = temp1
  enddo

  W=ZERO
  do i=1,n
    W=W+twkinv_tAkl(i)*twl(i)
  enddo

  CXCwl=ZERO
  CXCwk=ZERO
  do i=1,n
    temp = ZERO
    temp1 = ZERO
    do j=1,n 
      temp = temp + CXC(i,j)*twl(j)
      temp1 = temp1 + CXC(i,j)*twk(j)
    enddo
      CXCwl(i) = temp
      CXCwk(i) = temp1
  enddo

  WX = ZERO 
  do i=1,n 
    WX = WX + twk(i)*CXCwl(i)
  enddo


  if (p == q) then 
    gamma = inv_tAkl(p, p)
    gamma = ONE/sqrt(gamma)
    jijCwk = twkinv_tAkl(p)
    jijCwl = inv_tAkltwl(p)
    jijCXCwl = CXCwl(p)
    jijCXCwk = CXCwk(p)
    jijCXCjij = CXC(p, q)
  else
    gamma = inv_tAkl(p, p) + inv_tAkl(q, q) - inv_tAkl(p, q) -  inv_tAkl(q, p)
    gamma = ONE/sqrt(gamma)
    jijCwk = twkinv_tAkl(p) - twkinv_tAkl(q) 
    jijCwl = inv_tAkltwl(p) - inv_tAkltwl(q) 
    jijCXCwl = CXCwl(p) - CXCwl(q)
    jijCXCwk = CXCwk(p) - CXCwk(q)
    jijCXCjij = CXC(p, p) + CXC(q, q) - CXC(p, q) - CXC(q, p)
  endif
  Wij = jijCwk*jijCwl
  WijX = jijCwk*jijCXCwl 
  WXij = jijCXCwk*jijCwl
  WijXij = jijCwk*jijCXCjij*jijCwl
  trijX = jijCXCjij

  I1 = (THREE/TWO*trX*W + WX)*gamma 
  I2 = -(THREE/TWO*(trX*Wij + trijX*W) + WijX + WXij)*gamma**3/THREE 
  I3 = (THREE/TWO*trijX*Wij + WijXij)*gamma**5/FIVE

  commonFactor = Glob_Piraised3n2/(Glob_SqrtPi*det_tAkl*sqrt(det_tAkl))
  ME_rXr_over_rij_WkWl_real = (I1 + I2 + I3)*commonFactor


end function ME_rXr_over_rij_WkWl_real


function ME_over_rij_real(p,q,inv_tAkl,det_tAkl,tvk,tvl,twk,twl)

  real(wp)   ME_over_rij_real
  integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles
  !Arguments:
  integer       p, q
  real(wp)   tAl(nn,nn),tAk(nn,nn),inv_tAkl(nn,nn),det_tAkl
  real(wp)   tvk(nn),tvl(nn),twk(nn),twl(nn)

  !Local variables
  integer  i,j,k,n 
  real(wp)   tvkinv_tAkl(nn), twkinv_tAkl(nn), inv_tAkltvl(nn), inv_tAkltwl(nn)
  real(wp)   V,W,tV,tW,jijCvk,jijCvl,jijCwk,jijCwl,&
  Vij,Wij,tVij,tWij
  real(wp)  commonFactor, gamma, temp, temp1, temp2, temp3
  real(wp)  I1,I2,I3


  n = Glob_n
  tvkinv_tAkl = ZERO 
  twkinv_tAkl = ZERO 
  inv_tAkltvl = ZERO 
  inv_tAkltwl = ZERO 

  do i=1,n 
    temp = ZERO
    temp1 = ZERO 
    temp2 = ZERO 
    temp3 = ZERO
    do j=1,n 
      temp = temp + tvk(j)*inv_tAkl(j,i)
      temp1 = temp1 + twk(j)*inv_tAkl(j,i)
      temp2 = temp2 + inv_tAkl(i,j)*tvl(j)
      temp3 = temp3 + inv_tAkl(i,j)*twl(j)
    enddo
    tvkinv_tAkl(i) = temp
    twkinv_tAkl(i) = temp1 
    inv_tAkltvl(i) = temp2 
    inv_tAkltwl(i) = temp3 
  enddo

  V=ZERO
  W=ZERO
  tV=ZERO
  tW=ZERO
  do i=1,n
    V=V+tvkinv_tAkl(i)*tvl(i)
    W=W+twkinv_tAkl(i)*twl(i)
    tV=tV+tvkinv_tAkl(i)*twl(i)
    tW=tW+twkinv_tAkl(i)*tvl(i)
  enddo

  if (p==q) then
    gamma = inv_tAkl(p,p)
    gamma = ONE/sqrt(gamma)

    jijCvk = tvkinv_tAkl(p)
    jijCvl = inv_tAkltvl(p)
    jijCwk = twkinv_tAkl(p)
    jijCwl = inv_tAkltwl(p)
  else
    gamma = inv_tAkl(p,p) + inv_tAkl(q,q) - inv_tAkl(p,q) - inv_tAkl(q,p)
    gamma = ONE/sqrt(gamma)

    jijCvk = tvkinv_tAkl(p) - tvkinv_tAkl(q)
    jijCvl = inv_tAkltvl(p) - inv_tAkltvl(q)
    jijCwk = twkinv_tAkl(p) - twkinv_tAkl(q)
    jijCwl = inv_tAkltwl(p) - inv_tAkltwl(q)
  endif

  Vij = jijCvk * jijCvl
  Wij = jijCwk * jijCwl
  tVij = jijCvk * jijCwl
  tWij = jijCwk * jijCvl

  I1 = (V*W + tV*tW)*gamma 
  I2 = -(Vij*W + V*Wij + tVij*tW + tV*tWij)*gamma**3/THREE 
  I3 = (Vij*Wij + tVij*tWij)*gamma**5/FIVE 

  commonFactor = ONEHALF*Glob_Piraised3n2/(Glob_SqrtPi*det_tAkl*sqrt(det_tAkl))
  ME_over_rij_real = (I1 + I2 + I3)*commonFactor



end function ME_over_rij_real


function ME_rXr_rYr_over_rij_real(p,q,X,Y,inv_tAkl,det_tAkl,tvk,tvl,twk,twl)

  real(wp)   ME_rXr_rYr_over_rij_real
  integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles
  !Arguments:
  real(wp)   X(nn,nn),Y(nn,nn),inv_tAkl(nn,nn),det_tAkl
  real(wp)   tvk(nn),tvl(nn),twk(nn),twl(nn)
  integer       p,q



  !Local vars
  integer     :: i,j,k,n
  real(wp) ::   tvkinv_tAkl(nn), twkinv_tAkl(nn), inv_tAkltvl(nn), inv_tAkltwl(nn)
  real(wp) :: gamma, temp, temp1, temp2, temp3, &
  temp4, temp5, temp6, temp7, commonFactor
  real(wp) :: V, Vij, VX, VY, VYX, VXY, VijY, VYij, VijX, VXij, &
  VijYX, VYijX, VYXij, VijXY, VXijY, VXYij, VijYij, VijXij, &
  VijYijX, VijYXij, VYijXij,  VijXijY, VijXYij, VXijYij, &
  VijXijYij, VijYijXij
  real(wp) :: W, Wij, WX, WY, WYX, WXY, WijY, WYij, WijX, WXij, &
  WijYX, WYijX, WYXij, WijXY, WXijY, WXYij, WijYij, WijXij, &
  WijYijX, WijYXij, WYijXij, WijXijY, WijXYij, WXijYij, &
  WijXijYij, WijYijXij
  real(wp) :: tV, tVij, tVX, tVY, tVYX, tVXY, tVijY, tVYij, tVijX, tVXij, &
  tVijYX, tVYijX, tVYXij, tVijXY, tVXijY, tVXYij, tVijYij,tVijXij, &
  tVijYijX, tVijYXij, tVYijXij,  tVijXijY, tVijXYij, tVXijYij, &
  tVijXijYij, tVijYijXij
  real(wp) :: tW, tWij, tWX, tWY, tWYX, tWXY, tWijY, tWYij, tWijX, tWXij, &
  tWijYX, tWYijX, tWYXij, tWijXY, tWXijY, tWXYij, tWijYij, tWijXij, &
  tWijYijX, tWijYXij, tWYijXij, tWijXijY, tWijXYij, tWXijYij, &
  tWijXijYij, tWijYijXij
  real(wp) :: jijCvk, jijCvl, jijCwk, jijCwl, &
  jijCXCvk, jijCXCvl, jijCXCwk, jijCXCwl, jijCYCvk, jijCYCvl, jijCYCwk, jijCYCwl, &
  jijCYCXCvl, jijCYCXCwl, jijCYCXCvk, jijCYCXCwk, &
  jijCXCYCvl, jijCXCYCwl, jijCXCYCvk, jijCXCYCwk, &
  jijCYCjij, jijCXCjij, jijCYCXCjij, jijCXCYCjij
  real(wp) :: trX, trY, trYX, trYij, trXij, trijYX, trYijX, trijYijX
 
  real(wp) :: I11, I12, I13, I14, I15, I1, &
  I21, I22, I23, I24, I25, I2, &
  I31, I32, I33, I34, I35, I3, &
  I41, I42, I43, I44, I45, I4, &
  I51, I52, I53, I54, I55, I5

  real(wp) :: CX(nn, nn), XC(nn, nn), CXC(nn, nn), &
  CY(nn, nn), YC(nn, nn), CYC(nn, nn), &
  CXCYC(nn, nn), CYCXC(nn, nn)

  real(wp) :: CXCvl(nn), CXCwl(nn), CYCvl(nn), CYCwl(nn), &
  CXCvk(nn), CXCwk(nn), CYCvk(nn), CYCwk(nn), &
  CXCYCvl(nn), CXCYCwl(nn), CYCXCvl(nn), CYCXCwl(nn), &
  CXCYCvk(nn), CXCYCwk(nn), CYCXCvk(nn), CYCXCwk(nn)

  n = Glob_n

  CX = ZERO
  do i=1,n
    do j=1,n
      CX = CX + inv_tAkl(i,j)*X(j,i) 
    enddo
  enddo

  XC = ZERO
  do i=1,n
    do j=1,n
      temp = ZERO
      do k=1,n
        temp = temp + X(i,k)*inv_tAkl(k,j)
      enddo
      XC(i,j) = temp
    enddo
  enddo

  CXC = ZERO
  do i=1,n
    do j=1,n
      temp = ZERO
      do k=1,n
        temp = temp + inv_tAkl(i,k)*XC(k,j)
      enddo
      CXC(i,j) = temp
    enddo
  enddo

  CY = ZERO
  do i=1,n
    do j=1,n
      CY = CY + inv_tAkl(i,j)*Y(j,i) 
    enddo
  enddo

  YC = ZERO
  do i=1,n
    do j=1,n
      temp = ZERO
      do k=1,n
        temp = temp + Y(i,k)*inv_tAkl(k,j)
      enddo
      YC(i,j) = temp
    enddo
  enddo

  CYC = ZERO
  do i=1,n
    do j=1,n
      temp = ZERO
      do k=1,n
        temp = temp + inv_tAkl(i,k)*YC(k,j)
      enddo
      CYC(i,j) = temp
    enddo
  enddo

  CXCYC = ZERO
  do i=1,n
    do j=1,n
      temp = ZERO
      do k=1,n
        temp = temp + CXC(i,k)*YC(k,j)
      enddo
      CXCYC(i,j) = temp
    enddo
  enddo  

  CYCXC = ZERO
  do i=1,n
    do j=1,n
      temp = ZERO
      do k=1,n
        temp = temp + CYC(i,k)*XC(k,j)
      enddo
     CYCXC(i,j) = temp
    enddo
  enddo  


  trX = ZERO
  do i=1,n
    do k=1,n
       trX = trX + inv_tAkl(i,k)*X(k,i)
    enddo
  enddo
  
  trY = ZERO
  do i=1,n
    do k=1,n
       trY = trY + inv_tAkl(i,k)*Y(k,i)
    enddo
  enddo

  trYX = ZERO
  do i=1,n
    do k=1,n
       trYX = trYX + CYC(i,k)*X(k,i)
    enddo
  enddo

  tvkinv_tAkl = ZERO 
  twkinv_tAkl = ZERO 
  inv_tAkltvl = ZERO 
  inv_tAkltwl = ZERO 
  do i=1,n 
    temp = ZERO
    temp1 = ZERO 
    temp2 = ZERO 
    temp3 = ZERO
    do j=1,n 
      temp = temp + tvk(j)*inv_tAkl(j,i)
      temp1 = temp1 + twk(j)*inv_tAkl(j,i)
      temp2 = temp2 + inv_tAkl(i,j)*tvl(j)
      temp3 = temp3 + inv_tAkl(i,j)*twl(j)
    enddo
    tvkinv_tAkl(i) = temp
    twkinv_tAkl(i) = temp1 
    inv_tAkltvl(i) = temp2 
    inv_tAkltwl(i) = temp3 
  enddo

  V=ZERO
  W=ZERO
  tV=ZERO
  tW=ZERO
  do i=1,n
    V=V+tvkinv_tAkl(i)*tvl(i)
    W=W+twkinv_tAkl(i)*twl(i)
    tV=tV+tvkinv_tAkl(i)*twl(i)
    tW=tW+twkinv_tAkl(i)*tvl(i)
  enddo

  CXCvl=ZERO
  CYCvl=ZERO
  CXCwl=ZERO
  CYCwl=ZERO

  CXCvk=ZERO
  CYCvk=ZERO
  CXCwk=ZERO
  CYCwk=ZERO
  do i=1,n
    temp = ZERO
    temp1 = ZERO
    temp2 = ZERO
    temp3 = ZERO
    temp4 = ZERO
    temp5 = ZERO
    temp6 = ZERO 
    temp7 = ZERO
    do j=1,n 
      temp = temp + CXC(i,j)*tvl(j)
      temp1 = temp1 + CYC(i,j)*tvl(j)
      temp2 = temp2 + CXC(i,j)*twl(j)
      temp3 = temp3 + CYC(i,j)*twl(j)

      temp4 = temp4 + CXC(i,j)*tvk(j)
      temp5 = temp5 + CYC(i,j)*tvk(j)
      temp6 = temp6 + CXC(i,j)*twk(j)
      temp7 = temp7 + CYC(i,j)*twk(j)
    enddo
    CXCvl(i) = temp
    CYCvl(i) = temp1
    CXCwl(i) = temp2
    CYCwl(i) = temp3

    CXCvk(i) = temp4
    CYCvk(i) = temp5 
    CXCwk(i) = temp6
    CYCwk(i) = temp7
  enddo

  CXCYCvl=ZERO
  CYCXCvl=ZERO
  CXCYCwl=ZERO
  CYCXCwl=ZERO

  CXCYCvk=ZERO
  CYCXCvk=ZERO
  CXCYCwk=ZERO
  CYCXCwk=ZERO
  do i=1,n
    temp = ZERO
    temp1 = ZERO
    temp2 = ZERO
    temp3 = ZERO
    temp4 = ZERO
    temp5 = ZERO
    temp6 = ZERO 
    temp7 = ZERO
    do j=1,n 
      temp = temp + CXCYC(i,j)*tvl(j)
      temp1 = temp1 + CYCXC(i,j)*tvl(j)
      temp2 = temp2 + CXCYC(i,j)*twl(j)
      temp3 = temp3 + CYCXC(i,j)*twl(j)

      temp4 = temp4 + CXCYC(i,j)*tvk(j)
      temp5 = temp5 + CYCXC(i,j)*tvk(j)
      temp6 = temp6 + CXCYC(i,j)*twk(j)
      temp7 = temp7 + CYCXC(i,j)*twk(j)
    enddo
    CXCYCvl(i) = temp
    CYCXCvl(i) = temp1
    CXCYCwl(i) = temp2
    CYCXCwl(i) = temp3

    CXCYCvk(i) = temp4
    CYCXCvk(i) = temp5 
    CXCYCwk(i) = temp6
    CYCXCwk(i) = temp7
  enddo


  VX = ZERO 
  VY = ZERO
  VXY = ZERO 
  VYX = ZERO

  tVX = ZERO 
  tVY = ZERO
  tVXY = ZERO 
  tVYX = ZERO

  WX = ZERO 
  WY = ZERO
  WXY = ZERO 
  WYX = ZERO

  tWX = ZERO 
  tWY = ZERO
  tWXY = ZERO 
  tWYX = ZERO
  do i=1,n 
    VX = VX + tvk(i)*CXCvl(i)
    VY = VY + tvk(i)*CYCvl(i)
    VXY = VXY + tvk(i)*CXCYCvl(i)
    VYX = VYX + tvk(i)*CYCXCvl(i)

    tVX = tVX + tvk(i)*CXCwl(i)
    tVY = tVY + tvk(i)*CYCwl(i)
    tVXY = tVXY + tvk(i)*CXCYCwl(i)
    tVYX = tVYX + tvk(i)*CYCXCwl(i)

    WX = WX + twk(i)*CXCwl(i)
    WY = WY + twk(i)*CYCwl(i)
    WXY = WXY + twk(i)*CXCYCwl(i)
    WYX = WYX + twk(i)*CYCXCwl(i)

    tWX = tWX + twk(i)*CXCvl(i)
    tWY = tWY + twk(i)*CYCvl(i)
    tWXY = tWXY + twk(i)*CXCYCvl(i)
    tWYX = tWYX + twk(i)*CYCXCvl(i)

  enddo


   if (p == q) then
    gamma = inv_tAkl(p, p)
    gamma = ONE/sqrt(gamma)
    
    jijCvk = tvkinv_tAkl(p)
    jijCvl = inv_tAkltvl(p)
    jijCXCvl = CXCvl(p)
    jijCXCvk = CXCvk(p)
    jijCYCvl = CYCvl(p)
    jijCYCvk = CYCvk(p)
    jijCYCXCvl = CYCXCvl(p)
    jijCXCYCvl = CXCYCvl(p)
    jijCYCXCvk = CYCXCvk(p)
    jijCXCYCvk = CXCYCvk(p)

    jijCwk = twkinv_tAkl(p)
    jijCwl = inv_tAkltwl(p)
    jijCXCwl = CXCwl(p)
    jijCXCwk = CXCwk(p)
    jijCYCwl = CYCwl(p)
    jijCYCwk = CYCwk(p)
    jijCYCXCwl = CYCXCwl(p)
    jijCXCYCwl = CXCYCwl(p)
    jijCYCXCwk = CYCXCwk(p)
    jijCXCYCwk = CXCYCwk(p)

    jijCYCjij = CYC(p, q)
    jijCXCjij = CXC(p, q)
    jijCYCXCjij = CYCXC(p, q)
    jijCXCYCjij = CXCYC(p, q)

  else
    gamma = inv_tAkl(p, p) + inv_tAkl(q, q) - inv_tAkl(p, q) -  inv_tAkl(q, p)
    gamma = ONE/sqrt(gamma)
    
    jijCvk = tvkinv_tAkl(p) - tvkinv_tAkl(q)
    jijCvl = inv_tAkltvl(p) - inv_tAkltvl(q)
    jijCXCvl = CXCvl(p) - CXCvl(q) 
    jijCXCvk = CXCvk(p) - CXCvk(q)
    jijCYCvl = CYCvl(p) - CYCvl(q)
    jijCYCvk = CYCvk(p) - CYCvk(q)
    jijCYCXCvl = CYCXCvl(p) - CYCXCvl(q)
    jijCXCYCvl = CXCYCvl(p) - CXCYCvl(q)
    jijCYCXCvk = CYCXCvk(p) - CYCXCvk(q)
    jijCXCYCvk = CXCYCvk(p) - CXCYCvk(q)

    jijCwk = twkinv_tAkl(p) - twkinv_tAkl(q)
    jijCwl = inv_tAkltwl(p) - inv_tAkltwl(q) 
    jijCXCwl = CXCwl(p) - CXCwl(q)
    jijCXCwk = CXCwk(p) - CXCwk(q)
    jijCYCwl = CYCwl(p) - CYCwl(q)
    jijCYCwk = CYCwk(p) - CYCwk(q)
    jijCYCXCwl = CYCXCwl(p) - CYCXCwl(q)
    jijCXCYCwl = CXCYCwl(p) - CXCYCwl(q)
    jijCYCXCwk = CYCXCwk(p) - CYCXCwk(q)
    jijCXCYCwk = CXCYCwk(p) - CXCYCwk(q)

    jijCYCjij = CYC(p, p) + CYC(q, q) - CYC(p, q) - CYC(q, p)
    jijCXCjij = CXC(p, p) + CXC(q, q) - CXC(p, q) - CXC(q, p)
    jijCYCXCjij = CYCXC(p, p) + CYCXC(q, q) - CYCXC(p, q) - CYCXC(q, p)
    jijCXCYCjij = CXCYC(p, p) + CXCYC(q, q) - CXCYC(p, q) - CXCYC(q, p)

  endif

  Vij = jijCvk*jijCvl
  VijY = jijCvk*jijCYCvl
  VYij = jijCYCvk*jijCvl
  VijX = jijCvk*jijCXCvl
  VXij = jijCXCvk*jijCvl
  VijYX = jijCvk*jijCYCXCvl
  VYijX = jijCYCvk*jijCXCvl
  VYXij = jijCXCYCvk*jijCvl
  VijXY = jijCvk*jijCXCYCvl
  VXijY = jijCXCvk*jijCYCvl
  VXYij = jijCYCXCvk*jijCvl
  VijYij = jijCvk*jijCYCjij*jijCvl
  VijXij = jijCvk*jijCXCjij*jijCvl
  VijYijX = jijCvk*jijCYCjij*jijCXCvl
  VijYXij = jijCvk*jijCYCXCjij*jijCvl
  VYijXij = jijCYCvk*jijCXCjij*jijCvl
  VijXijY = jijCvk*jijCXCjij*jijCYCvl
  VijXYij = jijCvk*jijCXCYCjij*jijCvl
  VXijYij = jijCXCvk*jijCYCjij*jijCvl
  VijYijXij = jijCvk*jijCYCjij*jijCXCjij*jijCvl 
  VijXijYij = jijCvk*jijCXCjij*jijCYCjij*jijCvl

  tVij = jijCvk*jijCwl
  tVijY = jijCvk*jijCYCwl
  tVYij = jijCYCvk*jijCwl
  tVijX = jijCvk*jijCXCwl
  tVXij = jijCXCvk*jijCwl
  tVijYX = jijCvk*jijCYCXCwl
  tVYijX = jijCYCvk*jijCXCwl
  tVYXij = jijCXCYCvk*jijCwl
  tVijXY = jijCvk*jijCXCYCwl
  tVXijY = jijCXCvk*jijCYCwl
  tVXYij = jijCYCXCvk*jijCwl
  tVijYij = jijCvk*jijCYCjij*jijCwl
  tVijXij = jijCvk*jijCXCjij*jijCwl
  tVijYijX = jijCvk*jijCYCjij*jijCXCwl
  tVijYXij = jijCvk*jijCYCXCjij*jijCwl
  tVYijXij = jijCYCvk*jijCXCjij*jijCwl
  tVijXijY = jijCvk*jijCXCjij*jijCYCwl
  tVijXYij = jijCvk*jijCXCYCjij*jijCwl
  tVXijYij = jijCXCvk*jijCYCjij*jijCwl
  tVijYijXij = jijCvk*jijCYCjij*jijCXCjij*jijCwl 
  tVijXijYij = jijCvk*jijCXCjij*jijCYCjij*jijCwl

  Wij = jijCwk*jijCwl
  WijY = jijCwk*jijCYCwl
  WYij = jijCYCwk*jijCwl
  WijX = jijCwk*jijCXCwl
  WXij = jijCXCwk*jijCwl
  WijYX = jijCwk*jijCYCXCwl
  WYijX = jijCYCwk*jijCXCwl
  WYXij = jijCXCYCwk*jijCwl
  WijXY = jijCwk*jijCXCYCwl
  WXijY = jijCXCwk*jijCYCwl
  WXYij = jijCYCXCwk*jijCwl
  WijYij = jijCwk*jijCYCjij*jijCwl
  WijXij = jijCwk*jijCXCjij*jijCwl
  WijYijX = jijCwk*jijCYCjij*jijCXCwl
  WijYXij = jijCwk*jijCYCXCjij*jijCwl
  WYijXij = jijCYCwk*jijCXCjij*jijCwl
  WijXijY = jijCwk*jijCXCjij*jijCYCwl
  WijXYij = jijCwk*jijCXCYCjij*jijCwl
  WXijYij = jijCXCwk*jijCYCjij*jijCwl
  WijYijXij = jijCwk*jijCYCjij*jijCXCjij*jijCwl 
  WijXijYij = jijCwk*jijCXCjij*jijCYCjij*jijCwl

  tWij = jijCwk*jijCvl
  tWijY = jijCwk*jijCYCvl
  tWYij = jijCYCwk*jijCvl
  tWijX = jijCwk*jijCXCvl
  tWXij = jijCXCwk*jijCvl
  tWijYX = jijCwk*jijCYCXCvl
  tWYijX = jijCYCwk*jijCXCvl
  tWYXij = jijCXCYCwk*jijCvl
  tWijXY = jijCwk*jijCXCYCvl
  tWXijY = jijCXCwk*jijCYCvl
  tWXYij = jijCYCXCwk*jijCvl
  tWijYij = jijCwk*jijCYCjij*jijCvl
  tWijXij = jijCwk*jijCXCjij*jijCvl
  tWijYijX = jijCwk*jijCYCjij*jijCXCvl
  tWijYXij = jijCwk*jijCYCXCjij*jijCvl
  tWYijXij = jijCYCwk*jijCXCjij*jijCvl
  tWijXijY = jijCwk*jijCXCjij*jijCYCvl
  tWijXYij = jijCwk*jijCXCYCjij*jijCvl
  tWXijYij = jijCXCwk*jijCYCjij*jijCvl
  tWijYijXij = jijCwk*jijCYCjij*jijCXCjij*jijCvl 
  tWijXijYij = jijCwk*jijCXCjij*jijCYCjij*jijCvl
  
  trXij = jijCXCjij
  trYij = jijCYCjij
  trijYX = jijCYCXCjij
  trYijX = jijCXCYCjij
  trijYijX = jijCYCjij*jijCXCjij



  I11 = NINE/FOUR*trX*trY*(V*W + tV*tW)
  I12 = THREEHALF*trYX*(V*W + tV*tW)
  I13 = THREEHALF*trX*(VY*W + V*WY + tVY*tW + tV*tWY)
  I14 = THREEHALF*trY*(VX*W + V*WX + tVX*tW + tV*tWX)
  I15 = VYX*W + VXY*W + VX*WY + V*WYX + V*WXY + VY*WX + &
  tVYX*tW + tVXY*tW + tVX*tWY + tV*tWYX + tV*tWXY + tVY*tWX
  I1 = (I11 + I12 + I13 + I14 + I15)*gamma

  I21 = NINE/FOUR*(&
  (trYij*trX + trY*trXij)*(V*W + tV*tW) + trX*trY*(Vij*W + V*Wij + tVij*tW + tV*tWij ))
  I22 = THREEHALF*(&
    (trijYX + trYijX)*(V*W + tV*tW) + trYX*(Vij*W + V*Wij + tVij*tW + tV*tWij));
  I23 = THREEHALF*trX*((VijY + VYij)*W + VY*Wij + V*(WijY + WYij) + Vij*WY + &
  (tVijY + tVYij)*tW + tVY*tWij + tV*(tWijY + tWYij) + tVij*tWY) + &
   THREEHALF*trXij*(VY*W + V*WY + tVY*tW + tV*tWY)
  I24 = THREEHALF*trY*((VijX + VXij)*W + VX*Wij + V*(WijX + WXij) + Vij*WX + &
  (tVijX + tVXij)*tW + tVX*tWij + tV*(tWijX + tWXij) + tVij*tWX) + &
   THREEHALF*trYij*(VX*W + V*WX + tVX*tW + tV*tWX)
  I25 = (VijYX + VYijX + VYXij)*W + VYX*Wij + &
   (VijXY + VXijY + VXYij)*W + VXY*Wij + &
   (VijX + VXij)*WY + VX*(WijY + WYij) + &
   (WijYX + WYijX + WYXij)*V + WYX*Vij + &
   (WijXY + WXijY + WXYij)*V + WXY*Vij + &
   (WijX + WXij)*VY + WX*(VijY + VYij) + &
   (tVijYX + tVYijX + tVYXij)*tW + tVYX*tWij + &
   (tVijXY + tVXijY + tVXYij)*tW + tVXY*tWij + &
   (tVijX + tVXij)*tWY + tVX*(tWijY + tWYij) + &
   (tWijYX + tWYijX + tWYXij)*tV + tWYX*tVij + &
   (tWijXY + tWXijY + tWXYij)*tV + tWXY*tVij + &
   (tWijX + tWXij)*tVY + tWX*(tVijY + tVYij) 
  I2 = -(I21 + I22 + I23 + I24 + I25)*(gamma**3)/THREE


  I31 = NINE/FOUR*(&
    trYij*trXij*(V*W + tV*tW)+(trYij*trX+trY*trXij)*(Vij*W+V*Wij+tVij*tW+tV*tWij)+&
    trX*trY*(Vij*Wij + tVij*tWij))
  I32 = THREEHALF*(&
  trijYijX*(V*W + tV*tW)+(trijYX+trYijX)*(Vij*W+V*Wij+tVij*tW+tV*tWij)+&
  trYX*(Vij*Wij + tVij*tWij))
  I33 = THREEHALF*(&
  trXij*((VijY+VYij)*W+VY*Wij+V*(WijY+WYij)+Vij*WY + &
  (tVijY+tVYij)*tW+tVY*tWij+tV*(tWijY+tWYij)+tVij*tWY)+&
  trX*(VijYij*W+(VijY+VYij)*Wij+V*WijYij+Vij*(WijY+WYij)+&
  tVijYij*tW+(tVijY+tVYij)*tWij+tV*tWijYij+tVij*(tWijY+tWYij)));
  I34 = THREEHALF*(&
  trYij*((VijX+VXij)*W+VX*Wij+V*(WijX+WXij)+Vij*WX + &
  (tVijX+tVXij)*tW+tVX*tWij+tV*(tWijX+tWXij)+tVij*tWX)+&
  trY*(VijXij*W+(VijX+VXij)*Wij+V*WijXij+Vij*(WijX+WXij)+&
  tVijXij*tW+(tVijX+tVXij)*tWij+tV*tWijXij+tVij*(tWijX+tWXij)));
  I35 =  (VijYijX + VijYXij + VYijXij)*W + (VijYX + VYijX + VYXij)*Wij + &
  (VijXijY + VijXYij + VXijYij)*W + (VijXY + VXijY + VXYij)*Wij + &
  VijXij*WY + VX*WijYij + (VijX + VXij)*(WijY + WYij) + &
  (WijYijX + WijYXij + WYijXij)*V + (WijYX + WYijX + WYXij)*Vij + &
  (WijXijY + WijXYij + WXijYij)*V + (WijXY + WXijY + WXYij)*Vij + &
  WijXij*VY + WX*VijYij + (WijX + WXij)*(VijY + VYij) + &
  (tVijYijX + tVijYXij + tVYijXij)*tW + (tVijYX + tVYijX + tVYXij)*tWij + &
  (tVijXijY + tVijXYij + tVXijYij)*tW + (tVijXY + tVXijY + tVXYij)*tWij + &
  tVijXij*tWY + tVX*tWijYij + (tVijX + tVXij)*(tWijY + tWYij) + &
  (tWijYijX + tWijYXij + tWYijXij)*tV + (tWijYX + tWYijX + tWYXij)*tVij + &
  (tWijXijY + tWijXYij + tWXijYij)*tV + (tWijXY + tWXijY + tWXYij)*tVij + &
  tWijXij*tVY + tWX*tVijYij + (tWijX + tWXij)*(tVijY + tVYij)
  I3 = (I31 + I32 + I33 + I34 + I35)*gamma**5/FIVE

  I41 = NINE/FOUR*(&
    trYij*trXij*(Vij*W + V*Wij + tVij*tW + tV*tWij) + &
     (trYij*trX + trY*trXij)*(Vij*Wij + tVij*tWij));
  I42 = THREEHALF*(&
    trijYijX*(Vij*W + V*Wij + tVij*tW + tV*tWij) + &
     (trijYX + trYijX)*(Vij*Wij + tVij*tWij))
  I43 = THREEHALF*(&
    trXij*(VijYij*W + (VijY + VYij)*Wij + Vij*(WijY + WYij) + V*WijYij + &
     tVijYij*tW + (tVijY + tVYij)*tWij + tVij*(tWijY + tWYij) + tV*tWijYij) + &
     trX*(VijYij*Wij + Vij*WijYij + tVijYij*tWij + tVij*tWijYij))
  I44 = THREEHALF*(&
    trYij*(VijXij*W + (VijX + VXij)*Wij + Vij*(WijX + WXij) + V*WijXij + &
    tVijXij*tW + (tVijX + tVXij)*tWij + tVij*(tWijX + tWXij) + tV*tWijXij) + &
    trY*(VijXij*Wij + Vij*WijXij + tVijXij*tWij + tVij*tWijXij))
  I45 = VijYijXij*W + (VijYijX + VijYXij + VYijXij)*Wij + &
   VijXijYij*W + (VijXijY + VijXYij + VXijYij)*Wij + &
   VijXij*(WijY + WYij) + (VijX + VXij)*WijYij + &
   WijYijXij*V + (WijYijX + WijYXij + WYijXij)*Vij + &
   WijXijYij*V + (WijXijY + WijXYij + WXijYij)*Vij + &
   WijXij*(VijY + VYij) + (WijX + WXij)*VijYij + &
   tVijYijXij*tW + (tVijYijX + tVijYXij + tVYijXij)*tWij + &
   tVijXijYij*tW + (tVijXijY + tVijXYij + tVXijYij)*tWij + &
   tVijXij*(tWijY + tWYij) + (tVijX + tVXij)*tWijYij + &
   tWijYijXij*tV + (tWijYijX + tWijYXij + tWYijXij)*tVij + &
   tWijXijYij*tV + (tWijXijY + tWijXYij + tWXijYij)*tVij + &
   tWijXij*(tVijY + tVYij) + (tWijX + tWXij)*tVijYij
   I4 = -(I41 + I42 + I43 + I44 + I45)*gamma**7/SEVEN

   I51 = NINE/FOUR*trYij*trXij*(Vij*Wij + tVij*tWij)
   I52 = THREEHALF*trijYijX*(Vij*Wij + tVij*tWij)
   I53 = THREEHALF*trXij*(VijYij*Wij + Vij*WijYij + tVijYij*tWij + tVij*tWijYij)
   I54 = THREEHALF*trYij*(VijXij*Wij + Vij*WijXij + tVijXij*tWij + tVij*tWijXij)
   I55 = VijYijXij*Wij + VijXijYij*Wij + VijXij*WijYij + &
    WijYijXij*Vij + WijXijYij*Vij + WijXij*VijYij + &
    tVijYijXij*tWij + tVijXijYij*tWij + tVijXij*tWijYij + &
    tWijYijXij*tVij + tWijXijYij*tVij + tWijXij*tVijYij 
   I5 = (I51 + I52 + I53 + I53 + I55)*gamma**9/NINE

  commonFactor = Glob_Piraised3n2/(TWO*Glob_SqrtPi*det_tAkl*sqrt(det_tAkl))
  ME_rXr_rYr_over_rij_real = (I1 + I2 + I3 + I4 + I5)*commonFactor

  end function ME_rXr_rYr_over_rij_real

  function ME_rXr_over_rij_real(p,q,Xs,inv_tAkl,det_tAkl,tvk,tvl,twk,twl)
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
    real(wp)   ME_rXr_over_rij_real

    integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles
    !Arguments:
    integer       p,q
    real(wp)   Xs(nn,nn),inv_tAkl(nn,nn),det_tAkl
    real(wp)   tvk(nn),tvl(nn),twk(nn),twl(nn)

    !Local variables:
    integer       c,s,n,k,i,j
    real(wp)   tvkinv_tAkl(nn), twkinv_tAkl(nn), inv_tAkltvl(nn), inv_tAkltwl(nn)
    real(wp)  commonFactor, gamma, temp, temp1, temp2, temp3

    !Vars for Q-part
    real(wp) :: trXs, Xij, V, Vij, VX, VijX, VXij, VijXij, &
                   W, Wij, WX, WijX, WXij, WijXij, &
                   jijAvk, jijAvl, jijAwk, jijAwl, &
                   jijAXsAvk, jijAXsAvl, jijAXsAwk, jijAXsAwl
    real(wp) :: tV, tVij, tVX, tVijX, tVXij, tVijXij, &
                   tW, tWij, tWX, tWijX, tWXij, tWijXij
    real(wp) :: I11, I12, I13, I1, &
                   I21, I22, I23, I2, &
                   I31, I32, I33, I3, &
                   I41, I42, I43, I4
    real(wp) :: XAl(nn, nn), AXsA(nn, nn), XsA(nn, nn)
    real(wp) :: AXsA_Vl(nn), AXsA_Wl(nn), Vk_AXsA(nn), Wk_AXsA(nn)
    real(wp) :: Qans

    n = Glob_n
    !Build Xs matrix
    !Find trA
    trXs = ZERO
    do i=1,n
      do j=1,n
        trXs = trXs + inv_tAkl(i,j)*Xs(j,i)
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

    tvkinv_tAkl = ZERO
    twkinv_tAkl = ZERO
    inv_tAkltvl = ZERO
    inv_tAkltwl = ZERO

    do i=1,n
      temp = ZERO
      temp1 = ZERO
      temp2 = ZERO
      temp3 = ZERO
      do j=1,n
        temp = temp + tvk(j)*inv_tAkl(j,i)
        temp1 = temp1 + twk(j)*inv_tAkl(j,i)
        temp2 = temp2 + inv_tAkl(i,j)*tvl(j)
        temp3 = temp3 + inv_tAkl(i,j)*twl(j)
      enddo
      tvkinv_tAkl(i) = temp
      twkinv_tAkl(i) = temp1
      inv_tAkltvl(i) = temp2
      inv_tAkltwl(i) = temp3
    enddo

    V=ZERO
    W=ZERO
    tV=ZERO
    tW=ZERO
    do i=1,n
      V=V+tvkinv_tAkl(i)*tvl(i)
      W=W+twkinv_tAkl(i)*twl(i)
      tV=tV+tvkinv_tAkl(i)*twl(i)
      tW=tW+twkinv_tAkl(i)*tvl(i)
    enddo

    do i=1,n
      temp = ZERO
      temp1 = ZERO
      temp2 = ZERO
      temp3 = ZERO
      do j=1,n
        temp = temp + AXsA(i,j)*tvl(j)
        temp1 = temp1 + AXsA(i,j)*twl(j)
        temp2 = temp2 + tvk(j)*AXsA(j,i)
        temp3 = temp3 + twk(j)*AXsA(j,i)
      enddo
      AXsA_Vl(i) = temp
      AXsA_Wl(i) = temp1
      Vk_AXsA(i) = temp2
      Wk_AXsA(i) = temp3
    enddo

    VX = ZERO
    WX = ZERO
    tVX = ZERO
    tWX = ZERO
    do i=1,n
      VX = VX + tvk(i)*AXsA_vl(i)
      WX = WX + twk(i)*AXsA_wl(i)
      tVX = tVX + tvk(i)*AXsA_wl(i)
      tWX = tWX + twk(i)*AXsA_vl(i)
    enddo
  !!! END Q part !!!!

    if (p==q) then
      !Common part
      gamma = inv_tAkl(p,p)
      gamma = ONE/sqrt(gamma)
      !Q part
      Xij = AXsA(p, p)
      jijAvk = tvkinv_tAkl(p)
      jijAvl = inv_tAkltvl(p)
      jijAwk = twkinv_tAkl(p)
      jijAwl = inv_tAkltwl(p)
      jijAXsAvl = AXsA_Vl(p)
      jijAXsAvk = Vk_AXsA(p)
      jijAXsAwl = AXsA_Wl(p)
      jijAXsAwk = Wk_AXsA(p)
    else
      !Common part
      gamma = inv_tAkl(p,p) + inv_tAkl(q,q) - inv_tAkl(p,q) - inv_tAkl(q,p)
      gamma = ONE/sqrt(gamma)
      !Q part
      Xij = AXsA(p,p) + AXsA(q,q) - AXsA(p,q) - AXsA(q,p)
      jijAvk = tvkinv_tAkl(p) - tvkinv_tAkl(q)
      jijAvl = inv_tAkltvl(p) - inv_tAkltvl(q)
      jijAwk = twkinv_tAkl(p) - twkinv_tAkl(q)
      jijAwl = inv_tAkltwl(p) - inv_tAkltwl(q)
      jijAXsAvl = AXsA_Vl(p) - AXsA_Vl(q)
      jijAXsAvk = Vk_AXsA(p) - Vk_AXsA(q)
      jijAXsAwl = AXsA_Wl(p) - AXsA_Wl(q)
      jijAXsAwk = Wk_AXsA(p) - Wk_AXsA(q)
    endif

    commonFactor = ONEHALF*Glob_PiRaised3n2/(sqrt(Glob_Pi)*det_tAkl*sqrt(det_tAkl))

    !Q-part
    Vij = jijAvk * jijAvl
    VijX = jijAvk * jijAXsAvl
    VXij = jijAvl * jijAXsAvk
    VijXij = jijAvk * Xij * jijAvl

    Wij = jijAwk * jijAwl
    WijX = jijAwk * jijAXsAwl
    WXij = jijAwl * jijAXsAwk
    WijXij = jijAwk * Xij * jijAwl

    tVij = jijAvk * jijAwl
    tVijX = jijAvk * jijAXsAwl
    tVXij = jijAwl * jijAXsAvk
    tVijXij = jijAvk * Xij * jijAwl

    tWij = jijAwk * jijAvl
    tWijX = jijAwk * jijAXsAvl
    tWXij = jijAvl * jijAXsAwk
    tWijXij = jijAwk * Xij * jijAvl

    I11 = trXs*V*W + trXs*tV*tW
    I12 = VX*W + tVX*tW
    I13 = V*WX + tV*tWX
    I1 = (THREEHALF*I11 + I12 + I13)*gamma

    I21 = (Xij*V*W + trXs*Vij*W + trXs*V*Wij) + (Xij*tV*tW + trXs*tVij*tW + trXs*tV*tWij)
    I22 = (VijX*W + VXij*W + VX*Wij) + (tVijX*tW + tVXij*tW + tVX*tWij)
    I23 = (V*WijX + V*WXij + Vij*WX) + (tV*tWijX + tV*tWXij + tVij*tWX)
    I2 = -(THREEHALF*I21 + I22 + I23)*gamma**3/THREE

    I31 = (Xij*Vij*W + Xij*V*Wij + trXs*Vij*Wij) + (Xij*tVij*tW + Xij*tV*tWij + trXs*tVij*tWij)
    I32 = (VijXij*W + VijX*Wij + VXij*Wij) + (tVijXij*tW + tVijX*tWij + tVXij*tWij)
    I33 = (V*WijXij + Vij*WijX + Vij*WXij) + (tV*tWijXij + tVij*tWijX + tVij*tWXij)
    I3 = (THREEHALF*I31 + I32 + I33)*gamma**5/FIVE

    I41 = Xij*Vij*Wij + Xij*tVij*tWij
    I42 = VijXij*Wij + tVijXij*tWij
    I43 = Vij*WijXij +  tVij*tWijXij
    I4 = -(THREEHALF*I41 + I42 + I43)*gamma**7/SEVEN

    Qans = (I1 + I2 + I3 + I4)*commonFactor

    ME_rXr_over_rij_real = Qans

  end function ME_rXr_over_rij_real


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

    do i = 1,n
      do j = i+1,n
        t=ONEHALF*(W(j,i)+W(i,j))
        W(j,i) = t
        W(i,j) = t
      end do
    end do

  end subroutine symmetrize_matrix

  function SG_ME_rXr_rYr_over_rij(i,j,X,Y,inv_tAkl,det_tAkl)
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
    real(wp)   t_XYJ,t_YXJ,det_tAkl

    n=Glob_n
!Form Aj=inv_tAkl*ji        j/=i
!     Aj=inv_tAkl*(ji-jj)   j/=i
!Remember that Jii=ji*ji' and Jij=(ji-jj)*(ji-jj)'
    if (i==j) then
      do p=1,n
        Aj(p)=inv_tAkl(p,i)
      enddo
    else
      do p=1,n
        Aj(p)=inv_tAkl(p,i)-inv_tAkl(p,j)
      enddo
    endif

!Compute AjX'=Aj'*X
!    and AjY'=Aj'*Y
    do p=1,n
      temp1=ZERO
      temp2=ZERO
      do q=1,n
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
    do p=1,n
      do q=1,n
        temp1=ZERO
        temp2=ZERO
        do s=1,n
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
    do p=1,n
      temp3=ZERO
      temp4=ZERO
      do q=1,n
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

    do p=1,n
      t_XJ=t_XJ+AjX(p)*Aj(p)
      t_YJ=t_YJ+AjY(p)*Aj(p)
      t_XYJ=t_XYJ+AjX(p)*AYAj(p)
    enddo

!Compute t_YXJ=tr[inv_tAkl*Y*inv_tAkl*X*inv_tAkl*Jij]
    t_YXJ=t_XYJ

    temp1=Glob_PiRaised3n2/(Glob_SqrtPi*det_tAkl**(THREEHALF))
    temp3=1/t_J
    SG_ME_rXr_rYr_over_rij=THREE*temp1*temp3*sqrt(temp3)*(  &
                            THREEHALF*t_J*t_X*t_Y - ONEHALF*(t_Y*t_XJ + t_X*t_YJ) + &
                            t_J*t_XY  - ONETHIRD*(t_XYJ + t_YXJ) + ONEHALF*temp3*t_XJ*t_YJ&
                            )

  end function SG_ME_rXr_rYr_over_rij

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function SG_ME_rXr_rYr_rZr_over_rij(i,j,X,Y,Z,inv_tAkl,det_tAkl)
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
    real(wp)   SG_ME_rXr_rYr_rZr_over_rij
    integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles
!Arguments:
    real(wp)   X(nn,nn),Y(nn,nn),Z(nn,nn),inv_tAkl(nn,nn)
    integer       i,j
!Local variables:
    integer       p,q,s,n
    real(wp)   temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,temp9,temp10,temp11,temp12
    real(wp)   AX(nn,nn),AY(nn,nn),AZ(nn,nn)
    real(wp)   Aj(nn),AjX(nn),AjY(nn),AjZ(nn),AXAj(nn),AYAj(nn),AZAj(nn)
    real(wp)   AZAY(nn,nn),AYAZ(nn,nn),AZAX(nn,nn),AXAZ(nn,nn),AXAY(nn,nn),AYAX(nn,nn)
    real(wp)   AZAYAj(nn),AYAZAj(nn),AZAXAj(nn),AXAZAj(nn),AXAYAj(nn),AYAXAj(nn)
    real(wp)   t_J,t_X,t_Y,t_Z
    real(wp)   t_XJ,t_YJ,t_ZJ
    real(wp)   t_XY,t_ZY,t_ZX,t_YX,t_ZYX,t_YZX
    real(wp)   t_XYJ,t_YXJ,t_ZYJ,t_YZJ,t_XZJ,t_ZXJ
    real(wp)   t_ZYXJ,t_YZXJ,t_YXZJ,t_XYZJ,t_ZXYJ,t_XZYJ
    real(wp)   det_tAkl,term1,term2,term3,term4,term5,term6,term7,term8

    n=Glob_n
!Form Aj=inv_tAkl*ji        j/=i
!     Aj=inv_tAkl*(ji-jj)   j/=i
!Remember that Jii=ji*ji' and Jij=(ji-jj)*(ji-jj)'
    if (i==j) then
      do p=1,n
        Aj(p)=inv_tAkl(p,i)
      enddo
    else
      do p=1,n
        Aj(p)=inv_tAkl(p,i)-inv_tAkl(p,j)
      enddo
    endif

!Compute AjX'=Aj'*X
!    and AjY'=Aj'*Y
    do p=1,n
      temp1=ZERO
      temp2=ZERO
      temp3=ZERO
      do q=1,n
        temp1=temp1+Aj(q)*X(q,p)
        temp2=temp2+Aj(q)*Y(q,p)
        temp3=temp3+Aj(q)*Z(q,p)
      enddo
      AjX(p)=temp1
      AjY(p)=temp2
      AjZ(p)=temp3
    enddo

!Compute AX=inv_tAkl*X  t_X=tr[inv_tAkl*X]
!        AY=inv_tAkl*Y  t_Y=tr[inv_tAkl*Y]
    t_X=ZERO
    t_Y=ZERO
    t_Z=ZERO
    do p=1,n
      do q=1,n
        temp1=ZERO
        temp2=ZERO
        temp3=ZERO
        do s=1,n
          temp1=temp1+inv_tAkl(s,q)*X(p,s)
          temp2=temp2+inv_tAkl(s,q)*Y(p,s)
          temp3=temp3+inv_tAkl(s,q)*Z(p,s)
        enddo
        AX(q,p)=temp1
        AY(q,p)=temp2
        AZ(q,p)=temp3
      enddo
      t_X=t_X+AX(p,p)
      t_Y=t_Y+AY(p,p)
      t_Z=t_Z+AZ(p,p)
    enddo

    do p=1,n
      do q=1,n
        temp1=ZERO
        temp2=ZERO
        temp3=ZERO
        temp4=ZERO
        temp5=ZERO
        temp6=ZERO
        do s=1,n
          temp1=temp1+AZ(s,q)*AY(p,s)
          temp2=temp2+AY(s,q)*AZ(p,s)
          temp3=temp3+AZ(s,q)*AX(p,s)
          temp4=temp4+AX(s,q)*AZ(p,s)
          temp5=temp5+AX(s,q)*AY(p,s)
          temp6=temp6+AY(s,q)*AX(p,s)
        enddo
        AZAY(q,p)=temp1
        AYAZ(q,p)=temp2
        AZAX(q,p)=temp3
        AXAZ(q,p)=temp4
        AXAY(q,p)=temp5
        AYAX(q,p)=temp6
      enddo
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
    t_ZX=ZERO
    t_ZY=ZERO
    t_ZYX=ZERO
    t_YZX=ZERO
    do p=1,n
      temp3=ZERO
      temp4=ZERO
      temp5=ZERO
      temp6=ZERO
      temp7=ZERO
      temp8=ZERO
      temp9=ZERO
      temp10=ZERO
      temp11=ZERO
      do q=1,n
        t_XY=t_XY+AX(p,q)*AY(q,p)
        t_ZX=t_ZX+AZ(p,q)*AX(q,p)
        t_ZY=t_ZY+AZ(p,q)*AY(q,p)
        t_ZYX=t_ZYX+AZAY(p,q)*AX(q,p)
        t_YZX=t_YZX+AYAZ(p,q)*AX(q,p)
        temp3=temp3+AX(p,q)*Aj(q)
        temp4=temp4+AY(p,q)*Aj(q)
        temp5=temp5+AZ(p,q)*Aj(q)
        temp6=temp6+AZAY(p,q)*Aj(q)
        temp7=temp7+AYAZ(p,q)*Aj(q)
        temp8=temp8+AZAX(p,q)*Aj(q)
        temp9=temp9+AXAZ(p,q)*Aj(q)
        temp10=temp10+AXAY(p,q)*Aj(q)
        temp11=temp11+AYAX(p,q)*Aj(q)
      enddo
      AXAj(p)=temp3
      AYAj(p)=temp4
      AZAj(p)=temp5
      AZAYAj(p)=temp6
      AYAZAj(p)=temp7
      AZAXAj(p)=temp8
      AXAZAj(p)=temp9
      AXAYAj(p)=temp10
      AYAXAj(p)=temp11
    enddo

!Compute
!t_XJ=tr[inv_tAkl*X*inv_tAkl*Jij]=AjX'*Aj
!t_YJ=tr[inv_tAkl*Y*inv_tAkl*Jij]=AjY'*Aj
!t_XYJ=tr[inv_tAkl*X*inv_tAkl*Y*inv_tAkl*Jij]=AjX'*AYAj

    t_XJ=ZERO
    t_YJ=ZERO
    t_ZJ=ZERO
    t_XYJ=ZERO
    t_ZYJ=ZERO
    t_ZXJ=ZERO
    t_ZYXJ=ZERO
    t_YZXJ=ZERO
    t_ZXYJ=ZERO
    t_XZYJ=ZERO
    t_YXZJ=ZERO
    t_XYZJ=ZERO
    do p=1,n
      t_XJ=t_XJ+AjX(p)*Aj(p)
      t_YJ=t_YJ+AjY(p)*Aj(p)
      t_ZJ=t_ZJ+AjZ(p)*Aj(p)
      t_XYJ=t_XYJ+AjX(p)*AYAj(p)
      t_ZYJ=t_ZYJ+AjZ(p)*AYAj(p)
      t_ZXJ=t_ZXJ+AjZ(p)*AXAj(p)
      t_ZYXJ=t_ZYXJ+AjZ(p)*AYAXAj(p)
      t_YZXJ=t_YZXJ+AjY(p)*AZAXAj(p)
      t_ZXYJ=t_ZXYJ+AjZ(p)*AXAYAj(p)
      t_XZYJ=t_XZYJ+AjX(p)*AZAYAj(p)
      t_YXZJ=t_YXZJ+AjY(p)*AXAZAj(p)
      t_XYZJ=t_XYZJ+AjX(p)*AYAZAj(p)
    enddo

!Compute t_YXJ=tr[inv_tAkl*Y*inv_tAkl*X*inv_tAkl*Jij]
    t_YXJ=t_XYJ
    t_YZJ=t_ZYJ
    t_XZJ=t_ZXJ

    temp1=Glob_PiRaised3n2/(Glob_SqrtPi*det_tAkl**(THREEHALF))
    temp3=1/t_J
    term1=THREE*temp1*temp3*sqrt(temp3)*THREEHALF*t_Z*(  &
           THREEHALF*t_J*t_X*t_Y - ONEHALF*(t_Y*t_XJ + t_X*t_YJ) + &
           t_J*t_XY  - ONETHIRD*(t_XYJ + t_YXJ) + ONEHALF*temp3*t_XJ*t_YJ&
           )
    term2=-THREEHALF*t_ZY*(3*t_X/sqrt(t_J)-t_XJ/(t_J*sqrt(t_J)))
    term3=THREEHALF*t_Y*(THREEHALF*t_X*t_ZJ/(t_J*sqrt(t_J))-3*t_ZX/sqrt(t_J)&
                         -THREEHALF*t_XJ*t_ZJ/(sqrt(t_J)*t_J*t_J)+(t_ZXJ+t_XZJ)/(t_J*sqrt(t_J)))
    term4= THREEHALF*t_YX*t_ZJ/(t_J*sqrt(t_J))
    term5=-3*(t_ZYX+t_YZX)/sqrt(t_J)-9*t_X*t_YJ*t_ZJ/(4*t_J*t_J*sqrt(t_J))
    term6=THREEHALF*(t_ZX*t_YJ+t_X*(t_ZYJ+t_YZJ))/(t_J*sqrt(t_J))+15*t_XJ*t_YJ*t_ZJ/(4*t_J*t_J*t_J*sqrt(t_J))
    term7=-THREEHALF*(t_YJ*(t_ZXJ+t_XZJ)+t_XJ*(t_ZYJ+t_YZJ)+t_ZJ*(t_YXJ+t_XYJ))/(t_J*t_J*sqrt(t_J))
    term8=(t_ZYXJ+t_YZXJ+t_YXZJ+t_XYZJ+t_ZXYJ+t_XZYJ)/(t_J*sqrt(t_J))

    SG_ME_rXr_rYr_rZr_over_rij=term1-temp1*(term2+term3+term4+term5+term6+term7+term8)
  end function SG_ME_rXr_rYr_rZr_over_rij

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function ME_rXr_over_rij(i,j,X,inv_tAkl,det_tAkl,tvk,tvl,tbk,tbl)
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
    real(wp)   X(nn,nn),inv_tAkl(nn,nn),det_tAkl
    integer       i,j,tvk(nn),tvl(nn),tbk(nn),tbl(nn)

    real(wp)   inv_tAkltvl(nn),inv_tAkltbl(nn),inv_tAkltbk(nn)
    real(wp)   tvkinv_tAkl(nn),tbkinv_tAkl(nn),tvlinv_tAkl(nn)
    real(wp)   tau3,tau33,tau333,tau334,m1,m3,m
!Local variables:
    integer       p,q,n
    real(wp)   Aj(nn),AjX(nn)
    real(wp)   t_J,t_X,t_XJ
    real(wp)   t_JV1,t_XV1,t_JXV1,t_XJV1
    real(wp)   t_JV2,t_XV2,t_JXV2,t_XJV2
    real(wp)   t_JV5,t_XV5,t_JXV5,t_XJV5
    real(wp)   t_JV6,t_XV6,t_JXV6,t_XJV6
    real(wp)   Ajtvl,Ajtbl
    real(wp)    temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8
    real(wp)    temp11,temp22,temp33,mu,mX,mXJ,u,Xmu,muXJ

    n=Glob_n

!Compute inv_tAkltvl = inv_tAkl * tvl
    do p=1,n
      temp1=ZERO
      temp2=ZERO
      temp3=ZERO
      do q=1,n
        temp1=temp1+inv_tAkl(q,p)*tvl(q)
        temp2=temp2+inv_tAkl(q,p)*tbl(q)
        temp3=temp3+inv_tAkl(q,p)*tbk(q)
      enddo
      inv_tAkltvl(p)=temp1
      inv_tAkltbl(p)=temp2
      inv_tAkltbk(p)=temp3
    enddo

!Compute tvkinv_tAkl=tvk'*inv_tAkl
    do p=1,n
      temp1=ZERO
      temp2=ZERO
      temp3=ZERO
      do q=1,n
        temp1=temp1+tvk(q)*inv_tAkl(q,p)
        temp2=temp2+tbk(q)*inv_tAkl(q,p)
        temp3=temp3+tvl(q)*inv_tAkl(q,p)
      enddo
      tvkinv_tAkl(p)=temp1
      tbkinv_tAkl(p)=temp2
      tvlinv_tAkl(p)=temp3
    enddo

!Compute tau3=tvkinv_tAkl*tvl
    tau3=ZERO
    tau33=ZERO
    tau333=ZERO
    tau334=ZERO
    do p=1,n
      tau3=tau3+tvkinv_tAkl(p)*tvl(p)
      tau33=tau33+tbkinv_tAkl(p)*tbl(p)
      tau333=tau333+tvkinv_tAkl(p)*tbl(p)
      tau334=tau334+tbkinv_tAkl(p)*tvl(p)
    enddo
    m1=tau3*tau33
    m3=tau333*tau334
    m=m1+m3  !look formula 40 in document

!Form Aj=inv_tAkl*ji        j/=i
!     Aj=inv_tAkl*(ji-jj)   j/=i
!Remember that Jii=ji*ji' and Jij=(ji-jj)*(ji-jj)'
    if (i==j) then
      do p=1,n
        Aj(p)=inv_tAkl(p,i)
      enddo
    else
      do p=1,n
        Aj(p)=inv_tAkl(p,i)-inv_tAkl(p,j)
      enddo
    endif

!Compute AjX'=Aj'*X
    do p=1,n
      temp1=ZERO
      do q=1,n
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
    Ajtbl=ZERO
    t_XJ=ZERO
    temp1=ZERO
    temp2=ZERO
    temp11=ZERO
    temp22=ZERO
    do p=1,n
      Ajtvl=Ajtvl+Aj(p)*tvl(p)
      Ajtbl=Ajtbl+Aj(p)*tbl(p)
      t_XJ=t_XJ+AjX(p)*Aj(p)
      temp1=temp1+tvk(p)*Aj(p)
      temp2=temp2+AjX(p)*inv_tAkltvl(p)
      temp11=temp11+tbk(p)*Aj(p)
      temp22=temp22+AjX(p)*inv_tAkltbl(p)
    enddo
    t_JV1=temp11*Ajtbl
    t_JXV1=temp11*temp22
    t_JV2=temp1*Ajtvl
    t_JXV2=temp1*temp2
    t_JV5=temp11*Ajtvl
    t_JXV5=temp11*temp2
    t_JV6=temp1*Ajtbl
    t_JXV6=temp1*temp22
!Compute t_X=tr[inv_tAkl*X]
!        t_XV=tr[inv_tAkl*X*inv_tAkl*tvl*tvk']=tvkinv_tAkl'*X*inv_tAkltvl
!        t_XJV=tr[inv_tAkl*X*inv_tAkl*Jij*inv_tAkl*tvl*tvk']=tvkinv_tAkl'*X*Aj*Ajtvl
    t_X=ZERO
    t_XV1=ZERO
    t_XV2=ZERO
    t_XV5=ZERO
    t_XV6=ZERO
    temp2=ZERO
    temp22=ZERO
    do p=1,n
      temp1=ZERO
      temp11=ZERO
      do q=1,n
        t_X=t_X+inv_tAkl(q,p)*X(q,p)
        temp1=temp1+tvkinv_tAkl(q)*X(q,p)
        temp11=temp11+tbkinv_tAkl(q)*X(q,p)
      enddo
      t_XV1=t_XV1+temp11*inv_tAkltbl(p)
      t_XV2=t_XV2+temp1*inv_tAkltvl(p)
      t_XV5=t_XV5+temp11*inv_tAkltvl(p)
      t_XV6=t_XV6+temp1*inv_tAkltbl(p)
      temp2=temp2+temp1*Aj(p)
      temp22=temp22+temp11*Aj(p)
    enddo
    t_XJV1=temp22*Ajtbl
    t_XJV2=temp2*Ajtvl
    t_XJV5=temp22*Ajtvl
    t_XJV6=temp2*Ajtbl

    mu=tau3*t_JV1+tau33*t_JV2+tau333*t_JV5+tau334*t_JV6
    mX=tau3*t_XV1+tau33*t_XV2+tau333*t_XV5+tau334*t_XV6
    u=t_JV1*t_JV2+t_JV5*t_JV6
    mXJ=tau3*(t_XJV1+t_JXV1)+tau33*(t_XJV2+t_JXV2)+tau333*(t_XJV5+t_JXV5)+tau334*(t_XJV6+t_JXV6)
    muXJ=t_JV2*(t_XJV1+t_JXV1)+t_JV1*(t_XJV2+t_JXV2)+t_JV6*(t_XJV5+t_JXV5)+t_JV5*(t_XJV6+t_JXV6)
    Xmu=t_XV2*t_JV1+t_XV1*t_JV2+t_XV6*t_JV5+t_XV5*t_JV6

    temp1=Glob_PiRaised3n2/(TWO*Glob_SqrtPi*det_tAkl**(THREEHALF))
    temp2=+THREEHALF*t_X*(m+ONEFIFTH*u/(t_J*t_J)-ONETHIRD*mu/t_J)/sqrt(t_J)
    temp3=-m*t_XJ/(TWO*t_J*sqrt(t_J))
    temp4=+mX/sqrt(t_J)
    temp5=+mu*t_XJ/(TWO*t_J*t_J*sqrt(t_J))
    temp6=-(mXJ+Xmu)/(THREE*t_J*sqrt(t_J))
    temp7=-u*t_XJ/(TWO*t_J*t_J*t_J*sqrt(t_J))
    temp8=+ONEFIFTH*muXJ/(t_J*t_J*sqrt(t_J))
    ME_rXr_over_rij=temp1*(temp2+temp3+temp4+temp5+temp6+temp7+temp8)

  end function ME_rXr_over_rij

  function ME_over_rij(i,j,inv_tAkl,det_tAkl,tvk,tvl,tbk,tbl)
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
    real(wp)   ME_over_rij
    integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles
!Arguments:
    real(wp)   inv_tAkl(nn,nn),det_tAkl
    integer       i,j,tvk(nn),tvl(nn),tbk(nn),tbl(nn)

    real(wp)   inv_tAkltvl(nn),inv_tAkltbl(nn),inv_tAkltbk(nn)
    real(wp)   tvkinv_tAkl(nn),tbkinv_tAkl(nn),tvlinv_tAkl(nn)
    real(wp)   tau3,tau33,tau333,tau334,m1,m3,m
!Local variables:
    integer       p,q,n
    real(wp)   Aj(nn),AjX(nn),t_JV1,t_JV2,t_JV5,t_JV6,t_J,Ajtvl,Ajtbl
    real(wp)    temp1,temp2,temp3,temp11,temp22,temp33,mu,u

    n=Glob_n

!Compute inv_tAkltvl = inv_tAkl * tvl
    do p=1,n
      temp1=ZERO
      temp2=ZERO
      temp3=ZERO
      do q=1,n
        temp1=temp1+inv_tAkl(q,p)*tvl(q)
        temp2=temp2+inv_tAkl(q,p)*tbl(q)
        temp3=temp3+inv_tAkl(q,p)*tbk(q)
      enddo
      inv_tAkltvl(p)=temp1
      inv_tAkltbl(p)=temp2
      inv_tAkltbk(p)=temp3
    enddo

!Compute tvkinv_tAkl=tvk'*inv_tAkl
    do p=1,n
      temp1=ZERO
      temp2=ZERO
      temp3=ZERO
      do q=1,n
        temp1=temp1+tvk(q)*inv_tAkl(q,p)
        temp2=temp2+tbk(q)*inv_tAkl(q,p)
        temp3=temp3+tvl(q)*inv_tAkl(q,p)
      enddo
      tvkinv_tAkl(p)=temp1
      tbkinv_tAkl(p)=temp2
      tvlinv_tAkl(p)=temp3
    enddo

!Compute tau3=tvkinv_tAkl*tvl
    tau3=ZERO
    tau33=ZERO
    tau333=ZERO
    tau334=ZERO
    do p=1,n
      tau3=tau3+tvkinv_tAkl(p)*tvl(p)
      tau33=tau33+tbkinv_tAkl(p)*tbl(p)
      tau333=tau333+tvkinv_tAkl(p)*tbl(p)
      tau334=tau334+tbkinv_tAkl(p)*tvl(p)
    enddo
    m1=tau3*tau33
    m3=tau333*tau334
    m=m1+m3  !look formula 40 in document

!Form Aj=inv_tAkl*ji        j/=i
!     Aj=inv_tAkl*(ji-jj)   j/=i
!Remember that Jii=ji*ji' and Jij=(ji-jj)*(ji-jj)'
    if (i==j) then
      do p=1,n
        Aj(p)=inv_tAkl(p,i)
      enddo
    else
      do p=1,n
        Aj(p)=inv_tAkl(p,i)-inv_tAkl(p,j)
      enddo
    endif

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
    Ajtbl=ZERO
    temp1=ZERO
    temp11=ZERO
    do p=1,n
      Ajtvl=Ajtvl+Aj(p)*tvl(p)
      Ajtbl=Ajtbl+Aj(p)*tbl(p)
      temp1=temp1+tvk(p)*Aj(p)
      temp11=temp11+tbk(p)*Aj(p)
    enddo
    t_JV1=temp11*Ajtbl
    t_JV2=temp1*Ajtvl
    t_JV5=temp11*Ajtvl
    t_JV6=temp1*Ajtbl

    mu=tau3*t_JV1+tau33*t_JV2+tau333*t_JV5+tau334*t_JV6
    u=t_JV1*t_JV2+t_JV5*t_JV6

    temp1=Glob_PiRaised3n2/(TWO*Glob_SqrtPi*det_tAkl**(THREEHALF))

    ME_over_rij=temp1*(m+ONEFIFTH*u/(t_J*t_J)-ONETHIRD*mu/t_J)/sqrt(t_J)

  end function ME_over_rij

  function ME_d_X_over_rij_d(p,q,X,tAk,tAl,inv_tAkl,det_tAkl,tvk,tvl,twk,twl, &
                             tvkinv_tAkl, twkinv_tAkl, inv_tAkltvl, inv_tAkltwl)

    real(wp)   ME_d_X_over_rij_d
    integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles
    !Arguments:
    real(wp)   X(nn,nn),tAl(nn,nn),tAk(nn,nn),inv_tAkl(nn,nn),det_tAkl, &
      tvkinv_tAkl(nn), twkinv_tAkl(nn), inv_tAkltvl(nn), inv_tAkltwl(nn)
    integer       p,q,tvk(nn),tvl(nn),twk(nn),twl(nn)

    !Local variables:
    integer       c,s,n,k,i,j
    real(wp)   tvkXtAl(nn),tbkXtAl(nn)
    real(wp)   tAkX(nn,nn),tAkXtAl(nn,nn),tAkXtvl(nn),tAkXtbl(nn),XtAl(nn,nn)

    real(wp) :: commonFactor, gamma, temp, temp1, temp2, temp3

    !Vars for Q-part
    real(wp) :: trXs, Xij, V, Vij, VX, VijX, VXij, VijXij, &
                   W, Wij, WX, WijX, WXij, WijXij, &
                   jijAvk, jijAvl, jijAwk, jijAwl, &
                   jijAXsAvk, jijAXsAvl, jijAXsAwk, jijAXsAwl
    real(wp) :: tV, tVij, tVX, tVijX, tVXij, tVijXij, &
                   tW, tWij, tWX, tWijX, tWXij, tWijXij
    real(wp) :: I11, I12, I13, I1, &
                   I21, I22, I23, I2, &
                   I31, I32, I33, I3, &
                   I41, I42, I43, I4
    real(wp) :: XAl(nn, nn), AXsA(nn, nn), Xs(nn, nn), XsA(nn, nn)
    real(wp) :: AXsA_Vl(nn), AXsA_Wl(nn), Vk_AXsA(nn), Wk_AXsA(nn)
    real(wp) :: Qans
    !Vars for RVk part
    real(wp) :: AlA(nn,nn), XAlA(nn,nn), XAlA_Vl(nn), XAlA_Wl(nn), &
                   Vk_XAlA(nn), Wk_XAlA(nn), VkXAlAjij, VkXAlAVl, VkXAlAWl
    real(wp) :: RVk, RVk1, RVk2, RVk3
    !Vars for RWk part
    real(wp) :: WkXAlAVl, WkXAlAWl, WkXAlAjij, RWk, RWk1, RWk2, RWk3
    !Vars for RVl part
    real(wp) :: AkX(nn,nn), AAkX(nn,nn),  AAkX_Vl(nn), AAkX_Wl(nn), &
                   VkAAkXVl, WkAAkXWl, VkAAkXWl, WkAAkXVl, &
                   jijAAkXVl, jijAAkXWl
    real(wp) :: RVl, RVl1, RVl2, RVl3
    !Vars for RWl part
    real(wp) :: RWl, RWl1, RWl2, RWl3
    !Vars for D2 terms
    real(wp) :: X_Vl(nn), X_Wl(nn), &
                   VkXVl, VkXWl, WkXWl, WkXVl, DTwo1, DTwo2, DTwo3, DTwo4, DTwo

    n = Glob_n

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
          temp = temp + tAk(i,k) * XAl(k,j)
        enddo
        Xs(i,j) = temp
      enddo
    enddo

    !Symmetrize MS
    do i = 1,n
      do j = i+1,n
        temp=ONEHALF*(Xs(j,i)+Xs(i,j))
        Xs(j,i) = temp
        Xs(i,j) = temp
      enddo
    enddo

    !Find trA
    trXs = ZERO
    do i=1,n
      do j=1,n
        trXs = trXs + inv_tAkl(i,j)*Xs(j,i)
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

    V=ZERO
    W=ZERO
    tV=ZERO
    tW=ZERO
    do i=1,n
      V=V+tvkinv_tAkl(i)*tvl(i)
      W=W+twkinv_tAkl(i)*twl(i)
      tV=tV+tvkinv_tAkl(i)*twl(i)
      tW=tW+twkinv_tAkl(i)*tvl(i)
    enddo

    do i=1,n
      temp = ZERO
      temp1 = ZERO
      temp2 = ZERO
      temp3 = ZERO
      do j=1,n
        temp = temp + AXsA(i,j)*tvl(j)
        temp1 = temp1 + AXsA(i,j)*twl(j)
        temp2 = temp2 + tvk(j)*AXsA(j,i)
        temp3 = temp3 + twk(j)*AXsA(j,i)
      enddo
      AXsA_Vl(i) = temp
      AXsA_Wl(i) = temp1
      Vk_AXsA(i) = temp2
      Wk_AXsA(i) = temp3
    enddo

    VX = ZERO
    WX = ZERO
    tVX = ZERO
    tWX = ZERO
    do i=1,n
      VX = VX + tvk(i)*AXsA_vl(i)
      WX = WX + twk(i)*AXsA_wl(i)
      tVX = tVX + tvk(i)*AXsA_wl(i)
      tWX = tWX + twk(i)*AXsA_vl(i)
    enddo
  !!! END Q part !!!!

  !!!!! RVk & RWk part of M-matelem !!!!
    AlA = ZERO
    do i=1,n
      do j=1,n
        temp = ZERO
        do k=1,n
          temp = temp + tAl(i,k)*inv_tAkl(k,j)
        enddo
        AlA(i,j) = temp
      enddo
    enddo

    XAlA = ZERO
    do i=1,n
      do j=1,n
        temp = ZERO
        do k=1,n
          temp = temp + X(i,k)*AlA(k,j)
        enddo
        XAlA(i,j) = temp
      enddo
    enddo

    do i=1,n
      temp = ZERO
      temp1 = ZERO
      temp2 = ZERO
      temp3 = ZERO
      do j=1,n
        temp = temp + XAlA(i,j)*tvl(j)
        temp1 = temp1 + XAlA(i,j)*twl(j)
        temp2 = temp2 + tvk(j)*XAlA(j,i)
        temp3 = temp3 + twk(j)*XAlA(j,i)
      enddo
      XAlA_Vl(i) = temp
      XAlA_Wl(i) = temp1
      Vk_XAlA(i) = temp2
      Wk_XAlA(i) = temp3
    enddo

    VkXAlAVl = ZERO
    VkXAlAWl = ZERO
    WkXAlAWl = ZERO
    WkXAlAVl = ZERO
    do i=1,n
      VkXAlAVl = VkXAlAVl + tvk(i)*XAlA_Vl(i)
      VkXAlAWl = VkXAlAWl + tvk(i)*XAlA_Wl(i)
      WkXAlAWl = WkXAlAWl + twk(i)*XAlA_Wl(i)
      WkXAlAVl = WkXAlAVl + twk(i)*XAlA_Vl(i)
    enddo
  !!!!! END RVk & RWk part of M-matelem !!!!

  !!!!! RVl part of M-matelem !!!!
    AkX = ZERO
    do i=1,n
      do j=1,n
        temp = ZERO
        do k=1,n
          temp = temp + tAk(i,k)*X(k,j)
        enddo
        AkX(i,j) = temp
      enddo
    enddo

    AAkX = ZERO
    do i=1,n
      do j=1,n
        temp = ZERO
        do k=1,n
          temp = temp + inv_tAkl(i,k)*AkX(k,j)
        enddo
        AAkX(i,j) = temp
      enddo
    enddo

    do i=1,n
      temp = ZERO
      temp1 = ZERO
      do j=1,n
        temp = temp + AAkX(i,j)*tvl(j)
        temp1 = temp1 + AAkX(i,j)*twl(j)
      enddo
      AAkX_Vl(i) = temp
      AAkX_Wl(i) = temp1
    enddo

    VkAAkXVl = ZERO
    VkAAkXWl = ZERO
    WkAAkXWl = ZERO
    WkAAkXVl = ZERO
    do i=1,n
      VkAAkXVl =  VkAAkXVl + tvk(i)*AAkX_Vl(i)
      VkAAkXWl = VkAAkXWl + tvk(i)*AAkX_Wl(i)
      WkAAkXWl = WkAAkXWl + twk(i)*AAkX_Wl(i)
      WkAAkXVl = WkAAkXVl + twk(i)*AAkX_Vl(i)
    enddo
  !!!!! END RVl part of M-matelem !!!!

  !!!!! D2 terms !!!!
    do i=1,n
      temp = ZERO
      temp1 = ZERO
      do j=1,n
        temp = temp + X(i,j)*tvl(j)
        temp1 = temp1 + X(i,j)*twl(j)
      enddo
      X_Vl(i) = temp
      X_Wl(i) = temp1
    enddo

    VkXVl = ZERO
    VkXWl = ZERO
    WkXWl = ZERO
    WkXVl = ZERO
    do i=1,n
      VkXVl = VkXVl + tvk(i)*X_Vl(i)
      VkXWl = VkXWl + tvk(i)*X_Wl(i)
      WkXWl = WkXWl + tWk(i)*X_Wl(i)
      WkXVl = WkXVl + tWk(i)*X_Vl(i)
    enddo

  !!!!! END D2 terms !!!!

    if (p==q) then
      !Common part
      gamma = inv_tAkl(p,p)
      gamma = ONE/sqrt(gamma)
      !Q part
      Xij = AXsA(p, p)
      jijAvk = tvkinv_tAkl(p)
      jijAvl = inv_tAkltvl(p)
      jijAwk = twkinv_tAkl(p)
      jijAwl = inv_tAkltwl(p)
      jijAXsAvl = AXsA_Vl(p)
      jijAXsAvk = Vk_AXsA(p)
      jijAXsAwl = AXsA_Wl(p)
      jijAXsAwk = Wk_AXsA(p)
      !RVk part
      VkXAlAjij = Vk_XAlA(p)
      !RWk part
      WkXAlAjij = Wk_XAlA(p)
      ! RVl part
      jijAAkXvl = AAkX_Vl(p)
      jijAAkXWl = AAkX_Wl(p)
    else
      !Common part
      gamma = inv_tAkl(p,p) + inv_tAkl(q,q) - inv_tAkl(p,q) - inv_tAkl(q,p)
      gamma = ONE/sqrt(gamma)
      !Q part
      Xij = AXsA(p,p) + AXsA(q,q) - AXsA(p,q) - AXsA(q,p)
      jijAvk = tvkinv_tAkl(p) - tvkinv_tAkl(q)
      jijAvl = inv_tAkltvl(p) - inv_tAkltvl(q)
      jijAwk = twkinv_tAkl(p) - twkinv_tAkl(q)
      jijAwl = inv_tAkltwl(p) - inv_tAkltwl(q)
      jijAXsAvl = AXsA_Vl(p) - AXsA_Vl(q)
      jijAXsAvk = Vk_AXsA(p) - Vk_AXsA(q)
      jijAXsAwl = AXsA_Wl(p) - AXsA_Wl(q)
      jijAXsAwk = Wk_AXsA(p) - Wk_AXsA(q)
      !RVk part
      VkXAlAjij = Vk_XAlA(p) - Vk_XAlA(q)
      !RWk part
      WkXAlAjij = Wk_XAlA(p) - Wk_XAlA(q)
      ! RVl part
      jijAAkXvl = AAkX_Vl(p) - AAkX_Vl(q)
      jijAAkXWl = AAkX_Wl(p) - AAkX_Wl(q)
    endif

    commonFactor = Glob_PiRaised3n2/(sqrt(Glob_Pi)*det_tAkl*sqrt(det_tAkl))

    !Q-part
    Vij = jijAvk * jijAvl
    VijX = jijAvk * jijAXsAvl
    VXij = jijAvl * jijAXsAvk
    VijXij = jijAvk * Xij * jijAvl

    Wij = jijAwk * jijAwl
    WijX = jijAwk * jijAXsAwl
    WXij = jijAwl * jijAXsAwk
    WijXij = jijAwk * Xij * jijAwl

    tVij = jijAvk * jijAwl
    tVijX = jijAvk * jijAXsAwl
    tVXij = jijAwl * jijAXsAvk
    tVijXij = jijAvk * Xij * jijAwl

    tWij = jijAwk * jijAvl
    tWijX = jijAwk * jijAXsAvl
    tWXij = jijAvl * jijAXsAwk
    tWijXij = jijAwk * Xij * jijAvl

    I11 = trXs*V*W + trXs*tV*tW
    I12 = VX*W + tVX*tW
    I13 = V*WX + tV*tWX
    I1 = (THREEHALF*I11 + I12 + I13)*gamma

    I21 = (Xij*V*W + trXs*Vij*W + trXs*V*Wij) + (Xij*tV*tW + trXs*tVij*tW + trXs*tV*tWij)
    I22 = (VijX*W + VXij*W + VX*Wij) + (tVijX*tW + tVXij*tW + tVX*tWij)
    I23 = (V*WijX + V*WXij + Vij*WX) + (tV*tWijX + tV*tWXij + tVij*tWX)
    I2 = -(THREEHALF*I21 + I22 + I23)*gamma**3/THREE

    I31 = (Xij*Vij*W + Xij*V*Wij + trXs*Vij*Wij) + (Xij*tVij*tW + Xij*tV*tWij + trXs*tVij*tWij)
    I32 = (VijXij*W + VijX*Wij + VXij*Wij) + (tVijXij*tW + tVijX*tWij + tVXij*tWij)
    I33 = (V*WijXij + Vij*WijX + Vij*WXij) + (tV*tWijXij + tVij*tWijX + tVij*tWXij)
    I3 = (THREEHALF*I31 + I32 + I33)*gamma**5/FIVE

    I41 = Xij*Vij*Wij + Xij*tVij*tWij
    I42 = VijXij*Wij + tVijXij*tWij
    I43 = Vij*WijXij +  tVij*tWijXij
    I4 = -(THREEHALF*I41 + I42 + I43)*gamma**7/SEVEN

    Qans = (I1 + I2 + I3 + I4)*TWO*commonFactor

    !RVk-part of M-matelem
    RVk1 = gamma*(VkXAlAVl*W +VkXAlAWl*tW)
    RVk2 = -gamma**3/THREE*(VkXAlAjij*jijAvl*W + VkXAlAVl*Wij + &
                            VkXAlAjij*jijAWl*tW + VkXAlAWl*tWij)
    RVk3 = gamma**5/FIVE*(VkXAlAjij*jijAvl*Wij + VkXAlAjij*jijAWl*tWij)
    RVk = -(RVk1 + RVk2 + RVk3)*commonFactor

    !RWk-part of M-matelem
    RWk1 = gamma*(V*WkXAlAWl + tV*WkXAlAVl)
    RWk2 = -gamma**3/THREE*(Vij*WkXAlAWl + V*WkXAlAjij*jijAWl + &
                            tVij*WkXAlAVl + tV*WkXAlAjij*jijAVl)
    RWk3 = gamma**5/FIVE*(Vij*WkXAlAjij*jijAWl + tVij*WkXAlAjij*jijAVl)
    RWk = -(RWk1 + RWk2 + RWk3)*commonFactor

    !RVl-part of M-matelem
    RVl1 = gamma*(VkAAkXVl*W + tV*WkAAkXVl)
    RVl2 = -gamma**3/THREE*(jijAVk*jijAAkXVl*W + VkAAkXVl*Wij + &
                            jijAVk*jijAWl*WkAAkXVl +  tV*jijAWk*jijAAkXVl)
    RVl3 = gamma**5/FIVE*(jijAVk*jijAAkXVl*Wij + &
                          tVij*jijAWk*jijAAkXVl)

    RVl = -(RVl1 + RVl2 + RVl3)*Glob_PiRaised3n2/(sqrt(Glob_Pi)*det_tAkl*sqrt(det_tAkl))

    !RWl-part of M-matelem
    RWl1 = gamma*(V*WkAAkXWl + VkAAkXWl*tW)
    RWl2 = -gamma**3/THREE*(Vij*WkAAkXWl + V*jijAWk*jijAAkXWl + &
                            jijAVk*jijAAkXWl*tW + VkAAkXWl*tWij)
    RWl3 = gamma**5/FIVE*(Vij*jijAWk*jijAAkXWl + jijAVk*jijAAkXWl*tWij)
    RWl = -(RWl1 + RWl2 + RWl3)*Glob_PiRaised3n2/(sqrt(Glob_Pi)*det_tAkl*sqrt(det_tAkl))
    !END RWl-part of M-matelem !

    !D2 part
    DTwo1 = VkXWl*(gamma*tW - ONE/THREE*(gamma**3) * jijAwk * jijAVl)
    DTwo2 = WkXVl*(gamma*tV - ONE/THREE*(gamma**3) * jijAVk * jijAWl)
    DTwo3 = VkXVl*(gamma*W - ONE/THREE*(gamma**3) * jijAWk * jijAWl)
    DTwo4 = WkXWl*(gamma*V - ONE/THREE*(gamma**3) * jijAVk * jijAVl)
    DTwo = (DTwo1 + DTwo2 + DTwo3 + DTwo4)*commonFactor

    ME_d_X_over_rij_d = Qans + RVk + RVl + RWk + RWl + DTwo

  end function ME_d_X_over_rij_d

  function ME_rXr_rYr_over_rij(i,j,X,Y,inv_tAkl,det_tAkl,tvk,tvl,tbk,tbl)
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
    real(wp)   ME_rXr_rYr_over_rij
    integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles
!Arguments:
    real(wp)   X(nn,nn),Y(nn,nn),inv_tAkl(nn,nn),det_tAkl
    integer       i,j,tvk(nn),tvl(nn),tbk(nn),tbl(nn)

    real(wp)   inv_tAkltvl(nn),inv_tAkltbl(nn),inv_tAkltbk(nn)
    real(wp)   tvkinv_tAkl(nn),tbkinv_tAkl(nn),tvlinv_tAkl(nn)
    real(wp)   tau3,tau33,tau333,tau334,m1,m3,m
!Local variables:
    integer       p,q,n,s,k
    real(wp)   Aj(nn),AjX(nn),AjY(nn)
    real(wp)   t_J,t_X,t_XJ,t_Y,t_YJ,t_XYJ
    real(wp)   t_JV1,t_XV1,t_JXV1,t_XJV1,t_YV1,t_JYV1,t_YJV1,t_XYV1,t_YXV1,t_XYJV1,t_YXJV1,t_XJYV1,t_YJXV1,t_JXYV1,t_JYXV1
    real(wp)   t_JV2,t_XV2,t_JXV2,t_XJV2,t_YV2,t_JYV2,t_YJV2,t_XYV2,t_YXV2,t_XYJV2,t_YXJV2,t_XJYV2,t_YJXV2,t_JXYV2,t_JYXV2
    real(wp)   t_JV5,t_XV5,t_JXV5,t_XJV5,t_YV5,t_JYV5,t_YJV5,t_XYV5,t_YXV5,t_XYJV5,t_YXJV5,t_XJYV5,t_YJXV5,t_JXYV5,t_JYXV5
    real(wp)   t_JV6,t_XV6,t_JXV6,t_XJV6,t_YV6,t_JYV6,t_YJV6,t_XYV6,t_YXV6,t_XYJV6,t_YXJV6,t_XJYV6,t_YJXV6,t_JXYV6,t_JYXV6
    real(wp)    Ajtvl,Ajtbl,tvkAj,tbkAj,prod
    real(wp)    temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,temp9,temp44,temp77,temp88,temp99
    real(wp)    temp11,temp22,temp33,mu,mX,mXJ,u,Xmu,muXJ,mY,mYJ,Ymu,muYJ,t_XYJV11,t_XYJV22,t_XYJV55,t_XYJV66,YX,mXYJ
    real(wp)    temp2Y,temp22Y,t_XY,temp55,temp66,AXAj(nn),AYAj(nn),muYX,muXYJ,XYJ,YXJ,t_YXJ,XJYJ,mYX
    real(wp)    tvkinv_tAklX(nn),tbkinv_tAklX(nn),tvkinv_tAklY(nn),tbkinv_tAklY(nn),AY(nn,nn),AX(nn,nn)
    real(wp)    AXinv_tAkltvl(nn),AXinv_tAkltbl(nn),AYinv_tAkltvl(nn),AYinv_tAkltbl(nn)
    real(wp)    term1,term2,term3,term4,term5,term6,term7,term8,term9,term10,term11,term12,term13,temp1Y,temp11Y

    n=Glob_n

!Compute inv_tAkltvl = inv_tAkl * tvl
    do p=1,n
      temp1=ZERO
      temp2=ZERO
      temp3=ZERO
      do q=1,n
        temp1=temp1+inv_tAkl(q,p)*tvl(q)
        temp2=temp2+inv_tAkl(q,p)*tbl(q)
        temp3=temp3+inv_tAkl(q,p)*tbk(q)
      enddo
      inv_tAkltvl(p)=temp1
      inv_tAkltbl(p)=temp2
      inv_tAkltbk(p)=temp3
    enddo

!Compute tvkinv_tAkl=tvk'*inv_tAkl
    do p=1,n
      temp1=ZERO
      temp2=ZERO
      temp3=ZERO
      do q=1,n
        temp1=temp1+tvk(q)*inv_tAkl(q,p)
        temp2=temp2+tbk(q)*inv_tAkl(q,p)
        temp3=temp3+tvl(q)*inv_tAkl(q,p)
      enddo
      tvkinv_tAkl(p)=temp1
      tbkinv_tAkl(p)=temp2
      tvlinv_tAkl(p)=temp3
    enddo

!Compute tau3=tvkinv_tAkl*tvl
    tau3=ZERO
    tau33=ZERO
    tau333=ZERO
    tau334=ZERO
    do p=1,n
      tau3=tau3+tvkinv_tAkl(p)*tvl(p)
      tau33=tau33+tbkinv_tAkl(p)*tbl(p)
      tau333=tau333+tvkinv_tAkl(p)*tbl(p)
      tau334=tau334+tbkinv_tAkl(p)*tvl(p)
    enddo
    m1=tau3*tau33
    m3=tau333*tau334
    m=m1+m3  !look formula 40 in document

!Form Aj=inv_tAkl*ji        j/=i
!     Aj=inv_tAkl*(ji-jj)   j/=i
!Remember that Jii=ji*ji' and Jij=(ji-jj)*(ji-jj)'
    if (i==j) then
      do p=1,n
        Aj(p)=inv_tAkl(p,i)
      enddo
    else
      do p=1,n
        Aj(p)=inv_tAkl(p,i)-inv_tAkl(p,j)
      enddo
    endif

!Compute AjX'=Aj'*X
    do p=1,n
      temp1=ZERO
      temp2=ZERO
      do q=1,n
        temp1=temp1+Aj(q)*X(q,p)
        temp2=temp2+Aj(q)*Y(q,p)
      enddo
      AjX(p)=temp1
      AjY(p)=temp2
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
    Ajtbl=ZERO
    t_XJ=ZERO
    t_YJ=ZERO
    temp1=ZERO
    temp2=ZERO
    temp2Y=ZERO
    temp11=ZERO
    temp22=ZERO
    temp22Y=ZERO
    do p=1,n
      Ajtvl=Ajtvl+Aj(p)*tvl(p)
      Ajtbl=Ajtbl+Aj(p)*tbl(p)
      t_XJ=t_XJ+AjX(p)*Aj(p)
      t_YJ=t_YJ+AjY(p)*Aj(p)
      temp1=temp1+tvk(p)*Aj(p)
      temp2=temp2+AjX(p)*inv_tAkltvl(p)
      temp2Y=temp2Y+AjY(p)*inv_tAkltvl(p)
      temp11=temp11+tbk(p)*Aj(p)
      temp22=temp22+AjX(p)*inv_tAkltbl(p)
      temp22Y=temp22Y+AjY(p)*inv_tAkltbl(p)
    enddo
    t_JV1=temp11*Ajtbl
    t_JXV1=temp11*temp22
    t_JYV1=temp11*temp22Y
    t_JV2=temp1*Ajtvl
    t_JXV2=temp1*temp2
    t_JYV2=temp1*temp2Y
    t_JV5=temp11*Ajtvl
    t_JXV5=temp11*temp2
    t_JYV5=temp11*temp2Y
    t_JV6=temp1*Ajtbl
    t_JXV6=temp1*temp22
    t_JYV6=temp1*temp22Y
!Compute t_X=tr[inv_tAkl*X]
!        t_XV=tr[inv_tAkl*X*inv_tAkl*tvl*tvk']=tvkinv_tAkl'*X*inv_tAkltvl
!        t_XJV=tr[inv_tAkl*X*inv_tAkl*Jij*inv_tAkl*tvl*tvk']=tvkinv_tAkl'*X*Aj*Ajtvl
    do p=1,n
      do q=1,n
        temp1=ZERO
        temp2=ZERO
        do s=1,n
          temp1=temp1+inv_tAkl(s,q)*X(p,s)
          temp2=temp2+inv_tAkl(s,q)*Y(p,s)
        enddo
        AX(q,p)=temp1
        AY(q,p)=temp2
      enddo
    enddo

    t_X=ZERO
    t_Y=ZERO
    t_XV1=ZERO
    t_XV2=ZERO
    t_XV5=ZERO
    t_XV6=ZERO
    t_YV1=ZERO
    t_YV2=ZERO
    t_YV5=ZERO
    t_YV6=ZERO
    temp2=ZERO
    temp22=ZERO
    temp2Y=ZERO
    temp22Y=ZERO
    do p=1,n
      temp1=ZERO
      temp11=ZERO
      temp1Y=ZERO
      temp11Y=ZERO
      do q=1,n
        t_X=t_X+inv_tAkl(q,p)*X(q,p)
        t_Y=t_Y+inv_tAkl(q,p)*Y(q,p)
        temp1=temp1+tvkinv_tAkl(q)*X(q,p)
        temp11=temp11+tbkinv_tAkl(q)*X(q,p)
        temp1Y=temp1Y+tvkinv_tAkl(q)*Y(q,p)
        temp11Y=temp11Y+tbkinv_tAkl(q)*Y(q,p)
      enddo
      t_XV1=t_XV1+temp11*inv_tAkltbl(p)
      t_XV2=t_XV2+temp1*inv_tAkltvl(p)
      t_XV5=t_XV5+temp11*inv_tAkltvl(p)
      t_XV6=t_XV6+temp1*inv_tAkltbl(p)
      t_YV1=t_YV1+temp11Y*inv_tAkltbl(p)
      t_YV2=t_YV2+temp1Y*inv_tAkltvl(p)
      t_YV5=t_YV5+temp11Y*inv_tAkltvl(p)
      t_YV6=t_YV6+temp1Y*inv_tAkltbl(p)
      temp2=temp2+temp1*Aj(p)
      temp22=temp22+temp11*Aj(p)
      temp2Y=temp2Y+temp1Y*Aj(p)
      temp22Y=temp22Y+temp11Y*Aj(p)
    enddo
    t_XJV1=temp22*Ajtbl
    t_XJV2=temp2*Ajtvl
    t_XJV5=temp22*Ajtvl
    t_XJV6=temp2*Ajtbl
    t_YJV1=temp22Y*Ajtbl
    t_YJV2=temp2Y*Ajtvl
    t_YJV5=temp22Y*Ajtvl
    t_YJV6=temp2Y*Ajtbl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    t_XY=ZERO
    do p=1,n
      temp1=ZERO
      temp2=ZERO
      temp11=ZERO
      temp22=ZERO
      temp3=ZERO
      temp4=ZERO
      temp5=ZERO
      temp6=ZERO
      temp55=ZERO
      temp66=ZERO
      do q=1,n
        t_XY=t_XY+AX(p,q)*AY(q,p)
        temp1=temp1+AX(p,q)*inv_tAkltvl(q)
        temp2=temp2+AY(p,q)*inv_tAkltvl(q)
        temp11=temp11+AX(p,q)*inv_tAkltbl(q)
        temp22=temp22+AY(p,q)*inv_tAkltbl(q)
        temp3=temp3+AX(p,q)*Aj(q)
        temp4=temp4+AY(p,q)*Aj(q)
        temp5=temp5+tvkinv_tAkl(q)*X(q,p)
        temp6=temp6+tvkinv_tAkl(q)*Y(q,p)
        temp55=temp55+tbkinv_tAkl(q)*X(q,p)
        temp66=temp66+tbkinv_tAkl(q)*Y(q,p)
      enddo
      AXinv_tAkltvl(p)=temp1
      AYinv_tAkltvl(p)=temp2
      AXinv_tAkltbl(p)=temp11
      AYinv_tAkltbl(p)=temp22
      AXAj(p)=temp3
      AYAj(p)=temp4
      tvkinv_tAklX(p)=temp5
      tvkinv_tAklY(p)=temp6
      tbkinv_tAklX(p)=temp55
      tbkinv_tAklY(p)=temp66
    enddo

    tvkAj=ZERO
    tbkAj=ZERO
    temp1=ZERO
    temp11=ZERO
    temp2=ZERO
    temp22=ZERO
    temp3=ZERO
    temp33=ZERO
    t_XYJ=ZERO
    t_YXJ=ZERO
    t_XYV1=ZERO
    t_YXV1=ZERO
    t_XYV2=ZERO
    t_YXV2=ZERO
    t_XYV5=ZERO
    t_YXV5=ZERO
    t_XYV6=ZERO
    t_YXV6=ZERO
    temp4=ZERO
    temp5=ZERO
    temp44=ZERO
    temp55=ZERO
    temp6=ZERO
    temp7=ZERO
    temp8=ZERO
    temp9=ZERO
    temp66=ZERO
    temp77=ZERO
    temp88=ZERO
    temp99=ZERO
    do p=1,n
      tvkAj=tvkAj+tvk(p)*Aj(p)
      tbkAj=tbkAj+tbk(p)*Aj(p)
!  temp1=temp1+tvk(p)*Aj(p)
!  temp11=temp11+tbk(p)*Aj(p)
      temp2=temp2+AjX(p)*inv_tAkltvl(p)
      temp22=temp22+AjX(p)*inv_tAkltbl(p)
      temp3=temp3+AjY(p)*inv_tAkltvl(p)
      temp33=temp33+AjY(p)*inv_tAkltbl(p)
      t_XYJ=t_XYJ+AjX(p)*AYAj(p)
      t_YXJ=t_YXJ+AjY(p)*AXAj(p)
      t_XYV1=t_XYV1+tbkinv_tAklX(p)*AYinv_tAkltbl(p)
      t_XYV2=t_XYV2+tvkinv_tAklX(p)*AYinv_tAkltvl(p)
      t_XYV5=t_XYV5+tbkinv_tAklX(p)*AYinv_tAkltvl(p)
      t_XYV6=t_XYV6+tvkinv_tAklX(p)*AYinv_tAkltbl(p)
      t_YXV1=t_YXV1+tbkinv_tAklY(p)*AXinv_tAkltbl(p)
      t_YXV2=t_YXV2+tvkinv_tAklY(p)*AXinv_tAkltvl(p)
      t_YXV5=t_YXV5+tbkinv_tAklY(p)*AXinv_tAkltvl(p)
      t_YXV6=t_YXV6+tvkinv_tAklY(p)*AXinv_tAkltbl(p)
      temp4=temp4+tvkinv_tAklX(p)*Aj(p)
      temp5=temp5+tvkinv_tAklY(p)*Aj(p)
      temp44=temp44+tbkinv_tAklX(p)*Aj(p)
      temp55=temp55+tbkinv_tAklY(p)*Aj(p)
      temp6=temp6+tvkinv_tAklX(p)*AYAj(p)
      temp7=temp7+tvkinv_tAklY(p)*AXAj(p)
      temp8=temp8+AjX(p)*AYinv_tAkltvl(p)
      temp9=temp9+AjY(p)*AXinv_tAkltvl(p)
      temp66=temp66+tbkinv_tAklX(p)*AYAj(p)
      temp77=temp77+tbkinv_tAklY(p)*AXAj(p)
      temp88=temp88+AjX(p)*AYinv_tAkltbl(p)
      temp99=temp99+AjY(p)*AXinv_tAkltbl(p)
    enddo
    t_XYJV1=temp66*Ajtbl
    t_XYJV2=temp6*Ajtvl
    t_XYJV5=temp66*Ajtvl
    t_XYJV6=temp6*Ajtbl
    t_YXJV1=temp77*Ajtbl
    t_YXJV2=temp7*Ajtvl
    t_YXJV5=temp77*Ajtvl
    t_YXJV6=temp7*Ajtbl
    t_XJYV1=temp44*temp33
    t_XJYV2=temp4*temp3
    t_XJYV5=temp44*temp3
    t_XJYV6=temp4*temp33
    t_YJXV1=temp55*temp22
    t_YJXV2=temp5*temp2
    t_YJXV5=temp55*temp2
    t_YJXV6=temp5*temp22
    t_JXYV1=tbkAj*temp88
    t_JXYV2=tvkAj*temp8
    t_JXYV5=tbkAj*temp8
    t_JXYV6=tvkAj*temp88
    t_JYXV1=tbkAj*temp99
    t_JYXV2=tvkAj*temp9
    t_JYXV5=tbkAj*temp9
    t_JYXV6=tvkAj*temp99
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    t_XYJV11=(t_XYJV1+t_YXJV1+t_XJYV1+t_YJXV1+t_JXYV1+t_JYXV1)
    t_XYJV22=(t_XYJV2+t_YXJV2+t_XJYV2+t_YJXV2+t_JXYV2+t_JYXV2)
    t_XYJV55=(t_XYJV5+t_YXJV5+t_XJYV5+t_YJXV5+t_JXYV5+t_JYXV5)
    t_XYJV66=(t_XYJV6+t_YXJV6+t_XJYV6+t_YJXV6+t_JXYV6+t_JYXV6)
!!!!!!!!!!!!!!!!!!!!!!
    mu=tau3*t_JV1+tau33*t_JV2+tau333*t_JV5+tau334*t_JV6
    mX=tau3*t_XV1+tau33*t_XV2+tau333*t_XV5+tau334*t_XV6
    mY=tau3*t_YV1+tau33*t_YV2+tau333*t_YV5+tau334*t_YV6
    mYX=tau3*t_YXV1+tau33*t_YXV2+tau333*t_YXV5+tau334*t_YXV6
    u=t_JV1*t_JV2+t_JV5*t_JV6
    mXJ=tau3*(t_XJV1+t_JXV1)+tau33*(t_XJV2+t_JXV2)+tau333*(t_XJV5+t_JXV5)+tau334*(t_XJV6+t_JXV6)
    mYJ=tau3*(t_YJV1+t_JYV1)+tau33*(t_YJV2+t_JYV2)+tau333*(t_YJV5+t_JYV5)+tau334*(t_YJV6+t_JYV6)
    muXJ=t_JV2*(t_XJV1+t_JXV1)+t_JV1*(t_XJV2+t_JXV2)+t_JV6*(t_XJV5+t_JXV5)+t_JV5*(t_XJV6+t_JXV6)
    muYJ=t_JV2*(t_YJV1+t_JYV1)+t_JV1*(t_YJV2+t_JYV2)+t_JV6*(t_YJV5+t_JYV5)+t_JV5*(t_YJV6+t_JYV6)
    Xmu=t_XV2*t_JV1+t_XV1*t_JV2+t_XV6*t_JV5+t_XV5*t_JV6
    Ymu=t_YV2*t_JV1+t_YV1*t_JV2+t_YV6*t_JV5+t_YV5*t_JV6

    mXYJ=tau3*t_XYJV11+tau33*t_XYJV22+tau333*t_XYJV55+tau334*t_XYJV66
    YX=t_YV2*t_XV1+t_YV1*t_XV2+t_YV6*t_XV5+t_YV5*t_XV6
    muYX=t_JV1*t_YXV2+t_JV2*t_YXV1+t_JV6*t_YXV5+t_JV5*t_YXV6
    muXYJ=t_JV1*t_XYJV22+t_JV2*t_XYJV11+t_JV6*t_XYJV55+t_JV5*t_XYJV66
    XYJ=t_XV1*t_YJV2+t_XV2*t_YJV1+t_XV5*t_YJV6+t_XV6*t_YJV5
    YXJ=t_YV1*t_XJV2+t_YV2*t_XJV1+t_YV5*t_XJV6+t_YV6*t_XJV5
XJYJ=(t_XJV1+t_JXV1)*(t_YJV2+t_JYV2)+(t_XJV2+t_JXV2)*(t_YJV1+t_JYV1)+(t_XJV5+t_JXV5)*(t_YJV6+t_JYV6)+(t_XJV6+t_JXV6)*(t_YJV5+t_JYV5)

    temp1=Glob_PiRaised3n2/(TWO*Glob_SqrtPi*det_tAkl**(THREEHALF))
    temp2=+THREEHALF*t_X*(m+ONEFIFTH*u/(t_J*t_J)-ONETHIRD*mu/t_J)/sqrt(t_J)
    temp3=-m*t_XJ/(TWO*t_J*sqrt(t_J))
    temp4=+mX/sqrt(t_J)
    temp5=+mu*t_XJ/(TWO*t_J*t_J*sqrt(t_J))
    temp6=-(mXJ+Xmu)/(THREE*t_J*sqrt(t_J))
    temp7=-u*t_XJ/(TWO*t_J*t_J*t_J*sqrt(t_J))
    temp8=+ONEFIFTH*muXJ/(t_J*t_J*sqrt(t_J))

    term1=THREEHALF*temp1*t_Y*(temp2+temp3+temp4+temp5+temp6+temp7+temp8)
    term2=-THREEHALF*t_XY*(m-ONETHIRD*mu/t_J+ONEFIFTH*u/(t_J*t_J))/sqrt(t_J)
    term3=+THREEHALF*t_X*(t_YJ*m/(TWO*t_J)-mY-t_YJ*mu/(TWO*t_J*t_J) &
                          +ONETHIRD*(Ymu+mYJ)/t_J+t_YJ*u/(TWO*t_J**3)-ONEFIFTH*muYJ/(t_J*t_J))/sqrt(t_J)
    term4=-3*m*t_XJ*t_YJ/(4*t_J*t_J*sqrt(t_J))
    term5=+(mY*t_XJ+m*(t_YXJ+t_XYJ)+mX*t_YJ)/(2*t_J*sqrt(t_J))
    term6=-(YX+mYX)/sqrt(t_J)
    term7=+5*mu*t_XJ*t_YJ/(4*sqrt(t_J)*t_J**3)
    term8=-(mu*(t_YXJ+t_XYJ)+t_XJ*(Ymu+mYJ)+t_YJ*(Xmu+mXJ))/(2*sqrt(t_J)*t_J**2)
    term9=+ONETHIRD*(muYX+XYJ+YXJ+mXYJ)/(sqrt(t_J)*t_J)
    term10=-7*t_YJ*t_XJ*u/(4*sqrt(t_J)*t_J**4)
    term11=+(u*(t_YXJ+t_XYJ)+t_XJ*muYJ+t_YJ*muXJ)/(2*sqrt(t_J)*t_J**3)
    term12=-ONEFIFTH*(muXYJ+XJYJ)/(sqrt(t_J)*t_J**2)

    term13=temp1*(term2+term3+term4+term5+term6+term7+term8+term9+term10+term11+term12)

    ME_rXr_rYr_over_rij=term1-term13

  end function ME_rXr_rYr_over_rij

  function rPr_rQr(P,Q,tvk,tbk,inv_tAkl,det_tAkl,tau3,tau33,tau333,tau334,tvkinv_tAkl,tbkinv_tAkl,inv_tAkltvl,inv_tAkltbl)
    real(wp)   rPr_rQr
!arguments
    integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles
    integer       tvk(nn),tbk(nn)
    real(wp)   P(nn,nn),Q(nn,nn),inv_tAkl(nn,nn),tau3,tau33,tau333,tau334,det_tAkl
    real(wp)   inv_tAkltvl(nn),inv_tAkltbl(nn),tvkinv_tAkl(nn),tbkinv_tAkl(nn)

!declaration
    integer       i,j,n,k
    real(wp)   temp1,temp2,temp3,temp4,P1,P2,P5,P6,Q1,Q2,Q5,Q6,trP,trQ,trQP
    real(wp)   PQ1,PQ2,PQ5,PQ6,QP1,QP2,QP5,QP6
    real(wp)   inv_tAklP(nn,nn),inv_tAklQ(nn,nn),inv_tAklQP(nn,nn)
    real(wp)   tvkinv_tAklP(nn),tbkinv_tAklP(nn),tvkinv_tAklQ(nn),tbkinv_tAklQ(nn)
    real(wp)   Pgamma,Qgamma,PQgamma,QPgamma
    n=Glob_n

!Doing multiplication inv_tAkltAlM=inv_tAkltAl*M
    do i=1,n
      do j=1,n
        temp1=ZERO
        temp2=ZERO
        do k=1,n
          temp1=temp1+inv_tAkl(j,k)*P(k,i)
          temp2=temp2+inv_tAkl(j,k)*Q(k,i)
        enddo
        inv_tAklP(j,i)=temp1
        inv_tAklQ(j,i)=temp2
      enddo
    enddo

    do i=1,n
      do j=1,n
        temp1=ZERO
        do k=1,n
          temp1=temp1+inv_tAklQ(j,k)*inv_tAklP(k,i)
        enddo
        inv_tAklQP(j,i)=temp1
      enddo
    enddo

    trQP=ZERO
    trP=ZERO
    trQ=ZERO
    do i=1,n
      trQP=trQP+inv_tAklQP(i,i)
    enddo

    do i=1,n
      do j=1,n
        trP=trP+inv_tAkl(i,j)*P(j,i)
        trQ=trQ+inv_tAkl(i,j)*Q(j,i)
      enddo
    enddo

    Q1=ZERO
    Q2=ZERO
    Q5=ZERO
    Q6=ZERO
    P1=ZERO
    P2=ZERO
    P5=ZERO
    P6=ZERO
    do i=1,n
      temp1=ZERO
      temp2=ZERO
      temp3=ZERO
      temp4=ZERO
      do j=1,n
        temp1=temp1+tbkinv_tAkl(j)*P(j,i)
        temp2=temp2+tvkinv_tAkl(j)*P(j,i)
        temp3=temp3+tbkinv_tAkl(j)*Q(j,i)
        temp4=temp4+tvkinv_tAkl(j)*Q(j,i)
      enddo
      P1=P1+temp1*inv_tAkltbl(i)
      P2=P2+temp2*inv_tAkltvl(i)
      P5=P5+temp1*inv_tAkltvl(i)
      P6=P6+temp2*inv_tAkltbl(i)
      Q1=Q1+temp3*inv_tAkltbl(i)
      Q2=Q2+temp4*inv_tAkltvl(i)
      Q5=Q5+temp3*inv_tAkltvl(i)
      Q6=Q6+temp4*inv_tAkltbl(i)
    enddo

    do i=1,n
      temp1=ZERO
      temp2=ZERO
      temp3=ZERO
      temp4=ZERO
      do j=1,n
        temp1=temp1+tvk(j)*inv_tAklP(j,i)
        temp2=temp2+tbk(j)*inv_tAklP(j,i)
        temp3=temp3+tvk(j)*inv_tAklQ(j,i)
        temp4=temp4+tbk(j)*inv_tAklQ(j,i)
      enddo
      tvkinv_tAklP(i)=temp1
      tbkinv_tAklP(i)=temp2
      tvkinv_tAklQ(i)=temp3
      tbkinv_tAklQ(i)=temp4
    enddo

    PQ1=ZERO
    PQ2=ZERO
    PQ5=ZERO
    PQ6=ZERO
    QP1=ZERO
    QP2=ZERO
    QP5=ZERO
    QP6=ZERO
    do i=1,n
      temp1=ZERO
      temp2=ZERO
      temp3=ZERO
      temp4=ZERO
      do j=1,n
        temp1=temp1+tbkinv_tAklP(j)*inv_tAklQ(j,i)
        temp2=temp2+tvkinv_tAklP(j)*inv_tAklQ(j,i)
        temp3=temp3+tbkinv_tAklQ(j)*inv_tAklP(j,i)
        temp4=temp4+tvkinv_tAklQ(j)*inv_tAklP(j,i)
      enddo
      PQ1=PQ1+temp1*inv_tAkltbl(i)
      PQ2=PQ2+temp2*inv_tAkltvl(i)
      PQ5=PQ5+temp1*inv_tAkltvl(i)
      PQ6=PQ6+temp2*inv_tAkltbl(i)
      QP1=QP1+temp3*inv_tAkltbl(i)
      QP2=QP2+temp4*inv_tAkltvl(i)
      QP5=QP5+temp3*inv_tAkltvl(i)
      QP6=QP6+temp4*inv_tAkltbl(i)
    enddo

    temp1=Glob_PiRaised3n2/(FOUR*det_tAkl**(THREEHALF))
    temp2=THREEHALF*(tau3*tau33+tau333*tau334)*(trQP+THREEHALF*trQ*trP)
    temp3=Q1*P2+Q2*P1+Q5*P6+Q6*P5
    Pgamma=P1*tau3+P2*tau33+P5*tau333+P6*tau334
    Qgamma=Q1*tau3+Q2*tau33+Q5*tau333+Q6*tau334
    PQgamma=PQ1*tau3+PQ2*tau33+PQ5*tau333+PQ6*tau334
    QPgamma=QP1*tau3+QP2*tau33+QP5*tau333+QP6*tau334

    rPr_rQr= temp1*(THREEHALF*trQ*Pgamma+THREEHALF*trP*Qgamma+PQgamma+QPgamma+temp2+temp3)

  end function rPr_rQr

   function dXddYd(X,Y,tvk,tbk,tvl,tbl,tAl,tAk,inv_tAkl,det_tAkl,tau3,tau33,tau333,tau334,inv_tAkltAl,inv_tAkltAk)
    real(wp)   dXddYd
!arguments
    integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles
    integer       tvk(nn),tbk(nn),tvl(nn),tbl(nn)
    real(wp)   X(nn,nn),Y(nn,nn),inv_tAkl(nn,nn),tau3,tau33,tau333,tau334,det_tAkl
    real(wp)   tAk(nn,nn),tAl(nn,nn)
    real(wp)   inv_tAkltvl(nn),inv_tAkltbl(nn),tvkinv_tAkl(nn),tbkinv_tAkl(nn)
    real(wp)   inv_tAkltAl(nn,nn),inv_tAkltAk(nn,nn)
!declaration
    integer       i,j,n,k
    real(wp)   P(nn,nn),Q(nn,nn),trAkX,trAlY
    real(wp)   temp1,temp2,temp3,temp4,temp33,temp44,P1,P2,P5,P6,Q1,Q2,Q5,Q6,trP,trQ,trQP
    real(wp)   PQ1,PQ2,PQ5,PQ6,QP1,QP2,QP5,QP6
    real(wp)   inv_tAklP(nn,nn),inv_tAklQ(nn,nn),inv_tAklQP(nn,nn),tvkXtAk(nn),tbkXtAk(nn)
    real(wp)   tvkinv_tAklP(nn),tbkinv_tAklP(nn),tvkinv_tAklQ(nn),tbkinv_tAklQ(nn)
    real(wp)   XtAk(nn,nn),YtAl(nn,nn)
    real(wp)   Pgamma,Qgamma,PQgamma,QPgamma,big_term1,big_term2,big_term3,big_term4
    real(wp)   term1,term2,term3,term4,term5,term6,term7,term8,prod,gamma
    real(wp)   term9,term10,term11,term12,term13,term14,term15,term16
    real(wp)   term17,term18,term19,term20,term21,term22,term23,term24,term25
    real(wp)   tAkX(nn,nn),tAlY(nn,nn),Ytvl(nn),Ytbl(nn),Xtbk(nn),tvkX(nn),tbkX(nn),tvlY(nn)
    real(wp)   term3_1,term3_2,term4_1,term4_2,term8_1,term8_2,term9_1,term9_2,term11_1,term11_2,term16_1,term16_2
    real(wp)   term12_1,term12_2,term17_1,term17_2,term13_1,term14_2,term18_2,term19_1,temp5,temp6

    n=Glob_n
    prod=Glob_PiRaised3n2/(FOUR*det_tAkl**(THREEHALF))
    gamma=tau3*tau33+tau333*tau334

    do i=1,n
      temp1=ZERO
      temp2=ZERO
      do j=1,n
        temp1=temp1+inv_tAkl(i,j)*tvl(j)
        temp2=temp2+inv_tAkl(i,j)*tbl(j)
      enddo
      inv_tAkltvl(i)=temp1
      inv_tAkltbl(i)=temp2
    enddo

!Compute tvkinv_tAkl=tvk'*inv_tAkl
    do i=1,n
      temp1=ZERO
      temp2=ZERO
      do j=1,n
        temp1=temp1+tvk(j)*inv_tAkl(j,i)
        temp2=temp2+tbk(j)*inv_tAkl(j,i)
      enddo
      tvkinv_tAkl(i)=temp1
      tbkinv_tAkl(i)=temp2
    enddo
    !carefully look at the i,j indices !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i=1,n
      temp1=ZERO
      temp2=ZERO
      temp3=ZERO
      do j=1,n
        temp1=temp1+Y(i,j)*tvl(j)
        temp2=temp2+Y(i,j)*tbl(j)
        temp3=temp3+X(i,j)*tbk(j)
      enddo
      Ytvl(i)=temp1
      Ytbl(i)=temp2
      Xtbk(i)=temp3
    enddo

    do i=1,n
      temp1=ZERO
      temp2=ZERO
      temp3=ZERO
      do j=1,n
        temp1=temp1+tvk(j)*X(j,i)
        temp2=temp2+tbk(j)*X(j,i)
        temp3=temp3+tvl(j)*Y(j,i)
      enddo
      tvkX(i)=temp1
      tbkX(i)=temp2
      tvlY(i)=temp3
    enddo

    trAkX=ZERO
    trAlY=ZERO
    do i=1,n
      do j=1,n
        trAkX=trAkX+tAk(i,j)*X(j,i)
        trAlY=trAlY+tAl(i,j)*Y(j,i)
      enddo
    enddo

    do i=1,n
      do j=1,n
        temp1=ZERO
        temp2=ZERO
        temp3=ZERO
        temp4=ZERO
        do k=1,n
          temp1=temp1+tAk(j,k)*X(k,i)
          temp2=temp2+tAl(j,k)*Y(k,i)
          temp3=temp3+X(j,k)*tAk(k,i)
          temp4=temp4+Y(j,k)*tAl(k,i)
        enddo
        tAkX(j,i)=temp1
        tAlY(j,i)=temp2
        XtAk(j,i)=temp3
        YtAl(j,i)=temp4
      enddo
    enddo

    do i=1,n
      do j=1,n
        temp1=ZERO
        temp2=ZERO
        do k=1,n
          temp1=temp1+tAkX(j,k)*tAk(k,i)
          temp2=temp2+tAlY(j,k)*tAl(k,i)
        enddo
        P(j,i)=temp1
        Q(j,i)=temp2
      enddo
    enddo

    do i=1,n
      do j=1,n
        temp1=ZERO
        temp2=ZERO
        do k=1,n
          temp1=temp1+inv_tAkl(j,k)*P(k,i)
          temp2=temp2+inv_tAkl(j,k)*Q(k,i)
        enddo
        inv_tAklP(j,i)=temp1
        inv_tAklQ(j,i)=temp2
      enddo
    enddo

    trP=ZERO
    trQ=ZERO
    do i=1,n
      do j=1,n
        trP=trP+inv_tAkl(i,j)*P(j,i)
        trQ=trQ+inv_tAkl(i,j)*Q(j,i)
      enddo
    enddo

    Q1=ZERO
    Q2=ZERO
    Q5=ZERO
    Q6=ZERO
    P1=ZERO
    P2=ZERO
    P5=ZERO
    P6=ZERO
    do i=1,n
      temp1=ZERO
      temp2=ZERO
      temp3=ZERO
      temp4=ZERO
      do j=1,n
        temp1=temp1+tbkinv_tAkl(j)*P(j,i)
        temp2=temp2+tvkinv_tAkl(j)*P(j,i)
        temp3=temp3+tbkinv_tAkl(j)*Q(j,i)
        temp4=temp4+tvkinv_tAkl(j)*Q(j,i)
      enddo
      P1=P1+temp1*inv_tAkltbl(i)
      P2=P2+temp2*inv_tAkltvl(i)
      P5=P5+temp1*inv_tAkltvl(i)
      P6=P6+temp2*inv_tAkltbl(i)
      Q1=Q1+temp3*inv_tAkltbl(i)
      Q2=Q2+temp4*inv_tAkltvl(i)
      Q5=Q5+temp3*inv_tAkltvl(i)
      Q6=Q6+temp4*inv_tAkltbl(i)
    enddo
    Qgamma=Q1*tau3+Q2*tau33+Q5*tau333+Q6*tau334
    Pgamma=P1*tau3+P2*tau33+P5*tau333+P6*tau334

    term3_1=ZERO
    term3_2=ZERO
    term4_1=ZERO
    term4_2=ZERO
    term11_1=ZERO
    term11_2=ZERO
    term16_1=ZERO
    term16_2=ZERO
    do i=1,n
      temp1=ZERO
      temp2=ZERO
      temp3=ZERO
      temp4=ZERO
      temp33=ZERO
      temp44=ZERO
      do j=1,n
        temp1=temp1+tvkinv_tAkl(j)*tAl(j,i)
        temp2=temp2+tbkinv_tAkl(j)*tAl(j,i)
        temp3=temp3+tbkinv_tAkl(j)*tAl(j,i)
        temp4=temp4+tvkinv_tAkl(j)*tAl(j,i)
        temp33=temp33+tvkX(j)*tAk(j,i)
        temp44=temp44+tbkX(j)*tAk(j,i)
      enddo
      term3_1=term3_1+temp1*Ytvl(i)
      term3_2=term3_2+temp2*Ytvl(i)
      term4_1=term4_1+temp3*Ytbl(i)
      term4_2=term4_2+temp4*Ytbl(i)
      term11_1=term11_1+temp33*inv_tAkltvl(i)
      term11_2=term11_2+temp33*inv_tAkltbl(i)
      term16_1=term16_1+temp44*inv_tAkltbl(i)
      term16_2=term16_2+temp44*inv_tAkltvl(i)
    enddo

    do i=1,n
      temp1=ZERO
      temp2=ZERO
      temp3=ZERO
      temp4=ZERO
      do j=1,n
        temp1=temp1+tvk(j)*inv_tAklP(j,i)
        temp2=temp2+tbk(j)*inv_tAklP(j,i)
        temp3=temp3+tvk(j)*inv_tAklQ(j,i)
        temp4=temp4+tbk(j)*inv_tAklQ(j,i)
      enddo
      tvkinv_tAklP(i)=temp1
      tbkinv_tAklP(i)=temp2
      tvkinv_tAklQ(i)=temp3
      tbkinv_tAklQ(i)=temp4
    enddo
    do i=1,n
      temp1=ZERO
      temp2=ZERO
      do j=1,n
        temp1=temp1+tvkX(j)*tAk(j,i)
        temp2=temp2+tbkX(j)*tAk(j,i)
      enddo
      tvkXtAk(i)=temp1
      tbkXtAk(i)=temp2
    enddo
    term8_1=ZERO
    term8_2=ZERO
    term9_1=ZERO
    term9_2=ZERO
    term12_1=ZERO
    term12_2=ZERO
    term17_1=ZERO
    term17_2=ZERO
    term13_1=ZERO
    term14_2=ZERO
    term18_2=ZERO
    term19_1=ZERO
    do i=1,n
      temp1=ZERO
      temp2=ZERO
      temp3=ZERO
      temp4=ZERO
      temp5=ZERO
      temp6=ZERO
      do j=1,n
        temp1=temp1+tvkinv_tAklP(j)*inv_tAkltAl(j,i)
        temp2=temp2+tbkinv_tAklP(j)*inv_tAkltAl(j,i)
        temp3=temp3+tvkXtAk(j)*inv_tAklQ(j,i)
        temp4=temp4+tbkXtAk(j)*inv_tAklQ(j,i)
        temp5=temp5+tvkXtAk(j)*inv_tAkltAl(j,i)
        temp6=temp6+tbkXtAk(j)*inv_tAkltAl(j,i)
      enddo
      term8_1=term8_1+temp1*Ytvl(i)
      term8_2=term8_2+temp2*Ytvl(i)
      term9_1=term9_1+temp2*Ytbl(i)
      term9_2=term9_2+temp1*Ytbl(i)
      term12_1=term12_1+temp3*inv_tAkltvl(i)
      term12_2=term12_2+temp3*inv_tAkltbl(i)
      term17_1=term17_1+temp4*inv_tAkltbl(i)
      term17_2=term17_2+temp4*inv_tAkltvl(i)
      term13_1=term13_1+temp5*Ytvl(i)
      term14_2=term14_2+temp5*Ytbl(i)
      term18_2=term18_2+temp6*Ytvl(i)
      term19_1=term19_1+temp6*Ytbl(i)
    enddo
    term1=36*trAkX*trAlY*prod*gamma
    term2=-24*trAkX*prod*(THREEHALF*trQ*gamma + Qgamma)
    term3=24*trAkX*prod*(tau33*term3_1+tau333*term3_2)
    term4=24*trAkX*prod*(tau3*term4_1+tau334*term4_2)
    term6=-24*trAlY*prod*(THREEHALF*trP*gamma+Pgamma)
    term7=16*rPr_rQr(P,Q,tvk,tbk,inv_tAkl,det_tAkl,tau3,tau33,tau333,tau334,tvkinv_tAkl,tbkinv_tAkl,inv_tAkltvl,inv_tAkltbl)
    term8=-24*trP*prod*(tau33*term3_1+tau333*term3_2)-16*prod*(tau33*term8_1+P1*term3_1+P6*term3_2+tau333*term8_2)
    term9=-24*trP*prod*(tau3*term4_1+tau334*term4_2)-16*prod*(P2*term4_1+tau3*term9_1+tau334*term9_2+P5*term4_2)
    term11=24*trAlY*prod*(tau33*term11_1+tau334*term11_2)
    term12=-24*trQ*prod*(tau33*term11_1+tau334*term11_2)-16*prod*(tau33*term12_1+Q1*term11_1+tau334*term12_2+Q5*term11_2)
    term13=16*prod*(tau33*term13_1+term11_2*term3_2)
    term14=16*prod*(term11_1*term4_1+tau334*term14_2)
    term16=24*trAlY*prod*(tau3*term16_1+tau333*term16_2)
    term17=-24*trQ*prod*(tau3*term16_1+tau333*term16_2)&
            -16*prod*(Q2*term16_1+tau3*term17_1+Q6*term16_2+tau333*term17_2)
    term18=16*prod*(term3_1*term16_1+tau333*term18_2)
    term19=16*prod*(tau3*term19_1+term4_2*term16_2)

    big_term1=term1+term2+term3+term4+term6
    big_term2=term7+term8+term9+term11
    big_term3=term12+term13+term14
    big_term4=term16+term17+term18+term19
    dXddYd= big_term1+big_term2+big_term3+big_term4

  end function dXddYd
  function ME_over_rij_tbk(i,j,inv_tAkl,det_tAkl,tvk,tvl,tbk,tbl)
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
    real(wp)   ME_over_rij_tbk
    integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles
!Arguments:
    real(wp)   inv_tAkl(nn,nn),det_tAkl
    integer       i,j,tvl(nn),tvk(nn),tbl(nn)

    real(wp)   inv_tAkltvl(nn),inv_tAkltbl(nn),inv_tAkltbk(nn)
    real(wp)   tvkinv_tAkl(nn),tbkinv_tAkl(nn),tvlinv_tAkl(nn)
    real(wp)   tau3,tau33,tau333,tau334,m1,m3,m,tbk(nn)
!Local variables:
    integer       p,q,n
    real(wp)   Aj(nn),AjX(nn),t_JV1,t_JV2,t_JV5,t_JV6,t_J,Ajtvl,Ajtbl
    real(wp)    temp1,temp2,temp3,temp11,temp22,temp33,mu,u

    n=Glob_n

!Compute inv_tAkltvl = inv_tAkl * tvl
    do p=1,n
      temp1=ZERO
      temp2=ZERO
      temp3=ZERO
      do q=1,n
        temp1=temp1+inv_tAkl(q,p)*tvl(q)
        temp2=temp2+inv_tAkl(q,p)*tbl(q)
        temp3=temp3+inv_tAkl(q,p)*tbk(q)
      enddo
      inv_tAkltvl(p)=temp1
      inv_tAkltbl(p)=temp2
      inv_tAkltbk(p)=temp3
    enddo

!Compute tvkinv_tAkl=tvk'*inv_tAkl
    do p=1,n
      temp1=ZERO
      temp2=ZERO
      temp3=ZERO
      do q=1,n
        temp1=temp1+tvk(q)*inv_tAkl(q,p)
        temp2=temp2+tbk(q)*inv_tAkl(q,p)
        temp3=temp3+tvl(q)*inv_tAkl(q,p)
      enddo
      tvkinv_tAkl(p)=temp1
      tbkinv_tAkl(p)=temp2
      tvlinv_tAkl(p)=temp3
    enddo

!Compute tau3=tvkinv_tAkl*tvl
    tau3=ZERO
    tau33=ZERO
    tau333=ZERO
    tau334=ZERO
    do p=1,n
      tau3=tau3+tvkinv_tAkl(p)*tvl(p)
      tau33=tau33+tbkinv_tAkl(p)*tbl(p)
      tau333=tau333+tvkinv_tAkl(p)*tbl(p)
      tau334=tau334+tbkinv_tAkl(p)*tvl(p)
    enddo
    m1=tau3*tau33
    m3=tau333*tau334
    m=m1+m3  !look formula 40 in document

!Form Aj=inv_tAkl*ji        j/=i
!     Aj=inv_tAkl*(ji-jj)   j/=i
!Remember that Jii=ji*ji' and Jij=(ji-jj)*(ji-jj)'
    if (i==j) then
      do p=1,n
        Aj(p)=inv_tAkl(p,i)
      enddo
    else
      do p=1,n
        Aj(p)=inv_tAkl(p,i)-inv_tAkl(p,j)
      enddo
    endif

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
    Ajtbl=ZERO
    temp1=ZERO
    temp11=ZERO
    do p=1,n
      Ajtvl=Ajtvl+Aj(p)*tvl(p)
      Ajtbl=Ajtbl+Aj(p)*tbl(p)
      temp1=temp1+tvk(p)*Aj(p)
      temp11=temp11+tbk(p)*Aj(p)
    enddo
    t_JV1=temp11*Ajtbl
    t_JV2=temp1*Ajtvl
    t_JV5=temp11*Ajtvl
    t_JV6=temp1*Ajtbl

    mu=tau3*t_JV1+tau33*t_JV2+tau333*t_JV5+tau334*t_JV6
    u=t_JV1*t_JV2+t_JV5*t_JV6

    temp1=Glob_PiRaised3n2/(TWO*Glob_SqrtPi*det_tAkl**(THREEHALF))

    ME_over_rij_tbk=temp1*(m+ONEFIFTH*u/(t_J*t_J)-ONETHIRD*mu/t_J)/sqrt(t_J)

  end function ME_over_rij_tbk

  function ME_over_rij_tvk(i,j,inv_tAkl,det_tAkl,tvk,tvl,tbk,tbl)
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
    real(wp)   ME_over_rij_tvk
    integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles
!Arguments:
    real(wp)   inv_tAkl(nn,nn),det_tAkl
    integer       i,j,tvl(nn),tbk(nn),tbl(nn)

    real(wp)   inv_tAkltvl(nn),inv_tAkltbl(nn),inv_tAkltbk(nn)
    real(wp)   tvkinv_tAkl(nn),tbkinv_tAkl(nn),tvlinv_tAkl(nn)
    real(wp)   tau3,tau33,tau333,tau334,m1,m3,m,tvk(nn)
!Local variables:
    integer       p,q,n
    real(wp)   Aj(nn),AjX(nn),t_JV1,t_JV2,t_JV5,t_JV6,t_J,Ajtvl,Ajtbl
    real(wp)    temp1,temp2,temp3,temp11,temp22,temp33,mu,u

    n=Glob_n

!Compute inv_tAkltvl = inv_tAkl * tvl
    do p=1,n
      temp1=ZERO
      temp2=ZERO
      temp3=ZERO
      do q=1,n
        temp1=temp1+inv_tAkl(q,p)*tvl(q)
        temp2=temp2+inv_tAkl(q,p)*tbl(q)
        temp3=temp3+inv_tAkl(q,p)*tbk(q)
      enddo
      inv_tAkltvl(p)=temp1
      inv_tAkltbl(p)=temp2
      inv_tAkltbk(p)=temp3
    enddo

!Compute tvkinv_tAkl=tvk'*inv_tAkl
    do p=1,n
      temp1=ZERO
      temp2=ZERO
      temp3=ZERO
      do q=1,n
        temp1=temp1+tvk(q)*inv_tAkl(q,p)
        temp2=temp2+tbk(q)*inv_tAkl(q,p)
        temp3=temp3+tvl(q)*inv_tAkl(q,p)
      enddo
      tvkinv_tAkl(p)=temp1
      tbkinv_tAkl(p)=temp2
      tvlinv_tAkl(p)=temp3
    enddo

!Compute tau3=tvkinv_tAkl*tvl
    tau3=ZERO
    tau33=ZERO
    tau333=ZERO
    tau334=ZERO
    do p=1,n
      tau3=tau3+tvkinv_tAkl(p)*tvl(p)
      tau33=tau33+tbkinv_tAkl(p)*tbl(p)
      tau333=tau333+tvkinv_tAkl(p)*tbl(p)
      tau334=tau334+tbkinv_tAkl(p)*tvl(p)
    enddo
    m1=tau3*tau33
    m3=tau333*tau334
    m=m1+m3  !look formula 40 in document

!Form Aj=inv_tAkl*ji        j/=i
!     Aj=inv_tAkl*(ji-jj)   j/=i
!Remember that Jii=ji*ji' and Jij=(ji-jj)*(ji-jj)'
    if (i==j) then
      do p=1,n
        Aj(p)=inv_tAkl(p,i)
      enddo
    else
      do p=1,n
        Aj(p)=inv_tAkl(p,i)-inv_tAkl(p,j)
      enddo
    endif

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
    Ajtbl=ZERO
    temp1=ZERO
    temp11=ZERO
    do p=1,n
      Ajtvl=Ajtvl+Aj(p)*tvl(p)
      Ajtbl=Ajtbl+Aj(p)*tbl(p)
      temp1=temp1+tvk(p)*Aj(p)
      temp11=temp11+tbk(p)*Aj(p)
    enddo
    t_JV1=temp11*Ajtbl
    t_JV2=temp1*Ajtvl
    t_JV5=temp11*Ajtvl
    t_JV6=temp1*Ajtbl

    mu=tau3*t_JV1+tau33*t_JV2+tau333*t_JV5+tau334*t_JV6
    u=t_JV1*t_JV2+t_JV5*t_JV6

    temp1=Glob_PiRaised3n2/(TWO*Glob_SqrtPi*det_tAkl**(THREEHALF))

    ME_over_rij_tvk=temp1*(m+ONEFIFTH*u/(t_J*t_J)-ONETHIRD*mu/t_J)/sqrt(t_J)

  end function ME_over_rij_tvk

  function ME_over_rij_tbl(i,j,inv_tAkl,det_tAkl,tvk,tvl,tbk,tbl)
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
    real(wp)   ME_over_rij_tbl
    integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles
!Arguments:
    real(wp)   inv_tAkl(nn,nn),det_tAkl
    integer       i,j,tvk(nn),tvl(nn),tbk(nn)

    real(wp)   inv_tAkltvl(nn),inv_tAkltbl(nn),inv_tAkltbk(nn)
    real(wp)   tvkinv_tAkl(nn),tbkinv_tAkl(nn),tvlinv_tAkl(nn)
    real(wp)   tau3,tau33,tau333,tau334,m1,m3,m,tbl(nn)
!Local variables:
    integer       p,q,n
    real(wp)   Aj(nn),AjX(nn),t_JV1,t_JV2,t_JV5,t_JV6,t_J,Ajtvl,Ajtbl
    real(wp)    temp1,temp2,temp3,temp11,temp22,temp33,mu,u

    n=Glob_n

!Compute inv_tAkltvl = inv_tAkl * tvl
    do p=1,n
      temp1=ZERO
      temp2=ZERO
      temp3=ZERO
      do q=1,n
        temp1=temp1+inv_tAkl(q,p)*tvl(q)
        temp2=temp2+inv_tAkl(q,p)*tbl(q)
        temp3=temp3+inv_tAkl(q,p)*tbk(q)
      enddo
      inv_tAkltvl(p)=temp1
      inv_tAkltbl(p)=temp2
      inv_tAkltbk(p)=temp3
    enddo

!Compute tvkinv_tAkl=tvk'*inv_tAkl
    do p=1,n
      temp1=ZERO
      temp2=ZERO
      temp3=ZERO
      do q=1,n
        temp1=temp1+tvk(q)*inv_tAkl(q,p)
        temp2=temp2+tbk(q)*inv_tAkl(q,p)
        temp3=temp3+tvl(q)*inv_tAkl(q,p)
      enddo
      tvkinv_tAkl(p)=temp1
      tbkinv_tAkl(p)=temp2
      tvlinv_tAkl(p)=temp3
    enddo

!Compute tau3=tvkinv_tAkl*tvl
    tau3=ZERO
    tau33=ZERO
    tau333=ZERO
    tau334=ZERO
    do p=1,n
      tau3=tau3+tvkinv_tAkl(p)*tvl(p)
      tau33=tau33+tbkinv_tAkl(p)*tbl(p)
      tau333=tau333+tvkinv_tAkl(p)*tbl(p)
      tau334=tau334+tbkinv_tAkl(p)*tvl(p)
    enddo
    m1=tau3*tau33
    m3=tau333*tau334
    m=m1+m3  !look formula 40 in document

!Form Aj=inv_tAkl*ji        j/=i
!     Aj=inv_tAkl*(ji-jj)   j/=i
!Remember that Jii=ji*ji' and Jij=(ji-jj)*(ji-jj)'
    if (i==j) then
      do p=1,n
        Aj(p)=inv_tAkl(p,i)
      enddo
    else
      do p=1,n
        Aj(p)=inv_tAkl(p,i)-inv_tAkl(p,j)
      enddo
    endif

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
    Ajtbl=ZERO
    temp1=ZERO
    temp11=ZERO
    do p=1,n
      Ajtvl=Ajtvl+Aj(p)*tvl(p)
      Ajtbl=Ajtbl+Aj(p)*tbl(p)
      temp1=temp1+tvk(p)*Aj(p)
      temp11=temp11+tbk(p)*Aj(p)
    enddo
    t_JV1=temp11*Ajtbl
    t_JV2=temp1*Ajtvl
    t_JV5=temp11*Ajtvl
    t_JV6=temp1*Ajtbl

    mu=tau3*t_JV1+tau33*t_JV2+tau333*t_JV5+tau334*t_JV6
    u=t_JV1*t_JV2+t_JV5*t_JV6

    temp1=Glob_PiRaised3n2/(TWO*Glob_SqrtPi*det_tAkl**(THREEHALF))

    ME_over_rij_tbl=temp1*(m+ONEFIFTH*u/(t_J*t_J)-ONETHIRD*mu/t_J)/sqrt(t_J)

  end function ME_over_rij_tbl

  function ME_over_rij_tvl(i,j,inv_tAkl,det_tAkl,tvk,tvl,tbk,tbl)
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
    real(wp)   ME_over_rij_tvl
    integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles
!Arguments:
    real(wp)   inv_tAkl(nn,nn),det_tAkl
    integer       i,j,tvk(nn),tbk(nn),tbl(nn)

    real(wp)   inv_tAkltvl(nn),inv_tAkltbl(nn),inv_tAkltbk(nn)
    real(wp)   tvkinv_tAkl(nn),tbkinv_tAkl(nn),tvlinv_tAkl(nn)
    real(wp)   tau3,tau33,tau333,tau334,m1,m3,m,tvl(nn)
!Local variables:
    integer       p,q,n,k
    real(wp)   Aj(nn),AjX(nn),t_JV1,t_JV2,t_JV5,t_JV6,t_J,Ajtvl,Ajtbl
    real(wp)    temp1,temp2,temp3,temp11,temp22,temp33,mu,u

    n=Glob_n

!Compute inv_tAkltvl = inv_tAkl * tvl
    do p=1,n
      temp1=ZERO
      temp2=ZERO
      temp3=ZERO
      do q=1,n
        temp1=temp1+inv_tAkl(q,p)*tvl(q)
        temp2=temp2+inv_tAkl(q,p)*tbl(q)
        temp3=temp3+inv_tAkl(q,p)*tbk(q)
      enddo
      inv_tAkltvl(p)=temp1
      inv_tAkltbl(p)=temp2
      inv_tAkltbk(p)=temp3
    enddo

!Compute tvkinv_tAkl=tvk'*inv_tAkl
    do p=1,n
      temp1=ZERO
      temp2=ZERO
      temp3=ZERO
      do q=1,n
        temp1=temp1+tvk(q)*inv_tAkl(q,p)
        temp2=temp2+tbk(q)*inv_tAkl(q,p)
        temp3=temp3+tvl(q)*inv_tAkl(q,p)
      enddo
      tvkinv_tAkl(p)=temp1
      tbkinv_tAkl(p)=temp2
      tvlinv_tAkl(p)=temp3
    enddo

!Compute tau3=tvkinv_tAkl*tvl
    tau3=ZERO
    tau33=ZERO
    tau333=ZERO
    tau334=ZERO
    do p=1,n
      tau3=tau3+tvkinv_tAkl(p)*tvl(p)
      tau33=tau33+tbkinv_tAkl(p)*tbl(p)
      tau333=tau333+tvkinv_tAkl(p)*tbl(p)
      tau334=tau334+tbkinv_tAkl(p)*tvl(p)
    enddo
    m1=tau3*tau33
    m3=tau333*tau334
    m=m1+m3  !look formula 40 in document

!Form Aj=inv_tAkl*ji        j/=i
!     Aj=inv_tAkl*(ji-jj)   j/=i
!Remember that Jii=ji*ji' and Jij=(ji-jj)*(ji-jj)'
    if (i==j) then
      do p=1,n
        Aj(p)=inv_tAkl(p,i)
      enddo
    else
      do p=1,n
        Aj(p)=inv_tAkl(p,i)-inv_tAkl(p,j)
      enddo
    endif

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
    Ajtbl=ZERO
    temp1=ZERO
    temp11=ZERO
    do p=1,n
      Ajtvl=Ajtvl+Aj(p)*tvl(p)
      Ajtbl=Ajtbl+Aj(p)*tbl(p)
      temp1=temp1+tvk(p)*Aj(p)
      temp11=temp11+tbk(p)*Aj(p)
    enddo
    t_JV1=temp11*Ajtbl
    t_JV2=temp1*Ajtvl
    t_JV5=temp11*Ajtvl
    t_JV6=temp1*Ajtbl

    mu=tau3*t_JV1+tau33*t_JV2+tau333*t_JV5+tau334*t_JV6
    u=t_JV1*t_JV2+t_JV5*t_JV6

    temp1=Glob_PiRaised3n2/(TWO*Glob_SqrtPi*det_tAkl**(THREEHALF))

    ME_over_rij_tvl=temp1*(m+ONEFIFTH*u/(t_J*t_J)-ONETHIRD*mu/t_J)/sqrt(t_J)

  end function ME_over_rij_tvl

  function ME_dXd(X,tvk,tvl,inv_tAkltvl,inv_tAkl,tAk,tAl,inv_tAkltAl,Skl,tau3)
    real(wp)   ME_dXd
    integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles
    real(wp)   X(nn,nn),tAk(nn,nn),tAl(nn,nn),inv_tAkl(nn,nn),inv_tAkltAl(nn,nn)
    integer       i,j,n,k,tvk(nn),tvl(nn)
    real(wp)   inv_tAkltAlX(nn,nn),inv_tAkltAlXtAk(nn,nn),tvkinv_tAkltAlX(nn),inv_tAkltvl(nn)
    real(wp)   temp1, Skl,tau,tau1,tau2,tau3
    n=Glob_n
    do i=1,n
      do j=1,n
        temp1=ZERO
        do k=1,n
          temp1=temp1+inv_tAkltAl(j,k)*X(k,i)
        enddo
        inv_tAkltAlX(j,i)=temp1
      enddo
    enddo
    tau1=ZERO
    do i=1,n
      temp1=ZERO
      do k=1,n
        temp1=temp1+inv_tAkltAlX(i,k)*tAk(k,i)
      enddo
      tau1=tau1+temp1
    enddo
    do i=1,n
      temp1=ZERO
      do j=1,n
        temp1=temp1+tvk(j)*inv_tAkltAlX(j,i)
      enddo
      tvkinv_tAkltAlX(i)=temp1
    enddo
    tau2=ZERO
    do i=1,n
      temp1=ZERO
      do j=1,n
        temp1=temp1+tvkinv_tAkltAlX(j)*tAk(j,i)
      enddo
      tau2=tau2+temp1*inv_tAkltvl(i)
    enddo
    ME_dXd=Skl*(SIX*tau1+FOUR*tau2/tau3)
  end function ME_dXd

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

  subroutine spinPreCalc(n, nFactorial, parityFactor, SSFmassChargeCoefficient, SSNCmassChargeCoefficient, &
                         SOmassChargeCoefficient, AMMmassChargeCoefficient, AMMFinmassChargeCoefficient, &
                         ketMatrix, spatialYoung0, spatialYoung1, &
                         SSNCspinME, SiMinusME, SiPlusME, SziME, spinFreeME, SiSjME)
    use spinStuff
    implicit none

    !input vars:
    character(len = maxLen), intent(in) :: spatialYoung0, spatialYoung1
    integer, intent(in) :: n, nFactorial

    !output vars:
    real(wp), dimension(nFactorial), intent(out) :: parityFactor
    real(wp), dimension(n, n, 4), intent(out) :: SOmassChargeCoefficient, AMMmassChargeCoefficient, &
                                                    AMMFinmassChargeCoefficient
    real(wp), dimension(n, n), intent(out) :: SSNCmassChargeCoefficient, SSFmassChargeCoefficient
    real(wp), dimension(n, n, nFactorial), intent(out) :: ketMatrix
    real(wp), dimension(n, n, nFactorial), intent(out) :: SSNCspinME, SiSjME
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
    SSFmassChargeCoefficient = ZERO

    do i = 1, n
      do j = 1, n
        SSFmassChargeCoefficient(i, j) = -Glob_PseudoCharge(i) * Glob_PseudoCharge(j) / &
                                         (Glob_Mass(i + 1) * Glob_Mass(j + 1)) * EIGHT * Glob_Pi / THREE
      enddo
    enddo

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

    !Singlet wf
    call getSpinFunction(n, nFactorial, mySpatialYoung1, allPermutations, parities, &
                         finalSpinFunction0, primitives0, numberOfPrimitives0)

    !Triplet wf
    call getSpinFunction(n, nFactorial, mySpatialYoung0, allPermutations, parities, &
                         finalSpinFunction1, primitives1, numberOfPrimitives1)

    call getSpinOpMeanValues(n, nFactorial, allPermutations, finalSpinFunction0, finalSpinFunction1, primitives0, primitives1, &
                             numberOfPrimitives0, numberOfPrimitives1, spinFreeME, SSNCspinME, SiMinusME, SiPlusME, SziME, SiSjME)

    ketMatrix = ZERO
    do i = 1, nFactorial
      do k = 1, n
        do l = 1, n
          !note the transposition here
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
                                         SSNCkl, SO1kl, SO2kl, AMM1kl, AMM2kl, AMM1finkl, AMM2finkl)
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

    real(wp), intent(out)  :: SO1kl, SO2kl, AMM1finkl, AMM2finkl, AMM1kl, AMM2kl, SSNCkl
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
    real(wp) :: commonFactor, gamma, gamma_diag, jiAlAklinvVl, &
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

    !common factor (ONEHALF - for consistent normalization with Skl)
    if (Glob_selectTransition == 1) commonFactor = Glob_PiRaised3n2 / (Glob_SqrtPi * det_tAkl * sqrt(det_tAkl))
    if (Glob_selectTransition == 2) commonFactor = SQRT(THREE) / TWO * Glob_PiRaised3n2 / (Glob_SqrtPi * det_tAkl * sqrt(det_tAkl))

    SO1kl = ZERO
    SO2kl = ZERO

    AMM1kl = ZERO
    AMM2kl = ZERO

    AMM1finkl = ZERO
    AMM2finkl = ZERO

    SSNCkl = ZERO

    WkAklinvWl = inv_tAkl(pmm_k, pmm_l)
    WkAklinvVl = inv_tAkl(pmm_k, pm_l)
    VkAklinvWl = inv_tAkl(pm_k, pmm_l)
    VkAklinvVl = inv_tAkl(pm_k, pm_l)

    do indexI = 1, n

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

      temp1 = temp1010 + temp0101 + temp1001 + temp0110

      SO1kl = SO1kl + SOspinME(indexI) * SOmassChargeCoefficient(indexI, indexI, 1) * temp1
      AMM1kl = AMM1kl + SOspinME(indexI) * AMMmassChargeCoefficient(indexI, indexI, 1) * temp1
      AMM1finkl = AMM1finkl + SOspinME(indexI) * AMMFinmassChargeCoefficient(indexI, indexI, 1) * temp1

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

        temp1 = temp1010 + temp0101 + temp1001 + temp0110

        SO2kl = SO2kl + SOspinME(indexI) * SOmassChargeCoefficient(indexI, indexI, 2) * temp1
        AMM2finkl = AMM2finkl + SOspinME(indexI) * AMMFinmassChargeCoefficient(indexI, indexI, 2) * temp1

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

        temp1 = temp1010 + temp0101 + temp1001 + temp0110

        SO2kl = SO2kl + SOspinME(indexI) * SOmassChargeCoefficient(indexI, indexJ, 3) * temp1
        AMM2kl  = AMM2kl + SOspinME(indexI) * AMMmassChargeCoefficient(indexI, indexJ, 3) * temp1
        AMM2finkl  = AMM2finkl + SOspinME(indexI) * AMMfinmassChargeCoefficient(indexI, indexJ, 3) * temp1

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

        temp1 = temp1010 + temp0101 + temp1001 + temp0110

        SO2kl = SO2kl + SOspinME(indexI) * SOmassChargeCoefficient(indexI, indexJ, 4) * temp1
        AMM2kl = AMM2kl + SOspinME(indexI) * AMMmassChargeCoefficient(indexI, indexJ, 4) * temp1
        AMM2finkl = AMM2finkl + SOspinME(indexI) * AMMfinmassChargeCoefficient(indexI, indexJ, 4) * temp1

        if (Glob_selectTransition == 2 .or. indexJ <= indexI .or. abs(SSNCspinME(indexI, indexJ)) < localEps ) cycle
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

        temp1 = ONETHIRD * (temp1010 + temp1001 + temp0110 + temp0101) !ONETHIRD - from common factor and spin part (1/sqrt(6))

        SSNCkl = SSNCkl + SSNCspinME(indexI, indexJ) * SSNCmassChargeCoefficient(indexI, indexJ) * temp1

      enddo ! indexJ cycle
    enddo ! indexI cycle

    SSNCkl = SSNCkl * commonFactor
    SO1kl = SO1kl * commonFactor
    SO2kl = SO2kl * commonFactor
    AMM1kl = AMM1kl * commonFactor
    AMM2kl = AMM2kl * commonFactor
    AMM1finkl = AMM1finkl * commonFactor
    AMM2finkl = AMM2finkl * commonFactor

  end subroutine spinDependentMatrixElements

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

!function SG_ME_rXr_over_rij(i,j,X,inv_tAkl,det_tAkl)
!!function SG_ME_rXr_over_rij computes the following matrix element:
!!<\tilde psi_k| (r' X r)/r_ij |\tilde psi_l>
!!where psi_k = exp[-(r' Ak r)] is a Simple Gaussian wavefunction
!!Here X is a some real symmetric matrix. If matrix X is not symmetric
!!then user needs to symmetrize it before calling this function.
!!Input:
!!            X  :: n x n real matrix
!!      inv_tAkl :: n x n real matrix where the inverse of tAk+tAl is stored
!!           t_V :: scalar, t_V = tau3 = tr[inv_tAkl*tvl*tvk']
!!           Skl :: scalar, overlap Skl=<\tilde phi_k|\tilde phi_l>
!real(wp)   SG_ME_rXr_over_rij
!integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles
!!Arguments:
!real(wp)   X(nn,nn),inv_tAkl(nn,nn),det_tAkl
!integer       i,j
!!Local variables:
!integer       p,q,n
!real(wp)   temp1,temp2,temp3
!real(wp)   Aj(nn),AjX(nn)
!real(wp)   t_J,t_X,t_XJ
!
!n=Glob_n
!!Form Aj=inv_tAkl*ji        j/=i
!!     Aj=inv_tAkl*(ji-jj)   j/=i
!!Remember that Jii=ji*ji' and Jij=(ji-jj)*(ji-jj)'
!if (i==j) then
! do p=1,n
!   Aj(p)=inv_tAkl(p,i)
! enddo
!else
! do p=1,n
!   Aj(p)=inv_tAkl(p,i)-inv_tAkl(p,j)
! enddo
!endif
!
!!Compute AjX'=Aj'*X
!do p=1,n
!  temp1=ZERO
!  do q=1,n
!    temp1=temp1+Aj(q)*X(q,p)
!  enddo
!  AjX(p)=temp1
!enddo
!
!!Compute t_J=tr[inv_tAkl*Jij]
!if (i==j) then
!  t_J=inv_tAkl(i,i)
!else
!  t_J=inv_tAkl(i,i)+inv_tAkl(j,j)-inv_tAkl(j,i)-inv_tAkl(j,i)
!endif
!
!!Compute t_XJ=tr[inv_tAkl*X*inv_tAkl*Jij]=AjX'*Aj
!t_XJ=ZERO
!do p=1,n
!  t_XJ=t_XJ+AjX(p)*Aj(p)
!enddo
!
!!Compute t_X=tr[inv_tAkl*X]
!t_X=ZERO
!do p=1,n
!  do q=1,n
!    t_X=t_X+inv_tAkl(q,p)*X(q,p)
!  enddo
!enddo
!
!temp1=Glob_PiRaised3n2/(Glob_SqrtPi*det_tAkl**(THREEHALF))
!temp3=1/t_J
!SG_ME_rXr_over_rij=temp1*temp3*sqrt(temp3)*(THREE*t_J*t_X - t_XJ)
!
!end function SG_ME_rXr_over_rij
!function ME_d_X_over_rij_d1(i,j,X,tAl,inv_tAkl,det_tAkl,tvk,tvl,tbk,tbl)
!
!real(wp)   ME_d_X_over_rij_d1
!integer,parameter :: nn=Glob_AllowedNumOfPseudoParticles
!!Arguments:
!real(wp)   X(nn,nn),tAl(nn,nn),inv_tAkl(nn,nn),det_tAkl
!integer       i,j,tvk(nn),tvl(nn),tbk(nn),tbl(nn)
!
!!Local variables:
!integer       c,s,p,q,n,k
!real(wp)   Xtvl(nn),Xtbl(nn),tAlX(nn,nn),tAlXtAl(nn,nn),tAlXtvl(nn),tAlXtbl(nn)
!real(wp)   temp1,temp2,temp3,temp4,t_tAlX
!
!
!!Doing multiplication Xtvl=X*tvl,tvkX=tvk*X
!n=Glob_n
!do s=1,n
!  temp1=ZERO
!  temp2=ZERO
!  do c=1,n
!        temp1=temp1+X(c,s)*tvl(c)
!        temp2=temp2+X(c,s)*tbl(c)
!  enddo
!  Xtvl(s)=temp1
!  Xtbl(s)=temp2
!enddo
!
!do s=1,n
!  do c=1,n
!    temp1=ZERO
!    do k=1,n
!      temp1=temp1+tAl(c,k)*X(k,s)
!    enddo
!    tAlX(c,s)=temp1
!  enddo
!enddo
!
!!Doing multiplication tAlXtAk=tAlX*tAk,tAkXtAl=tAkX*tAl
!do s=1,n
!  do c=1,n
!    temp1=ZERO
!    do k=1,n
!      temp1=temp1+tAlX(c,k)*tAl(k,s)
!    enddo
!    tAlXtAl(c,s)=temp1
!  enddo
!enddo
!
!do s=1,n
!  temp1=ZERO
!  temp2=ZERO
!  do c=1,n
!        temp1=temp1+tAlX(c,s)*tvl(c)
!        temp2=temp2+tAlX(c,s)*tbl(c)
!  enddo
!  tAlXtvl(s)=temp1
!  tAlXtbl(s)=temp2
!enddo
!
!t_tAlX=ZERO
!do p=1,n
!  do q=1,n
!    t_tAlX=t_tAlX+tAl(q,p)*X(q,p)
!  enddo
!enddo
!
!temp1=-6*t_tAlX*ME_over_rij(i,j,inv_tAkl,det_tAkl,tvk,tvl,tbk,tbl)
!temp2=4*ME_rXr_over_rij(i,j,tAlXtAl,inv_tAkl,det_tAkl,tvk,tvl,tbk,tbl)
!temp3=-4*ME_over_rij_tvl(i,j,inv_tAkl,det_tAkl,tvk,tAlXtvl,tbk,tbl)
!temp4=-4*ME_over_rij_tbl(i,j,inv_tAkl,det_tAkl,tvk,tvl,tbk,tAlXtbl)
!ME_d_X_over_rij_d1=-(temp1+temp2+temp3+temp4)
!end function ME_d_X_over_rij_d1

end module matelem
