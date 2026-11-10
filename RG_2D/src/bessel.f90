module bessel
! High-precision working-kind versions of the SLATEC DBESK, DBSKNU, and DASYIK routines.
!
! Enables evaluation of modified Bessel functions of the second kind K_nu(x)
! for real order nu and positive argument x with working kind precision.

  use globvars
  use bessel_coeff
  implicit none

contains

  subroutine DBESK(X,FNU,KODE,N,Y,NZ,F1MinusOne)
  ! High-precision working-kind version of the SLATEC DBESK routine.
  !
  ! ECGPACK modifications:
  !   * floating-point arithmetic uses the kind exported by the local wp_def;
  !   * integer-order starting values are also obtained from DBSKNU, avoiding
  !     the separate order-zero/order-one Chebyshev routine dependency chain;
  !   * the optional F1MinusOne output evaluates
  !       sqrt(2*X/pi)*exp(X)*K_1(X)-1
  !     directly from its large-2*X expansion when needed, avoiding
  !     cancellation when the result is close to zero;
  !   * XERMSG calls are replaced by descriptive ERROR STOP statements.
  ! Derived from the public-domain SLATEC 4.1 source (author D. E. Amos, SNLA),
  ! https://netlib.org/slatec/; see THIRD-PARTY-NOTICES.md.

    integer, intent(in) :: KODE,N
    integer, intent(out) :: NZ
    real(wp), intent(in) :: X,FNU
    real(wp), intent(out), dimension(N) :: Y
    real(wp), intent(out), optional :: F1MinusOne

    integer :: I,J,K,MZ,NB,ND,NN,NUD,NULIM(2)
    integer, external :: I1MACH
    real(wp) :: CN,DNU,ELIM,ETX,FLGIK,FN,FNN,GLN,GNU,RTZ,S,S1,S2,T,TM
    real(wp) :: TAU,TERM,TRX,W(2),XLIM,ZN
    real(wp), external :: D1MACH

    ! The original limits 35 and 70 were selected for about 14 decimal digits.
    ! Retain them for double precision and increase them conservatively when the
    ! requested working precision carries more digits.
    if (precision(ONE).le.15) then
      NULIM = [35,70]
    else
      NULIM = [3*precision(ONE),6*precision(ONE)]
    endif
    NN = -I1MACH(15)
    ELIM = log(TEN)*(NN*D1MACH(5)-THREE)
    XLIM = D1MACH(1)*1.0E+3_wp
    if (KODE.lt.1 .or. KODE.gt.2) error stop 'DBESK: KODE must be 1 or 2'
    if (FNU.lt.ZERO) error stop 'DBESK: FNU must be nonnegative'
    if (X.le.ZERO) error stop 'DBESK: X must be greater than zero'
    if (N.lt.1) error stop 'DBESK: N must be greater than zero'
    if (X.lt.XLIM) error stop 'DBESK: X is too small for the working kind'
    if (present(F1MinusOne)) then
      if (FNU.ne.ONE .or. KODE.ne.2) &
        error stop 'DBESK: F1MinusOne requires FNU=1 and KODE=2'
    endif
    ETX = KODE-1

    ND = N
    NZ = 0
    NUD = int(FNU)
    DNU = FNU-NUD
    GNU = FNU
    NN = min(2,ND)
    FN = FNU+N-1
    FNN = FN
    if (FN.lt.TWO) goto 150

    ZN = X/FN
    if (ZN.eq.ZERO) goto 320
    RTZ = sqrt(ONE+ZN*ZN)
    GLN = log((ONE+RTZ)/ZN)
    T = RTZ*(ONE-ETX)+ETX/(ZN+RTZ)
    CN = -FN*(T-GLN)
    if (CN.gt.ELIM) goto 320
    if (NUD.lt.NULIM(NN)) goto 30
    if (NN.eq.1) goto 20
  10 continue
    FN = GNU
    ZN = X/FN
    RTZ = sqrt(ONE+ZN*ZN)
    GLN = log((ONE+RTZ)/ZN)
    T = RTZ*(ONE-ETX)+ETX/(ZN+RTZ)
    CN = -FN*(T-GLN)
  20 continue
    if (CN.lt.-ELIM) goto 230

    FLGIK = -ONE
    call DASYIK(X,GNU,KODE,FLGIK,RTZ,CN,NN,Y)
    if (NN.eq.1) goto 240
    TRX = TWO/X
    TM = (GNU+GNU+TWO)/X
    goto 130

  30 continue
    if (KODE.ne.2 .and. X.gt.ELIM) goto 230
    NB = 2
    if (NUD.eq.0 .and. ND.eq.1) NB = 1
    call DBSKNU(X,DNU,KODE,NB,W,NZ)
    S1 = W(1)
    if (NB.eq.1) goto 120
    S2 = W(2)
    TRX = TWO/X
    TM = (DNU+DNU+TWO)/X
    if (ND.eq.1) NUD = NUD-1
    if (NUD.gt.0) goto 100
    if (ND.gt.1) goto 120
    S1 = S2
    goto 120
  100 continue
    do I=1,NUD
      S = S2
      S2 = TM*S2+S1
      S1 = S
      TM = TM+TRX
    enddo
    if (ND.eq.1) S1 = S2
  120 continue
    Y(1) = S1
    if (ND.eq.1) goto 240
    Y(2) = S2
  130 continue
    if (ND.eq.2) goto 240
    do I=3,ND
      Y(I) = TM*Y(I-1)+Y(I-2)
      TM = TM+TRX
    enddo
    goto 240

  150 continue
    if (KODE.ne.2 .and. X.gt.ELIM) goto 230
    if (FN.gt.ONE) then
      if (-FN*(log(X)-log(TWO)).gt.ELIM) goto 320
    endif
    call DBSKNU(X,FNU,KODE,ND,Y,MZ)
    goto 240

  230 continue
    NUD = NUD+1
    ND = ND-1
    if (ND.eq.0) goto 240
    NN = min(2,ND)
    GNU = GNU+ONE
    if (FNN.lt.TWO) goto 230
    if (NUD.lt.NULIM(NN)) goto 230
    goto 10
  240 continue
    NZ = N-ND
    if (NZ.ne.0) then
      if (ND.ne.0) then
        do I=1,ND
          J = N-I+1
          K = ND-I+1
          Y(J) = Y(K)
        enddo
      endif
      Y(1:NZ) = ZERO
    endif

    if (present(F1MinusOne)) then
      TAU = TWO*X
      if (TAU.gt.50.0_wp) then
        F1MinusOne = ZERO
        TERM = ONE
        do K=1,50
          TERM = TERM*(FOUR-(TWO*real(K,wp)-ONE)**2) &
                     /(FOUR*real(K,wp)*TAU)
          F1MinusOne = F1MinusOne+TERM
          if (abs(TERM).le.epsilon(ONE)*abs(F1MinusOne)) exit
        enddo
      else
        F1MinusOne = sqrt(TAU/Glob_Pi)*Y(1)-ONE
      endif
    endif
    return

  320 continue
    error stop 'DBESK: overflow; FNU or N is too large, or X is too small'
  end subroutine DBESK


  subroutine DBSKNU(X,FNU,KODE,N,Y,NZ)
  ! High-precision working-kind version of the SLATEC DBSKNU routine.
  !
  ! ECGPACK modifications:
  !   * floating-point arithmetic uses the kind exported by the local wp_def;
  !   * the standard generic GAMMA intrinsic replaces SLATEC DGAMMA;
  !   * XERMSG calls are replaced by descriptive ERROR STOP statements.
  !   * precision-dependent convergence limits replace double-precision limits.
  ! Derived from the public-domain SLATEC 4.1 source (author D. E. Amos, SNLA),
  ! https://netlib.org/slatec/; see THIRD-PARTY-NOTICES.md.

    integer, intent(in) :: KODE,N
    integer, intent(out) :: NZ
    real(wp), intent(in) :: X,FNU
    real(wp), intent(out), dimension(N) :: Y

    integer :: I,IFLAG,INU,J,K,KK,KODED,NN
    integer, external :: I1MACH
    real(wp) :: A(512),AK,A1,A2,B(512),BK,CK,CKABS,COEF,CX,DK,DNU,DNU2
    real(wp) :: ELIM,ETEST,EX,F,FC,FHS,FK,FKS,FLRX,FMU,G1,G2,P,PT,P1,P2,Q
    real(wp) :: PREVABS,RX,S,SMU,SQK,ST,S1,S2,TM,TOL,T1,T2,X2
    real(wp), external :: D1MACH
    real(wp), parameter :: X1=2.0_wp
    real(wp), parameter :: RTHPI=sqrt(0.5_wp*Glob_Pi)

    KK = -I1MACH(15)
    ELIM = log(TEN)*(KK*D1MACH(5)-3.0_wp)
    TOL = epsilon(1.0_wp)
    ! The optimally truncated large-X expansion has an error of order EXP(-2X).
    ! Move its lower boundary outward as the requested precision increases.
    X2 = max(17.0_wp,-ONEHALF*log(TOL)+TWO)
    if (X.le.ZERO) error stop 'DBSKNU: X must be greater than zero'
    if (FNU.lt.ZERO) error stop 'DBSKNU: FNU must be nonnegative'
    if (KODE.lt.1 .or. KODE.gt.2) error stop 'DBSKNU: KODE must be 1 or 2'
    if (N.lt.1) error stop 'DBSKNU: N must be greater than zero'
    NZ = 0
    IFLAG = 0
    KODED = KODE
    RX = TWO/X
    INU = int(FNU+ONEHALF)
    DNU = FNU-INU
    if (abs(DNU).eq.ONEHALF) goto 120
    DNU2 = ZERO
    if (abs(DNU).ge.TOL) DNU2 = DNU*DNU
    if (X.gt.X1) goto 120

    A1 = ONE-DNU
    A2 = ONE+DNU
    T1 = ONE/gamma(A1)
    T2 = ONE/gamma(A2)
    if (abs(DNU).gt.ONEHALF) goto 40
    S = Bessel_RecipGammaC(1)
    AK = ONE
    do K=2,Bessel_RecipGammaCount
      AK = AK*DNU2
      TM = Bessel_RecipGammaC(K)*AK
      S = S+TM
      if (abs(TM).le.TOL*abs(S)) exit
    enddo
    G1 = -S
    goto 50
  40 continue
    G1 = (T1-T2)/(DNU+DNU)
  50 continue
    G2 = (T1+T2)*ONEHALF
    SMU = ONE
    FC = ONE
    FLRX = log(RX)
    FMU = DNU*FLRX
    if (DNU.ne.ZERO) then
      FC = DNU*Glob_Pi
      FC = FC/sin(FC)
      if (FMU.ne.ZERO) SMU = sinh(FMU)/FMU
    endif
    F = FC*(G1*cosh(FMU)+G2*FLRX*SMU)
    FC = exp(FMU)
    P = ONEHALF*FC/T2
    Q = ONEHALF/(FC*T1)
    AK = ONE
    CK = ONE
    BK = ONE
    S1 = F
    S2 = P
    if (INU.gt.0 .or. N.gt.1) goto 90
    if (X.lt.TOL) goto 80
    CX = X*X*ONEFOURTH
  70 continue
    F = (AK*F+P+Q)/(BK-DNU2)
    P = P/(AK-DNU)
    Q = Q/(AK+DNU)
    CK = CK*CX/AK
    T1 = CK*F
    S1 = S1+T1
    BK = BK+AK+AK+ONE
    AK = AK+ONE
    S = abs(T1)/(ONE+abs(S1))
    if (S.gt.TOL) goto 70
  80 continue
    Y(1) = S1
    if (KODED.eq.1) return
    Y(1) = S1*exp(X)
    return
  90 continue
    if (X.lt.TOL) goto 110
    CX = X*X*ONEFOURTH
  100 continue
    F = (AK*F+P+Q)/(BK-DNU2)
    P = P/(AK-DNU)
    Q = Q/(AK+DNU)
    CK = CK*CX/AK
    T1 = CK*F
    S1 = S1+T1
    T2 = CK*(P-AK*F)
    S2 = S2+T2
    BK = BK+AK+AK+ONE
    AK = AK+ONE
    S = abs(T1)/(ONE+abs(S1))+abs(T2)/(ONE+abs(S2))
    if (S.gt.TOL) goto 100
  110 continue
    S2 = S2*RX
    if (KODED.eq.1) goto 170
    F = exp(X)
    S1 = S1*F
    S2 = S2*F
    goto 170

  120 continue
    COEF = RTHPI/sqrt(X)
    if (KODED.eq.2) goto 130
    if (X.gt.ELIM) goto 330
    COEF = COEF*exp(-X)
  130 continue
    if (abs(DNU).eq.ONEHALF) goto 340
    if (X.gt.X2) goto 280

    ETEST = cos(Glob_Pi*DNU)/(Glob_Pi*X*TOL)
    FKS = ONE
    FHS = ONEFOURTH
    FK = ZERO
    CK = X+X+TWO
    P1 = ZERO
    P2 = ONE
    K = 0
  140 continue
    K = K+1
    if (K.gt.size(A)) error stop 'DBSKNU: Miller workspace exhausted'
    FK = FK+ONE
    AK = (FHS-DNU2)/(FKS+FK)
    BK = CK/(FK+ONE)
    PT = P2
    P2 = BK*P2-AK*P1
    P1 = PT
    A(K) = AK
    B(K) = BK
    CK = CK+TWO
    FKS = FKS+FK+FK+ONE
    FHS = FHS+FK+FK
    if (ETEST.gt.FK*P1) goto 140
    KK = K
    S = ONE
    P1 = ZERO
    P2 = ONE
    do I=1,K
      PT = P2
      P2 = (B(KK)*P2-P1)/A(KK)
      P1 = PT
      S = S+P2
      KK = KK-1
    enddo
    S1 = COEF*(P2/S)
    if (INU.gt.0 .or. N.gt.1) goto 160
    goto 200
  160 continue
    S2 = S1*(X+DNU+ONEHALF-P1/P2)/X

  170 continue
    CK = (DNU+DNU+TWO)/X
    if (N.eq.1) INU = INU-1
    if (INU.gt.0) goto 180
    if (N.gt.1) goto 200
    S1 = S2
    goto 200
  180 continue
    do I=1,INU
      ST = S2
      S2 = CK*S2+S1
      S1 = ST
      CK = CK+RX
    enddo
    if (N.eq.1) S1 = S2
  200 continue
    if (IFLAG.eq.1) goto 220
    Y(1) = S1
    if (N.eq.1) return
    Y(2) = S2
    if (N.eq.2) return
    do I=3,N
      Y(I) = CK*Y(I-1)+Y(I-2)
      CK = CK+RX
    enddo
    return

  220 continue
    S = -X+log(S1)
    Y(1) = ZERO
    NZ = 1
    if (S.lt.-ELIM) goto 230
    Y(1) = exp(S)
    NZ = 0
  230 continue
    if (N.eq.1) return
    S = -X+log(S2)
    Y(2) = ZERO
    NZ = NZ+1
    if (S.lt.-ELIM) goto 240
    NZ = NZ-1
    Y(2) = exp(S)
  240 continue
    if (N.eq.2) return
    KK = 2
    if (NZ.lt.2) goto 260
    do I=3,N
      KK = I
      ST = S2
      S2 = CK*S2+S1
      S1 = ST
      CK = CK+RX
      S = -X+log(S2)
      NZ = NZ+1
      Y(I) = ZERO
      if (S.lt.-ELIM) cycle
      Y(I) = exp(S)
      NZ = NZ-1
      goto 260
    enddo
    return
  260 continue
    if (KK.eq.N) return
    S2 = S2*CK+S1
    CK = CK+RX
    KK = KK+1
    Y(KK) = exp(-X+log(S2))
    if (KK.eq.N) return
    KK = KK+1
    do I=KK,N
      Y(I) = CK*Y(I-1)+Y(I-2)
      CK = CK+RX
    enddo
    return

  280 continue
    NN = 2
    if (INU.eq.0 .and. N.eq.1) NN = 1
    DNU2 = DNU+DNU
    FMU = ZERO
    if (abs(DNU2).ge.TOL) FMU = DNU2*DNU2
    EX = X*EIGHT
    S2 = ZERO
    do K=1,NN
      S1 = S2
      S = ONE
      AK = ZERO
      CK = ONE
      SQK = ONE
      DK = EX
      PREVABS = huge(ONE)
      do J=1,200
        CK = CK*(FMU-SQK)/DK
        CKABS = abs(CK)
        if (CKABS.gt.PREVABS) exit
        S = S+CK
        DK = DK+EX
        AK = AK+EIGHT
        SQK = SQK+AK
        if (CKABS.le.TOL*max(ONE,abs(S))) exit
        PREVABS = CKABS
      enddo
      S2 = S*COEF
      FMU = FMU+EIGHT*DNU+FOUR
    enddo
    if (NN.gt.1) goto 170
    S1 = S2
    goto 200
  330 continue
    KODED = 2
    IFLAG = 1
    goto 120
  340 continue
    S1 = COEF
    S2 = COEF
    goto 170
  end subroutine DBSKNU


  subroutine DASYIK(X,FNU,KODE,FLGIK,RA,ARG,IN,Y)
  ! High-precision working-kind version of the SLATEC DASYIK routine.
  !
  ! The Debye table contains exact-rational coefficients through U_24, emitted
  ! with 50 significant decimal digits.  Derived from the public-domain SLATEC 4.1
  ! source (author D. E. Amos, SNLA); see THIRD-PARTY-NOTICES.md.

    integer, intent(in) :: IN,KODE
    real(wp), intent(in) :: X,FNU,FLGIK
    real(wp), intent(inout) :: RA,ARG
    real(wp), intent(out), dimension(IN) :: Y

    integer :: J,JN,K,KK,L,LEND
    real(wp) :: AK,AP,COEF,ETX,FN,GLN,S1,S2,T,TOL,T2,Z
    real(wp), parameter :: CON(2) = [sqrt(ONEHALF/Glob_Pi), &
                                    sqrt(ONEHALF*Glob_Pi)]
    TOL = epsilon(ONE)
    FN = FNU
    Z = (THREE-FLGIK)/TWO
    KK = int(Z)
    do JN=1,IN
      if (JN.ne.1) then
        FN = FN-FLGIK
        Z = X/FN
        RA = sqrt(ONE+Z*Z)
        GLN = log((ONE+RA)/Z)
        ETX = KODE-1
        T = RA*(ONE-ETX)+ETX/(Z+RA)
        ARG = FN*(T-GLN)*FLGIK
      endif
      COEF = exp(ARG)
      T = ONE/RA
      T2 = T*T
      T = sign(T/FN,FLGIK)
      S2 = ONE
      AP = ONE
      do K=1,Bessel_DebyeMaxOrder
        L = 1+(K-1)*(K+2)/2
        LEND = L+K
        S1 = Bessel_DebyeC(L)
        do J=L+1,LEND
          S1 = S1*T2+Bessel_DebyeC(J)
        enddo
        AP = AP*T
        AK = AP*S1
        S2 = S2+AK
        if (max(abs(AK),abs(AP)).le.TOL*max(ONE,abs(S2))) exit
      enddo
      Y(JN) = S2*COEF*sqrt(abs(T))*CON(KK)
    enddo
  end subroutine DASYIK

end module bessel
