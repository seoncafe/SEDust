module fixed_orientation_amplitude
   !! Amplitude matrix of a nonspherical particle in a fixed orientation.
   !!
   !! AMPL and its helper VIGAMPL are Mishchenko's routines of ampld.lp.f
   !! ("Calculation of the amplitude matrix for a nonspherical particle in
   !! a fixed orientation", Appl. Opt. 39, 1026-1031, 2000), migrated from
   !! fixed to free form.  Their names are kept unchanged so that they stay
   !! checkable against the original source.
   !!
   !! Differences from the original, all structural:
   !!
   !!   - The converged T-matrix is read from the caller-owned workspace
   !!     component WORK%FULL_TSTORE instead of COMMON /TMAT/.  Slot k of
   !!     its fourth dimension holds, in order, Re T11, Re T12, Re T21,
   !!     Re T22, Im T11, Im T12, Im T21, Im T22 -- the storage order the
   !!     numerical core writes and the expansion routine reads back.
   !!   - The original reported an out-of-range angle through WRITE and
   !!     STOP.  A library must do neither, so the condition returns
   !!     IERR = 1 with a zero amplitude matrix.  Three further status codes
   !!     replace what used to be undefined behavior: a workspace holding no
   !!     converged T-matrix, one whose T-matrix has since been replaced, and
   !!     a truncation order outside the stored block.
   !!   - NMAX and the wavelength are not arguments.  They are read from the
   !!     SCATMAT the evaluation returned and from the workspace itself, so
   !!     the multipole sum and the 1/DK scaling cannot be given values that
   !!     belong to some other particle.
   !!   - Shared DO termination labels are written as nested DO/END DO.
   !!   - The terminal reporting block of the original AMPL (the WRITE and
   !!     PRINT statements echoing the angles and the four amplitude-matrix
   !!     elements, with their FORMAT labels) is absent, so the routine is
   !!     silent inside a large node loop.
   !!
   !! Every expression, operation order, array rank and extent, and working
   !! precision is that of the original.
   !!
   !! Workspace contract: AMPL reads the T-matrix that the numerical core
   !! left in the workspace.  The same workspace must therefore have
   !! completed a successful scattering-matrix evaluation for this particle
   !! first, and SCATMAT must be the one that evaluation returned.  Both the
   !! truncation order and the wavelength come from there, and the generation
   !! tag the two carry catches a T-matrix replaced in between.  The expansion
   !! routine only reads that storage, so it survives the evaluation intact;
   !! nothing else in the library writes it.
   use tmatrix_kinds, only: wp
   use tmatrix_types, only: tmatrix_workspace_t, tmatrix_scatmat_t
   implicit real(wp) (a-h, o-z)
   private
   public :: ampl

   integer, parameter :: npn4 = 80, npn6 = npn4 + 1

   !! Fourth-dimension slots of WORK%FULL_TSTORE.
   integer, parameter :: slot_rt11 = 1, slot_rt12 = 2, slot_rt21 = 3, slot_rt22 = 4
   integer, parameter :: slot_it11 = 5, slot_it12 = 6, slot_it21 = 7, slot_it22 = 8

contains

!   CALCULATION OF THE AMPLITUDE MATRIX
!
!   AMPL returns the 2x2 amplitude matrix (VV, VH, HV, HH), each carrying
!   the dimension of length (the 1/DK, DK = 2 pi / DLAM, scaling is applied
!   just before the R/R1 rotations below).  TL, TL1 are the incidence and
!   scattering polar angles (deg); PL, PL1 the corresponding azimuths (deg);
!   ALPHA, BETA the Euler angles orienting the particle symmetry axis (deg).
!   SCATMAT is what tmatrix_eval_scattering_matrix returned for the particle
!   whose T-matrix WORK holds; the multipole truncation NMAX and the
!   wavelength DLAM are taken from it and from WORK, never from the caller.
!
!   IERR = 0 success
!          1 an angular parameter is outside its allowable range
!          2 the workspace holds no converged T-matrix
!          3 NMAX outside 1..NPN4, i.e. outside the stored T-matrix block
!          4 SCATMAT belongs to a T-matrix a later evaluation replaced

      SUBROUTINE AMPL (WORK,SCATMAT,TL,TL1,PL,PL1,ALPHA,BETA, &
&                      VV,VH,HV,HH,IERR)
      IMPLICIT real(wp) (A-B,D-H,O-Z), complex(wp) (C)
      TYPE(tmatrix_workspace_t), INTENT(IN) :: WORK
      TYPE(tmatrix_scatmat_t), INTENT(IN) :: SCATMAT
      real(wp), INTENT(IN) :: TL,TL1,PL,PL1,ALPHA,BETA
      complex(wp), INTENT(OUT) :: VV,VH,HV,HH
      INTEGER, INTENT(OUT) :: IERR
      INTEGER NMAX
      real(wp) DLAM
      real(wp) AL(3,2),AL1(3,2),AP(2,3),AP1(2,3),B(3,3), &
&              R(2,2),R1(2,2),C(3,2),CA,CB,CT,CP,CTP,CPP,CT1,CP1, &
&              CTP1,CPP1
      real(wp) DV1(NPN6),DV2(NPN6),DV01(NPN6),DV02(NPN6)
      complex(wp) CAL(NPN4,NPN4)

      IERR=0
      VV=(0D0,0D0)
      VH=(0D0,0D0)
      HV=(0D0,0D0)
      HH=(0D0,0D0)

      IF (.NOT.ALLOCATED(WORK%FULL_TSTORE).OR..NOT.WORK%TM_FILLED) THEN
          IERR=2
          RETURN
      ENDIF
      IF (SCATMAT%STATE_TAG.NE.WORK%TM_TAG) THEN
          IERR=4
          RETURN
      ENDIF
      NMAX=SCATMAT%NMAX_TM
      DLAM=WORK%TM_LAMBDA_UM
      IF (NMAX.LT.1.OR.NMAX.GT.NPN4) THEN
          IERR=3
          RETURN
      ENDIF
      IF (ALPHA.LT.0D0.OR.ALPHA.GT.360D0.OR. &
&         BETA.LT.0D0.OR.BETA.GT.180D0.OR. &
&         TL.LT.0D0.OR.TL.GT.180D0.OR. &
&         TL1.LT.0D0.OR.TL1.GT.180D0.OR. &
&         PL.LT.0D0.OR.PL.GT.360D0.OR. &
&         PL1.LT.0D0.OR.PL1.GT.360D0) THEN
          IERR=1
          RETURN
      ENDIF
      PIN=DACOS(-1D0)
      PIN2=PIN*0.5D0
      PI=PIN/180D0
      ALPH=ALPHA*PI
      BET=BETA*PI
      THETL=TL*PI
      PHIL=PL*PI
      THETL1=TL1*PI
      PHIL1=PL1*PI

      EPS=1D-7
      IF (THETL.LT.PIN2) THETL=THETL+EPS
      IF (THETL.GT.PIN2) THETL=THETL-EPS
      IF (THETL1.LT.PIN2) THETL1=THETL1+EPS
      IF (THETL1.GT.PIN2) THETL1=THETL1-EPS
      IF (PHIL.LT.PIN) PHIL=PHIL+EPS
      IF (PHIL.GT.PIN) PHIL=PHIL-EPS
      IF (PHIL1.LT.PIN) PHIL1=PHIL1+EPS
      IF (PHIL1.GT.PIN) PHIL1=PHIL1-EPS
      IF (BET.LE.PIN2.AND.PIN2-BET.LE.EPS) BET=BET-EPS
      IF (BET.GT.PIN2.AND.BET-PIN2.LE.EPS) BET=BET+EPS

!_____________COMPUTE THETP, PHIP, THETP1, AND PHIP1, EQS. (8), (19), AND (20)

      CB=DCOS(BET)
      SB=DSIN(BET)
      CT=DCOS(THETL)
      ST=DSIN(THETL)
      CP=DCOS(PHIL-ALPH)
      SP=DSIN(PHIL-ALPH)
      CTP=CT*CB+ST*SB*CP
      THETP=DACOS(CTP)
      CPP=CB*ST*CP-SB*CT
      SPP=ST*SP
      PHIP=DATAN(SPP/CPP)
      IF (PHIP.GT.0D0.AND.SP.LT.0D0) PHIP=PHIP+PIN
      IF (PHIP.LT.0D0.AND.SP.GT.0D0) PHIP=PHIP+PIN
      IF (PHIP.LT.0D0) PHIP=PHIP+2D0*PIN

      CT1=DCOS(THETL1)
      ST1=DSIN(THETL1)
      CP1=DCOS(PHIL1-ALPH)
      SP1=DSIN(PHIL1-ALPH)
      CTP1=CT1*CB+ST1*SB*CP1
      THETP1=DACOS(CTP1)
      CPP1=CB*ST1*CP1-SB*CT1
      SPP1=ST1*SP1
      PHIP1=DATAN(SPP1/CPP1)
      IF (PHIP1.GT.0D0.AND.SP1.LT.0D0) PHIP1=PHIP1+PIN
      IF (PHIP1.LT.0D0.AND.SP1.GT.0D0) PHIP1=PHIP1+PIN
      IF (PHIP1.LT.0D0) PHIP1=PHIP1+2D0*PIN

!____________COMPUTE MATRIX BETA, EQ. (21)

      CA=DCOS(ALPH)
      SA=DSIN(ALPH)
      B(1,1)=CA*CB
      B(1,2)=SA*CB
      B(1,3)=-SB
      B(2,1)=-SA
      B(2,2)=CA
      B(2,3)=0D0
      B(3,1)=CA*SB
      B(3,2)=SA*SB
      B(3,3)=CB

!____________COMPUTE MATRICES AL AND AL1, EQ. (14)

      CP=DCOS(PHIL)
      SP=DSIN(PHIL)
      CP1=DCOS(PHIL1)
      SP1=DSIN(PHIL1)
      AL(1,1)=CT*CP
      AL(1,2)=-SP
      AL(2,1)=CT*SP
      AL(2,2)=CP
      AL(3,1)=-ST
      AL(3,2)=0D0
      AL1(1,1)=CT1*CP1
      AL1(1,2)=-SP1
      AL1(2,1)=CT1*SP1
      AL1(2,2)=CP1
      AL1(3,1)=-ST1
      AL1(3,2)=0D0

!____________COMPUTE MATRICES AP^(-1) AND AP1^(-1), EQ. (15)

      CT=CTP
      ST=DSIN(THETP)
      CP=DCOS(PHIP)
      SP=DSIN(PHIP)
      CT1=CTP1
      ST1=DSIN(THETP1)
      CP1=DCOS(PHIP1)
      SP1=DSIN(PHIP1)
      AP(1,1)=CT*CP
      AP(1,2)=CT*SP
      AP(1,3)=-ST
      AP(2,1)=-SP
      AP(2,2)=CP
      AP(2,3)=0D0
      AP1(1,1)=CT1*CP1
      AP1(1,2)=CT1*SP1
      AP1(1,3)=-ST1
      AP1(2,1)=-SP1
      AP1(2,2)=CP1
      AP1(2,3)=0D0

!____________COMPUTE MATRICES R AND R^(-1), EQ. (13)
      DO I=1,3
         DO J=1,2
            X=0D0
            DO K=1,3
               X=X+B(I,K)*AL(K,J)
            ENDDO
            C(I,J)=X
         ENDDO
      ENDDO
      DO I=1,2
         DO J=1,2
            X=0D0
            DO K=1,3
               X=X+AP(I,K)*C(K,J)
            ENDDO
            R(I,J)=X
         ENDDO
      ENDDO
      DO I=1,3
         DO J=1,2
            X=0D0
            DO K=1,3
               X=X+B(I,K)*AL1(K,J)
            ENDDO
            C(I,J)=X
         ENDDO
      ENDDO
      DO I=1,2
         DO J=1,2
            X=0D0
            DO K=1,3
               X=X+AP1(I,K)*C(K,J)
            ENDDO
            R1(I,J)=X
         ENDDO
      ENDDO
      D=1D0/(R1(1,1)*R1(2,2)-R1(1,2)*R1(2,1))
      X=R1(1,1)
      R1(1,1)=R1(2,2)*D
      R1(1,2)=-R1(1,2)*D
      R1(2,1)=-R1(2,1)*D
      R1(2,2)=X*D

      CI=(0D0,1D0)
      DO NN=1,NMAX
         DO N=1,NMAX
            CN=CI**(NN-N-1)
            DNN=DFLOAT((2*N+1)*(2*NN+1))
            DNN=DNN/DFLOAT( N*NN*(N+1)*(NN+1) )
            RN=DSQRT(DNN)
            CAL(N,NN)=CN*RN
         ENDDO
      ENDDO
      DCTH0=CTP
      DCTH=CTP1
      PH=PHIP1-PHIP
      VV=(0D0,0D0)
      VH=(0D0,0D0)
      HV=(0D0,0D0)
      HH=(0D0,0D0)
      DO M=0,NMAX
         M1=M+1
         NMIN=MAX(M,1)
         CALL VIGAMPL (DCTH, NMAX, M, DV1, DV2)
         CALL VIGAMPL (DCTH0, NMAX, M, DV01, DV02)
         FC=2D0*DCOS(M*PH)
         FS=2D0*DSIN(M*PH)
         DO NN=NMIN,NMAX
            DV1NN=M*DV01(NN)
            DV2NN=DV02(NN)
            DO N=NMIN,NMAX
               DV1N=M*DV1(N)
               DV2N=DV2(N)

               CT11=DCMPLX(WORK%FULL_TSTORE(M1,N,NN,SLOT_RT11), &
&                          WORK%FULL_TSTORE(M1,N,NN,SLOT_IT11))
               CT22=DCMPLX(WORK%FULL_TSTORE(M1,N,NN,SLOT_RT22), &
&                          WORK%FULL_TSTORE(M1,N,NN,SLOT_IT22))

               IF (M.EQ.0) THEN

                  CN=CAL(N,NN)*DV2N*DV2NN

                  VV=VV+CN*CT22
                  HH=HH+CN*CT11

                 ELSE

                  CT12=DCMPLX(WORK%FULL_TSTORE(M1,N,NN,SLOT_RT12), &
&                             WORK%FULL_TSTORE(M1,N,NN,SLOT_IT12))
                  CT21=DCMPLX(WORK%FULL_TSTORE(M1,N,NN,SLOT_RT21), &
&                             WORK%FULL_TSTORE(M1,N,NN,SLOT_IT21))

                  CN1=CAL(N,NN)*FC
                  CN2=CAL(N,NN)*FS

                  D11=DV1N*DV1NN
                  D12=DV1N*DV2NN
                  D21=DV2N*DV1NN
                  D22=DV2N*DV2NN

                  VV=VV+(CT11*D11+CT21*D21 &
&                       +CT12*D12+CT22*D22)*CN1

                  VH=VH+(CT11*D12+CT21*D22 &
&                       +CT12*D11+CT22*D21)*CN2

                  HV=HV-(CT11*D21+CT21*D11 &
&                       +CT12*D22+CT22*D12)*CN2

                  HH=HH+(CT11*D22+CT21*D12 &
&                       +CT12*D21+CT22*D11)*CN1
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      DK=2D0*PIN/DLAM
      VV=VV/DK
      VH=VH/DK
      HV=HV/DK
      HH=HH/DK
      CVV=VV*R(1,1)+VH*R(2,1)
      CVH=VV*R(1,2)+VH*R(2,2)
      CHV=HV*R(1,1)+HH*R(2,1)
      CHH=HV*R(1,2)+HH*R(2,2)
      VV=R1(1,1)*CVV+R1(1,2)*CHV
      VH=R1(1,1)*CVH+R1(1,2)*CHH
      HV=R1(2,1)*CVV+R1(2,2)*CHV
      HH=R1(2,1)*CVH+R1(2,2)*CHH

      RETURN
      END

!*****************************************************************
!
!     Calculation of the functions
!     DV1(N)=dvig(0,m,n,arccos x)/sin(arccos x)
!     and
!     DV2(N)=[d/d(arccos x)] dvig(0,m,n,arccos x)
!     1.LE.N.LE.NMAX
!     0.LE.X.LE.1

      SUBROUTINE VIGAMPL (X, NMAX, M, DV1, DV2)
      IMPLICIT real(wp) (A-H,O-Z)
      real(wp), INTENT(IN) :: X
      INTEGER, INTENT(IN) :: NMAX, M
      real(wp), INTENT(OUT) :: DV1(NPN6), DV2(NPN6)
      DO N=1,NMAX
         DV1(N)=0D0
         DV2(N)=0D0
      ENDDO
      DX=DABS(X)
      IF (DABS(1D0-DX).LE.1D-10) GO TO 100
      A=1D0
      QS=DSQRT(1D0-X*X)
      QS1=1D0/QS
      DSI=QS1
      IF (M.NE.0) GO TO 20
      D1=1D0
      D2=X
      DO N=1,NMAX
         QN=DFLOAT(N)
         QN1=DFLOAT(N+1)
         QN2=DFLOAT(2*N+1)
         D3=(QN2*X*D2-QN*D1)/QN1
         DER=QS1*(QN1*QN/QN2)*(-D1+D3)
         DV1(N)=D2*DSI
         DV2(N)=DER
         D1=D2
         D2=D3
      ENDDO
      RETURN
   20 QMM=DFLOAT(M*M)
      DO I=1,M
         I2=I*2
         A=A*DSQRT(DFLOAT(I2-1)/DFLOAT(I2))*QS
      ENDDO
      D1=0D0
      D2=A
      DO N=M,NMAX
         QN=DFLOAT(N)
         QN2=DFLOAT(2*N+1)
         QN1=DFLOAT(N+1)
         QNM=DSQRT(QN*QN-QMM)
         QNM1=DSQRT(QN1*QN1-QMM)
         D3=(QN2*X*D2-QNM*D1)/QNM1
         DER=QS1*(-QN1*QNM*D1+QN*QNM1*D3)/QN2
         DV1(N)=D2*DSI
         DV2(N)=DER
         D1=D2
         D2=D3
      ENDDO
      RETURN
  100 IF (M.NE.1) RETURN
      DO N=1,NMAX
         DN=DFLOAT(N*(N+1))
         DN=0.5D0*DSQRT(DN)
         IF (X.LT.0D0) DN=DN*(-1)**(N+1)
         DV1(N)=DN
         IF (X.LT.0D0) DN=-DN
         DV2(N)=DN
      ENDDO
      RETURN
      END

end module fixed_orientation_amplitude
