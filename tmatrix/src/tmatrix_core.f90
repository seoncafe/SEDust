module tmatrix_core
   !! Free-form, caller-owned direct-storage T-matrix numerical core.
   !!
   !! This is a mechanical fixed-to-free-form migration of
   !! tmd_one_full_direct.f.  Formulae, storage ranks, loop order, and
   !! working precision are intentionally unchanged.  The retained
   !! implicit mapping matches the original F77 source; a later,
   !! separately-regressed pass may make every local declaration explicit.
   use, intrinsic :: iso_fortran_env, only: real32
   use tmatrix_kinds, only: wp
   use tmatrix_types, only: tmatrix_workspace_t
   implicit real(wp) (a-h, o-z)
   private
   public :: tmd_one_full_direct
   !! Gauss-Legendre nodes and weights of the angular quadrature this core
   !! integrates over.  Exposed so a caller that has to integrate over the
   !! same angular grid uses the identical nodes, bit for bit, instead of
   !! carrying a second copy of the rule.  tmatrix_api re-exports it as GAUSS,
   !! Mishchenko's original name for it.
   public :: gauss_full_direct

   integer, parameter :: npn1 = 100, npng1 = 300, npng2 = 2*npng1
   integer, parameter :: npn2 = 2*npn1, npl = npn2 + 1, npn3 = npn1 + 1
   integer, parameter :: npn4 = 80, npn5 = 2*npn4, npn6 = npn4 + 1, npl1 = npn5 + 1

   interface
      subroutine zgetrf(m, n, a, lda, ipiv, info)
         import :: wp
         integer, intent(in) :: m, n, lda
         complex(wp), intent(inout) :: a(lda, *)
         integer, intent(out) :: ipiv(*), info
      end subroutine zgetrf
      subroutine zgetri(n, a, lda, ipiv, work, lwork, info)
         import :: wp
         integer, intent(in) :: n, lda, lwork
         complex(wp), intent(inout) :: a(lda, *), work(*)
         integer, intent(in) :: ipiv(*)
         integer, intent(out) :: info
      end subroutine zgetri
   end interface

contains
!***********************************************************************
!  FULL direct workspace backend.
!  All mutable workspace is supplied by the caller.
!  TMAT, FAC, and SS are explicit caller-owned workspace arrays.
!
!  WS carries the largest scratch arrays of TT, TMATR0/TMATR, and GSP.  It is
!  named WS rather than WORK because TT_FULL_DIRECT already declares a local
!  array WORK, and Fortran names are case insensitive.
!***********************************************************************
      SUBROUTINE TMD_ONE_FULL_DIRECT(AXI, LAM, MRR, MRI, EPS, NP, &
&                         DDELT, NDGS, &
&                         QEXT, QSCA, WALB, ASYMM, IERR, LAPACK_INFO, &
&                         TR1,TI1,QR,QI,RGQR,RGQI, &
&                         J,Y,JR,JI,DJ,DY,DJR,DJI,B1R,B1I,B2R,B2I, &
&                         TSTORE,TMATR_WORK,D1G,D2G,D3G,D4G,D5RG,D5IG, &
&                         FAC,SSIGN,WS, &
&                         AL1_OUT,AL2_OUT,AL3_OUT,AL4_OUT,BE1_OUT,BE2_OUT, &
&                         LMAX_OUT,NMAX_TM_OUT)
      IMPLICIT real(wp) (A-H,O-Z)
      TYPE(tmatrix_workspace_t), INTENT(INOUT) :: WS
      real(wp)  AXI, LAM, MRR, MRI, EPS, DDELT
      real(wp)  QEXT, QSCA, WALB, ASYMM
      INTEGER NP, NDGS, IERR, LAPACK_INFO
      real(wp)  X(NPNG2),W(NPNG2),S(NPNG2),SS(NPNG2), &
&              AN(NPN1),R(NPNG2),DR(NPNG2), &
&              DDR(NPNG2),DRR(NPNG2),DRI(NPNG2),ANN(NPN1,NPN1)
      real(wp)  TR1(NPN2,NPN2),TI1(NPN2,NPN2), &
&              QR(NPN2,NPN2),QI(NPN2,NPN2),RGQR(NPN2,NPN2), &
&              RGQI(NPN2,NPN2), &
&              J(NPNG2,NPN1),Y(NPNG2,NPN1),JR(NPNG2,NPN1), &
&              JI(NPNG2,NPN1),DJ(NPNG2,NPN1),DY(NPNG2,NPN1), &
&              DJR(NPNG2,NPN1),DJI(NPNG2,NPN1), &
&              AL1(NPL),AL2(NPL),AL3(NPL),AL4(NPL),BE1(NPL),BE2(NPL)
      real(real32) B1R(NPL1,NPL1,NPN4),B1I(NPL1,NPL1,NPN4), &
&             B2R(NPL1,NPL1,NPN4),B2I(NPL1,NPL1,NPN4)
      real(real32) TSTORE(NPN6,NPN4,NPN4,8), &
&             D1G(NPL1,NPN4,NPN4),D2G(NPL1,NPN4,NPN4), &
&             D3G(NPL1,NPN4,NPN4),D4G(NPL1,NPN4,NPN4), &
&             D5RG(NPL1,NPN4,NPN4),D5IG(NPL1,NPN4,NPN4)
      real(wp) TMATR_WORK(NPN1,NPN1,16),FAC(900),SSIGN(900)
!  Optional random-orientation scattering-matrix expansion.  AL1..BE2 and
!  LMAX are already computed by the GSP call below; NMAX_TM is the
!  multipole truncation NMAX1 to which TSTORE is filled.  Every numerical
!  statement of this routine is independent of these arguments, so a
!  caller that omits them gets exactly the cross-section-only evaluation.
      real(wp), INTENT(OUT), OPTIONAL :: AL1_OUT(NPL),AL2_OUT(NPL), &
&              AL3_OUT(NPL),AL4_OUT(NPL),BE1_OUT(NPL),BE2_OUT(NPL)
      INTEGER, INTENT(OUT), OPTIONAL :: LMAX_OUT, NMAX_TM_OUT

      P    = DACOS(-1D0)
      IERR = 0
      LAPACK_INFO = 0
      QEXT  = 0D0
      QSCA  = 0D0
      WALB  = 0D0
      ASYMM = 0D0

!  Define the optional outputs up front so that an early error return
!  leaves a well-defined (all-zero, LMAX = 0) expansion rather than
!  uninitialized storage.
      IF (PRESENT(AL1_OUT)) AL1_OUT = 0D0
      IF (PRESENT(AL2_OUT)) AL2_OUT = 0D0
      IF (PRESENT(AL3_OUT)) AL3_OUT = 0D0
      IF (PRESENT(AL4_OUT)) AL4_OUT = 0D0
      IF (PRESENT(BE1_OUT)) BE1_OUT = 0D0
      IF (PRESENT(BE2_OUT)) BE2_OUT = 0D0
      IF (PRESENT(LMAX_OUT)) LMAX_OUT = 0
      IF (PRESENT(NMAX_TM_OUT)) NMAX_TM_OUT = 0

!  NCHECK setup (unchanged from original main program).
      NCHECK = 0
      IF (NP.EQ.-1.OR.NP.EQ.-2) NCHECK = 1
      IF (NP.GT.0.AND.(-1)**NP.EQ.1) NCHECK = 1

!  Original main program did `DDELT = 0.1D0*DDELT` directly on the input
!  variable. As a wrapper that is now an argument, do this on a local
!  copy DDELT_LOC instead — caller may pass a Fortran PARAMETER (read-only
!  storage), in which case writing to DDELT segfaults with SIGBUS.
      DDELT_LOC = 0.1D0*DDELT

!  Single particle, RAT = 1 (volume-equivalent sphere): A = AXI.
      A    = AXI
      XEV  = 2D0*P*A/LAM
      IXXX = XEV + 4.05D0*XEV**0.333333D0
      INM1 = MAX0(4,IXXX)
      IF (INM1.GE.NPN1) THEN
         IERR = 1
         RETURN
      ENDIF

      QEXT1 = 0D0
      QSCA1 = 0D0
      DO 50 NMA = INM1, NPN1
         NMAX   = NMA
         MMAX   = 1
         NGAUSS = NMAX*NDGS
         IF (NGAUSS.GT.NPNG1) THEN
            IERR = 2
            RETURN
         ENDIF
         CALL CONST_FULL_DIRECT(NGAUSS,NMAX,MMAX,P,X,W,AN,ANN,S,SS,NP,EPS)
         CALL VARY_FULL_DIRECT(LAM,MRR,MRI,A,EPS,NP,NGAUSS,X,P,PPI,PIR,PII,R, &
&                   DR,DDR,DRR,DRI,NMAX,IERR,J,Y,JR,JI,DJ,DY,DJR,DJI)
         IF (IERR.NE.0) RETURN
         CALL TMATR0_FULL_DIRECT(NGAUSS,X,W,AN,ANN,S,SS,PPI,PIR,PII,R,DR, &
&                     DDR,DRR,DRI,NMAX,NCHECK,IERR,LAPACK_INFO, &
&                     TR1,TI1,QR,QI,RGQR,RGQI,J,Y,JR,JI,DJ,DY,DJR,DJI, &
&                     TMATR_WORK,WS)
         IF (IERR.NE.0) RETURN
         QEXT = 0D0
         QSCA = 0D0
         DO N = 1, NMAX
            N1     = N + NMAX
            TR1NN  = TR1(N,N)
            TI1NN  = TI1(N,N)
            TR1NN1 = TR1(N1,N1)
            TI1NN1 = TI1(N1,N1)
            DN1    = DFLOAT(2*N+1)
            QSCA   = QSCA + DN1*(TR1NN*TR1NN + TI1NN*TI1NN &
&                                + TR1NN1*TR1NN1 + TI1NN1*TI1NN1)
            QEXT   = QEXT + (TR1NN + TR1NN1)*DN1
         ENDDO
         DSCA  = DABS((QSCA1-QSCA)/QSCA)
         DEXT  = DABS((QEXT1-QEXT)/QEXT)
         QEXT1 = QEXT
         QSCA1 = QSCA
         NMIN  = DFLOAT(NMAX)/2D0 + 1D0
         DO 10 N = NMIN, NMAX
            N1     = N + NMAX
            TR1NN  = TR1(N,N)
            TI1NN  = TI1(N,N)
            TR1NN1 = TR1(N1,N1)
            TI1NN1 = TI1(N1,N1)
            DN1    = DFLOAT(2*N+1)
            DQSCA  = DN1*(TR1NN*TR1NN + TI1NN*TI1NN &
&                         + TR1NN1*TR1NN1 + TI1NN1*TI1NN1)
            DQEXT  = (TR1NN + TR1NN1)*DN1
            DQSCA  = DABS(DQSCA/QSCA)
            DQEXT  = DABS(DQEXT/QEXT)
            NMAX1  = N
            IF (DQSCA.LE.DDELT_LOC.AND.DQEXT.LE.DDELT_LOC) GO TO 12
   10    CONTINUE
   12    CONTINUE
         IF (DSCA.LE.DDELT_LOC.AND.DEXT.LE.DDELT_LOC) GO TO 55
         IF (NMA.EQ.NPN1) THEN
            IERR = 3
            RETURN
         ENDIF
   50 CONTINUE

   55 NNNGGG = NGAUSS + 1
      MMAX   = NMAX1
      DO 150 NGAUS = NNNGGG, NPNG1
         NGAUSS = NGAUS
         CALL CONST_FULL_DIRECT(NGAUSS,NMAX,MMAX,P,X,W,AN,ANN,S,SS,NP,EPS)
         CALL VARY_FULL_DIRECT(LAM,MRR,MRI,A,EPS,NP,NGAUSS,X,P,PPI,PIR,PII,R, &
&                   DR,DDR,DRR,DRI,NMAX,IERR,J,Y,JR,JI,DJ,DY,DJR,DJI)
         IF (IERR.NE.0) RETURN
         CALL TMATR0_FULL_DIRECT(NGAUSS,X,W,AN,ANN,S,SS,PPI,PIR,PII,R,DR, &
&                     DDR,DRR,DRI,NMAX,NCHECK,IERR,LAPACK_INFO, &
&                     TR1,TI1,QR,QI,RGQR,RGQI,J,Y,JR,JI,DJ,DY,DJR,DJI, &
&                     TMATR_WORK,WS)
         IF (IERR.NE.0) RETURN
         QEXT = 0D0
         QSCA = 0D0
         DO 104 N = 1, NMAX
            N1     = N + NMAX
            TR1NN  = TR1(N,N)
            TI1NN  = TI1(N,N)
            TR1NN1 = TR1(N1,N1)
            TI1NN1 = TI1(N1,N1)
            DN1    = DFLOAT(2*N+1)
            QSCA   = QSCA + DN1*(TR1NN*TR1NN + TI1NN*TI1NN &
&                                + TR1NN1*TR1NN1 + TI1NN1*TI1NN1)
            QEXT   = QEXT + (TR1NN + TR1NN1)*DN1
  104    CONTINUE
         DSCA  = DABS((QSCA1-QSCA)/QSCA)
         DEXT  = DABS((QEXT1-QEXT)/QEXT)
         QEXT1 = QEXT
         QSCA1 = QSCA
         IF (DSCA.LE.DDELT_LOC.AND.DEXT.LE.DDELT_LOC) GO TO 155
  150 CONTINUE
      IERR = 5
  155 CONTINUE
      QSCA = 0D0
      QEXT = 0D0
      NNM  = NMAX*2
      DO 204 N = 1, NNM
         QEXT = QEXT + TR1(N,N)
  204 CONTINUE
      IF (NMAX1.GT.NPN4) THEN
         IERR = 4
         RETURN
      ENDIF
      DO 213 N2 = 1, NMAX1
         NN2 = N2 + NMAX
         DO 213 N1 = 1, NMAX1
            NN1 = N1 + NMAX
            ZZ1 = TR1(N1,N2)
            TSTORE(1,N1,N2,1) = ZZ1
            ZZ2 = TI1(N1,N2)
            TSTORE(1,N1,N2,5) = ZZ2
            ZZ3 = TR1(N1,NN2)
            TSTORE(1,N1,N2,2) = ZZ3
            ZZ4 = TI1(N1,NN2)
            TSTORE(1,N1,N2,6) = ZZ4
            ZZ5 = TR1(NN1,N2)
            TSTORE(1,N1,N2,3) = ZZ5
            ZZ6 = TI1(NN1,N2)
            TSTORE(1,N1,N2,7) = ZZ6
            ZZ7 = TR1(NN1,NN2)
            TSTORE(1,N1,N2,4) = ZZ7
            ZZ8 = TI1(NN1,NN2)
            TSTORE(1,N1,N2,8) = ZZ8
            QSCA = QSCA + ZZ1*ZZ1 + ZZ2*ZZ2 + ZZ3*ZZ3 + ZZ4*ZZ4 &
&                        + ZZ5*ZZ5 + ZZ6*ZZ6 + ZZ7*ZZ7 + ZZ8*ZZ8
  213 CONTINUE
      DO 220 M = 1, NMAX1
         CALL TMATR_FULL_DIRECT(M,NGAUSS,X,W,AN,ANN,S,SS,PPI,PIR,PII,R,DR, &
&                    DDR,DRR,DRI,NMAX,NCHECK,IERR,LAPACK_INFO, &
&                    TR1,TI1,QR,QI,RGQR,RGQI,J,Y,JR,JI,DJ,DY,DJR,DJI, &
&                    TMATR_WORK,WS)
         IF (IERR.NE.0) RETURN
         NM  = NMAX  - M + 1
         NM1 = NMAX1 - M + 1
         M1  = M + 1
         QSC = 0D0
         DO 214 N2 = 1, NM1
            NN2 = N2 + M  - 1
            N22 = N2 + NM
            DO 214 N1 = 1, NM1
               NN1 = N1 + M  - 1
               N11 = N1 + NM
               ZZ1 = TR1(N1,N2)
               TSTORE(M1,NN1,NN2,1) = ZZ1
               ZZ2 = TI1(N1,N2)
               TSTORE(M1,NN1,NN2,5) = ZZ2
               ZZ3 = TR1(N1,N22)
               TSTORE(M1,NN1,NN2,2) = ZZ3
               ZZ4 = TI1(N1,N22)
               TSTORE(M1,NN1,NN2,6) = ZZ4
               ZZ5 = TR1(N11,N2)
               TSTORE(M1,NN1,NN2,3) = ZZ5
               ZZ6 = TI1(N11,N2)
               TSTORE(M1,NN1,NN2,7) = ZZ6
               ZZ7 = TR1(N11,N22)
               TSTORE(M1,NN1,NN2,4) = ZZ7
               ZZ8 = TI1(N11,N22)
               TSTORE(M1,NN1,NN2,8) = ZZ8
               QSC = QSC + (ZZ1*ZZ1 + ZZ2*ZZ2 + ZZ3*ZZ3 + ZZ4*ZZ4 &
&                          + ZZ5*ZZ5 + ZZ6*ZZ6 + ZZ7*ZZ7 + ZZ8*ZZ8)*2D0
  214    CONTINUE
         NNM = 2*NM
         QXT = 0D0
         DO 215 N = 1, NNM
            QXT = QXT + TR1(N,N)*2D0
  215    CONTINUE
         QSCA = QSCA + QSC
         QEXT = QEXT + QXT
  220 CONTINUE
      COEFF1 = LAM*LAM*0.5D0/P
      CSCA = QSCA*COEFF1
      CEXT = -QEXT*COEFF1
      CALL GSP_FULL_DIRECT(NMAX1,CSCA,LAM,AL1,AL2,AL3,AL4,BE1,BE2,LMAX, &
&                          IERR,B1R,B1I,B2R,B2I,TSTORE,D1G,D2G,D3G,D4G, &
&                          D5RG,D5IG,FAC,SSIGN,WS)
      IF (IERR.NE.0) RETURN

!  Single particle: ALPH/BET reduce to AL/BE without quadrature
!  weighting, hence WALB = CSCA/CEXT and ASYMM = AL1(2)/3.
      WALB  = CSCA/CEXT
      ASYMM = AL1(2)/3D0

!  Q convention: Q = C / (pi * AXI**2). AXI is in microns; CSCA, CEXT
!  carry the LAM**2 factor from COEFF1 above and inherit the units of
!  LAM**2 (microns**2), so Q is dimensionless.
      QEXT = CEXT/(P*AXI*AXI)
      QSCA = CSCA/(P*AXI*AXI)

!  GSP writes elements 1..LMAX+1 only, so the zeroed tail set on entry is
!  what the caller sees beyond the truncation order.
      IF (PRESENT(AL1_OUT)) AL1_OUT(1:LMAX+1) = AL1(1:LMAX+1)
      IF (PRESENT(AL2_OUT)) AL2_OUT(1:LMAX+1) = AL2(1:LMAX+1)
      IF (PRESENT(AL3_OUT)) AL3_OUT(1:LMAX+1) = AL3(1:LMAX+1)
      IF (PRESENT(AL4_OUT)) AL4_OUT(1:LMAX+1) = AL4(1:LMAX+1)
      IF (PRESENT(BE1_OUT)) BE1_OUT(1:LMAX+1) = BE1(1:LMAX+1)
      IF (PRESENT(BE2_OUT)) BE2_OUT(1:LMAX+1) = BE2(1:LMAX+1)
      IF (PRESENT(LMAX_OUT)) LMAX_OUT = LMAX
!  NMAX1 is the multipole truncation to which TSTORE has just been filled:
!  the DO 213 and DO 220 loops above run to N = NMAX1 and M = NMAX1.  It is
!  the NMAX a fixed-orientation amplitude evaluation must iterate over.
      IF (PRESENT(NMAX_TM_OUT)) NMAX_TM_OUT = NMAX1

      RETURN
      END

      SUBROUTINE VARY_FULL_DIRECT (LAM,MRR,MRI,A,EPS,NP,NGAUSS,X,P,PPI,PIR,PII, &
&                       R,DR,DDR,DRR,DRI,NMAX,IERR,J,Y,JR,JI,DJ,DY,DJR,DJI)
      IMPLICIT real(wp) (A-H,O-Z)
      real(wp)  X(NPNG2),R(NPNG2),DR(NPNG2),MRR,MRI,LAM, &
&              Z(NPNG2),ZR(NPNG2),ZI(NPNG2), &
&              J(NPNG2,NPN1),Y(NPNG2,NPN1),JR(NPNG2,NPN1), &
&              JI(NPNG2,NPN1),DJ(NPNG2,NPN1), &
&              DJR(NPNG2,NPN1),DJI(NPNG2,NPN1),DDR(NPNG2), &
&              DRR(NPNG2),DRI(NPNG2), &
&              DY(NPNG2,NPN1)
      IERR=0
      NG=NGAUSS*2
      IF (NP.EQ.-1) CALL RSP1_FULL_DIRECT(X,NG,NGAUSS,A,EPS,NP,R,DR)
      IF (NP.GE.0) CALL RSP2_FULL_DIRECT(X,NG,A,EPS,NP,R,DR)
      IF (NP.EQ.-2) CALL RSP3_FULL_DIRECT(X,NG,NGAUSS,A,EPS,R,DR)
      PI=P*2D0/LAM
      PPI=PI*PI
      PIR=PPI*MRR
      PII=PPI*MRI
      V=1D0/(MRR*MRR+MRI*MRI)
      PRR=MRR*V
      PRI=-MRI*V
      TA=0D0
      DO 10 I=1,NG
           VV=DSQRT(R(I))
           V=VV*PI
           TA=MAX(TA,V)
           VV=1D0/V
           DDR(I)=VV
           DRR(I)=PRR*VV
           DRI(I)=PRI*VV
           V1=V*MRR
           V2=V*MRI
           Z(I)=V
           ZR(I)=V1
           ZI(I)=V2
   10 CONTINUE
      IF (NMAX.GT.NPN1) THEN
         IERR=6
         RETURN
      ENDIF
      TB=TA*DSQRT(MRR*MRR+MRI*MRI)
      TB=DMAX1(TB,DFLOAT(NMAX))
      NNMAX1=1.2D0*DSQRT(DMAX1(TA,DFLOAT(NMAX)))+3D0
      NNMAX2=(TB+4D0*(TB**0.33333D0)+1.2D0*DSQRT(TB))
      NNMAX2=NNMAX2-NMAX+5
      CALL BESS_FULL_DIRECT(Z,ZR,ZI,NG,NMAX,NNMAX1,NNMAX2,J,Y,JR,JI,DJ,DY,DJR,DJI)
      RETURN
      END

!**********************************************************************

      SUBROUTINE BESS_FULL_DIRECT (X,XR,XI,NG,NMAX,NNMAX1,NNMAX2,J,Y,JR,JI,DJ,DY,DJR,DJI)
      IMPLICIT real(wp) (A-H,O-Z)
      real(wp) X(NG),XR(NG),XI(NG), &
&              J(NPNG2,NPN1),Y(NPNG2,NPN1),JR(NPNG2,NPN1), &
&              JI(NPNG2,NPN1),DJ(NPNG2,NPN1),DY(NPNG2,NPN1), &
&              DJR(NPNG2,NPN1),DJI(NPNG2,NPN1), &
&              AJ(NPN1),AY(NPN1),AJR(NPN1),AJI(NPN1), &
&              ADJ(NPN1),ADY(NPN1),ADJR(NPN1), &
&              ADJI(NPN1)

      DO 10 I=1,NG
           XX=X(I)
           CALL RJB_FULL_DIRECT(XX,AJ,ADJ,NMAX,NNMAX1)
           CALL RYB_FULL_DIRECT(XX,AY,ADY,NMAX)
           YR=XR(I)
           YI=XI(I)
           CALL CJB_FULL_DIRECT(YR,YI,AJR,AJI,ADJR,ADJI,NMAX,NNMAX2)
           DO 10 N=1,NMAX
                J(I,N)=AJ(N)
                Y(I,N)=AY(N)
                JR(I,N)=AJR(N)
                JI(I,N)=AJI(N)
                DJ(I,N)=ADJ(N)
                DY(I,N)=ADY(N)
                DJR(I,N)=ADJR(N)
                DJI(I,N)=ADJI(N)
   10 CONTINUE
      RETURN
      END

!**********************************************************************


!  Inner subroutines below copied verbatim from tmd.lp.f.
!***********************************************************************

      SUBROUTINE TMATR0_FULL_DIRECT (NGAUSS,X,W,AN,ANN,S,SS,PPI,PIR,PII,R,DR,DDR, &
&                        DRR,DRI,NMAX,NCHECK,IERR,LAPACK_INFO,TR1,TI1,QR,QI,RGQR,RGQI, &
&                        J,Y,JR,JI,DJ,DY,DJR,DJI,TMATR_WORK,WS)
      IMPLICIT real(wp) (A-H,O-Z)
!  D1 and D2 are WS%TMATR_D1 and WS%TMATR_D2.  TMATR0 and TMATR are never
!  active at the same time, so they share the one pair.
      TYPE(tmatrix_workspace_t), INTENT(INOUT) :: WS
      real(wp)  X(NPNG2),W(NPNG2),AN(NPN1),S(NPNG2),SS(NPNG2), &
&              R(NPNG2),DR(NPNG2),SIG(NPN2), &
&              J(NPNG2,NPN1),Y(NPNG2,NPN1), &
&              JR(NPNG2,NPN1),JI(NPNG2,NPN1),DJ(NPNG2,NPN1), &
&              DY(NPNG2,NPN1),DJR(NPNG2,NPN1), &
&              DJI(NPNG2,NPN1),DDR(NPNG2),DRR(NPNG2), &
&              DRI(NPNG2),DS(NPNG2),DSS(NPNG2),RR(NPNG2), &
&              DV1(NPN1),DV2(NPN1)

      real(wp)  ANN(NPN1,NPN1),TMATR_WORK(NPN1,NPN1,16)
      real(wp) TR1(NPN2,NPN2),TI1(NPN2,NPN2), &
&             QR(NPN2,NPN2),QI(NPN2,NPN2), &
&             RGQR(NPN2,NPN2),RGQI(NPN2,NPN2)
      MM1=1
      NNMAX=NMAX+NMAX
      NG=2*NGAUSS
      NGSS=NG
      FACTOR=1D0
      IF (NCHECK.EQ.1) THEN
            NGSS=NGAUSS
            FACTOR=2D0
         ELSE
            CONTINUE
      ENDIF
      SI=1D0
      DO 5 N=1,NNMAX
           SI=-SI
           SIG(N)=SI
    5 CONTINUE
   20 DO 25 I=1,NGAUSS
         I1=NGAUSS+I
         I2=NGAUSS-I+1
         CALL VIG_FULL_DIRECT( X(I1), NMAX, 0, DV1, DV2)
         DO 25 N=1,NMAX
            SI=SIG(N)
            DD1=DV1(N)
            DD2=DV2(N)
            WS%TMATR_D1(I1,N)=DD1
            WS%TMATR_D2(I1,N)=DD2
            WS%TMATR_D1(I2,N)=DD1*SI
            WS%TMATR_D2(I2,N)=-DD2*SI
   25 CONTINUE
   30 DO 40 I=1,NGSS
           RR(I)=W(I)*R(I)
   40 CONTINUE

      DO 300  N1=MM1,NMAX
           AN1=AN(N1)
           DO 300 N2=MM1,NMAX
                AN2=AN(N2)
                AR12=0D0
                AR21=0D0
                AI12=0D0
                AI21=0D0
                GR12=0D0
                GR21=0D0
                GI12=0D0
                GI21=0D0
                IF (NCHECK.EQ.1.AND.SIG(N1+N2).LT.0D0) GO TO 205

                DO 200 I=1,NGSS
                    D1N1=WS%TMATR_D1(I,N1)
                    D2N1=WS%TMATR_D2(I,N1)
                    D1N2=WS%TMATR_D1(I,N2)
                    D2N2=WS%TMATR_D2(I,N2)
                    A12=D1N1*D2N2
                    A21=D2N1*D1N2
                    A22=D2N1*D2N2
                    AA1=A12+A21

                    QJ1=J(I,N1)
                    QY1=Y(I,N1)
                    QJR2=JR(I,N2)
                    QJI2=JI(I,N2)
                    QDJR2=DJR(I,N2)
                    QDJI2=DJI(I,N2)
                    QDJ1=DJ(I,N1)
                    QDY1=DY(I,N1)

                    C1R=QJR2*QJ1
                    C1I=QJI2*QJ1
                    B1R=C1R-QJI2*QY1
                    B1I=C1I+QJR2*QY1

                    C2R=QJR2*QDJ1
                    C2I=QJI2*QDJ1
                    B2R=C2R-QJI2*QDY1
                    B2I=C2I+QJR2*QDY1

                    DDRI=DDR(I)
                    C3R=DDRI*C1R
                    C3I=DDRI*C1I
                    B3R=DDRI*B1R
                    B3I=DDRI*B1I

                    C4R=QDJR2*QJ1
                    C4I=QDJI2*QJ1
                    B4R=C4R-QDJI2*QY1
                    B4I=C4I+QDJR2*QY1

                    DRRI=DRR(I)
                    DRII=DRI(I)
                    C5R=C1R*DRRI-C1I*DRII
                    C5I=C1I*DRRI+C1R*DRII
                    B5R=B1R*DRRI-B1I*DRII
                    B5I=B1I*DRRI+B1R*DRII

                    URI=DR(I)
                    RRI=RR(I)

                    F1=RRI*A22
                    F2=RRI*URI*AN1*A12
                    AR12=AR12+F1*B2R+F2*B3R
                    AI12=AI12+F1*B2I+F2*B3I
                    GR12=GR12+F1*C2R+F2*C3R
                    GI12=GI12+F1*C2I+F2*C3I

                    F2=RRI*URI*AN2*A21
                    AR21=AR21+F1*B4R+F2*B5R
                    AI21=AI21+F1*B4I+F2*B5I
                    GR21=GR21+F1*C4R+F2*C5R
                    GI21=GI21+F1*C4I+F2*C5I
  200           CONTINUE

  205           AN12=ANN(N1,N2)*FACTOR
                TMATR_WORK(N1,N2,2)=AR12*AN12
                TMATR_WORK(N1,N2,3)=AR21*AN12
                TMATR_WORK(N1,N2,6)=AI12*AN12
                TMATR_WORK(N1,N2,7)=AI21*AN12
                TMATR_WORK(N1,N2,10)=GR12*AN12
                TMATR_WORK(N1,N2,11)=GR21*AN12
                TMATR_WORK(N1,N2,14)=GI12*AN12
                TMATR_WORK(N1,N2,15)=GI21*AN12
  300 CONTINUE

      TPIR=PIR
      TPII=PII
      TPPI=PPI

      NM=NMAX
      DO 310 N1=MM1,NMAX
           K1=N1-MM1+1
           KK1=K1+NM
           DO 310 N2=MM1,NMAX
                K2=N2-MM1+1
                KK2=K2+NM

                TAR12= TMATR_WORK(N1,N2,6)
                TAI12=-TMATR_WORK(N1,N2,2)
                TGR12= TMATR_WORK(N1,N2,14)
                TGI12=-TMATR_WORK(N1,N2,10)

                TAR21=-TMATR_WORK(N1,N2,7)
                TAI21= TMATR_WORK(N1,N2,3)
                TGR21=-TMATR_WORK(N1,N2,15)
                TGI21= TMATR_WORK(N1,N2,11)

                QR(K1,K2)=TPIR*TAR21-TPII*TAI21+TPPI*TAR12
                QI(K1,K2)=TPIR*TAI21+TPII*TAR21+TPPI*TAI12
                RGQR(K1,K2)=TPIR*TGR21-TPII*TGI21+TPPI*TGR12
                RGQI(K1,K2)=TPIR*TGI21+TPII*TGR21+TPPI*TGI12

                QR(K1,KK2)=0D0
                QI(K1,KK2)=0D0
                RGQR(K1,KK2)=0D0
                RGQI(K1,KK2)=0D0

                QR(KK1,K2)=0D0
                QI(KK1,K2)=0D0
                RGQR(KK1,K2)=0D0
                RGQI(KK1,K2)=0D0

                QR(KK1,KK2)=TPIR*TAR12-TPII*TAI12+TPPI*TAR21
                QI(KK1,KK2)=TPIR*TAI12+TPII*TAR12+TPPI*TAI21
                RGQR(KK1,KK2)=TPIR*TGR12-TPII*TGI12+TPPI*TGR21
                RGQI(KK1,KK2)=TPIR*TGI12+TPII*TGR12+TPPI*TGI21
  310 CONTINUE

      CALL TT_FULL_DIRECT(NMAX,NCHECK,IERR,LAPACK_INFO,TR1,TI1,QR,QI,RGQR,RGQI,WS)
      RETURN
      END

      SUBROUTINE TMATR_FULL_DIRECT (M,NGAUSS,X,W,AN,ANN,S,SS,PPI,PIR,PII,R,DR,DDR, &
&                        DRR,DRI,NMAX,NCHECK,IERR,LAPACK_INFO,TR1,TI1,QR,QI,RGQR,RGQI, &
&                        J,Y,JR,JI,DJ,DY,DJR,DJI,TMATR_WORK,WS)
      IMPLICIT real(wp) (A-H,O-Z)
!  D1 and D2 are WS%TMATR_D1 and WS%TMATR_D2, shared with TMATR0.
      TYPE(tmatrix_workspace_t), INTENT(INOUT) :: WS
      real(wp)  X(NPNG2),W(NPNG2),AN(NPN1),S(NPNG2),SS(NPNG2), &
&              R(NPNG2),DR(NPNG2),SIG(NPN2), &
&              J(NPNG2,NPN1),Y(NPNG2,NPN1), &
&              JR(NPNG2,NPN1),JI(NPNG2,NPN1),DJ(NPNG2,NPN1), &
&              DY(NPNG2,NPN1),DJR(NPNG2,NPN1), &
&              DJI(NPNG2,NPN1),DDR(NPNG2),DRR(NPNG2), &
&              DRI(NPNG2),DS(NPNG2),DSS(NPNG2),RR(NPNG2), &
&              DV1(NPN1),DV2(NPN1)

      real(wp)  ANN(NPN1,NPN1),TMATR_WORK(NPN1,NPN1,16)
      real(wp) TR1(NPN2,NPN2),TI1(NPN2,NPN2), &
&             QR(NPN2,NPN2),QI(NPN2,NPN2), &
&             RGQR(NPN2,NPN2),RGQI(NPN2,NPN2)
      MM1=M
      QM=DFLOAT(M)
      QMM=QM*QM
      NG=2*NGAUSS
      NGSS=NG
      FACTOR=1D0
      IF (NCHECK.EQ.1) THEN
            NGSS=NGAUSS
            FACTOR=2D0
         ELSE
            CONTINUE
      ENDIF
      SI=1D0
      NM=NMAX+NMAX
      DO 5 N=1,NM
           SI=-SI
           SIG(N)=SI
    5 CONTINUE
   20 DO 25 I=1,NGAUSS
         I1=NGAUSS+I
         I2=NGAUSS-I+1
         CALL VIG_FULL_DIRECT(X(I1),NMAX,M,DV1,DV2)
         DO 25 N=1,NMAX
            SI=SIG(N)
            DD1=DV1(N)
            DD2=DV2(N)
            WS%TMATR_D1(I1,N)=DD1
            WS%TMATR_D2(I1,N)=DD2
            WS%TMATR_D1(I2,N)=DD1*SI
            WS%TMATR_D2(I2,N)=-DD2*SI
   25 CONTINUE
   30 DO 40 I=1,NGSS
           WR=W(I)*R(I)
           DS(I)=S(I)*QM*WR
           DSS(I)=SS(I)*QMM
           RR(I)=WR
   40 CONTINUE

      DO 300  N1=MM1,NMAX
           AN1=AN(N1)
           DO 300 N2=MM1,NMAX
                AN2=AN(N2)
                AR11=0D0
                AR12=0D0
                AR21=0D0
                AR22=0D0
                AI11=0D0
                AI12=0D0
                AI21=0D0
                AI22=0D0
                GR11=0D0
                GR12=0D0
                GR21=0D0
                GR22=0D0
                GI11=0D0
                GI12=0D0
                GI21=0D0
                GI22=0D0
                SI=SIG(N1+N2)

                DO 200 I=1,NGSS
                    D1N1=WS%TMATR_D1(I,N1)
                    D2N1=WS%TMATR_D2(I,N1)
                    D1N2=WS%TMATR_D1(I,N2)
                    D2N2=WS%TMATR_D2(I,N2)
                    A11=D1N1*D1N2
                    A12=D1N1*D2N2
                    A21=D2N1*D1N2
                    A22=D2N1*D2N2
                    AA1=A12+A21
                    AA2=A11*DSS(I)+A22
                    QJ1=J(I,N1)
                    QY1=Y(I,N1)
                    QJR2=JR(I,N2)
                    QJI2=JI(I,N2)
                    QDJR2=DJR(I,N2)
                    QDJI2=DJI(I,N2)
                    QDJ1=DJ(I,N1)
                    QDY1=DY(I,N1)

                    C1R=QJR2*QJ1
                    C1I=QJI2*QJ1
                    B1R=C1R-QJI2*QY1
                    B1I=C1I+QJR2*QY1

                    C2R=QJR2*QDJ1
                    C2I=QJI2*QDJ1
                    B2R=C2R-QJI2*QDY1
                    B2I=C2I+QJR2*QDY1

                    DDRI=DDR(I)
                    C3R=DDRI*C1R
                    C3I=DDRI*C1I
                    B3R=DDRI*B1R
                    B3I=DDRI*B1I

                    C4R=QDJR2*QJ1
                    C4I=QDJI2*QJ1
                    B4R=C4R-QDJI2*QY1
                    B4I=C4I+QDJR2*QY1

                    DRRI=DRR(I)
                    DRII=DRI(I)
                    C5R=C1R*DRRI-C1I*DRII
                    C5I=C1I*DRRI+C1R*DRII
                    B5R=B1R*DRRI-B1I*DRII
                    B5I=B1I*DRRI+B1R*DRII

                    C6R=QDJR2*QDJ1
                    C6I=QDJI2*QDJ1
                    B6R=C6R-QDJI2*QDY1
                    B6I=C6I+QDJR2*QDY1

                    C7R=C4R*DDRI
                    C7I=C4I*DDRI
                    B7R=B4R*DDRI
                    B7I=B4I*DDRI

                    C8R=C2R*DRRI-C2I*DRII
                    C8I=C2I*DRRI+C2R*DRII
                    B8R=B2R*DRRI-B2I*DRII
                    B8I=B2I*DRRI+B2R*DRII

                    URI=DR(I)
                    DSI=DS(I)
                    DSSI=DSS(I)
                    RRI=RR(I)

                    IF (NCHECK.EQ.1.AND.SI.GT.0D0) GO TO 150

                    E1=DSI*AA1
                    AR11=AR11+E1*B1R
                    AI11=AI11+E1*B1I
                    GR11=GR11+E1*C1R
                    GI11=GI11+E1*C1I
                    IF (NCHECK.EQ.1) GO TO 160

  150               F1=RRI*AA2
                    F2=RRI*URI*AN1*A12
                    AR12=AR12+F1*B2R+F2*B3R
                    AI12=AI12+F1*B2I+F2*B3I
                    GR12=GR12+F1*C2R+F2*C3R
                    GI12=GI12+F1*C2I+F2*C3I

                    F2=RRI*URI*AN2*A21
                    AR21=AR21+F1*B4R+F2*B5R
                    AI21=AI21+F1*B4I+F2*B5I
                    GR21=GR21+F1*C4R+F2*C5R
                    GI21=GI21+F1*C4I+F2*C5I
                    IF (NCHECK.EQ.1) GO TO 200

  160               E2=DSI*URI*A11
                    E3=E2*AN2
                    E2=E2*AN1
                    AR22=AR22+E1*B6R+E2*B7R+E3*B8R
                    AI22=AI22+E1*B6I+E2*B7I+E3*B8I
                    GR22=GR22+E1*C6R+E2*C7R+E3*C8R
                    GI22=GI22+E1*C6I+E2*C7I+E3*C8I
  200           CONTINUE
                AN12=ANN(N1,N2)*FACTOR
                TMATR_WORK(N1,N2,1)=AR11*AN12
                TMATR_WORK(N1,N2,2)=AR12*AN12
                TMATR_WORK(N1,N2,3)=AR21*AN12
                TMATR_WORK(N1,N2,4)=AR22*AN12
                TMATR_WORK(N1,N2,5)=AI11*AN12
                TMATR_WORK(N1,N2,6)=AI12*AN12
                TMATR_WORK(N1,N2,7)=AI21*AN12
                TMATR_WORK(N1,N2,8)=AI22*AN12
                TMATR_WORK(N1,N2,9)=GR11*AN12
                TMATR_WORK(N1,N2,10)=GR12*AN12
                TMATR_WORK(N1,N2,11)=GR21*AN12
                TMATR_WORK(N1,N2,12)=GR22*AN12
                TMATR_WORK(N1,N2,13)=GI11*AN12
                TMATR_WORK(N1,N2,14)=GI12*AN12
                TMATR_WORK(N1,N2,15)=GI21*AN12
                TMATR_WORK(N1,N2,16)=GI22*AN12

  300 CONTINUE
      TPIR=PIR
      TPII=PII
      TPPI=PPI
      NM=NMAX-MM1+1
      DO 310 N1=MM1,NMAX
           K1=N1-MM1+1
           KK1=K1+NM
           DO 310 N2=MM1,NMAX
                K2=N2-MM1+1
                KK2=K2+NM

                TAR11=-TMATR_WORK(N1,N2,1)
                TAI11=-TMATR_WORK(N1,N2,5)
                TGR11=-TMATR_WORK(N1,N2,9)
                TGI11=-TMATR_WORK(N1,N2,13)

                TAR12= TMATR_WORK(N1,N2,6)
                TAI12=-TMATR_WORK(N1,N2,2)
                TGR12= TMATR_WORK(N1,N2,14)
                TGI12=-TMATR_WORK(N1,N2,10)

                TAR21=-TMATR_WORK(N1,N2,7)
                TAI21= TMATR_WORK(N1,N2,3)
                TGR21=-TMATR_WORK(N1,N2,15)
                TGI21= TMATR_WORK(N1,N2,11)

                TAR22=-TMATR_WORK(N1,N2,4)
                TAI22=-TMATR_WORK(N1,N2,8)
                TGR22=-TMATR_WORK(N1,N2,12)
                TGI22=-TMATR_WORK(N1,N2,16)

                QR(K1,K2)=TPIR*TAR21-TPII*TAI21+TPPI*TAR12
                QI(K1,K2)=TPIR*TAI21+TPII*TAR21+TPPI*TAI12
                RGQR(K1,K2)=TPIR*TGR21-TPII*TGI21+TPPI*TGR12
                RGQI(K1,K2)=TPIR*TGI21+TPII*TGR21+TPPI*TGI12

                QR(K1,KK2)=TPIR*TAR11-TPII*TAI11+TPPI*TAR22
                QI(K1,KK2)=TPIR*TAI11+TPII*TAR11+TPPI*TAI22
                RGQR(K1,KK2)=TPIR*TGR11-TPII*TGI11+TPPI*TGR22
                RGQI(K1,KK2)=TPIR*TGI11+TPII*TGR11+TPPI*TGI22

                QR(KK1,K2)=TPIR*TAR22-TPII*TAI22+TPPI*TAR11
                QI(KK1,K2)=TPIR*TAI22+TPII*TAR22+TPPI*TAI11
                RGQR(KK1,K2)=TPIR*TGR22-TPII*TGI22+TPPI*TGR11
                RGQI(KK1,K2)=TPIR*TGI22+TPII*TGR22+TPPI*TGI11

                QR(KK1,KK2)=TPIR*TAR12-TPII*TAI12+TPPI*TAR21
                QI(KK1,KK2)=TPIR*TAI12+TPII*TAR12+TPPI*TAI21
                RGQR(KK1,KK2)=TPIR*TGR12-TPII*TGI12+TPPI*TGR21
                RGQI(KK1,KK2)=TPIR*TGI12+TPII*TGR12+TPPI*TGI21
  310 CONTINUE

      CALL TT_FULL_DIRECT(NM,NCHECK,IERR,LAPACK_INFO,TR1,TI1,QR,QI,RGQR,RGQI,WS)

      RETURN
      END

      SUBROUTINE GSP_FULL_DIRECT(NMAX,CSCA,LAM,ALF1,ALF2,ALF3,ALF4,BET1,BET2,LMAX, &
&                     IERR,B1R,B1I,B2R,B2I,TSTORE,D1,D2,D3,D4,D5R,D5I, &
&                     FAC,SSIGN,WS)
      IMPLICIT real(wp) (A-B,D-H,O-Z),complex(wp) (C)
!  TR1, TR2, TI1, TI2, G1, G2, FR, FI, and FF are WS%GSP_TR1 ... WS%GSP_FF.
      TYPE(tmatrix_workspace_t), INTENT(INOUT) :: WS
      real(wp) LAM,FAC(900),SSIGN(900)
      real(wp)  CSCA,SSI(NPL),SSJ(NPN1), &
&              ALF1(NPL),ALF2(NPL),ALF3(NPL), &
&              ALF4(NPL),BET1(NPL),BET2(NPL), &
&              AR1(NPN4),AR2(NPN4),AI1(NPN4),AI2(NPN4)
      real(real32) B1R(NPL1,NPL1,NPN4),B1I(NPL1,NPL1,NPN4), &
&             B2R(NPL1,NPL1,NPN4),B2I(NPL1,NPL1,NPN4), &
&             D1(NPL1,NPN4,NPN4),D2(NPL1,NPN4,NPN4), &
&             D3(NPL1,NPN4,NPN4),D4(NPL1,NPN4,NPN4), &
&             D5R(NPL1,NPN4,NPN4),D5I(NPL1,NPN4,NPN4), &
&             TSTORE(NPN6,NPN4,NPN4,8)
      complex(wp) CIM(NPN1)

      IERR=0
      CALL FACT_FULL_DIRECT(FAC)
      CALL SIGNUM_FULL_DIRECT(SSIGN)
      LMAX=2*NMAX
      L1MAX=LMAX+1
      CI=(0D0,1D0)
      CIM(1)=CI
      DO 2 I=2,NMAX
         CIM(I)=CIM(I-1)*CI
    2 CONTINUE
      SSI(1)=1D0
      DO 3 I=1,LMAX
         I1=I+1
         SI=DFLOAT(2*I+1)
         SSI(I1)=SI
         IF(I.LE.NMAX) SSJ(I)=DSQRT(SI)
    3 CONTINUE
      CI=-CI
      DO 5 I=1,NMAX
         SI=SSJ(I)
         CCI=CIM(I)
         DO 4 J=1,NMAX
            SJ=1D0/SSJ(J)
            CCJ=CIM(J)*SJ/CCI
            WS%GSP_FR(J,I)=CCJ
            WS%GSP_FI(J,I)=CCJ*CI
            WS%GSP_FF(J,I)=SI*SJ
    4    CONTINUE
    5 CONTINUE
      NMAX1=NMAX+1

! *****  CALCULATION OF THE ARRAYS B1 AND B2  *****

      K1=1
      K2=0
      K3=0
      K4=1
      K5=1
      K6=2

 3300 FORMAT (' B1 AND B2')
      DO 100 N=1,NMAX

! *****  CALCULATION OF THE ARRAYS T1 AND T2  *****


         DO 10 NN=1,NMAX
            M1MAX=MIN0(N,NN)+1
            DO 6 M1=1,M1MAX
               M=M1-1
               L1=NPN6+M
               TT1=TSTORE(M1,N,NN,1)
               TT2=TSTORE(M1,N,NN,2)
               TT3=TSTORE(M1,N,NN,3)
               TT4=TSTORE(M1,N,NN,4)
               TT5=TSTORE(M1,N,NN,5)
               TT6=TSTORE(M1,N,NN,6)
               TT7=TSTORE(M1,N,NN,7)
               TT8=TSTORE(M1,N,NN,8)
               T1=TT1+TT2
               T2=TT3+TT4
               T3=TT5+TT6
               T4=TT7+TT8
               WS%GSP_TR1(L1,NN)=T1+T2
               WS%GSP_TR2(L1,NN)=T1-T2
               WS%GSP_TI1(L1,NN)=T3+T4
               WS%GSP_TI2(L1,NN)=T3-T4
               IF(M.EQ.0) GO TO 6
               L1=NPN6-M
               T1=TT1-TT2
               T2=TT3-TT4
               T3=TT5-TT6
               T4=TT7-TT8
               WS%GSP_TR1(L1,NN)=T1-T2
               WS%GSP_TR2(L1,NN)=T1+T2
               WS%GSP_TI1(L1,NN)=T3-T4
               WS%GSP_TI2(L1,NN)=T3+T4
    6       CONTINUE
   10    CONTINUE

!  *****  END OF THE CALCULATION OF THE ARRAYS T1 AND T2  *****

         NN1MAX=NMAX1+N
         DO 40 NN1=1,NN1MAX
            N1=NN1-1

!  *****  CALCULATION OF THE ARRAYS A1 AND A2  *****

            CALL CCG_FULL_DIRECT(N,N1,NMAX,K1,K2,WS%GSP_G1,IERR,FAC,SSIGN)
            IF (IERR.NE.0) RETURN
            NNMAX=MIN0(NMAX,N1+N)
            NNMIN=MAX0(1,IABS(N-N1))
            KN=N+NN1
            DO 15 NN=NNMIN,NNMAX
               NNN=NN+1
               SIG=SSIGN(KN+NN)
               M1MAX=MIN0(N,NN)+NPN6
               AAR1=0D0
               AAR2=0D0
               AAI1=0D0
               AAI2=0D0
               DO 13 M1=NPN6,M1MAX
                  M=M1-NPN6
                  SSS=WS%GSP_G1(M1,NNN)
                  RR1=WS%GSP_TR1(M1,NN)
                  RI1=WS%GSP_TI1(M1,NN)
                  RR2=WS%GSP_TR2(M1,NN)
                  RI2=WS%GSP_TI2(M1,NN)
                  IF(M.EQ.0) GO TO 12
                  M2=NPN6-M
                  RR1=RR1+WS%GSP_TR1(M2,NN)*SIG
                  RI1=RI1+WS%GSP_TI1(M2,NN)*SIG
                  RR2=RR2+WS%GSP_TR2(M2,NN)*SIG
                  RI2=RI2+WS%GSP_TI2(M2,NN)*SIG
   12             AAR1=AAR1+SSS*RR1
                  AAI1=AAI1+SSS*RI1
                  AAR2=AAR2+SSS*RR2
                  AAI2=AAI2+SSS*RI2
   13          CONTINUE
               XR=WS%GSP_FR(NN,N)
               XI=WS%GSP_FI(NN,N)
               AR1(NN)=AAR1*XR-AAI1*XI
               AI1(NN)=AAR1*XI+AAI1*XR
               AR2(NN)=AAR2*XR-AAI2*XI
               AI2(NN)=AAR2*XI+AAI2*XR
   15       CONTINUE

!  *****  END OF THE CALCULATION OF THE ARRAYS A1 AND A2 ****

            CALL CCG_FULL_DIRECT(N,N1,NMAX,K3,K4,WS%GSP_G2,IERR,FAC,SSIGN)
            IF (IERR.NE.0) RETURN
            M1=MAX0(-N1+1,-N)
            M2=MIN0(N1+1,N)
            M1MAX=M2+NPN6
            M1MIN=M1+NPN6
            DO 30 M1=M1MIN,M1MAX
               BBR1=0D0
               BBI1=0D0
               BBR2=0D0
               BBI2=0D0
               DO 25 NN=NNMIN,NNMAX
                  NNN=NN+1
                  SSS=WS%GSP_G2(M1,NNN)
                  BBR1=BBR1+SSS*AR1(NN)
                  BBI1=BBI1+SSS*AI1(NN)
                  BBR2=BBR2+SSS*AR2(NN)
                  BBI2=BBI2+SSS*AI2(NN)
   25          CONTINUE
               B1R(NN1,M1,N)=BBR1
               B1I(NN1,M1,N)=BBI1
               B2R(NN1,M1,N)=BBR2
               B2I(NN1,M1,N)=BBI2
   30       CONTINUE
   40    CONTINUE
  100 CONTINUE

!  *****  END OF THE CALCULATION OF THE ARRAYS B1 AND B2 ****

!  *****  CALCULATION OF THE ARRAYS D1,D2,D3,D4, AND D5  *****

 3301 FORMAT(' D1, D2, ...')
      DO 200 N=1,NMAX
         DO 190 NN=1,NMAX
            M1=MIN0(N,NN)
            M1MAX=NPN6+M1
            M1MIN=NPN6-M1
            NN1MAX=NMAX1+MIN0(N,NN)
            DO 180 M1=M1MIN,M1MAX
               M=M1-NPN6
               NN1MIN=IABS(M-1)+1
               DD1=0D0
               DD2=0D0
               DO 150 NN1=NN1MIN,NN1MAX
                  XX=SSI(NN1)
                  X1=B1R(NN1,M1,N)
                  X2=B1I(NN1,M1,N)
                  X3=B1R(NN1,M1,NN)
                  X4=B1I(NN1,M1,NN)
                  X5=B2R(NN1,M1,N)
                  X6=B2I(NN1,M1,N)
                  X7=B2R(NN1,M1,NN)
                  X8=B2I(NN1,M1,NN)
                  DD1=DD1+XX*(X1*X3+X2*X4)
                  DD2=DD2+XX*(X5*X7+X6*X8)
  150          CONTINUE
               D1(M1,NN,N)=DD1
               D2(M1,NN,N)=DD2
  180       CONTINUE
            MMAX=MIN0(N,NN+2)
            MMIN=MAX0(-N,-NN+2)
            M1MAX=NPN6+MMAX
            M1MIN=NPN6+MMIN
            DO 186 M1=M1MIN,M1MAX
               M=M1-NPN6
               NN1MIN=IABS(M-1)+1
               DD3=0D0
               DD4=0D0
               DD5R=0D0
               DD5I=0D0
               M2=-M+2+NPN6
               DO 183 NN1=NN1MIN,NN1MAX
                  XX=SSI(NN1)
                  X1=B1R(NN1,M1,N)
                  X2=B1I(NN1,M1,N)
                  X3=B2R(NN1,M1,N)
                  X4=B2I(NN1,M1,N)
                  X5=B1R(NN1,M2,NN)
                  X6=B1I(NN1,M2,NN)
                  X7=B2R(NN1,M2,NN)
                  X8=B2I(NN1,M2,NN)
                  DD3=DD3+XX*(X1*X5+X2*X6)
                  DD4=DD4+XX*(X3*X7+X4*X8)
                  DD5R=DD5R+XX*(X3*X5+X4*X6)
                  DD5I=DD5I+XX*(X4*X5-X3*X6)
  183          CONTINUE
               D3(M1,NN,N)=DD3
               D4(M1,NN,N)=DD4
               D5R(M1,NN,N)=DD5R
               D5I(M1,NN,N)=DD5I
  186       CONTINUE
  190    CONTINUE
  200 CONTINUE

!  *****  END OF THE CALCULATION OF THE D-ARRAYS *****

!  *****  CALCULATION OF THE EXPANSION COEFFICIENTS *****

 3303 FORMAT (' G1, G2, ...')

      DK=LAM*LAM/(4D0*CSCA*DACOS(-1D0))
      DO 300 L1=1,L1MAX
         G1L=0D0
         G2L=0D0
         G3L=0D0
         G4L=0D0
         G5LR=0D0
         G5LI=0D0
         L=L1-1
         SL=SSI(L1)*DK
         DO 290 N=1,NMAX
            NNMIN=MAX0(1,IABS(N-L))
            NNMAX=MIN0(NMAX,N+L)
            IF(NNMAX.LT.NNMIN) GO TO 290
            CALL CCG_FULL_DIRECT(N,L,NMAX,K1,K2,WS%GSP_G1,IERR,FAC,SSIGN)
            IF (IERR.NE.0) RETURN
            IF(L.GE.2) CALL CCG_FULL_DIRECT(N,L,NMAX,K5,K6,WS%GSP_G2,IERR,FAC,SSIGN)
            IF (IERR.NE.0) RETURN
            NL=N+L
            DO 280  NN=NNMIN,NNMAX
               NNN=NN+1
               MMAX=MIN0(N,NN)
               M1MIN=NPN6-MMAX
               M1MAX=NPN6+MMAX
               SI=SSIGN(NL+NNN)
               DM1=0D0
               DM2=0D0
               DO 270 M1=M1MIN,M1MAX
                  M=M1-NPN6
                  IF(M.GE.0) SSS1=WS%GSP_G1(M1,NNN)
                  IF(M.LT.0) SSS1=WS%GSP_G1(NPN6-M,NNN)*SI
                  DM1=DM1+SSS1*D1(M1,NN,N)
                  DM2=DM2+SSS1*D2(M1,NN,N)
  270          CONTINUE
               FFN=WS%GSP_FF(NN,N)
               SSS=WS%GSP_G1(NPN6+1,NNN)*FFN
               G1L=G1L+SSS*DM1
               G2L=G2L+SSS*DM2*SI
               IF(L.LT.2) GO TO 280
               DM3=0D0
               DM4=0D0
               DM5R=0D0
               DM5I=0D0
               MMAX=MIN0(N,NN+2)
               MMIN=MAX0(-N,-NN+2)
               M1MAX=NPN6+MMAX
               M1MIN=NPN6+MMIN
               DO 275 M1=M1MIN,M1MAX
                  M=M1-NPN6
                  SSS1=WS%GSP_G2(NPN6-M,NNN)
                  DM3=DM3+SSS1*D3(M1,NN,N)
                  DM4=DM4+SSS1*D4(M1,NN,N)
                  DM5R=DM5R+SSS1*D5R(M1,NN,N)
                  DM5I=DM5I+SSS1*D5I(M1,NN,N)
  275          CONTINUE
               G5LR=G5LR-SSS*DM5R
               G5LI=G5LI-SSS*DM5I
               SSS=WS%GSP_G2(NPN4,NNN)*FFN
               G3L=G3L+SSS*DM3
               G4L=G4L+SSS*DM4*SI
  280       CONTINUE
  290    CONTINUE
         G1L=G1L*SL
         G2L=G2L*SL
         G3L=G3L*SL
         G4L=G4L*SL
         G5LR=G5LR*SL
         G5LI=G5LI*SL
         ALF1(L1)=G1L+G2L
         ALF2(L1)=G3L+G4L
         ALF3(L1)=G3L-G4L
         ALF4(L1)=G1L-G2L
         BET1(L1)=G5LR*2D0
         BET2(L1)=G5LI*2D0
         LMAX=L
         IF(DABS(G1L).LT.1D-6) GO TO 500
  300 CONTINUE
  500 CONTINUE
      RETURN
      END

!****************************************************************

!   CALCULATION OF THE QUANTITIES F(N+1)=0.5*LN(N!)
!   0.LE.N.LE.899


!*****************************************************************


!**********************************************************************


!**********************************************************************

      SUBROUTINE TT_FULL_DIRECT(NMAX,NCHECK,IERR,LAPACK_INFO,TR1,TI1,QR,QI,RGQR,RGQI,WS)
      IMPLICIT real(wp) (A-H,O-Z)
!  ZQ is WS%TT_ZQ; ZW and IPIV are small and stay on the stack.  The B, WORK,
!  F, A, C, D, E, and IPVT arrays of the original Gaussian-elimination
!  inversion are gone: this routine inverts with ZGETRF/ZGETRI and never
!  referenced them.
      TYPE(tmatrix_workspace_t), INTENT(INOUT) :: WS
      real(wp) TR1(NPN2,NPN2),TI1(NPN2,NPN2), &
&             QR(NPN2,NPN2),QI(NPN2,NPN2), &
&             RGQR(NPN2,NPN2),RGQI(NPN2,NPN2)
      complex(wp) ZW(NPN2)
      INTEGER IPIV(NPN2)
      IERR=0
      LAPACK_INFO=0
      NDIM=NPN2
      NNMAX=2*NMAX

!     Matrix inversion from LAPACK

      DO I=1,NNMAX
         DO J=1,NNMAX
            WS%TT_ZQ(I,J)=DCMPLX(QR(I,J),QI(I,J))
         ENDDO
      ENDDO
      INFO=0
      CALL ZGETRF(NNMAX,NNMAX,WS%TT_ZQ,NPN2,IPIV,INFO)
      IF (INFO.NE.0) THEN
         IERR=8
         LAPACK_INFO=INFO
         RETURN
      ENDIF
      CALL ZGETRI(NNMAX,WS%TT_ZQ,NPN2,IPIV,ZW,NPN2,INFO)
      IF (INFO.NE.0) THEN
         IERR=9
         LAPACK_INFO=INFO
         RETURN
      ENDIF
      DO I=1,NNMAX
         DO J=1,NNMAX
            TR=0D0
            TI=0D0
          DO K=1,NNMAX
                 ARR=RGQR(I,K)
                 ARI=RGQI(I,K)
                 AR=WS%TT_ZQ(K,J)
                 AI=DIMAG(WS%TT_ZQ(K,J))
                 TR=TR-ARR*AR+ARI*AI
                 TI=TI-ARR*AI-ARI*AR
            ENDDO
          TR1(I,J)=TR
            TI1(I,J)=TI
         ENDDO
      ENDDO
      RETURN
      END

!********************************************************************
!                                                                   *
!   CALCULATION OF THE EXPANSION COEFFICIENTS FOR (I,Q,U,V) -       *
!   REPRESENTATION.                                                 *
!                                                                   *
!   INPUT PARAMETERS:                                               *
!                                                                   *
!      LAM - WAVELENGTH OF LIGHT                                    *
!      CSCA - SCATTERING CROSS SECTION                              *
!      TR AND TI - ELEMENTS OF THE T-MATRIX. TRANSFERRED THROUGH    *
!                  explicit caller-owned T-matrix storage            *
!      NMAX - DIMENSION OF T(M)-MATRICES                            *
!                                                                   *
!   OUTPUT INFORTMATION:                                            *
!                                                                   *
!      ALF1,...,ALF4,BET1,BET2 - EXPANSION COEFFICIENTS             *
!      LMAX - NUMBER OF COEFFICIENTS MINUS 1                        *
!                                                                   *
!********************************************************************


!  Inner subroutines below copied verbatim from tmd.lp.f.
!***********************************************************************
      SUBROUTINE FACT_FULL_DIRECT(FAC)
      real(wp) FAC(900)

      FAC(1)=0D0
      FAC(2)=0D0
      DO 2 I=3,900
         I1=I-1
         FAC(I)=FAC(I1)+0.5D0*DLOG(DFLOAT(I1))
    2 CONTINUE
      RETURN
      END

!************************************************************

!   CALCULATION OF THE ARRAY SSIGN(N+1)=SIGN(N)
!   0.LE.N.LE.899

      SUBROUTINE SIGNUM_FULL_DIRECT(SSIGN)
      real(wp) SSIGN(900)

      SSIGN(1)=1D0
      DO 2 N=2,899
         SSIGN(N)=-SSIGN(N-1)
    2 CONTINUE
      RETURN
      END

!******************************************************************
!
!   CALCULATION OF CLEBSCH-GORDAN COEFFICIENTS
!   (N,M:N1,M1/NN,MM)
!   FOR GIVEN N AND N1. M1=MM-M, INDEX MM IS FOUND FROM M AS
!   MM=M*K1+K2
!
!   INPUT PARAMETERS :  N,N1,NMAX,K1,K2
!                               N.LE.NMAX
!                               N.GE.1
!                               N1.GE.0
!                               N1.LE.N+NMAX
!   OUTPUT PARAMETERS : GG(M+NPN6,NN+1) - ARRAY OF THE CORRESPONDING
!                                       COEFFICIENTS
!                               /M/.LE.N
!                               /M1/=/M*(K1-1)+K2/.LE.N1
!                               NN.LE.MIN(N+N1,NMAX)
!                               NN.GE.MAX(/MM/,/N-N1/)
!   IF K1=1 AND K2=0, THEN 0.LE.M.LE.N


      SUBROUTINE CCG_FULL_DIRECT(N,N1,NMAX,K1,K2,GG,IERR,FAC,SSIGN)
      IMPLICIT real(wp) (A-H,O-Z)
      real(wp) GG(NPL1,NPN6),CD(0:NPN5),CU(0:NPN5),FAC(900),SSIGN(900)
      IERR=0
      IF (NMAX.LE.NPN4 .AND. 0.LE.N1 .AND. N1.LE.NMAX+N .AND. N.GE.1 .AND. N.LE.NMAX) GO TO 1
      IERR=7
      RETURN
    1 NNF=MIN0(N+N1,NMAX)
      MIN=NPN6-N
      MF=NPN6+N
      IF(K1.EQ.1.AND.K2.EQ.0) MIN=NPN6
      DO 100 MIND=MIN,MF
         M=MIND-NPN6
         MM=M*K1+K2
         M1=MM-M
         IF(IABS(M1).GT.N1) GO TO 90
         NNL=MAX0(IABS(MM),IABS(N-N1))
         IF(NNL.GT.NNF) GO TO 90
         NNU=N+N1
         NNM=(NNU+NNL)*0.5D0
         IF (NNU.EQ.NNL) NNM=NNL
         CALL CCGIN_FULL_DIRECT(N,N1,M,MM,C,IERR,FAC,SSIGN)
         IF (IERR.NE.0) RETURN
         CU(NNL)=C
         IF (NNL.EQ.NNF) GO TO 50
         C2=0D0
         C1=C
         DO 7 NN=NNL+1,MIN0(NNM,NNF)
            A=DFLOAT((NN+MM)*(NN-MM)*(N1-N+NN))
            A=A*DFLOAT((N-N1+NN)*(N+N1-NN+1)*(N+N1+NN+1))
            A=DFLOAT(4*NN*NN)/A
            A=A*DFLOAT((2*NN+1)*(2*NN-1))
            A=DSQRT(A)
            B=0.5D0*DFLOAT(M-M1)
            D=0D0
            IF(NN.EQ.1) GO TO 5
            B=DFLOAT(2*NN*(NN-1))
            B=DFLOAT((2*M-MM)*NN*(NN-1)-MM*N*(N+1)+ &
&                     MM*N1*(N1+1))/B
            D=DFLOAT(4*(NN-1)*(NN-1))
            D=D*DFLOAT((2*NN-3)*(2*NN-1))
            D=DFLOAT((NN-MM-1)*(NN+MM-1)*(N1-N+NN-1))/D
            D=D*DFLOAT((N-N1+NN-1)*(N+N1-NN+2)*(N+N1+NN))
            D=DSQRT(D)
    5       C=A*(B*C1-D*C2)
            C2=C1
            C1=C
            CU(NN)=C
    7    CONTINUE
         IF (NNF.LE.NNM) GO TO 50
         CALL DIRECT_FULL_DIRECT(N,M,N1,M1,NNU,MM,C,FAC)
         CD(NNU)=C
         IF (NNU.EQ.NNM+1) GO TO 50
         C2=0D0
         C1=C
         DO 12 NN=NNU-1,NNM+1,-1
            A=DFLOAT((NN-MM+1)*(NN+MM+1)*(N1-N+NN+1))
            A=A*DFLOAT((N-N1+NN+1)*(N+N1-NN)*(N+N1+NN+2))
            A=DFLOAT(4*(NN+1)*(NN+1))/A
            A=A*DFLOAT((2*NN+1)*(2*NN+3))
            A=DSQRT(A)
            B=DFLOAT(2*(NN+2)*(NN+1))
            B=DFLOAT((2*M-MM)*(NN+2)*(NN+1)-MM*N*(N+1) &
&                     +MM*N1*(N1+1))/B
            D=DFLOAT(4*(NN+2)*(NN+2))
            D=D*DFLOAT((2*NN+5)*(2*NN+3))
            D=DFLOAT((NN+MM+2)*(NN-MM+2)*(N1-N+NN+2))/D
            D=D*DFLOAT((N-N1+NN+2)*(N+N1-NN-1)*(N+N1+NN+3))
            D=DSQRT(D)
            C=A*(B*C1-D*C2)
            C2=C1
            C1=C
            CD(NN)=C
   12    CONTINUE
   50    DO 9 NN=NNL,NNF
            IF (NN.LE.NNM) GG(MIND,NN+1)=CU(NN)
            IF (NN.GT.NNM) GG(MIND,NN+1)=CD(NN)
    9    CONTINUE
   90    CONTINUE
  100 CONTINUE
      RETURN
      END

!*********************************************************************

      SUBROUTINE DIRECT_FULL_DIRECT(N,M,N1,M1,NN,MM,C,FAC)
      IMPLICIT real(wp) (A-H,O-Z)
      real(wp) FAC(900)

      C=FAC(2*N+1)+FAC(2*N1+1)+FAC(N+N1+M+M1+1)+FAC(N+N1-M-M1+1)
      C=C-FAC(2*(N+N1)+1)-FAC(N+M+1)-FAC(N-M+1)-FAC(N1+M1+1)-FAC(N1-M1+1)
      C=DEXP(C)
      RETURN
      END

!*********************************************************************
!
!   CALCULATION OF THE CLEBCSH-GORDAN COEFFICIENTS
!   G=(N,M:N1,MM-M/NN,MM)
!   FOR GIVEN N,N1,M,MM, WHERE NN=MAX(/MM/,/N-N1/)
!                               /M/.LE.N
!                               /MM-M/.LE.N1
!                               /MM/.LE.N+N1

      SUBROUTINE CCGIN_FULL_DIRECT(N,N1,M,MM,G,IERR,FAC,SSIGN)
      IMPLICIT real(wp) (A-H,O-Z)
      real(wp) FAC(900),SSIGN(900)


      IERR=0
      M1=MM-M
      IF (N.GE.IABS(M) .AND. N1.GE.IABS(M1) .AND. IABS(MM).LE.(N+N1)) GO TO 1
      IERR=7
      RETURN
    1 IF (IABS(MM).GT.IABS(N-N1)) GO TO 100
      L1=N
      L2=N1
      L3=M
      IF(N1.LE.N) GO TO 50
      K=N
      N=N1
      N1=K
      K=M
      M=M1
      M1=K
   50 N2=N*2
      M2=M*2
      N12=N1*2
      M12=M1*2
      G=SSIGN(N1+M1+1) &
&       *DEXP(FAC(N+M+1)+FAC(N-M+1)+FAC(N12+1)+FAC(N2-N12+2)-FAC(N2+2) &
&             -FAC(N1+M1+1)-FAC(N1-M1+1)-FAC(N-N1+MM+1)-FAC(N-N1-MM+1))
      N=L1
      N1=L2
      M=L3
      RETURN
  100 A=1D0
      L1=M
      L2=MM
      IF(MM.GE.0) GO TO 150
      MM=-MM
      M=-M
      M1=-M1
      A=SSIGN(MM+N+N1+1)
  150 G=A*SSIGN(N+M+1) &
&         *DEXP(FAC(2*MM+2)+FAC(N+N1-MM+1)+FAC(N+M+1)+FAC(N1+M1+1) &
&              -FAC(N+N1+MM+2)-FAC(N-N1+MM+1)-FAC(-N+N1+MM+1)-FAC(N-M+1) &
&              -FAC(N1-M1+1))
      M=L1
      MM=L2
      RETURN
      END

!*****************************************************************

      SUBROUTINE CONST_FULL_DIRECT (NGAUSS,NMAX,MMAX,P,X,W,AN,ANN,S,SS,NP,EPS)
      IMPLICIT real(wp) (A-H,O-Z)
      real(wp) X(NPNG2),W(NPNG2),X1(NPNG1),W1(NPNG1), &
&              X2(NPNG1),W2(NPNG1), &
&              S(NPNG2),SS(NPNG2), &
&              AN(NPN1),ANN(NPN1,NPN1),DD(NPN1)

      DO 10 N=1,NMAX
           NN=N*(N+1)
           AN(N)=DFLOAT(NN)
           D=DSQRT(DFLOAT(2*N+1)/DFLOAT(NN))
           DD(N)=D
           DO 10 N1=1,N
                DDD=D*DD(N1)*0.5D0
                ANN(N,N1)=DDD
                ANN(N1,N)=DDD
   10 CONTINUE
      NG=2*NGAUSS
      IF (NP.EQ.-2) GO  TO 11
      CALL GAUSS_FULL_DIRECT(NG,0,0,X,W)
      GO TO 19
   11 NG1=DFLOAT(NGAUSS)/2D0
      NG2=NGAUSS-NG1
      XX=-DCOS(DATAN(EPS))
      CALL GAUSS_FULL_DIRECT(NG1,0,0,X1,W1)
      CALL GAUSS_FULL_DIRECT(NG2,0,0,X2,W2)
      DO 12 I=1,NG1
         W(I)=0.5D0*(XX+1D0)*W1(I)
         X(I)=0.5D0*(XX+1D0)*X1(I)+0.5D0*(XX-1D0)
   12 CONTINUE
      DO 14 I=1,NG2
         W(I+NG1)=-0.5D0*XX*W2(I)
         X(I+NG1)=-0.5D0*XX*X2(I)+0.5D0*XX
   14 CONTINUE
      DO 16 I=1,NGAUSS
         W(NG-I+1)=W(I)
         X(NG-I+1)=-X(I)
   16 CONTINUE
   19 DO 20 I=1,NGAUSS
           Y=X(I)
           Y=1D0/(1D0-Y*Y)
           SS(I)=Y
           SS(NG-I+1)=Y
           Y=DSQRT(Y)
           S(I)=Y
           S(NG-I+1)=Y
   20 CONTINUE
      RETURN
      END

!**********************************************************************
      SUBROUTINE RSP1_FULL_DIRECT (X,NG,NGAUSS,REV,EPS,NP,R,DR)
      IMPLICIT real(wp) (A-H,O-Z)
      real(wp) X(NG),R(NG),DR(NG)
      A=REV*EPS**(1D0/3D0)
      AA=A*A
      EE=EPS*EPS
      EE1=EE-1D0
      DO 50 I=1,NGAUSS
          C=X(I)
          CC=C*C
          SS=1D0-CC
          S=DSQRT(SS)
          RR=1D0/(SS+EE*CC)
          R(I)=AA*RR
          R(NG-I+1)=R(I)
          DR(I)=RR*C*S*EE1
          DR(NG-I+1)=-DR(I)
   50 CONTINUE
      RETURN
      END

!**********************************************************************

      SUBROUTINE RSP2_FULL_DIRECT (X,NG,REV,EPS,N,R,DR)
      IMPLICIT real(wp) (A-H,O-Z)
      real(wp) X(NG),R(NG),DR(NG)
      DNP=DFLOAT(N)
      DN=DNP*DNP
      DN4=DN*4D0
      EP=EPS*EPS
      A=1D0+1.5D0*EP*(DN4-2D0)/(DN4-1D0)
      I=(DNP+0.1D0)*0.5D0
      I=2*I
      IF (I.EQ.N) A=A-3D0*EPS*(1D0+0.25D0*EP)/ &
&                    (DN-1D0)-0.25D0*EP*EPS/(9D0*DN-1D0)
      R0=REV*A**(-1D0/3D0)
      DO 50 I=1,NG
         XI=DACOS(X(I))*DNP
         RI=R0*(1D0+EPS*DCOS(XI))
         R(I)=RI*RI
         DR(I)=-R0*EPS*DNP*DSIN(XI)/RI
   50 CONTINUE
      RETURN
      END

!**********************************************************************

      SUBROUTINE RSP3_FULL_DIRECT (X,NG,NGAUSS,REV,EPS,R,DR)
      IMPLICIT real(wp) (A-H,O-Z)
      real(wp) X(NG),R(NG),DR(NG)
      H=REV*( (2D0/(3D0*EPS*EPS))**(1D0/3D0) )
      A=H*EPS
      DO 50 I=1,NGAUSS
         CO=-X(I)
         SI=DSQRT(1D0-CO*CO)
         IF (SI/CO.GT.A/H) GO TO 20
         RAD=H/CO
         RTHET=H*SI/(CO*CO)
         GO TO 30
   20    RAD=A/SI
         RTHET=-A*CO/(SI*SI)
   30    R(I)=RAD*RAD
         R(NG-I+1)=R(I)
         DR(I)=-RTHET/RAD
         DR(NG-I+1)=-DR(I)
   50 CONTINUE
      RETURN
      END

!************************************************************************
      SUBROUTINE RJB_FULL_DIRECT(X,Y,U,NMAX,NNMAX)
      IMPLICIT real(wp) (A-H,O-Z)
      real(wp) Y(NMAX),U(NMAX),Z(800)
      L=NMAX+NNMAX
      XX=1D0/X
      Z(L)=1D0/(DFLOAT(2*L+1)*XX)
      L1=L-1
      DO 5 I=1,L1
         I1=L-I
         Z(I1)=1D0/(DFLOAT(2*I1+1)*XX-Z(I1+1))
    5 CONTINUE
      Z0=1D0/(XX-Z(1))
      Y0=Z0*DCOS(X)*XX
      Y1=Y0*Z(1)
      U(1)=Y0-Y1*XX
      Y(1)=Y1
      DO 10 I=2,NMAX
         YI1=Y(I-1)
         YI=YI1*Z(I)
         U(I)=YI1-DFLOAT(I)*YI*XX
         Y(I)=YI
   10 CONTINUE
      RETURN
      END

!**********************************************************************

      SUBROUTINE RYB_FULL_DIRECT(X,Y,V,NMAX)
      IMPLICIT real(wp) (A-H,O-Z)
      real(wp) Y(NMAX),V(NMAX)
      C=DCOS(X)
      S=DSIN(X)
      X1=1D0/X
      X2=X1*X1
      X3=X2*X1
      Y1=-C*X2-S*X1
      Y(1)=Y1
      Y(2)=(-3D0*X3+X1)*C-3D0*X2*S
      NMAX1=NMAX-1
      DO 5 I=2,NMAX1
    5     Y(I+1)=DFLOAT(2*I+1)*X1*Y(I)-Y(I-1)
      V(1)=-X1*(C+Y1)
      DO 10 I=2,NMAX
  10       V(I)=Y(I-1)-DFLOAT(I)*X1*Y(I)
      RETURN
      END

!**********************************************************************
!                                                                     *
!   CALCULATION OF SPHERICAL BESSEL FUNCTIONS OF THE FIRST KIND       *
!   J=JR+I*JI OF COMPLEX ARGUMENT X=XR+I*XI OF ORDERS FROM 1 TO NMAX  *
!   BY USING BACKWARD RECURSION. PARAMETR NNMAX DETERMINES NUMERICAL  *
!   ACCURACY. U=UR+I*UI - FUNCTION (1/X)(D/DX)(X*J(X))                *
!                                                                     *
!**********************************************************************

      SUBROUTINE CJB_FULL_DIRECT (XR,XI,YR,YI,UR,UI,NMAX,NNMAX)
      IMPLICIT real(wp) (A-H,O-Z)
      real(wp) YR(NMAX),YI(NMAX),UR(NMAX),UI(NMAX)
      real(wp) CYR(NPN1),CYI(NPN1),CZR(1200),CZI(1200), &
&             CUR(NPN1),CUI(NPN1)
      L=NMAX+NNMAX
      XRXI=1D0/(XR*XR+XI*XI)
      CXXR=XR*XRXI
      CXXI=-XI*XRXI
      QF=1D0/DFLOAT(2*L+1)
      CZR(L)=XR*QF
      CZI(L)=XI*QF
      L1=L-1
      DO I=1,L1
         I1=L-I
         QF=DFLOAT(2*I1+1)
         AR=QF*CXXR-CZR(I1+1)
         AI=QF*CXXI-CZI(I1+1)
         ARI=1D0/(AR*AR+AI*AI)
         CZR(I1)=AR*ARI
         CZI(I1)=-AI*ARI
      ENDDO
      AR=CXXR-CZR(1)
      AI=CXXI-CZI(1)
      ARI=1D0/(AR*AR+AI*AI)
      CZ0R=AR*ARI
      CZ0I=-AI*ARI
      CR=DCOS(XR)*DCOSH(XI)
      CI=-DSIN(XR)*DSINH(XI)
      AR=CZ0R*CR-CZ0I*CI
      AI=CZ0I*CR+CZ0R*CI
      CY0R=AR*CXXR-AI*CXXI
      CY0I=AI*CXXR+AR*CXXI
      CY1R=CY0R*CZR(1)-CY0I*CZI(1)
      CY1I=CY0I*CZR(1)+CY0R*CZI(1)
      AR=CY1R*CXXR-CY1I*CXXI
      AI=CY1I*CXXR+CY1R*CXXI
      CU1R=CY0R-AR
      CU1I=CY0I-AI
      CYR(1)=CY1R
      CYI(1)=CY1I
      CUR(1)=CU1R
      CUI(1)=CU1I
      YR(1)=CY1R
      YI(1)=CY1I
      UR(1)=CU1R
      UI(1)=CU1I
      DO I=2,NMAX
         QI=DFLOAT(I)
         CYI1R=CYR(I-1)
         CYI1I=CYI(I-1)
         CYIR=CYI1R*CZR(I)-CYI1I*CZI(I)
         CYII=CYI1I*CZR(I)+CYI1R*CZI(I)
         AR=CYIR*CXXR-CYII*CXXI
         AI=CYII*CXXR+CYIR*CXXI
         CUIR=CYI1R-QI*AR
         CUII=CYI1I-QI*AI
         CYR(I)=CYIR
         CYI(I)=CYII
         CUR(I)=CUIR
         CUI(I)=CUII
         YR(I)=CYIR
         YI(I)=CYII
         UR(I)=CUIR
         UI(I)=CUII
      ENDDO
      RETURN
      END

!**********************************************************************
      SUBROUTINE VIG_FULL_DIRECT (X, NMAX, M, DV1, DV2)
      IMPLICIT real(wp) (A-H,O-Z)
      real(wp) DV1(NPN1),DV2(NPN1)

      A=1D0
      QS=DSQRT(1D0-X*X)
      QS1=1D0/QS
      DO N=1,NMAX
         DV1(N)=0D0
         DV2(N)=0D0
      ENDDO
      IF (M.NE.0) GO TO 20
      D1=1D0
      D2=X
      DO N=1,NMAX
         QN=DFLOAT(N)
         QN1=DFLOAT(N+1)
         QN2=DFLOAT(2*N+1)
         D3=(QN2*X*D2-QN*D1)/QN1
         DER=QS1*(QN1*QN/QN2)*(-D1+D3)
         DV1(N)=D2
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
         DV1(N)=D2
         DV2(N)=DER
         D1=D2
         D2=D3
      ENDDO
      RETURN
      END

!**********************************************************************
!    CALCULATION OF POINTS AND WEIGHTS OF GAUSSIAN QUADRATURE         *
!    FORMULA. IF IND1 = 0 - ON INTERVAL (-1,1), IF IND1 = 1 - ON      *
!    INTERVAL  (0,1). IF  IND2 = 1 RESULTS ARE PRINTED.               *
!    N - NUMBER OF POINTS                                             *
!    Z - DIVISION POINTS                                              *
!    W - WEIGHTS                                                      *
!**********************************************************************
!    IND2 selected the diagnostic print, which this library does not
!    carry; the argument is kept so the signature stays Mishchenko's.
!**********************************************************************
      SUBROUTINE GAUSS_FULL_DIRECT ( N,IND1,IND2,Z,W )
      IMPLICIT real(wp) (A-H,P-Z)
      real(wp) Z(N),W(N)
      A=1D0
      B=2D0
      C=3D0
      IND=MOD(N,2)
      K=N/2+IND
      F=DFLOAT(N)
      DO 100 I=1,K
          M=N+1-I
          IF(I.EQ.1) X=A-B/((F+A)*F)
          IF(I.EQ.2) X=(Z(N)-A)*4D0+Z(N)
          IF(I.EQ.3) X=(Z(N-1)-Z(N))*1.6D0+Z(N-1)
          IF(I.GT.3) X=(Z(M+1)-Z(M+2))*C+Z(M+3)
          IF(I.EQ.K.AND.IND.EQ.1) X=0D0
          NITER=0
          CHECK=1D-16
   10     PB=1D0
          NITER=NITER+1
          IF (NITER.LE.100) GO TO 15
          CHECK=CHECK*10D0
   15     PC=X
          DJ=A
          DO 20 J=2,N
              DJ=DJ+A
              PA=PB
              PB=PC
   20         PC=X*PB+(X*PB-PA)*(DJ-A)/DJ
          PA=A/((PB-X*PC)*F)
          PB=PA*PC*(A-X*X)
          X=X-PB
          IF(DABS(PB).GT.CHECK*DABS(X)) GO TO 10
          Z(M)=X
          W(M)=PA*PA*(A-X*X)
          IF(IND1.EQ.0) W(M)=B*W(M)
          IF(I.EQ.K.AND.IND.EQ.1) GO TO 100
          Z(I)=-Z(M)
          W(I)=W(M)
  100 CONTINUE
      IF(IND1.EQ.0) GO TO 140
      DO 120 I=1,N
  120     Z(I)=(A+Z(I))/B
  140 CONTINUE
      RETURN
      END

!****************************************************************
end module tmatrix_core
