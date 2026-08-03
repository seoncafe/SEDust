module scattering_matrix_expansion
   !! Angular scattering matrix reconstructed from its generalized-spherical-
   !! function expansion, and the van der Mee & Hovenier necessary conditions
   !! on that expansion.
   !!
   !! Both routines are pure array arithmetic on a set of expansion
   !! coefficients: they hold no state, need no workspace, and may be called
   !! concurrently.  They are silent and return status instead of printing,
   !! which is what separates them from Mishchenko's MATR and HOVENR of
   !! ampld.lp.f / tmd.lp.f:
   !!
   !!   - SCATMAT_FROM_MOMENTS uses the angular expansion of MATR, with
   !!     identical recurrences and identical normalization, but writes the
   !!     six scattering-matrix elements into caller-supplied arrays instead
   !!     of PRINTing them.
   !!   - VDM_HOVENIER_TEST applies the same necessary conditions as HOVENR,
   !!     but returns the verdict in KONTR instead of PRINTing it, and
   !!     accumulates it over all L.  HOVENR resets its verdict at the top of
   !!     the L loop, so its message reflects only the last L.
   !!
   !! Shared DO termination labels are written as nested DO/END DO.  No
   !! expression, operation order, array rank or extent, or working
   !! precision differs from the original.
   use tmatrix_kinds, only: wp
   use tmatrix_types, only: tmatrix_expansion_size
   implicit real(wp) (a-h, o-z)
   private
   public :: scatmat_from_moments, vdm_hovenier_test

   integer, parameter :: npl = tmatrix_expansion_size

contains

!****************************************************************
!
!    A1,...,B2 - expansion coefficients, dimensioned NPL
!    LMAX      - number of coefficients minus 1
!    NPNA      - number of scattering angles; angle i (i = 1..NPNA) is
!                180*(I-1)/(NPNA-1) degrees, returned in THETA
!    THETA     - scattering angle [degrees], dimension NPNA
!    F11,...,F34 - scattering-matrix elements, dimension NPNA.
!                Unnormalized in the same sense as MATR: F11 carries the
!                normalization of A1 (F11 = 1 isotropic when A1(1) = 1),
!                and F22, F33, F44, F12, F34 are NOT divided by F11.
!                The degree of linear polarization for unpolarized
!                incident light is -F12/F11.

      SUBROUTINE SCATMAT_FROM_MOMENTS(A1,A2,A3,A4,B1,B2,LMAX,NPNA, &
&                         THETA,F11,F22,F33,F44,F12,F34)
      IMPLICIT real(wp) (A-H,O-Z)
      INTEGER, INTENT(IN) :: LMAX,NPNA
      real(wp), INTENT(IN) :: A1(NPL),A2(NPL),A3(NPL),A4(NPL),B1(NPL),B2(NPL)
      real(wp), INTENT(OUT) :: THETA(NPNA),F11(NPNA),F22(NPNA),F33(NPNA), &
&                              F44(NPNA),F12(NPNA),F34(NPNA)
      N=NPNA
      DN=1D0/DFLOAT(N-1)
      DA=DACOS(-1D0)*DN
      DB=180D0*DN
      L1MAX=LMAX+1
      TB=-DB
      TAA=-DA
      D6=DSQRT(6D0)*0.25D0
      DO I1=1,N
         TAA=TAA+DA
         TB=TB+DB
         U=DCOS(TAA)
         FF11=0D0
         F2=0D0
         F3=0D0
         FF44=0D0
         FF12=0D0
         FF34=0D0
         P1=0D0
         P2=0D0
         P3=0D0
         P4=0D0
         PP1=1D0
         PP2=0.25D0*(1D0+U)*(1D0+U)
         PP3=0.25D0*(1D0-U)*(1D0-U)
         PP4=D6*(U*U-1D0)
         DO L1=1,L1MAX
            L=L1-1
            DL=DFLOAT(L)
            DL1=DFLOAT(L1)
            FF11=FF11+A1(L1)*PP1
            FF44=FF44+A4(L1)*PP1
            IF(L.EQ.LMAX) GO TO 350
            PL1=DFLOAT(2*L+1)
            P=(PL1*U*PP1-DL*P1)/DL1
            P1=PP1
            PP1=P
  350       IF(L.LT.2) GO TO 400
            F2=F2+(A2(L1)+A3(L1))*PP2
            F3=F3+(A2(L1)-A3(L1))*PP3
            FF12=FF12+B1(L1)*PP4
            FF34=FF34+B2(L1)*PP4
            IF(L.EQ.LMAX) GO TO 400
            PL2=DFLOAT(L*L1)*U
            PL3=DFLOAT(L1*(L*L-4))
            PL4=1D0/DFLOAT(L*(L1*L1-4))
            P=(PL1*(PL2-4D0)*PP2-PL3*P2)*PL4
            P2=PP2
            PP2=P
            P=(PL1*(PL2+4D0)*PP3-PL3*P3)*PL4
            P3=PP3
            PP3=P
            P=(PL1*U*PP4-DSQRT(DFLOAT(L*L-4))*P4)/DSQRT(DFLOAT(L1*L1-4))
            P4=PP4
            PP4=P
  400       CONTINUE
         ENDDO
         THETA(I1)=TB
         F11(I1)=FF11
         F22(I1)=(F2+F3)*0.5D0
         F33(I1)=(F2-F3)*0.5D0
         F44(I1)=FF44
         F12(I1)=FF12
         F34(I1)=FF34
      ENDDO
      RETURN
      END

!****************************************************************
!
!    L1        - number of coefficients (LMAX+1)
!    A1,...,B2 - expansion coefficients, dimension at least L1
!    KONTR     - 1 test satisfied, 2 test violated
!    LVIOL     - lowest L at which the test fails, -1 if none

      SUBROUTINE VDM_HOVENIER_TEST(L1,A1,A2,A3,A4,B1,B2,KONTR,LVIOL)
      IMPLICIT real(wp) (A-H,O-Z)
      INTEGER, INTENT(IN) :: L1
      real(wp), INTENT(IN) :: A1(L1),A2(L1),A3(L1),A4(L1),B1(L1),B2(L1)
      INTEGER, INTENT(OUT) :: KONTR,LVIOL
      INTEGER KL
      KONTR=1
      LVIOL=-1
      DO L=1,L1
         KL=1
         LL=L-1
         DL=DFLOAT(LL)*2D0+1D0
         DDL=DL*0.48D0
         AA1=A1(L)
         AA2=A2(L)
         AA3=A3(L)
         AA4=A4(L)
         BB1=B1(L)
         BB2=B2(L)
         IF(LL.GE.1.AND.DABS(AA1).GE.DL) KL=2
         IF(DABS(AA2).GE.DL) KL=2
         IF(DABS(AA3).GE.DL) KL=2
         IF(DABS(AA4).GE.DL) KL=2
         IF(DABS(BB1).GE.DDL) KL=2
         IF(DABS(BB2).GE.DDL) KL=2
         C=-0.1D0
         DO I=1,11
            C=C+0.1D0
            CC=C*C
            C1=CC*BB2*BB2
            C2=C*AA4
            C3=C*AA3
            IF((DL-C*AA1)*(DL-C*AA2)-CC*BB1*BB1.LE.-1D-4) KL=2
            IF((DL-C2)*(DL-C3)+C1.LE.-1D-4) KL=2
            IF((DL+C2)*(DL-C3)-C1.LE.-1D-4) KL=2
            IF((DL-C2)*(DL+C3)-C1.LE.-1D-4) KL=2
         ENDDO
         IF(KL.EQ.2) THEN
            KONTR=2
            IF(LVIOL.LT.0) LVIOL=LL
         ENDIF
      ENDDO
      RETURN
      END

end module scattering_matrix_expansion
