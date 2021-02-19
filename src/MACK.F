      !#########################################################################
      !>
      !> \brief MACK METHOD OF SOLVING ASSIGNMENT PROBLEM
      !>
      !> \details THIS SUBROUTINE SOLVES THE FULL DENSITY LINEAR ASSIGNMENT PROBLEM
      !>          ACCORDING TO Банди Б. Основы линейного программирования:
      !>          Пер. с англ. - М.:Радио м связь, 1989, стр 113-123
      !> \param[in]   IFINDMAX    1-FIND MAX, 0-FIND MIN
      !> \param[in]   N           ORDER OF MATRIX ASSIGNCOST
      !> \param[in]   IASSIGNCOST COST MATRIX
      !> \param[out]  RES         RESULT OF ASSIGNMENT (INDEX OF COL PER ROW)
      !> \param[out]  Z           SUMM OF ASSIGNED COSTS
      !> \n
      SUBROUTINE MACK(IFINDMAX,N,ASSIGNCOST,RES,Z)
      !INTEGER, PARAMETER :: NSEL_LIMIT = 2048  ! MAX NUMBER OF ASSIGNMENT PROBLEM RESULT ARRAY
      !INTEGER, PARAMETER :: NCLU_LIMIT = 512  ! MAX NUMBER OF ASSIGNMENT PROBLEM RESULT ARRAY
      !-------------------------------------------------------------------------
      ! LIST OF FORMAL PARAMETERS
      INTEGER IFINDMAX        ! 1-FIND MAX, 0-FIND MIN
      INTEGER N               ! ORDER OF MATRIX ASSIGNCOST
      REAL ASSIGNCOST(N,N)    ! COST SQUARE MATRIX
      INTEGER RES(N)          ! RESULT OF ASSIGNMENT (INDEX OF COL PER ROW)
      REAL Z                  ! SUMM OF ASSIGNED COSTS
      !-------------------------------------------------------------------------

C      REAL P(512,512)     ! CLUSTER (VALUES NOT SAVED)
C      INTEGER(4) JV(2048) ! OUTPUTTED ARRAY OF OPTIMAL APPOINTMENTS
C      !-------------------------------------------------------------------------
C      REAL MA(512),MB(512)
C      INTEGER(4) IP(512),LR(512),IMM(512)
C      INTEGER(4) NM(512),IC(512,512)
C      INTEGER(4) JR(512),JM(512),JK(512)
      !INTEGER,PARAMETER::NORD=1000
      REAL P(N,N)     ! CLUSTER (VALUES NOT SAVED)
      INTEGER(4) JV(2048) ! OUTPUTTED ARRAY OF OPTIMAL APPOINTMENTS
      !-------------------------------------------------------------------------
      REAL MA(2000),MB(2000)
      INTEGER(4) IP(2000),LR(2000),IMM(2000)
      INTEGER(4) NM(2000),IC(2000,2000)
      INTEGER(4) JR(2000),JM(2000),JK(2000)
      !INTEGER IFINDMAX

C      REAL P(64,64)     ! CLUSTER (VALUES NOT SAVED)
C      INTEGER(4) JV(2048) ! OUTPUTTED ARRAY OF OPTIMAL APPOINTMENTS
C      !-------------------------------------------------------------------------
C      REAL MA(64),MB(64)
C      INTEGER(4) IP(64),LR(64),IMM(64)
C      INTEGER(4) NM(64),IC(64,64)
C      INTEGER(4) JR(64),JM(64),JK(64)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! SUBROITINE MAY RANDOMLY FAIL IF NOT DO THIS ?!
      RES=0
      IP=0
      LR=0
      IMM=0
      NM=0
      IC=0
      NCLU=0
      JR=0
      JM=0
      JK=0
      JV=0
      !
      JC=0
      IR=0
      L=0

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !PRINT *,"MACK N=",N
      DO 22 I=1,N
      NM(I)=0
      DO 22 J=1,N
C     *******  MAXIMUM FINDING (*-1) *******
      !P(J,I)=-P(J,I)
      IF(IFINDMAX.EQ.1)THEN
        P(J,I)=-ASSIGNCOST(J,I) ! FIND MAX
      ELSEIF(IFINDMAX.EQ.0)THEN
        P(J,I)=ASSIGNCOST(J,I) ! FIND MIN
      ELSE
        RES=0
        RETURN
      ENDIF
   22 CONTINUE
      RMA=1.E10
      DO 36 I=1,N
      RIM=RMA
      DO 31 J=1,N
      RIZ=P(I,J)
      IF(RIZ.GT.RIM) GOTO 31
      RIM=RIZ
      L=J
   31 CONTINUE
      NM(L)=NM(L)+1
      K=NM(L)
      IC(L,K)=I
      MA(I)=RIM
      JR(I)=L
   36 CONTINUE
   40 J=0
   41 J=J+1
      IF(J.GT.N) GOTO 150
      IF(NM(J).LT.2) GOTO 41
      JU1=NM(J)
      DO 47 I=1,N
      IP(I)=IC(J,I)
   47 CONTINUE
      NC=1
      LR(1)=J
      JK(J)=1
      MB(J)=0
   52 RIV=RMA
      DO 64 K=1,JU1
      I=IP(K)
      DO 63 JD=1,N
      IF(JK(JD).EQ.1) GOTO 63
      RIZ=P(I,JD)-MA(I)
      IF(RIZ.GT.RIV) GOTO 63
      RIV=RIZ
      JC=JD
      IR=I
   63 CONTINUE
   64 CONTINUE
      DO 67 JX=1,NC
      ILR=LR(JX)
      MB(ILR)=MB(ILR)+RIV
   67 CONTINUE
      DO 70 K=1,JU1
      IIP=IP(K)
      MA(IIP)=MA(IIP)+RIV
   70 CONTINUE
      MB(JC)=0
      JK(JC)=1
      NC=NC+1
      LR(NC)=JC
      IMM(JC)=IR
      JM(IR)=JC
      JY=NM(JC)
      IF(JY.EQ.0) GOTO 84
      DO 82 JX=1,JY
      JU1=JU1+1
      IP(JU1)=IC(JC,JX)
   82 CONTINUE
      GOTO 52
   84 CONTINUE
      DO 89 JX=1,NC
      LS=LR(JX)
      JK(LS)=0
      DO 90 I=1,N
      P(I,LS)=P(I,LS)+MB(LS)
   90 CONTINUE
   89 CONTINUE
      NM(JC)=1
      IC(JC,1)=IR
   93 JP=JR(IR)
      JR(IR)=JC
      IW=0
      JQ=NM(JP)
      DO 102 IL=1,JQ
      IZ=IC(JP,IL)
      IF(IZ.EQ.IR) GOTO 102
      IW=IW+1
      IC(JP,IW)=IC(JP,IL)
  102 CONTINUE
      IF(JQ.GT.1) GOTO 108
      IR=IMM(JP)
      JC=JP
      IC(JP,JQ)=IR
      GOTO 93
  108 NM(JP)=JQ-1
      GOTO 40
  150 CONTINUE
      DO 153 J=1,N
      ICJ=IC(J,1)
      JV(ICJ)=J
  153 CONTINUE
      ! FORM RESULT:
      Z=0.0
      DO I=1,N
        RES(I)=JV(I)
        Z=Z+ASSIGNCOST(I,JV(I))
      ENDDO
      RETURN
      ENDSUBROUTINE
