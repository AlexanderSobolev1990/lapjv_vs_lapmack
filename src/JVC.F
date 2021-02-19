      !#########################################################################
      !>
      !> \brief JONKER VOLGENANT CASTANON METHOD OF SOLVING ASSIGNMENT PROBLEM
      !>
      !> \details THIS SUBROUTINE SOLVES THE FULL DENSITY LINEAR ASSIGNMENT PROBLEM
      !>          ACCORDING TO:
      !>       1) "A Shortest Augmenting Path Algorithm for Dense and Sparse Linear
      !>          Assignment Problems," Computing 38, 325-340, 1987 by
      !>          R. Jonker and A. Volgenant, University of Amsterdam.
      !>
      !>       2) https://github.com/yongyanghz/LAPJV-algorithm-c
      !>
      !>       3) https://www.mathworks.com/matlabcentral/fileexchange/26836-
      !>          lapjv-jonker-volgenant-algorithm-for-linear-assignment-problem-v3-0
      !>
      !> \remark      FIXED TOLERANCE PROBLEMS
      !> \param[in]   IFINDMAX    1-FIND MAX, 0-FIND MIN
      !> \param[in]   N           ORDER OF MATRIX ASSIGNCOST
      !> \param[in]   ASSIGNCOST  COST MATRIX
      !> \param[in]   ABSMAXCOST  ABSOLUTE VALUE OF MAXIMUM COST (INCLUDE EMPTY COSTS)
      !> \param[in]   RESOLUTION  TWO REALS, DIFFER LESS THEN RESOLUTION, SUPPOSED TO BE EQUAL
      !> \param[out]  RES         RESULT OF ASSIGNMENT (INDEX OF COL PER ROW)
      !> \param[out]  Z           SUMM OF ASSIGNED COSTS
      !> \n
      SUBROUTINE JVC(IFINDMAX,N,ASSIGNCOST,ABSMAXCOST,RESOLUTION,RES,Z)
      !-------------------------------------------------------------------------
      ! LIST OF FORMAL PARAMETERS
      INTEGER IFINDMAX        ! 1-FIND MAX, 0-FIND MIN
      INTEGER N               ! ORDER OF MATRIX ASSIGNCOST
      REAL ASSIGNCOST(N,N)    ! COST SQUARE MATRIX
      REAL ABSMAXCOST         ! ABSOLUTE VALUE OF MAXIMUM COST (INCLUDE EMPTY COSTS)
      REAL RESOLUTION         ! TWO REALS, DIFFER LESS THEN RESOLUTION, SUPPOSED TO BE EQUAL
      INTEGER RES(N)          ! RESULT OF ASSIGNMENT (INDEX OF COL PER ROW)
      REAL Z                  ! SUMM OF ASSIGNED COSTS
      !-------------------------------------------------------------------------
      INTEGER ROWSOL(N),COLSOL(N)
      REAL U(N),V(N)
      INTEGER I,IMIN,NUMFREE,PRVNUMFREE,F,I0,K,FREEROW,PRED(N),FREE1(N),
     & J,J1,J2,ENDOFPATH,LAST,LOW,UP,LOOPCNT,COLLIST(N),MATCHES(N)
      REAL DMIN,H,UMIN,USUBMIN,V2,D(N)
      LOGICAL UNASSIGNEDFOUND
      REAL COST(N,N)

      I=0;IMIN=0;NUMFREE=0;PRVNUMFREE=0;F=0;I0=0;K=0;FREEROW=0;PRED=0
      FREE1=0;J=0;J1=0;J2=0;ENDOFPATH=0;LAST=0;LOW=0;UP=0;LOOPCNT=0
      COLLIST=0;MATCHES=0
      DMIN=0.;H=0.;UMIN=0.;USUBMIN=0.;V2=0.;D=0.
      UNASSIGNEDFOUND=.FALSE.
      RES=0
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! MAXIMUM FINDING (*-1) - COMMENT THIS BLOCK IF MINIMUM IS NEEDED
      DO I=1,N
        DO J=1,N
          IF(IFINDMAX.EQ.1)THEN
            COST(J,I)=-ASSIGNCOST(J,I) ! FIND MAX
          ELSEIF(IFINDMAX.EQ.0)THEN
            COST(J,I)=ASSIGNCOST(J,I) ! FIND MIN
          ELSE
            RES=0
            RETURN
          ENDIF
        ENDDO
      ENDDO
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! COLUMN REDUCTION
      DO J=N,1,-1 ! REVERSE ORDER GIVES BETTER RESULTS
        DMIN=COST(J,1)
        IMIN=1
        DO I=2,N
          IF(COST(J,I).LT.DMIN)THEN
            DMIN=COST(J,I)
            IMIN=I
          ENDIF
        ENDDO
        V(J)=DMIN
        MATCHES(IMIN)=MATCHES(IMIN)+1
        IF(MATCHES(IMIN).EQ.1)THEN
          ! INIT ASSIGNMENT IF MINIMUM ROW ASSIGNED FOR FIRST TIME
          ROWSOL(IMIN)=J
          COLSOL(J)=IMIN
        ELSE IF(V(J).LT.V(ROWSOL(IMIN)))THEN
          J1=ROWSOL(IMIN)
          ROWSOL(IMIN)=J
          COLSOL(J)=IMIN
          COLSOL(J1)=-1!0!
        ELSE
          COLSOL(J)=-1!0! ! ROW ALREADY ASSIGNED, COLUMN NOT ASSIGNED
        ENDIF
      ENDDO

      ! REDUCTION TRANSFER
      NUMFREE=0
      DO I=1,N
        IF(MATCHES(I).EQ.0)THEN ! FILL LIST OF UNASSIGNED 'FREE' ROWS
          NUMFREE=NUMFREE+1
          FREE1(NUMFREE)=I
        ELSE IF(MATCHES(I).EQ.1)THEN ! TRANSFER REDUCTION FROM ROWS THAT ARE ASSIGNED ONCE
          J1=ROWSOL(I)
          DMIN=ABSMAXCOST
          DO J=1,N
            IF(J.NE.J1)THEN
              IF((COST(J,I)-V(J)).LT.(DMIN+RESOLUTION))THEN
                DMIN=(COST(J,I)-V(J))
              ENDIF
            ENDIF
          ENDDO
          V(J1)=V(J1)-DMIN
        ENDIF
      ENDDO

      ! AUGMENTING ROW REDUCTION
      DO LOOPCNT=1,2
        ! SCAN ALL FREE ROWS
        ! IN SOME CASES, A FREE ROW MAY BE REPLACED WITH ANOTHER ONE TO BE SCANNED NEXT
        K=1
        PRVNUMFREE=NUMFREE
        NUMFREE=0
        DO WHILE(K.LE.PRVNUMFREE)
          I=FREE1(K)
          K=K+1
          ! FIND MINIMUM AND SECOND MINIMUM REDUCED COST OVER COLUMNS
          UMIN=(COST(1,I)-V(1))
          J1=1
          USUBMIN=ABSMAXCOST!BIG
          DO J=2,N
            H=(COST(J,I)-V(J))
            IF(H.LT.USUBMIN)THEN
              IF((H.GT.UMIN).OR.(ABS(H-UMIN).LT.RESOLUTION))THEN
                USUBMIN=H
                J2=J
              ELSE
                USUBMIN=UMIN
                UMIN=H
                J2=J1
                J1=J
              ENDIF
            ENDIF
          ENDDO
          I0=COLSOL(J1)
          IF((USUBMIN-UMIN).GT.RESOLUTION)THEN
            ! CHANGE THE REDUCTION OF THE MINIMUM COLUMN TO INCREASE THE MINIMUM
            ! REDUCED COST IN THE ROW TO THE SUBMINIMUM
            V(J1)=V(J1)-(USUBMIN-UMIN)
          ELSE ! MINIMUM AND SUBMINIMUM EQUAL
            IF(I0.GT.-1)THEN ! MINIMUM COLUMN J1 IS ASSIGNED
              ! SWAP COLUMNS J1 AND J2, AS J2 MAY BE UNASSIGNED
              J1=J2
              I0=COLSOL(J2)
            ENDIF
          ENDIF
          ! (RE-)ASSIGN I TO J1, POSSIBLY DE-ASSIGNING AN I0
          ROWSOL(I)=J1
          COLSOL(J1)=I
          IF(I0.GT.-1)THEN
            IF((USUBMIN-UMIN).GT.RESOLUTION)THEN
              ! PUT IN CURRENT K, AND GO BACK TO THAT K
              ! CONTINUE AUGMENTING PATH I - J1 WITH I0
              K=K-1
              FREE1(K)=I0
            ELSE
              ! NO FURTHER AUGMENTING REDUCTION POSSIBLE
              ! STORE I0 IN LIST OF FREE ROWS FOR NEXT PHASE
              NUMFREE=NUMFREE+1
              FREE1(NUMFREE)=I0
            ENDIF
          ENDIF
        ENDDO
      ENDDO

      ! AUGMENT SOLUTION FOR EACH FREE ROW
      DO F=1,NUMFREE
        FREEROW=FREE1(F) ! START ROW OF AUGMENTING PATH
        ! DIJKSTRA SHORTEST PATH ALGORITHM
        ! RUNS UNTIL UNASSIGNED COLUMN ADDED TO SHORTEST PATH TREE
        DO J=1,N !DO J=N,1,-1
          D(J)=COST(J,FREEROW)-V(J)
          PRED(J)=FREEROW
          COLLIST(J)=J ! INIT COLUMN LIST
        ENDDO
        ! COLUMNS IN 0..LOW-1 ARE READY, NOW NONE.
        ! COLUMNS IN LOW..UP-1 ARE TO BE SCANNED FOR CURRENT MINIMUM, NOW NONE.
        ! COLUMNS IN UP..DIM-1 ARE TO BE CONSIDERED LATER TO FIND NEW MINIMUM,
        ! AT THIS STAGE THE LIST SIMPLY CONTAINS ALL COLUMNS
        LOW=1
        UP=1
        UNASSIGNEDFOUND=.FALSE.
        DO WHILE(.NOT.UNASSIGNEDFOUND)
          IF(UP.EQ.LOW)THEN ! NO MORE COLUMNS TO BE SCANNED FOR CURRENT MINIMUM
            LAST=LOW-1
            ! SCAN COLUMNS FOR UP..DIM-1 TO FIND ALL INDICES FOR WHICH NEW MINIMUM OCCURS
            ! STORE THESE INDICES BETWEEN LOW..UP-1 (INCREASING UP)
            DMIN=D(COLLIST(UP))
            UP=UP+1
            DO K=UP,N
              J=COLLIST(K)
              H=D(J)
              IF((H.LT.DMIN).OR.(ABS(H-DMIN).LT.RESOLUTION))THEN
                IF(H.LT.DMIN)THEN ! NEW MINIMUM
                  UP=LOW
                  DMIN=H
                ENDIF
                ! NEW INDEX WITH SAME MINIMUM, PUT ON UNDEX UP, AND EXTEND LIST
                COLLIST(K)=COLLIST(UP)
                COLLIST(UP)=J
                UP=UP+1
              ENDIF
            ENDDO
            ! CHECK IF ANY OF THE MINIMUM COLUMNS HAPPENS TO BE UNASSIGNED
            ! IF SO, WE HAVE AN AUGMENTING PATH RIGHT AWAY
            DO K=LOW,(UP-1)
              IF(COLSOL(COLLIST(K)).LT.1)THEN !
                ENDOFPATH=COLLIST(K)
                UNASSIGNEDFOUND=.TRUE.
                EXIT ! BREAK DO LOOP
              ENDIF
            ENDDO
          ENDIF
          IF(.NOT.UNASSIGNEDFOUND)THEN
            ! UPDATE 'DISTANCES' BETWEEN FREEROW AND ALL UNSCANNED COLUMNS, VIA NEXT SCANNED COLUMN
            J1=COLLIST(LOW)
            LOW=LOW+1
            I=COLSOL(J1)
            H=(COST(J1,I)-V(J1)-DMIN)
            DO K=UP,N
              J=COLLIST(K)
              V2=(COST(J,I)-V(J)-H)
              IF(V2.LT.D(J))THEN
                PRED(J)=I
                IF(ABS(V2-DMIN).LT.RESOLUTION)THEN ! NEW COLUMN FOUND AT SAME MINIMUM VALUE!IF(ABS(V2-DMIN).LT.EPS)THEN ! NEW COLUMN FOUND AT SAME MINIMUM VALUE
                  IF(COLSOL(J).LT.0)THEN
                    ! IF UNASSIGNED, SHORTEST AUGMENTING PATH IS COMPLETE
                    ENDOFPATH=J
                    UNASSIGNEDFOUND=.TRUE.
                    EXIT ! BREAK DO LOOP
                  ELSE
                    ! ELSE ADD TO LIST TO BE SCANNED RIGHT AWAY
                    COLLIST(K)=COLLIST(UP)
                    COLLIST(UP)=J
                    UP=UP+1
                  ENDIF
                ENDIF
                D(J)=V2
              ENDIF
            ENDDO
          ENDIF
        ENDDO

        ! UPDATE COLUMN PRICES
        DO K=1,LAST
          J1=COLLIST(K)
          V(J1)=V(J1)+D(J1)-DMIN
        ENDDO

        ! RESET ROW AND COLUMN ASSIGNMENTS ALONG THE ALTERNATING PATH
        DO WHILE(I.NE.FREEROW)
          I=PRED(ENDOFPATH)
          COLSOL(ENDOFPATH)=I
          J1=ENDOFPATH
          ENDOFPATH=ROWSOL(I)
          ROWSOL(I)=J1
        ENDDO
      ENDDO

      ! CALCULATE OPTIMAL COST
      Z=0.0
      DO I=1,N ! DO I=N,1,-1
        J=ROWSOL(I)
        U(I)=(COST(J,I)-V(J))
        Z=Z+ASSIGNCOST(J,I)
        RES(I)=COLSOL(I) ! FORM RESULT
      ENDDO

      ENDSUBROUTINE
