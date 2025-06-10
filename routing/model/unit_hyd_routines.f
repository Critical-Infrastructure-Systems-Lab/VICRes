
c     SUBROUTINES RELATED TO UNIT HYDROGRAPH
C************************************************************************************************************************************************************************************
c     Calculate impulse response function for grid cells
c     using equation (15) from Lohmann, et al. (1996)  Tellus article
c ************************************************************************************************************************************************************************************
      SUBROUTINE MAKE_UHM(UH,VELO,DIFF,XMASK,NCOL,NROW,LE,DT,IROW,ICOL)
      CHARACTER*50 RCSID
      REAL    UH(NCOL,NROW,LE), INTE
      INTEGER I, J, K, LE
      REAL    T, DT, POT, VELO(NCOL,NROW), DIFF(NCOL,NROW)
      REAL    XMASK(NCOL,NROW), H, PI
      PI = 4.0*atan(1.0)
      DO I = 1, ICOL
          DO J =1, IROW
            T = 0.0
            DO K = 1,LE
                T = T + DT
                IF (VELO(I,J) .GT. 0.0) THEN
                  POT = ((VELO(I,J)*T-XMASK(I,J))**2.0)/
     &                  (4.0*DIFF(I,J)*T)
                  
                  IF (POT .GT. 69.0) THEN
                      H = 0.0
                      ELSE
                      H = 1.0/(2.0*SQRT(PI*DIFF(I,J))) *
     &                   XMASK(I,J)/(T**1.5) * EXP(-POT)
                  END IF
                ELSE
                  H = 0.0
                END IF
                UH(I,J,K) = H
            END DO
          END DO
      END DO
      DO I = 1, ICOL
          DO J = 1, IROW
            INTE = 0.0
            DO K = 1,LE
                INTE = INTE + UH(I,J,K)
            END DO
            IF (INTE .GT. 0.0) THEN
                DO K = 1,LE
                  UH(I,J,K) = UH(I,J,K)/INTE
                END DO
            END IF
          END DO
      END DO
      RETURN
      END

C************************************************************************************************************************************************************************************
c Making GRID Unit Hydrograoph
C************************************************************************************************************************************************************************************
      SUBROUTINE MAKE_GRID_UH
     & (DIREC, NOB, UH_DAY, TMAX, PI, PJ, LE, UH_DAILY, KE,
     &  CATCHIJ, UHM, FR, PMAX, NCOL, NROW, UH_BOX,
     &  UH_S, UH_STRING, NAME5,NORESERVOIRS, RESER,
     &	RESERNAME,RES_DIRECT,RESORDER,FINAL,NRESER_MAX,UH_PATH,CLEN,STEPBYSTEP)
      IMPLICIT NONE
      INTEGER UH_DAY, TMAX, PMAX, KE, NRESER_MAX
      INTEGER NOB(NRESER_MAX)
      INTEGER PI, PJ, LE, NCOL, NROW, FINAL,CLEN
      INTEGER DIREC(NCOL,NROW,2)
      REAL    UH_DAILY(PMAX,UH_DAY)
      INTEGER CATCHIJ(PMAX,2,NRESER_MAX)
      REAL    UHM(NCOL,NROW,LE), FR(TMAX,2)
      REAL    UH_BOX(PMAX,KE)
      REAL    UH_S(PMAX,KE+UH_DAY-1,NRESER_MAX)
      INTEGER N, I, J, K, L, T, II, JJ, TT, U
      INTEGER 	RESER(NCOL,NROW)
   	  REAL    SUM
   	  INTEGER stat
   	  INTEGER 	RES_DIRECT(NRESER_MAX,3)
   	  CHARACTER*200 UH_STRING,UH_PATH       !new, AW
   	  CHARACTER(len=*)  NAME5
   	  INTEGER NORESERVOIRS,RESORDER, RESERNAME
      LOGICAL      STEPBYSTEP
   	  RESORDER = NORESERVOIRS
      
   	  DO I = 1, (NORESERVOIRS-1)    ! loop over reservoirs only [NORESERVOIRS is the last element which corresponds to the station under consideration]
        IF (RESERNAME .EQ. RES_DIRECT(I,1)) THEN
			    RESORDER = I
        END IF
   	  END DO
      
      If (FINAL .EQ. 1) THEN    ! this is the case of the station under consideration (not a real reservoir). From SEARCH_WHOLECATCHMENT; NORESERVOIRS=NORESERVOIRS+1 after looping over all reservoirs, i.e., the last reservoir is the station
        RESORDER = NORESERVOIRS
      END IF
      
      IF (UH_STRING(1:4) .ne. 'NONE') THEN       ! read UH_S grid, not make it
        print*, 'reading UH_S grid from file: ',UH_PATH(1:CLEN)//trim(ADJUSTL(NAME5))//'.uh_s'
        !open(98, file=UH_STRING, status='old')
        open(98,file=UH_PATH(1:CLEN)//trim(ADJUSTL(NAME5))//'.uh_s',status='old')
        DO N = 1,NOB(RESORDER)
          READ(98, *) (UH_S(N,K,RESORDER), K = 1,KE+UH_DAY-1)
        END DO
      ELSE		
        print*, 'Making UH_S grid...',' '//trim(NAME5)//'.uh_s'		         
        open(98, iostat=stat, file = UH_PATH(1:CLEN)//trim(ADJUSTL(NAME5))//'.uh_s', status='replace')    ! If the file is successfully opened, IOSTAT will contain the value 0. If the file exists: it is deleted and replaced. If it does not exist: it is created.
        IF (stat .NE. 0) THEN ! If the first open fails (e.g., due to file system permission issues or invalid path), this IF block executes.
          print*, 'Error opening file for UH_S grid. Status:', stat
        ENDIF  

		    DO N = 1, NOB(RESORDER)
          !print*, 'grid cell', N,' out of', NOB(RESORDER)
          DO K = 1,UH_DAY 
            UH_DAILY(N,K) = 0.0
          END DO
          I = CATCHIJ(N,1,RESORDER)
          J = CATCHIJ(N,2,RESORDER)
          DO K = 1,24
            FR(K,1) = 1.0 / 24.0
            FR(K,2) = 0.0
          END DO
          DO K = 25,TMAX
            FR(K,1) = 0.0
            FR(K,2) = 0.0
          END DO
 100      CONTINUE 
         
          IF ((I .NE. PI) .OR. (J .NE. PJ)) THEN
            DO T = 1, TMAX
              DO L = 1, LE
                IF ((T-L) .GT. 0) THEN
                  FR(T,2) = FR(T,2) + FR(T-L,1)*UHM(I,J,L)
                END IF
              END DO
            END DO
            II = DIREC(I,J,1)
            JJ = DIREC(I,J,2)
            I = II
            J = JJ
            DO T = 1, TMAX
              FR(T,1) = FR(T,2)
              FR(T,2) = 0.0
            END DO
          END IF
          
          IF ((I .NE. PI) .OR. (J .NE. PJ)) THEN
            GOTO 100
          END IF
          
          DO T = 1,TMAX
            TT = (T+23)/24
            UH_DAILY(N,TT) = UH_DAILY(N,TT) + FR(T,1)
          END DO
        END DO
        
        DO N = 1,NOB(RESORDER)
          DO K = 1, KE+UH_DAY-1
            UH_S(N,K,RESORDER) = 0.0
          END DO
        END DO
        
        DO N = 1,NOB(RESORDER)
          DO K = 1,KE
            DO U =1,UH_DAY
              UH_S(N,K+U-1,RESORDER) = UH_S(N,K+U-1,RESORDER)+UH_BOX(N,K) * UH_DAILY(N,U)
            END DO
          END DO
          SUM = 0
          DO K = 1,KE+UH_DAY-1
            SUM = SUM + UH_S(N,K,RESORDER)
          END DO
          DO K = 1,KE+UH_DAY-1
            UH_S(N,K,RESORDER) = UH_S(N,K,RESORDER) / SUM
          END DO
        END DO
        
c   write out the grid for future reference...
        DO N = 1,NOB(RESORDER)
          WRITE(98, *) (UH_S(N,K,RESORDER), K = 1,KE+UH_DAY-1)
        END DO
      ENDIF  ! end of if-conditoinal from line 91 (IF (UH_STRING(1:4) .ne. 'NONE') THEN)
      close(98)
      RETURN
      END



C************************************************************************************************************************************************************************************
c Find RES Grid (I,J)
C************************************************************************************************************************************************************************************
       
        SUBROUTINE FIND_RESIJ(RESER, NROW, NCOL, RESID, I_RES,J_RES)
        IMPLICIT NONE
        INTEGER NROW, NCOL
        INTEGER RESER(NCOL,NROW), RESID
        INTEGER I_RES, J_RES
        INTEGER I, J

        IF (RESID==0) THEN
                print*, 'Outlet station or sink cell with no downstream routing..'
        ELSE

                I_RES = -1
                J_RES = -1
        
                DO J = 1, NROW
                        DO I = 1, NCOL
                                IF (RESER(I, J) == RESID) THEN
                                        I_RES = I
                                        J_RES = J
        
                                        RETURN
                                END IF
                        END DO
                END DO
        END IF
        END



C ************************************************************************************************************************************************************************************
c Making Reservoir to Reservoir (or Station) Unit Hydrograoph
C ************************************************************************************************************************************************************************************
        SUBROUTINE MAKE_RESERVOIR_UH
     & (DIREC, GRID_COUNT, UH_DAY, TMAX, PI, PJ, LE, UH_DAILY, KE,
     &  CATCHIJ, UHM, FR, PMAX, NCOL, NROW, UH_BOX,
     &  UH_R, UH_STRING, NAME5,NORESERVOIRS, RESER,
     &	RESID,DS_RESID,RES_DIRECT,RESORDER,FINAL,NRESER_MAX,UH_PATH,CLEN,STEPBYSTEP)
        IMPLICIT NONE
        INTEGER UH_DAY, TMAX, PMAX, KE, NRESER_MAX
        INTEGER PI, PJ, LE, NCOL, NROW, FINAL,CLEN
        INTEGER DIREC(NCOL,NROW,2)
        REAL    UH_DAILY(1,UH_DAY)
        INTEGER CATCHIJ(PMAX,2,NRESER_MAX)
        REAL    UHM(NCOL,NROW,LE), FR(TMAX,2)
        REAL    UH_BOX(PMAX,KE)
        REAL    UH_R(1,KE+UH_DAY-1)  ! UH_R: Unit Hydrograph to Reservoir
        INTEGER I, J, K, L, T, II, JJ, TT, U,I_DS,J_DS,GRID_COUNT
        INTEGER 	RESER(NCOL,NROW)
        REAL    SUM
        INTEGER stat,M
        INTEGER RES_DIRECT(NRESER_MAX,3)
        CHARACTER*100 UH_STRING,UH_PATH,UH_NAME       !new, AW
        CHARACTER*10  NAME5,RESID_STR,DS_RESID_STR,US_RESID_STR
        INTEGER NORESERVOIRS,RESORDER,RESID,DS_RESID
        LOGICAL STEPBYSTEP,RES_IN_BASIN

                   
10      FORMAT(I10) 
        WRITE(RESID_STR,10) RESID
        DS_RESID=RES_DIRECT(RESORDER,2)
        WRITE(DS_RESID_STR,10) DS_RESID  
        ! check if the downstream reservoir is part of the current sub basin (not downstream of the outlet station)
        RES_IN_BASIN = .FALSE.
        DO M = 1, NORESERVOIRS
            IF (RES_DIRECT(M,1) == DS_RESID) THEN
                RES_IN_BASIN = .TRUE.
                EXIT
            END IF
        END DO

        IF (.NOT. RES_IN_BASIN) THEN
            PRINT *, 'Downstream Reservoir: ',trim(adjustl(DS_RESID_STR)), ' is outside the sub-basin (downstream of current station)'
        END IF


        IF (DS_RESID>0 .AND. RES_IN_BASIN) THEN
          CALL FIND_RESIJ(RESER,NROW,NCOL,DS_RESID,I_DS,J_DS)
          UH_NAME = 'RES'//trim(adjustl(RESID_STR))//"_"//'RES'//trim(adjustl(DS_RESID_STR))
        ELSE
          DS_RESID=0
          I_DS=PI
          J_DS=PJ       
          UH_NAME ='RES'//trim(adjustl(RESID_STR))//"_"//trim(adjustl(NAME5))
        END IF
        print*, 'Making UH_R grid...',' '//trim(UH_NAME)//'.uh_r' 
        open(88, iostat=stat, file = UH_PATH(1:CLEN)//trim(ADJUSTL(UH_NAME))//'.uh_r', status='replace')    ! If the file is successfully opened, IOSTAT will contain the value 0. If the file exists: it is deleted and replaced. If it does not exist: it is created.
        IF (stat .NE. 0) THEN ! If the first open fails (e.g., due to file system permission issues or invalid path), this block executes.
          print*, 'Error opening file for UH_R grid. Status:', stat
        END IF  
        
        DO K = 1,UH_DAY 
            UH_DAILY(1,K) = 0.0
        END DO
        I = CATCHIJ(GRID_COUNT,1,NORESERVOIRS)   ! use NORESERVOIRS since the GRID_COUNT is the grid number (N) in the catchment of the station NOB(I) (and not the grid number (N) in the catchment of the reservoir NOB(RESORDER)
        J = CATCHIJ(GRID_COUNT,2,NORESERVOIRS)
        DO K = 1,24
          FR(K,1) = 1.0 / 24.0
          FR(K,2) = 0.0
        END DO
        DO K = 25,TMAX
          FR(K,1) = 0.0
          FR(K,2) = 0.0
        END DO
        
  100      CONTINUE 
        
        IF ((I .NE. I_DS) .OR. (J .NE. J_DS)) THEN
          DO T = 1, TMAX
            DO L = 1, LE
              IF ((T-L) .GT. 0) THEN
                FR(T,2) = FR(T,2) + FR(T-L,1)*UHM(I,J,L)
                
              END IF
            END DO
          END DO
          II = DIREC(I,J,1)
          JJ = DIREC(I,J,2)
          I = II
          J = JJ
          DO T = 1, TMAX
            FR(T,1) = FR(T,2)
            FR(T,2) = 0.0
          END DO
        END IF
        
        IF ((I .NE. I_DS) .OR. (J .NE. J_DS)) THEN
          GOTO 100
        END IF
        
        DO T = 1,TMAX
          TT = (T+23)/24
          UH_DAILY(1,TT) = UH_DAILY(1,TT) + FR(T,1)
          
        END DO

        DO K = 1, KE+UH_DAY-1
          UH_R(1,K) = 0.0
        END DO

        DO K = 1,KE
          DO U =1,UH_DAY
            UH_R(1,K+U-1) = UH_R(1,K+U-1)+UH_BOX(GRID_COUNT,K) * UH_DAILY(1,U)
            
          END DO
        END DO
        
        SUM = 0
        DO K = 1,KE+UH_DAY-1
          SUM = SUM + UH_R(1,K)
        END DO
        
        DO K = 1,KE+UH_DAY-1
          UH_R(1,K) = UH_R(1,K) / SUM
        END DO
        
c   write out the grid for future reference...
            WRITE(88, *) (UH_R(1,K), K = 1,KE+UH_DAY-1)
            close(88)
        
        RETURN
        END
          
