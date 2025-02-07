c     SUBROUTINES RELATED TO WRITING RESULTS TO FILES

C************************************************************************************************************************************************************************************
C     WRITE DAILY DATA AT THE BASIN OUTLET
C************************************************************************************************************************************************************************************   
      SUBROUTINE WRITE_DATA
     & (FLOW, DAYS, NAME5, FACTOR_SUM, OUTPATH,IDAY,IMONTH,IYEAR,
     &  VOL, FLOWIN, FLOWOUT, HHO, HYDROPOWER,HTK, RESER, NCOL,NROW,
     &  ICOL, IROW, RPATH, NORESERVOIRS,RES_DIRECT,RES_EVAPORATION,NO_OF_BOX, NRESER_MAX,STEPBYSTEP,NDAY_SIM,SIM_YEAR, SIM_MON, SIM_DAY)
      IMPLICIT NONE
c     Declare variables
c     RCS ID STRING
      INTEGER       NRESER_MAX,NDAY_SIM,SIM_YEAR, SIM_MON, SIM_DAY
      INTEGER       IROW, ICOL
      INTEGER       NO_OF_BOX(NRESER_MAX)
      INTEGER       DAYS
      INTEGER       NROW, NCOL
      INTEGER       IDAY(DAYS), IMONTH(DAYS), IYEAR(DAYS)
      INTEGER       I, J, K, CLEN
      INTEGER       RES_DIRECT(NRESER_MAX,3)
      INTEGER       RESER(NCOL,NROW)
      INTEGER       H(NCOL,NROW)
      INTEGER       RULE, IRRIGATION
      INTEGER       NORESERVOIRS, NAM
      REAL          VOL(NRESER_MAX,DAYS),FLOWIN(NRESER_MAX,DAYS),FLOWOUT(NRESER_MAX,DAYS)
      REAL          HHO(NRESER_MAX,DAYS),HYDROPOWER(NRESER_MAX,DAYS),HTK(NRESER_MAX,DAYS)
      REAL          SEEPAGE, INFILTRATION
      REAL          RES_EVAPORATION(NRESER_MAX,DAYS)
      REAL          HRESERMAX, HRESERMIN
      REAL          QRESER,VRESER,HPHATDIEN,QSPILL,QTUR,VDEAD,VINI
      REAL          X1,X2,X3
      REAL          FLOW(DAYS)
      REAL          FACTOR_SUM
      CHARACTER*201 TEMPRPATH, ABC
      CHARACTER*72  RPATH, IPATH
      CHARACTER*20  CHUOI
      CHARACTER*10   NAME5
      CHARACTER*72  OUTPATH, RESERNAME
      LOGICAL   STEPBYSTEP

c     Subroutine main body   
      IF (STEPBYSTEP) THEN
            print*, 'writing data for step-by-step version...'
            CLEN = INDEX(OUTPATH,' ')-1
            IF (NDAY_SIM.EQ.1) THEN
                  OPEN(30, FILE = OUTPATH(1:CLEN)//'/STEPBYSTEP/'//trim(ADJUSTL(NAME5))//'.day')
                  OPEN(31, FILE = OUTPATH(1:CLEN)//'/STEPBYSTEP/'//trim(ADJUSTL(NAME5))//'.day_mm')
                  WRITE(30,*) SIM_YEAR, SIM_MON, SIM_DAY,FLOW(NDAY_SIM)    
                  WRITE(31,*) SIM_YEAR, SIM_MON, SIM_DAY,FLOW(NDAY_SIM) / FACTOR_SUM
            ELSE
                  OPEN(30, FILE = OUTPATH(1:CLEN)//'/STEPBYSTEP/'//trim(ADJUSTL(NAME5))//'.day',status='old', position='append',ERR=9004)
                  OPEN(31, FILE = OUTPATH(1:CLEN)//'/STEPBYSTEP/'//trim(ADJUSTL(NAME5))//'.day_mm',status='old', position='append',ERR=9004)
                  WRITE(30,*) SIM_YEAR, SIM_MON, SIM_DAY,FLOW(NDAY_SIM)    
                  WRITE(31,*) SIM_YEAR, SIM_MON, SIM_DAY,FLOW(NDAY_SIM) / FACTOR_SUM 
            END IF 
            CLOSE(30)
            CLOSE(31)
      ELSE      
            print*, 'writing data for all steps...',NAME5                
            CLEN = INDEX(OUTPATH,' ')-1
            OPEN(30, FILE = OUTPATH(1:CLEN)//trim(ADJUSTL(NAME5))//'.day')
            OPEN(31, FILE = OUTPATH(1:CLEN)//trim(ADJUSTL(NAME5))//'.day_mm')
            DO I = 1,DAYS
                  WRITE(30,*) IYEAR(I),IMONTH(I),IDAY(I),FLOW(I)    
                  WRITE(31,*) IYEAR(I),IMONTH(I),IDAY(I),FLOW(I) / FACTOR_SUM
            END DO
            CLOSE(30)
            CLOSE(31)
      END IF
      RETURN
      STOP     
9004  WRITE(0,*) 'ERROR in Opening Output File; Please Check if Previous Time Step is Simulated'        
      END



C************************************************************************************************************************************************************************************
C     WRITE RESERVOIRS DATA AT THE BASIN OUTLET
C************************************************************************************************************************************************************************************
c     Write reservoirs output of the last station
      SUBROUTINE WRITE_RES
     & (FLOW, DAYS, NAME5, FACTOR_SUM, OUTPATH,IDAY,IMONTH,IYEAR,
     &  VOL, FLOWIN, FLOWOUT, HHO, HYDROPOWER,HTK, RESER, NCOL,NROW,
     &  ICOL, IROW, RPATH, NORESERVOIRS,RES_DIRECT,RES_EVAPORATION,NO_OF_BOX,VRESER, NRESER_MAX,STEPBYSTEP,NDAY_SIM,SIM_YEAR, SIM_MON, SIM_DAY)
      IMPLICIT NONE
c     Declare variables
c     RCS ID STRING
      INTEGER       NSTATIONS_MAX, NRESER_MAX,NDAY_SIM,STEPS,SIM_YEAR, SIM_MON, SIM_DAY
      INTEGER       IROW, ICOL
      INTEGER       NO_OF_BOX(NRESER_MAX)
      INTEGER       DAYS
      INTEGER       NROW, NCOL
      INTEGER       IDAY(DAYS), IMONTH(DAYS), IYEAR(DAYS)
      INTEGER       I, J, K, CLEN
      INTEGER       RES_DIRECT(NRESER_MAX,3)
      INTEGER       RESER(NCOL,NROW)
      INTEGER       H(NCOL,NROW)
      INTEGER       RULE, IRRIGATION
      INTEGER       NORESERVOIRS, NAM
      REAL          VOL(NRESER_MAX,DAYS),FLOWIN(NRESER_MAX,DAYS),FLOWOUT(NRESER_MAX,DAYS), VRESER(NRESER_MAX,DAYS)
      REAL          HHO(NRESER_MAX,DAYS),HYDROPOWER(NRESER_MAX,DAYS),HTK(NRESER_MAX,DAYS)
      REAL          SEEPAGE, INFILTRATION
      REAL          RES_EVAPORATION(NRESER_MAX,DAYS)
      REAL          HRESERMAX, HRESERMIN
      REAL          QRESER,VRESERE,HPHATDIEN,QSPILL,QTUR,VDEAD,VINI
      REAL          X1,X2,X3
      REAL          FLOW(DAYS)
      REAL          FACTOR_SUM
      CHARACTER*201 TEMPRPATH, ABC
      CHARACTER*72  RPATH, IPATH
      CHARACTER*20  CHUOI
      CHARACTER*5   NAME5
      CHARACTER*72  OUTPATH, RESERNAME
      LOGICAL   STEPBYSTEP

c     Subroutine main body
      CLEN = INDEX(OUTPATH,' ')-1
      DO I=1, NORESERVOIRS-1
            WRITE(CHUOI,*) RES_DIRECT(I,1)
            TEMPRPATH = trim(RPATH)//"res"//trim(ADJUSTL(CHUOI))//".txt"
            OPEN(25, FILE = TEMPRPATH,FORM = 'FORMATTED',STATUS='OLD')
            READ(25,*)
            READ(25,*) HRESERMAX, HRESERMIN, VRESERE, VDEAD, HPHATDIEN,QRESER, NAM, VINI, RESERNAME
            READ(25,*)
            READ(25,*) SEEPAGE, INFILTRATION
            READ(25,*)
            READ(25,*) IRRIGATION
            READ(25,*) IPATH
            READ(25,*)
            READ(25,*) RULE
            CLOSE(25)
            WRITE(CHUOI,*) RES_DIRECT(I,1)
            print*, 'Exporting data for Reservoir ',RESERNAME
            IF (STEPBYSTEP) THEN
                  print*, 'writing reservoir data for step-by-step version...'
                  STEPS=1         
                  IF (NDAY_SIM.EQ.1) THEN
                        ! Create a new file for step 1 and write the header
                        OPEN(40, FILE = OUTPATH(1:CLEN)//'/STEPBYSTEP/'//'reservoir_'
     &                 //trim(ADJUSTL(CHUOI))//'_'//trim(ADJUSTL(NAME5))//'.day')
                        IF ((RULE .EQ. 1) .OR. (RULE .EQ. 2) .OR. (RULE .EQ. 3) .OR. (RULE .EQ. 5)) THEN
                              WRITE(40,*) 'Year Month Day Volume_1000cm FullSupplyVolume WaterLevel_m Qinflow_cms '
     &                     //'Qspilways_cms Qturbine_cms Energy_MW'
                        ELSE
                              WRITE(40,*) 'Year Month Day Volume_1000cm FullSupplyVolume WaterLevel_m Qinflow_cms '
     &                        //'Qspilways_cms Qoutflow_cms Energy_MW'
                        ENDIF            
                  ELSE
                        ! Write to old file created from step 1 
                        OPEN(40, FILE = OUTPATH(1:CLEN)//'/STEPBYSTEP/'//'reservoir_'
     &                 //trim(ADJUSTL(CHUOI))//'_'//trim(ADJUSTL(NAME5))//'.day',status='old', position='append',ERR=9004)
                  ENDIF
            ELSE  ! Not a Step-By-Step Mode  
                  STEPS=DAYS    
                  OPEN(40, FILE = OUTPATH(1:CLEN)//'reservoir_'
     &            //trim(ADJUSTL(CHUOI))//'_'//trim(ADJUSTL(NAME5))//'.day')
                  IF ((RULE .EQ. 1) .OR. (RULE .EQ. 2) .OR. (RULE .EQ. 3) .OR. (RULE .EQ. 5)) THEN
                        WRITE(40,*) 'Year Month Day Volume_1000cm FullSupplyVolume WaterLevel_m Qinflow_cms '
     &                     //'Qspilways_cms Qturbine_cms Energy_MW'
                  ELSE
                        WRITE(40,*) 'Year Month Day Volume_1000cm FullSupplyVolume WaterLevel_m Qinflow_cms '
     &                        //'Qspilways_cms Qoutflow_cms Energy_MW'
                  ENDIF            
            ENDIF
            ! Write to output file: reservoir data (level,inflow,storage,outflow, hydropower)
            IF ((RULE .EQ. 1) .OR. (RULE .EQ. 2) .OR. (RULE .EQ. 3) .OR. (RULE .EQ. 5)) THEN
                  DO K = 1, STEPS
                        IF (FLOWOUT(I,K)>QRESER) THEN
                              QSPILL = FLOWOUT(I,K) - QRESER
                              QTUR = QRESER
                        ELSE
                              QSPILL = 0
                              QTUR = FLOWOUT(I,K)
                        END IF
                        IF (VOL(I,K)<0) THEN
                              VOL(I,K) = 0
                        END IF
                        IF (HHO(I,K)<0) THEN
                              HHO(I,K) = 0
                        END IF
                        IF (FLOWIN(I,K)<0) THEN
                              FLOWIN(I,K) = 0
                        END IF
                        IF (HYDROPOWER(I,K)<0) THEN
                              HYDROPOWER(I,K) = 0
                        END IF
                        
                        IF (STEPBYSTEP) THEN
                              WRITE(40,*) SIM_YEAR, SIM_MON, SIM_DAY, VOL(I,K), VRESER(I,K), HHO(I,K), FLOWIN(I,K), QSPILL, QTUR, HYDROPOWER(I,K) !This function also writes the maximum storage capacity (VRESER), writing hydropower output in the 7th column (important for optimization.py)
                        ELSE
                              WRITE(40,*) IYEAR(K),IMONTH(K),IDAY(K), VOL(I,K), VRESER(I,K), HHO(I,K), FLOWIN(I,K), QSPILL, QTUR, HYDROPOWER(I,K) !This function also writes the maximum storage capacity (VRESER), writing hydropower output in the 7th column (important for optimization.py)
                        ENDIF       
                  END DO
            ELSE  !Operation Strategy 4
                  QSPILL = 0        
                  DO K = 1, STEPS
                        IF (VOL(I,K)<0) THEN
                              VOL(I,K) = 0
                        END IF
                        IF (HHO(I,K)<0) THEN
                              HHO(I,K) = 0
                        END IF
                        IF (FLOWIN(I,K)<0) THEN
                              FLOWIN(I,K) = 0
                        END IF
                        IF (HYDROPOWER(I,K)<0) THEN
                              HYDROPOWER(I,K) = 0
                        END IF

                        IF (STEPBYSTEP) THEN
                              WRITE(40,*) SIM_YEAR, SIM_MON, SIM_DAY, VOL(I,K), VRESER(I,K), HHO(I,K), FLOWIN(I,K), QSPILL, FLOWOUT(I,K), HYDROPOWER(I,K)
                        ELSE
                              WRITE(40,*) IYEAR(K),IMONTH(K),IDAY(K), VOL(I,K), VRESER(I,K), HHO(I,K), FLOWIN(I,K), QSPILL, FLOWOUT(I,K), HYDROPOWER(I,K)
                        ENDIF
                  END DO
            END IF
            CLOSE(40)
      END DO
      RETURN
      STOP     
9004  WRITE(0,*) 'ERROR in Opening Output File; Please Check if Previous Time Step is Simulated'        
      END




C************************************************************************************************************************************************************************************
C     WRITE MONTHLY/YEAR DATA AT THE BASIN OUTLET
C************************************************************************************************************************************************************************************
      SUBROUTINE WRITE_MONTH
     & (DAYS_IN_MONTH,START_YEAR, STOP_YEAR, FIRST_YEAR, LAST_YEAR,
     &  START_MO, STOP_MO, FIRST_MO, LAST_MO,
     &  NAME5, DAYS, FLOW, FACTOR_SUM, MONTHLY, MONTHLY_mm,
     &  YEARLY, YEARLY_mm,OUTPATH,NDAY,IMONTH,IYEAR,MO,YR,NMONTHS,NYR)
      IMPLICIT NONE
      INTEGER DAYS_IN_MONTH(12)
      INTEGER NYR
      INTEGER START_YEAR, STOP_YEAR, FIRST_YEAR, LAST_YEAR
      INTEGER START_MO, STOP_MO, FIRST_MO, LAST_MO   !AWW-092700
      INTEGER DAYS,NDAY,NMONTHS
      INTEGER IMONTH(DAYS),IYEAR(DAYS)
      INTEGER SKIPTO, STOPAT
      INTEGER OLDMO
      INTEGER I, MONTH, YEAR, DAY_IN_MONTH
      INTEGER M, MCOUNT(12)     !AWW-092700
      INTEGER MO(12*NYR),YR(12*NYR)
c     INTEGER MNTH_INDX
      REAL    FLOW(DAYS)
      REAL    FACTOR_SUM
      REAL    MONTHLY(12*(STOP_YEAR-START_YEAR+1))
      REAL    MONTHLY_mm(12*(STOP_YEAR-START_YEAR+1))
      REAL    YEARLY(12)
      REAL    YEARLY_mm(12)
      CHARACTER*5  NAME5
      CHARACTER*72 OUTPATH, TMPPTH
c     concatenate output string
      I=INDEX(OUTPATH,' ')-1
C      OUTPATH(I:I+4)=NAME5
C      I=I+4
      OPEN(40, FILE = OUTPATH(1:I)//NAME5//'.month')
      OPEN(41, FILE = OUTPATH(1:I)//NAME5//'.month_mm')
      OPEN(42, FILE = OUTPATH(1:I)//NAME5//'.year')
      OPEN(43, FILE = OUTPATH(1:I)//NAME5//'.year_mm')
      OPEN(77, FILE = OUTPATH(1:I)//NAME5//'.end_of_month')
c     iniitalize monthly averages
      DO I = 1, 12*(STOP_YEAR-START_YEAR+1)
            MONTHLY(I) = 0.0
            MONTHLY_mm(I) = 0.0
      END DO
c     Average flows for each month in simulation
      M=1
      OLDMO=MO(1)
      DO I = 1, NDAY
            IF(IMONTH(I).ne.OLDMO) THEN
                  M=M+1
                  OLDMO=IMONTH(I)
            ENDIF
            MONTHLY(M)=MONTHLY(M)+FLOW(I)/DAYS_IN_MONTH(IMONTH(I))
            MONTHLY_mm(M) = MONTHLY_mm(M) + FLOW(I)/FACTOR_SUM
      END DO
C     writing monthly averages
      DO I = 1,12
            YEARLY(I) = 0.0
            YEARLY_mm(I) = 0.0
            MCOUNT(I) = 0
      END DO
c     Find months in time series to start and stop writing data
c     Note array starts at 1 regardless of actual month number
      SKIPTO = (FIRST_YEAR-START_YEAR)*12+(FIRST_MO-START_MO)+1
      STOPAT = NMONTHS-((STOP_YEAR-LAST_YEAR)*12+(STOP_MO-LAST_MO))
      DO I=SKIPTO,STOPAT
            WRITE(40,*) YR(I),MO(I), MONTHLY(I)
            WRITE(41,*) YR(I),MO(I), MONTHLY_mm(I)
            YEARLY(MOD(I-1,12)+1) = YEARLY(MOD(I-1,12)+1) + MONTHLY(I)
            YEARLY_mm(MOD(I-1,12)+1) =YEARLY_mm(MOD(I-1,12)+1) + MONTHLY_mm(I)
            MCOUNT(MO(I)) = MCOUNT(MO(I))+1
      END DO
      DO I = 1, 12
            IF(MCOUNT(I) .GT. 0) THEN
                  WRITE(42,*) I, YEARLY(I)/MCOUNT(I)
                  WRITE(43,*) I, YEARLY_mm(I)/MCOUNT(I)
            ELSE
                  WRITE(42,*) I, '  0'
                  WRITE(43,*) I, '  0'
            END IF
      END DO
      CLOSE(40)
      CLOSE(41)
      CLOSE(42)
      CLOSE(43)
      CLOSE(77)
      RETURN
      END
C     END OF FILE
C************************************************************************************************************************************************************************************
