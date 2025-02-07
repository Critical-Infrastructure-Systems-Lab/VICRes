      PROGRAM rout
c     Routing algorithm developed by D. Lohmann.
c     Code first developed by the University of Washington. See WA Hydrology Homepage for more details.
c     Modified by the Resilient Water Systems Group/Singapore University of Technology and Design to account for reservoir operations.
c     Reservoir presentation and operations are incorporated by the following steps
c     1. We first calculate the time lag from each cell to a reservoir / the basin outlet
c     2. We then calculate the flow to each reservoirs and the outlet
c     3. Track the reservoir network
c     4. Reservoir operations are modelled based on rule curves (1)/ operating rules (2)/ predefined time-series data (3)
      IMPLICIT NONE

c************************************************************************************************************************************************************************************
C     Declare variables
c************************************************************************************************************************************************************************************
c     RCS ID STRING
      CHARACTER*50 RCSID
c     DATA RCSID/"$Id: read_routines.f,v1.0 2019/04/17$"/
      INTEGER   IARGC
      INTEGER   isaleap
      EXTERNAL  isaleap
      INTEGER PMAX
c------------------------------------------------------------------------------------------------------------------------------------------------------------------            
C     USER DEFINED INPUTS
c------------------------------------------------------------------------------------------------------------------------------------------------------------------            
C     Change dimensions here:
C     1. NROW and NCOL should be larger than the grid
C     2. NYR should equal run length yrs+1
C     3. NSTATIONS_MAX and NRESER_MAX are the maximum number of stations and reservoirs that can be considered (change if you need to consider more of them)
C     4. PMAX is the total number of cells that can be routed to a specific station (need to change it for large basins)
      INTEGER NROW, NCOL, DAYS, NYR, NSTATIONS_MAX, NRESER_MAX
      PARAMETER (NROW = 400, NCOL = 300, NSTATIONS_MAX=10, NRESER_MAX=150)
      PARAMETER (NYR = 50)
      PARAMETER (PMAX = 14000)      ! Changed to accommodate number of cells in large basins (e.g., Mekong basin) 
      
C     -------------- No changes after here -------------------------------------------------------
C     -------------- Note: NRESER_MAX is the maximum number of reservoirs. Change if more reservoirs added
c------------------------------------------------------------------------------------------------------------------------------------------------------------------            
C     Unit-hydrograph parameters
c------------------------------------------------------------------------------------------------------------------------------------------------------------------            
      INTEGER   KE, LE, TMAX, UH_DAY
      REAL      DT
      PARAMETER (DAYS=NYR*366)
      PARAMETER (KE   = 12)
      PARAMETER (LE   = 48)
      PARAMETER (DT   = 3600.0)
      PARAMETER (UH_DAY = 96)
      PARAMETER (TMAX = UH_DAY*24)
      REAL      UH_BOX(PMAX,KE),UHM(NCOL,NROW,LE)
      REAL      UH_S(PMAX,KE+UH_DAY-1,NRESER_MAX)
      REAL      UH_DAILY(PMAX,UH_DAY)
      REAL      FR(TMAX,2)
      CHARACTER*80 UH_STRING(NSTATIONS_MAX)
c------------------------------------------------------------------------------------------------------------------------------------------------------------------
C     Routing
c------------------------------------------------------------------------------------------------------------------------------------------------------------------
      REAL      AREA(NSTATIONS_MAX)
      REAL      FACTOR_SUM
      REAL      XC,YC,SIZE
      REAL      FDUM
      REAL      VELO(NCOL,NROW), DIFF(NCOL,NROW)
      REAL      XMASK(NCOL,NROW), FRACTION(NCOL,NROW)
      REAL      BASE(DAYS,NSTATIONS_MAX), RUNO(DAYS,NSTATIONS_MAX), FLOW(DAYS,NSTATIONS_MAX)
      INTEGER   DIREC(NCOL,NROW,2)
      INTEGER   IDAY(DAYS), IMONTH(DAYS), IYEAR(DAYS)
      INTEGER   MO(12*NYR),YR(12*NYR)
      INTEGER   NO_OF_BOX(NRESER_MAX,NSTATIONS_MAX)
      INTEGER   NO_OF_BOX_STA(NSTATIONS_MAX)
      INTEGER   CATCHIJ(PMAX,2,NRESER_MAX,NSTATIONS_MAX)
      INTEGER   H(NCOL,NROW)
      INTEGER   PI(NSTATIONS_MAX),PJ(NSTATIONS_MAX),NR(NSTATIONS_MAX)
      INTEGER   PII,PJJ
      INTEGER   IROW,ICOL
      INTEGER   LP,M,Y,J,I,K,pp,uu,rr
      INTEGER   DPREC,FINAL
      INTEGER   NDAY,NDAY_SIM, CLEN
      INTEGER   NMONTHS
      LOGICAL   TORF
      LOGICAL   STEPBYSTEP, WRITE_OUTPUT
c------------------------------------------------------------------------------------------------------------------------------------------------------------------
C     Reservoir parameters
c------------------------------------------------------------------------------------------------------------------------------------------------------------------
      INTEGER   RESER(NCOL,NROW)          ! storing reservoir_location.txt file 
      INTEGER   NRESER                    ! actual number of reservoirs as read from the reservoirlocation.txt file 
      INTEGER   RESORDER(NSTATIONS_MAX,NRESER_MAX)  ! order of reservoir from US to DS as routing to the station under consideration 
      REAL      VRESER(NRESER_MAX,NSTATIONS_MAX,DAYS),VOL(NRESER_MAX,DAYS,NSTATIONS_MAX),FLOWIN(NRESER_MAX,DAYS),FLOWOUT(NRESER_MAX,DAYS,NSTATIONS_MAX)
      REAL      HHO(NRESER_MAX,DAYS),LUONGDIEN(NRESER_MAX,DAYS),HTK(NRESER_MAX,DAYS)
      INTEGER   NORESERVOIRS(NSTATIONS_MAX)     ! number of reservoirs contributing to the station under consideration (must be less than or equal NRESER)
      INTEGER   NSTATIONS                       ! actual number of stations as read from the stations.txt file      
      INTEGER   N,D
      INTEGER   RES_DIRECT(NRESER_MAX,3,NSTATIONS_MAX)    ! Reservoir Direction Array has 3 Layers: #Layer 1 stores the reservoir id (from reservoir location file) if the reservoir contributes to the station (from SEARCH_WHOLECATCHMENT in init_routines.f). #Layer 2 stores the downstream reservoir contributing to the reservoir under consideration. #Layer 3 stores the number of cells (NO_OF_BOX) contributes to the reservoir under consideration   
      REAL      RES_EVAPORATION(NRESER_MAX,DAYS)
      REAL      TRVLTIME(NRESER_MAX)
c------------------------------------------------------------------------------------------------------------------------------------------------------------------      
C     Filename
c------------------------------------------------------------------------------------------------------------------------------------------------------------------
      CHARACTER*21 NAME(NSTATIONS_MAX)
      CHARACTER*21  NAMERS51
      CHARACTER*21  NAMERS5(NSTATIONS_MAX,NRESER_MAX)
      CHARACTER*20 NAMERS(NSTATIONS_MAX,NRESER_MAX)
      CHARACTER*10  NAME5(NSTATIONS_MAX)
      CHARACTER*72 FILE_INPUT,FILENAME,RFILENAME
      CHARACTER*72 INPATH,OUTPATH,RPATH
c------------------------------------------------------------------------------------------------------------------------------------------------------------------
C     Variables for monthly means
c------------------------------------------------------------------------------------------------------------------------------------------------------------------
      INTEGER   DAYS_IN_MONTH(12)
      DATA      DAYS_IN_MONTH /31,28,31,30,31,30,31,31,30,31,30,31/
      INTEGER   START_YEAR,STOP_YEAR,FIRST_YEAR,LAST_YEAR
      INTEGER   START_MO,STOP_MO,FIRST_MO,LAST_MO
      INTEGER   SIM_YEAR,SIM_MON,SIM_DAY,START_YEAR_SBS,START_MON_SBS,START_DAY_SBS
      INTEGER   PREV_YEAR,PREV_MON,PREV_DAY
      REAL      MONTHLY(12*NYR)
      REAL      MONTHLY_mm(12*NYR)
      REAL      YEARLY(12)
      REAL      YEARLY_mm(12)
      NORESERVOIRS = 0    ! index of reservoir 
      NSTATIONS = 0
      
C************************************************************************************************************************************************************************************
C     OPEN NECESSARY FILES
C************************************************************************************************************************************************************************************
c     Process commandline args
      IF(IARGC().NE.1)THEN
           PRINT*, 'USAGE:  rout <infile>'
           STOP
      ENDIF
      CALL GETARG(1,FILE_INPUT)
      OPEN(1,FILE=FILE_INPUT,STATUS='OLD',ERR=9001)
      READ(1,'(//A)') FILENAME
      CALL READ_DIREC(DIREC,NCOL,NROW,H,XC,YC,SIZE,
     $     FILENAME,IROW,ICOL)
c     Process velocity file
      READ(1,*)
      READ(1,*) TORF
      IF(TORF)THEN
           READ(1,'(A)') FILENAME
           CALL READ_VELO(VELO,NCOL,NROW,FILENAME,IROW,ICOL)
      ELSE
           READ(1,*) FDUM
           CALL INIT_ARRAY(VELO,NCOL,NROW,FDUM)
      ENDIF
c     Process diffusion file
      READ(1,*)
      READ(1,*)TORF
      IF(TORF)THEN
          READ(1,'(A)') FILENAME
          CALL READ_DIFF(DIFF,NCOL,NROW,FILENAME,IROW,ICOL)
      ELSE
          READ(1,*) FDUM
          CALL INIT_ARRAY(DIFF,NCOL,NROW,FDUM)
      ENDIF
c     Process xmask file
      READ(1,*)
      READ(1,*)TORF
      IF(TORF)THEN
          READ(1,'(A)') FILENAME
          CALL READ_XMASK(XMASK,NCOL,NROW,FILENAME,IROW,ICOL)
      ELSE
          READ(1,*) FDUM
          CALL INIT_ARRAY(XMASK,NCOL,NROW,FDUM)
      ENDIF
c     Read fraction file
      READ(1,*)
      READ(1,*)TORF
      IF(TORF)THEN
          READ(1,'(A)') FILENAME
          CALL READ_FRACTION(FRACTION,NCOL,NROW,FILENAME,IROW,ICOL)
      ELSE
          READ(1,*) FDUM
          CALL INIT_ARRAY(FRACTION,NCOL,NROW,FDUM)
      ENDIF
c     Read station file
      READ(1,'(/A)')FILENAME
      OPEN(10,FILE=FILENAME)
c     Read input path and precision of VIC filenames
      READ(1,'(/A)')INPATH
      READ(1,*)DPREC
c     Read output pathname
      READ(1,'(/A)')OUTPATH
c     Read input path of reservoir information
      READ(1,'(/A)') RFILENAME
      READ(1,'(/A)') RPATH
c     Read input file name of reservoir locations
      CALL READ_RESE(RESER,ICOL,IROW,NCOL,NROW,RFILENAME,NRESER)
c     Number of days to process
      READ(1,*)
c     Start and end year/month from VIC simulation
      READ(1,*) START_YEAR, START_MO, STOP_YEAR, STOP_MO
c     Calculate number of days & months in simulation
      M=START_MO
      Y=START_YEAR
      NMONTHS = 0
      NDAY=0
      DO J=START_MO,12*(STOP_YEAR-START_YEAR)+STOP_MO
        IF(M.EQ.2) THEN
           LP=isaleap(Y)
        ELSE
           LP=0
        ENDIF
        NDAY = NDAY+DAYS_IN_MONTH(M)+LP
        NMONTHS = NMONTHS + 1
        MO(NMONTHS) = M
        YR(NMONTHS) = Y
        M = M + 1
        IF (M .GT. 12) THEN
            M = 1
            Y  = Y + 1
        ENDIF
      END DO
      IF(NDAY.GT.DAYS) THEN
         PRINT*, 'IN ROUT.F RESET DAYS TO ', NDAY
         STOP
      ENDIF
      PRINT*,'NDAY = ',NDAY, ' NMONTHS = ',NMONTHS
c     Read start and end year/month for writing output
      READ(1,*) FIRST_YEAR, FIRST_MO, LAST_YEAR, LAST_MO
c     Read uh file
      READ(1,'(/A)')FILENAME
c     Read simulation mode: Step-By-Step vs All Steps
      READ(1,*)
      READ(1,*) TORF
      IF(TORF)THEN
         STEPBYSTEP=.TRUE.
      ELSE
         STEPBYSTEP=.FALSE.
      ENDIF
c     VIC simulation day for StepByStep version
      READ(1,*)
      READ(1,*) START_YEAR_SBS, START_MON_SBS, START_DAY_SBS,SIM_YEAR, SIM_MON, SIM_DAY

c     Calculate number of simulation days
      M=START_MON_SBS
      Y=START_YEAR_SBS
      D=START_DAY_SBS
      NDAY_SIM=0
      IF (SIM_YEAR.GT.START_YEAR_SBS .OR.  SIM_MON.GT.START_MON_SBS) THEN
         DO J=START_MON_SBS,12*(SIM_YEAR-START_YEAR_SBS)+SIM_MON-1
            IF(M.EQ.2) THEN
               LP=isaleap(Y)
            ELSE
               LP=0
            ENDIF
            NDAY_SIM = NDAY_SIM+DAYS_IN_MONTH(M)+LP
            M = M + 1
            IF (M .GT. 12) THEN
               M = 1
               Y  = Y + 1
            ENDIF
         END DO
      ENDIF
      !Calculate number of days in last month
       NDAY_SIM=NDAY_SIM+SIM_DAY
      !Define the previous day
      IF (SIM_DAY.EQ.1 .AND. SIM_MON.GT.1)  THEN
           IF(SIM_MON.EQ.2) THEN
               LP=isaleap(Y)
           ELSE
               LP=0
           ENDIF
          PREV_DAY=DAYS_IN_MONTH(SIM_MON-1)+LP
          PREV_MON=SIM_MON-1
          PREV_YEAR=SIM_YEAR
      ElSEIF (SIM_DAY.EQ.1 .AND. SIM_MON.EQ.1) THEN
          PREV_DAY=31
          PREV_MON=12
          PREV_YEAR=SIM_YEAR-1
      ElSE
          PREV_DAY=SIM_DAY-1
          PREV_MON=SIM_MON
          PREV_YEAR=SIM_YEAR
      ENDIF
c     Check if writing outputs or not (important for step-by-step version when using it in the coupling environment and there are iterations for convergence) 
      READ(1,*)
      READ(1,*) TORF
      IF(TORF)THEN
            WRITE_OUTPUT=.TRUE.
      ELSE
            WRITE_OUTPUT=.FALSE.
      ENDIF



C************************************************************************************************************************************************************************************
C     START MODELLING
C************************************************************************************************************************************************************************************

c------------------------------------------------------------------------------------------------------------------------------------------------------------------            
c     STEP 1: Calculate impulse response function for grid cells
c------------------------------------------------------------------------------------------------------------------------------------------------------------------            
      CALL MAKE_UHM(UHM,VELO,DIFF,XMASK,NCOL,NROW,LE,DT,IROW,ICOL)
c------------------------------------------------------------------------------------------------------------------------------------------------------------------            
c     STEP 2: Read Station Info
c------------------------------------------------------------------------------------------------------------------------------------------------------------------                  
      I=1   ! loop over required stations
 100  CONTINUE
      READ(10,*,END=110) NR(I),NAME(I),PI(I),PJ(I),AREA(I)   ! Read stations.txt
      READ(10,'(A80)',END=110) UH_STRING(I)
      IF (NR(I) .EQ. 1) THEN           ! only if the station is 1 ==> Active
            print*,'--------------------'
            PRINT*, 'Routing station: ', trim(NAME(I)), ' [Active]' ,'-->','Unit Hydrograph File:',UH_STRING(I) 
            PRINT*, 'Location:',PI(I),',', PJ(I)
            !WRITE(*,'(I2,2X,A,I4,I4,G12.6)') NR(I), NAME(I), PI(I), PJ(I)
            print*,'--------------------'
            PI(I)=ICOL+1-PI(I)						!note, the arrays are flipped left to right
            NAME5(I) = NAME(I)
            NSTATIONS=NSTATIONS+1                                 !At the end of the loop NSTATIONS will be equal to the actual number of considered stations

c------------------------------------------------------------------------------------------------------------------------------------------------------------------            
c     STEP 3: Searching for Cells/Reservoirs contributing to Station
c------------------------------------------------------------------------------------------------------------------------------------------------------------------            
c     Look for cells, contributing to the station                
            PII=PI(I)
            PJJ=PJ(I)
            print*, 'Searching Catchment...'
            CALL SEARCH_WHOLECATCHMENT
     &           (PI(I),PJ(I),DIREC,NCOL,NROW,
     &            NO_OF_BOX(:,I),CATCHIJ(:,:,:,I),PMAX,IROW,ICOL,
     &            NORESERVOIRS(I), RES_DIRECT(:,:,I), RESER,NRESER_MAX)
      
            ! Outputs from SEARCH_WHOLECATCHMENT:
            ! # NO_OF_BOX(:,I): number of cells contributing to station
            ! # CATCHIJ(:,:,:,I): the location (i,j) of the contributing cells
            ! # NORESERVOIRS(I): number of reservoirs contributing to station
            ! # RES_DIRECT(:,:,I): the first layer of RES_DIRECT (:,1,I) storing the index (or the numbering) of contributing reservoirs (as read from the reservoirlocation.txt)          

c------------------------------------------------------------------------------------------------------------------------------------------------------------------            
c     STEP 4: Reading Pre-Defined UH File (path in the configuration file)
c------------------------------------------------------------------------------------------------------------------------------------------------------------------   
            print*, 'Reading grid unit hydrograph (UH)...'
c           Read a pre-defined UH grid
            CALL READ_GRID_UH
     &           (UH_BOX,KE,PMAX,NO_OF_BOX(:,I), CATCHIJ(:,:,:,I),FILENAME,NORESERVOIRS(I),NRESER_MAX)    ! reading UH file from ../../RoutingSetup/UH.all (defined in the configuration.txt)

c------------------------------------------------------------------------------------------------------------------------------------------------------------------            
c     STEP 5: Making Grid UH only for Reservoir Catchments (cells contributing to reservoirs and then to station under consideration)
c------------------------------------------------------------------------------------------------------------------------------------------------------------------ 
            IF (STEPBYSTEP .AND. (NDAY_SIM.GT.1)) THEN
                  print*, 'Running STEPBYSTEP Version...'
                  print*, 'No need to remake grid UH after step 1 (only define name of reservoir-station (NAMERS5))...'
                  D=1    ! iterate only when the cell is a reservoir
                  DO N = 1,NO_OF_BOX(NORESERVOIRS(I),I) ! loop over cells contributing to station
                        IF ((RESER(CATCHIJ(N,1,NORESERVOIRS(I),I),CATCHIJ(N,2,NORESERVOIRS(I),I)).GT.0) .AND. (RESER(CATCHIJ(N,1,NORESERVOIRS(I),I),CATCHIJ(N,2,NORESERVOIRS(I),I)).NE.9999)) THEN	
                              WRITE(NAMERS(I,D),*) RESER(CATCHIJ(N,1,NORESERVOIRS(I),I),CATCHIJ(N,2,NORESERVOIRS(I),I))
                              NAMERS5(I,RESORDER(I,D)) = 'RES'//trim(adjustl(NAMERS(I,D)))//"_"//trim(adjustl(NAME5(I)))
                              D=D+1  
                        ENDIF            
                  END DO
            ELSE      
                  print*, 'Making Grid UH for Reservoir Catchments...'
                  D=1  ! iterate only when the cell is a reservoir
                  CLEN = INDEX(OUTPATH,' ')-1
                  DO N = 1,NO_OF_BOX(NORESERVOIRS(I),I)    ! loop over cells contributing to station
                        IF ((RESER(CATCHIJ(N,1,NORESERVOIRS(I),I),CATCHIJ(N,2,NORESERVOIRS(I),I)).GT.0) .AND. (RESER(CATCHIJ(N,1,NORESERVOIRS(I),I),CATCHIJ(N,2,NORESERVOIRS(I),I)).NE.9999)) THEN				
c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx            
c     STEP 5-1: Searching for Cells contributing to Reservoirs
c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                              WRITE(NAMERS(I,D),*) RESER(CATCHIJ(N,1,NORESERVOIRS(I),I),CATCHIJ(N,2,NORESERVOIRS(I),I))
                              NAMERS51 = 'RES'//trim(adjustl(NAMERS(I,D)))//"_"//trim(adjustl(NAME5(I)))
                              CALL SEARCH_CATCHMENTRS(CATCHIJ(N,1,NORESERVOIRS(I),I),
     &                        CATCHIJ(N,2,NORESERVOIRS(I),I),DIREC,NCOL,NROW,NO_OF_BOX(:,I),CATCHIJ(:,:,:,I),PMAX,
     &                        IROW,ICOL,NORESERVOIRS(I),RES_DIRECT(:,:,I),RESER,0,SIZE,VELO,PI(I),PJ(I),NRESER_MAX)
c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 
c     STEP 5-2: Making Grid UH only for Reservoir Catchments
c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                              CALL MAKE_GRID_UH
     &                        (DIREC, NO_OF_BOX(:,I), UH_DAY, TMAX, PI(I), PJ(I), LE, UH_DAILY, KE,
     &                        CATCHIJ(:,:,:,I),UHM, FR, PMAX, NCOL, NROW, UH_BOX, UH_S,
     &                        UH_STRING(I),NAMERS51,NORESERVOIRS(I),RESER,
     &                        RESER(CATCHIJ(N,1,NORESERVOIRS(I),I),CATCHIJ(N,2,NORESERVOIRS(I),I)),
     &                        RES_DIRECT(:,:,I),RESORDER(I,D),0,NRESER_MAX,OUTPATH,CLEN,STEPBYSTEP)
                              
                              NAMERS5(I,RESORDER(I,D)) = 'RES'//trim(adjustl(NAMERS(I,D)))//"_"//trim(adjustl(NAME5(I)))
                              D=D+1
                        END IF               
                  END DO
c------------------------------------------------------------------------------------------------------------------------------------------------------------------            
c     STEP 6: Making Grid UH for the rest of the basin (cells other than those contributing to reservoirs) to the station under consideration--> FINAL=1
c------------------------------------------------------------------------------------------------------------------------------------------------------------------ 
                  print*, 'Making Grid UH for Non-Reservoir Catchments...'
                  CLEN = INDEX(OUTPATH,' ')-1
c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx            
c     STEP 6-1: Searching for Cells contributing to Station (but not to reservoirs)
c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                  CALL SEARCH_CATCHMENTRS(PI(I),PJ(I),
     &            DIREC,NCOL,NROW,NO_OF_BOX(:,I),CATCHIJ(:,:,:,I),PMAX,
     &            IROW,ICOL,NORESERVOIRS(I),RES_DIRECT(:,:,I),RESER,1,SIZE,VELO,PI(I),PJ(I),NRESER_MAX)
c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 
c     STEP 6-2: Making Grid UH only for Non-Reservoir Catchments
c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx            
                  CALL MAKE_GRID_UH
     &                 (DIREC, NO_OF_BOX(:,I), UH_DAY, TMAX, PI(I), PJ(I), LE, UH_DAILY, KE,
     &                 CATCHIJ(:,:,:,I),UHM, FR, PMAX, NCOL, NROW, UH_BOX, UH_S,
     &                 UH_STRING(I),NAME5(I),NORESERVOIRS(I),RESER,NORESERVOIRS(I),RES_DIRECT(:,:,I),RESORDER(I,D),1,NRESER_MAX,OUTPATH,CLEN,STEPBYSTEP)
            ENDIF    ! end of if conditional for stepbystep [Line 308]                              
      I=I+1          ! increment to next station
      ENDIF          ! end of if conditional for active cell [Line 272] 
      GOTO 100       ! Go back to line 100 (Step 2) for reading next station info                     
 110  CONTINUE

c------------------------------------------------------------------------------------------------------------------------------------------------------------------            
c     STEP 7: Flow generation for the required station (Convolution Method)
c------------------------------------------------------------------------------------------------------------------------------------------------------------------  
      print*, 'Making Convolution...'
      CALL MAKE_CONVOLUTIONRS
     &     (RPATH,RESER,NCOL, NROW, NO_OF_BOX, PMAX, DAYS,
     &     CATCHIJ, BASE, RUNO, FLOW, KE, UH_DAY, FRACTION,
     &     FACTOR_SUM,XC,YC,SIZE,DPREC,INPATH,ICOL,NDAY,
     &     IDAY,IMONTH,IYEAR, START_YEAR, START_MO, MO, YR, NYR, VOL,
     &     FLOWIN, FLOWOUT, HHO, LUONGDIEN,HTK,DIREC,IROW,
     &     PI,PJ,NORESERVOIRS,RES_DIRECT,RES_EVAPORATION,TRVLTIME,RESORDER,NAMERS5,NAME5,NSTATIONS_MAX,UH_S,VRESER, NRESER_MAX, NSTATIONS,STEPBYSTEP,NDAY_SIM,OUTPATH)   
c------------------------------------------------------------------------------------------------------------------------------------------------------------------            
c     STEP 8: Writing Data into Output Files
c------------------------------------------------------------------------------------------------------------------------------------------------------------------  
c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 
c     STEP 8-1: Write Flow of all Stations in one output file (OUTPUT.day)
c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx           
      IF (STEPBYSTEP .AND. WRITE_OUTPUT) THEN
            print*, 'writing data for step-by-step...',NDAY_SIM
            CLEN = INDEX(OUTPATH,' ')-1     
            IF (NDAY_SIM.EQ.1) THEN
                  ! Create a new file for step 1 and write the header
                  OPEN(80, FILE = OUTPATH(1:CLEN)//'/STEPBYSTEP/'//'OUTPUT.day')
                  OPEN(81, FILE = OUTPATH(1:CLEN)//'/STEPBYSTEP/'//'OUTPUT.day_mm')
                  WRITE(80,*) '       YEAR','        MONTH','        DAY ',NAME(1:NSTATIONS)
                  WRITE(81,*) '       YEAR','        MONTH','        DAY ',NAME(1:NSTATIONS)
            ELSE
                  ! Write to old file created from step 1 
                  OPEN(80, FILE = OUTPATH(1:CLEN)//'/STEPBYSTEP/'//'OUTPUT.day',status='old', position='append',ERR=9004)
                  OPEN(81, FILE = OUTPATH(1:CLEN)//'/STEPBYSTEP/'//'OUTPUT.day_mm',status='old', position='append',ERR=9004)
            END IF
            ! Write to output file: date and flow values
            WRITE(80,*) SIM_YEAR, SIM_MON, SIM_DAY, FLOW(NDAY_SIM, 1:NSTATIONS)
            WRITE(81,*) SIM_YEAR, SIM_MON, SIM_DAY, (FLOW(NDAY_SIM,J) / FACTOR_SUM,J = 1, NSTATIONS)
            CLOSE(80)
            CLOSE(81)
      ELSE     ! Not a Step-By-Step Mode           
            print*, 'writing data...'
            CLEN = INDEX(OUTPATH,' ')-1
            OPEN(80, FILE = OUTPATH(1:CLEN)//'OUTPUT.day')
            OPEN(81, FILE = OUTPATH(1:CLEN)//'OUTPUT.day_mm')
            WRITE(80,*) '       YEAR','        MONTH','        DAY ',NAME(1:NSTATIONS)
            WRITE(81,*) '       YEAR','        MONTH','        DAY ',NAME(1:NSTATIONS)
            DO I = 1,NDAY
                  WRITE(80,*) IYEAR(I), IMONTH(I), IDAY(I), FLOW(I, 1:NSTATIONS)
                  WRITE(81,*) IYEAR(I), IMONTH(I), IDAY(I), (FLOW(I,J)/FACTOR_SUM,J = 1, NSTATIONS)
            END DO
            CLOSE(80)
            CLOSE(81)
      END IF  ! end of if conditional for STEPBYSTEP (Line 389

! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 
!     STEP 8-2: Write outputs for each Station/Reservoir
! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 
      IF (WRITE_OUTPUT) THEN   ! only if writing output is true in the configuration file
      Do J=1,NSTATIONS
c     Writing discharge data (only for the last station considered in stations.txt)  ! Bruno 
            print*, 'writing data for station: ',NAME(J)
            CALL WRITE_DATA
     &           (FLOW(:,J), NDAY, NAME(J), FACTOR_SUM, OUTPATH,IDAY,IMONTH,IYEAR,
     &           VOL(:,:,J),FLOWIN, FLOWOUT(:,:,J), HHO, LUONGDIEN,HTK,RESER,NCOL,NROW,
     &           ICOL,IROW, RPATH, NORESERVOIRS(J),RES_DIRECT(:,:,J),RES_EVAPORATION,NO_OF_BOX(:,J),NRESER_MAX,STEPBYSTEP,NDAY_SIM,SIM_YEAR, SIM_MON, SIM_DAY)
c     Writing reservoirs data (only for the last station considered in stations.txt: if it corresponds to the basin outlet, it will automatically write the output for all the reservoirs)
            print*, 'writing data for reservoirs upstream of station: ',NAME(J)
            CALL WRITE_RES
     &           (FLOW(:,J), NDAY, NAME(J), FACTOR_SUM, OUTPATH,IDAY,IMONTH,IYEAR,
     &           VOL(:,:,J),FLOWIN, FLOWOUT(:,:,J), HHO, LUONGDIEN,HTK,RESER,NCOL,NROW,
     &           ICOL,IROW, RPATH, NORESERVOIRS(J),RES_DIRECT(:,:,J),RES_EVAPORATION,NO_OF_BOX(:,J),VRESER(:,J,:),NRESER_MAX,STEPBYSTEP,NDAY_SIM,SIM_YEAR, SIM_MON, SIM_DAY)
! c     Writing monthly output data
            IF (.not. STEPBYSTEP) THEN
                  CALL WRITE_MONTH
     &                 (DAYS_IN_MONTH,START_YEAR, STOP_YEAR, FIRST_YEAR, LAST_YEAR, START_MO, STOP_MO, FIRST_MO, LAST_MO,
     &                 NAME(J), DAYS, FLOW(:,J), FACTOR_SUM, MONTHLY, MONTHLY_mm,
     &                 YEARLY, YEARLY_mm, OUTPATH, NDAY, IMONTH, IYEAR, MO, YR, NMONTHS, NYR)
            ENDIF
      END DO
      ENDIF
      STOP     
 9001 WRITE(*,*) 'CANNOT OPEN: ', FILE_INPUT
 9004 WRITE(0,*) 'ERROR in Opening Output File; Please Check if Previous Time Step is Simulated'     
      END
C************************************************************************************************************************************************************************************
c USER DEFINED FUNCTIONS
C************************************************************************************************************************************************************************************
c     FUNCTION  ISALEAP
      integer function isaleap( iyr )
c     return 1 if a leap yr else 0
      if( (mod(iyr,4) .eq. 0 .and. mod(iyr,100) .ne.0)
     $                       .or. mod(iyr,400) .eq. 0) then
         isaleap = 1
      else
         isaleap = 0
      endif
      end
C     END OF FILE