c       Contain main subroutines of the routing component and reservoir operation
c       Build based on the old make_convolution file 
c       Note: currently the maximum number of reservoirs is NRESER_MAX. Modify the variable declaration parts to increase (if needed)
c       28 Dec 2020 - Correct the bug related to the dead storage - Thank Mr. Arivoli and his collegues

C************************************************************************************************************************************************************************************
C       Read diffusion file
C************************************************************************************************************************************************************************************
        SUBROUTINE READ_DIFF(DIFF,NCOL,NROW,FILENAME,IROW, ICOL)
c       Declare variables
        INTEGER NCOL,NROW,IROW,ICOL,I,J
        REAL DIFF(NCOL,NROW)
        CHARACTER*72 FILENAME
c       Subroutine main body
        print*,'READING DIFFUSIVITY FROM FILE:', FILENAME    
        OPEN(10, FILE = FILENAME,FORM = 'FORMATTED',STATUS='OLD',ERR=9001)
        DO I = 1,6													! Read file, skip header
            READ(10,*)
        END DO
        DO J = IROW,1,-1
            READ(10,*) (DIFF(I,J), I=ICOL,1,-1)
        END DO
        CLOSE(10)
        RETURN
 9001   WRITE(*,*)'CANNOT OPEN INPUT FILE IN READ_DIFF',FILENAME
        END

C************************************************************************************************************************************************************************************
C       Read fraction file
C************************************************************************************************************************************************************************************
        SUBROUTINE READ_FRACTION(FRACTION,NCOL,NROW,FILENAME,IROW,ICOL)
c       Declare variables
        INTEGER NCOL,NROW,ICOL,IROW,I,J
        REAL FRACTION(NCOL,NROW)
        CHARACTER*72 FILENAME
c       Subroutine main body
        OPEN(22, FILE = FILENAME,FORM = 'FORMATTED',STATUS='OLD',ERR=9001)
        DO I = 1,6													! Read file, skip header 
            READ(22,*)
        END DO
        DO J = IROW,1,-1
            READ(22,*) (FRACTION(I,J), I=ICOL,1,-1)
        END DO
        CLOSE(22)
        RETURN
 9001   WRITE(*,*) 'CANNOT OPEN INPUT FILE IN READ_FRACTION',FILENAME
        STOP
        END

C************************************************************************************************************************************************************************************
C       Read GRID_UH
C************************************************************************************************************************************************************************************
        SUBROUTINE READ_GRID_UH
     &    (UH_BOX,KE,PMAX,NOB,CATCHIJ,FILENAME,NORESERVOIRS, NRESER_MAX)
        IMPLICIT NONE
c       Declare variables
        INTEGER KE, PMAX, NRESER_MAX
        INTEGER NOB(NRESER_MAX)
        INTEGER CATCHIJ(PMAX,2,NRESER_MAX)
        INTEGER N,K,NORESERVOIRS
        REAL    UH_BOX(PMAX,KE)
        REAL    JUNK
        CHARACTER*72 FILENAME
c       Subroutine main body
        DO N = 1,NOB(NORESERVOIRS)
            OPEN(14,FILE = FILENAME,FORM = 'FORMATTED',
     $         STATUS='OLD', ERR=9001)
            DO K = 1,KE
               READ(14,*) JUNK, UH_BOX(N,K)
            END DO
            CLOSE(14)
        END DO
        RETURN
 9001   WRITE(*,*) 'CANNOT OPEN INPUT FILE IN GRID_UH',FILENAME
        STOP
        END

C************************************************************************************************************************************************************************************
C       Read VELOCITY FILE
C************************************************************************************************************************************************************************************
        SUBROUTINE READ_VELO(VELO,NCOL,NROW,FILENAME,
     $      IROW,ICOL)
c       Declare variables
        INTEGER NCOL,NROW,IROW,ICOL,I,J
        REAL VELO(NCOL,NROW)
        CHARACTER*72 FILENAME
c       Subroutine main body
        print*,'READING VELOCITY FROM FILE:', FILENAME    
        OPEN(10, FILE = FILENAME,FORM = 'FORMATTED',STATUS='OLD', ERR=9001)
        DO I = 1,6												! Read file, skip header 
            READ(10,*)
        END DO
        DO J = IROW,1,-1
            READ(10,*) (VELO(I,J), I=ICOL,1,-1)
        END DO
        CLOSE(10)
        RETURN
 9001   WRITE(*,*) 'CANNOT OPEN INPUT FILE IN  READ_VELO',FILENAME
        STOP
        END

C************************************************************************************************************************************************************************************
C       Read XMASK FILE
C************************************************************************************************************************************************************************************
        SUBROUTINE READ_XMASK(XMASK,NCOL,NROW,FILENAME,
     $      IROW,ICOL)
c       Declare variables
        INTEGER NCOL,NROW,ICOL,IROW,I,J
        REAL XMASK(NCOL,NROW)
        CHARACTER*72, FILENAME
c       Subroutine main body
        OPEN(10, FILE = FILENAME,FORM = 'FORMATTED',STATUS='OLD',ERR=9001)
        DO I = 1,6												! Read file, skip header 
            READ(10,*)
        END DO
        DO J = IROW,1,-1
            READ(10,*, END=20) (XMASK(I,J), I=ICOL,1,-1)
        END DO
        CLOSE(10)
        RETURN
 20     WRITE(*,*) 'REACHED END OF XMASK:  LAST ROW', j
 9001   WRITE(*,*) 'CANNOT OPEN INPUT FILE IN READ_XMASK'
        STOP
        END

C************************************************************************************************************************************************************************************
C       Read flow direction file
C************************************************************************************************************************************************************************************
        SUBROUTINE READ_DIREC(DIREC,NCOL,NROW,H,XC,YC,SIZE,FILENAME,IROW,ICOL)
        IMPLICIT NONE
c       Declare variables
        INTEGER NCOL,NROW,I,J,IROW,ICOL,IMISS
        INTEGER DIREC(NCOL,NROW,2)       ! direction array has 2 dimensions: 1 stores to ith index of downstream cell and 2 stores the jth index of downstream cell. Example: direction=1 (north) means direc(,,1)=i and direc(,,2)=j+1
        INTEGER H(NCOL,NROW)
        REAL XC, YC, SIZE
        CHARACTER*72 FILENAME
        CHARACTER*14 CDUM
c       Subroutine main body
        OPEN(10, FILE = FILENAME, FORM = 'FORMATTED',STATUS='OLD',ERR=9001)
        READ(10,*) CDUM, ICOL											! number of columms (flow direction matrix)
        READ(10,*) CDUM, IROW											! number of rows
        READ(10,*) CDUM, XC
        READ(10,*) CDUM, YC
        READ(10,*) CDUM, SIZE											! cell size
        READ(10,*) CDUM, IMISS
        IF(IROW.GT.NROW .OR. ICOL.GT.NCOL)THEN
            WRITE(*,*) 'Incorrect dimensions:'
            WRITE(*,*) 'Reset nrow and ncol in main to;',irow, icol
            STOP
        ENDIF
        DO J = IROW,1,-1
            
            READ(10,*) (H(I,J), I=ICOL,1,-1)
        END DO
        CLOSE(10)
C       Convert flow direction from 1 - 8 into grid-based cells
        DO I = 1, ICOL												! convert the flow direction matrix
            DO J = 1,IROW
                IF (H(I,J) .EQ. 0) THEN
                    DIREC(I,J,1) = 0
                    DIREC(I,J,2) = 0
                ELSE IF (H(I,J) .EQ. 1) THEN
                    DIREC(I,J,1) = I
                    DIREC(I,J,2) = J+1
                ELSE IF (H(I,J) .EQ. 2) THEN
                    DIREC(I,J,1) = I-1
                    DIREC(I,J,2) = J+1
                ELSE IF (H(I,J) .EQ. 3) THEN
                    DIREC(I,J,1) = I-1
                    DIREC(I,J,2) = J
                ELSE IF (H(I,J) .EQ. 4) THEN
                    DIREC(I,J,1) = I-1
                    DIREC(I,J,2) = J-1
                ELSE IF (H(I,J) .EQ. 5) THEN
                    DIREC(I,J,1) = I
                    DIREC(I,J,2) = J-1
                ELSE IF (H(I,J) .EQ. 6) THEN
                    DIREC(I,J,1) = I+1
                    DIREC(I,J,2) = J-1
                ELSE IF (H(I,J) .EQ. 7) THEN
                    DIREC(I,J,1) = I+1
                    DIREC(I,J,2) = J
                ELSE IF (H(I,J) .EQ. 8) THEN
                    DIREC(I,J,1) = I+1
                    DIREC(I,J,2) = J+1
                END IF
            END DO
        END DO
        RETURN
 9001   WRITE(*,*) 'CANNOT OPEN INPUT FILE IN READ_DIREC',FILENAME
        STOP
        END

C************************************************************************************************************************************************************************************
C       Estimate flows at reservoirs/outets based on the unit hydrograph and linearlized Saint Vernant methods
C************************************************************************************************************************************************************************************
        SUBROUTINE MAKE_CONVOLUTION
     &      (RPATH,RESER,NCOL, NROW, NO_OF_BOX, PMAX, DAYS, CATCHIJ,
     &      BASE, RUNO, FLOW, KE, UH_S, UH_DAY, FRACTION, FACTOR_SUM,
     &      XC, YC, SIZE, DPREC, INPATH,ICOL,NDAY,IDAY,IMONTH,IYEAR, START_YEAR, START_MO,
     &      MO, YR, NYR, VOL, FLOWIN, FLOWOUT, HHO, RESFLOWS, NO, resn, RES_EVAPORATION, NORESERVOIRS, NRESER_MAX)

        IMPLICIT NONE
c       Declare variables
        INTEGER     N, NO, resn, I, J, DAYS, NDAY, II, JJ, K, SODONG, NRESER_MAX
        INTEGER     NCOL,NROW,ICOL,PMAX,KE,UH_DAY
        INTEGER     NO_OF_BOX(NRESER_MAX)
        INTEGER     CATCHIJ(PMAX,2,NRESER_MAX)
        INTEGER     NYR, START_YEAR, START_MO
        INTEGER     RESER(NCOL,NROW)
        INTEGER     IDAY(DAYS), IMONTH(DAYS), IYEAR(DAYS)
        INTEGER     MO(12*NYR),YR(12*NYR)
        INTEGER     MISS_NUM
        INTEGER     CURRENTYEAR
        INTEGER     MONTH_OF_YEAR
        INTEGER     DPREC, CLEN
        INTEGER     OPEYEAR, NORESERVOIRS
        REAL        RES_EVAPORATION(NRESER_MAX,DAYS)
        REAL        UH_S(PMAX,KE+UH_DAY-1,NRESER_MAX)
        REAL        BASE(DAYS), RUNO(DAYS), FLOW(DAYS), AIRTEMP(DAYS), WINDS(DAYS), RHD(DAYS)
        REAL        FRACTION(NCOL,NROW)
        REAL        VRESER(NRESER_MAX), QRESER(NRESER_MAX)
        REAL        HRESERMAX(NRESER_MAX), HRESERMIN(NRESER_MAX)
        REAL        V0, FLOWINN
        REAL        PI, RERD, FACTOR, FACTOR_SUM
        REAL        DESIGNWL, CURRENTWL
        REAL        JLOC, ILOC, EVA_FACTOR
        REAL        XC, YC, SIZE
        REAL        AREA, AREA_SUM
        REAL        VOL(NRESER_MAX,DAYS), FLOWIN(NRESER_MAX,DAYS), FLOWOUT(NRESER_MAX,DAYS)
        REAL        HHO(NRESER_MAX,DAYS)
        REAL        K_CONST
        REAL        DUM1,DUM2,DUM3,DUM4,DUM5,DUM6,DUM7,DUM8,DUM9,DUM10,DUM11,DUM12,DUM13
        REAL        RESFLOWS(NRESER_MAX,DAYS)
        REAL        nN, Ra, deltagamma
        REAL        potentialevapo(DAYS)
        REAL        Ramatrix(12, 13)
        REAL        weightingmatrix(45), eamatrix(43), deltaTa4matrix(23)
        REAL        temp, ea, deltaTa4, ed, wind
        REAL        vic_eva_vege(DAYS), vic_trans_vege(DAYS), vic_eva_soil(DAYS)
        LOGICAL      TORF
        CHARACTER*20 RESNO
        CHARACTER*72 RPATH
        CHARACTER*100 TEMPRPATH
        CHARACTER*20 LOC
        CHARACTER*72 INPATH
        CHARACTER*5 NAMES
        PARAMETER   (RERD  = 6371229.0)  											!radius of earth in meters
c       Subroutine main body
C       Ra matrix (mm day-1) for the Penman formula
        Ramatrix = reshape((/1.4, 3.6, 7.0, 11.1, 14.6, 16.4, 15.6, 12.6, 8.5, 4.7, 2.0, 0.9,
     &                     3.7, 6.0, 9.2, 12.7, 15.5, 16.6, 16.1, 13.7, 10.4, 7.1, 4.4, 3.1,
     &                     6.2, 8.4, 11.1, 13.8, 15.9, 16.7, 16.3, 14.7, 12.1, 9.3, 6.8, 5.6,
     &                     8.1, 10.5, 12.8, 14.7, 16.1, 16.5, 16.2, 15.2, 13.5, 11.2, 9.1, 7.9,
     &                     10.8, 12.4, 14.0, 15.2, 15.7, 15.8, 15.8, 15.4, 14.4, 12.9, 11.3, 10.4,
     &                     12.8, 13.9, 14.8, 15.2, 15.0, 14.8, 14.9, 15.0, 14.8, 14.2, 13.1, 12.5,
     &                     14.6, 15.0, 15.2, 14.7, 13.9, 13.4, 13.6, 14.3, 14.9, 15.0, 14.6, 14.3,
     &                     15.9, 15.7, 15.1, 13.9, 12.5, 11.7, 12.0, 13.1, 14.4, 15.4, 15.7, 15.8,
     &                     16.8, 16.0, 14.5, 12.5, 10.7, 9.7, 10.1, 11.6, 13.6, 15.3, 16.4, 16.9,
     &                     17.2, 15.8, 13.5, 10.9, 8.6, 7.5, 7.9, 9.7, 12.3, 14.8, 16.7, 17.5,
     &                     17.3, 15.1, 12.2, 8.9, 6.4, 5.2, 5.6, 7.6, 10.7, 13.8, 16.5, 17.8,
     &                     16.9, 14.1, 10.4, 6.7, 4.1, 2.9, 3.4, 5.4, 8.7, 12.5, 16.0, 17.6,
     &                     16.5, 12.6, 8.3, 4.3, 1.8, 0.9, 1.3, 3.1, 6.5, 10.8, 15.1, 17.5 /),
     &                     shape(Ramatrix))	 
!       Weighting factor delta/gamma matrix for the Penman formula
        weightingmatrix = (/0.69, 0.71, 0.73, 0.76, 0.78, 0.80, 0.83, 0.85, 0.88, 0.91, 0.94, 0.97,
     &                     1.00, 1.03, 1.06, 1.09, 1.13, 1.16, 1.19, 1.23, 1.26, 1.30, 1.34, 1.38, 
     &                     1.42, 1.47, 1.51, 1.56, 1.60, 1.65, 1.69, 1.74, 1.79, 1.84, 1.89, 1.94,
     &                     1.99, 2.04, 2.11, 2.17, 2.23, 2.28, 2.35, 2.41, 2.48/)
!       Ea matrix
        eamatrix = (/4.58, 4.75, 4.93, 5.11, 5.30, 5.49, 5.69, 5.89, 6.10, 6.32, 6.54, 6.77, 7.02, 7.26,
     &                     7.52, 7.79, 8.05, 8.33, 8.62, 8.91, 9.21, 9.52, 9.84, 10.18, 10.52, 10.87, 11.23, 
     &                     11.61, 11.99, 12.38, 12.79, 13.13, 13.63, 14.08, 14.53, 15.00, 15.49, 15.97, 16.47,
     &                     17.53, 18.68, 19.80, 21.10/)
!       DeltaTa4matrix for the Penman formula
        deltaTa4matrix = (/11.0, 11.2, 11.4, 11.5, 11.7, 11.9, 12.0, 12.2, 12.3, 12.5, 12.7, 12.9, 13.1,
     &                     13.3, 13.5, 13.7, 13.9, 14.0, 14.2, 14.4, 14.6, 14.8, 15.0/)
!       MISS_NUM is the number of grid cell output files not found
        MISS_NUM=0
!       *** 0 <= K_CONST = 1.0
! ***   K_CONST smaller 1.0 makes it a simple linear storage
        K_CONST = 1.0
        PI = ATAN(1.0) * 4.0
        AREA_SUM   = 0.0
        FACTOR_SUM = 0.0
        OPEYEAR = 0
        CURRENTYEAR = START_YEAR
        DO I = 1,NDAY
            FLOW(I) = 0.0
            RES_EVAPORATION(NO,I) = 0.0
        END DO
!       Look for starting year
        IF (NO .NE. NORESERVOIRS) THEN
            WRITE(RESNO,*) resn
            TEMPRPATH = trim(RPATH)//"res"//trim(ADJUSTL(RESNO))//".txt"			! only calculate open surface water evaporation after the commision year
            OPEN(26, FILE = TEMPRPATH,FORM = 'FORMATTED', STATUS='OLD',ERR=9002)
            READ(26,*)
            READ(26,*) DUM1,DUM2,DUM3,DUM4,DUM5,DUM6, OPEYEAR
            CLOSE(26)
        END IF
        ! Added to avoid using values from previous time step when having headwater cell (no upstream cells and so no convolution needed) [HD July 2024]
        IF (NO_OF_BOX(NO).EQ.0) THEN
            print*,'Headwater Cell or Souce Cell: No Upstream Cells....'
            print*, 'No Convolution Needed...'
            DO I = 1,NDAY
                RUNO(I) = 0.0
                BASE(I) = 0.0
                ! Adjust Runoff and Basefloww with cell values
                FLOW(I) = FLOW(I)+(BASE(I-J+1)+RUNO(I-J+1))
                RESFLOWS(NO,I) = 0
                FLOWIN(NO,I)=0
                FLOWOUT(NO,I)=0
            END DO
        ENDIF

!       Calculate the area of each cell        
        DO N = 1,NO_OF_BOX(NO) 
            DO I = 1,NDAY
                RUNO(I) = 0.0
                BASE(I) = 0.0
            END DO
            II = CATCHIJ(N,1,NO)
            JJ = CATCHIJ(N,2,NO)
!       the grid has been flipped left to right
!       find the revised cooordinates
            ILOC=XC + (ICOL-II)*SIZE + SIZE/2.0
            JLOC=YC + JJ*SIZE - SIZE/2.0
!       CONVERSIONFACTOR for mm/day to ft**3/sec
        AREA =  RERD**2*ABS(SIZE)*PI/180*										! give area of box in square meters
     &          ABS(SIN((JLOC-SIZE/2.0)*PI/180)-
     $          SIN((JLOC+SIZE/2.0)*PI/180))
        AREA_SUM = AREA_SUM + AREA
        FACTOR = FRACTION(II,JJ)*35.315*AREA/(86400.0*1000.0)  				!convert to sq.mi. by cell fract (original) - later, we convert back to SI unit
        FACTOR_SUM = FACTOR_SUM + FACTOR
        call create_vic_names(jloc,iloc,loc,clen,dprec)
        INQUIRE(FILE=INPATH(1:(INDEX(INPATH,' ')-1))//LOC(1:CLEN),
     $      EXIST=TORF)  
        
        IF (torf) THEN
            OPEN(20,FILE=INPATH(1:(INDEX(INPATH,' ')-1))//
     $      LOC(1:CLEN),
     $      STATUS='OLD',ERR=9001)
            DO I = 1,NDAY 
                READ(20,*,END=9001,ERR=9001) IYEAR(I), IMONTH(I), IDAY(I),
     &              DUM1, DUM2, RUNO(I), BASE(I), DUM3, DUM4, DUM5, DUM6, DUM7, vic_eva_vege(I),
     &              vic_trans_vege(I), vic_eva_soil(I), DUM8, DUM9, DUM10, DUM11, DUM12, RHD(I), DUM13,
     &              AIRTEMP(I),WINDS(I) 												! (IMPORTANT) modify here accordingly depending on your selected outputs from rainfall-runoff
!       check to be sure dates in VIC file start at same time specified in input file
                IF(I.eq.1) THEN
                    IF(IYEAR(I).ne.YR(1) . or. IMONTH(I).ne.MO(1)) THEN
                        print*, 'VIC output file does not match specified period in input file'
                        stop
                    ENDIF
                ENDIF
            END DO
        ELSE
            print*, INPATH(1:(INDEX(INPATH,' ')-1))//LOC(1:CLEN),' NOT FOUND, INSERTING ZEROS'
            miss_num = miss_num+1
            DO i=1,nday
                IYEAR(I)=9999
                IMONTH(I)=99
                IDAY(I)=99
                runo(i)=0
                base(i)=0
                potentialevapo(i) = 0
            END DO
        ENDIF
        IF (RESER(II,JJ) .EQ. 9999) THEN
        ! Calculate evaporation via water surface according to Penman formua (see Penman, 1948)
        ! Note: the vic rainfall-runoff uses equations in Handbook of Hydrology (Penman-Monteith evapotranspiration)
        ! E0 = (delta/gamma*H+Ea)/(delta/gamma+1)
        ! H = (1-r)Rin-R0
        ! (1-r)Rin = 0.95Ra(0.18*0.55n/N)
        ! n/N is represented by nN which is the ration between actual sunshine hours and possible sunshine hours
            nN = 0.6 															! modify if needed
            DO I = 1,NDAY
                ! Ra is the solar radian (depending on location, see J.S.G. McCulloch, 1965. E. African Forest. J. (3) 286-295 for appropiate values
                MONTH_OF_YEAR = INT(MOD(NDAY,365)/30) + START_MO	! appropiate (ignore 29,31, not so important in this case)
                CURRENTYEAR = START_YEAR + INT(MOD(NDAY,365))
                IF (CURRENTYEAR>=OPEYEAR) THEN
                    IF (MONTH_OF_YEAR>12) THEN
                        MONTH_OF_YEAR = 12
                    END IF
                    IF (INT(JLOC/10)>60) THEN
                        Ra = Ramatrix(MONTH_OF_YEAR,1)
                    ELSE IF ((INT(JLOC/10)<=60) .AND. (INT(ILOC/10)>=0)) THEN
                        Ra = Ramatrix(MONTH_OF_YEAR,(7-INT(JLOC/10)))
                    ELSE IF (INT(JLOC/10)<-60) THEN
                        Ra = Ramatrix(MONTH_OF_YEAR,13)
                    ELSE 
                        Ra = Ramatrix(MONTH_OF_YEAR,(7-INT(JLOC/10)))
                    END IF
                    ! delta/gamma is depedent on temperature
                    temp = AIRTEMP(I)        										! (degree celcius)
                    wind = WINDS(I) * 24 * 3600 / 1609       						! (miles/day)
                    IF (temp<0) THEN
                        deltagamma = 0.69											! could be modifed (using intepolation)
                    ELSE IF (temp>22) THEN
                        deltagamma = 2.48											! could be modifed (using intepolation)
                    ELSE
                        deltagamma = weightingmatrix(1+int(temp*2))
                    END IF
                    ! calculate the saturation vapor pressure ea
                    IF (temp<0) THEN
                        ea = 4.4													! could be modifed (using intepolation)
                    ELSE if (temp>22) THEN
                        ea = 21.1													! could be modifed (using intepolation)
                    ELSE
                        ea = eamatrix(1+int(temp*2))
                    END IF
                    ! calculate deltaTa4
                    IF (temp<-1) THEN												! minimum value in the table is -1o Celcius
                       deltaTa4 = 11.0												! could be modifed (using intepolation)
                    ELSE IF (temp>21) THEN 											! maximum value in the table is 21o Celcius
                       deltaTa4 = 15.0												! could be modifed (using intepolation)
                    ELSE
                       deltaTa4 = deltaTa4matrix(int(temp)+2)
                    END IF
                    ! actual vapor pressure
                    ed = ea * RHD(I) / 100
                    ! potential evaporation 
                    potentialevapo(I) = (deltagamma*(0.95*(0.18+0.55*nN)*Ra
     &                 - deltaTa4*(0.56-0.09*sqrt(ed))*(0.10+0.90*nN))
     &                 + 0.35*(0.5+wind/100)*(ea-ed))/(deltagamma+1)
                    ! 0.95 = 1.0 - 0.05 (albedo for water)
                    ! convert potential to actual evaporation by multiplying a factor of K (modify if needed)
                    EVA_FACTOR = 0.9
                    RES_EVAPORATION(NO,I) = RES_EVAPORATION(NO,I)+potentialevapo(I)*EVA_FACTOR	! this is the total water loss due to evaporation in a reservoir (for checking only)
                END IF
            END DO
        ELSE
              potentialevapo=0  ! We are not in a reservoir cell, we do not have evaporation
        END IF  
        
        ! Calculate baseflow and runoff
        DO I = 1,NDAY
            RUNO(I) = RUNO(I) * FACTOR * 0.028 			 ! 0.028 convert from cfs into cms
            BASE(I) = BASE(I) * FACTOR * 0.028
        END DO

        DO I = 1,NDAY
            DO J = 1,KE+UH_DAY-1
                IF ((I-J+1) .GE. 1) THEN
                        FLOW(I) = FLOW(I)+UH_S(N,J,NO)*(BASE(I-J+1)+RUNO(I-J+1))
                END IF   
            END DO 
            ! allow negative values - evaporation is higher than inflow
            IF (potentialevapo(I)*EVA_FACTOR>(vic_eva_vege(I)+vic_eva_soil(I)+vic_trans_vege(I))) THEN
                ! also consider the water losses due to evapotranspiration calculated by VIC rainfall-runoff
                FLOW(I) = FLOW(I)-potentialevapo(I)*FACTOR*0.028/1000*EVA_FACTOR 
     &          + (vic_eva_soil(I)+vic_eva_vege(I)+vic_trans_vege(I))*FACTOR*0.028/1000
            END IF
            IF (potentialevapo(I)<0.0001) THEN          ! if too small, set to 0, also avoid negative values
                potentialevapo(I) = 0
            END IF
            RESFLOWS(NO,I) = FLOW(I)
           
               
        END DO
        CLOSE(20)
        END DO
       
        IF (MISS_NUM>0) THEN
            print*, MISS_NUM, ' files not found, zero runoff/baseflow used'
        ENDIF
        RETURN
 9001   WRITE(*,*) 'Error reading time-series data, ',
     $              'insufficient data or missing input file',
     $               INPATH(1:INDEX(INPATH,' ')-1)//LOC(1:CLEN)
 9002   WRITE(*,*) 'Error in reading reservoir data'
        END

C************************************************************************************************************************************************************************************
C       Determine the distance/time from the downstream cell of a reservoir to its cascade reservoir
C************************************************************************************************************************************************************************************
        SUBROUTINE CALCULATE_DISTANCE(SI,SJ,ICOL,IROW,NCOL,NROW,DIREC,TRVLTIME,
     &  RESORDER,RESER,FINAL,NORESERVOIRS,SIZE_CELL,VELO,FI,FJ,NRESER_MAX)
        IMPLICIT NONE
c       Declare variables  
        INTEGER SI,SJ,II,JJ,COUNTCELL,FINAL,NCOL,NROW,NORESERVOIRS,FI,FJ,NRESER_MAX
        INTEGER DIREC(NCOL,NROW,2)
        INTEGER RESER(NCOL,NROW)
        INTEGER I,J,ICOL,IROW,III,JJJ,RESORDER
        REAL    VELO(NCOL,NROW)
        REAL    TRVLTIME(NRESER_MAX)
        REAL    SIZE_CELL
c       Subroutine main body
        II = SI
        JJ = SJ
        COUNTCELL = 0
 400    CONTINUE
        IF ((II .GT. ICOL) .OR. (II .LT.1) .OR.
     &      (JJ .GT. IROW) .OR. (JJ .LT.1)) THEN
            GOTO 410
        END IF
        IF (((RESER(II,JJ)>0) .AND. (RESER(II,JJ) .NE. 9999))
     &     .OR. (II .EQ. FI) .AND. (JJ .EQ. FJ) .AND. (FINAL .EQ. 1)) THEN
            GOTO 410
        ELSE
            IF ((DIREC(II,JJ,1).NE.0) .AND.    										!check if the current
     &         (DIREC(II,JJ,2) .NE.0)) THEN   											!ii,jj cell routes down
                III = DIREC(II,JJ,1)         										!to the subbasin outlet
                JJJ = DIREC(II,JJ,2)         										!point, following the
                II  = III                    										!direction of direc(,)
                JJ  = JJJ                    										!from each cell
                COUNTCELL = COUNTCELL + 1											!this is for checking purposes, could be deleted
                IF (FINAL .EQ. 0) THEN
                    TRVLTIME(RESORDER) = TRVLTIME(RESORDER) + SIZE_CELL*1000*116/VELO(SI,SJ)	! +1 for the current cell (116km = 1 degree = 116000m; ignore the distortion of cells in the upper and lower lattitude)
                  ELSE 
                    TRVLTIME(NORESERVOIRS) = TRVLTIME(NORESERVOIRS) + SIZE_CELL*1000*116/VELO(SI,SJ)			! +1 for the current cell (116km = 1 degree = 116000m; ignore the distortion of cells in the upper and lower lattitude)
                END IF
                GOTO 400
            END IF                           											!if you get there,
        END IF                                										!no_of_box increments
 410    CONTINUE       
        TRVLTIME(RESORDER) = TRVLTIME(RESORDER)/24/3600								! convert to traveling time in number of days
        RETURN
        END

C************************************************************************************************************************************************************************************
C       Determine which cells contribute to flows at each reservoirs
C************************************************************************************************************************************************************************************
        SUBROUTINE SEARCH_CATCHMENTRS
     &     (PI,PJ,DIREC,NCOL,NROW,NO_OF_BOX,CATCHIJ,PMAX,
     $      IROW,ICOL,NORESERVOIRS,RES_DIRECT,RESER,FINAL,SIZE_CELL,VELO,FI,FJ,NRESER_MAX,TRVLTIME)
        IMPLICIT NONE
c       Declare variables
        INTEGER PI,PJ,I,J,NCOL,NROW,PMAX,ICOL,IROW,NRESER_MAX
        INTEGER II, JJ, III, JJJ, K,FI,FJ
        INTEGER DIREC(NCOL,NROW,2)
        INTEGER NO_OF_BOX(NRESER_MAX)
        INTEGER CATCHIJ(PMAX,2,NRESER_MAX)
        INTEGER RESER(NCOL,NROW)
        INTEGER RES_DIRECT(NRESER_MAX,3)
        INTEGER NORESERVOIRS,RESORDER
        INTEGER FINAL
        REAL    SIZE_CELL
        REAL    TRVLTIME(NRESER_MAX)
        REAL    VELO(NCOL,NROW)
c       Subroutine main body
        RESORDER = NORESERVOIRS   ! station is not at reservoir location (station is downstream of all its reservoirs = NORESERVOIRS)
        DO I = 1, (NORESERVOIRS-1)    ! loop over reservoirs only [the last NORESERVOIRS is the station under consideration downstream of reservoirs]
            IF (RESER(PI,PJ).EQ. RES_DIRECT(I,1)) THEN   ! PI,PJ here is the reservoir location [PI=CATCHIJ(N,1,NORESERVOIRS(I),I),PJ=CATCHIJ(N,2,NORESERVOIRS(I),I)]
               RESORDER = I      ! find the order of the reservoir 
            END IF
        END DO
        
        NO_OF_BOX(RESORDER)=0
        If (FINAL .EQ. 1) THEN    ! this is the case of the station under consideration (not a real reservoir). From SEARCH_WHOLECATCHMENT; NORESERVOIRS=NORESERVOIRS+1 after looping over all reservoirs, i.e., the last reservoir is the station
            RESORDER = NORESERVOIRS
        END IF
        ! This part adds cell to bags (each bag represents one reservoir)
        DO I = 1, ICOL
            DO J = 1, IROW
                II = I
                JJ = J
 300            CONTINUE 
                IF ((II .GT. ICOL) .OR. (II .LT.1) .OR.
     &              (JJ .GT. IROW) .OR. (JJ .LT.1)) THEN
                        GOTO 320
                END IF
                IF ((II .EQ.  PI) .AND. (JJ .EQ. PJ)) THEN
                    NO_OF_BOX(RESORDER) = NO_OF_BOX(RESORDER) + 1
                    CATCHIJ(NO_OF_BOX(RESORDER),1,RESORDER) = I
                    CATCHIJ(NO_OF_BOX(RESORDER),2,RESORDER) = J
                    !WRITE(*,*) 'Cell (',I,',',J,')---> Res ID ',RES_DIRECT(RESORDER,1)				! show which cells go to which reservoir
                    IF (FINAL .EQ. 0) THEN
                        GOTO 310 
                    ELSE
                        GOTO 320
                    END IF
                ELSE IF ((FINAL .EQ. 1) .AND. (RESER(II,JJ)>0)) THEN   ! reservoir over the grid and we target outlet station, then skip routing this grid cell to the station, why?? This will reduce number of upstream cells and so change unit hydrograph (because reservoirs are already dealed separtaely and no need to duplicate their cells)
                    GOTO 320
                ELSE IF ((RESER(II,JJ)>0) .AND.(RESER(II,JJ) .NE. 9999)) THEN   ! reservoir over the grid and we target a reservoir grid, then check routing this reservoir to downstream cell
                    GOTO 310
                ELSE
                    IF ((DIREC(II,JJ,1).NE.0) .AND.    								!check if the current
     &                  (DIREC(II,JJ,2) .NE.0)) THEN   								!ii,jj cell routes down
                        III = DIREC(II,JJ,1)         								!to the subbasin outlet
                        JJJ = DIREC(II,JJ,2)         								!point, following the
                        II  = III                    								!direction of direc(,)
                        JJ  = JJJ                    								!from each cell
                        GOTO 300
                    END IF   														!if you get there,
                END IF                                								!no_of_box increments
 310            CONTINUE                              								!and you try another  
                IF ((I .GE. 2) .AND. (I .LT. (ICOL-1)).AND. (J .GE. 2)				! this part to track which cell a reservoir releases to
     &              .AND. (J .LT. (IROW-1))) THEN
                    
                    IF ((RESER(I-1,J) .GT. 0) .AND. (DIREC(I-1,J,1) .EQ. I)
     &                  .AND. (DIREC(I-1,J,2) .EQ. J) .AND. (RESER(I-1,J) .NE. 9999)) THEN
                            DO K = 1, NORESERVOIRS
                                IF (RESER(I-1,J) .EQ. RES_DIRECT(K,1)) THEN
                                    RES_DIRECT(K,3) =  NO_OF_BOX(RESORDER)
                                    CALL CALCULATE_DISTANCE(I,J,ICOL,IROW,NCOL,NROW,DIREC,TRVLTIME,
     &                                                      RESORDER,RESER,FINAL,NORESERVOIRS,SIZE_CELL,VELO,FI,FJ,NRESER_MAX)
                                END IF
                            END DO
                    END IF

                    IF ((RESER(I+1,J) .GT.0) .AND. (DIREC(I+1,J,1) .EQ. I)
     &                  .AND. (DIREC(I+1,J,2) .EQ. J) .AND. (RESER(I+1,J) .NE.9999)) THEN
                        DO K = 1, NORESERVOIRS
                            IF (RESER(I+1,J) .EQ. RES_DIRECT(K,1)) THEN
                                RES_DIRECT(K,3) =  NO_OF_BOX(RESORDER)
                                CALL CALCULATE_DISTANCE(I,J,ICOL,IROW,NCOL,NROW,DIREC,TRVLTIME,
     &                                                  RESORDER,RESER,FINAL,NORESERVOIRS,SIZE_CELL,VELO,FI,FJ,NRESER_MAX)
                            END IF
                        END DO
                    END IF

                    IF ((RESER(I,J-1) .GT.0) .AND. (DIREC(I,J-1,1) .EQ. I)
     &                  .AND. (DIREC(I,J-1,2) .EQ. J) .AND. (RESER(I,J-1) .NE.9999)) THEN
                        DO K = 1, NORESERVOIRS
                            IF (RESER(I,J-1) .EQ. RES_DIRECT(K,1)) THEN
                                RES_DIRECT(K,3) =  NO_OF_BOX(RESORDER)
                                CALL CALCULATE_DISTANCE(I,J,ICOL,IROW,NCOL,NROW,DIREC,TRVLTIME,
     &                                                  RESORDER,RESER,FINAL,NORESERVOIRS,SIZE_CELL,VELO,FI,FJ,NRESER_MAX)
                            END IF
                        END DO
                    END IF

                    IF ((RESER(I,J+1) .GT.0) .AND. (DIREC(I,J+1,1) .EQ. I)
     &                  .AND. (DIREC(I,J+1,2) .EQ. J) .AND. (RESER(I,J+1) .NE.9999)) THEN
                        DO K = 1, NORESERVOIRS
                            IF (RESER(I,J+1) .EQ. RES_DIRECT(K,1)) THEN
                                RES_DIRECT(K,3) =  NO_OF_BOX(RESORDER)
                                CALL CALCULATE_DISTANCE(I,J,ICOL,IROW,NCOL,NROW,DIREC,TRVLTIME,
     &                                                  RESORDER,RESER,FINAL,NORESERVOIRS,SIZE_CELL,VELO,FI,FJ,NRESER_MAX)
                            END IF
                        END DO
                    END IF

                    IF ((RESER(I-1,J-1) .GT.0) .AND. (DIREC(I-1,J-1,1) .EQ. I)
     &                  .AND. (DIREC(I-1,J-1,2) .EQ. J) .AND. (RESER(I-1,J-1) .NE.9999)) THEN
                        DO K = 1, NORESERVOIRS
                            IF (RESER(I-1,J-1) .EQ. RES_DIRECT(K,1)) THEN
                                RES_DIRECT(K,3) =  NO_OF_BOX(RESORDER)
                                CALL CALCULATE_DISTANCE(I,J,ICOL,IROW,NCOL,NROW,DIREC,TRVLTIME,
     &                                                  RESORDER,RESER,FINAL,NORESERVOIRS,SIZE_CELL,VELO,FI,FJ,NRESER_MAX)
                            END IF
                        END DO
                    END IF

                    IF ((RESER(I-1,J+1) .GT.0) .AND. (DIREC(I-1,J+1,1) .EQ. I)
     &                  .AND. (DIREC(I-1,J+1,2) .EQ. J) .AND. (RESER(I-1,J+1) .NE.9999)) THEN
                        DO K = 1, NORESERVOIRS
                            IF (RESER(I-1,J+1) .EQ. RES_DIRECT(K,1)) THEN
                                RES_DIRECT(K,3) =  NO_OF_BOX(RESORDER)	
                                CALL CALCULATE_DISTANCE(I,J,ICOL,IROW,NCOL,NROW,DIREC,TRVLTIME,
     &                                                  RESORDER,RESER,FINAL,NORESERVOIRS,SIZE_CELL,VELO,FI,FJ,NRESER_MAX)
                            END IF
                        END DO
                    END IF

                    IF ((RESER(I+1,J-1) .GT.0) .AND. (DIREC(I+1,J-1,1) .EQ. I)
     &                  .AND. (DIREC(I+1,J-1,2) .EQ. J) .AND. (RESER(I+1,J-1) .NE.9999)) THEN
                        DO K = 1, NORESERVOIRS
                            IF (RESER(I+1,J-1) .EQ. RES_DIRECT(K,1)) THEN
                                RES_DIRECT(K,3) =  NO_OF_BOX(RESORDER)
                                CALL CALCULATE_DISTANCE(I,J,ICOL,IROW,NCOL,NROW,DIREC,TRVLTIME,
     &                                                  RESORDER,RESER,FINAL,NORESERVOIRS,SIZE_CELL,VELO,FI,FJ,NRESER_MAX)
                            END IF
                        END DO
                    END IF

                    IF ((RESER(I+1,J+1) .GT.0) .AND. (DIREC(I+1,J+1,1) .EQ. I)
     &                  .AND. (DIREC(I+1,J+1,2) .EQ. J) .AND. (RESER(I+1,J+1) .NE.9999)) THEN
                        DO K = 1, NORESERVOIRS
                            IF (RESER(I+1,J+1) .EQ. RES_DIRECT(K,1)) THEN
                                RES_DIRECT(K,3) =  NO_OF_BOX(RESORDER)
                                CALL CALCULATE_DISTANCE(I,J,ICOL,IROW,NCOL,NROW,DIREC,TRVLTIME,
     &                                                  RESORDER,RESER,FINAL,NORESERVOIRS,SIZE_CELL,VELO,FI,FJ,NRESER_MAX)
                            END IF
                        END DO
                    END IF
                END IF
320             CONTINUE
            END DO
        END DO
        
        ! Find downstream reservoirs
        II = PI
        JJ = PJ
500     CONTINUE				! this part is to locate the sequence of reservoirs
        IF ((II .GT. ICOL) .OR. (II .LT.1) .OR.
     &      (JJ .GT. IROW) .OR. (JJ .LT.1)) THEN
            RES_DIRECT(RESORDER,2) = 0
            GOTO 510
        END IF
        IF ((RESER(II,JJ)>0) .AND. (RESER(II,JJ) .NE. 9999)     
     &     .AND. ((II .NE. PI) .OR. (JJ .NE. PJ))) THEN
            RES_DIRECT(RESORDER,2) = RESER(II,JJ)
            GOTO 510
        ELSE
            IF ((DIREC(II,JJ,1).NE.0) .AND.    										!check if the current
     &          (DIREC(II,JJ,2) .NE.0)) THEN   											!ii,jj cell routes down
                III = DIREC(II,JJ,1)         										!to the subbasin outlet
                JJJ = DIREC(II,JJ,2)         										!point, following the
                II  = III                    										!direction of direc(,)
                JJ  = JJJ                    										!from each cell
                GOTO 500
            END IF                           											!if you get there,
        END IF                                										!no_of_box increments
510     CONTINUE 
        !WRITE(*,*) 'Res ',RES_DIRECT(RESORDER,1),'---> Res ',RES_DIRECT(RESORDER,2)
        WRITE(*,*) 'Number of grid cells upstream of present reservoir',no_of_box(RESORDER)
        RETURN
        END

C************************************************************************************************************************************************************************************
C     Calling Convolution and Considering Reservoir Rule Curves/Opearting Rules
C************************************************************************************************************************************************************************************

        SUBROUTINE MAKE_CONVOLUTIONRS
     & (RPATH, RESER, NCOL, NROW, NO_OF_BOX, PMAX, DAYS, CATCHIJ,
     &  BASE, RUNO, FLOW, KE, UH_DAY, FRACTION, FACTOR_SUM,
     &  XC, YC, SIZE, DPREC, INPATH,ICOL,NDAY,IDAY,IMONTH,IYEAR, START_YEAR,START_MO,
     &  MO, YR, NYR, VOL, FLOWIN, FLOWOUT,FLOWOUT_TURB, HHO, ENERGYPRO, HTK, DIREC, IROW,
     &  XNE, YNE, NORESERVOIRS, RES_DIRECT, RES_EVAPORATION, TRVLTIME,RESORDER,NAMERS5,NAME5,NSTATIONS_MAX,UH_SS,VRESER,NRESER_MAX, NSTATIONS,STEPBYSTEP,SIM_YEAR,NDAY_SIM,COUPLER_ITERATION,OUTPATH,UH_PATH,NF_PATH,NF_EXIST)
        IMPLICIT NONE
c       Declare variables
        INTEGER     N, DAYS, NDAY,NDAY_SIM, SIM_YEAR,START_YEAR,START_MO,NSTATIONS_MAX,NRESER_MAX, NSTATIONS,PREV_DAY,PREV_MON,PREV_YEAR,STEPS,UH_CLEN,NF_CLEN,OUT_CLEN
        INTEGER     NCOL,NROW,ICOL,PMAX,KE,UH_DAY,KK
        INTEGER     CATCHIJ(PMAX,2,NRESER_MAX,NSTATIONS_MAX)
        INTEGER     NYR,COUPLER_ITERATION
        INTEGER     RESER(NCOL,NROW)
        INTEGER     IDAY(DAYS), IMONTH(DAYS), IYEAR(DAYS)
        INTEGER     MO(12*NYR),YR(12*NYR)
        INTEGER     DIREC(NCOL,NROW,2)
        INTEGER     RES_DIRECT(NRESER_MAX,3,NSTATIONS_MAX)
        INTEGER     YEAR(NRESER_MAX,DAYS,NSTATIONS_MAX),MONTH(NRESER_MAX,DAYS,NSTATIONS_MAX),DAY(NRESER_MAX,DAYS,NSTATIONS_MAX)
        INTEGER     OPEYEAR(NRESER_MAX,NSTATIONS_MAX)
        INTEGER     IROW, CURRENTYEAR
        INTEGER     NO_OF_BOX(NRESER_MAX,NSTATIONS_MAX)
        INTEGER     XNE(NSTATIONS_MAX), YNE(NSTATIONS_MAX)
        INTEGER     NORESERVOIRS(NSTATIONS_MAX)
        INTEGER     I, J, K, L, CAL_MONTH, W, x, II,RESID,DS_RESID
        INTEGER     NO_OF_ROW(NRESER_MAX,NSTATIONS_MAX)
        INTEGER     RULE(NRESER_MAX,NSTATIONS_MAX), WITHIRRIGATION(NRESER_MAX,NSTATIONS_MAX),BTHMET_ACTIVE(NRESER_MAX,NSTATIONS_MAX),HYDROPOWER_GENERATOR(NRESER_MAX,NSTATIONS_MAX)
        INTEGER     STARTDAY(NRESER_MAX,NSTATIONS_MAX)
        INTEGER     DPREC, RESORDER(NSTATIONS_MAX,NRESER_MAX), NOOFROW(NSTATIONS_MAX), CRTDATE,RESORD(NSTATIONS_MAX)
        REAL        UH_SS(PMAX,KE+UH_DAY-1,NRESER_MAX),UH_RR(1,KE+UH_DAY-1)
        REAL        TRVLTIME(NRESER_MAX,NSTATIONS_MAX)
        REAL        BASE(DAYS,NSTATIONS_MAX), RUNO(DAYS,NSTATIONS_MAX)
        REAL        FLOW(DAYS,NSTATIONS_MAX)
        REAL        FRACTION(NCOL,NROW)
        REAL        PI, RERD, FACTOR, FACTOR_SUM
        REAL        CURRENTWL(NSTATIONS_MAX), DESIGNWL(NSTATIONS_MAX)
        REAL        RES_EVAPORATION(NRESER_MAX,DAYS)
        REAL        XC, YC, SIZE, TEMP
        REAL        VOL(NRESER_MAX,DAYS,NSTATIONS_MAX), FLOWIN(NRESER_MAX,DAYS,NSTATIONS_MAX), FLOWOUT(NRESER_MAX,DAYS,NSTATIONS_MAX),FLOWOUT_TURB(NRESER_MAX,DAYS,NSTATIONS_MAX),Qmin(NRESER_MAX,DAYS,NSTATIONS_MAX)
        REAL        HHO(NRESER_MAX,DAYS,NSTATIONS_MAX), ENERGYPRO(NRESER_MAX,DAYS,NSTATIONS_MAX), HTK(NRESER_MAX,DAYS,NSTATIONS_MAX)
        REAL        HMAX(NRESER_MAX,NSTATIONS_MAX), HMIN(NRESER_MAX,NSTATIONS_MAX)
        REAL        RESFLOWS(NRESER_MAX,DAYS,NSTATIONS_MAX) 
        REAL        HRESERMAX(NRESER_MAX,NSTATIONS_MAX), HRESERMIN(NRESER_MAX,NSTATIONS_MAX), VINITIAL(NRESER_MAX,NSTATIONS_MAX), H0(NRESER_MAX,NSTATIONS_MAX)
        REAL        QRESER(NRESER_MAX,NSTATIONS_MAX),VRESER(NRESER_MAX,NSTATIONS_MAX,DAYS),REALHEAD(NRESER_MAX,NSTATIONS_MAX),VRESERTHAT(NRESER_MAX,NSTATIONS_MAX),VDEAD(NRESER_MAX,NSTATIONS_MAX)
        REAL        QRESERTHAT(NRESER_MAX,NSTATIONS_MAX), HYDRAUHEAD(NRESER_MAX,NSTATIONS_MAX)
        REAL        OP1(NRESER_MAX,2,NSTATIONS_MAX) ! rule curve
        REAL        OP4(NRESER_MAX,DAYS,NSTATIONS_MAX),OP6(NRESER_MAX,DAYS,NSTATIONS_MAX) ! pre-defined time series data
        REAL        X1(NRESER_MAX,NSTATIONS_MAX), X2(NRESER_MAX,NSTATIONS_MAX), X3(NRESER_MAX,NSTATIONS_MAX), X4(NRESER_MAX,NSTATIONS_MAX)	! operating rule
        REAL        BTHMET_A(NRESER_MAX,NSTATIONS_MAX),BTHMET_B(NRESER_MAX,NSTATIONS_MAX),BTHMET_C(NRESER_MAX,NSTATIONS_MAX),BTHMET_D(NRESER_MAX,NSTATIONS_MAX)
        REAL        SEEPAGE(NRESER_MAX,NSTATIONS_MAX), INFILTRATION(NRESER_MAX,NSTATIONS_MAX), Demand(NRESER_MAX,NSTATIONS_MAX)
        REAL        TEMPO
        REAL        RESLV(NRESER_MAX,12,NSTATIONS_MAX), OP5X1(NRESER_MAX,12,NSTATIONS_MAX), OP5X2(NRESER_MAX,12,NSTATIONS_MAX), OP5X3(NRESER_MAX,12,NSTATIONS_MAX), OP5X4(NRESER_MAX,12,NSTATIONS_MAX), DEMAND5(NRESER_MAX,12,NSTATIONS_MAX)
        REAL        IRRIGATION(NRESER_MAX,DAYS),ENV_FLOW(NRESER_MAX,NSTATIONS_MAX)
        REAL        KFACTOR(NSTATIONS_MAX)
        CHARACTER*10 WTEMP, W_STR,NDAY_SIM_STR,NDAY_SIM_STR_NEXT,NDAY_SIM_STR_PREV
        CHARACTER*72 RESFLOWSTRING,RESFLOWSTRING_SIM,RESID_STR,DS_RESID_STR
        CHARACTER*72 RPATH
        CHARACTER*72 IPATH(NSTATIONS_MAX)
        CHARACTER*100 TEMPRPATH(NSTATIONS_MAX)
        CHARACTER*100 PATHRES3,PATHVOL3
        CHARACTER*200 INPATH,OUTPATH,UH_PATH,NF_PATH,UH_NAME
        CHARACTER*20 CHUOI(NSTATIONS_MAX)
        CHARACTER*21  NAMERS5(NSTATIONS_MAX,NRESER_MAX)
        CHARACTER*10  NAME5(NSTATIONS_MAX)
        CHARACTER*100 FILENAME
        CHARACTER*100 storage_filename,previous_storage_filename,previous_natflow_filename
        CHARACTER*100 res_filename
        LOGICAL      STEPBYSTEP
        CHARACTER(LEN=300) :: check_file
        LOGICAL :: file_exists,NF_EXIST,RES_IN_BASIN
        
c       Subroutine main body    

c------------------------------------------------------------------------------------------------------------------------------------------------------------------            
c     STEP 1: Reading unit hydrograph files (UH_S) for reservoirs (NORESERVOIRS) and station (W)
c------------------------------------------------------------------------------------------------------------------------------------------------------------------            
c       Look for reservoirs sequences
        IF (STEPBYSTEP) THEN
            CURRENTYEAR = SIM_YEAR
        ELSE
            CURRENTYEAR = START_YEAR   ! This needs to be updated to account for year of simulation instead of start_year (otherwise we need to make sure that all dams are set to be operating before start_year)
        ENDIF
        
        DO W=1,NSTATIONS
            print*,'Number of Stations=', NSTATIONS
            print*,'Station:',W
            print*, 'Number of Reservoirs', NORESERVOIRS(W)
c------------------------------------------------------------------------------------------------------------------------------------------------------------------            
c     STEP 2: Making Convolution to Produce Naturalized Flow 
c------------------------------------------------------------------------------------------------------------------------------------------------------------------                        
            OUT_CLEN=INDEX(OUTPATH,' ')-1    ! Character Length for Output Path (from configuration.txt)  
            NF_CLEN=INDEX(NF_PATH,' ')-1 
            UH_CLEN=INDEX(UH_PATH,' ')-1    ! Character Length for UH_Path
            WRITE(WTEMP, 10) W
            write(NDAY_SIM_STR,10) NDAY_SIM    
10          FORMAT(I4)
            RESFLOWSTRING = trim(trim('NF'//trim(ADJUSTL(WTEMP)))//'.txt')
            print*,'Checking whether to make convolution or not.....'
            IF (STEPBYSTEP .AND. ((NDAY_SIM.GT.1) .OR. (COUPLER_ITERATION.GT.1))) THEN
                print*, 'No need to rerun MAKE_CONVOLUTION for Step-By-Step Version'  ! Note that natuarlized flow from previous run has to be prepared for the same number of days (NDAY), i.e., same range of start/end year/month of simulation 
                print*, 'Read updated naturalized flow from previous run'
                !READ RESFLOWS
                RESFLOWSTRING_SIM = trim(trim('NF'//trim(ADJUSTL(WTEMP)))//'_STEP'//trim(adjustl(NDAY_SIM_STR))//'.txt')
                open(72, FILE=NF_PATH(1:NF_CLEN)//'/stepbystep/'//RESFLOWSTRING_SIM, status='old')
                DO N = 1, NORESERVOIRS(W)
                    READ(72, *) (RESFLOWS(N,K,W), K = 1, NDAY)
                ENDDO
                !close(72,status='delete')
                close(72)
            ELSEIF (STEPBYSTEP .AND. (NDAY_SIM.EQ.1) .AND. (NF_EXIST)) THEN
                print*, 'No need to run MAKE_CONVOLUTION for Step-By-Step Version'
                print*, 'Read updated naturalized flow from previous run'
                !READ RESFLOWS
                RESFLOWSTRING_SIM = trim(trim('NF'//trim(ADJUSTL(WTEMP)))//'_STEP'//trim(adjustl(NDAY_SIM_STR))//'.txt')
                open(72, FILE=NF_PATH(1:NF_CLEN)//'/stepbystep/'//RESFLOWSTRING_SIM, status='old')
                DO N = 1, NORESERVOIRS(W)
                    READ(72, *) (RESFLOWS(N,K,W), K = 1, NDAY)
                ENDDO
                
            ELSEIF (.NOT. STEPBYSTEP .AND. (NF_EXIST)) THEN ! Note that natuarlized flow from previous run has to be prepared for the same number of days (NDAY), i.e., same range of start/end year/month of simulation 
                print*, 'No need to run MAKE_CONVOLUTION'
                print*, 'Read updated naturalized flow from previous run'
                !READ RESFLOWS
                RESFLOWSTRING = trim(trim('NF'//trim(ADJUSTL(WTEMP)))//'.txt')
                open(72, FILE = NF_PATH(1:NF_CLEN)//RESFLOWSTRING, status='old')
                DO N = 1, NORESERVOIRS(W)
                    READ(72, *) (RESFLOWS(N,K,W), K = 1, NDAY)
                ENDDO
       
                
            ELSE ! Not STEPBYSTEP Version or STEPBYSTEP with step=1 & NF_EXIST=FALSE
                print*,'reading grid from files uh_s' 
                ! Reading UH_S files for reservoirs
                IF (STEPBYSTEP) THEN
                    DO J = 1,(NORESERVOIRS(W)-1)
                        OPEN(98,file = UH_PATH(1:UH_CLEN)//trim(adjustl(NAMERS5(W,RESORDER(W,J))))//'.uh_s', status='old')
                        DO N = 1,NO_OF_BOX(RESORDER(W,J),W)
                            READ(98, *) (UH_SS(N,K,RESORDER(W,J)), K = 1,KE+UH_DAY-1)  
                        ENDDO
                        CLOSE(98)
                    ENDDO
                    ! Reading UH_S files for last reservoir, which is the station under consideration
                    J=NORESERVOIRS(W)
                    OPEN(97,file = UH_PATH(1:UH_CLEN)//trim(adjustl(NAME5(W)))//'.uh_s', status='old')
                    DO N= 1,NO_OF_BOX(RESORDER(W,J),W)
                        READ(97, *) (UH_SS(N,K,RESORDER(W,J)), K = 1,KE+UH_DAY-1)
                    ENDDO
                    CLOSE(97)
                ELSE ! Not STEPBYSTEP Version
                    DO J = 1,(NORESERVOIRS(W)-1)
                        OPEN(98,file = UH_PATH(1:UH_CLEN)//trim(adjustl(NAMERS5(W,RESORDER(W,J))))//'.uh_s', status='old')
                        DO N = 1,NO_OF_BOX(RESORDER(W,J),W)
                            READ(98, *) (UH_SS(N,K,RESORDER(W,J)), K = 1,KE+UH_DAY-1)  
                        ENDDO
                        CLOSE(98)
                    ENDDO    
                    ! Reading UH_S files for station
                    J=NORESERVOIRS(W)
                    OPEN(97,file = UH_PATH(1:UH_CLEN)//trim(adjustl(NAME5(W)))//'.uh_s', status='old')
                    DO N= 1,NO_OF_BOX(RESORDER(W,J),W)
                        READ(97, *) (UH_SS(N,K,RESORDER(W,J)), K = 1,KE+UH_DAY-1)
                    ENDDO
                    CLOSE(97)
                ENDIF
                print*, 'Making Convolution to produce naturalized flow ...'
                RESFLOWS(:,:,W)=0
                DO N = 1,NORESERVOIRS(W)
                    print*,N, 'Working on reservoir no...',RES_DIRECT(N,1,W)
                    !RESFLOWS(N,:,W)=0
                    CALL MAKE_CONVOLUTION
     &                   (RPATH,RESER,NCOL, NROW, NO_OF_BOX(:,W), PMAX, DAYS,
     &                   CATCHIJ(:,:,:,W), BASE(:,W), RUNO(:,W), FLOW(:,W), KE, UH_SS,UH_DAY,FRACTION,
     &                   FACTOR_SUM,XC,YC,SIZE,DPREC,INPATH,ICOL,NDAY,
     &                   IDAY,IMONTH,IYEAR, START_YEAR, START_MO, MO, YR, NYR, VOL(:,:,W),
     &                   FLOWIN(:,:,W), FLOWOUT(:,:,W), HHO(:,:,W), RESFLOWS(:,:,W),N,RES_DIRECT(N,1,W),
     &                   RES_EVAPORATION, NORESERVOIRS(W),NRESER_MAX)    
                ENDDO
c               !SAVE RESFLOWS TO READ THEM LATER    
                IF (STEPBYSTEP) THEN
                    RESFLOWSTRING= trim(trim('NF'//trim(ADJUSTL(WTEMP)))//'_STEP'//trim(adjustl(NDAY_SIM_STR))//'.txt')
                    open(71, FILE = NF_PATH(1:NF_CLEN)//'/stepbystep/'//RESFLOWSTRING, status='unknown')
                ELSE
                    RESFLOWSTRING = trim(trim('NF'//trim(ADJUSTL(WTEMP)))//'.txt')
                    open(71, FILE = NF_PATH(1:NF_CLEN)//RESFLOWSTRING, status='unknown')
                ENDIF
                DO N = 1,NORESERVOIRS(W)
                    WRITE(71, *) (RESFLOWS(N,K,W), K = 1,NDAY)
                ENDDO
                close(71)
c               !READ RESFLOWS
                IF (STEPBYSTEP) THEN
                    open(72, FILE=NF_PATH(1:NF_CLEN)//'/stepbystep/'//RESFLOWSTRING, status='old',ERR=9004)
                ELSE
                    open(72, FILE=NF_PATH(1:NF_CLEN)//RESFLOWSTRING, status='old',ERR=9004)
                ENDIF
                DO N = 1, NORESERVOIRS(W)
                    READ(72, *) (RESFLOWS(N,K,W), K = 1, NDAY)
                ENDDO
                close(72)             
            ENDIF
            
c------------------------------------------------------------------------------------------------------------------------------------------------------------------            
c     STEP 3: Initiate Reservoir Variables (Fluxes and States) and Read Reservoir Parameters
c------------------------------------------------------------------------------------------------------------------------------------------------------------------                                    
            IF (STEPBYSTEP) THEN
               ! I = 1,2            ! only one step to run and next step to save I+1 step variables
                DO I = 1, 2
                    DO J = 1, NORESERVOIRS(W)-1
                        FLOWIN(J,I,W) = 0
                        FLOWOUT(J,I,W) = 0
                        Qmin(J,I,W) = 0
                        VOL(J,I,W) = 0
                        ENERGYPRO(J,I,W) = 0
                        HTK(J,I,W) = 0
                    END DO
                END DO
                IF ((NDAY_SIM.GT.1) .OR. (COUPLER_ITERATION.GT.1)) THEN
                    write(W_STR,10) W ! convert integer to string to be able to concatenate in filename with strings
                    write(NDAY_SIM_STR,10) NDAY_SIM
                    storage_filename='RESERVOIR_VOL_STATION'//trim(adjustl(W_STR))//'_STEP'//trim(adjustl(NDAY_SIM_STR))//'.txt'
                    open(99,FILE = OUTPATH(1:OUT_CLEN)//'/stepbystep/'//storage_filename,STATUS='OLD',ERR=9004)
                    READ(99, *) (VOL(J,1,W), J = 1, NORESERVOIRS(W))
                    close(99)
                END IF  
            ELSE ! Not STEPBYSTEP Version 
                DO I = 1, NDAY
                    DO J = 1, NORESERVOIRS(W)-1
                        FLOWIN(J,I,W) = 0
                        FLOWOUT(J,I,W) = 0
                        Qmin(J,I,W) = 0
                        VOL(J,I,W) = 0
                        ENERGYPRO(J,I,W) = 0
                        HTK(J,I,W) = 0
                    END DO
                END DO
            ENDIF
c------------------------------------------------------------------------------------------------------------------------------------------------------------------            
c     STEP 4: Read Reservoir Parameters
c------------------------------------------------------------------------------------------------------------------------------------------------------------------                         
!           Read reservoir characteristics/parameters from input file (/Reservoirs/res1.txt, /Reservoirs/res2.txt,...)
            print*, 'Reading Reservoir Parameters.....'
            DO  J  = 1, NORESERVOIRS(W)-1
                WRITE(CHUOI(W),*) RES_DIRECT(J,1,W)
                TEMPRPATH(W) = trim(RPATH)//"res"//trim(ADJUSTL(CHUOI(W)))//".txt"
                OPEN(25, FILE = TEMPRPATH(W),FORM = 'FORMATTED',
     &          STATUS='OLD',ERR=9002)
                READ(25,*)  
                READ(25,*) HRESERMAX(J,W),HRESERMIN(J,W),VRESERTHAT(J,W),VDEAD(J,W),
     &          HYDRAUHEAD(J,W), QRESERTHAT(J,W), OPEYEAR(J,W), VINITIAL(J,W)
                READ(25,*)
                READ(25,*) SEEPAGE(J,W), INFILTRATION(J,W)
                READ(25,*)
                READ(25,*) WITHIRRIGATION(J,W)
                READ(25,'(A)') IPATH(W)
                READ(25,*)
                ! convert units of storage to m3
                !VRESERTHAT(J,W)=VRESERTHAT(J,W)*1000
                !VDEAD(J,W)=VDEAD(J,W)*1000
                !VINITIAL(J,W)=VINITIAL(J,W)*1000
                IF (WITHIRRIGATION(J,W) .EQ. 0) THEN
                    DO L = 1, DAYS-1
                        IRRIGATION(J,L) = 0
                    END DO
                ELSE IF (WITHIRRIGATION(J,W) .EQ. 1) THEN
                    OPEN(27, FILE = IPATH(W),FORM = 'FORMATTED',STATUS='OLD',ERR=9003)
                    READ(27,*)
                    READ(27,*) NOOFROW(W)
                    READ(27,*)
                    DO L = 1, NOOFROW(W)
                        READ(27,*) IRRIGATION(J,L)
                    END DO
                    CLOSE(27)
                END IF
                
                READ(25,*) HYDROPOWER_GENERATOR(J,W)  ! Read the HYDROPOWER_GENERATOR of reservoir operation; constrained by hydropower or not (HYDROPOWER_GENERATOR=1 -> Generators ON / HYDROPOWER_GENERATOR=0 -> Disconnect turbines from generators)
                READ(25,*)
                READ(25,*) ENV_FLOW(J,W)  ! Read the percentage of environmental flow (% of maximum flow)
                READ(25,*)
                READ(25,*) RULE(J,W)
                IF (RULE(J,W) .EQ. 1) THEN ! simplified rule curve
                    READ(25,*) HMAX(J,W), HMIN(J,W), OP1(J,1,W),OP1(J,2,W)
                    IF (HMAX(J,W) .LT. HMIN(J,W)) THEN
                        TEMP=HMAX(J,W)
                        HMAX(J,W)=HMIN(J,W)
                        HMIN(J,W)=TEMP
                    ENDIF
                ELSE IF (RULE(J,W) .EQ. 2) THEN ! rule curve
                    READ(25,*) RESLV(J,1,W), RESLV(J,2,W), RESLV(J,3,W), RESLV(J,4,W), RESLV(J,5,W), RESLV(J,6,W),
     &              RESLV(J,7,W), RESLV(J,8,W), RESLV(J,9,W), RESLV(J,10,W), RESLV(J,11,W), RESLV(J,12,W)
                ELSE IF (RULE(J,W) .EQ. 3) THEN ! operating rule
                    READ(25,*) Demand(J,W), X1(J,W), X2(J,W), X3(J,W), X4(J,W)
                ELSE IF (RULE(J,W) .EQ. 4) THEN! predefined time-series of release
                    READ(25,'(A)')FILENAME
                    PATHRES3 = trim(trim(ADJUSTL(FILENAME)))
		            print*,'Predefined Time Series of Release:', PATHRES3
                    OPEN(26, FILE = PATHRES3,FORM = 'FORMATTED',ERR=9002)
                    READ(26,*) NO_OF_ROW(J,W) 						! Read the number of rows of reservoir J th
                    DO L = 1, NO_OF_ROW(J,W)
                        READ(26,*) YEAR(J,L,W), MONTH(J,L,W), DAY(J,L,W), OP4(J,L,W)
                        IF ((YEAR(J,L,W) .EQ. START_YEAR) .AND. (MONTH(J,L,W) .EQ. START_MO) .AND. (DAY(J,L,W) .EQ. 1)) THEN
                            STARTDAY(J,W) = L-1							! If the time-series is shorter than the simulation period (release = 0); Changed to "L-1" instead of "L" so STARDAY=0 if time series starts with simulation period
                        END IF
                    END DO
                    CLOSE(26)
                ELSE IF (RULE(J,W) .EQ. 5) THEN
                    DO L = 1, 12
                        READ(25,*) DEMAND5(J,L,W), OP5X1(J,L,W), OP5X2(J,L,W), OP5X3(J,L,W), OP5X4(J,L,W)
                    END DO
                ELSE IF (RULE(J,W) .EQ. 6) THEN  !predefined time-series of storage
                    READ(25,'(A)')FILENAME
                    PATHRES3 = trim(trim(ADJUSTL(FILENAME)))
		            print*,'Predefined Time Series of Storage:', PATHRES3
                    OPEN(27, FILE = PATHRES3,FORM = 'FORMATTED',ERR=9002)
                    READ(27,*) NO_OF_ROW(J,W) 						! Read the number of rows of reservoir J th
                    DO L = 1, NO_OF_ROW(J,W)
                        READ(27,*) YEAR(J,L,W), MONTH(J,L,W), DAY(J,L,W), OP6(J,L,W)
                        IF ((YEAR(J,L,W) .EQ. START_YEAR) .AND. (MONTH(J,L,W) .EQ. START_MO) .AND. (DAY(J,L,W) .EQ. 1)) THEN
                            STARTDAY(J,W) = L-1							! If the time-series is shorter than the simulation period (storage = 0); Changed to "L-1" instead of "L" so STARDAY=0 if time series starts with simulation period
                        END IF
                    END DO
                    CLOSE(27)
                END IF
                ! Read Bathymetry Parameters
                READ(25,*)
                READ(25,*) BTHMET_ACTIVE(J,W)
                !######################## Calculate the bottom level of reservoir [H0] [Volume =Zero]
                IF (BTHMET_ACTIVE(J,W).EQ.1) THEN
                    READ(25,*) BTHMET_A(J,W), BTHMET_B(J,W),BTHMET_C(J,W),BTHMET_D(J,W)
                    H0(J,W)=BTHMET_A(J,W)/(1+BTHMET_D(J,W))   ! Volume =Zero in the Fitting Equation
                ELSE
                    KFACTOR(W) = SQRT(VDEAD(J,W)/VRESERTHAT(J,W))							! calculate the bottom level of reservoir
                    H0(J,W) = (KFACTOR(W)*HRESERMAX(J,W) - HRESERMIN(J,W))/(KFACTOR(W)-1)		! pyramid shape reservoir
                    !H0(J,W) = (HRESERMAX(J,W)-HRESERMIN(J,W)*(VRESERTHAT(J,W)+VDEAD(J,W))/VDEAD(J,W))/(1-(VRESERTHAT(J,W)+VDEAD(J,W))/VDEAD(J,W))    !This actually performs better, but we will fix the problem of H0 with reservoirs bathimetry
                ENDIF
                CLOSE(25)
                 
                IF (CURRENTYEAR>=OPEYEAR(J,W)) THEN
                    IF (STEPBYSTEP .AND. ((NDAY_SIM.GT.1) .OR. (COUPLER_ITERATION.GT.1))) THEN   ! CURRENTYEAR=START_YEAR; do we need to change this condition for the case when NDAY_SIM>1 but CURRENTYEAR>=OPEYEAR, then we need to assign VOL=Vinitial?
                        print*,'Read the storage volume from the previous time step'
                    ELSE
                        VOL(J,1,W) = VINITIAL(J,W)
                    ENDIF 
                ELSE 
                    VOL(J,1,W) = 0.001
                END IF
            END DO   ! End of J Loop for NORESERVOIRS [Line 916]
        END DO        ! End of W Loop for NSTATIONS [Line 790]
c------------------------------------------------------------------------------------------------------------------------------------------------------------------            
c     STEP 5: Reservoir Modeling Considering Operating Curve
c------------------------------------------------------------------------------------------------------------------------------------------------------------------             
        print*, 'Reservoir Modeling Considering Operating Curve.....'
        IF (STEPBYSTEP) THEN
           STEPS = 1            ! only one step to run
        ELSE
           STEPS=NDAY
        ENDIF
        !LOOP OVER TIMESTEPS
        DO I = 1, STEPS
            IF (STEPBYSTEP) THEN
                II=NDAY_SIM    ! II=index of simulated day 
            ELSE
                II=I    ! II=same as the steps (NDAY)
            ENDIF   
            !LOOP OVER STATIONS AND  RESERVOIRS FOR OPERATING CURVE  
            DO W=1, NSTATIONS
                DO J = 1, NORESERVOIRS(W)-1
c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx            
c     STEP 5-1: Check if the reservoir has started operation (reaching the year of construction) 
c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                    CURRENTYEAR = START_YEAR + INT(II/365)		! approximate, does not consider leap years
                    IF (CURRENTYEAR>=OPEYEAR(J,W)) THEN
                        VRESER(J,W,I) = (VRESERTHAT(J,W)+VDEAD(J,W))  
                        REALHEAD(J,W) = HYDRAUHEAD(J,W)
                        QRESER(J,W) = QRESERTHAT(J,W)    
                    ELSE
                        GOTO 125
                    END IF
                    ! Consider Inflow = Naturalized Flow
                    FLOWIN(J,I,W) = RESFLOWS(J,II,W)     ! change here I (step) to II (simulated day)
                    ! Considering Environmental Flow (Minimum Flow = 10% of mean flow or ~ 5% of the Qmax=QRESER)
                    Qmin(J,I,W)=ENV_FLOW(J,W) * QRESERTHAT(J,W)
c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx            
c     STEP 5-2: Reservoir Operation using Strategy 1 or 2 
c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                    !Calculate the designed water level
                    !Note RULE = 1: simplified rule curve - 2: rule curve - 3: operating rules: - 4 pre-defined time-series data - 5 12 month operating rule
                    CRTDATE = INT(1.0* mod(II,365)+(START_MO-1)*30)						! approximate
                    IF ((RULE(J,W) .EQ. 1) .or. (RULE(J,W) .EQ. 2)) THEN   ! (Options 1 and 2: rule curves)
                        IF (CURRENTYEAR<OPEYEAR(J,W)) THEN
                            FLOWOUT(J,I,W) = FLOWIN(J,I,W)
                            VOL(J,I+1,W) = VOL(J,I,W)
                            GOTO 123
                        END IF
                        IF (RULE(J,W) .EQ. 1) THEN
                            IF (OP1(J,1,W)>OP1(J,2,W)) THEN
                            ! Caculate target water level
                                IF ((CRTDATE .GT. OP1(J,2,W)) .and. (CRTDATE .LT. OP1(J,1,W))) THEN     ! OP1 are the time T1 and T2 at which predefined maximum (HMAX) and minimum (HMIN) levels are reached
                                    DESIGNWL(W)=(CRTDATE-OP1(J,2,W))/(OP1(J,1,W)-OP1(J,2,W))*(HMAX(J,W)-HMIN(J,W))
                                ELSE IF (CRTDATE .GE. OP1(J,1,W)) THEN 
                                    DESIGNWL(W)=(HMAX(J,W)-HMIN(J,W))
     &                                  -(CRTDATE-OP1(J,1,W))/(365-OP1(J,1,W)+OP1(J,2,W))*
     &                                  (HMAX(J,W)-HMIN(J,W))
                                ELSE
                                    DESIGNWL(W)=(HMAX(J,W)-HMIN(J,W))
     &                              -(CRTDATE+365-OP1(J,1,W))/(365-OP1(J,1,W)+OP1(J,2,W))*
     &                              (HMAX(J,W)-HMIN(J,W))
                                ENDIF
                            ELSE
                                IF ((CRTDATE .GT. OP1(J,1,W)) .and. (CRTDATE .LT. OP1(J,2,W))) THEN
                                    DESIGNWL(W)=(HMAX(J,W)-HMIN(J,W))-(CRTDATE-OP1(J,1,W))/(OP1(J,2,W)-OP1(J,1,W))*(HMAX(J,W)-HMIN(J,W))
                                ELSE IF (CRTDATE .GE. OP1(J,2,W)) THEN 
                                    DESIGNWL(W)=(HMAX(J,W)-HMIN(J,W))/(365-OP1(J,2,W)+OP1(J,1,W))*(CRTDATE-OP1(J,2,W))
                                ELSE 
                                    DESIGNWL(W)=(HMAX(J,W)-HMIN(J,W))/(365-OP1(J,2,W)+OP1(J,1,W))*
     &                              (CRTDATE-OP1(J,1,W))+(HMAX(J,W)-HMIN(J,W))
                                ENDIF
                            END IF
                            DESIGNWL(W) = DESIGNWL(W) + (HMIN(J,W) - H0(J,W))
                        ELSE   ! Operation Rule #2
                            DESIGNWL(W) = RESLV(J,CAL_MONTH(CRTDATE),W)
                            DESIGNWL(W) = DESIGNWL(W) - H0(J,W)
                        ENDIF
                        HTK(J,I,W) = DESIGNWL(W) + H0(J,W)													! water head
                        !Calculate the Water Level using Pyramid Shape of bathymetry or provided Parameters
                        IF (BTHMET_ACTIVE(J,W).EQ.1) THEN
                            CURRENTWL(W)=BTHMET_A(J,W)/(EXP(-1*VOL(J,I,W)/1000*BTHMET_B(J,W))+VOL(J,I,W)/1000*BTHMET_C(J,W)+BTHMET_D(J,W))   ! Fitting Equation by Shanti [convert volume from 1000M3 to MCM]
                        ELSE
                            CURRENTWL(W) = VOL(J,I,W) * (HRESERMAX(J,W)-H0(J,W))/VRESER(J,W,I)
                        ENDIF
                        ! Check Zones/Cases of Rule Curve 
                        IF (CURRENTWL(W)>=DESIGNWL(W)) THEN											! Zone 3
                            IF ((VOL(J,I,W)+FLOWIN(J,I,W)*24*3.6 -QRESER(J,W)*24*3.6)						! Case 2
     &                              >(DESIGNWL(W)* VRESER(J,W,I)) /(HRESERMAX(J,W)-H0(J,W))) THEN
                                VOL(J,I+1,W) = VOL(J,I,W) + FLOWIN(J,I,W)*24*3.6-QRESER(J,W)*24*3.6
                                FLOWOUT(J,I,W) = QRESER(J,W)
                                !VOL(J,I+1,W)=(DESIGNWL(W)*VRESER(J,W,I))/(HRESERMAX(J,W)-H0(J,W))  ! from Shanti Version (no need for this line: bug fixed by Bruno)
                            ELSE																	! Case 1
                                VOL(J,I+1,W)=(DESIGNWL(W)*VRESER(J,W,I))/(HRESERMAX(J,W)-H0(J,W))
                                FLOWOUT(J,I,W)=(VOL(J,I,W)-VOL(J,I+1,W))/24/3.6+ FLOWIN(J,I,W)
                            END IF
                            GOTO 123   ! Done with the Operating Curve
                        ELSE																		! Zone 2
                            IF ((VOL(J,I,W)+FLOWIN(J,I,W)*24*3.6)>((DESIGNWL(W)* VRESER(J,W,I))				! Case 2
     &                              /(HRESERMAX(J,W)-H0(J,W)))) THEN
                                VOL(J,I+1,W)=(DESIGNWL(W) * VRESER(J,W,I))/(HRESERMAX(J,W)-H0(J,W))
                                FLOWOUT(J,I,W)=FLOWIN(J,I,W)-(VOL(J,I+1,W)-VOL(J,I,W))/24/3.6
                                IF (FLOWOUT(J,I,W)>QRESER(J,W)) THEN
                                    VOL(J,I+1,W)=VOL(J,I+1,W)+(FLOWOUT(J,I,W)-QRESER(J,W))*24*3.6
                                    FLOWOUT(J,I,W) = QRESER(J,W)
                                END IF
                                GOTO 123
                            ELSE																	! Case 1 (this case covers Zone 1 + Zone 2 case 1)
                                VOL(J,I+1,W) = VOL(J,I,W) + FLOWIN(J,I,W)*24*3.6
                                FLOWOUT(J,I,W) = Qmin(J,I,W)
                                GOTO 123
                            END IF
                        END IF
                        
c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx            
c     STEP 5-3: Reservoir Operation using Strategy 3 
c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx                    
                    ELSE IF (RULE(J,W) .EQ. 3) THEN ! opearting rule (Option 3)
                        ! Note x1 and x4 in radian (0 to pi/2), not degree
                        IF (VOL(J,I,W) < VDEAD(J,W)) THEN ! below dead water level
                            FLOWOUT(J,I,W) = 0																					! case 1
                        ELSE IF (VOL(J,I,W) .LE. X2(J,W)) THEN ! hedging
                            IF ((VOL(J,I,W)-VDEAD(J,W)+FLOWIN(J,I,W)*24*3.6) .LE. (Demand(J,W)+(VOL(J,I,W)-X2(J,W))*tan(X1(J,W))*24*3.6)) THEN		! discharge more than the water amount in reservoir (case 2)
                                FLOWOUT(J,I,W) =(VOL(J,I,W)-VDEAD(J,W))/24/3.6+FLOWIN(J,I,W)											! discharge all of water in the reservoir
                            ELSE	! (case 3)
                                FLOWOUT(J,I,W) = Demand(J,W)+(VOL(J,I,W)-X2(J,W))*tan(X1(J,W))
                            END IF
                        ELSE IF (VOL(J,I,W).GE. X3(J,W)) THEN 	! spilling
                            IF ((Demand(J,W) + (VOL(J,I,W)-X3(J,W))*tan(X4(J,W)))*24*3.6>(VOL(J,I,W))-VDEAD(J,W)) THEN		! discharge more than the water in reservoir (just to make sure)
                                FLOWOUT(J,I,W) = (VOL(J,I,W)-VDEAD(J,W))/24/3.6										! discharge all of water in the reservoir
                            ELSE IF ((Demand(J,W) + (VOL(J,I,W)-X3(J,W))*tan(X4(J,W)))>QRESER(J,W)) THEN
                                FLOWOUT(J,I,W) = QRESER(J,W)
                            ELSE	 ! case 5
                                FLOWOUT(J,I,W) = Demand(J,W) + (VOL(J,I,W)-X3(J,W))*tan(X4(J,W))
                            END IF
                        ELSE ! releasing
                            IF (Demand(J,W)*24*3.6>(VOL(J,I,W)-VDEAD(J,W))) THEN
                                FLOWOUT(J,I,W) = (VOL(J,I,W)-VDEAD(J,W))/24/3.6
                            ELSE
                                FLOWOUT(J,I,W) = Demand(J,W)			! case 4
                            END IF
                        END IF
                        IF (FLOWOUT(J,I,W)>QRESER(J,W)) THEN			! double check just in case users chose unrealistic x1 and x4
                            FLOWOUT(J,I,W) = QRESER(J,W)
                        END IF
                        VOL(J,I+1,W) = VOL(J,I,W) + (FLOWIN(J,I,W)-FLOWOUT(J,I,W))*24*3.6
                        GOTO 123
c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx            
c     STEP 5-4: Reservoir Operation using Strategy 4
c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                    ELSE IF (RULE(J,W) .EQ. 4) THEN! pre-defined time series (Option 4)
                        IF ((OP4(J,I+STARTDAY(J,W),W) .GT. QRESER(J,W)) .AND. (VOL(J,I,W) .LT. VRESER(J,W,I))) THEN
                            OP4(J,I+STARTDAY(J,W),W) = QRESER(J,W)
                        END IF
                        IF (OP4(J,I+STARTDAY(J,W),W)*24*3.6 .GT.((VOL(J,I,W)-VDEAD(J,W)))+FLOWIN(J,I,W)*24*3.6) THEN
                            FLOWOUT(J,I,W) =(VOL(J,I,W)-VDEAD(J,W))/24/3.6 + FLOWIN(J,I,W)    ! Maximum flow that can be released (inflow + storage)
                        ELSE
                            FLOWOUT(J,I,W) = OP4(J,I+STARTDAY(J,W),W)
                        END IF
                        IF (FLOWOUT(J,I,W)<0) THEN
                            FLOWOUT(J,I,W)=0
                        END IF
                        IF (HYDROPOWER_GENERATOR(J,W) .EQ. 1) THEN
                            FLOWOUT_TURB(J,I,W)=FLOWOUT(J,I,W)   ! forcing the turbine release to be equal to the provided release in operation strategy 4 (otherwise flowout is going through changes due to changes in storage; if storage is greater than scap)
                        ELSE
                            FLOWOUT_TURB(J,I,W)=0.0
                        END IF    
                        VOL(J,I+1,W) = VOL(J,I,W) + (FLOWIN(J,I,W)-FLOWOUT(J,I,W))*24*3.6
                        GOTO 123
c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx            
c     STEP 5-5: Reservoir Operation using Strategy 5
c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                    ELSE IF (RULE(J,W) .EQ. 5) THEN ! note: this option is similar to OP3 but for a periodic demand
                    ! Note x1 and x4 in radian (0 to pi/2), not degree, this part can be shorthen
                        X1(J,W) = OP5X1(J,CAL_MONTH(CRTDATE),W)
                        X2(J,W) = OP5X2(J,CAL_MONTH(CRTDATE),W)
                        X3(J,W) = OP5X3(J,CAL_MONTH(CRTDATE),W)
                        X4(J,W) = OP5X4(J,CAL_MONTH(CRTDATE),W)
                        Demand(J,W) = DEMAND5(J,CAL_MONTH(CRTDATE),W)
                        IF (VOL(J,I,W) < VDEAD(J,W)) THEN ! below dead water level
                            FLOWOUT(J,I,W) = 0																					! case 1
                        ELSE IF (VOL(J,I,W) .LE. X2(J,W)) THEN ! hedging
                            IF ((VOL(J,I,W)-VDEAD(J,W)+FLOWIN(J,I,W)*24*3.6) .LE. (Demand(J,W)+(VOL(J,I,W)-X2(J,W))*tan(X1(J,W))*24*3.6)) THEN		! discharge more than the water amount in reservoir (case 2)                    
                                FLOWOUT(J,I,W) =(VOL(J,I,W)-VDEAD(J,W))/24/3.6+FLOWIN(J,I,W)											! discharge all of water in the reservoir
                            ELSE	! (case 3)
                                FLOWOUT(J,I,W) = Demand(J,W)+(VOL(J,I,W)-X2(J,W))*tan(X1(J,W))
                            END IF
                        ELSE IF (VOL(J,I,W).GE. X3(J,W)) THEN 	! spilling
                            IF ((Demand(J,W) + (VOL(J,I,W)-X3(J,W))*tan(X4(J,W)))*24*3.6>(VOL(J,I,W))-VDEAD(J,W)) THEN		! discharge more than the water in reservoir (just to make sure)
                                FLOWOUT(J,I,W) = (VOL(J,I,W)-VDEAD(J,W))/24/3.6										! discharge all of water in the reservoir
                            ELSE IF ((Demand(J,W) + (VOL(J,I,W)-X3(J,W))*tan(X4(J,W)))>QRESER(J,W)) THEN
                                FLOWOUT(J,I,W) = QRESER(J,W)
                            ELSE	 ! case 5
                                FLOWOUT(J,I,W) = Demand(J,W) + (VOL(J,I,W)-X3(J,W))*tan(X4(J,W))
                            END IF
                        ELSE ! releasing
                            IF (Demand(J,W)*24*3.6>(VOL(J,I,W)-VDEAD(J,W))) THEN
                                FLOWOUT(J,I,W) = (VOL(J,I,W)-VDEAD(J,W))/24/3.6
                            ELSE
                                FLOWOUT(J,I,W) = Demand(J,W)			! case 4
                            END IF
                        END IF
                        IF (FLOWOUT(J,I,W)>QRESER(J,W)) THEN			! double check just in case users chose unrealistic x1 and x4
                            FLOWOUT(J,I,W) = QRESER(J,W)
                        END IF
                        VOL(J,I+1,W) = VOL(J,I,W) + (FLOWIN(J,I,W)-FLOWOUT(J,I,W))*24*3.6
                        GOTO 123
c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx            
c     STEP 5-6: Reservoir Operation using Strategy 6
c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                    ELSE IF (RULE(J,W) .EQ. 6) THEN
                        IF (FLOWIN(J,I,W)-(OP6(J,I+STARTDAY(J,W),W)-VOL(J,I,W))/24/3.6 >0) THEN
                            VOL(J,I+1,W) = OP6(J,I+STARTDAY(J,W),W)
                            FLOWOUT(J,I,W) = FLOWIN(J,I,W)-(VOL(J,I+1,W)-VOL(J,I,W))/24/3.6
                            
                        ELSE
                            FLOWOUT(J,I,W) = Qmin(J,I,W)
                            VOL(J,I+1,W) = VOL(J,I,W) + (FLOWIN(J,I,W)-FLOWOUT(J,I,W))*24*3.6
                        END IF
                        GOTO 123
 123                END IF   ! End of Loop for Operation Strategies (Step 5)
c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx            
c     STEP 5-7: Calculate Reservoir Variables considering Irrigation and Seepage
c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                    ! Check if there are any negative values
                    IF (ENERGYPRO(J,I,W)<0) THEN
                        ENERGYPRO(J,I,W)=0
                    END IF
                    IF (FLOWOUT(J,I,W)<0) THEN
                        FLOWOUT(J,I,W)=Qmin(J,I,W)
                    END IF
                    IF (VOL(J,I,W)<0)THEN ! Not allow dropping below the minimum water level (mostly due to evaporation)
                        VOL(J,I,W)=0
                    END IF
c                   Remote water for irrigation
                    IF (FLOWOUT(J,I,W)>=IRRIGATION(J,II)) THEN
                        FLOWOUT(J,I,W) = FLOWOUT(J,I,W) - IRRIGATION(J,II)
                    ELSE
                        FLOWOUT(J,I,W) = Qmin(J,I,W)
                    END IF
                    
                    
                    ! Calculate energy production
                    IF (VOL(J,I+1,W)<=VRESER(J,W,I)) THEN
                        IF (RULE(J,W) .EQ. 4) THEN  ! Force the hydropower to follow the input release (even if VOL>VRESER)
                            ENERGYPRO(J,I,W) = 0.9 * 9.81 * FLOWOUT_TURB(J,I,W)    
                        ELSE    
                            ENERGYPRO(J,I,W) = 0.9 * 9.81 * FLOWOUT(J,I,W)
                        END IF 
                    ELSE
                        IF (RULE(J,W) .EQ. 4) THEN  ! Force the hydropower to follow the input release (even if VOL>VRESER)
                            ENERGYPRO(J,I,W) = 0.9 * 9.81 * FLOWOUT_TURB(J,I,W)   
                        ELSE 
                            ENERGYPRO(J,I,W) = 0.9 * 9.81 * QRESER(J,W)
                        END IF    
                    END IF
                    
                    ! Update water losses due to seepage and infiltration
                    ! Infiltration is permanent losses + water seepage is added to outflow
                    ! Note seepage occurs until there is no water left (considering dead volume also)
                
                    IF (VOL(J,I+1,W) - (SEEPAGE(J,W)+INFILTRATION(J,W))*24*3.6 .GT. 0) THEN
                        VOL(J,I+1,W) = VOL(J,I+1,W) - (SEEPAGE(J,W)+INFILTRATION(J,W))*24*3.6
                    !Note that seepage does not contribute to energy production, so we add here
                        FLOWOUT(J,I,W) = FLOWOUT(J,I,W) + SEEPAGE(J,W)
                    ELSE
                        VOL(J,I+1,W) = 0
                    END IF
                    
c                   ! Check the neccesity to spill water
                    IF (VOL(J,I+1,W)>VRESER(J,W,I)) THEN
                        FLOWOUT(J,I,W) = FLOWOUT(J,I,W)+(VOL(J,I+1,W)-VRESER(J,W,I))/24/3.6
                        VOL(J,I+1,W) = VRESER(J,W,I)
                    END IF
                    HHO(J,I,W)=VOL(J,I,W)/VRESER(J,W,I)*(HRESERMAX(J,W)-H0(J,W))+H0(J,W)
                    HHO(J,I+1,W)=VOL(J,I+1,W)/VRESER(J,W,I)* (HRESERMAX(J,W)-H0(J,W))+H0(J,W)
                    ! Note: hydraulic head calculated from the maximum water level
                    IF (HYDROPOWER_GENERATOR(J,W) .EQ. 0) THEN
                        ENERGYPRO(J,I,W) = 0     ! If the turbine is not connected to generator (generator shutdown)
                    ELSE
                        ENERGYPRO(J,I,W) = ENERGYPRO(J,I,W) *((HHO(J,I,W)+HHO(J,I+1,W))/2-(HRESERMAX(J,W)-REALHEAD(J,W)))/1000			! this part is for hydropower production estimation, ignore if work with irrigation reservoirs   
                    END IF 
c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx            
c     STEP 6: Propagate water to the downstream reservoir, considering the time lag
c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                    
 125                IF (CURRENTYEAR<OPEYEAR(J,W)) THEN
                        FLOWIN(J,I,W) = RESFLOWS(J,II,W)
                        FLOWOUT(J,I,W) = FLOWIN(J,I,W)
                        VOL(J,I,W) =0
                        ENERGYPRO(J,I,W) = 0
                        HTK(J,I,W) = 0
                        HHO(J,I,W) = 0
                        VRESER(J,W,I) = 0 
                        REALHEAD(J,W) = 0
                        QRESER(J,W) = 0
                        FLOWOUT_TURB(J,I,W)=0.0  
                    ELSE
                        IF (FLOWOUT(J,I,W) .LE. Qmin(J,I,W)) THEN
                            FLOWOUT(J,I,W) = Qmin(J,I,W)
                        END IF   
                    END IF
                    
                    ! Reading UH_R files for reservoirs
                    RESID=RES_DIRECT(J,1,W)
                    DS_RESID=RES_DIRECT(J,2,W)
                    WRITE(RESID_STR,10) RESID
                    WRITE(DS_RESID_STR,10) DS_RESID
                    ! check if the downstream reservoir is part of the current sub basin (not downstream of the outlet station)
                    RES_IN_BASIN = .FALSE.
                    DO K= 1, NORESERVOIRS(W)
                        IF (RES_DIRECT(K,1,W) == DS_RESID) THEN
                            RES_IN_BASIN = .TRUE.
                            EXIT
                        END IF
                    END DO

                   

                    IF (DS_RESID>0 .AND. RES_IN_BASIN) THEN
                        UH_NAME = 'RES'//trim(adjustl(RESID_STR))//"_"//'RES'//trim(adjustl(DS_RESID_STR))
                        OPEN(22,file = UH_PATH(1:UH_CLEN)//trim(ADJUSTL(UH_NAME))//'.uh_r', status='old')
                        READ(22, *) (UH_RR(1,K), K = 1,KE+UH_DAY-1)  
                        CLOSE(22)
                    ELSE
                        UH_NAME = 'RES'//trim(adjustl(RESID_STR))//"_"//trim(adjustl(NAME5(W)))  
                        OPEN(22,file = UH_PATH(1:UH_CLEN)//trim(ADJUSTL(UH_NAME))//'.uh_r', status='old')
                        READ(22, *) (UH_RR(1,K), K = 1,KE+UH_DAY-1)  
                        CLOSE(22)
                    END IF        
                    
                    ! Find the index of the downstream reservoir --> RESORD(W)
                    IF (DS_RESID > 0 .AND. RES_IN_BASIN) THEN
                        DO K = 1, NORESERVOIRS(W)
                            IF (RES_DIRECT(K,1,W) == DS_RESID) THEN
                                RESORD(W) = K
                                EXIT
                            END IF
                        END DO
                    ELSE
                        RESORD(W) = NORESERVOIRS(W)  ! Route to station
                    END IF
                    ! Added by HD-May25 to account for UH between reservoirs
                    ! RESFLOWS(NORESERVOIRS(W),II,W) is the reservoir at the downstream reservoir.
                    ! Update the flow for the downstream reservoir for the next time step (so the new upstream flow can be accounted for in the next step; this is because the J loop is happening here from downstream to upstream)
                    DO KK = 1,KE+UH_DAY-1
                        IF ((II-KK+1) .GE. 1) THEN 
                            RESFLOWS(RESORD(W), II+1, W) = RESFLOWS(RESORD(W), II+1, W) + UH_RR(1, KK) * FLOWOUT(J, I-KK+1, W)
                        END IF
                    END DO

                    ! Option to consider uniform UH with direct flow of runoff using the Travel Time
                    !IF (RES_DIRECT(J,2,W) .EQ. 0) THEN  ! No downstream reservoirs
                       ! RESFLOWS(NORESERVOIRS(W),II+1+INT(TRVLTIME(J,W)),W)= RESFLOWS(NORESERVOIRS(W),II+1+INT(TRVLTIME(J,W)),W) + FLOWOUT(J,I,W)     
                    !ELSE
                        !RESFLOWS(RESORD(W),II+1+INT(TRVLTIME(J,W)),W) =RESFLOWS(RESORD(W),II+1+INT(TRVLTIME(J,W)),W) + FLOWOUT(J,I,W)
                    !END IF          
                END DO     !!! END OF J LOOP OVER RESERVOIRS of STATION [W] 
            
                IF (RESFLOWS(NORESERVOIRS(W),II,W) .LT. 0) THEN
                    RESFLOWS(NORESERVOIRS(W),II,W) = 0
                END IF
                ! Update Flow array (naturalized flow) to include regulated flow from reservoirs
                FLOW(II,W) = RESFLOWS(NORESERVOIRS(W),II,W)   ! change I(step) to II(simulated_day) 
                
c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx            
c     STEP 7: Save Reservoir States/Fluxes for STEPBYSTEP Version
c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                IF (STEPBYSTEP) THEN
                    write(W_STR,10) W ! convert integer to string to be able to concatenate in filename with strings
                    write(NDAY_SIM_STR,10) (NDAY_SIM)
                    write(NDAY_SIM_STR_NEXT,10) (NDAY_SIM+1)
                    write(NDAY_SIM_STR_PREV,10) (NDAY_SIM-1)
c                   !SAVE RESFLOWS TO READ THEM LATER for StepByStep Version 
                    RESFLOWSTRING_SIM = trim(trim('NF'//trim(ADJUSTL(W_STR)))//'_STEP'//trim(adjustl(NDAY_SIM_STR_NEXT))//'.txt')
                    open(73, file = NF_PATH(1:NF_CLEN)//'/stepbystep/'//RESFLOWSTRING_SIM)
                    DO N = 1,NORESERVOIRS(W)
                        WRITE(73, *) (RESFLOWS(N,K,W), K = 1,NDAY)
                    END DO
                    close(73)

                    ! DELETE PREVIOUS NATURALIZED FLOW FILE [only keep the first step for future simulations]
                    IF (NDAY_SIM > 2) then
                        previous_natflow_filename = trim(trim('NF'//trim(ADJUSTL(W_STR)))//'_STEP'//trim(adjustl(NDAY_SIM_STR_PREV))//'.txt')
                        inquire(file=NF_PATH(1:NF_CLEN)//'/stepbystep/'//previous_natflow_filename, exist=file_exists)
                        IF (file_exists) then
                            open(74, file=NF_PATH(1:NF_CLEN)//'/stepbystep/'//previous_natflow_filename, status='old')
                            close(74, status='delete')
                        end IF
                    END IF

                    !SAVE Storage TO READ THEM LATER for StepByStep Version
                    IF (NDAY_SIM .EQ. 1)  THEN   ! Save for the first timestep so it can be used in case of iteartions of conevrgence
                        storage_filename='RESERVOIR_VOL_STATION'//trim(adjustl(W_STR))//'_STEP'//trim(adjustl(NDAY_SIM_STR))//'.txt'     ! needed for convergence iterations
                        open(88, file=OUTPATH(1:OUT_CLEN)//'/stepbystep/'//storage_filename)
                        DO J = 1, NORESERVOIRS(W)
                            WRITE(88, *) (VOL(J,I,W))
                        END DO
                        close(88)
                        storage_filename='RESERVOIR_VOL_STATION'//trim(adjustl(W_STR))//'_STEP'//trim(adjustl(NDAY_SIM_STR_NEXT))//'.txt'     ! needed for next time step
                        open(87, file=OUTPATH(1:OUT_CLEN)//'/stepbystep/'//storage_filename)
                        DO J = 1, NORESERVOIRS(W)
                            WRITE(87, *) (VOL(J,I+1,W))
                        END DO
                        close(87)
                    ELSE
                        storage_filename='RESERVOIR_VOL_STATION'//trim(adjustl(W_STR))//'_STEP'//trim(adjustl(NDAY_SIM_STR_NEXT))//'.txt'     ! needed for next time step
                        open(87, file=OUTPATH(1:OUT_CLEN)//'/stepbystep/'//storage_filename)
                        DO J = 1, NORESERVOIRS(W)
                            WRITE(87, *) (VOL(J,I+1,W))
                        END DO
                        close(87)
                    ENDIF
                    ! DELETE PREVIOUS STORAGE FILE
                    IF (NDAY_SIM > 1) then
                        previous_storage_filename = 'RESERVOIR_VOL_STATION'//trim(adjustl(W_STR))//'_STEP'//trim(adjustl(NDAY_SIM_STR_PREV))//'.txt'
                        inquire(file=OUTPATH(1:OUT_CLEN)//'/stepbystep/'//previous_storage_filename, exist=file_exists)
                        IF (file_exists) then
                            open(89, file=OUTPATH(1:OUT_CLEN)//'/stepbystep/'//previous_storage_filename, status='old')
                            close(89, status='delete')
                        END IF
                    END IF
                END IF
            ENDDO       !!! END OF W LOOP OVER STATIONS               
        ! At the end of this loop, we will have the discharge for each station under consideration at time i, in FLOW(I,:)
        ENDDO    !!! END OF I LOOP OVER STEPS (STEPS=1 For SBS Version) or (STEPS=NDAY For All Steps)
        RETURN
 9001   WRITE(*,*) 'Error reading UH data'
 9002   WRITE(*,*) 'Error in reading reservoir data'
 9003   WRITE(*,*) 'Error in reading irrigation data'
 9004   WRITE(0,*) 'ERROR in Reading Storage Volume for Current Time Step; Please Check if Previous Time Step is Simulated'     
        END
C************************************************************************************************************************************************************************************
C       Convert from day to month
C************************************************************************************************************************************************************************************
        INTEGER FUNCTION CAL_MONTH(CRTDATE)
        INTEGER CRTDATE
        IF (CRTDATE .LE. 31) THEN !jan
            CAL_MONTH = 1
        ELSE IF (CRTDATE .LE. 59) THEN ! feb
            CAL_MONTH = 2
        ELSE IF (CRTDATE .LE. 90) THEN ! mar
            CAL_MONTH = 3
        ELSE IF (CRTDATE .LE. 120) THEN ! apr
            CAL_MONTH = 4
        ELSE IF (CRTDATE .LE. 151) THEN ! may
            CAL_MONTH = 5
        ELSE IF (CRTDATE .LE. 181) THEN ! jun
            CAL_MONTH = 6
        ELSE IF (CRTDATE .LE. 212) THEN ! jul
            CAL_MONTH = 7
        ELSE IF (CRTDATE .LE. 243) THEN ! aug
            CAL_MONTH = 8
        ELSE IF (CRTDATE .LE. 273) THEN ! sep
            CAL_MONTH = 9
        ELSE IF (CRTDATE .LE. 304) THEN ! oct
            CAL_MONTH = 10
        ELSE IF (CRTDATE .LE. 334) THEN ! nov
            CAL_MONTH = 11
        ELSE ! dec
            CAL_MONTH = 12
        END IF
        RETURN
        END
