c       Contain main subroutines of the routing component and reservoir operation
c       Build based on the old make_convolution file 
c       Note: currently the maximum number of reservoirs is NUMRES. Modify the variable declaration parts to increase (if needed)
c       28 Dec 2020 - Correct the bug related to the dead storage - Thank Mr. Arivoli and his collegues

C************************************************************************************************************************************************************************************
C       Read diffusion file

        SUBROUTINE READ_DIFF(DIFF,NCOL,NROW,FILENAME,IROW, ICOL)
c       Declare variables
        INTEGER NCOL,NROW,IROW,ICOL,I,J
        REAL DIFF(NCOL,NROW)
        CHARACTER*72 FILENAME
c       Subroutine main body
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
     &    (UH_BOX,KE,PMAX,NOB,CATCHIJ,FILENAME,NORESERVOIRS, NUMRES)
        IMPLICIT NONE
c       Declare variables
        INTEGER KE, PMAX, NUMRES
        INTEGER NOB(NUMRES)
        INTEGER CATCHIJ(PMAX,2,NUMRES)
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
        INTEGER DIREC(NCOL,NROW,2)
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
     &      MO, YR, NYR, VOL, FLOWIN, FLOWOUT, HHO, RESFLOWS, NO, resn, RES_EVAPORATION, NO_STAS, NUMRES)
        IMPLICIT NONE
c       Declare variables
        INTEGER     N, NO, resn, I, J, DAYS, NDAY, II, JJ, K, SODONG, NUMRES
        INTEGER     NCOL,NROW,ICOL,PMAX,KE,UH_DAY
        INTEGER     NO_OF_BOX(NUMRES)
        INTEGER     CATCHIJ(PMAX,2,NUMRES)
        INTEGER     NYR, START_YEAR, START_MO
        INTEGER     RESER(NCOL,NROW)
        INTEGER     IDAY(DAYS), IMONTH(DAYS), IYEAR(DAYS)
        INTEGER     MO(12*NYR),YR(12*NYR)
        INTEGER     MISS_NUM
        INTEGER     CURRENTYEAR
        INTEGER     MONTH_OF_YEAR
        INTEGER     DPREC, CLEN
        INTEGER     OPEYEAR, NO_STAS
        REAL        RES_EVAPORATION(NUMRES,DAYS)
        REAL        UH_S(PMAX,KE+UH_DAY-1,NUMRES)
        REAL        BASE(DAYS), RUNO(DAYS), FLOW(DAYS), AIRTEMP(DAYS), WINDS(DAYS), RHD(DAYS)
        REAL        FRACTION(NCOL,NROW)
        REAL        VRESER(NUMRES), QRESER(NUMRES)
        REAL        HRESERMAX(NUMRES), HRESERMIN(NUMRES)
        REAL        V0, FLOWINN
        REAL        PI, RERD, FACTOR, FACTOR_SUM
        REAL        DESIGNWL, CURRENTWL
        REAL        JLOC, ILOC, EVA_FACTOR
        REAL        XC, YC, SIZE
        REAL        AREA, AREA_SUM
        REAL        VOL(NUMRES,DAYS), FLOWIN(NUMRES,DAYS), FLOWOUT(NUMRES,DAYS)
        REAL        HHO(NUMRES,DAYS)
        REAL        K_CONST
        REAL        DUM1,DUM2,DUM3,DUM4,DUM5,DUM6,DUM7,DUM8,DUM9,DUM10,DUM11,DUM12,DUM13
        REAL        RESFLOWS(NUMRES,DAYS)
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
c       Weighting factor delta/gamma matrix for the Penman formula
        weightingmatrix = (/0.69, 0.71, 0.73, 0.76, 0.78, 0.80, 0.83, 0.85, 0.88, 0.91, 0.94, 0.97,
     &                     1.00, 1.03, 1.06, 1.09, 1.13, 1.16, 1.19, 1.23, 1.26, 1.30, 1.34, 1.38, 
     &                     1.42, 1.47, 1.51, 1.56, 1.60, 1.65, 1.69, 1.74, 1.79, 1.84, 1.89, 1.94,
     &                     1.99, 2.04, 2.11, 2.17, 2.23, 2.28, 2.35, 2.41, 2.48/)
C       Ea matrix
        eamatrix = (/4.58, 4.75, 4.93, 5.11, 5.30, 5.49, 5.69, 5.89, 6.10, 6.32, 6.54, 6.77, 7.02, 7.26,
     &                     7.52, 7.79, 8.05, 8.33, 8.62, 8.91, 9.21, 9.52, 9.84, 10.18, 10.52, 10.87, 11.23, 
     &                     11.61, 11.99, 12.38, 12.79, 13.13, 13.63, 14.08, 14.53, 15.00, 15.49, 15.97, 16.47,
     &                     17.53, 18.68, 19.80, 21.10/)
C       DeltaTa4matrix for the Penman formula
        deltaTa4matrix = (/11.0, 11.2, 11.4, 11.5, 11.7, 11.9, 12.0, 12.2, 12.3, 12.5, 12.7, 12.9, 13.1,
     &                     13.3, 13.5, 13.7, 13.9, 14.0, 14.2, 14.4, 14.6, 14.8, 15.0/)
C       MISS_NUM is the number of grid cell output files not found
        MISS_NUM=0
C       *** 0 <= K_CONST = 1.0
C ***   K_CONST smaller 1.0 makes it a simple linear storage
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
C       Look for starting year
        IF (NO .NE. NO_STAS) THEN
            WRITE(RESNO,*) resn
            TEMPRPATH = trim(RPATH)//"res"//trim(ADJUSTL(RESNO))//".txt"			! only calculate open surface water evaporation after the commision year
            OPEN(26, FILE = TEMPRPATH,FORM = 'FORMATTED', STATUS='OLD',ERR=9002)
            READ(26,*)
            READ(26,*) DUM1,DUM2,DUM3,DUM4,DUM5,DUM6, OPEYEAR
            CLOSE(26)
        END IF


        
C       Calculate the area of each cell        
        DO N = 1,NO_OF_BOX(NO) 
            DO I = 1,NDAY
                RUNO(I) = 0.0
                BASE(I) = 0.0
            END DO
            II = CATCHIJ(N,1,NO)
            JJ = CATCHIJ(N,2,NO)
c       the grid has been flipped left to right
c       find the revised cooordinates
            ILOC=XC + (ICOL-II)*SIZE + SIZE/2.0
            JLOC=YC + JJ*SIZE - SIZE/2.0
C       CONVERSIONFACTOR for mm/day to ft**3/sec
        AREA =  RERD**2*ABS(SIZE)*PI/180*										! give area of box in square meters
     &          ABS(SIN((JLOC-SIZE/2.0)*PI/180)-
     $          SIN((JLOC+SIZE/2.0)*PI/180))
        AREA_SUM = AREA_SUM + AREA
c       WRITE(*,*) N, ILOC, JLOC
        FACTOR = FRACTION(II,JJ)*35.315*AREA/(86400.0*1000.0)  				!convert to sq.mi. by cell fract (original) - later, we convert back to SI unit
        FACTOR_SUM = FACTOR_SUM + FACTOR
        call create_vic_names(jloc,iloc,loc,clen,dprec)
        INQUIRE(FILE=INPATH(1:(INDEX(INPATH,' ')-1))//LOC(1:CLEN),
     $      EXIST=TORF)
        IF(torf)THEN
        OPEN(20,FILE=INPATH(1:(INDEX(INPATH,' ')-1))//
     $      LOC(1:CLEN),
     $      STATUS='OLD',ERR=9001)
        DO I = 1,NDAY 
            READ(20,*,END=9001,ERR=9001) IYEAR(I), IMONTH(I), IDAY(I),
     &      DUM1, DUM2, RUNO(I), BASE(I), DUM3, DUM4, DUM5, DUM6, DUM7, vic_eva_vege(I),
     &      vic_trans_vege(I), vic_eva_soil(I), DUM8, DUM9, DUM10, DUM11, DUM12, RHD(I), DUM13,
     &      AIRTEMP(I),WINDS(I) 												! (IMPORTANT) modify here accordingly depending on your selected outputs from rainfall-runoff
c       check to be sure dates in VIC file start at same time specified
c       in input file
            IF(I.eq.1) THEN
                IF(IYEAR(I).ne.YR(1) . or. IMONTH(I).ne.MO(1)) THEN
                    print*, 'VIC output file does not match specified '
                    print*, 'period in input file.'
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
        if(MISS_NUM>0) then
            print*, MISS_NUM, ' files not found, zero runoff/baseflow used'
        end if
        RETURN
 9001   WRITE(*,*) 'Error reading time-series data, ',
     $     'insufficient data or missing input file',
     $     INPATH(1:INDEX(INPATH,' ')-1)//LOC(1:CLEN)
 9002   WRITE(*,*) 'Error in reading reservoir data'
        END

C************************************************************************************************************************************************************************************
C       Determine the distance/time from the downstream cell of a reservoir to its cascade reservoir
C************************************************************************************************************************************************************************************
        SUBROUTINE CALCULATE_DISTANCE(SI,SJ,ICOL,IROW,NCOL,NROW,DIREC,TRVLTIME,
     &  RESORDER,RESER,FINAL,NORESERVOIRS,SIZE_CELL,VELO,FI,FJ,NUMRES)
        IMPLICIT NONE
c       Declare variables  
        INTEGER SI,SJ,II,JJ,COUNTCELL,FINAL,NCOL,NROW,NORESERVOIRS,FI,FJ,NUMRES
        INTEGER DIREC(NCOL,NROW,2)
        INTEGER RESER(NCOL,NROW)
        INTEGER I,J,ICOL,IROW,III,JJJ,RESORDER
        REAL    VELO(NCOL,NROW)
        REAL    TRVLTIME(NUMRES)
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
     $     IROW,ICOL,NORESERVOIRS,RES_DIRECT,RESER,FINAL,SIZE_CELL,VELO,FI,FJ,NUMRES)
        IMPLICIT NONE
c       Declare variables
        INTEGER PI,PJ,I,J,NCOL,NROW,PMAX,ICOL,IROW,NUMRES
        INTEGER II, JJ, III, JJJ, K,FI,FJ
        INTEGER DIREC(NCOL,NROW,2)
        INTEGER NO_OF_BOX(NUMRES)
        INTEGER CATCHIJ(PMAX,2,NUMRES)
        INTEGER RESER(NCOL,NROW)
        INTEGER RES_DIRECT(NUMRES,3)
        INTEGER NORESERVOIRS,RESORDER
        INTEGER FINAL
        REAL    SIZE_CELL
        REAL    TRVLTIME(NUMRES)
        REAL    VELO(NCOL,NROW)
c       Subroutine main body
        RESORDER = NORESERVOIRS
        DO I = 1, NORESERVOIRS
            IF (RESER(PI,PJ).EQ. RES_DIRECT(I,1)) THEN
               RESORDER = I
            END IF
        END DO
        NO_OF_BOX(RESORDER)=0
        If (FINAL .EQ. 1) THEN
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
                WRITE(*,*) 'Cell (',I,',',J,')---> Res ID ',RES_DIRECT(RESORDER,1)				! show which cells go to which reservoir
                IF (FINAL .EQ. 0) THEN
                    GOTO 310 
                ELSE
                    GOTO 320
                END IF
            ELSE IF ((FINAL .EQ. 1) .AND. (RESER(II,JJ)>0)) THEN
                GOTO 320
            ELSE IF ((RESER(II,JJ)>0) .AND.(RESER(II,JJ) .NE. 9999)) THEN
                GOTO 310
            ELSE
                IF ((DIREC(II,JJ,1).NE.0) .AND.    								!check if the current
     &              (DIREC(II,JJ,2) .NE.0)) THEN   								!ii,jj cell routes down
                    III = DIREC(II,JJ,1)         								!to the subbasin outlet
                    JJJ = DIREC(II,JJ,2)         								!point, following the
                    II  = III                    								!direction of direc(,)
                    JJ  = JJJ                    								!from each cell
                    GOTO 300
                END IF   														!if you get there,
            END IF                                								!no_of_box increments
 310        CONTINUE                              								!and you try another
            IF ((I .GE. 2) .AND. (I .LT. (ICOL-1)).AND. (J .GE. 2)				! this part to track which cell a reservoir releases to
     &          .AND. (J .LT. (IROW-1))) THEN
                IF ((RESER(I-1,J) .GT. 0) .AND. (DIREC(I-1,J,1) .EQ. I)
     &             .AND. (DIREC(I-1,J,2) .EQ. J) .AND. (RESER(I-1,J) .NE. 9999)) THEN
                    DO K = 1, NORESERVOIRS
                        IF (RESER(I-1,J) .EQ. RES_DIRECT(K,1)) THEN
                            RES_DIRECT(K,3) =  NO_OF_BOX(RESORDER)
                            CALL CALCULATE_DISTANCE(I,J,ICOL,IROW,NCOL,NROW,DIREC,TRVLTIME,
     &                            RESORDER,RESER,FINAL,NORESERVOIRS,SIZE_CELL,VELO,FI,FJ,NUMRES)
                        END IF
                    END DO
                END IF
                IF ((RESER(I+1,J) .GT.0) .AND. (DIREC(I+1,J,1) .EQ. I)
     &              .AND. (DIREC(I+1,J,2) .EQ. J) .AND. (RESER(I+1,J) .NE.9999)) THEN
                    DO K = 1, NORESERVOIRS
                        IF (RESER(I+1,J) .EQ. RES_DIRECT(K,1)) THEN
                            RES_DIRECT(K,3) =  NO_OF_BOX(RESORDER)
                            CALL CALCULATE_DISTANCE(I,J,ICOL,IROW,NCOL,NROW,DIREC,TRVLTIME,
     &                      RESORDER,RESER,FINAL,NORESERVOIRS,SIZE_CELL,VELO,FI,FJ,NUMRES)
                        END IF
                    END DO
                END IF
                IF ((RESER(I,J-1) .GT.0) .AND. (DIREC(I,J-1,1) .EQ. I)
     &            .AND. (DIREC(I,J-1,2) .EQ. J) .AND. (RESER(I,J-1) .NE.9999)) THEN
                    DO K = 1, NORESERVOIRS
                        IF (RESER(I,J-1) .EQ. RES_DIRECT(K,1)) THEN
                            RES_DIRECT(K,3) =  NO_OF_BOX(RESORDER)
                            CALL CALCULATE_DISTANCE(I,J,ICOL,IROW,NCOL,NROW,DIREC,TRVLTIME,
     &                      RESORDER,RESER,FINAL,NORESERVOIRS,SIZE_CELL,VELO,FI,FJ,NUMRES)
                        END IF
                    END DO
                END IF
                IF ((RESER(I,J+1) .GT.0) .AND. (DIREC(I,J+1,1) .EQ. I)
     &            .AND. (DIREC(I,J+1,2) .EQ. J) .AND. (RESER(I,J+1) .NE.9999)) THEN
                    DO K = 1, NORESERVOIRS
                        IF (RESER(I,J+1) .EQ. RES_DIRECT(K,1)) THEN
                            RES_DIRECT(K,3) =  NO_OF_BOX(RESORDER)
                            CALL CALCULATE_DISTANCE(I,J,ICOL,IROW,NCOL,NROW,DIREC,TRVLTIME,
     &                      RESORDER,RESER,FINAL,NORESERVOIRS,SIZE_CELL,VELO,FI,FJ,NUMRES)
                        END IF
                    END DO
                END IF
                IF ((RESER(I-1,J-1) .GT.0) .AND. (DIREC(I-1,J-1,1) .EQ. I)
     &          .AND. (DIREC(I-1,J-1,2) .EQ. J) .AND. (RESER(I-1,J-1) .NE.9999)) THEN
                    DO K = 1, NORESERVOIRS
                        IF (RESER(I-1,J-1) .EQ. RES_DIRECT(K,1)) THEN
                            RES_DIRECT(K,3) =  NO_OF_BOX(RESORDER)
                            CALL CALCULATE_DISTANCE(I,J,ICOL,IROW,NCOL,NROW,DIREC,TRVLTIME,
     &                      RESORDER,RESER,FINAL,NORESERVOIRS,SIZE_CELL,VELO,FI,FJ,NUMRES)
                        END IF
                    END DO
                END IF
                IF ((RESER(I-1,J+1) .GT.0) .AND. (DIREC(I-1,J+1,1) .EQ. I)
     &          .AND. (DIREC(I-1,J+1,2) .EQ. J) .AND. (RESER(I-1,J+1) .NE.9999)) THEN
                    DO K = 1, NORESERVOIRS
                        IF (RESER(I-1,J+1) .EQ. RES_DIRECT(K,1)) THEN
                            RES_DIRECT(K,3) =  NO_OF_BOX(RESORDER)	
                            CALL CALCULATE_DISTANCE(I,J,ICOL,IROW,NCOL,NROW,DIREC,TRVLTIME,
     &                      RESORDER,RESER,FINAL,NORESERVOIRS,SIZE_CELL,VELO,FI,FJ,NUMRES)
                        END IF
                    END DO
                END IF
                IF ((RESER(I+1,J-1) .GT.0) .AND. (DIREC(I+1,J-1,1) .EQ. I)
     &          .AND. (DIREC(I+1,J-1,2) .EQ. J) .AND. (RESER(I+1,J-1) .NE.9999)) THEN
                    DO K = 1, NORESERVOIRS
                        IF (RESER(I+1,J-1) .EQ. RES_DIRECT(K,1)) THEN
                            RES_DIRECT(K,3) =  NO_OF_BOX(RESORDER)
                            CALL CALCULATE_DISTANCE(I,J,ICOL,IROW,NCOL,NROW,DIREC,TRVLTIME,
     &                      RESORDER,RESER,FINAL,NORESERVOIRS,SIZE_CELL,VELO,FI,FJ,NUMRES)
                        END IF
                    END DO
                END IF
                IF ((RESER(I+1,J+1) .GT.0) .AND. (DIREC(I+1,J+1,1) .EQ. I)
     &           .AND. (DIREC(I+1,J+1,2) .EQ. J) .AND. (RESER(I+1,J+1) .NE.9999)) THEN
                    DO K = 1, NORESERVOIRS
                        IF (RESER(I+1,J+1) .EQ. RES_DIRECT(K,1)) THEN
                            RES_DIRECT(K,3) =  NO_OF_BOX(RESORDER)
                            CALL CALCULATE_DISTANCE(I,J,ICOL,IROW,NCOL,NROW,DIREC,TRVLTIME,
     &                      RESORDER,RESER,FINAL,NORESERVOIRS,SIZE_CELL,VELO,FI,FJ,NUMRES)
                        END IF
                    END DO
                END IF
            END IF
320         CONTINUE
        END DO
        END DO
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
        WRITE(*,*) 'Res ',RES_DIRECT(RESORDER,1),'---> Res ',RES_DIRECT(RESORDER,2)
        WRITE(*,*) 'Number of grid cells upstream of present reservoir',
     $      no_of_box(RESORDER)
        RETURN
        END

C************************************************************************************************************************************************************************************
C     Consider rule curves and opearting rules
C************************************************************************************************************************************************************************************

        SUBROUTINE MAKE_CONVOLUTIONRS
     & (RPATH, RESER, NCOL, NROW, NO_OF_BOX, PMAX, DAYS, CATCHIJ,
     &  BASE, RUNO, FLOW, KE, UH_DAY, FRACTION, FACTOR_SUM,
     &  XC, YC, SIZE, DPREC, INPATH,ICOL,NDAY,IDAY,IMONTH,IYEAR, START_YEAR,START_MO,
     &  MO, YR, NYR, VOL, FLOWIN, FLOWOUT, HHO, ENERGYPRO, HTK, DIREC, IROW,
     &  XNE, YNE, NO_STAS, RES_DIRECT, RES_EVAPORATION, TRVLTIME,RESORDER,NAMERS5,NAME5,NREACH,UH_SS,VRESER,NUMRES, NSTATIONS)
        IMPLICIT NONE
c       Declare variables
        INTEGER     N, DAYS, NDAY, START_YEAR,START_MO,NREACH,NUMRES, NSTATIONS
        INTEGER     NCOL,NROW,ICOL,PMAX,KE,UH_DAY
        INTEGER     CATCHIJ(PMAX,2,NUMRES,NREACH)
        INTEGER     NYR
        INTEGER     RESER(NCOL,NROW)
        INTEGER     IDAY(DAYS), IMONTH(DAYS), IYEAR(DAYS)
        INTEGER     MO(12*NYR),YR(12*NYR)
        INTEGER     DIREC(NCOL,NROW,2)
        INTEGER     RES_DIRECT(NUMRES,3,NREACH)
        INTEGER     YEAR(NUMRES,DAYS),MONTH(NUMRES,DAYS),DAY(NUMRES,DAYS)
        INTEGER     OPEYEAR(NUMRES,NREACH)
        INTEGER     IROW, CURRENTYEAR
        INTEGER     NO_OF_BOX(NUMRES,NREACH)
        INTEGER     XNE(NREACH), YNE(NREACH)
        INTEGER     NO_STAS(NREACH)
        INTEGER     I, J, K, L, CAL_MONTH, W, x
        INTEGER     NO_OF_ROW(NUMRES,NREACH)
        INTEGER     RULE(NUMRES,NREACH), WITHIRRIGATION(NUMRES,NREACH)
        INTEGER     STARTDAY(NUMRES,NREACH)
        INTEGER     DPREC, RESORDER(NREACH,NUMRES), NOOFROW(NREACH), CRTDATE,RESORD(NREACH)
        REAL        UH_SS(PMAX,KE+UH_DAY-1,NUMRES)
        REAL        TRVLTIME(NUMRES)
        REAL        BASE(DAYS,NREACH), RUNO(DAYS,NREACH)
        REAL        FLOW(DAYS,NREACH)
        REAL        FRACTION(NCOL,NROW)
        REAL        PI, RERD, FACTOR, FACTOR_SUM
        REAL        CURRENTWL(NREACH), DESIGNWL(NREACH)
        REAL        RES_EVAPORATION(NUMRES,DAYS)
        REAL        XC, YC, SIZE, TEMP
        REAL        VOL(NUMRES,DAYS,NREACH), FLOWIN(NUMRES,DAYS), FLOWOUT(NUMRES,DAYS,NREACH)
        REAL        HHO(NUMRES,DAYS), ENERGYPRO(NUMRES,DAYS), HTK(NUMRES,DAYS)
        REAL        HMAX(NUMRES,NREACH), HMIN(NUMRES,NREACH)
        REAL        RESFLOWS(NUMRES,DAYS,NREACH) 
        REAL        HRESERMAX(NUMRES,NREACH), HRESERMIN(NUMRES,NREACH), VINITIAL(NUMRES,NREACH), H0(NUMRES,NREACH)
        REAL        QRESER(NUMRES,NREACH),VRESER(NUMRES,NREACH,DAYS),REALHEAD(NUMRES,NREACH),VRESERTHAT(NUMRES,NREACH),VDEAD(NUMRES,NREACH)
        REAL        QRESERTHAT(NUMRES,NREACH), HYDRAUHEAD(NUMRES,NREACH)
        REAL        OP1(NUMRES,2,NREACH) ! rule curve
        REAL        OP2(NUMRES,DAYS) ! pre-defined time series data
        REAL        X1(NUMRES,NREACH), X2(NUMRES,NREACH), X3(NUMRES,NREACH), X4(NUMRES,NREACH)	! operating rule
        REAL        SEEPAGE(NUMRES,NREACH), INFILTRATION(NUMRES,NREACH), Demand(NUMRES,NREACH)
        REAL        TEMPO
        REAL        RESLV(NUMRES,12,NREACH), OP5X1(NUMRES,12,NREACH), OP5X2(NUMRES,12,NREACH), OP5X3(NUMRES,12,NREACH), OP5X4(NUMRES,12,NREACH), DEMAND5(NUMRES,12,NREACH)
        REAL        IRRIGATION(NUMRES,DAYS)
        REAL        KFACTOR(NREACH)
        CHARACTER*10 WTEMP
        CHARACTER*72 RESFLOWSTRING
        CHARACTER*72 RPATH
        CHARACTER*72 IPATH(NREACH)
        CHARACTER*100 TEMPRPATH(NREACH)
        CHARACTER*100 PATHRES3
        CHARACTER*72 INPATH
        CHARACTER*20 CHUOI(NREACH)
        CHARACTER*7  NAMERS5(NREACH,NUMRES)
        CHARACTER*7  NAME5(NREACH)
        CHARACTER*100 FILENAME

c       Subroutine main body    

c       Look for reservoirs sequences
        CURRENTYEAR = START_YEAR
        print*,'reading grid from files uh_s'
        
        DO W=1,NSTATIONS
            DO J = 1,NO_STAS(W)-1
                  
                  OPEN(98, file = trim(adjustl(NAMERS5(W,RESORDER(W,J))))//'.uh_s', status='old')
                  DO N = 1,NO_OF_BOX(RESORDER(W,J),W)
                        PRINT*,NO_OF_BOX(RESORDER(W,J),W)
                        READ(98, *) (UH_SS(N,K,RESORDER(W,J)), K = 1,KE+UH_DAY-1)  
                  ENDDO
                  CLOSE(98) 
            ENDDO
            J=NO_STAS(W)
            print*, W, J,RESORDER(W,J), NAME5(W)
            OPEN(98, file = trim(adjustl(NAME5(W)))//'.uh_s', status='old')
            DO N= 1,NO_OF_BOX(RESORDER(W,J),W)
                READ(98, *) (UH_SS(N,K,RESORDER(W,J)), K = 1,KE+UH_DAY-1)
            ENDDO
            CLOSE(98)
            DO N = 1,NO_STAS(W)
                print*, 'Working on reservoir no...',RES_DIRECT(N,1,W)
                CALL MAKE_CONVOLUTION
     &              (RPATH,RESER,NCOL, NROW, NO_OF_BOX(:,W), PMAX, DAYS,
     &              CATCHIJ(:,:,:,W), BASE(:,W), RUNO(:,W), FLOW(:,W), KE, UH_SS, UH_DAY, FRACTION,
     &              FACTOR_SUM,XC,YC,SIZE,DPREC,INPATH,ICOL,NDAY,
     &              IDAY,IMONTH,IYEAR, START_YEAR, START_MO, MO, YR, NYR, VOL(:,:,W),
     &              FLOWIN, FLOWOUT(:,:,W), HHO, RESFLOWS(:,:,W),N,RES_DIRECT(N,1,W),RES_EVAPORATION, NO_STAS(W),NUMRES)
                
            END DO
            WRITE(WTEMP, 10) W
 10         FORMAT(I4)
            RESFLOWSTRING = trim(trim('RF'//trim(ADJUSTL(WTEMP)))//'.txt')
C           SAVE RESFLOWS TO READ THEM LATER
        
           open(71, file=RESFLOWSTRING, status='unknown')
           DO N = 1,NO_STAS(W)
                WRITE(71, *) (RESFLOWS(N,K,W), K = 1,NDAY)
           END DO
           close(71)
c           READ RESFLOWS
            open(72, file=RESFLOWSTRING, status='old')

            DO N = 1, NO_STAS(W)
                READ(72, *) (RESFLOWS(N,K,W), K = 1, NDAY)
            END DO


            close(72)                
            !print*,RESFLOWS(20,1,W)
c       Initiate reservoir parameters
            DO I = 1, NDAY
                DO J = 1, NO_STAS(W)-1
                    FLOWIN(J,I) = 0
                    FLOWOUT(J,I,W) = 0
                    VOL(J,I,W) = 0
                    ENERGYPRO(J,I) = 0
                    HTK(J,I) = 0
                END DO
            END DO
            DO  J  = 1, NO_STAS(W)-1
                WRITE(CHUOI(W),*) RES_DIRECT(J,1,W)
                TEMPRPATH(W) = trim(RPATH)//"res"//trim(ADJUSTL(CHUOI(W)))//".txt"
                OPEN(25, FILE = TEMPRPATH(W),FORM = 'FORMATTED',
     &          STATUS='OLD',ERR=9002)
                READ(25,*)  
                READ(25,*) HRESERMAX(J,W),HRESERMIN(J,W),VRESERTHAT(J,W),VDEAD(J,W),
     &          HYDRAUHEAD(J,W), QRESERTHAT(J,W), OPEYEAR(J,W), VINITIAL(J,W)
                KFACTOR(W) = SQRT(VDEAD(J,W)/VRESERTHAT(J,W))							! calculate the bottom level of reservoir
                H0(J,W) = (KFACTOR(W)*HRESERMAX(J,W) - HRESERMIN(J,W))/(KFACTOR(W)-1)		! pyramid shape reservoir
c                H0(J,W) = (HRESERMAX(J,W)-HRESERMIN(J,W)*(VRESERTHAT(J,W)+VDEAD(J,W))/VDEAD(J,W))/(1-(VRESERTHAT(J,W)+VDEAD(J,W))/VDEAD(J,W))    !This actually performs better, but we will wix the problem of H0 with reservoirs bathimetry
                READ(25,*)
                READ(25,*) SEEPAGE(J,W), INFILTRATION(J,W)
                READ(25,*)
                READ(25,*) WITHIRRIGATION(J,W)
                READ(25,'(A)') IPATH(W)
                READ(25,*)
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
                READ(25,*) RULE(J,W)
                IF (RULE(J,W) .EQ. 1) THEN ! simplified rule curve
                    READ(25,*) HMAX(J,W), HMIN(J,W), OP1(J,1,W),OP1(J,2,W)
                    IF (HMAX(J,W) .LT. HMIN(J,W)) THEN
                        TEMP=HMAX(J,W)
                        HMAX(J,W)=HMIN(J,W)
                        HMIN(J,W)=TEMP
                    END IF
                ELSE IF (RULE(J,W) .EQ. 2) THEN ! rule curve
                    READ(25,*) RESLV(J,1,W), RESLV(J,2,W), RESLV(J,3,W), RESLV(J,4,W), RESLV(J,5,W), RESLV(J,6,W),
     &              RESLV(J,7,W), RESLV(J,8,W), RESLV(J,9,W), RESLV(J,10,W), RESLV(J,11,W), RESLV(J,12,W)
                ELSE IF (RULE(J,W) .EQ. 3) THEN ! operating rule
                    READ(25,*) Demand(J,W), X1(J,W), X2(J,W), X3(J,W), X4(J,W)
                ELSE IF (RULE(J,W) .EQ. 4) THEN! predefined time-series
                    READ(25,'(/A)')FILENAME
                    PATHRES3 = trim(ADJUSTL(FILENAME))
                    OPEN(26, FILE = PATHRES3,FORM = 'FORMATTED',STATUS='OLD',ERR=9002)
                    READ(26,*) NO_OF_ROW(J,W) 						! Read the number of rows of reservoir J th
                    DO L = 1, NO_OF_ROW(J,W)
                        READ(26,*) YEAR(J,L), MONTH(J,L), DAY(J,L), OP2(J,L)
                        IF ((YEAR(J,L) .EQ. START_YEAR) .AND. (MONTH(J,L) .EQ. START_MO) .AND. (DAY(J,L) .EQ. 1)) THEN
                            STARTDAY(J,W) = L							! If the time-series is shorter than the simulation period (release = 0)
                        END IF
                    END DO
                    CLOSE(26)
                ELSE IF (RULE(J,W) .EQ. 5) THEN
                    DO L = 1, 12
                        READ(25,*) DEMAND5(J,L,W), OP5X1(J,L,W), OP5X2(J,L,W), OP5X3(J,L,W), OP5X4(J,L,W)
                    END DO
                END IF
                CLOSE(25)
                IF (CURRENTYEAR>=OPEYEAR(J,W)) THEN
                    VOL(J,1,W) = VINITIAL(J,W)
                ELSE 
                    VOL(J,1,W) = 0.001
                END IF
            END DO
        ENDDO
        DO I = 1, NDAY
            DO W=1, NSTATIONS
                DO J = 1, NO_STAS(W)-1
                    CURRENTYEAR = START_YEAR + INT(I/365)		! approximate, does not consider leap years
                IF (CURRENTYEAR>=OPEYEAR(J,W)) THEN
                    VRESER(J,W,I) = VRESERTHAT(J,W)+VDEAD(J,W)
                    REALHEAD(J,W) = HYDRAUHEAD(J,W)
                    QRESER(J,W) = QRESERTHAT(J,W)
                ELSE
                    VRESER(J,W,I) = 0.001
                    REALHEAD(J,W) = 0.001
                    QRESER(J,W) = 999999
                END IF
                FLOWIN(J,I) = RESFLOWS(J,I,W)
c               Calculate the designed water level
c               Note RULE = 1: simplified rule curve - 2: rule curve - 3: operating rules: - 4 pre-defined time-series data - 5 12 month operating rule
                CRTDATE = INT(1.0* mod(I,365)+(START_MO-1)*30)						! approximate
                IF ((RULE(J,W) .EQ. 1) .or. (RULE(J,W) .EQ. 2)) THEN   ! (Options 1 and 2: rule curves)
                    IF (CURRENTYEAR<OPEYEAR(J,W)) THEN
                        FLOWOUT(J,I,W) = FLOWIN(J,I)
                        VOL(J,I+1,W) = VOL(J,I,W)
                        GOTO 123
                    END IF
                    IF (RULE(J,W) .EQ. 1) THEN
                        IF (OP1(J,1,W)>OP1(J,2,W)) THEN
                    ! Caculate target water level
                        IF ((CRTDATE .GT. OP1(J,2,W)) .and. (CRTDATE .LT. OP1(J,1,W))) THEN
                           DESIGNWL(W)=(CRTDATE-OP1(J,2,W))/(OP1(J,1,W)-OP1(J,2,W))
     &                      *(HMAX(J,W)-HMIN(J,W))
                        ELSE IF (CRTDATE .GE. OP1(J,1,W)) THEN 
                           DESIGNWL(W)=(HMAX(J,W)-HMIN(J,W))
     &                     -(CRTDATE-OP1(J,1,W))/(365-OP1(J,1,W)+OP1(J,2,W))*
     &                     (HMAX(J,W)-HMIN(J,W))
                        ELSE
                            DESIGNWL(W)=(HMAX(J,W)-HMIN(J,W))
     &                          -(CRTDATE+365-OP1(J,1,W))/(365-OP1(J,1,W)+OP1(J,2,W))*
     &                          (HMAX(J,W)-HMIN(J,W))
                        END IF
                    ELSE
                        IF ((CRTDATE .GT. OP1(J,1,W)) .and. (CRTDATE .LT. OP1(J,2,W))) THEN
                           DESIGNWL(W)=(HMAX(J,W)-HMIN(J,W))-(CRTDATE-OP1(J,1,W))/(OP1(J,2,W)-OP1(J,1,W))
     &                      *(HMAX(J,W)-HMIN(J,W))
                        ELSE IF (CRTDATE .GE. OP1(J,2,W)) THEN 
                           DESIGNWL(W)=(HMAX(J,W)-HMIN(J,W))/(365-OP1(J,2,W)+OP1(J,1,W))*(CRTDATE-OP1(J,2,W))
                        ELSE 
                            DESIGNWL(W)=(HMAX(J,W)-HMIN(J,W))/(365-OP1(J,2,W)+OP1(J,1,W))*
     &                      (CRTDATE-OP1(J,1,W))+(HMAX(J,W)-HMIN(J,W))
                        END IF

                    END IF
                    DESIGNWL(W) = DESIGNWL(W) + (HMIN(J,W) - H0(J,W))
                    ELSE
                        DESIGNWL(W) = RESLV(J,CAL_MONTH(CRTDATE),W)
                        DESIGNWL(W) = DESIGNWL(W) - H0(J,W)
                    END IF
                    HTK(J,I) = DESIGNWL(W) + H0(J,W)													! water head
                    CURRENTWL(W) = VOL(J,I,W) * (HRESERMAX(J,W)-H0(J,W))/VRESER(J,W,I)
                    IF (CURRENTWL(W)>=DESIGNWL(W)) THEN											! Zone 3
                        IF ((VOL(J,I,W)+FLOWIN(J,I)*24*3.6 -QRESER(J,W)*24*3.6)						! Case 2
     &                  >(DESIGNWL(W)* VRESER(J,W,I)) /(HRESERMAX(J,W)-H0(J,W))) THEN
                            VOL(J,I+1,W) = VOL(J,I,W) + FLOWIN(J,I)*24*3.6-QRESER(J,W)*24*3.6
                            FLOWOUT(J,I,W) = QRESER(J,W)
                            
                        ELSE																	! Case 1
                            VOL(J,I+1,W)=(DESIGNWL(W)*VRESER(J,W,I))/(HRESERMAX(J,W)-H0(J,W))
                            FLOWOUT(J,I,W)=(VOL(J,I,W)-VOL(J,I+1,W))/24/3.6 + FLOWIN(J,I)
                        END IF
                        
                        GOTO 123
                    ELSE																		! Zone 2
                        IF ((VOL(J,I,W)+FLOWIN(J,I)*24*3.6)>((DESIGNWL(W)* VRESER(J,W,I))				! Case 2
     &                  /(HRESERMAX(J,W)-H0(J,W)))) THEN
                            VOL(J,I+1,W)=(DESIGNWL(W) * VRESER(J,W,I))/(HRESERMAX(J,W)-H0(J,W))
                            FLOWOUT(J,I,W)=FLOWIN(J,I)-(VOL(J,I+1,W)-VOL(J,I,W))/24/3.6
                            IF (FLOWOUT(J,I,W)>QRESER(J,W)) THEN
                                VOL(J,I+1,W)=VOL(J,I+1,W)+(FLOWOUT(J,I,W)-QRESER(J,W))*24*3.6
                                FLOWOUT(J,I,W) = QRESER(J,W)
                            END IF
                            GOTO 123
                        ELSE																	! Case 1 (this case covers Zone 1 + Zone 2 case 1)
                            VOL(J,I+1,W) = VOL(J,I,W) + FLOWIN(J,I)*24*3.6
                            FLOWOUT(J,I,W) = 0
                            GOTO 123
                        END IF
                    END IF
                ELSE IF (RULE(J,W) .EQ. 3) THEN ! opearting rule (Option 3)
        ! Note x1 and x4 in radian (0 to pi/2), not degree
                    IF (VOL(J,I,W) < VDEAD(J,W)) THEN ! below dead water level
                        FLOWOUT(J,I,W) = 0																					! case 1
                    ELSE IF (VOL(J,I,W) .LE. X2(J,W)) THEN ! hedging
                        IF ((VOL(J,I,W)-VDEAD(J,W)+FLOWIN(J,I)*24*3.6) .LE. (Demand(J,W)+(VOL(J,I,W)-X2(J,W))*tan(X1(J,W))*24*3.6)) THEN		! discharge more than the water amount in reservoir (case 2)
                            FLOWOUT(J,I,W) = (VOL(J,I,W)-VDEAD(J,W))/24/3.6+FLOWIN(J,I)											! discharge all of water in the reservoir
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
                    VOL(J,I+1,W) = VOL(J,I,W) + (FLOWIN(J,I)-FLOWOUT(J,I,W))*24*3.6
                    GOTO 123
                ELSE IF (RULE(J,W) .EQ. 4) THEN! pre-defined time series (Option 4)
                    IF ((OP2(J,I+STARTDAY(J,W)) .GT. QRESER(J,W)) .AND. (VOL(J,I,W) .LT. VRESER(J,W,I))) THEN
                        OP2(J,I+STARTDAY(J,W)) = QRESER(J,W)
                    END IF
                    IF (OP2(J,I+STARTDAY(J,W))*24*3.6 .GT. ((VOL(J,I,W)-VDEAD(J,W)))+FLOWIN(J,I)*24*3.6) THEN
                        FLOWOUT(J,I,W) = (VOL(J,I,W)-VDEAD(J,W))/24/3.6 + FLOWIN(J,I)
                    ELSE
                        FLOWOUT(J,I,W) = OP2(J,I+STARTDAY(J,W))
                    END IF
                    VOL(J,I+1,W) = VOL(J,I,W) + (FLOWIN(J,I)-FLOWOUT(J,I,W))*24*3.6
                    GOTO 123
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
                        IF ((VOL(J,I,W)-VDEAD(J,W)+FLOWIN(J,I)*24*3.6) .LE. (Demand(J,W)+(VOL(J,I,W)-X2(J,W))*tan(X1(J,W))*24*3.6)) THEN		! discharge more than the water amount in reservoir (case 2)                    
                            FLOWOUT(J,I,W) = (VOL(J,I,W)-VDEAD(J,W))/24/3.6+FLOWIN(J,I)											! discharge all of water in the reservoir
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
                    VOL(J,I+1,W) = VOL(J,I,W) + (FLOWIN(J,I)-FLOWOUT(J,I,W))*24*3.6
                    GOTO 123
 123            END IF
        ! Check if there are any negative values
                IF (ENERGYPRO(J,I)<0) THEN
                    ENERGYPRO(J,I)=0
                END IF
                IF (FLOWOUT(J,I,W)<0) THEN
                    FLOWOUT(J,I,W)=0
                END IF
                IF (VOL(J,I,W)<0)THEN ! Not allow dropping below the minimum water level (mostly due to evaporation)
                    VOL(J,I,W)=0
                END IF
c               Remote water for irrigation
                IF (FLOWOUT(J,I,W)>=IRRIGATION(J,I)) THEN
                    FLOWOUT(J,I,W) = FLOWOUT(J,I,W) - IRRIGATION(J,I)
                ELSE
                    FLOWOUT(J,I,W) = 0
                END IF
        ! Calculate energy production
                IF (VOL(J,I+1,W)<=VRESER(J,W,I)) THEN
                    ENERGYPRO(J,I) = 0.9 * 9.81 * FLOWOUT(J,I,W)
                ELSE
                    ENERGYPRO(J,I) = 0.9 * 9.81 * QRESER(J,W)
                END IF
        !IF (J .EQ. 3) THEN
        !    WRITE(*,*) RES_DIRECT(J,1),' - ',ENERGYPRO(J,I),'  ',HHO(J,I),'  ',HHO(J,I+1),'  ',FLOWOUT(J,I),'  ',HRESERMAX(J),'  ',REALHEAD(J)
        !END IF
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
                
c               Check the neccesity to spill water
                IF (VOL(J,I+1,W)>VRESER(J,W,I)) THEN
                    FLOWOUT(J,I,W) = FLOWOUT(J,I,W)+(VOL(J,I+1,W)-VRESER(J,W,I))/24/3.6
                    VOL(J,I+1,W) = VRESER(J,W,I)
                END IF
                HHO(J,I)=VOL(J,I,W)/VRESER(J,W,I)*(HRESERMAX(J,W)-H0(J,W))+H0(J,W)
                HHO(J,I+1)=VOL(J,I+1,W)/VRESER(J,W,I)* (HRESERMAX(J,W)-H0(J,W))+H0(J,W)
                ! Note: hydraulic head calculated from the maximum water level
                ENERGYPRO(J,I) = ENERGYPRO(J,I) *
     &          (( HHO(J,I)+HHO(J,I+1))/2-(HRESERMAX(J,W)-REALHEAD(J,W)))/1000			! this part is for hydropower production estimation, ignore if work with irrigation reservoirs
                ! Propagate water to the downstream reservoir, considering the time lag
                RESORD(W) = NO_STAS(W)
                DO K = 1, NO_STAS(W)
                    IF ((RES_DIRECT(K,1,W) .EQ. RES_DIRECT(J,2,W))) THEN
                        RESORD(W) = K
                    END IF
                END DO
                IF (FLOWOUT(J,I,W) .LT. 0) THEN
                    FLOWOUT(J,I,W) = 0
                END IF
                IF (RES_DIRECT(J,2,W) .EQ. 0) THEN
                    RESFLOWS(NO_STAS(W),I+1+INT(TRVLTIME(J)),W) = RESFLOWS(NO_STAS(W),I+1+INT(TRVLTIME(J)),W) + FLOWOUT(J,I,W)
                ELSE
                    RESFLOWS(RESORD(W),I+1+INT(TRVLTIME(J)),W) = RESFLOWS(RESORD(W),I+1+INT(TRVLTIME(J)),W) + FLOWOUT(J,I,W)
                END IF
      !   WRITE(*,*) '-------------------------------------------------'
      !   WRITE(*,*) 'Hdesign ',HTK(J,I)
      !   WRITE(*,*) 'CRTDATE',CRTDATE,'Reservoir No.',J,' FLOWIN ',FLOWIN(J,I),
     &!  'FLOWOUT',FLOWOUT(J,I),' VOL ',VOL(J,I),' ELEC', ENERGYPRO(J,I)
                END DO
                IF (RESFLOWS(NO_STAS(W),I,W) .LT. 0) THEN
                     RESFLOWS(NO_STAS(W),I,W) = 0
                END IF
                FLOW(I,W) = RESFLOWS(NO_STAS(W),I,W)
            ENDDO
            ! At the end of this loop, we will have the discharge for each considered stations at time i, in FLOW(I,:)
        ENDDO
        RETURN
 9001 WRITE(*,*) 'Error reading UH data'
 9002 WRITE(*,*) 'Error in reading reservoir data'
 9003 WRITE(*,*) 'Error in reading irrigation data'
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
