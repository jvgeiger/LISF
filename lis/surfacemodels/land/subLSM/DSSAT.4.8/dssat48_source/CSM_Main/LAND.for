C=======================================================================
C COPYRIGHT 1998-2021 
C                     DSSAT Foundation
C                     University of Florida, Gainesville, Florida
C                     International Fertilizer Development Center
C                     
C ALL RIGHTS RESERVED
C=======================================================================
C=======================================================================
C  LAND UNIT Module. G.Hoogenboom, J.W.Jones, C.Porter
C-----------------------------------------------------------------------
C  Land Unit Module.  Provides the interface between soil, weather
C  and crops.  Based on the original CROPGRO routine
C=======================================================================
C  REVISION       HISTORY
C  12/01/2001 CHP Written.
C  12/12/2001 GH  Rename to Land
!  10/24/2005 CHP Put weather variables in constructed variable. 
!  02/28/2006 CHP Rename Alt_Plant to Plant, move call to CROPGRO there
!  03/03/2006 CHP Added tillage (A.Andales & WDBatchelor).
!  03/21/2006 CHP Added mulch effects
!  10/31/2007 CHP Added simple K model.

C=======================================================================
      SUBROUTINE LAND(CONTROL, ISWITCH, 
     &                YRPLT, MDATE, YREND, nest, t)
      
C-----------------------------------------------------------------------
      USE ModuleDefs    
      USE ModuleData, only : PUT !Pang: for PUT CONTROL
      USE FloodModule      
      USE CsvOutput   ! VSH 
      USE dssat48_lsmMod !Pang: Use memory from LIS
      IMPLICIT NONE
      SAVE
C------------------------------------------------
C     LIS-DSSAT
C     Pang 2023.09.29  
C------------------------------------------------
      INTEGER nest, t
C-----------------------------------------------------------------------
C     Crop, Experiment, Command line Variables
C-----------------------------------------------------------------------
      CHARACTER*2  CROP
      CHARACTER*6  ERRKEY
      PARAMETER   (ERRKEY = 'LAND  ')
      CHARACTER*8  MODEL
      CHARACTER*120 FILEIO
      
C-----------------------------------------------------------------------
C     Date / Timing / Sequencing Variables
C-----------------------------------------------------------------------
      INTEGER      DYNAMIC, YRSIM, YRDOY

C-----------------------------------------------------------------------
C     Input and Output Handling
C-----------------------------------------------------------------------
      CHARACTER*1  IDETS, IPLTI
      CHARACTER*78 MSG(2)

C-----------------------------------------------------------------------
C     Weather module Variables
C-----------------------------------------------------------------------
      TYPE (WeatherType)  WEATHER
      CHARACTER*12 FILEW
C-----------------------------------------------------------------------
C     Soil Processes Module Variables 
C-----------------------------------------------------------------------
      REAL SNOW, WINF
      REAL, DIMENSION(NL) :: NH4_plant, NO3_plant, SPi_Avail, SKi_Avail
      REAL, DIMENSION(NL) :: ST, UPPM, SW, SWDELTS, UPFLOW
      TYPE (SoilType) SOILPROP    !type defined in ModuleDefs
      TYPE (FloodWatType) FLOODWAT
      TYPE (FloodNType)   FloodN
      TYPE (MulchType)    MULCH
!     Needed for ORYZA-Rice
      REAL, DIMENSION(0:NL) :: SomLitC
      REAL, DIMENSION(0:NL,NELEM) :: SomLitE

C-----------------------------------------------------------------------
C     Soil - Plant - Atmosphere Module Variables
C-----------------------------------------------------------------------
      REAL EO, EOP, ES, SRFTEMP, TRWUP
      REAL SWDELTU(NL), SWDELTX(NL) !, RWU(NL)
!     Needed for CaneGro_SA
      REAL EOS, EP, TRWU
!     Calculated by ORYZA-Rice
      REAL UH2O(NL)
!     Needed for SALUS
      REAL RWU(NL)

C-----------------------------------------------------------------------
C     PLANT Module Variables
C-----------------------------------------------------------------------
      INTEGER MDATE
      INTEGER STGDOY(20)
      REAL CANHT, EORATIO, NSTRES, PORMIN, PSTRES1, RWUMX
      REAL XHLAI, XLAI
      REAL KSEVAP, KTRANS
      REAL, Dimension(NL) :: PUptake, RLV, FracRts, UNH4, UNO3, KUptake
      Type (ResidueType) HARVRES  !type defined in ModuleDefs
      Type (ResidueType) SENESCE  
      
C-----------------------------------------------------------------------
C     Operations Management Module Variables 
C-----------------------------------------------------------------------
      TYPE (TillType) TILLVALS
      INTEGER YREND, YRPLT
      REAL IRRAMT
      REAL, DIMENSION(2) :: HARVFRAC   !Harvest & byproduct fractions
      TYPE (FertType) FERTDATA         !Fertilizer application
      TYPE (OrgMatAppType)OMAData      !Organic matter application

C-----------------------------------------------------------------------
!!     Temporary timer function
!!     Date / time variables
!      INTEGER DATE_TIME(8)
!!      date_time(1)  The 4-digit year  
!!      date_time(2)  The month of the year  
!!      date_time(3)  The day of the month  
!!      date_time(4)  The time difference with respect to Coordinated Universal Time (UTC) in minutes  
!!      date_time(5)  The hour of the day (range 0 to 23) - local time  
!!      date_time(6)  The minutes of the hour (range 0 to 59) - local time  
!!      date_time(7)  The seconds of the minute (range 0 to 59) - local time  
!!      date_time(8)  The milliseconds of the second (range 0 to 999) - local time  
!      REAL TIME0, TIME1, TIME_START DELTA_TIME
C-----------------------------------------------------------------------

C     Define constructed variable types based on definitions in
C     ModuleDefs.for.
      TYPE (ControlType) CONTROL
      TYPE (SwitchType)  ISWITCH

      CALL PUT(CONTROL) !PL: Put the current CRONTROL to the memory for GET in modules
      CALL PUT(ISWITCH)
C     Transfer values from constructed data types into local variables.
      CROP    = CONTROL % CROP
      DYNAMIC = CONTROL % DYNAMIC
      FILEIO  = CONTROL % FILEIO
      MODEL   = CONTROL % MODEL
      YRDOY   = CONTROL % YRDOY
      YRSIM   = CONTROL % YRSIM

      IPLTI   = ISWITCH % IPLTI

      !---- Pang: 2023.09.26 -------------------------------------------
      !---- SOIL -------------------------------------------------------
      IDETS   = ISWITCH % IDETS !Used for Daily Loop

      SOILPROP = dssat48_struc(nest)%SOILPROP(t) !Pang: Obtain SOILPROP from memory
      MULCH = dssat48_struc(nest)%MULCH(t) !Pang: Obtain MULCH from memory
      FLOODWAT = dssat48_struc(nest)%FLOODWAT(t) !Pang: Obtain FLOODWAT from memory
      FLOODN = dssat48_struc(nest)%FLOODN(t) !Pang: Obtain FLOODN from memory
      SW = dssat48_struc(nest)%dssat48(t)%SW
      SNOW = dssat48_struc(nest)%dssat48(t)%SNOW
      SomLitC = dssat48_struc(nest)%dssat48(t)%SomLitC
      SKi_Avail = dssat48_struc(nest)%dssat48(t)%SKi_Avail
      !----- SOIL ADDED ------------------------------------------------
      SomLitE = dssat48_struc(nest)%dssat48(t)%SomLitE
      NH4_plant = dssat48_struc(nest)%dssat48(t)%NH4_plant
      NO3_Plant = dssat48_struc(nest)%dssat48(t)%NO3_Plant
      UPPM = dssat48_struc(nest)%dssat48(t)%UPPM
      SPi_AVAIL = dssat48_struc(nest)%dssat48(t)%SPi_AVAIL
      SWDELTS = dssat48_struc(nest)%dssat48(t)%SWDELTS
      WINF = dssat48_struc(nest)%dssat48(t)%WINF
      !PRINT*, 'SOILPROP bg LAND: ', SOILPROP
      !PRINT*, 'FLOODWAT bg LAND: ', FLOODWAT
      !PRINT*, 'FLOODN bg LAND: ', FLOODN
      !PRINT*, 'MULCH bg LAND: ', MULCH
      !PRINT*, 'SW bg LAND: ', SW
      !PRINT*, 'SNOW bg LAND: ', SNOW
      !PRINT*, 'SomLitC bg LAND: ', SomLitC
      !PRINT*, 'Ski_Avail bg LAND: ', SKi_Avail
      !---- SPAM -------------------------------------------------------
      !---- Pang: 2024.01.12 -------------------------------------------
      ES = dssat48_struc(nest)%dssat48(t)%ES
      ST = dssat48_struc(nest)%dssat48(t)%ST
      SWDELTX = dssat48_struc(nest)%dssat48(t)%SWDELTX
      UPFLOW = dssat48_struc(nest)%dssat48(t)%UPFLOW
      !---- Pang: 2024.02.05 -------------------------------------------
      EO = dssat48_struc(nest)%dssat48(t)%EO
      EOP = dssat48_struc(nest)%dssat48(t)%EOP
      EOS = dssat48_struc(nest)%dssat48(t)%EOS
      EP = dssat48_struc(nest)%dssat48(t)%EP
      SRFTEMP = dssat48_struc(nest)%dssat48(t)%SRFTEMP
      TRWU = dssat48_struc(nest)%dssat48(t)%TRWU
      TRWUP = dssat48_struc(nest)%dssat48(t)%TRWUP
      RWU = dssat48_struc(nest)%dssat48(t)%RWU
      SWDELTU = dssat48_struc(nest)%dssat48(t)%SWDELTU
      !---- Pang: 2024.01.17 -------------------------------------------
      !---- PLANT ------------------------------------------------------
      HARVRES = dssat48_struc(nest)%HARVRES(t)
      SENESCE = dssat48_struc(nest)%SENESCE(t)

      CANHT = dssat48_struc(nest)%dssat48(t)%CANHT
      EORATIO = dssat48_struc(nest)%dssat48(t)%EORATIO
      NSTRES = dssat48_struc(nest)%dssat48(t)%NSTRES
      PSTRES1 = dssat48_struc(nest)%dssat48(t)%PSTRES1
      PORMIN = dssat48_struc(nest)%dssat48(t)%PORMIN
      RWUMX = dssat48_struc(nest)%dssat48(t)%RWUMX
      KSEVAP = dssat48_struc(nest)%dssat48(t)%KSEVAP
      KTRANS = dssat48_struc(nest)%dssat48(t)%KTRANS
      KUptake = dssat48_struc(nest)%dssat48(t)%KUptake
      PUptake = dssat48_struc(nest)%dssat48(t)%PUptake
      RLV = dssat48_struc(nest)%dssat48(t)%RLV
      FracRts = dssat48_struc(nest)%dssat48(t)%FracRts
      UNH4 = dssat48_struc(nest)%dssat48(t)%UNH4
      UNO3 = dssat48_struc(nest)%dssat48(t)%UNO3
      STGDOY = dssat48_struc(nest)%dssat48(t)%STGDOY
      XHLAI = dssat48_struc(nest)%dssat48(t)%XHLAI
      XLAI = dssat48_struc(nest)%dssat48(t)%XLAI
      UH2O = dssat48_struc(nest)%dssat48(t)%UH2O

      !---- Pang: 2024.01.17 -------------------------------------------
      !---- MGMTOPS ------------------------------------------------------

      FERTDATA = dssat48_struc(nest)%FERTDATA(t)
      OMADATA = dssat48_struc(nest)%OMADATA(t)
      TILLVALS = dssat48_struc(nest)%TILLVALS(t)
      HARVFRAC = dssat48_struc(nest)%dssat48(t)%HARVFRAC
      IRRAMT = dssat48_struc(nest)%dssat48(t)%IRRAMT

C***********************************************************************
C***********************************************************************
C     Run Initialization - Called once per simulation
C***********************************************************************
      IF (DYNAMIC .EQ. RUNINIT) THEN
C-----------------------------------------------------------------------
!!     Temporary timer function
!      !Get initial time
!      CALL DATE_AND_TIME (VALUES=DATE_TIME)
!!     Convert time to seconds
!      TIME0 = DATE_TIME(7) 
!     &      + DATE_TIME(8) / 1000.  
!     &      + DATE_TIME(6) * 60.  
!     &      + DATE_TIME(5) * 3600.
!      TIME_START = TIME0
C-----------------------------------------------------------------------
C     Read switches from FILEIO
C-----------------------------------------------------------------------
      !Pang: Read INP file. LUNIO & N_ELEMS are Updated After IPIBS
      !      while ISWITCH is same before and after
      CALL IPIBS (CONTROL, ISWITCH, CROP, IDETS, MODEL)
C-----------------------------------------------------------------------
C     Read input parameters for weather routines
C-----------------------------------------------------------------------
      CALL WEATHR(CONTROL, ISWITCH, nest, t, WEATHER, YREND)  !JE
C-----------------------------------------------------------------------
C     Read initial soil data 
C-----------------------------------------------------------------------
      CALL SOIL(CONTROL, ISWITCH, nest, t,               !Pang2023.09.19 
     &    ES, FERTDATA, FracRts, HARVRES, IRRAMT,         !Input
     &    KTRANS, KUptake, OMAData, PUptake, RLV,         !Input
     &    SENESCE, ST, SWDELTX,TILLVALS, UNH4, UNO3,      !Input
     &    WEATHER, XHLAI,                                 !Input
     &    FLOODN, FLOODWAT, MULCH, UPFLOW,                !I/O
     &    NH4_plant, NO3_plant, SKi_AVAIL, SNOW,          !Output
     &    SPi_AVAIL, SOILPROP, SomLitC, SomLitE,          !Output
     &    SW, SWDELTS, SWDELTU, UPPM, WINF, YREND)        !Output
C-----------------------------------------------------------------------
C     Read initial soil-plant-atmosphere data
C-----------------------------------------------------------------------
      CALL SPAM(CONTROL, ISWITCH,nest, t,
     &    CANHT, EORATIO, KSEVAP, KTRANS, MULCH,          !Input
     &    PSTRES1, PORMIN, RLV, RWUMX, SOILPROP, SW,      !Input
     &    SWDELTS, UH2O, WEATHER, WINF, XHLAI, XLAI,      !Input
     &    FLOODWAT, SWDELTU,                              !I/O
     &    EO, EOP, EOS, EP, ES, RWU, SRFTEMP, ST,         !Output
     &    SWDELTX, TRWU, TRWUP, UPFLOW)                   !Output

C-----------------------------------------------------------------------
C     Read initial plant module data
C-----------------------------------------------------------------------
      CALL PLANT(CONTROL, ISWITCH, nest, t,              !Pang2024.01.18
     &    EO, EOP, EOS, EP, ES, FLOODWAT, HARVFRAC,       !Input
     &    NH4_plant, NO3_plant, SKi_Avail, SomLitC, SomLitE, !Input
     &    SPi_AVAIL, SNOW, SOILPROP, SRFTEMP, ST, SW,     !Input
     &    TRWU, TRWUP, UPPM, WEATHER, YREND, YRPLT,       !Input
     &    IRRAMT,                                         !Input
     &    FLOODN,                                         !I/O
     &    CANHT, EORATIO, HARVRES, KSEVAP, KTRANS,        !Output
     &    KUptake, MDATE, NSTRES, PSTRES1,                !Output
     &    PUptake, PORMIN, RLV, RWUMX, SENESCE,           !Output
     &    STGDOY, FracRts, UH2O, UNH4, UNO3, XHLAI, XLAI) !Output

C-----------------------------------------------------------------------
C     Initialize summary.out information
C-----------------------------------------------------------------------
      !CALL OPSUM (CONTROL, ISWITCH, YRPLT) !Pang 2024.05.17: turn off to reduce I/O

C*********************************************************************** 
C*********************************************************************** 
C     SEASONAL INITIALIZATION
C*********************************************************************** 
      ELSEIF (DYNAMIC .EQ. SEASINIT) THEN
C-----------------------------------------------------------------------
C     Call WEATHR for initialization - reads first day of weather
C     data for use in soil N and soil temp initialization.
C-----------------------------------------------------------------------
      CALL WEATHR(CONTROL, ISWITCH, nest, t, WEATHER, YREND)

C-----------------------------------------------------------------------
C     Set planting date, adjust operations dates for seasonal or 
C     sequenced runs.
C-----------------------------------------------------------------------
      CALL MGMTOPS(CONTROL, ISWITCH, nest, t, !Pang 2024.03.11
     &    FLOODWAT, HARVRES, SOILPROP, ST,                !Input 
     &    STGDOY, SW, WEATHER,                            !Input
     &    YREND, FERTDATA, HARVFRAC, IRRAMT,              !Output
     &    MDATE, OMADATA, TILLVALS, YRPLT)                !Output

C-----------------------------------------------------------------------
      IF (YRPLT < YRSIM .AND. CROP /= 'FA' .AND.
     &    INDEX('AF', IPLTI) == 0) THEN
          CALL ERROR(ERRKEY,2,' ',0)
      ENDIF

C-----------------------------------------------------------------------
C     Seasonal initialization for soil processes
C-----------------------------------------------------------------------
      !SOILPROP = dssat48_struc(nest)%SOILPROP(t) !Pang: Obtain SOILPROP from mem
      CALL SOIL(CONTROL, ISWITCH, nest, t,                 !Pang2023.09.19
     &    ES, FERTDATA, FracRts, HARVRES, IRRAMT,         !Input
     &    KTRANS, KUptake, OMAData, PUptake, RLV,         !Input
     &    SENESCE, ST, SWDELTX,TILLVALS, UNH4, UNO3,      !Input
     &    WEATHER, XHLAI,                                 !Input
     &    FLOODN, FLOODWAT, MULCH, UPFLOW,                !I/O
     &    NH4_plant, NO3_plant, SKi_AVAIL, SNOW,          !Output
     &    SPi_AVAIL, SOILPROP, SomLitC, SomLitE,          !Output
     &    SW, SWDELTS, SWDELTU, UPPM, WINF, YREND)        !Output
C-----------------------------------------------------------------------
C     Seasonal initialization for soil-plant-atmosphere processes
!     chp moved this before SOIL, so soil temp is available 
!     update 2020-12-04 - order makes no difference
C-----------------------------------------------------------------------
      CALL SPAM(CONTROL, ISWITCH,nest, t,
     &    CANHT, EORATIO, KSEVAP, KTRANS, MULCH,          !Input
     &    PSTRES1, PORMIN, RLV, RWUMX, SOILPROP, SW,      !Input
     &    SWDELTS, UH2O, WEATHER, WINF, XHLAI, XLAI,      !Input
     &    FLOODWAT, SWDELTU,                              !I/O
     &    EO, EOP, EOS, EP, ES, RWU, SRFTEMP, ST,         !Output
     &    SWDELTX, TRWU, TRWUP, UPFLOW)                   !Output
!C-----------------------------------------------------------------------
!C     Seasonal initialization for soil processes
!C-----------------------------------------------------------------------
!      CALL SOIL(CONTROL, ISWITCH, 
!     &    ES, FERTDATA, HARVRES, IRRAMT, KTRANS,          !Input
!     &    KUptake, OMAData, PUptake, SENESCE, ST,         !Input
!     &    FracRts, SWDELTX,TILLVALS, UNH4, UNO3, UPFLOW,  !Input
!     &    WEATHER, XHLAI, FLOODN, FLOODWAT, MULCH,        !I/O
!     &    NH4_plant, NO3_plant, SKi_AVAIL, SNOW,          !Output
!     &    SPi_AVAIL, SOILPROP, SomLitC, SomLitE,          !Output
!     &    SW, SWDELTS, SWDELTU, UPPM, WINF, YREND)        !Output

C-----------------------------------------------------------------------
C     Initialize PLANT routines (including phenology and pest)
C-----------------------------------------------------------------------
      CALL PLANT(CONTROL, ISWITCH, nest, t,              !Pang2024.01.18 
     &    EO, EOP, EOS, EP, ES, FLOODWAT, HARVFRAC,       !Input
     &    NH4_plant, NO3_plant, SKi_Avail, SomLitC, SomLitE, !Input
     &    SPi_AVAIL, SNOW, SOILPROP, SRFTEMP, ST, SW,     !Input
     &    TRWU, TRWUP, UPPM, WEATHER, YREND, YRPLT,       !Input
     &    IRRAMT,                                         !Input
     &    FLOODN,                                         !I/O
     &    CANHT, EORATIO, HARVRES, KSEVAP, KTRANS,        !Output
     &    KUptake, MDATE, NSTRES, PSTRES1,                !Output
     &    PUptake, PORMIN, RLV, RWUMX, SENESCE,           !Output
     &    STGDOY, FracRts, UH2O, UNH4, UNO3, XHLAI, XLAI) !Output
C-----------------------------------------------------------------------
C     Initialize summary output file - possible output from 
C     various modules.
C-----------------------------------------------------------------------
      IF (IDETS .EQ. 'Y' .OR. IDETS .EQ. 'A') THEN
       ! CALL OPSUM (CONTROL, ISWITCH, YRPLT) !Pang 2024.05.17: turn off to reduce I/O
      ENDIF
C***********************************************************************
C***********************************************************************
C     DAILY RATE CALCULATIONS
C***********************************************************************
      ELSE IF (DYNAMIC .EQ. RATE) THEN
C-----------------------------------------------------------------------
C     Call WEATHER Subroutine to input weather data and to
C     calculate hourly radiation and air temperature values
C     Note: First day of weather has already been read by 
C       initialization call to WEATHR.
C-----------------------------------------------------------------------
      CALL WEATHR(CONTROL, ISWITCH, nest, t, WEATHER, YREND)

C-----------------------------------------------------------------------
C     Call Operations Management module to determine today's 
C     applications of irrigation, tillage, etc.
C-----------------------------------------------------------------------
      CALL MGMTOPS(CONTROL, ISWITCH,nest, t, !Pang 2024.03.11 
     &    FLOODWAT, HARVRES, SOILPROP, ST,                !Input 
     &    STGDOY, SW, WEATHER,                            !Input
     &    YREND, FERTDATA, HARVFRAC, IRRAMT,              !Output
     &    MDATE, OMADATA, TILLVALS, YRPLT)                !Output

C-----------------------------------------------------------------------
C     Call Soil processes module to determine today's rates of 
C     change of soil properties.
C-----------------------------------------------------------------------
      CALL SOIL(CONTROL, ISWITCH, nest, t,                !Pang2023.09.19
     &    ES, FERTDATA, FracRts, HARVRES, IRRAMT,         !Input
     &    KTRANS, KUptake, OMAData, PUptake, RLV,         !Input
     &    SENESCE, ST, SWDELTX,TILLVALS, UNH4, UNO3,      !Input
     &    WEATHER, XHLAI,                                 !Input
     &    FLOODN, FLOODWAT, MULCH, UPFLOW,                !I/O
     &    NH4_plant, NO3_plant, SKi_AVAIL, SNOW,          !Output
     &    SPi_AVAIL, SOILPROP, SomLitC, SomLitE,          !Output
     &    SW, SWDELTS, SWDELTU, UPPM, WINF, YREND)        !Output

C-----------------------------------------------------------------------
C     Call Soil-plant-atmosphere module to determine today's
C     rates of evapotranspiration.
C-----------------------------------------------------------------------
      CALL SPAM(CONTROL, ISWITCH,nest, t,
     &    CANHT, EORATIO, KSEVAP, KTRANS, MULCH,          !Input
     &    PSTRES1, PORMIN, RLV, RWUMX, SOILPROP, SW,      !Input
     &    SWDELTS, UH2O, WEATHER, WINF, XHLAI, XLAI,      !Input
     &    FLOODWAT, SWDELTU,                              !I/O
     &    EO, EOP, EOS, EP, ES, RWU, SRFTEMP, ST,         !Output
     &    SWDELTX, TRWU, TRWUP, UPFLOW)                   !Output

C-----------------------------------------------------------------------
C     Call PLANT Subroutine to calculate crop growth and
C     development rates.
C     Skip plant growth and development routines for fallow runs
C-----------------------------------------------------------------------
!      IF (CROP .NE. 'FA' .AND. 
!     &    YRDOY .GE. YRPLT .AND. YRPLT .NE. -99) THEN
        CALL PLANT(CONTROL, ISWITCH,  nest, t,              !Pang2024.01.18
     &    EO, EOP, EOS, EP, ES, FLOODWAT, HARVFRAC,       !Input
     &    NH4_plant, NO3_plant, SKi_Avail, SomLitC, SomLitE, !Input
     &    SPi_AVAIL, SNOW, SOILPROP, SRFTEMP, ST, SW,     !Input
     &    TRWU, TRWUP, UPPM, WEATHER, YREND, YRPLT,       !Input
     &    IRRAMT,                                         !Input
     &    FLOODN,                                         !I/O
     &    CANHT, EORATIO, HARVRES, KSEVAP, KTRANS,        !Output
     &    KUptake, MDATE, NSTRES, PSTRES1,                !Output
     &    PUptake, PORMIN, RLV, RWUMX, SENESCE,           !Output
     &    STGDOY, FracRts, UH2O, UNH4, UNO3, XHLAI, XLAI) !Output
!      ENDIF

C***********************************************************************
C     DAILY INTEGRATION 
C***********************************************************************
      ELSE IF (DYNAMIC .EQ. INTEGR) THEN
C***********************************************************************
C     Integrate soil state variables
C-----------------------------------------------------------------------
      CALL SOIL(CONTROL, ISWITCH, nest, t,                !Pang2023.09.19
     &    ES, FERTDATA, FracRts, HARVRES, IRRAMT,         !Input
     &    KTRANS, KUptake, OMAData, PUptake, RLV,         !Input
     &    SENESCE, ST, SWDELTX,TILLVALS, UNH4, UNO3,      !Input
     &    WEATHER, XHLAI,                                 !Input
     &    FLOODN, FLOODWAT, MULCH, UPFLOW,                !I/O
     &    NH4_plant, NO3_plant, SKi_AVAIL, SNOW,          !Output
     &    SPi_AVAIL, SOILPROP, SomLitC, SomLitE,          !Output
     &    SW, SWDELTS, SWDELTU, UPPM, WINF, YREND)        !Output

C-----------------------------------------------------------------------
C     Compute cumulative totals for soil-plant-atmosphere processes
C-----------------------------------------------------------------------
      CALL SPAM(CONTROL, ISWITCH,nest, t,
     &    CANHT, EORATIO, KSEVAP, KTRANS, MULCH,          !Input
     &    PSTRES1, PORMIN, RLV, RWUMX, SOILPROP, SW,      !Input
     &    SWDELTS, UH2O, WEATHER, WINF, XHLAI, XLAI,      !Input
     &    FLOODWAT, SWDELTU,                              !I/O
     &    EO, EOP, EOS, EP, ES, RWU, SRFTEMP, ST,         !Output
     &    SWDELTX, TRWU, TRWUP, UPFLOW)                   !Output

C-----------------------------------------------------------------------
C     Call Plant module to integrate daily plant processes and update
C     plant state variables.
C-----------------------------------------------------------------------
      IF (CROP .NE. 'FA' .AND. 
     &        YRDOY .GE. YRPLT .AND. YRPLT .NE. -99) THEN
        CALL PLANT(CONTROL, ISWITCH,  nest, t,              !Pang2024.01.18
     &    EO, EOP, EOS, EP, ES, FLOODWAT, HARVFRAC,       !Input
     &    NH4_plant, NO3_plant, SKi_Avail, SomLitC, SomLitE, !Input
     &    SPi_AVAIL, SNOW, SOILPROP, SRFTEMP, ST, SW,     !Input
     &    TRWU, TRWUP, UPPM, WEATHER, YREND, YRPLT,       !Input
     &    IRRAMT,                                         !Input
     &    FLOODN,                                         !I/O
     &    CANHT, EORATIO, HARVRES, KSEVAP, KTRANS,        !Output
     &    KUptake, MDATE, NSTRES, PSTRES1,                !Output
     &    PUptake, PORMIN, RLV, RWUMX, SENESCE,           !Output
     &    STGDOY, FracRts, UH2O, UNH4, UNO3, XHLAI, XLAI) !Output
      ENDIF

C-----------------------------------------------------------------------
C     Call Operations Management module to check for harvest end, 
C     accumulate variables.
C-----------------------------------------------------------------------
      CALL MGMTOPS(CONTROL, ISWITCH, nest, t, !Pang 2024.03.11
     &    FLOODWAT, HARVRES, SOILPROP, ST,                !Input 
     &    STGDOY, SW, WEATHER,                            !Input
     &    YREND, FERTDATA, HARVFRAC, IRRAMT,              !Output
     &    MDATE, OMADATA, TILLVALS, YRPLT)                !Output

C***********************************************************************
C***********************************************************************
C     Daily Output
C***********************************************************************
      ELSE IF (DYNAMIC .EQ. OUTPUT) THEN

      CALL WEATHR(CONTROL, ISWITCH, nest, t, WEATHER, YREND)

        CALL SOIL(CONTROL, ISWITCH, nest, t,              !Pang2023.09.19
     &    ES, FERTDATA, FracRts, HARVRES, IRRAMT,         !Input
     &    KTRANS, KUptake, OMAData, PUptake, RLV,         !Input
     &    SENESCE, ST, SWDELTX,TILLVALS, UNH4, UNO3,      !Input
     &    WEATHER, XHLAI,                                 !Input
     &    FLOODN, FLOODWAT, MULCH, UPFLOW,                !I/O
     &    NH4_plant, NO3_plant, SKi_AVAIL, SNOW,          !Output
     &    SPi_AVAIL, SOILPROP, SomLitC, SomLitE,          !Output
     &    SW, SWDELTS, SWDELTU, UPPM, WINF, YREND)        !Output

        CALL SPAM(CONTROL, ISWITCH,nest, t,
     &    CANHT, EORATIO, KSEVAP, KTRANS, MULCH,          !Input
     &    PSTRES1, PORMIN, RLV, RWUMX, SOILPROP, SW,      !Input
     &    SWDELTS, UH2O, WEATHER, WINF, XHLAI, XLAI,      !Input
     &    FLOODWAT, SWDELTU,                              !I/O
     &    EO, EOP, EOS, EP, ES, RWU, SRFTEMP, ST,         !Output
     &    SWDELTX, TRWU, TRWUP, UPFLOW)                   !Output

C-----------------------------------------------------------------------
C     Call plant module for daily printout.
C-----------------------------------------------------------------------
        IF (CROP .NE. 'FA') THEN
          CALL PLANT(CONTROL, ISWITCH, nest, t,           !Pang2024.01.18
     &    EO, EOP, EOS, EP, ES, FLOODWAT, HARVFRAC,       !Input
     &    NH4_plant, NO3_plant, SKi_Avail, SomLitC, SomLitE, !Input
     &    SPi_AVAIL, SNOW, SOILPROP, SRFTEMP, ST, SW,     !Input
     &    TRWU, TRWUP, UPPM, WEATHER, YREND, YRPLT,       !Input
     &    IRRAMT,                                         !Input
     &    FLOODN,                                         !I/O
     &    CANHT, EORATIO, HARVRES, KSEVAP, KTRANS,        !Output
     &    KUptake, MDATE, NSTRES, PSTRES1,                !Output
     &    PUptake, PORMIN, RLV, RWUMX, SENESCE,           !Output
     &    STGDOY, FracRts, UH2O, UNH4, UNO3, XHLAI, XLAI) !Output
        ENDIF

        CALL MGMTOPS(CONTROL, ISWITCH, nest, t,  !Pang 2024.03.11 
     &    FLOODWAT, HARVRES, SOILPROP, ST,                !Input 
     &    STGDOY, SW, WEATHER,                            !Input
     &    YREND, FERTDATA, HARVFRAC, IRRAMT,              !Output
     &    MDATE, OMADATA, TILLVALS, YRPLT)                !Output

C*********************************************************************** 
C***********************************************************************
C     Seasonal Output
C*********************************************************************** 
      ELSE IF (DYNAMIC .EQ. SEASEND) THEN

C     Call WEATHER module to close current weather file 
      CALL WEATHR(CONTROL, ISWITCH, nest, t, WEATHER, YREND)

C     Print seasonal summaries and close files.
      CALL SOIL(CONTROL, ISWITCH, nest, t,               !Pang2023.09.19
     &    ES, FERTDATA, FracRts, HARVRES, IRRAMT,         !Input
     &    KTRANS, KUptake, OMAData, PUptake, RLV,         !Input
     &    SENESCE, ST, SWDELTX,TILLVALS, UNH4, UNO3,      !Input
     &    WEATHER, XHLAI,                                 !Input
     &    FLOODN, FLOODWAT, MULCH, UPFLOW,                !I/O
     &    NH4_plant, NO3_plant, SKi_AVAIL, SNOW,          !Output
     &    SPi_AVAIL, SOILPROP, SomLitC, SomLitE,          !Output
     &    SW, SWDELTS, SWDELTU, UPPM, WINF, YREND)        !Output

      CALL SPAM(CONTROL, ISWITCH,nest,t,
     &    CANHT, EORATIO, KSEVAP, KTRANS, MULCH,          !Input
     &    PSTRES1, PORMIN, RLV, RWUMX, SOILPROP, SW,      !Input
     &    SWDELTS, UH2O, WEATHER, WINF, XHLAI, XLAI,      !Input
     &    FLOODWAT, SWDELTU,                              !I/O
     &    EO, EOP, EOS, EP, ES, RWU, SRFTEMP, ST,         !Output
     &    SWDELTX, TRWU, TRWUP, UPFLOW)                   !Output

      CALL PLANT(CONTROL, ISWITCH, nest, t,              !Pang2024.01.18 
     &    EO, EOP, EOS, EP, ES, FLOODWAT, HARVFRAC,       !Input
     &    NH4_plant, NO3_plant, SKi_Avail, SomLitC, SomLitE, !Input
     &    SPi_AVAIL, SNOW, SOILPROP, SRFTEMP, ST, SW,     !Input
     &    TRWU, TRWUP, UPPM, WEATHER, YREND, YRPLT,       !Input
     &    IRRAMT,                                         !Input
     &    FLOODN,                                         !I/O
     &    CANHT, EORATIO, HARVRES, KSEVAP, KTRANS,        !Output
     &    KUptake, MDATE, NSTRES, PSTRES1,                !Output
     &    PUptake, PORMIN, RLV, RWUMX, SENESCE,           !Output
     &    STGDOY, FracRts, UH2O, UNH4, UNO3, XHLAI, XLAI) !Output

!     Call management operations module for seasonal printout.
      CALL MGMTOPS(CONTROL, ISWITCH, nest, t,  !Pang 2024.03.11 
     &    FLOODWAT, HARVRES, SOILPROP, ST,                !Input 
     &    STGDOY, SW, WEATHER,                            !Input
     &    YREND, FERTDATA, HARVFRAC, IRRAMT,              !Output
     &    MDATE, OMADATA, TILLVALS, YRPLT)                !Output

C-----------------------------------------------------------------------
C     Seasonal Output
C     Call end of season and summary output subroutines
C-----------------------------------------------------------------------
      !CALL OPSUM (CONTROL, ISWITCH, YRPLT) !Pang 2024.05.17: turn off to reduce I/O

!!     Temporary timer function
!      CALL DATE_AND_TIME (VALUES=DATE_TIME)
!      
!!     Convert time to seconds
!      TIME1 = DATE_TIME(7) 
!     &      + DATE_TIME(8) / 1000.  
!     &      + DATE_TIME(6) * 60.  
!     &      + DATE_TIME(5) * 3600.
!      DELTA_TIME = TIME1 - TIME0
!      WRITE(200,'(1X,"RUN ",I3,3X,F10.3)') RUN, DELTA_TIME
!      TIME0 = TIME1

      IF (CONTROL % ERRCODE > 0) THEN
        WRITE(MSG(1),'(A,I8)') "End of run ", CONTROL % RUN
        WRITE(MSG(2),'("Simulation ended with error code ",I3)') 
     &      CONTROL % ERRCODE
        CALL WARNING(2,'ENDRUN',MSG)
        CALL INFO(2,'ENDRUN',MSG)
      ELSE
        WRITE(MSG(1),'(A,I8)') "Normal end of run ", CONTROL % RUN
        CALL WARNING(0,'ENDRUN',MSG)
        CALL INFO(1,'ENDRUN',MSG)
      ENDIF
      
!     VSH
      if (SOILPROP % NLAYR > maxnlayers ) then
         maxnlayers = SOILPROP % NLAYR
      end if 
C*********************************************************************** 
C***********************************************************************
C     End of Run
C*********************************************************************** 
      ELSE IF (DYNAMIC .EQ. ENDRUN) THEN
        CALL SOIL(CONTROL, ISWITCH, nest, t,             !Pang2023.09.19
     &    ES, FERTDATA, FracRts, HARVRES, IRRAMT,         !Input
     &    KTRANS, KUptake, OMAData, PUptake, RLV,         !Input
     &    SENESCE, ST, SWDELTX,TILLVALS, UNH4, UNO3,      !Input
     &    WEATHER, XHLAI,                                 !Input
     &    FLOODN, FLOODWAT, MULCH, UPFLOW,                !I/O
     &    NH4_plant, NO3_plant, SKi_AVAIL, SNOW,          !Output
     &    SPi_AVAIL, SOILPROP, SomLitC, SomLitE,          !Output
     &    SW, SWDELTS, SWDELTU, UPPM, WINF, YREND)        !Output

!!     Temporary timer function
!      CALL DATE_AND_TIME (VALUES=DATE_TIME)
!      
!!     Convert time to seconds
!      TIME1 = DATE_TIME(7) 
!     &      + DATE_TIME(8) / 1000.  
!     &      + DATE_TIME(6) * 60.  
!     &      + DATE_TIME(5) * 3600.
!      DELTA_TIME = TIME1 - TIME_START
!      WRITE(200,'(/," Total Time",F10.3)') RUN, DELTA_TIME

!      VSH CSV outputs
       IF (ISWITCH % FMOPT == 'C') THEN
          CALL CsvOutputs(CONTROL % MODEL(1:5), CONTROL % N_ELEMS,
     & maxnlayers)
        END IF 
!***********************************************************************
!***********************************************************************
!     END OF DYNAMIC IF CONSTRUCT
!***********************************************************************
      ENDIF
!------ Pang: 2023.09.27 --------------------------------------------------
!------ SOIL --------------------------------------------------------------
      dssat48_struc(nest)%SOILPROP(t) = SOILPROP !Assign updated var to mem
      dssat48_struc(nest)%MULCH(t) = MULCH
      dssat48_struc(nest)%FLOODWAT(t) = FLOODWAT
      dssat48_struc(nest)%FLOODN(t) = FLOODN
      dssat48_struc(nest)%dssat48(t)%SW = SW
      dssat48_struc(nest)%dssat48(t)%SNOW = SNOW
      dssat48_struc(nest)%dssat48(t)%SomLitC = SomLitC
      dssat48_struc(nest)%dssat48(t)%SKi_Avail = SKi_Avail
!----- SOIL ADDED ------------------------------------------------
      dssat48_struc(nest)%dssat48(t)%SomLitE = SomLitE
      dssat48_struc(nest)%dssat48(t)%NH4_plant = NH4_plant
      dssat48_struc(nest)%dssat48(t)%NO3_Plant = NO3_Plant
      dssat48_struc(nest)%dssat48(t)%UPPM = UPPM
      dssat48_struc(nest)%dssat48(t)%SPi_AVAIL = SPi_AVAIL
      dssat48_struc(nest)%dssat48(t)%SWDELTS = SWDELTS
      dssat48_struc(nest)%dssat48(t)%WINF = WINF

      !PRINT*, 'SOILPROP af LAND: ', SOILPROP
      !PRINT*, 'FLOODWAT af LAND: ', FLOODWAT
      !PRINT*, 'FLOODN af LAND: ', FLOODN
      !PRINT*, 'MULCH af LAND: ', MULCH
      !PRINT*, 'SW af LAND: ', SW
      !PRINT*, 'SNOW af LAND: ', SNOW
      !PRINT*, 'SomLitC af LAND: ', SomLitC
      !PRINT*, 'Ski_Avail af LAND: ', SKi_Avail
!------ SPAM --------------------------------------------------------------
!------ Pang: 2024.01.12 --------------------------------------------------
      dssat48_struc(nest)%dssat48(t)%ES = ES
      dssat48_struc(nest)%dssat48(t)%ST = ST
      dssat48_struc(nest)%dssat48(t)%SWDELTX = SWDELTX
      dssat48_struc(nest)%dssat48(t)%UPFLOW = UPFLOW
!---- Pang: 2024.02.05 -------------------------------------------
      dssat48_struc(nest)%dssat48(t)%EO = EO
      dssat48_struc(nest)%dssat48(t)%EOP = EOP
      dssat48_struc(nest)%dssat48(t)%EOS = EOS
      dssat48_struc(nest)%dssat48(t)%EP = EP
      dssat48_struc(nest)%dssat48(t)%SRFTEMP =SRFTEMP
      dssat48_struc(nest)%dssat48(t)%TRWU = TRWU
      dssat48_struc(nest)%dssat48(t)%TRWUP = TRWUP
      dssat48_struc(nest)%dssat48(t)%RWU = RWU
      dssat48_struc(nest)%dssat48(t)%SWDELTU = SWDELTU

!------ Pang: 2024.01.12 --------------------------------------------------
!------ PLANT --------------------------------------------------------------
      dssat48_struc(nest)%HARVRES(t) = HARVRES
      dssat48_struc(nest)%SENESCE(t) = SENESCE

      dssat48_struc(nest)%dssat48(t)%CANHT = CANHT
      dssat48_struc(nest)%dssat48(t)%EORATIO = EORATIO
      dssat48_struc(nest)%dssat48(t)%NSTRES = NSTRES
      dssat48_struc(nest)%dssat48(t)%PSTRES1 = PSTRES1
      dssat48_struc(nest)%dssat48(t)%PORMIN = PORMIN
      dssat48_struc(nest)%dssat48(t)%RWUMX = RWUMX
      dssat48_struc(nest)%dssat48(t)%KSEVAP = KSEVAP
      dssat48_struc(nest)%dssat48(t)%KTRANS = KTRANS
      dssat48_struc(nest)%dssat48(t)%KUptake = KUptake
      dssat48_struc(nest)%dssat48(t)%PUptake = PUptake
      dssat48_struc(nest)%dssat48(t)%RLV = RLV
      dssat48_struc(nest)%dssat48(t)%FracRts = FracRts 
      dssat48_struc(nest)%dssat48(t)%UNH4 = UNH4
      dssat48_struc(nest)%dssat48(t)%UNO3 = UNO3
      dssat48_struc(nest)%dssat48(t)%STGDOY = STGDOY
      dssat48_struc(nest)%dssat48(t)%XHLAI = XHLAI
      dssat48_struc(nest)%dssat48(t)%XLAI = XLAI
      dssat48_struc(nest)%dssat48(t)%UH2O = UH2O
!---- Pang: 2024.01.17 -------------------------------------------
!---- MGMTOPS ------------------------------------------------------

      dssat48_struc(nest)%FERTDATA(t) = FERTDATA
      dssat48_struc(nest)%OMADATA(t) = OMADATA
      dssat48_struc(nest)%TILLVALS(t) = TILLVALS
      dssat48_struc(nest)%dssat48(t)%HARVFRAC = HARVFRAC
      dssat48_struc(nest)%dssat48(t)%IRRAMT = IRRAMT


      RETURN
      END SUBROUTINE LAND 
