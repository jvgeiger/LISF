C=======================================================================
C  COPYRIGHT 1998-2015 The University of Georgia, Griffin, Georgia
C                      University of Florida, Gainesville, Florida
C                      International Fertilzier Development Center
C                      University of Guelph, Guelph, Ontario
C                      Washington State University, Prosser, Washington
C  ALL RIGHTS RESERVED
C=======================================================================
C=======================================================================
C  PEST, Subroutine
C  Calculates pest damage.
C-----------------------------------------------------------------------
C  REVISION       HISTORY
C  02/23/1998 CHP Written based on Pest damage code in PLANT subroutine
C  01/12/1999 GH  Incorporated into CROPGRO
C  06/19/2001 GH  Added mowing option
C  01/28/2002 GH  Expanded number of pests to 100 (from 40)
C  05/09/2003 CHP Expanded number of pests to 200 (from 100)
C  06/18/2015 GH  Moved MaxPest to ModuleDefs
C-----------------------------------------------------------------------
C  Called by: PLANT
C  Calls:     ASMDM
C             IPPARM
C             IPPEST
C             IPPROG
C             LINDM
C             OPPEST
C             PESTCP
C             ROOTDM
C             SEEDDM
C             VEGDM
C=======================================================================
      SUBROUTINE PEST(CONTROL, ISWITCH,nest, t,       !Pang 2024.01.31 
     &    AREALF, CLW, CSW, LAGSD, LNGPEG, NR2, PGAVL,    !Input
     &    PHTIM, PLTPOP, RTWT, SLA, SLDOT, SOILPROP,      !Input
     &    SSDOT, STMWT, TOPWT, WLFDOT, WTLF, YRPLT,       !Input
     &    RLV, SDNO, SHELN, SWIDOT,                       !Input/Output
     &    VSTAGE, WSHIDT, WTSD, WTSHE,                    !Input/Output
     &    ASMDOT, DISLA, NPLTD, PPLTD,                    !Output
     &    SDDES, WLIDOT, WRIDOT, WSIDOT, SDWT)            !Output

!-----------------------------------------------------------------------
      USE ModuleDefs     !Definitions of constructed variable types, 
                         ! which contain control information, soil
                         ! parameters, hourly weather data.
      USE dssat48_lsmMod
      IMPLICIT NONE
      SAVE
!-----------------------------------------------------------------------
      CHARACTER*1   ISWDIS
      CHARACTER*5   PSTHD(6)
      CHARACTER*5   PID(MAXPEST)
      CHARACTER*5   PCPID(MAXPEST,6)
      CHARACTER*12  FILEP, FILET
      CHARACTER*30  FILEIO
      CHARACTER*80  PATHPE, PATHEX
      INTEGER nest, t
      INTEGER DYNAMIC, LUNIO, MULTI, PCN

      INTEGER TRTNUM, NR2
      INTEGER YRDOY, YRPLT, DAP, DAS
      INTEGER TIMDIF
      INTEGER PNO(6), POBS(6)
      INTEGER PCTID(MAXPEST)
      INTEGER IDAP(6,MAXPEST)

      REAL PHTHRS8, SLA
      REAL PL(6)
      REAL PHTIM(NCOHORTS)
      REAL PDCF1(MAXPEST,6)
      REAL YPL(6,MAXPEST)

C     Leaf Variables
      REAL TLFAD, PLFAD, TLFMD, PLFMD, PCLMT, PCLMA
      REAL CLFM, CLAI, LAIDOT, WLIDOT
      REAL WTLF, SLDOT
      REAL DISLA, DISLAP, AREALF
      REAL CLW, WLFDOT
C     Diseased Leaf Area
      REAL PDLA, TDLA
C     Stem Variables
      REAL WSTMD, PSTMD, PCSTMD, CSTEM
      REAL WSIDOT
      REAL STMWT
      REAL VSTAGE, PVSTGD, VSTGD
      REAL CSW, SSDOT
C     Root Variables
      REAL WRTMD, PRTMD, TRTLV, PRTLV, PRTLF, TRTLF, WRIDOT
      REAL RLVDOT, RLFDOT, CRLV, CRLF, CRTM
      REAL RTWT
      REAL RLV(NL)
      REAL PRLV
C     Seed Variables
      REAL NSDDS, NSDDL, NSDDM, PSDDS, PSDDL, PSDDM, WSDDS, WSDDL, WSDDM
      REAL CSDM, CSDN
      REAL SDDES(NCOHORTS), SDIDOT, SWIDOT
      REAL TSDNOS, TSDNOL, TSDNOM, TSDWTS, TSDWTL, TSDWTM
      REAL LAGSD, LNGPEG
      REAL WSDD,PSDD,SDWT
!     REAL SEEDNO
      REAL WTSD(NCOHORTS), SDNO(NCOHORTS)
C     Shell Variables
      REAL NSHDS, NSHDL, NSHDM, PSHDS, PSHDL, PSHDM, WSHDS, WSHDL, WSHDM
      REAL CSHM, CSHN
      REAL SHIDOT, WSHIDT
      REAL TSHNOS, TSHNOL, TSHNOM, TSHWTS, TSHWTL, TSHWTM
      REAL  WTSHE(NCOHORTS), SHELN(NCOHORTS)
C     Plant Variables
      REAL NPLTD, PPLTD, CPPLTD, TOPWT
      REAL PLTPOP
C     Photosynthesis Variables
      REAL TPSR, PPSR, CASM, ASMDOT
      REAL PGAVL
C** WDB 1/2022 Added Target LAI from remote sensing observations (TLAI)      
      REAL TLAI      
      
!     Redundant with SAVE stmt earlier
!      SAVE CASM, CRLV, CRLF, CRTM
!      SAVE CSDM, CSDN, CSHM, CSHN
!      SAVE CLAI, CLFM, CSTEM

!-----------------------------------------------------------------------
!     Define constructed variable types based on definitions in
!     ModuleDefs.for.

!     The variable "CONTROL" is of type "ControlType".
      TYPE (ControlType) CONTROL

!     The variable "SOILPROP" is of type "SoilType".
      TYPE (SoilType) SOILPROP

!     The variable "ISWITCH" is of type "SwitchType".
      TYPE (SwitchType) ISWITCH

!     Transfer values from constructed data types into local variables.
      DAS     = CONTROL % DAS
      DYNAMIC = CONTROL % DYNAMIC
      FILEIO  = CONTROL % FILEIO
      LUNIO   = CONTROL % LUNIO
      MULTI   = CONTROL % MULTI
      YRDOY   = CONTROL % YRDOY

      ISWDIS  = ISWITCH % ISWDIS
!-----------------------------------------------------------------------
!------ Obtain Vars From Memory ----------------------------------------
!------ Pang-Wei Liu 2024.01.31 ----------------------------------------
!------ IPPEST, IPPARM, IPPROG -----------------------------------------
      FILET = dssat48_struc(nest)%dssat48(t)%FILET
      PATHEX = dssat48_struc(nest)%dssat48(t)%PATHEX
      PHTHRS8 = dssat48_struc(nest)%dssat48(t)%PHTHRS8
      TRTNUM = dssat48_struc(nest)%dssat48(t)%TRTNUM
      PCPID = dssat48_struc(nest)%dssat48(t)%PCPID
      PCTID = dssat48_struc(nest)%dssat48(t)%PCTID
      PDCF1 = dssat48_struc(nest)%dssat48(t)%PDCF1
      PID = dssat48_struc(nest)%dssat48(t)%PID
      IDAP = dssat48_struc(nest)%dssat48(t)%IDAP
      PCN = dssat48_struc(nest)%dssat48(t)%PCN
      PNO= dssat48_struc(nest)%dssat48(t)%PNO
      YPL= dssat48_struc(nest)%dssat48(t)%YPL
!------ PESTCP ---------------------------------------------------------
      NSDDL = dssat48_struc(nest)%dssat48(t)%NSDDL
      NSDDM = dssat48_struc(nest)%dssat48(t)%NSDDM
      NSDDS = dssat48_struc(nest)%dssat48(t)%NSDDS
      NSHDL = dssat48_struc(nest)%dssat48(t)%NSHDL
      NSHDM = dssat48_struc(nest)%dssat48(t)%NSHDM
      NSHDS = dssat48_struc(nest)%dssat48(t)%NSHDS
      TLFAD = dssat48_struc(nest)%dssat48(t)%TLFAD
      TLFMD = dssat48_struc(nest)%dssat48(t)%TLFMD
      TRTLV = dssat48_struc(nest)%dssat48(t)%TRTLV
      WRTMD = dssat48_struc(nest)%dssat48(t)%WRTMD
      WSDDL= dssat48_struc(nest)%dssat48(t)%WSDDL
      WSDDM= dssat48_struc(nest)%dssat48(t)%WSDDM
      WSDDS= dssat48_struc(nest)%dssat48(t)%WSDDS
      WSHDL= dssat48_struc(nest)%dssat48(t)%WSHDL
      WSHDM= dssat48_struc(nest)%dssat48(t)%WSHDM
      WSHDS= dssat48_struc(nest)%dssat48(t)%WSHDS
      CPPLTD= dssat48_struc(nest)%dssat48(t)%CPPLTD
      PCLMA= dssat48_struc(nest)%dssat48(t)%PCLMA
      PCLMT= dssat48_struc(nest)%dssat48(t)%PCLMT
      PCSTMD= dssat48_struc(nest)%dssat48(t)%PCSTMD
      PDLA= dssat48_struc(nest)%dssat48(t)%PDLA
      PLFAD= dssat48_struc(nest)%dssat48(t)%PLFAD
      PLFMD= dssat48_struc(nest)%dssat48(t)%PLFMD
      PPSR= dssat48_struc(nest)%dssat48(t)%PPSR
      PRTLF= dssat48_struc(nest)%dssat48(t)%PRTLF
      PRTLV= dssat48_struc(nest)%dssat48(t)%PRTLV
      PRTMD= dssat48_struc(nest)%dssat48(t)%PRTMD
      PSDDL= dssat48_struc(nest)%dssat48(t)%PSDDL
      PSDDM= dssat48_struc(nest)%dssat48(t)%PSDDM
      PSDDS= dssat48_struc(nest)%dssat48(t)%PSDDS
      PSHDL= dssat48_struc(nest)%dssat48(t)%PSHDL
      PSHDM= dssat48_struc(nest)%dssat48(t)%PSHDM
      PSHDS= dssat48_struc(nest)%dssat48(t)%PSHDS
      PSTMD= dssat48_struc(nest)%dssat48(t)%PSTMD
      PVSTGD= dssat48_struc(nest)%dssat48(t)%PVSTGD
      TDLA= dssat48_struc(nest)%dssat48(t)%TDLA
      TPSR= dssat48_struc(nest)%dssat48(t)%TPSR
      TRTLF= dssat48_struc(nest)%dssat48(t)%TRTLF
      VSTGD= dssat48_struc(nest)%dssat48(t)%VSTGD
      WSTMD= dssat48_struc(nest)%dssat48(t)%WSTMD
      WSDD= dssat48_struc(nest)%dssat48(t)%WSDD
      PSDD= dssat48_struc(nest)%dssat48(t)%PSDD
      PRLV = dssat48_struc(nest)%dssat48(t)%PRLV
!------ ASMDM ---------------------------------------------------------
      CASM = dssat48_struc(nest)%dssat48(t)%CASM
!------ SEEDDM ---------------------------------------------------------
      CSDM = dssat48_struc(nest)%dssat48(t)%CSDM
      CSDN = dssat48_struc(nest)%dssat48(t)%CSDN
      CSHM = dssat48_struc(nest)%dssat48(t)%CSHM
      CSHN = dssat48_struc(nest)%dssat48(t)%CSHN
      SDIDOT = dssat48_struc(nest)%dssat48(t)%SDIDOT
      SHIDOT = dssat48_struc(nest)%dssat48(t)%SHIDOT
      TSDNOL = dssat48_struc(nest)%dssat48(t)%TSDNOL
      TSDNOM = dssat48_struc(nest)%dssat48(t)%TSDNOM
      TSDNOS = dssat48_struc(nest)%dssat48(t)%TSDNOS
      TSDWTL = dssat48_struc(nest)%dssat48(t)%TSDWTL
      TSDWTM = dssat48_struc(nest)%dssat48(t)%TSDWTM
      TSDWTS = dssat48_struc(nest)%dssat48(t)%TSDWTS
      TSHNOL = dssat48_struc(nest)%dssat48(t)%TSHNOL
      TSHNOM = dssat48_struc(nest)%dssat48(t)%TSHNOM
      TSHNOS = dssat48_struc(nest)%dssat48(t)%TSHNOS
      TSHWTL = dssat48_struc(nest)%dssat48(t)%TSHWTL
      TSHWTM = dssat48_struc(nest)%dssat48(t)%TSHWTM
      TSHWTS = dssat48_struc(nest)%dssat48(t)%TSHWTS
!------ VEGDM ---------------------------------------------------------
      CLAI = dssat48_struc(nest)%dssat48(t)%CLAI
      CLFM = dssat48_struc(nest)%dssat48(t)%CLFM
      CSTEM = dssat48_struc(nest)%dssat48(t)%CSTEM
      DISLAP = dssat48_struc(nest)%dssat48(t)%DISLAP
      LAIDOT = dssat48_struc(nest)%dssat48(t)%LAIDOT
!------ ROOTDM ---------------------------------------------------------
      CLAI = dssat48_struc(nest)%dssat48(t)%CLAI
      CRLF = dssat48_struc(nest)%dssat48(t)%CRLF
      CRLV = dssat48_struc(nest)%dssat48(t)%CRLV
      CRTM = dssat48_struc(nest)%dssat48(t)%CRTM
      RLFDOT = dssat48_struc(nest)%dssat48(t)%RLFDOT
      RLVDOT = dssat48_struc(nest)%dssat48(t)%RLVDOT
C***********************************************************************
C***********************************************************************
!     Run Initialization - Called once per simulation
C***********************************************************************
      IF (DYNAMIC .EQ. RUNINIT) THEN
C-----------------------------------------------------------------------
C     Call IPPEST to read data from FILEIO
C-----------------------------------------------------------------------
      CALL IPPEST(
     &    FILEIO, LUNIO,                                  !Input
     &    FILEP, FILET, PATHEX, PATHPE, PHTHRS8, TRTNUM)  !Output
!     Note: PHTHRS8 should be imported from plant routines unless
!       cultivar sections are standardized.  CHP 08/27/2003
C-----------------------------------------------------------------------
C     Subroutine IPPARM reads FILEP, the PEST progress file.
C-----------------------------------------------------------------------
      CALL IPPARM(
     &    FILEP, PATHPE, ISWDIS,                          !Input
     &    PCPID, PCTID, PDCF1, PID)                       !Output

C***********************************************************************
C***********************************************************************
!     Seasonal initialization - run once per season
C***********************************************************************
      ELSEIF (DYNAMIC .EQ. SEASINIT) THEN
C-----------------------------------------------------------------------
C     Subroutine IPPROG reads FILET, the pest time series file.
C-----------------------------------------------------------------------
      IF (ISWDIS .EQ. 'Y') THEN
        CALL IPPROG(CONTROL, 
     &    FILET, PATHEX, PID, YRPLT, TRTNUM,              !Input
     &    IDAP, PCN, PNO, POBS, PSTHD, YPL)               !Output
      ENDIF
C-----------------------------------------------------------------------
C  Initialize whole plant factors
C-----------------------------------------------------------------------
      CALL PESTCP(nest, t,                      !Pang 2024.01.01
     &    PCN, PCPID, PCTID, PDCF1,                       !Input
     &    PL, PLTPOP, PNO, RTWT, SLA, STMWT, TOPWT,        !Input
     &    TSDNOL, TSDNOM, TSDNOS, TSDWTL, TSDWTM, TSDWTS, !Input
     &    TSHNOL, TSHNOM, TSHNOS, TSHWTL, TSHWTM, TSHWTS, !Input
     &    TLAI, VSTAGE, WTLF,                             !Input
     &    NSDDL, NSDDM, NSDDS, NSHDL, NSHDM, NSHDS,       !Input/Output
     &    PPLTD, TLFAD, TLFMD, TRTLV,                     !Input/Output
     &    WRTMD, WSDDL, WSDDM, WSDDS,                     !Input/Output
     &    WSHDL, WSHDM, WSHDS,                            !Input/Output
     &    CPPLTD, NPLTD, PCLMA, PCLMT,                    !Output
     &    PCSTMD, PDLA, PLFAD, PLFMD, PPSR,               !Output
     &    PRTLF, PRTLV, PRTMD, PSDDL, PSDDM, PSDDS,       !Output
     &    PSHDL, PSHDM, PSHDS, PSTMD, PVSTGD,             !Output
     &    TDLA, TPSR, TRTLF, VSTGD, WSTMD,                !Output
     &    SEASINIT,WSDD,PSDD,PRLV)
C-----------------------------------------------------------------------
C  Initialize assimilate, seed, vegetative and root pest damage factors
C-----------------------------------------------------------------------
      CALL ASMDM(
     &    PGAVL, PPSR, TPSR,                              !Input
     &    ASMDOT, CASM,                                   !Output
     &    SEASINIT)                                       !Control

      CALL SEEDDM(nest, t,                      !Pang 2024.02.01
     &    DAS, LAGSD, LNGPEG, NR2, PHTIM, PHTHRS8,        !Input
     &    PSDDL, PSDDM, PSDDS, PSHDL, PSHDM, PSHDS,       !Input
     &    NSDDL, NSDDM, NSDDS, NSHDL, NSHDM, NSHDS,       !Input/Output
     &    SDNO, SHELN, SWIDOT,                            !Input/Output
     &    WSDDL, WSDDM, WSDDS,                            !Input/Output
     &    WSHDL, WSHDM, WSHDS, WTSD, WTSHE, WSHIDT,       !Input/Output
     &    CSDM, CSDN, CSHM, CSHN, SDDES, SDIDOT,          !Output
     &    SHIDOT, TSDNOL, TSDNOM, TSDNOS, TSDWTL,         !Output
     &    TSDWTM, TSDWTS, TSHNOL, TSHNOM, TSHNOS,         !Output
     &    TSHWTL, TSHWTM, TSHWTS,                         !Output
     &    SEASINIT,WSDD,PSDD,SDWT)                             !Control

      CALL VEGDM(SEASINIT,nest,t,               !Pang 2024.02.01
     &    AREALF, CLW, CSW, PCLMT, PCSTMD, PDLA, PLFAD,   !Input
     &    PLFMD, PSTMD, PVSTGD, SLA, SLDOT, SSDOT,        !Input
     &    STMWT, TDLA, TLAI,VSTGD, WLFDOT, WSTMD, WTLF,        !Input
     &    TLFAD, TLFMD, VSTAGE, WLIDOT,                   !Input/Output
     &    CLAI, CLFM, CSTEM, DISLA, DISLAP,               !Output
     &    LAIDOT, WSIDOT)                                 !Output

      CALL ROOTDM(nest, t,                       !Pang 2024.02.01
     &    PRLV,PRTLF, PRTLV, PRTMD, RTWT, SOILPROP, TRTLF,     !Input
     &    RLV, TRTLV, WRTMD,                              !Input/Output
     &    CRLF, CRLV, CRTM, RLFDOT, RLVDOT, WRIDOT,       !Output
     &    SEASINIT)                                       !Control

      CALL OPPEST(CONTROL, ISWITCH, 
     &    ASMDOT, CASM, CLAI, CLFM, CPPLTD, CRLF, CRLV,      
     &    CRTM, CSDM, CSDN, CSHM, CSHN, CSTEM, DISLA, DISLAP,   
     &    LAIDOT, PPLTD, RLFDOT, RLVDOT, SDIDOT, SHIDOT, 
     &    SWIDOT, WLIDOT, WRIDOT, WSIDOT, WSHIDT, YRPLT)     

!***********************************************************************
!***********************************************************************
!     Daily rate calculations
!***********************************************************************
      ELSEIF (DYNAMIC .EQ. RATE) THEN
C-----------------------------------------------------------------------
C     Interpolate between pest damage factors to get today's pest level.
C-----------------------------------------------------------------------
      IF (PCN .LE. 0) RETURN

      DAP   = MAX(0,TIMDIF(YRPLT,YRDOY))
      CALL LINDM(
     &    DAP, IDAP, PCN, YPL,                            !Input
     &    PL)                                             !Output

C-----------------------------------------------------------------------
C     Compute damage applied to each coupling point.
C-----------------------------------------------------------------------
      CALL PESTCP(nest, t,                      !Pang 2024.01.01
     &    PCN, PCPID, PCTID, PDCF1,                       !Input
     &    PL, PLTPOP, PNO, RTWT, SLA, STMWT, TOPWT,       !Input
     &    TSDNOL, TSDNOM, TSDNOS, TSDWTL, TSDWTM, TSDWTS, !Input
     &    TSHNOL, TSHNOM, TSHNOS, TSHWTL, TSHWTM, TSHWTS, !Input
     &    TLAI, VSTAGE, WTLF,                             !Input
     &    NSDDL, NSDDM, NSDDS, NSHDL, NSHDM, NSHDS,       !Input/Output
     &    PPLTD, TLFAD, TLFMD, TRTLV,                     !Input/Output
     &    WRTMD, WSDDL, WSDDM, WSDDS,                     !Input/Output
     &    WSHDL, WSHDM, WSHDS,                            !Input/Output
     &    CPPLTD, NPLTD, PCLMA, PCLMT,                    !Output
     &    PCSTMD, PDLA, PLFAD, PLFMD, PPSR,               !Output
     &    PRTLF, PRTLV, PRTMD, PSDDL, PSDDM, PSDDS,       !Output
     &    PSHDL, PSHDM, PSHDS, PSTMD, PVSTGD,             !Output
     &    TDLA, TPSR, TRTLF, VSTGD, WSTMD,                !Output
     &    RATE,WSDD,PSDD,PRLV)

      CALL SEEDDM(nest, t,                      !Pang 2024.02.01
     &    DAS, LAGSD, LNGPEG, NR2, PHTIM, PHTHRS8,        !Input
     &    PSDDL, PSDDM, PSDDS, PSHDL, PSHDM, PSHDS,       !Input
     &    NSDDL, NSDDM, NSDDS, NSHDL, NSHDM, NSHDS,       !Input/Output
     &    SDNO, SHELN, SWIDOT,                            !Input/Output
     &    WSDDL, WSDDM, WSDDS,                            !Input/Output
     &    WSHDL, WSHDM, WSHDS, WTSD, WTSHE, WSHIDT,       !Input/Output
     &    CSDM, CSDN, CSHM, CSHN, SDDES, SDIDOT,          !Output
     &    SHIDOT, TSDNOL, TSDNOM, TSDNOS, TSDWTL,         !Output
     &    TSDWTM, TSDWTS, TSHNOL, TSHNOM, TSHNOS,         !Output
     &    TSHWTL, TSHWTM, TSHWTS,                         !Output
     &    RATE, WSDD,PSDD,SDWT)                           !Control
C-----------------------------------------------------------------------
C     Call vegetative pest damage routine and compute damage rates
C-----------------------------------------------------------------------
      CALL VEGDM(RATE, nest, t,                     !Pang 2024.02.01
     &    AREALF, CLW, CSW, PCLMT, PCSTMD, PDLA, PLFAD,   !Input
     &    PLFMD, PSTMD, PVSTGD, SLA, SLDOT, SSDOT,        !Input
     &    STMWT, TDLA, TLAI,VSTGD, WLFDOT, WSTMD, WTLF,        !Input
     &    TLFAD, TLFMD, VSTAGE, WLIDOT,                   !Input/Output
     &    CLAI, CLFM, CSTEM, DISLA, DISLAP,               !Output
     &    LAIDOT, WSIDOT)                                 !Output
C-----------------------------------------------------------------------
C     Call root pest damage routine and compute damage rates
C-----------------------------------------------------------------------
      CALL ROOTDM(nest, t,                           !Pang 2024.02.01
     &    PRLV,PRTLF, PRTLV, PRTMD, RTWT, SOILPROP, TRTLF,     !Input
     &    RLV, TRTLV, WRTMD,                              !Input/Output
     &    CRLF, CRLV, CRTM, RLFDOT, RLVDOT, WRIDOT,       !Output
     &    RATE)                                           !Control

!***********************************************************************
!***********************************************************************
!     Daily integration
!***********************************************************************
      ELSEIF (DYNAMIC .EQ. INTEGR) THEN
C-----------------------------------------------------------------------
      IF (PCN .LE. 0) RETURN
C-----------------------------------------------------------------------
C     Call assimilative damage routine to update assimilative damage
!          variables.
C-----------------------------------------------------------------------
      CALL ASMDM(
     &    PGAVL, PPSR, TPSR,                              !Input
     &    ASMDOT, CASM,                                   !Output
     &    INTEGR)                                         !Control
C-----------------------------------------------------------------------
C     Call routine to apply damage to seed and shell
C-----------------------------------------------------------------------
      CALL SEEDDM(nest, t,                          !Pang 2024.02.01
     &    DAS, LAGSD, LNGPEG, NR2, PHTIM, PHTHRS8,        !Input
     &    PSDDL, PSDDM, PSDDS, PSHDL, PSHDM, PSHDS,       !Input
     &    NSDDL, NSDDM, NSDDS, NSHDL, NSHDM, NSHDS,       !Input/Output
     &    SDNO, SHELN, SWIDOT,                            !Input/Output
     &    WSDDL, WSDDM, WSDDS,                            !Input/Output
     &    WSHDL, WSHDM, WSHDS, WTSD, WTSHE, WSHIDT,       !Input/Output
     &    CSDM, CSDN, CSHM, CSHN, SDDES, SDIDOT,          !Output
     &    SHIDOT, TSDNOL, TSDNOM, TSDNOS, TSDWTL,         !Output
     &    TSDWTM, TSDWTS, TSHNOL, TSHNOM, TSHNOS,         !Output
     &    TSHWTL, TSHWTM, TSHWTS,                         !Output
     &    INTEGR, WSDD,PSDD,SDWT)                         !Control
C-----------------------------------------------------------------------
C     Call routine to apply damage to leaf and stem
C-----------------------------------------------------------------------
      CALL VEGDM(INTEGR, nest, t,                  !Pang 2024.02.01
     &    AREALF, CLW, CSW, PCLMT, PCSTMD, PDLA, PLFAD,   !Input
     &    PLFMD, PSTMD, PVSTGD, SLA, SLDOT, SSDOT,        !Input
     &    STMWT, TDLA, TLAI,VSTGD, WLFDOT, WSTMD, WTLF,        !Input
     &    TLFAD, TLFMD, VSTAGE, WLIDOT,                   !Input/Output
     &    CLAI, CLFM, CSTEM, DISLA, DISLAP,               !Output
     &    LAIDOT, WSIDOT)                                 !Output
C-----------------------------------------------------------------------
C     Call root pest damage routine and compute damage factors
C-----------------------------------------------------------------------
      CALL ROOTDM(nest, t,                         !Pang 2024.02.01
     &    PRLV, PRTLF, PRTLV, PRTMD, RTWT, SOILPROP, TRTLF,     !Input
     &    RLV, TRTLV, WRTMD,                              !Input/Output
     &    CRLF, CRLV, CRTM, RLFDOT, RLVDOT, WRIDOT,       !Output
     &    INTEGR)                                         !Control

!***********************************************************************
!***********************************************************************
!     OUTPUT/SEASEND
!***********************************************************************
      ELSEIF (DYNAMIC .EQ. OUTPUT .OR. DYNAMIC .EQ. SEASEND) THEN
C-----------------------------------------------------------------------
      CALL OPPEST(CONTROL, ISWITCH, 
     &    ASMDOT, CASM, CLAI, CLFM, CPPLTD, CRLF, CRLV,      
     &    CRTM, CSDM, CSDN, CSHM, CSHN, CSTEM, DISLA, DISLAP,   
     &    LAIDOT, PPLTD, RLFDOT, RLVDOT, SDIDOT, SHIDOT, 
     &    SWIDOT, WLIDOT, WRIDOT, WSIDOT, WSHIDT, YRPLT)     

!***********************************************************************
!***********************************************************************
!     END OF DYNAMIC IF CONSTRUCT
!***********************************************************************
      ENDIF
!***********************************************************************
!------ ASSIGN Vars To Memory ----------------------------------------
!------ Pang-Wei Liu 2024.01.31 ----------------------------------------
!------ IPPEST, IPPARM, IPPROG -----------------------------------------
      dssat48_struc(nest)%dssat48(t)%FILET = FILET
      dssat48_struc(nest)%dssat48(t)%PATHEX = PATHEX
      dssat48_struc(nest)%dssat48(t)%PHTHRS8 = PHTHRS8
      dssat48_struc(nest)%dssat48(t)%TRTNUM = TRTNUM
      dssat48_struc(nest)%dssat48(t)%PCPID = PCPID
      dssat48_struc(nest)%dssat48(t)%PCTID = PCTID
      dssat48_struc(nest)%dssat48(t)%PDCF1 = PDCF1
      dssat48_struc(nest)%dssat48(t)%PID = PID
      dssat48_struc(nest)%dssat48(t)%IDAP = IDAP
      dssat48_struc(nest)%dssat48(t)%PCN = PCN
      dssat48_struc(nest)%dssat48(t)%PNO = PNO
      dssat48_struc(nest)%dssat48(t)%YPL = YPL
!------ PESTCP ---------------------------------------------------------
      dssat48_struc(nest)%dssat48(t)%NSDDL = NSDDL
      dssat48_struc(nest)%dssat48(t)%NSDDM = NSDDM
      dssat48_struc(nest)%dssat48(t)%NSDDS = NSDDS
      dssat48_struc(nest)%dssat48(t)%NSHDL = NSHDL
      dssat48_struc(nest)%dssat48(t)%NSHDM = NSHDM
      dssat48_struc(nest)%dssat48(t)%NSHDS = NSHDS
      dssat48_struc(nest)%dssat48(t)%TLFAD = TLFAD
      dssat48_struc(nest)%dssat48(t)%TLFMD = TLFMD
      dssat48_struc(nest)%dssat48(t)%TRTLV = TRTLV
      dssat48_struc(nest)%dssat48(t)%WRTMD = WRTMD
      dssat48_struc(nest)%dssat48(t)%WSDDL = WSDDL
      dssat48_struc(nest)%dssat48(t)%WSDDM = WSDDM
      dssat48_struc(nest)%dssat48(t)%WSDDS = WSDDS
      dssat48_struc(nest)%dssat48(t)%WSHDL = WSHDL
      dssat48_struc(nest)%dssat48(t)%WSHDM = WSHDM
      dssat48_struc(nest)%dssat48(t)%WSHDS = WSHDS
      dssat48_struc(nest)%dssat48(t)%CPPLTD = CPPLTD
      dssat48_struc(nest)%dssat48(t)%PCLMA = PCLMA
      dssat48_struc(nest)%dssat48(t)%PCLMT = PCLMT
      dssat48_struc(nest)%dssat48(t)%PCSTMD = PCSTMD
      dssat48_struc(nest)%dssat48(t)%PDLA = PDLA
      dssat48_struc(nest)%dssat48(t)%PLFAD = PLFAD
      dssat48_struc(nest)%dssat48(t)%PLFMD = PLFMD
      dssat48_struc(nest)%dssat48(t)%PPSR = PPSR
      dssat48_struc(nest)%dssat48(t)%PRTLF = PRTLF
      dssat48_struc(nest)%dssat48(t)%PRTLV = PRTLV
      dssat48_struc(nest)%dssat48(t)%PRTMD = PRTMD
      dssat48_struc(nest)%dssat48(t)%PSDDL = PSDDL
      dssat48_struc(nest)%dssat48(t)%PSDDM = PSDDM
      dssat48_struc(nest)%dssat48(t)%PSDDS = PSDDS
      dssat48_struc(nest)%dssat48(t)%PSHDL = PSHDL
      dssat48_struc(nest)%dssat48(t)%PSHDM = PSHDM
      dssat48_struc(nest)%dssat48(t)%PSHDS = PSHDS
      dssat48_struc(nest)%dssat48(t)%PSTMD = PSTMD
      dssat48_struc(nest)%dssat48(t)%PVSTGD= PVSTGD
      dssat48_struc(nest)%dssat48(t)%TDLA = TDLA
      dssat48_struc(nest)%dssat48(t)%TPSR = TPSR
      dssat48_struc(nest)%dssat48(t)%TRTLF = TRTLF
      dssat48_struc(nest)%dssat48(t)%VSTGD = VSTGD
      dssat48_struc(nest)%dssat48(t)%WSTMD = WSTMD
      dssat48_struc(nest)%dssat48(t)%WSDD = WSDD
      dssat48_struc(nest)%dssat48(t)%PSDD = PSDD
      dssat48_struc(nest)%dssat48(t)%PRLV = PRLV
!------ ASMDM ---------------------------------------------------------
      dssat48_struc(nest)%dssat48(t)%CASM = CASM
!------ SEEDDM ---------------------------------------------------------
      dssat48_struc(nest)%dssat48(t)%CSDM = CSDM
      dssat48_struc(nest)%dssat48(t)%CSDN = CSDN
      dssat48_struc(nest)%dssat48(t)%CSHM = CSHM
      dssat48_struc(nest)%dssat48(t)%CSHN =CSHN
      dssat48_struc(nest)%dssat48(t)%SDIDOT = SDIDOT
      dssat48_struc(nest)%dssat48(t)%SHIDOT = SHIDOT
      dssat48_struc(nest)%dssat48(t)%TSDNOL = TSDNOL
      dssat48_struc(nest)%dssat48(t)%TSDNOM = TSDNOM
      dssat48_struc(nest)%dssat48(t)%TSDNOS = TSDNOS
      dssat48_struc(nest)%dssat48(t)%TSDWTL = TSDWTL
      dssat48_struc(nest)%dssat48(t)%TSDWTM = TSDWTM
      dssat48_struc(nest)%dssat48(t)%TSDWTS = TSDWTS
      dssat48_struc(nest)%dssat48(t)%TSHNOL = TSHNOL
      dssat48_struc(nest)%dssat48(t)%TSHNOM = TSHNOM
      dssat48_struc(nest)%dssat48(t)%TSHNOS = TSHNOS
      dssat48_struc(nest)%dssat48(t)%TSHWTL = TSHWTL
      dssat48_struc(nest)%dssat48(t)%TSHWTM = TSHWTM
      dssat48_struc(nest)%dssat48(t)%TSHWTS = TSHWTS
!------ VEGDM ---------------------------------------------------------
      dssat48_struc(nest)%dssat48(t)%CLAI = CLAI
      dssat48_struc(nest)%dssat48(t)%CLFM = CLFM
      dssat48_struc(nest)%dssat48(t)%CSTEM = CSTEM
      dssat48_struc(nest)%dssat48(t)%DISLAP = DISLAP
      dssat48_struc(nest)%dssat48(t)%LAIDOT = LAIDOT
!------ ROOTDM ---------------------------------------------------------
      dssat48_struc(nest)%dssat48(t)%CLAI = CLAI
      dssat48_struc(nest)%dssat48(t)%CRLF = CRLF
      dssat48_struc(nest)%dssat48(t)%CRLV = CRLV
      dssat48_struc(nest)%dssat48(t)%CRTM = CRTM
      dssat48_struc(nest)%dssat48(t)%RLFDOT = RLFDOT
      dssat48_struc(nest)%dssat48(t)%RLVDOT = RLVDOT
      END  ! SUBROUTINE PEST

!-----------------------------------------------------------------------
!     PEST Variable Definitions
!-----------------------------------------------------------------------
! AREALF     Area of leaves (one side) per unit ground area
!              (cm2[leaf] / m2[ground])
! ASMDOT     Daily assimilative damage (g[CH2O] /m2 / d)
! CASM       Cumulative assimilate damage (g[CH2O] / m2)
! CLAI       Cumulative leaf area index destroyed (m2/m2)
! CLFM       Cumulative leaf mass destroyed  (g/m2)
! CLW        Cumulative leaf growth (g[leaf]/m2)
! CPPLTD     Cumulative percent of plants destroyed (%)
! CRLF       Cumulative root length flux (cm root / cm2 ground)
! CRLV       Cumulative root length density (cm root / cm3 soil)
! CRTM       Cumulative root mass (g/m2)
! CSDM       Cumulative seed mass destroyed (g/m2)
! CSDN       Cumulative number of seeds destroyed (#/m2)
! CSHM       Cumulative shell mass destroyed (g/m2)
! CSHN       Cumulative number of shells destroyed (#/m2)
! CSTEM      Cumulative stem mass destroyed (g/m2)
! CSW        Cumulative stem growth (g[stem]/m2)
! DAP        Number of days after planting (d)
! DAS        Days after start of simulation (d)
! DISLA      Diseased leaf area (cm2[leaf]/m2[ground]/d)
! DISLAP     Percent diseased leaf area (%/d)
! DLAYR(L)   Soil thickness in layer L (cm)
! FILEIO     Filename for input file (e.g., IBSNAT35.INP) 
! FILEP      Filename for pest coefficient file 
! FILET      Pest time series file 
! IDAP(I,J)  Day of pest damage for pest I, observation J
!              (days after planting)
! LAGSD      Time required between shell growth and seed growth, per cohort
!              (Photo-thermal days)
! LAIDOT     Daily change in leaf area index due to pest damage (m2/m2/d)
! LNGPEG     Time between start of peg (full flower) and shell formation 
!              (for peanuts only).  Defines slow growth period.
!              (Photo-thermal days)
! LUNIO      Logical unit number for FILEIO 
! MULTI      Current simulation year (=1 for first or single simulation, 
!              =NYRS for last seasonal simulation) 
! NL         Maximum number of soil layers = 20 
! NLAYR      Actual number of soil layers 
! NPLTD      Number of plants destroyed (#/m2/d)
! NR2        Day when 50% of plants have one peg (peanuts only) (d)
! NSDDL      Daily number of large seeds damaged (#/m2/d)
! NSDDM      Daily number of mature seeds damaged (#/m2/d)
! NSDDS      Daily number of small seeds damaged (#/m2/d)
! NSHDL      Daily number of large shells damaged (#/m2/d)
! NSHDM      Daily number of mature shells damaged (#/m2/d)
! NSHDS      Daily number of small shells damaged (#/m2/d)
! OUTD       File name for pest damage output file (e.g., PEST.OUT) 
! PATHPE     Path name for pest information file 
! PCLMA      Percent observed leaf mass (WTLF) damage (%)
! PCLMT      Percent of total leaf mass (WTLF + senescence) destroyed (%)
! PCN        Number of pests in FILET 
! PCPID(I,J) Pest coupling points identification code for pest I, coupling 
!              point J 
! PCSTMD     Observed cumulative percentage stem mass damage (%)
! PCTID(I)   Pest damage characterization method for pest I: (1) absolute 
!              daily damage rate, (2) percent observed damage, (3) daily 
!              percent damage rate, (4) absolute daily damage rate w/ 
!              preference and competition. 
! PDCF1(I,J) Pest damage coefficient associated with pest I, coupling point 
!              J 
! PDLA       Percent diseased leaf area (%)
! PGAVL      Total available CH2O available for growth & respiration
!              (g[CH2O] / m2)
! PHTHRS8    Threshold time that must accumulate in phase 8 for the next 
!              stage to occur.  Equivalent to PHTHRS(8) in Subroutine 
!              PHENOLOG. (photothermal days)
! PHTIM      Cumulative photothermal time ages of seeds and shells 
! PID(I)     Pest identification header from FILEP (max 200) 
! PL(I)      Pest level from pest progress file (FILET) for pest I (varies)
! PLFAD      Daily percent leaf area damage (%/d)
! PLFMD      Percent leaf mass damage (%/d)
! PLTPOP     Plant population (# plants / m2)
! PNO(I)     Pest number from File P for pest I 
! POBS(I)    Number of observations for pest I 
! PPLTD      Percent plants destroyed  (%/m2/d)
! PPSR       Assimilate damage in daily percent photosynthesis reduction
!              (%/d)
! PRLV       Percent of root length volume destroyed, %
! PRTLF      Percent of root flux destroyed (%/d)
! PRTLV      Daily percent reduction in root length volume  (%/d)
! PRTMD      Daily percent root mass damage (%/d)
! PSDD       Percent of seed mass destroyed (%/d)
! PSDDL      Percent large seed mass damage (%/d)
! PSDDM      Percent mature seed mass damage (%/d)
! PSDDS      Percent small seed mass damage (%/d)
! PSHDL      Percent large shell mass damage (%/d)
! PSHDM      Percent mature shell mass damage (%/d)
! PSHDS      Percent small shell mass damage (%/d)
! PSTHD(I)   Column heading for each pest progress curve in FILET (max 6) 
! PSTMD      Daily percent stem mass damage (%)
! PVSTGD     Percent V-stage damage (%)
! RLFDOT     Daily root length flux damage (cm root / cm2 ground)
! RLV(L)     Root length density for soil layer L (cm[root] / cm3[soil])
! RLVDOT     Daily root length damage (cm root / cm3 soil)
! RTWT       Dry mass of root tissue, including C and N
!              (g[root] / m2[ground])
! SDDES(J)   Number of seeds destroyed today in cohort J when shells are 
!              not destroyed (#/m2/day)
! SDIDOT     Number of seeds destroyed on the current day (#/m2/d)
! SDNO(J)    Number of seeds for cohort J (#/m2)
! SDWT       Seed weight, g/m2
! SHELN(J)   Number of shells for cohort J (#/m2)
! SHIDOT     Number of shells destroyed on current day (#/m2/d)
! SLA        Specific leaf area (cm2[leaf] / m2[ground])
! SLDOT      Defoliation due to daily leaf senescence (g/m2/day)
! SSDOT      Daily senescence of petioles (g / m2 / d)
! STMWT      Dry mass of stem tissue, including C and N
!              (g[stem] / m2[ground)
! SWIDOT     Daily seed mass damage (g/m2/day)
! TDLA       Total diseased leaf area (cm2/m2)
! TIMDIF     Integer function which calculates the number of days between 
!              two Julian dates (da)
! TLFAD      Total leaf area damage (cm2/cm2/d)
! TLFMD      Total leaf mass damage (g/m2/day)
! TOPWT      Total above ground biomass (g/m2)
! TPSR       Daily absolute assimilate damage (g[CH2O]/m2/d)
! TRTLF      Total root flux destroyed (cm/cm2/d)
! TRTLV      Daily absolute reduction in root length volume  (cm/cm3/d)
! TRTNUM     Treatment number being simulated (from FILEX) 
! TSDNOL     Total number of large seeds (#/m2)
! TSDNOM     Total number of mature seeds (#/m2)
! TSDNOS     Total number of small seeds (#/m2)
! TSDWTL     Seed mass for large seeds (g/m2)
! TSDWTM     Seed mass for mature seeds (g/m2)
! TSDWTS     Seed mass for small seeds (g/m2)
! TSHNOL     Number shells with large seeds (#/m2)
! TSHNOM     Number shells with mature seeds (#/m2)
! TSHNOS     Number shells with small seeds (#/m2)
! TSHWTL     Shell mass with large seeds (g/m2)
! TSHWTM     Shell mass with mature seeds (g/m2)
! TSHWTS     Shell mass with small seeds (g/m2)
! VSTAGE     Number of nodes on main stem of plant (nodes)
! VSTGD      Absolute daily V-stage damage (nodes/day)
! WLFDOT     Leaf weight losses due to freezing (g[leaf]/m2-d)
! WLIDOT     Daily pest or freeze damage to leaf mass (g/m2/day)
! WRIDOT     Daily pest damage to root mass (g/m2/day)
! WRTMD      Daily absolute root mass reduction  (g/m2/day)
! WSDD       Daily weight of seeds destroyed (g/m2/d)
! WSDDL      Daily percent of large seed damaged (g/m2/day)
! WSDDM      Daily percent of mature seed damaged (g/m2/day)
! WSDDS      Daily percent of small seed damaged (g/m2/day)
! WSHDL      Daily mass of large shell damaged (g/m2/day)
! WSHDM      Daily mass of mature shell damaged (g/m2/day)
! WSHDS      Daily mass of small shell damaged (g/m2/day)
! WSHIDT     Weight of shell tissue consumed by pests today (g[shell]/m2-d)
! WSIDOT     Daily pest damage to stem mass (g/m2/day)
! WSTMD      Daily absolute stem damage (g/m2/day)
! WTLF       Dry mass of leaf tissue including C and N
!              (g[leaf] / m2[ground])
! WTSD(J)    Seed mass  for cohort J (g/m2)
! WTSHE(J)   Shell mass  for cohort J (g/m2)
! YPL(I,J)   Array for storage of pest data for pest or damage type I, 
!              observation J 
! YRDOY      Current day of simulation (YYDDD)
! YRPLT      Planting date (YYDDD)
!-----------------------------------------------------------------------
!     End Subroutine PEST
!-----------------------------------------------------------------------

