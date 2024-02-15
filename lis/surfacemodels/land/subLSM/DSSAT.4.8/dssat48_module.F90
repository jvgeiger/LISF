!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

module dssat48_module
!BOP
!
! !MODULE: dssat48_module
!
! !DESCRIPTION:
!  The code in this file provides a description of the
!  data structure containing the DSSAT model 1-d variables.
!  The variables specified in the data structure include:
!
!  \begin{description}      
!   \item[n]
!     nest id. unit: -
!
! !REVISION HISTORY:
!  11 May 2023: Pang-Wei Liu
!  31 Aug 2023: J. Erlingis; Add forcing
!EOP
    implicit none
    private
    type, public :: dssat48dec
        !-------------------------------------------------------------------------
        ! forcing
        !-------------------------------------------------------------------------
        real               :: tair, tmax, tmin       !from LIS forcing
        real               :: qair, swdown, lwdown
        real               :: uwind, vwind, wndspd
        real               :: psurf, rainf, snowf
        real               :: tdew, totprc
        real               :: forc_pres, forc_precip !forc_* is passed to DSSAT
        real               :: forc_tmax, forc_tmin
        real               :: forc_tdew, forc_swrad
        real               :: forc_wind      
        !-------------------------------------------------------------------------
        ! spatial parameter
        !-------------------------------------------------------------------------
        !real               :: CDL
        CHARACTER*10 SLNO   !SOIL ID Pang 2024.02.07
        real                :: lat, lon, elev
        !-------------------------------------------------------------------------
        ! multilevel spatial parameter
        !-------------------------------------------------------------------------
        !real, pointer      :: ALB(:) 
        !-------------------------------------------------------------------------
        ! state
        !-------------------------------------------------------------------------
        !real, pointer      :: SM(:)
        !real, pointer      :: SNOWRHO(:)
        !-------------------------------------------------------------------------
        ! output
        !-------------------------------------------------------------------------
        !real               :: SM
        !real               :: YIELD
        !real               :: LAI
        !-------------------------------------------------------------------------
        ! TIME CONTROL
        ! Pang 2023.09.29
        !-------------------------------------------------------------------------
         INTEGER            :: EXPNO, TRTALL
         INTEGER            :: NREPS, YRDOY_END
         CHARACTER*1        :: RNMODE
         CHARACTER*80       :: PATHEX
         CHARACTER*92       :: FILEX_P
         CHARACTER*120      :: FILECTL
        !------------------------------------------------------------------------
        !  LAND
        !------------------------------------------------------------------------
         INTEGER            :: YRPLT, MDATE, YREND        
        !-------------------------------------------------------------------------
        !  LAND - SOIL
        !  Pang 2023.09.29
        !------------------------------------------------------------------------
          REAL, DIMENSION(0:20) :: SomLitC
          REAL, DIMENSION(0:20,3) :: SomLitE
          REAL, DIMENSION(20) :: SKi_Avail, SW, NH4_plant, NO3_Plant, UPPM, SPi_AVAIL, SWDELTS
          REAL                :: SNOW, WINF
        !-------------------------------------------------------------------------
        !  SOIL
        !  Pang 2023.10.03
        !-------------------------------------------------------------------------
          REAL, DIMENSION(20) :: SomLit, SPi_Labile, NO3, NH4, DRN
          REAL, DIMENSION(0:20) :: SSOMC, LITC, newCO2
          REAL, DIMENSION(0:20,3) :: IMM, MNR
          REAL                  :: TDFC, DRAIN
          INTEGER               :: TDLNO
        !-------------------------------------------------------------------------
        !  SOIL - SOILDYN
        ! Pang 2023.09.18
        !-------------------------------------------------------------------------
          REAL :: CN_INIT, CRAIN, LCRAIN, SUMKE
          REAL, DIMENSION(20) :: BD_INIT, DLAYR_INIT, DS_INIT, DUL_INIT
          REAL, DIMENSION(20) :: LL_INIT, SWCN_INIT, SAT_INIT, SW_INIT
          REAL, DIMENSION(20) :: SomLit_INIT, SOM_PCT_init, BD_calc_init
          REAL, DIMENSION(20) :: BD_SOM, DLAYR_SOM, DS_SOM, DUL_SOM, LL_SOM
          REAL, DIMENSION(20) :: OC_INIT, TOTN_INIT, TotOrgN_init !Pang: May not need
          REAL, DIMENSION(0:20) :: KECHGE
          LOGICAL :: TILLED, SOILDYN_FIRST
        !-------------------------------------------------------------------------
        !  SOIL - WATBAL
        !  Pang 2023.10.11
        !-------------------------------------------------------------------------
          REAL :: ICWD_INIT !(IPWBAL modeul)
          REAL :: WB_CRAIN,TDRAIN, TRUNOF, TSWINI!TSW !(WBSUM)
          REAL :: WTDEP, TDFD, WATAVL
          REAL, DIMENSION(20) :: DLAYR_YEST, SWDELTT, SWDELTL
        !-------------------------------------------------------------------------
        !  SOIL - SoilOrg
        !  Pang 2023.10.05
        !-------------------------------------------------------------------------
          REAL :: SEN_AM, SEN_PRCEL, SEN_PRCHO, SEN_PRLIG,SEN_EXTFAC, SEN_WATFAC
          REAL :: SWEF, CUMSENWT, CUMSENN , CUMSENP, RDCHO, RDCEL, RDLIG, DMINR, DSNC
          REAL :: PRCEL, PRCHO, PRLIG, AM, EXTFAC, WATFAC
          REAL :: FOM(20), FON(20), FOP(20), FPOOL(20,3), SSOME(0:20,3)
          REAL :: CNRAT(20)
          REAL :: ACCCO2(0:1)
          CHARACTER(LEN=2) :: PREV_CROP
        !------------------------------------------------------------------------
        !  SOIL - CENTURY
        !  Pang 2023.10.10
        !------------------------------------------------------------------------
          !NELEM = 3
          REAL, DIMENSION(3) :: CEDAM, CES21T, CES23T, CES3M, CES3T
          REAL, DIMENSION(3) :: CES3X, CESTR, FRDAE
          REAL, DIMENSION(0:0,3) :: CES21I
          REAL, DIMENSION(0:1,3) :: CES1M, CES1T, CES1X, CES21M, CES21S
          REAL, DIMENSION(0:1,3) :: CES21X, CES23LM
          REAL, DIMENSION(0:1,3) :: CES23LX, CES23M, CES23X
          REAL, DIMENSION(0:1,3) :: CES2LI, CES2LM, CES2LS, CES2LX
          REAL, DIMENSION(0:1) :: CO2MET, DECMET, DECS1, DECSTR, LIGSTR
          REAL :: DECS2(1), DECS3(1), METABE(0:20,3)
          REAL :: SOM1E(0:20,3),SOM23E(20,3), SOM2E(20,3), SOM3E(20,3)
          REAL :: STRUCE(0:20,3)
          REAL :: CO2S2, CO2S3, CULMETQ, CULS1Q, CULS2Q, CULS3Q, CULSTRQ
          REAL :: FRMETI, FRMETS, RESDAX, SENESSUMC, SENESSUMN, SENESSUMP
          REAL, DIMENSION(0:20) :: FRLSTR, CO2S1, LIGC, METABC, SOM1C
          REAL, DIMENSION(0:20) :: STRUCC
          REAL, DIMENSION(20) :: S1S3, S2S3, SOM23C, SOM2C, SOM3C, TXS1
          INTEGER DISTURBENDMAX, DISTURBNUM
          REAL :: DISTURBDEPTHMAX
          INTEGER, DIMENSION(9000*3) :: DISTURBEND  !Do we really need 9000*3?
          REAL, DIMENSION(9000*3) :: DISTURBDEPTH
          REAL, DIMENSION(0:1,2) :: CO2STR(0:1,2)
        !------------------------------------------------------------------------
        !  SOIL - SoilOrg/CENTURY - MethaneDynamics
        !  Pang 2023.10.05
        !------------------------------------------------------------------------
          REAL :: Buffer(20,2), CH4Stored_Y, CumNewCO2, lamda_rho
          LOGICAL FirstTime
        !------------------------------------------------------------------------
        !  SOIL - SoilNi
        !  Pang 2023.10.12
        !------------------------------------------------------------------------
           LOGICAL :: IUON
           INTEGER :: LFD10
           REAL :: TNOM, CMINERN, CIMMOBN, WTNUP, CumSumFert
           REAL :: CLeach, CNTILEDR !(NFLUX)
           !REAL :: CNITRIFY, CN2Onitrif, CNOflux !(N2O_data)
           REAL :: TFNITY(20), SNH4(20), SNO3(20), UREA(20)
           REAL :: ALGFIX, CUMFNRO, TOTAML !(FLOOD_CHEM)
           REAL :: CNOX, TNOXD !(Denit)
        !------------------------------------------------------------------------
        !  SOIL - SoilPi
        !  Pang 2023.10.16
        !------------------------------------------------------------------------
           REAL, DIMENSION(20) :: FracRtsY
           REAL, DIMENSION(20) :: SPiLabRts, SPiLabNoRts
           REAL, DIMENSION(20) :: SPi_Active, SPiActRts, SPiActNoRts
           REAL, DIMENSION(20) :: SPi_Stable, SPiStaRts, SPiStaNoRts
           REAL, DIMENSION(20) :: SPi_Soluble, SPiSolRts, SPiSolNoRts
           REAL, DIMENSION(20) :: PiActRts, PiStaRts, FracPSoln
           REAL :: CMINERP, CIMMOBP, PERROR
           LOGICAL :: SOILPI_FIRST
        !------------------------------------------------------------------------
        !  SOIL - SoilKi
        !  Pang 2023.10.16
        !------------------------------------------------------------------------
           REAL, DIMENSION(20) :: SKi_Tot
        !-------------------------------------------------------------------------
        !  LAND - SPAM
        !  Pang 2024.01.12
        !------------------------------------------------------------------------
          REAL :: EO, EOP, EOS, EP, ES, SRFTEMP, TRWU, TRWUP 
          REAL, DIMENSION(20) :: RWU, ST, SWDELTX, UPFLOW, SWDELTU
        !------------------------------------------------------------------------
        !  In SPAM
        !  Pang 2024.02.05
        !-------------------------------------------------------------------------
           REAL :: EF, CEF, EM, CEM, CEO, CEP, CES, CET, EVAP, CEVAP, ES_LYR(20), ET0(24)
            !---------------------------------------------------------------------
            !  SPAM -> In ETPHOT.for (Pang 2024.02.12)
            !---------------------------------------------------------------------
               CHARACTER TYPPGL*3
               CHARACTER(len=2) PGPATH
               REAL AZIR, BETN, CEC, DLAYR2(20), DUL2(20), DULE, LFANGD(3), LL2(20), LLE, &
                    LWIDTH, PALBW, RCUTIC, SAT2(20), SCVIR, SCVP, &
                    XSW(20,3), YSCOND(20,3), YSHCAP(20,3), FNPGL(4), &
                    LMXREF, NSLOPE, QEREF, SLWREF, SLWSLO, XLMAXT(6), YLMAXT(6), &
                    CCNEFF, CICAD, CMXSF, CQESF, PCINPD, PGNOON, PCINPN, SLWSLN, &
                    SLWSHN, PNLSLN, PNLSHN, LMXSLN, LMXSHN, TSHR(20), TEMPN, TSRF(3), &
                    TSRFN(3), TAV, TAMP, COLDSTR
            !---------------------------------------------------------------------
            !  SPAM -> In STEMPi_EPIC.for (Pang 2024.02.12)
            !---------------------------------------------------------------------
               INTEGER WetDay(30), NDays
               REAL WFT, CV, BCV1, BCV2, BCV, X2_AVG
            !---------------------------------------------------------------------
            !  SPAM -> In STEMP.for (Pang 2024.02.12)
            !---------------------------------------------------------------------
               REAL, DIMENSION(20) :: SWI, DSI, DLI, DSMID
               REAL HDAY, TBD, TLL, TSW, TDL, CUMDPT, PESW, ABD, FX, DP, WW, &
                    ALBEDO, TMA(5), ATOT
            !---------------------------------------------------------------------
            !  SPAM -> In ROOTWU.for (Pang 2024.02.12)
            !---------------------------------------------------------------------
               REAL  SWCON2(20), ROOTWU_TSS(20)
            !---------------------------------------------------------------------
            !  SPAM -> In SOILEV.for (Pang 2024.02.12)
            !---------------------------------------------------------------------
               REAL SUMES1, SUMES2, T_prm
        !-------------------------------------------------------------------------
        !  LAND - Management
        !  Pang 2024.01.17
        !-------------------------------------------------------------------------
          REAL IRRAMT     
          REAL, DIMENSION(2) :: HARVFRAC
        !-------------------------------------------------------------------------
        !  LAND - PLANT
        !  Pang 2024.01.17
        !-------------------------------------------------------------------------
          REAL CANHT, EORATIO, NSTRES, PSTRES1, PORMIN, RWUMX
          REAL KSEVAP, KTRANS
          REAL, Dimension(20) :: KUptake, PUptake, RLV, FracRts, UNH4, UNO3
          INTEGER STGDOY(20)
          REAL XHLAI, XLAI
          REAL UH2O(20)

        !-------------------------------------------------------------------------
        !  PLANT
        !  Pang 2024.01.18
        !-------------------------------------------------------------------------
           REAL KCAN, KEP
           LOGICAL FixCanht
           !------------------------------------------------------------------------
           ! PLANT - CROPGRO
           ! Pang 2024.01.24
           !------------------------------------------------------------------------
           !---- PLANT - IPPLNT Module
             CHARACTER*1 DETACH
             CHARACTER*6 ECONO
             CHARACTER(LEN=92) FILECC, FILEGC
             INTEGER NOUTDO
             REAL CADPR1, CMOBMX, FRCNOD, FREEZ1, FREEZ2, KC_SLOPE, PCARSH, &
                  PCH2O, PLIPSH, PLIGSD, PLIGSH, PMINSD, PMINSH, POASD, POASH, PROLFI, &
                  PRORTI, PROSHI, PROSTI, R30C2, RCH2O, RES30C, RFIXN, RLIG, RLIP, RMIN, &
                  RNH4C, RNO3C, ROA, RPRO, RWUEP1, TTFIX
           !------------------------------------------------------------------------
           !---- CROPGRO - PHOTO------------------------------------------------------
              REAL AGEFAC, PG
              !---- In PHOTP.for
                 CHARACTER*3  TYPPGN, TYPPGT
                 REAL CCEFF, CCMAX, CCMP, FNPGN(4), FNPGT(4), LNREF, PARMAX, PHTHRS10, PHTMAX, ROWSPC, &
                      XPGSLW(15), YPGSLW(15), PGLFMX, CUMSTR
           !-----------------------------------------------------------------------
           !---- CROPGROW - PHENOL ----------------------------------
           !-----------------------------------------------------------------------

           !-----------------------------------------------------------------------
           ! PLANT - MZ_CERES
           ! Pang 2024.01.24
           !------------------------------------------------------------------------
            CHARACTER*10 STNAME(20)

           !-----------------------------------------------------------------------
           ! PLANT - MZ_CERES - MZ_PHENOL & MZ_IX_PHENOL
             INTEGER ISDATE, ISTAGE, YREMRG, CDAY
             REAL CUMDTT, DTT, EARS, GPP, SUMDTT, XNTI,TLNO,XSTAGE,RUE, P3, TSEN, SeedFrac, VegFrac 
             REAL PEAR, PSTM, GDDAE, Z2STAGE !For MZ_IX_PHENOL only
             !----- In  MZ_PHENOL.for & MZ_IX_PHENOL.for
                INTEGER PATHL, NDAS, L0, L
                REAL PLTPOP,SDEPTH,P1,P2,P5,G2, G3, PHINT, DSGT, DGET, SWCG, TBASE, TOPT,ROPT,P2O,DJTI,GDDE,DSGFT, &
                     DUMMY,TNSOIL,TMSOIL,TH,TEMPCX,TEMPCR,TDSOIL,SWSD,SNUP,SNDN,S1,RATEIN,PSKER,PDTT, &
                     P9, DLV, DEC, C1, ACOEF, DOPT

           !-----------------------------------------------------------------------
           ! PLANT - MZ_CERES - MZ_GROSUB & MZ_IX_GROSUB & SW_GROSUB
           ! Pang 2024.01.30
             INTEGER LEAFNO, RSTAGE
             REAL APTNUP, AREALF, CANNAA, CANWAA, CANWH, CARBO, GNUP, GPSM, GRNWT, GRORT, HI, HIP, &
                  PCNGRN, PCNL, PCNRT, PCNST, PCNVEG, PODNO, PConc_Root, PConc_Seed, &
                  PConc_Shel, PConc_Shut, PODWT, PSTRES2, PTF, RLWR, ROOTN, RTWT, &
                  RTWTO, SATFAC, SDWT, SEEDNO, SHELPC, SI1(6), SI3(6), SKERWT, SLA, STMWTO, STOVER,  &
                  STOVN, STOVWT, SUMP, SWFAC, TOPWT, TURFAC, VSTAGE, WTLF, WTNCAN, WTNLF, WTNSD, WTNST, &
                  WTNVEG, XGNP, XN, YIELD, KSTRES
              !----- In MZ_GROSUB.for
                  INTEGER CMAT, EMAT, ICOLD, ICSDUR
                  CHARACTER*1 ISWNIT, ISWPOT, ISWPHO, ISWWAT, IDETO 
                  !CHARACTER*92 FILECC
                  REAL PRFTC(4), RGFIL(4), PARSR, CO2X(10), CO2Y(10), FSLFW, FSLFN, FSLFP, &
                       SDSZ, RSGR, RSGRT, CARBOT, STMWTE, RTWTE, LFWTE, SEEDRVE, LEAFNOE, PLAE, TMNC, &
                       TANCE, RCNP, RANCE, CTCNP1, CTCNP2, PLIGLF, PLIGRT, ASGDD, BSGDD, &
                       BIOMAS, CANHT_POT, CUMDTTEG, CumLeafSenes, CumLeafSenesY, CumLfNSenes, CUMPH, &
                       EARWT, EP1, GRAINN, GRF, GNP, GROEAR, GROGRN, GROLF, GROSTM, PAR, LAI, &
                       LAIDOT, LFWT, LIFAC, MAXLAI, NPOOL, NPOOL1, NPOOL2, NSDR, NSINK,PCO2,PCNSD,PDWI, &
                       PLA, PLAG, PLAS, PRFT, RANC, RMNC, RNLAB, SEEDRV, SENLA, SLAN, SLFC, SLFN, SLFP, & 
                       SLFT, SLFW, SI2(6), SI4(6), SI5(6), SI6(6), Stem2Ear, STMWT, SUMRL, SUMEX, SWEXF, &
                       SWMAX, SWMIN, TANC, TAVGD, TCNP, TEMPM, TFAC, TI, TNLAB, TOTNUP, TRNU, TSS(20), &
                       VANC, VMNC, XANC, XLFWT, XNF, YIELDB, NDEF3, NFAC, PGRORT
                 !----- In MZ_GROSUB -> P_Ceres.for
                       CHARACTER*2 CROP
                       REAL SenSoilP, SenSurfP, CumSenSurfP, PestShut, PestRoot, PestShel, PestSeed, ShutMob, &
                            RootMob, ShelMob, PShut_kg, PShel_kg, PSeed_kg
                    !----- IN MZ_GROSUB -> P_Ceres -> RootSoilVol.for
                           REAL PlantPop, ROOTRAD
                           LOGICAL RootSoilVol_FIRST
                 !----- In MZ_GROSUB -> P_Ceres.for -> P_Plant
                       REAL PRoot_kg
                     !----- IN MZ_GROSUB -> P_Ceres -> IN P_Plant.for
                           REAL PConc_Shut_opt, PConc_Root_opt, PConc_Shel_opt, PConc_Seed_opt, &
                                PConc_Shut_min, PConc_Root_min, PConc_Shel_min, PConc_Seed_min, &
                                PConc_Plant, PPlant_kg, PSTRESS_RATIO, Shut_kg, Plant_kg
                     !----- P_IPPLNP
                           REAL, DIMENSION(3) :: N2Pmax, N2Pmin
                           REAL, DIMENSION(3) :: PCShutMin, PCLeafMin, PCStemMin, PCRootMin, PCShelMin, &
                                                 PCSeedMin, PCShutOpt, PCLeafOpt, PCStemOpt, PCRootOpt, &
                                                 PCShelOpt, PCSeedOpt
                           LOGICAL UseShoots
                           REAL :: FracPMobil, FracPUptake, SRATPHOTO, SRATPART
                    !----- P_Demand
                           REAL DeltPRoot, DeltPSeed, DeltPShel, DeltPShut, PRootDem, PSeedDem, PShelDem, &
                                PShutDem, PTotDem
                       !----- IN MZ_GROSUB -> P_Ceres -> P_Plant -> In PPlantSubs.for (P_Demand)
                           REAL PShutMobPool, PRootMobPool, PShelMobPool, PSeedMobPool
                    !----- P_Uptake
                           REAL N2P, PUptakeProf
                       !----- IN MZ_GROSUB -> P_Ceres -> P_Plant -> In P_Uptake.for
                                REAL WF_PHOS
           !-----------------------------------------------------------------------
           ! PLANT - MZ_CERES - PEST
             REAL SDNO(300), SHELN(300), WTSHE(300), SDDES(300), WTSD(300)   !300 is NCOHORTS in DSSAT
             REAL SWIDOT, WSHIDT, ASMDOT, DISLA, NPLTD, PPLTD, WLIDOT, WRIDOT, WSIDOT
              !----- In PEST.for
                  CHARACTER*5   PID(500), PCPID(500,6)   !500 is MAXPEST in DSSAT
                  CHARACTER*12  FILET
                  INTEGER TRTNUM, PCTID(500), IDAP(6,500), PCN, PNO(6)
                  REAL PHTHRS8, PDCF1(500,6), YPL(6,500)
              !----- PEST - PESTCP
                  REAL NSDDL, NSDDM, NSDDS, NSHDL, NSHDM, NSHDS, TLFAD, TLFMD, TRTLV, WRTMD, WSDDL, WSDDM, &
                       WSDDS, WSHDL, WSHDM, WSHDS, CPPLTD, PCLMA, PCLMT, PCSTMD, PDLA, PLFAD, PLFMD, PPSR, &
                       PRTLF, PRTLV, PRTMD, PSDDL, PSDDM, PSDDS, PSHDL, PSHDM, PSHDS, PSTMD, PVSTGD, TDLA, &
                       TPSR, TRTLF, VSTGD, WSTMD, WSDD, PSDD, PRLV
                 !----- In PEST -> PESTCP.for
                       REAL PL(6), LAIW, PLAIW, PDAM, FSM
              !----- PEST - ASMDM
                  REAL CASM
              !----- PEST - SEEDDM
                  REAL CSDM, CSDN, CSHM, CSHN, SDIDOT, SHIDOT, TSDNOL, TSDNOM, TSDNOS, TSDWTL, TSDWTM, & 
                       TSDWTS, TSHNOL, TSHNOM, TSHNOS, TSHWTL, TSHWTM, TSHWTS
                 !----- In PEST -> SEEDDM.for
                       REAL SDWDES(300),SDNDES(300),SHWDES(300),SHNDES(300),TSDNO,TSDWT,TSHNO,TSHWT  
              !----- PEST - VEGDM
                  REAL CLAI, CLFM, CSTEM, DISLAP
                 !----- In PEST -> VEGDM.for
                       REAL CLSEN, CLFRZ, LDAM, LAIDAM, SDAM
              !----- PEST - ROOTDM
                  REAL CRLF, CRLV, CRTM, RLFDOT, RLVDOT
                 !----- In PEST -> ROOTDM.for
                  REAL CUMDEP
           !-----------------------------------------------------------------------
           ! PLANT - MZ_CERES - MZ_ROOTGR
             REAL RTDEP
              !----- In MZ_ROOTS.for
                  REAL ESW(20),RLDF(20),RNLF,RNFAC,RLNEW
        !------------------------------------------------------------------------
        !---- ADDITIONAL CONTROL ------------------------------------------------
        !---- Pang 2023.09.19 ---------------------------------------------------
        !------------------------------------------------------------------------
          LOGICAL :: doseasinit
    end type dssat48dec
          !INTEGER, PARAMETER :: &
          !Dynamic variable values
          !  RUNINIT  = 1, &
          !  INIT     = 2, & !Will take the place of RUNINIT & SEASINIT
          !               !     (not fully implemented)
          !  SEASINIT = 2, &
          !  RATE     = 3, &
          !  EMERG    = 3, & !Used for some plant processes.  
          !  INTEGR   = 4, & 
          !  OUTPUT   = 5, & 
          !  SEASEND  = 6, &
          !  ENDRUN   = 7    
end module dssat48_module
