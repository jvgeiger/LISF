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
         INTEGER            :: NREPS
         CHARACTER*1        :: RNMODE
         !CHARACTER*80       :: PATHEX
        !------------------------------------------------------------------------
        !  CSM
        !------------------------------------------------------------------------
           LOGICAL :: DONE
           INTEGER :: YRDOY_END
        !------------------------------------------------------------------------
        !  LAND
        !------------------------------------------------------------------------
         INTEGER            :: YRPLT, MDATE, YREND        
        !------------------------------------------------------------------------
        ! WEATHER
        ! Pang 2024.03.14
        !------------------------------------------------------------------------
          INTEGER NEV
        !------------------------------------------------------------------------
        ! WEATHER - CO2VAL
        !------------------------------------------------------------------------
          INTEGER NVals, INDX
          REAL CO2BAS
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
          REAL :: WB_CRAIN,TDRAIN, TRUNOF, TSWINI, WB_TSW !(WBSUM)
          REAL :: WTDEP, TDFD, WATAVL
          REAL, DIMENSION(20) :: DLAYR_YEST, SWDELTT, SWDELTL
           !----------------------------------------------------------------------
           !  SOIL - WATBAL -> MULCHWATER.for ------------------------------------
           !----------------------------------------------------------------------
             REAL CNRAIN, NRAIN, MULWATADD, RESWATADD, FLOODY, SNOWY, MWY
           !----------------------------------------------------------------------
           !  SOIL - WATBAL -> WBAL.for ------------------------------------------
           !----------------------------------------------------------------------
             REAL CUMRESWATADD,  CUMWBAL, CUMMULEVAP, TSWY
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
          REAL, DIMENSION(3) :: CEDAM, CES21T, CES23T, CES3M, CES3T, CES3X, CESTR, FRDAE
          REAL, DIMENSION(3) :: TLITE, TMETABE, TSOM1E, TSOM2E, TSOM23E, TSOM3E, TSOME, TSTRUCE
          REAL, DIMENSION(0:0,3) :: CES21I
          REAL, DIMENSION(0:1,3) :: CES1M, CES1T, CES1X, CES21M, CES21S
          REAL, DIMENSION(0:1,3) :: CES21X, CES23LM
          REAL, DIMENSION(0:1,3) :: CES23LX, CES23M, CES23X
          REAL, DIMENSION(0:1,3) :: CES2LI, CES2LM, CES2LS, CES2LX
          REAL, DIMENSION(0:1) :: CO2MET, DECMET, DECS1, DECSTR, LIGSTR
          REAL :: DECS2(1), DECS3(1), METABE(0:20,3)
          REAL :: SOM1E(0:20,3),SOM23E(20,3), SOM2E(20,3), SOM3E(20,3), AMINRL(20,3), CES2(20,3), CES3(20,3)
          REAL :: STRUCE(0:20,3), ACCMNR(0:20,3), CES1(0:20,3), LITE(0:20,3)
          REAL :: CO2S2, CO2S3, CULMETQ, CULS1Q, CULS2Q, CULS3Q, CULSTRQ
          REAL :: FRMETI, FRMETS, RESDAX, SENESSUMC, SENESSUMN, SENESSUMP
          REAL, DIMENSION(0:20) :: FRLSTR, CO2S1, LIGC, METABC, SOM1C
          REAL, DIMENSION(0:20) :: STRUCC
          REAL, DIMENSION(20) :: S1S3, S2S3, SOM23C, SOM2C, SOM3C, TXS1
          INTEGER DISTURBENDMAX, DISTURBNUM, DNUM
          REAL :: DISTURBDEPTHMAX, TLITC, TMETABC, TSOM1C, TSOM2C, TSOM3C, TSOMC, TSTRUCC
          INTEGER, DIMENSION(9000*3) :: DISTURBEND, DEND  !Do we really need 9000*3?
          REAL, DIMENSION(9000*3) :: DISTURBDEPTH, DDEP
          REAL, DIMENSION(0:1,2) :: CO2STR(0:1,2)
          LOGICAL :: DOCULT(0:20), ADDMETABEFLAG, FRMETFLAG
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
           INTEGER :: IUOF, LFD10
           REAL :: TNOM, CMINERN, CIMMOBN, WTNUP, CumSumFert
           REAL :: CLeach, CNTILEDR, CNETMINRN, CNUPTAKE !(NFLUX)
           REAL :: ALI, nox_puls, TNH4, TNO3, TUREA, TN2OnitrifD, TNH4NO3
           !REAL :: CNITRIFY, CN2Onitrif, CNOflux !(N2O_data)
           REAL :: TFNITY(20), SNH4(20), SNO3(20), UREA(20), dD0(20), DENITRIF(20)
           REAL :: ALGFIX, CUMFNRO, TOTAML,BD1,TOTFLOODN !(FLOOD_CHEM)
           REAL :: CNOX, TNOXD !(Denit)
           !------------------------------------------------------------------------
           !  SOIL - SoilNi - FLOOD_CHEM
           !  Pang 2024.03.01
           !------------------------------------------------------------------------
              REAL :: FNI, BD2, OXFAC
              REAL :: CUMPERCN, CUMFNU, ALI_CHEM
              LOGICAL AKILL, IMMOBIL
              !------------------------------------------------------------------------
              !  SOIL - SoilNi - FLOOD_CHEM - FCHEM
              !  Pang 2024.03.04
              !------------------------------------------------------------------------
                 REAL :: FPI, YALGA, ALGACT_FCHEM, FTI, FPH
           !------------------------------------------------------------------------
           !  SOIL - SoilNi - OXLAYER
           !  Pang 2024.03.04
           !------------------------------------------------------------------------
              REAL :: OXNI, ALGACT_OXLAYR, OXPH
           !------------------------------------------------------------------------
           !  SOIL - SoilNi - Denit_DayCent
           !  Pang 2024.03.04
           !------------------------------------------------------------------------
              REAL :: A(4), min_nitrate
              INTEGER :: DD_layer
           !------------------------------------------------------------------------
           !  SOIL - SoilNi - Denit_Ceres
           !  Pang 2024.03.04
           !------------------------------------------------------------------------
              REAL :: DLAG(20)
           !------------------------------------------------------------------------
           !  SOIL - SoilNi - N2Oemit
           !  Pang 2024.03.04
           !------------------------------------------------------------------------
              REAL :: n2o_soil(20), no_soil(20), n2_soil(20)
           !------------------------------------------------------------------------
           !  SOIL - SoilNi - nox_pulse
           !  Pang 2024.03.04
           !------------------------------------------------------------------------
              REAL cumppt(0:15-1), pl_nox(0:2-1), pm(0:6-1), ph(0:13-1), mtplr(0:13-1)
              INTEGER npl, npm, nph, nppt, mptr, pflag
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
           REAL, DIMENSION(20) :: FracNoRts, PiLabRts, PiSolRts, PiActive
           REAL, DIMENSION(20) :: PiLabile, PiStable, SPi_Total
           REAL :: TMINERP, TIMMOBP, SPiAvlProf, SPiLabProf, SPiActProf, SPiStaProf, SPiSolProf
           REAL :: SPiTotProf, SPiLabRtsProf, SPiLabNoRtsProf, SPiSolRtsProf, SPiSolNoRtsProf
        !------------------------------------------------------------------------
        !  SOIL - SoilKi
        !  Pang 2023.10.16
        !------------------------------------------------------------------------
           REAL, DIMENSION(20) :: SKi_Tot, Ki_Avail, DLTSKiAvail
           REAL :: SKiAvlProf, SKiTotProf
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
            !  SPAM -> In STEMP_EPIC.for (Pang 2024.02.12)
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
        !  In MGMTOPS
        !  Pang 2024.03.06
        !-------------------------------------------------------------------------
          INTEGER HDATE(3), NHAR, TILLNO, NAP, TIL_IRR, TOTIR, NCHEM
        !-------------------------------------------------------------------------
           !-------------------------------------------------------------------------
           !  In MGMTOPS -> AUTPLT
           !  Pang 2024.03.08
           !-------------------------------------------------------------------------
             CHARACTER*1  PLME
             CHARACTER*2  CROP_MGM
             REAL PTX, PTTN, SWPLTL, SWPLTH, SWPLTD
             INTEGER PWDINF, PWDINL
           !-------------------------------------------------------------------------
           !  In MGMTOPS -> AUTHAR
           !  Pang 2024.03.08
           !-------------------------------------------------------------------------
             INTEGER HDLAY, HLATE, HSTG(3)
           !-------------------------------------------------------------------------
           !  In MGMTOPS -> Fert_Place
           !  Pang 2024.03.11
           !-------------------------------------------------------------------------
             CHARACTER*5  FERMET(9000), FERTYPE_CDE(9000) !NAPPL=9000
             INTEGER NFERT, NSR
             REAL AMTFER(3)    !NELEM = 3
             INTEGER FDAY(9000), FERTYP(9000), FERDEP(9000), ANFER(9000), APFER(9000)
             INTEGER AKFER(9000), NAPFER(3)
           !-------------------------------------------------------------------------
           !  In MGMTOPS -> OM_Place
           !  Pang 2024.03.11
           !-------------------------------------------------------------------------
             CHARACTER*5  RESTYP(9000)
             INTEGER NRESDL, NRESAP, RESDAY(9000), NAPRes
             REAL RESLIGNIN(9000), DRESMG, RESIDUE(9000), RESN(9000), RESP(9000), RIP(9000), RESDEP(9000)
           !-------------------------------------------------------------------------
           !  In MGMTOPS -> Tillage
           !  Pang 2024.03.11
           !-------------------------------------------------------------------------
             INTEGER NTIL, TILLDATE(9000), NLYRS(9000)
             REAL, DIMENSION(9000) :: CNP, RINP, MIXT, TDEP
             REAL, DIMENSION(9000, 20) :: DEP, BDP, SWCNP
           !-------------------------------------------------------------------------
           !  In MGMTOPS -> IRRIG
           !  Pang 2024.03.11
           !-------------------------------------------------------------------------
             LOGICAL SeasonalWL
             INTEGER NAPW, NTBL, NCOND, NPERC, NPUD, DaysSinceIrrig, AIRRCOD, IRON(20)
             INTEGER IFREQ(20), NGSIrrigs, IRINC, NDAYS_DRY
             REAL TOTEFFIRR, GSWatUsed, ACCUM_ET, DSOIL, THETAC, THETAU, AIRAMT, EFFIRR
             REAL IMDEP(20), ITHRL(20), ITHRU(20), IRAMT(20), IREFF(20), AVWATI(20), AVWAT
             INTEGER, DIMENSION(9000) :: IRRCOD, CONDAT, IBDAT, IIRRCV, IPDAT, JULAPL, JULWTB, PUDDAT
             REAL, DIMENSION(9000) :: AMIR, BUND, COND, IPERC, PWAT
              !-------------------------------------------------------------------------
              !  In MGMTOPS -> IRRIG -> Flood_Irrig
              !  Pang 2024.03.11
              !-------------------------------------------------------------------------
                 LOGICAL CONVERTED
                 REAL APWAT, PERMW
           !-------------------------------------------------------------------------
           !  In MGMTOPS -> PADDY_MGMT
           !  Pang 2024.03.12
           !-------------------------------------------------------------------------
             INTEGER NDRY
           !-------------------------------------------------------------------------
           !  In MGMTOPS -> OpMgmt
           !  Pang 2024.03.12
           !-------------------------------------------------------------------------
             LOGICAL DPRINT
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
             REAL CMINEP, CNODMN, NAVL, PGAVL, RO, RP, RPROAV, RSPNO3, RSPNH4 !SWFAC, TURFAC, ,KSTRES
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
           !---- CROPGRO - PHENOL -----------------------------------------------
           !-----------------------------------------------------------------------
              REAL    :: DRPP, DTX, DXR57, FRACDN, PHTHRS(20), RVSTGE
              REAL    :: TDUMX, TDUMX2 !VSTAGE SeedFrac, VegFrac
              INTEGER :: NDLEAF, NDSET, NR1, NR2, NR5, NR7, NVEG0 !RSTAGE
              INTEGER :: YRNR1, YRNR2, YRNR3, YRNR5, YRNR7 !YREMRG
              !----- In PHENOL.for -------------------------------------------------
                  !CHARACTER*1 PLME
                  CHARACTER*3 CTMP(20), DLTYP(20)
                  REAL :: ATEMP, CLDVAR, CLDVRR, CSDVAR, CSDVRR, EVMODC, NSENP(20), OPTBI, PSENP(20)
                  REAL :: SDAGE, SLOBI, THVAR, TRIFOL, TB(5), TO1(5), TO2(5), TM(5), WSENP(20), MNEMV1, MNFLLL
                  REAL :: TNTFAC, TNTFC2, FNSTR(20), FPSTR(20), FSW(20), FT(20), FUDAY(20), PHZACC(20) !SDEPTH
                  INTEGER NPRIOR(20), TSELC(20), JPEND, NDVST, NVALPH(20), NVEG1
                 !----- In PHENOL.for -> RSTAGES.for
                    INTEGER NR0, NR3
                    REAL VegTime, PROG(20), PHTEM
                 !----- In PHENOL.for -> VSTAGES.for
                    REAL VSTGED, VSTAGP
           !-----------------------------------------------------------------------
           !---- CROPGRO - DEMAND -----------------------------------------------
           !-----------------------------------------------------------------------
               REAL AGRSD1, AGRSD2, AGRVG, AGRVG2, CDMREP, F, FNINL, FNINR, FNINS, FNINSD, FRLF, FRRT, FRSTM, &
                    GDMSD, GRRAT1, NDMNEW, NDMOLD, NDMREP, NDMSDR, NDMTOT, NDMVEG, NMINEP, NMOBR, PHTIM(300), &
                    PNTIM(300), POTCAR, POTLIP, SDGR, TURADD, XFRT
              !----- In DEMAND.for ------------------------------------------------
                   CHARACTER*3 TYPSDT
                   REAL CARMIN, FINREF, FNSDT(4), FRLFF, FRLFMX, FRSTMF, LIPOPT, LIPTB, LNGSH, NMOBMX, &
                      NRCVR, NVSMOB, PROLFF, PRORTF, PROSTF, & !PROSTI, PLIGSD, PMINSD, POASD, PROLFI, PRORTI
                      SDLIP, SDPRO, SHLAG, SLAMAX, SLAMIN, & !RCH2O, RLIG, RLIP, RMIN, RNO3C, ROA, RPRO
                      SLAPAR, SLAREF, SLAVAR, SLOSUM, SIZELF, SIZREF, SRMAX, THRESH, TURSLA, VSSINK, XFRMAX, &
                      XFRUIT, XLEAF(25), XSLATM(10), XTRFAC(10), XVGROW(6), XXFTEM(10), &
                      YLEAF(25), YSLATM(10), YSTEM(25), YTRFAC(10), YVREF(6), YXFTEM(10)
                   REAL CDMSDR, GDMSDR, NDMSD, NDMSH, RPRPUN, TMPFAC, DUMFAC, FVEG, SLAMN, SLAMX, GROMAX, SIZRAT, YVGROW(6)
                   REAL ADDSHL, TURXFR, CDMSD
                   INTEGER NAGE
           !-----------------------------------------------------------------------
           !---- CROPGRO - INCOMP -----------------------------------------------
           !-----------------------------------------------------------------------
               REAL AGRLF, AGRNOD, AGRRT, AGRSH1, AGRSH2, AGRSTM, SDPROR
              !----- In INCOMP.for ------------------------------------------------
                   !REAL SDPRO, SDLIP, RNO3C, RCH2O, RLIP, RLIG, ROA, RMIN, PROLFI, PROSTI, PRORTI, PROSHI, 
                        !PCARSH, PLIGSH, PLIPSH,PLIGSD,POASH, POASD, PMINSH, PMINSD, PLIGLF,  PLIGRT
                   REAL SDPROS, PCARLF, PCARST, PCARRT, PCARSD, PCARNO, PLIPLF, PLIPST, PLIPRT, PLIPNO, &
                        PLIGST, PLIGNO, POALF, POAST, POART, POANO, PMINLF, PMINST, PMINRT, PMINNO 
           !-----------------------------------------------------------------------
           !---- CROPGRO - NUPTAKE -----------------------------------------------
           !-----------------------------------------------------------------------
               REAL TRNH4U, TRNO3U !TRNU
              !----- In NUPTAK.for ------------------------------------------------
                 REAL RTNO3, RTNH4
           !-----------------------------------------------------------------------
           !---- CROPGRO - NFIX ------------------------------------------------
           !-----------------------------------------------------------------------
               REAL CNOD, DWNOD, DWNODA, NDTH, NFIXN, NODGR, WTNFX, SENNOD(20)
              !----- In NFIX.for --------------------------------------------------
                  CHARACTER*3 TYPFXT,TYPNGT,TYPFXD,TYPFXW,TYPFXA
                  REAL EFINOC, EFNFIX, PRONOD, SNACTM, NODRGM, DWNODI, NDTHMX, CNODCR, FNFXT(4), FNNGT(4), FNFXD(4), &
                       FNFXW(4), FNFXA(4) !RFIXN
                  INTEGER SDWNOD
           !-----------------------------------------------------------------------
           !---- CROPGRO - PODS ------------------------------------------------
           !-----------------------------------------------------------------------
               REAL AGRSD3, LAGSD, LNGPEG, NGRSD, NGRSH, PCTMAT, SDVAR, SHVAR, WSDDTN, & !PODNO, SEEDNO
                    WSHDTN, WTABRT, WTSHMT, FLWN(300)
              !----- In PODS.for -------------------------------------------------- 
                 CHARACTER*3   TYPPDT
                 INTEGER NR2TIM, TRIGGR
                 REAL  SDPDVR, PODUR, SETMAX, RFLWAB, XMPAGE, DSWBAR, & !THRESH, PROLFF, PROSHI, SHLAG,
                       FNPDT(4), XSWBAR(10), YSWBAR(10), XSWFAC(10), YSWFAC(10) !LNGSH
                 REAL MNESPM !LNGPEG, LAGSD, SDVAR, SHVAR
                 REAL FNINSH, NAVPOD, PGNPOD, WTSHM, PGAVLR !RPRPUN
                 REAL ACCAGE, AFLW, CNSTRES, CPSTRES, FLWRDY, PODADD, SHMINE, TEMPOD, SUPDE(300), AVTEM(300)
                 REAL ANINSD, CUMSIG, RSD !PODCOMP
                 !----- In PODS.for->PODCOM ----------------------------------------------------
                    REAL PROMIN, PROMAX, THETA !RCH2O, RLIP, RLIG, ROA, RMIN, PLIGSD, POASD, PMINSD
                    REAL RATION, RATIOC
                 !----- In PODS.for->FreshWt.for ----------------------------------------------------
                    REAL AvgDMC, AvgDPW, AvgFPW, PodAge !SHELPC
           !-----------------------------------------------------------------------
           !---- CROPGRO - PODDET ------------------------------------------------
           !-----------------------------------------------------------------------
              REAL PODWTD
              !----- In PODDET.for -------------------------------------------------- 
                 REAL DWC, PR1DET, PR2DET, XP1DET, XP2DET, DTC(300), DAYS(300), WPODY(300)
                 !TB(5), TO1(5), TO2(5), TM(5)
           !-----------------------------------------------------------------------
           !---- CROPGRO - VEGGR ------------------------------------------------
           !----------------------------------------------------------------------- 
              REAL CADLF, CADST, CMINEA, CRUSLF, CRUSRT, CRUSSH, CRUSST, EXCESS, NADLF, &
                   NADRT, NADST, NGRLF, NGRRT, NGRST, TNLEAK, WLDOTN, WRDOTN, WSDOTN !CANWH
              !----- In VEGGR.for--------------------------------------------------
                 REAL PROLFG, PROSTG, PRORTG, CADSTF, ATOP !PROLFI, PROSTI, PRORTI, CMOBMX
                 REAL CUMTUR, FNINLG, FNINRG, FNINSG, PGLEFT, SUPPN, VGRDEM
                 !----- In VEGGR.for -> CANOPY.for ----------------------------------
                    REAL XVSHT(15),YVSHT(15),YVSWH(15),XHWTEM(10),YHWTEM(10),XHWPAR(10),YHWPAR(10), RWIDTH, RHGHT
           !-----------------------------------------------------------------------
           !---- CROPGRO - SENES ------------------------------------------------
           !-----------------------------------------------------------------------
              REAL SLDOT, SLNDOT, SSDOT, SSNDOT
              !----- In SENES.for ------------------------------------------------
                 REAL PORPT, SENRTE, SENRT2, SENDAY, ICMP, TCMP, SENMAX(4), SENPOR(4), XSENMX(4), XSTAGE_CROPGRO(4), &
                      RATTP, SWFCAB(5)
           !-----------------------------------------------------------------------
           !---- CROPGRO - ROOTS ------------------------------------------------
           !-----------------------------------------------------------------------
              REAL SENRT(20), SRDOT !SATFAC, RTDEP
              !----- In ROOTS.for -------------------------------------------------
                 REAL RFAC1, RLDSM, RTDEPI, RTEXF, RTSEN, RTSDF, RTWTMIN, XRTFAC(4), YRTFAC(4) 
                 REAL DEPMAX, RFAC2, RFAC3, CumRootMass, TRLV !SUMRL, SUMEX
           !-----------------------------------------------------------------------
           !---- CROPGRO - GROW ------------------------------------------------
           !-----------------------------------------------------------------------
              REAL WLFDOT, CLW, CSW, GROWTH, GRWRES, LAIMX, PCCSD, & !AREALF, CANNAA, CANWAA, BETN
                   PCLSD, PCNSH, & ! PCNSD, PLIGLF, PLIGRT, PLIGNO, PLIGST, PLTPOP, PCNL, PCNRT, PCNST
                   PUNCSD, PUNCTR, RHOL, RHOS, RNITP, SDNPL, SDRATE, SEEDNI, & !ROWSPC, POTWT, RTWT
                   SHELWT, SLAAD, TOTWT, WCRLF, WCRRT, WCRSH, WCRST, WNRLF, WNRRT, & !SLA, TOPWT, STMWT
                   WNRSH, WNRST, WTCO, WTLO, WTMAIN, WTNEW, WTNLA, WTNLO, WTNNA, & !WTNCAN, WTNLF, WTLF
                   WTNNAG, WTNNO, WTNNOD, WTNOO, WTNRA, WTNRO, WTNRT, WTNSA, WTNSDA, WTNSDO, & !WTNSD
                   WTNSH, WTNSHA, WTNSHO, WTNSO, WTRO, WTSDO, WTSHO, WTSO, XPOD !WTNUP, WTNST
                   !ShutMob, RootMob, ShelMob
              !----- In GROW.for --------------------------------------------------
                 !CHARACTER*1  PLME
                 CHARACTER*2  XPODF 
                 REAL ALPHL, ALPHR, ALPHS, ALPHSH, PROSHF, & !PLIGST, PLIGLF, PLIGRT, PLIGSH, PLIPNO, 
                      !ROWSPC, PROLFF, PROSTF, PRORTF, PRONOD, PROLFI, PROSTI, PRORTI, PLTPOP
                      SDWTPL, WTFSD, WTPSD, CPFLF, CPFSTM, CPFRT, CPFNOD, CPFSH1, CPFSD1, PCNMIN, & !SDPRO
                      ALFDOT, AREAH, NLDOT, NSDOT, NRDOT, NSDDOT, NSHDOT, NTOVR, RHOR, RHOSH, SDWTAM, &
                      SenWt(0:20), SenLig(0:20), SenE(0:20,3), TGROW, WSDDOT, WSHDOT, WTCSD, WTLSD, &
                      WTNMOB, WTNTOT, SDPDOT, PUNDOT, NLPEST, &
                      WLDOT, WSDOT, WRDOT, WNDOT, WPDOT, WRCLDT, WRCSDT, WRCRDT, WRCSHD !WLFDOT
           !-----------------------------------------------------------------------
           !---- CROPGRO - P_CGRO ------------------------------------------------
           !-----------------------------------------------------------------------
              !REAL PConc_Shut, PConc_Root, PConc_Shel, PConc_Seed, PStres2
              !----- In P_CGRO.for ------------------------------------------------
                 REAL Leaf_kg, Stem_kg, Root_kg, Shel_kg, Seed_kg !SenSoilP, SenSurfP, PestShut, PestRoot, PestShel, PestSeed
           !-----------------------------------------------------------------------
           !---- CROPGRO - MOBIL ------------------------------------------------
           !-----------------------------------------------------------------------
              REAL NMINEA, NRUSLF, NRUSRT, NRUSSH, NRUSST
           !-----------------------------------------------------------------------
           !---- CROPGRO - OPGROW ------------------------------------------------
           !-----------------------------------------------------------------------
              !----- In Opgrow.for ------------------------------------------------
                 CHARACTER*12 OUTG, OUTPN, OUTPC
                 INTEGER N_LYR 
                 REAL CUMSENSURF, CUMSENSOIL, CUMSENSURFN, CUMSENSOILN, SWF_AV, TUR_AV, NST_AV, EXW_AV, PS1_AV, PS2_AV, KST_AV
                 LOGICAL FIRST_OPGROW
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
                INTEGER IDURP !07/10/2024
                REAL SIND, SUMDTT_2 !07/10/2024
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
                       EARWT, EP1, GRAINN, GRF, GNP, GROEAR, GROGRN, GROLF, GROSTM, IPAR, LAI, &
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
             INTEGER COUNT
        !------------------------------------------------------------------------
        ! JE Tight Coupling Variables
          REAL, DIMENSION(20):: LIS_sm            !Dimensioned the same as SW
          REAL, DIMENSION(20):: LIS_sm_start
          REAL LIS_lai
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
