! Moved from GHG_mod.for to share with dssat48_lsmMod.F90.
      MODULE GHG_types_mod
      USE ModuleDefs, only : NL
      TYPE N2O_type
!            Daily        Cumul        Layer         
        REAL TNOXD,       CNOX,        DENITRIF(NL)  ![N] Denitrified
!       REAL TN2OdenitD,  CN2Odenit,   N2Odenit(NL)  !N2O[N] from denit
        REAL TN2OdenitD,  CN2Odenit,   N2Odenit(NL)  
        REAL TN2OnitrifD, CN2Onitrif,  N2ONitrif(NL) !N2O[N] from nitr

        REAL TN2D,        CN2,         N2flux(NL)    !N2[N] from denit
!                                      N2Oflux = N2Odenit + N2ONitrif
        REAL                           N2OFLUX(NL)   
        REAL TNOfluxD,    CNOflux,     NOflux(NL)    !NO flux from nitr

        REAL TNITRIFY,    CNITRIFY,    NITRIF(NL)    ![N] Nitrified 

        REAL N2_emitted,  CN2_emitted                !N2[N] emitted
        REAL N2O_emitted, CN2O_emitted               !N2O[N] emitted
        REAL NO_emitted,  CNO_emitted                !NO[N] emitted

        REAL TNGSoil  !N2, N2O, and NO in soil

        REAL, DIMENSION(NL) :: WFPS
      END TYPE N2O_type

      TYPE CH4_type
        REAL CH4Consumption, CH4Emission, CH4Leaching, CH4Stored
        REAL CumCH4Consumpt, CumCH4Emission, CumCH4Leaching
        REAL CO2emission, CumCO2Emission                      
      END TYPE CH4_type
      END MODULE GHG_types_mod
