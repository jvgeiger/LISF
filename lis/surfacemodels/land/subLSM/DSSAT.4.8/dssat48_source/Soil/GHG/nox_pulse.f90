!=======================================================================
! DayCent nox_pulse, Subroutine
!=======================================================================
!
!/*              Copyright 1993 Colorado State University                    */
!/*                      All Rights Reserved                                 */
!
!/*****************************************************************************
!**
!**  FILE:      nox_pulse.c
!**
!**  FUNCTION:  float nox_pulse() 
!**
!**  PURPOSE:   
!**  
!**  INPUTS:
!**    ppt  - daily precip (cm)
!**    snow - snow cover (cm SWE)
!**
!**  GLOBAL VARIABLES:
!**    None
!**
!**  LOCAL VARIABLES:
!**    cumppt[] - circular array holding precipitation values
!**    ii       - loop control variable
!**    indx     - current array index
!**    mptr     - starting position in mtplr circular array 
!**    mtplr[]  - circular array ??
!**    nph      - starting position in ph circular array
!**    npl      - starting position in pl circular array
!**    npm      - starting position in pm circular array
!**    nppt     - starting position in cumppt circular array
!**    pflag    - ??
!**    ph[]     - circular array ??
!**    PHDAYS   - ?? (13)
!**    pl[]     - circular array ??
!**    PLDAYS   - ?? (2)
!**    pm[]     - circular array ??
!**    PMDAYS   - ?? (6)
!**    PPTDAYS  - number of consecutive days to track precipitation values (15)
!**    retval   - return value, increase of NO due to moisture and rain >= 1.0
!**    sumppt   - sum of precipitation values in cumppt array
!**
!**  OUTPUTS:
!**    retval - increase of NO due to moisture and rain >=1.0
!**
!**  CALLED BY:
!**    trace_gas_model()
!**
!**  CALLS:
!**    max3 - return the maximum of three input values
!**
!*****************************************************************************/
      Subroutine nox_pulse (dynamic, rain, snow, nox_puls, nest, t) !Pang 2024.03.04
      use ModuleDefs
      use dssat48_lsmMod
      save

      integer, intent(in) :: dynamic
      integer nest, t
      real, intent(in) :: rain, snow
      real, intent(out) :: nox_puls

      integer, parameter ::   &
          PPTDAYS = 15,       &
          PLDAYS = 2,         &
          PMDAYS = 6,         &
          PHDAYS = 13 
      
      real cumppt(0:PPTDAYS-1)    !0:14
      real pl(0:PLDAYS-1)         !0:1         
      real pm(0:PMDAYS-1)         !0:5
      real ph(0:PHDAYS-1)         !0:12
      real mtplr(0:PHDAYS-1)      !0:12

      integer npl, npm, nph, nppt, mptr
      integer pflag

      real sumppt
      integer ii, indx
!------ Obtain Vars from Memory -----------------------------------
!------ Pang 2024.03.04 -------------------------------------------
       cumppt = dssat48_struc(nest)%dssat48(t)%cumppt 
       pl = dssat48_struc(nest)%dssat48(t)%pl_nox
       pm = dssat48_struc(nest)%dssat48(t)%pm
       ph = dssat48_struc(nest)%dssat48(t)%ph
       mtplr = dssat48_struc(nest)%dssat48(t)%mtplr
       npl = dssat48_struc(nest)%dssat48(t)%npl
       npm = dssat48_struc(nest)%dssat48(t)%npm
       nph = dssat48_struc(nest)%dssat48(t)%nph
       nppt = dssat48_struc(nest)%dssat48(t)%nppt
       mptr = dssat48_struc(nest)%dssat48(t)%mptr
       pflag = dssat48_struc(nest)%dssat48(t)%pflag

! =================================================================
      select case (dynamic)
! =================================================================
      case (runinit, seasinit)
! =================================================================

      cumppt = 0.0
      pl     = 1.0
      pm     = 1.0
      ph     = 1.0
      mtplr  = 1.0

      npl   = 0
      npm   = 0
      nph   = 0
      nppt  = 0
      mptr  = 0
      pflag = 0

      nox_puls = 1.0

!     write(334,'(" nppt   rain cumppt sumppt  npl     pl  npm     pm  nph     ph mptr  mtplr pflag nox_puls")')
!     write(334,'(i5,3f7.2,4(i5,f7.3),i6, f8.3)')  &
!        nppt, rain, cumppt(nppt), sumppt, npl, pl(npl), npm, pm(npm),   &
!        nph, ph(nph), mptr, mtplr(mptr), pflag, nox_puls

! =================================================================
      case (rate)
! =================================================================
      sumppt = 0.0
      cumppt(nppt) = rain
      do ii=1, PPTDAYS-1
        indx = Mod(nppt+ii,PPTDAYS)
        sumppt = sumppt + cumppt(indx)
      enddo

      if (snow > 0.0) then
        mtplr(mptr) = 0.0

      else if ((sumppt <= 1.0) .and. (rain > 0.1)) then
 !    /* initiate new pulse */
        if (rain < 0.5) then
          do ii=0, PLDAYS-1
            indx = MOD(npl+ii,PLDAYS)
            pl(indx) = 2.8 * exp(-0.805 * (ii+1))
          enddo
          pflag = 2

        else if ((rain >= 0.5) .and. (rain <= 1.5)) then
          do ii=0, PMDAYS-1
            indx = MOD(npm+ii,PMDAYS)
            pm(indx) = 3.67 * exp(-0.384 * (ii+1))
          enddo
          pflag = 6

        else
          do ii=0, PHDAYS-1
            indx = MOD(nph+ii,PHDAYS)
            ph(indx) = 4.615 * exp(-0.208 * (ii+1))
          enddo
          pflag = 13

        endif
        
        mtplr(mptr) = amax1(pl(npl), pm(npm), ph(nph))
        pflag = pflag - 1

      else if (pflag > 0) then
        mtplr(mptr) = amax1(pl(npl), pm(npm), ph(nph))
        pflag = pflag - 1

      else
        mtplr(mptr) = 1.0
      endif

      nox_puls = mtplr(mptr)

!     write(334,'(i5,3f7.2,4(i5,f7.3),i6,f8.3)')  &
!        nppt, rain, cumppt(nppt), sumppt, npl, pl(npl), npm, pm(npm),   &
!        nph, ph(nph), mptr, mtplr(mptr), pflag, nox_puls

      pl(npl) = 1.0
      pm(npm) = 1.0
      ph(nph) = 1.0
      
!     /* increment pointers in circular arrays */
      npl = MOD(npl+1, PLDAYS)
      npm = MOD(npm+1, PMDAYS)
      nph = MOD(nph+1, PHDAYS)
      nppt = MOD(nppt+1,PPTDAYS)
      mptr = MOD(mptr+1,PHDAYS)

! =================================================================
      end select
! =================================================================
!------ Save Vars to Memory -----------------------------------
!------ Pang 2024.03.04 -------------------------------------------
       dssat48_struc(nest)%dssat48(t)%cumppt = cumppt
       dssat48_struc(nest)%dssat48(t)%pl_nox = pl
       dssat48_struc(nest)%dssat48(t)%pm = pm
       dssat48_struc(nest)%dssat48(t)%ph = ph
       dssat48_struc(nest)%dssat48(t)%mtplr = mtplr
       dssat48_struc(nest)%dssat48(t)%npl = npl
       dssat48_struc(nest)%dssat48(t)%npm = npm
       dssat48_struc(nest)%dssat48(t)%nph = nph
       dssat48_struc(nest)%dssat48(t)%nppt = nppt
       dssat48_struc(nest)%dssat48(t)%mptr = mptr
       dssat48_struc(nest)%dssat48(t)%pflag = pflag
      return
      end subroutine nox_pulse
!=======================================================================

