!> \file
!> Central module for CTEM (biogeochemical)-related utilities
module ctemUtilities

  implicit none

  public :: genSortIndex
  public :: dayEndCTEMPreparation
  public :: accumulateForCTEM
  public :: ctemInit

contains

  ! --------------------------------------------------------------------------------------------------------------------
  !> \ingroup ctemutilities_genSortIndex
  
  !> Generate the sort index for correspondence between the CTEM pfts and the
  !! array of values in the parameter vectors (e.g. for 9 CTEM, the array is of
  !! size 12)
  !! @author V.Arora, J. Melton
  function genSortIndex ()

    use classicParams, only : ican, l2max, icc, nol2pfts

    implicit none

    integer :: genSortIndex(icc)

    integer :: icount, j, m, n

    icount = 0
    do j = 1,ican
      do m = 1,nol2pfts(j)
        n = (j - 1) * l2max + m
        icount = icount + 1
        genSortIndex(icount) = n
      end do ! loop 96
    end do ! loop 95

  end function genSortIndex
  
  ! --------------------------------------------------------------------------------------------------------------------
  !> \ingroup ctemutilities_genSortIndex
  
  !> Calculates the specific leaf area based on leaf life span.
  function calcSLA(peatlandType,sort)

    use classicParams, only : icc, lfespany, lfespany_peat, &
                              ilg, specsla, abszero

    implicit none

    real :: calcSLA
    character(8), intent(in) :: peatlandType !< Peatland type (else 'None')
    integer, intent(in) :: sort      !< index for correspondence between 9 pfts and the 12 values in parameters vectors

    !> The default is to use the specified SLA read in from the run parameters file. If desired that specsla
    !! can be set to 0.0 and the calculation below occurs. This calculation sets the SLA to be directly
    !! derived from the leaf life span.
    if (specsla(sort) > abszero) then
      calcSLA = specsla(sort)
    else 
      if (peatlandType == 'None') then ! Uplands
        calcSLA = 25.0 * (lfespany(sort) ** ( -0.50))
      else ! peatlands
        calcSLA = 25.0 * (lfespany_peat(sort) ** ( -0.50))
      end if 
    end if 

  end function calcSLA
  
  ! --------------------------------------------------------------------------------------------------------------------
  !> \ingroup ctemutilities_dayEndCTEMPreparation
  
  !> Prepare the CTEM input (physics) variables at the end of the day.
  !! @author V.Arora, J. Melton
  subroutine dayEndCTEMPreparation (vgat, ctem_tile, nml, nday)

    use classicParams, only : icc, ignd, ilg, nlat, nmos
    use ctemStateVars, only : veg_gat, ctem_tile_level

    implicit none

    integer, intent(in) :: nml       !< Counter representing number of mosaic tiles on modelled domain that are land
    integer, intent(in) :: nday      !< Number of short (physics) timesteps in one day. e.g., if physics timestep is 15 min this is 48.
    type(veg_gat(ilg, ignd)),         intent(inout) :: vgat
    type(ctem_tile_level(ilg, ignd)), intent(inout) :: ctem_tile
    
    integer :: i, j
    real :: fsstar_gat
    real :: flstar_gat

    associate( &
    altotcount_ctm    => vgat%altotcount_ctm,           & !< integer(:):! nlat   Counter used for calculating total albedo 
    anmossac_t        => ctem_tile%anmossac_t,          & !< real(:): 
    rmlmossac_t       => ctem_tile%rmlmossac_t,         & !< real(:): 
    gppmossac_t       => ctem_tile%gppmossac_t,         & !< real(:): 
    fsinacc_gat       => vgat%fsinacc_gat,              & !< real(:):!(ilg)     
    flutacc_gat       => vgat%flutacc_gat,              & !< real(:):!(ilg)     
    flinacc_gat       => vgat%flinacc_gat,              & !< real(:):!(ilg)     
    altotacc_gat      => vgat%altotacc_gat,             & !< real(:):!(ilg)    
    netrad_gat        => vgat%netrad_gat,               & !< real(:):!(ilg)      
    mossPresentgat    => vgat%mossPresent,             & !< character(8):(ilg) 
    tbaraccgat_t      => ctem_tile%tbaraccgat_t,        & !< real(:,:): 
    thliqacc_t        => ctem_tile%thliqacc_t,          & !< real(:,:): 
    thiceacc_t        => ctem_tile%thiceacc_t,          & !< real(:,:): 
    ancgvgac_t        => ctem_tile%ancgvgac_t,          & !< real(:,:): 
    rmlcgvga_t        => ctem_tile%rmlcgvga_t,          & !< real(:,:): 
    fsnowacc_t        => ctem_tile%fsnowacc_t,          & !< real(:): 
    taaccgat_t        => ctem_tile%taaccgat_t,          & !< real(:): 
    uvaccgat_t        => ctem_tile%uvaccgat_t,          & !< real(:): 
    vvaccgat_t        => ctem_tile%vvaccgat_t,          & !< real(:): 
    CFLUX_GAacc_t     => ctem_tile%CFLUX_GAacc_t,       & !< real(:): Daily accu. aerodynamic conductance (\f$[m day^{-1} ]\f$), inverse of boundary layer aerodynamic resistance (ra) 
    USTARBS_GAacc_t   => ctem_tile%USTARBS_GAacc_t,     & !< real(:): Daily accu. friction velocity to be used in nitrogen volatilization  \f$[m day^{-1} ]\f$ 
    ROFBacc_t         => ctem_tile%ROFBacc_t,           & !< real(:): Daily accu. base flow from bottom of soil column \f$[kg m^{-2} day^{-1} ]\f$ 
    NoonTempcrun          => ctem_tile%NoonTempcrun,          & !< real(:):
    NoonPreciprun        => ctem_tile%NoonPreciprun,      & !< real(:): 
    NoonPrecip        => ctem_tile%NoonPrecip,      & !< real(:): 
    QFCacc_t          => ctem_tile%QFCacc_t             & !< real(:,:): Daily accu. water removed from soil layers by transpiration \f$[kg m^{-2} day^{-1}]\f$ 
    )

    do i = 1,nml

      ! net radiation and precipitation estimates for ctem's bioclim

      uvaccgat_t(i) = uvaccgat_t(i)/real(nday)
      vvaccgat_t(i) = vvaccgat_t(i)/real(nday)
      fsinacc_gat(i) = fsinacc_gat(i)/real(nday)
      flinacc_gat(i) = flinacc_gat(i)/real(nday)
      flutacc_gat(i) = flutacc_gat(i)/real(nday)

      if (altotcount_ctm(i) > 0) then
        altotacc_gat(i) = altotacc_gat(i)/real(altotcount_ctm(i))
      else
        altotacc_gat(i) = 0.
      end if

      fsstar_gat = fsinacc_gat(i) * (1. - altotacc_gat(i))
      flstar_gat = flinacc_gat(i) - flutacc_gat(i)
      netrad_gat(i) = fsstar_gat + flstar_gat

      fsnowacc_t(i) = fsnowacc_t(i)/real(nday)
      taaccgat_t(i) = taaccgat_t(i)/real(nday)
      CFLUX_GAacc_t(i) = CFLUX_GAacc_t(i)/real(nday)
      USTARBS_GAacc_t(i) = USTARBS_GAacc_t(i)/real(nday)
      ROFBacc_t(i) = ROFBacc_t(i)/real(nday)

      do j = 1,ignd
        QFCacc_t(i,j) = QFCacc_t(i,j)/real(nday)
        tbaraccgat_t(i,j) = tbaraccgat_t(i,j)/real(nday)
        thliqacc_t(i,j) = thliqacc_t(i,j)/real(nday)
        thiceacc_t(i,j) = thiceacc_t(i,j)/real(nday)
      end do ! loop 831

      do j = 1,icc
        ancgvgac_t(i,j) = ancgvgac_t(i,j)/real(nday)
        rmlcgvga_t(i,j) = rmlcgvga_t(i,j)/real(nday)
      end do ! loop 832

      !     Daily average moss C fluxes-------------------\
      !     Capitulum biomass = 0.22 kg/m2 in hummock, 0.1 kg/m2 in lawn
      !     stem biomass = 1.65 kg/m2 in hummock, 0.77 kg/m2 in lawn (Bragazza et al.2004)
      !     the ratio between stem and capitulum = 7.5 and 7.7
      if (mossPresentgat(i) /= 'None') then                           
        anmossac_t(i) = anmossac_t(i)/real(nday)                  
        rmlmossac_t(i) = rmlmossac_t(i)/real(nday)
        gppmossac_t(i) = gppmossac_t(i)/real(nday)
      end if

      NoonTempcrun(i) = -9999
      NoonPrecip(i) = NoonPreciprun(i)      
      NoonPreciprun(i) = 0

    end do ! nml loop

    end associate
  end subroutine dayEndCTEMPreparation
  
  ! --------------------------------------------------------------------------------------------------------------------
  !> \ingroup ctemutilities_accumulateForCTEM
  
  !> Accumulate the CTEM input (physics) variables at the end of each physics timestep
  !! @author V.Arora, J. Melton
  subroutine accumulateForCTEM (class_gat, vgat, ctem_tile, nml, dofire)

    use classicParams,  only : icc, ignd, ilg, nlat, nmos, nbs, DELT, SBC
    use classStateVars, only : class_gather
    use ctemStateVars,  only : veg_gat, ctem_tile_level

    implicit none

    integer, intent(in) :: nml       !< Counter representing number of mosaic tiles on modelled domain that are land
    logical, intent(in) :: dofire

    type(class_gather(ilg, ignd, nbs)),         intent(inout) :: class_gat
    type(veg_gat(ilg, ignd)),                   intent(inout) :: vgat
    type(ctem_tile_level(ilg, ignd)),           intent(inout) :: ctem_tile
    
    integer :: i, j
    real :: fssgat
    real :: anmossgat, gppmossgat, rmlmossgat

    associate( &
    anmossac_t        => ctem_tile%anmossac_t,          & !< real(:): 
    rmlmossac_t       => ctem_tile%rmlmossac_t,         & !< real(:): 
    gppmossac_t       => ctem_tile%gppmossac_t,         & !< real(:): 
    altotcount_ctm    => vgat%altotcount_ctm,           & !< integer(:):! nlat   Counter used for calculating total albedo 
    fsinacc_gat       => vgat%fsinacc_gat,              & !< real(:):!(ilg)     
    flutacc_gat       => vgat%flutacc_gat,              & !< real(:):!(ilg)     
    flinacc_gat       => vgat%flinacc_gat,              & !< real(:):!(ilg)     
    altotacc_gat      => vgat%altotacc_gat,             & !< real(:):!(ilg)    
    netrad_gat        => vgat%netrad_gat,               & !< real(:):!(ilg)      
    preacc_gat        => vgat%preacc_gat,               & !< real(:):!(ilg)      
    mossPresentgat    => vgat%mossPresent,              & !< character(8)(ilg)
    tbaraccgat_t      => ctem_tile%tbaraccgat_t,        & !< real(:,:): 
    thliqacc_t        => ctem_tile%thliqacc_t,          & !< real(:,:): 
    thiceacc_t        => ctem_tile%thiceacc_t,          & !< real(:,:): 
    ancgvgac_t        => ctem_tile%ancgvgac_t,          & !< real(:,:): 
    rmlcgvga_t        => ctem_tile%rmlcgvga_t,          & !< real(:,:): 
    fsnowacc_t        => ctem_tile%fsnowacc_t,          & !< real(:): 
    taaccgat_t        => ctem_tile%taaccgat_t,          & !< real(:): 
    uvaccgat_t        => ctem_tile%uvaccgat_t,          & !< real(:): 
    vvaccgat_t        => ctem_tile%vvaccgat_t,          & !< real(:): 
    CFLUX_GAacc_t     => ctem_tile%CFLUX_GAacc_t,       & !< real(:): daily accu. boundary layer aerodynamic conductance \f$[m day^{-1} ]\f$ 
    USTARBS_GAacc_t   => ctem_tile%USTARBS_GAacc_t,     & !< real(:): daily accu. friction velocity \f$[m day^{-1} ]\f$ 
    ROFBacc_t         => ctem_tile%ROFBacc_t,           & !< real(:): daily avg. base flow from bottom of soil column \f$[kg m^{-2} s^{-1} ]\f$ 
    QFCacc_t          => ctem_tile%QFCacc_t,            & !< real(:,:): daily accu. water removed from soil layers by transpiration \f$[kg m^{-2} day^{-1}]\f$ 
    NoonTempc         => ctem_tile%NoonTempc,           & !< real(:): 
    NoonPreciprun     => ctem_tile%NoonPreciprun,       & !< real(:): 
    NoonWindtot       => ctem_tile%NoonWindtot,         & !< real(:): 
    NoonRhum          => ctem_tile%NoonRhum,            & !< real(:):
    NoonTempcrun      => ctem_tile%NoonTempcrun,        & !< real(:):
    ancsmoss          => vgat%ancsmoss,                 & !< real(:): 
    angsmoss          => vgat%angsmoss,                 & !< real(:): 
    ancmoss           => vgat%ancmoss,                  & !< real(:): 
    angmoss           => vgat%angmoss,                  & !< real(:): 
    rmlcsmoss         => vgat%rmlcsmoss,                & !< real(:): 
    rmlgsmoss         => vgat%rmlgsmoss,                & !< real(:): 
    rmlcmoss          => vgat%rmlcmoss,                 & !< real(:): 
    rmlgmoss          => vgat%rmlgmoss,                 & !< real(:): 
    FC                => class_gat%FC,                  & !< real(:): 
    FG                => class_gat%FG,                  & !< real(:): 
    FCS               => class_gat%FCS,                 & !< real(:): 
    FGS               => class_gat%FGS,                 & !< real(:): 
    FSIHGAT           => class_gat%FSIHGAT,             & !< real(:): Near-infrared radiation incident on horizontal surface \f$[W m^{-2} ]\f$ 
    FSVHGAT           => class_gat%FSVHGAT,             & !< real(:): Visible radiation incident on horizontal surface \f$[W m^{-2} ]\f$ 
    ALIRGAT           => class_gat%ALIRGAT,             & !< real(:): Diagnosed total near-infrared albedo of land surface [ ] 
    ALVSGAT           => class_gat%ALVSGAT,             & !< real(:): Diagnosed total visible albedo of land surface [ ] 
    GTGAT             => class_gat%GTGAT,               & !< real(:): Diagnosed effective surface black-body temperature [K] 
    FDLGAT            => class_gat%FDLGAT,              & !< real(:): Downwelling longwave radiation at bottom of atmosphere (i.e. incident on modelled land surface elements \f$[W m^{-2} ]\f$ 
    PREGAT            => class_gat%PREGAT,              & !< real(:): Surface precipitation rate \f$[kg m^{-2} s^{-1} ]\f$ 
    FSNOGAT           => class_gat%FSNOGAT,             & !< real(:): Diagnosed fractional snow coverage [ ] 
    TAGAT             => class_gat%TAGAT,               & !< real(:): Air temperature at reference height [K] 
    VLGAT             => class_gat%VLGAT,               & !< real(:): Meridional component of wind velocity \f$[m s^{-1} ]\f$ 
    ULGAT             => class_gat%ULGAT,               & !< real(:): Zonal component of wind velocity \f$[m s^{-1} ]\f$ 
    FSGGGAT           => class_gat%FSGGGAT,             & !< real(:): Diagnosed net shortwave radiation at soil surface \f$[W m^{-2} ]\f$ 
    TBARGAT           => class_gat%TBARGAT,             & !< real(r8), dimension(:,:): Temperature of soil layers [K] 
    FSGSGAT           => class_gat%FSGSGAT,             & !< real(:): Diagnosed net shortwave radiation at snow surface \f$[W m^{-2} ]\f$ 
    FSGVGAT           => class_gat%FSGVGAT,             & !< real(:): Diagnosed net shortwave radiation on vegetation canopy \f$[W m^{-2} ]\f$ 
    THICGAT           => class_gat%THICGAT,             & !< real(:,:): Volumetric frozen water content of soil layers \f$[m^3 m^{-3} ]\f$ 
    THLQGAT           => class_gat%THLQGAT,             & !< real(:,:): Volumetric liquid water content of soil layers \f$[m^3 m^{-3} ]\f$ 
    QFCGAT            => class_gat%QFCGAT,              & !< real(:,:): water removed from soil layers by transpiration \f$[kg m^{-2} s^{-1}]\f$ 
    QAGAT             => class_gat%QAGAT ,              &
    PRESGAT           => class_gat%PRESGAT,             & !< real(:,:):
    CFLUX_GAgat       => vgat%CFLUX_GA,                 & !< real(:): aerodynamic conductance (\f$[m day^{-1} ]\f$) 
    USTARBS_GA        => vgat%USTARBS_GA,               & !< real(:): grid average friction velocity to be used in nitrogen volatilization  \f$[m day^{-1} ]\f$ 
    ROFBGAT           => vgat%ROFB,                     & !< real(:): base flow from the bottom of soil column \f$[kg m^{-2} s^{-1} ]\f$ 
    ancsveggat        => vgat%ancsveg,                  & !< real(:,:): 
    ancgveggat        => vgat%ancgveg,                  & !< real(:,:): 
    rmlcsveggat       => vgat%rmlcsveg,                 & !< real(:,:): 
    rmlcgveggat       => vgat%rmlcgveg                  & !< real(:,:): 
    )

    do i = 1,nml

      fssgat = fsvhgat(i) + fsihgat(i)  ! total shortwave radiation
      fsinacc_gat(i) = fsinacc_gat(i) + fssgat 
      flinacc_gat(i) = flinacc_gat(i) + fdlgat(i)
      flutacc_gat(i) = flutacc_gat(i) + sbc * gtgat(i) ** 4
      preacc_gat(i) = preacc_gat(i) + pregat(i) * delt
      fsnowacc_t(i) = fsnowacc_t(i) + fsnogat(i)
      taaccgat_t(i) = taaccgat_t(i) + tagat(i)
      vvaccgat_t(i) = vvaccgat_t(i) + vlgat(i)
      uvaccgat_t(i) = uvaccgat_t(i) + ulgat(i)
      CFLUX_GAacc_t(i) = CFLUX_GAacc_t(i) + CFLUX_GAgat(i)
      USTARBS_GAacc_t(i) = USTARBS_GAacc_t(i) + USTARBS_GA(i)
      ROFBacc_t(i) = ROFBacc_t(i) + ROFBGAT(i)

      if (fssgat > 0.) then
        altotacc_gat(i) = altotacc_gat(i) + (fssgat - &
                          (FSGVGAT(I) + FSGSGAT(I) + FSGGGAT(I))) &
                          /fssgat
        altotcount_ctm = altotcount_ctm + 1
      end if

      do j = 1,ignd
        QFCacc_t(i,j) = QFCacc_t(i,j) + QFCGAT(i,j)
        tbaraccgat_t(i,j) = tbaraccgat_t(i,j) + tbargat(i,j)
        thliqacc_t(i,j) = thliqacc_t(i,j) + THLQGAT(i,j)
        thiceacc_t(i,j) = thiceacc_t(i,j) + THICGAT(i,j)
      end do ! loop 710

      do j = 1,icc
        ancgvgac_t(i,j) = ancgvgac_t(i,j) + (1. - fsnogat(i)) * ancgveggat(i,j) + fsnogat(i) * ancsveggat(i,j)
        rmlcgvga_t(i,j) = rmlcgvga_t(i,j) + (1. - fsnogat(i)) * rmlcgveggat(i,j) + fsnogat(i) * rmlcsveggat(i,j)
      end do ! loop 713

      !   Accumulate moss C fluxes to tile level then daily           
      if (mossPresentgat(i) /= 'None') then
        anmossgat = fcs(i) * ancsmoss(i) + fgs(i) * angsmoss(i) + fc(i) * ancmoss(i) + fg(i) * angmoss(i)
        rmlmossgat = fcs(i) * rmlcsmoss(i) + fgs(i) * rmlgsmoss(i) + fc(i) * rmlcmoss(i) + fg(i) * rmlgmoss(i)
        gppmossgat = anmossgat + rmlmossgat

        anmossac_t(i) = anmossac_t(i)   + anmossgat
        rmlmossac_t(i) = rmlmossac_t(i)  + rmlmossgat
        gppmossac_t(i) = gppmossac_t(i)  + gppmossgat

      end if

      !if fire is active update the noon climate variables for FWI calculations
      NoonPreciprun(i) = NoonPreciprun(i) + PREGAT(i) * delt
      NoonTempcrun(i) = max(TAGAT(i),NoonTempcrun(i))
      if (dofire .and. NoonTempcrun(i) == TAGAT(i)) then
        NoonTempc(i) = TAGAT(i) - 273.15
        NoonWindtot(i) = sqrt((ULGAT(i)**2)+(VLGAT(i)**2)) * 3.6
        NoonRhum(i) = ((QAGAT(i) * (PRESGAT(i)*0.01) / (0.378 * QAGAT(i) + 0.622))/(6.112 * exp((17.67 * NoonTempc(i))/(NoonTempc(i) + 243.5))))*100
        if (NoonRhum(i)> 100.0) NoonRhum(i)=100.0 !!! make sure RHUM <= 100.0
        if (NoonRhum(i)< 0.0) NoonRhum(i)=0.0 !!! make sure RHUM not < 0
      end if
    end do 

    end associate
  end subroutine accumulateForCTEM
  
  !
  ! --------------------------------------------------------------------------------------------------------------------
  !> \ingroup ctemutilities_ctemInit
  
  !> Find mosaic tile (grid) average vegetation biomass, litter mass, and soil c mass.
  !! Also initialize additional variables which are used by CTEM (biogeochemical processes).
  !! @author V.Arora, J. Melton, A. Asaadi
  subroutine ctemInit (class_rot, vrot, nltest, nmtest)

    use classicParams,  only : icc, ignd, iccp1, nlat, nmos, nbs
    use ctemStateVars,  only : veg_rot
    use classStateVars, only : class_rowtile
    use peatlandsMod,   only : peatStorage, peatDepth

    implicit none

    integer, intent(in) :: nltest
    integer, intent(in) :: nmtest
    type(class_rowtile(nlat, nmos, ignd, nbs)), intent(inout) :: class_rot
    type(veg_rot(nlat, nmos, ignd)),            intent(inout) :: vrot
    integer :: i, m, j, k

    associate( &
    peatlandTyperow   => vrot%peatlandType,             & !< character(8), dimension(:,:) :
    mossPresentrow    => vrot%mossPresent,              & !< character(8), dimension(:,:) :
    upMossSoilCrow    => vrot%upMossSoilC,              & !< real, dimension(:,:,:) : Moss soil C mass (non-peat), \f$kg C/m^2\f$ 
    slairow           => vrot%slai,                     & !< real, dimension(:,:,:) : 
    vgbiomasrow       => vrot%vgbiomas,                 & !< real, dimension(:,:) :
    gavglairow        => vrot%gavglai,                  & !< real, dimension(:,:) : 
    gavgltmsrow       => vrot%gavgltms,                 & !< real, dimension(:,:) : 
    gavgscmsrow       => vrot%gavgscms,                 & !< real, dimension(:,:) :
    lmaxtrow          => vrot%lmaxt,                    & !< real, dimension(:,:,:) : 
    smaxtrow          => vrot%smaxt,                    & !< real, dimension(:,:,:) : 
    rmaxtrow          => vrot%rmaxt,                    & !< real, dimension(:,:,:) : 
    fcanrot           => class_rot%fcanrot,             & !< real, dimension(:,:,:) : 
    gleafmasrow       => vrot%gleafmas,                 & !< real, dimension(:,:,:) :! 
    bleafmasrow       => vrot%bleafmas,                 & !< real, dimension(:,:,:) :! 
    stemmassrow       => vrot%stemmass,                 & !< real, dimension(:,:,:) :! 
    rootmassrow       => vrot%rootmass,                 & !< real, dimension(:,:,:) :! 
    litrmassrow       => vrot%litrmass,                 & !< real, dimension(:,:,:,:) : 
    soilcmasrow       => vrot%soilcmas,                 & !< real, dimension(:,:,:,:) : 
    fcancmxrow        => vrot%fcancmx,                  & !< real, dimension(:,:,:) : 
    peatdeprow        => vrot%peatdep,                  & !< real, dimension(:,:) : 
    litrmsmossrow     => vrot%litrmsmoss,               & !< real, dimension(:,:,:) : 
    Cmossmasrow       => vrot%Cmossmas,                 & !< real, dimension(:,:) : 
    peatSoilCrow      => vrot%peatSoilC                 & !< real, dimension(:,:) : Peat soil C pool (kgC/m2)
    )

    ! ------

    ! Initialize to zero:
    slairow(:,:,:) = 0.0     ! if bio2str is not called we need to initialize this to zero
    vgbiomasrow(:,:) = 0.0
    gavglairow(:,:) = 0.0
    gavgltmsrow(:,:) = 0.0
    gavgscmsrow(:,:) = 0.0
    lmaxtrow(:,:,:) = 0.0
    smaxtrow(:,:,:) = 0.0
    rmaxtrow(:,:,:) = 0.0

    do i = 1,nltest ! loop 115
      do m = 1,nmtest
        do j = 1,iccp1
          if (j < iccp1) then 
            vgbiomasrow(i,m) = vgbiomasrow(i,m) + fcancmxrow(i,m,j) * &
                               (gleafmasrow(i,m,j) + stemmassrow(i,m,j) + &
                               rootmassrow(i,m,j) + bleafmasrow(i,m,j))
             
            do k = 1,ignd
              gavgltmsrow(i,m)=gavgltmsrow(i,m)+fcancmxrow(i,m,j)* &
                  &                       litrmassrow(i,m,j,k)
              gavgscmsrow(i,m)=gavgscmsrow(i,m)+fcancmxrow(i,m,j)* &
                  &         soilcmasrow(i,m,j,k)
            end do ! ignd
                      
          else ! bareground or consider peatland pools  
                                
            if (peatlandTyperow(i,m) == 'None') then ! NON-peatland tile (upland)           
              ! Only add the litter and soil C since this is bareground.
              do k = 1,ignd
                gavgltmsrow(i,m)=gavgltmsrow(i,m)+ (1.0-sum(fcanrot(i,m,:)))*litrmassrow(i,m,j,k)
                gavgscmsrow(i,m)=gavgscmsrow(i,m)+ (1.0-sum(fcanrot(i,m,:)))*soilcmasrow(i,m,j,k)
              end do
            else ! peatlands
                ! Find the peat depth
                peatdeprow(i,m) = peatDepth(peatSoilCrow(i,m))
                
                gavgscmsrow(i,m) = gavgscmsrow(i,m) + peatSoilCrow(i,m)

            end if 
            
            if (mossPresentrow(i,m) /= 'None') then
              ! Add the moss litter mass.
              do k = 1, ignd
                gavgltmsrow(i,m) = gavgltmsrow(i,m) + litrmsmossrow(i,m,k)

                ! Add the moss derived upland soil C
                gavgscmsrow(i,m) = gavgscmsrow(i,m) + upMossSoilCrow(i,m,k)
              end do 

              ! Add the moss C mass 
              vgbiomasrow(i,m) = vgbiomasrow(i,m) + Cmossmasrow(i,m)

            end if   

          end if 
        end do ! loop 116
      end do
    end do ! loop 115

    end associate
  end subroutine ctemInit
  
  !> \namespace ctemutilities
  !! Central module for CTEM (biogeochemical)-related utilities
end module ctemUtilities
