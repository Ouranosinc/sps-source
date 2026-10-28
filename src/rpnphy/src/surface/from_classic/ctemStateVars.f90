!> \file
!> Contains the biogeochemistry-related variable type structures.
!! @author J. Melton
!! Variable types herein:
!!
!! 1. vrot (veg_rot) - CTEM's 'rot' vars
!! 2. vgat (veg_gat) - CTEM's 'gat' vars
!! 3. ctem_tile (ctem_tile_level) - CTEM's variables per tile

module ctemStateVars

  ! S.R.C.    April 2022 - This module was modified to statically allocate these variables thus removing allocCtemVars

  ! J. Melton Apr 2015

  use classicParams,  only : ican, icp1, icc, iccp2, iccp1

  implicit none

  public :: initRowVarsBioGeoChem     ! Initializes 'row' variables
  public :: resetMosaicAccum ! Resets physics accumulator variables (used as input to CTEM) after CTEM has been called
  public :: ctemdump ! Dumps all the ctemdump statevars for diagnostic purposes.

  !=================================================================================
  !> CTEM's 'rot' vars
  type veg_rot(nlat, nmos, ignd)
  
    integer, len :: nlat
    integer, len :: nmos
    integer, len :: ignd

    logical, dimension(nlat,nmos,icc) :: pftexist  !< logical array indicating pfts exist(t) or not(f)
    real,    dimension(nlat,nmos,icc) :: cc        !< colonization rate
    real,    dimension(nlat,nmos,icc) :: mm        !< mortality rate
    integer, dimension(nlat,nmos,icc) :: lfstatus  !< leaf phenology status
    integer, dimension(nlat,nmos,icc) :: pandays   !< days with positive net photosynthesis(an) for use in
    !< the phenology subroutine
    real, dimension(nlat,nmos,icc) :: gleafmas     !< green leaf mass for each of the ctem pfts, \f$kg c/m^2\f$
    real, dimension(nlat,nmos,icc) :: gleafmas_ns  !< non-structural green leaf mass for each of the ctem pfts, \f$kg c/m^2\f$
    real, dimension(nlat,nmos,icc) :: gleafmas_s   !< structural green leaf mass for each of the ctem pfts, \f$kg c/m^2\f$
    real, dimension(nlat,nmos,icc) :: leafns2s     !< carbon flux from non-structural to structural leaf pool, \f$kg c/m^2.day\f$
    real, dimension(nlat,nmos,icc) :: stemns2s     !< carbon flux from non-structural to structural stem pool, \f$kg c/m^2.day\f$
    real, dimension(nlat,nmos,icc) :: rootns2s     !< carbon flux from non-structural to structural root pool, \f$kg c/m^2.day\f$
    real, dimension(nlat,nmos,icc) :: re_alloc_s2l !< amount of nsc reallocated from stem to leaves during leafout
    real, dimension(nlat,nmos,icc) :: re_alloc_r2l !< amount of nsc reallocated from root to leaves during leafout
    real, dimension(nlat,nmos,icc) :: re_alloc_sr2l!< amount of nsc reallocated from stem and root to leaves during leafout
    real, dimension(nlat,nmos,icc) :: bleafmas     !< brown leaf mass for each of the ctem pfts, \f$kg c/m^2\f$
    real, dimension(nlat,nmos,icc) :: stemmass     !< stem mass for each of the ctem pfts, \f$kg c/m^2\f$
    real, dimension(nlat,nmos,icc) :: stemmass_ns  !< non-structural stem mass for each of the ctem pfts, \f$kg c/m^2\f$
    real, dimension(nlat,nmos,icc) :: stemmass_s   !< structural stem mass for each of the ctem pfts, \f$kg c/m^2\f$
    real, dimension(nlat,nmos,icc) :: rootmass     !< root mass for each of the ctem pfts, \f$kg c/m^2\f$
    real, dimension(nlat,nmos,icc) :: rootmass_ns  !< non-structural root mass for each of the ctem pfts, \f$kg c/m^2\f$
    real, dimension(nlat,nmos,icc) :: rootmass_s   !< structural root mass for each of the ctem pfts, \f$kg c/m^2\f$
    real, dimension(nlat,nmos,icc) :: pstemmass    !< stem mass from previous timestep, is value before fire. used by burntobare subroutine
    real, dimension(nlat,nmos,icc) :: pgleafmass   !< root mass from previous timestep, is value before fire. used by burntobare subroutine
    real, dimension(nlat,nmos,icc) :: fcancmx      !< max. fractional coverage of ctem's pfts, but this can be
    !< modified by land-use change,and competition between pfts
    real, dimension(nlat,nmos,icc) :: ngleafmas    !< green leaf nitrogen mass for each of the ctem pfts, \f$g N/m^2\f$
    real, dimension(nlat,nmos,icc) :: ngleafmas_ns !< non-structural green leaf nitrogen mass for each of the ctem pfts, \f$g N/m^2\f$
    real, dimension(nlat,nmos,icc) :: ngleafmas_s  !< structural green leaf nitrogen mass for each of the ctem pfts, \f$g N/m^2\f$
    real, dimension(nlat,nmos,icc) :: nbleafmas    !< brown leaf nitrogen mass for each of the ctem pfts, \f$g N/m^2\f$
    real, dimension(nlat,nmos,icc) :: nstemmass    !< stem nitrogen mass for each of the ctem pfts, \f$g N/m^2\f$
    real, dimension(nlat,nmos,icc) :: nstemmass_ns !<non-structural stem nitrogen mass for each of the ctem pfts, \f$g N/m^2\f$
    real, dimension(nlat,nmos,icc) :: nstemmass_s  !<structural stem nitrogen mass for each of the ctem pfts, \f$g N/m^2\f$
    real, dimension(nlat,nmos,icc) :: nrootmass    !<root nitrogen mass for each of the ctem pfts, \f$g N/m^2\f$
    real, dimension(nlat,nmos,icc) :: nrootmass_ns !<non-structural root nitrogen mass for each of the ctem pfts, \f$g N/m^2\f$
    real, dimension(nlat,nmos,icc) :: nrootmass_s  !<structural root nitrogen mass for each of the ctem pfts, \f$g N/m^2\f$

    real, dimension(nlat,nmos,icc) :: ailcg        !< Green LAI for CTEM's pfts
    real, dimension(nlat,nmos,icc) :: ailcgs       !< Green LAI for canopy over snow sub-area
    real, dimension(nlat,nmos,icc) :: fcancs       !< Fraction of canopy over snow for ctem's pfts
    real, dimension(nlat,nmos,icc) :: fcanc        !< Fractional coverage of carbon pfts, canopy over snow
    real, dimension(nlat,nmos,icc) :: co2i1cg      !< Intercellular CO2 conc for pfts for canopy over ground subarea(Pa) - for single/sunlit leaf
    real, dimension(nlat,nmos,icc) :: co2i1cs      !< Same as above but for shaded leaf(above being co2i1cg)
    real, dimension(nlat,nmos,icc) :: co2i2cg      !< Intercellular CO2 conc for pfts for canopy over snowsubarea(pa) - for single/sunlit leaf
    real, dimension(nlat,nmos,icc) :: co2i2cs      !< Same as above but for shaded leaf(above being co2i2cg)
    real, dimension(nlat,nmos,icc) :: ancsveg      !< Net photosynthetic rate for CTEM's pfts for canopy over snow subarea
    real, dimension(nlat,nmos,icc) :: ancgveg      !< Net photosynthetic rate for CTEM's pfts for canopy over ground subarea
    real, dimension(nlat,nmos,icc) :: rmlcsveg     !< Leaf respiration rate for CTEM' pfts forcanopy over snow subarea
    real, dimension(nlat,nmos,icc) :: rmlcgveg     !< Leaf respiration rate for CTEM' pfts forcanopy over ground subarea
    real, dimension(nlat,nmos,icc) :: slai         !< storage/imaginary lai for phenology purposes
    real, dimension(nlat,nmos,icc) :: ailcb        !< brown lai for ctem's 9 pfts. for now we assume only grasses can have brown lai
    real, dimension(nlat,nmos,icc) :: flhrloss     !< fall or harvest loss for deciduous trees and crops, respectively, \f$kg c/m^2\f$il1
    real, dimension(nlat,nmos,icc) :: flhrloss_ns     !< fall or harvest loss for deciduous trees and crops, respectively, \f$kg c/m^2\f$il1
    real, dimension(nlat,nmos,icc) :: flhrloss_s     !< fall or harvest loss for deciduous trees and crops, respectively, \f$kg c/m^2\f$il1
    real, dimension(nlat,nmos,icc) :: grwtheff     !< growth efficiency. change in biomass per year per unit max.
    !< lai(\f$kg c/m^2\f$)/(m2/m2),for use in mortality subroutine
    real, dimension(nlat,nmos,icc) :: lystmmas     !< stem mass at the end of last year
    real, dimension(nlat,nmos,icc) :: lyrotmas     !< root mass at the end of last year
    real, dimension(nlat,nmos,icc) :: tymaxlai     !< this year's maximum lai
    real, dimension(nlat,nmos,icc) :: stmhrlos     !< stem harvest loss for crops, \f$kg c/m^2\f$
    real, dimension(nlat,nmos,icc) :: vgbiomas_veg !< vegetation biomass for each pft
    real, dimension(nlat,nmos,icc) :: emit_co2     !< carbon dioxide
    real, dimension(nlat,nmos,icc) :: emit_co      !< carbon monoxide
    real, dimension(nlat,nmos,icc) :: emit_ch4     !< methane
    real, dimension(nlat,nmos,icc) :: emit_nmhc    !< non-methane hydrocarbons
    real, dimension(nlat,nmos,icc) :: emit_h2      !< hydrogen gas
    real, dimension(nlat,nmos,icc) :: emit_nox     !< nitrogen oxides
    real, dimension(nlat,nmos,icc) :: emit_n2o     !< nitrous oxide
    real, dimension(nlat,nmos,icc) :: emit_nh3     !< ammonia (kg <species> $m^{-2}$$s^{-1}$)
    real, dimension(nlat,nmos,icc) :: emit_pm25    !< particulate matter less than 2.5 um in diameter
    real, dimension(nlat,nmos,icc) :: emit_tpm     !< total particulate matter
    real, dimension(nlat,nmos,icc) :: emit_tc      !< total carbon
    real, dimension(nlat,nmos,icc) :: emit_bc      !< black carbon
    real, dimension(nlat,nmos,icc) :: emit_oc      !< organic carbon
    real, dimension(nlat,nmos,icc) :: emit_so2     !< sulfur dioxide
    real, dimension(nlat)          :: elcgrow      !< prognostic cloud-to-ground lightning \f$flashes/km^2/year\f$
    real, dimension(nlat,nmos,icc) :: burnvegf     !< per PFT fraction burned of that PFT's area
    real, dimension(nlat,nmos,icc) :: smfuncveg    !<
    real, dimension(nlat,nmos,icc) :: bterm        !< biomass term for fire probabilty calc
    real, dimension(nlat,nmos,icc) :: drgtstrs     !< biomass term for fire probabilty calc
    real, dimension(nlat,nmos,ignd) :: betadrgt    !< biomass term for fire probabilty calc
    real, dimension(nlat,nmos,icc) :: fFireCVeg    !< Live biomass combusted by fire, \f$kg c/m^2\f$
    real, dimension(nlat,nmos,icc) :: fFireLitter  !< above ground necromass combusted by fire (considers layer 1 litter as surface), \f$kg c/m^2\f$
    real, dimension(nlat,nmos,icc) :: fFireCSoil   !< soil carbon combusted by fire (considers soil c + litter(2:ignd)), \f$kg c/m^2\f$

    real :: SRHrow(nlat,nmos)
    real :: isirow(nlat,nmos)        !< The initial spead index (ISI) used by FWI fire, unitless
    real :: buirow(nlat,nmos)        !< The build up index (BUI) used by FWI fire, unitless
    real :: fwirow(nlat,nmos)        !< The fire weather index (FWI) used by FWI fire, unitless
    real :: ffmcrow(nlat,nmos)       !< The fine fuel moisture code (FFMC) used by FWI fire, unitless
    real :: dmcrow(nlat,nmos)        !< The duff moisture code (DMC) used by FWI fire, unitless
    real :: dcrow(nlat,nmos)         !< The drought code (DC) used by FWI fire, unitless

    real, dimension(nlat,nmos,icc) :: mterm        !< moisture term for fire probabilty calc
    real, dimension(nlat,nmos,icc) :: bmasveg      !< total(gleaf + stem + root) biomass for each ctem pft, \f$kg c/m^2\f$
    real, dimension(nlat,nmos,icc) :: veghght      !< vegetation height(meters)
    real, dimension(nlat,nmos,icc) :: rootdpth     !< 99% soil rooting depth(meters)
    !< both veghght & rootdpth can be used as diagnostics to see
    !< how vegetation grows above and below ground, respectively
    real, dimension(nlat,nmos,icc) :: tltrleaf     !< total leaf litter fall rate(u-mol co2/m2.sec)
    real, dimension(nlat,nmos,icc) :: tltrstem     !< total stem litter fall rate(u-mol co2/m2.sec)
    real, dimension(nlat,nmos,icc) :: tltrroot     !< total root litter fall rate(u-mol co2/m2.sec)
    real, dimension(nlat,nmos,icc) :: leaflitr     !< leaf litter fall rate(u-mol co2/m2.sec). this leaf litter
    !< does not include litter generated due to mortality/fire
    real, dimension(nlat,nmos,icc) :: roottemp     !< root temperature, k
    real, dimension(nlat,nmos,icc) :: afrleaf      !< allocation fraction for leaves
    real, dimension(nlat,nmos,icc) :: afrstem      !< allocation fraction for stem
    real, dimension(nlat,nmos,icc) :: afrroot      !< allocation fraction for root
    real, dimension(nlat,nmos,icc) :: wtstatus     !< soil water status used for calculating allocation fractions
    real, dimension(nlat,nmos,icc) :: ltstatus     !< light status used for calculating allocation fractions
    real, dimension(nlat,nmos,icc) :: gppveg       !< ! gross primary productity for each pft
    real, dimension(nlat,nmos,icc) :: vcmax0       !< max. photosynthetic rate at the top of canopy(\f$(mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(nlat,nmos,icc) :: avesmfunc    !< soil mositure scalar on photosnthesis (uniless)
    real, dimension(nlat,nmos,icc) :: nppveg       !< npp for individual pfts, u-mol co2/m2.sec
    real, dimension(nlat,nmos,icc) :: autoresveg   !<
    real, dimension(nlat,nmos,icc) :: rmlvegacc    !<
    real, dimension(nlat,nmos,icc) :: rmsveg       !< stem maintenance resp. rate for each pft
    real, dimension(nlat,nmos,icc) :: rmrveg       !< root maintenance resp. rate for each pft
    real, dimension(nlat,nmos,icc) :: rgveg        !< growth resp. rate for each pft
    real, dimension(nlat,nmos,icc) :: litrfallveg  !< litter fall in for each pft(\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(nlat,nmos,icc) :: rothrlos     !< root death as crops are harvested, \f$kg c/m^2\f$
    real, dimension(nlat,nmos,icc) :: pfcancmx     !< previous year's fractional coverages of pfts
    real, dimension(nlat,nmos,icc) :: nfcancmx     !< next year's fractional coverages of pfts
    real, dimension(nlat,nmos,icc) :: anveg        !< net photosynthesis rate for each pft
    real, dimension(nlat,nmos,icc) :: rmlveg       !< leaf maintenance resp. rate for each pft

    real, dimension(nlat,nmos,iccp1) :: bnf_free         !< free-living biological nitrogen fixation(\f$g N/m^2day\f$)
    real, dimension(nlat,nmos,icc) :: bnf_ant          !< anthropogenic biological nitrogen fixation(\f$g N/m^2day\f$)
    real, dimension(nlat,nmos,icc) :: bnf_nat          !< natural biological nitrogen fixation(\f$g N/m^2day\f$)
    real, dimension(nlat,nmos,iccp1) :: bnf_tot          !< total biological nitrogen fixation(\f$g N/m^2day\f$)
    real, dimension(nlat,nmos)       :: bnf_moss       !< moss-related biological nitrogen fixation(\f$g N/m^2day\f$)
    real, dimension(nlat,nmos,icc) :: nstress          !< N stress
    real, dimension(nlat,nmos,iccp1) :: nitrifveg              !< nitrification flux from nh4_mass to no3_mass pool(\f$g N/m^2 day\f$) for individual PFTs + bareground
    real, dimension(nlat,nmos,iccp1) :: no_nitveg              !< NO loss through nitrification(\f$g N/m^2 day\f$) for individual PFTs + bareground
    real, dimension(nlat,nmos,iccp1) :: no_denitveg            !< NO loss through denitrification(\f$g N/m^2 day\f$) for individual PFTs + bareground
    real, dimension(nlat,nmos,iccp1) :: no_nitdenitveg         !< total NO loss from denitrification and nitrification(\f$g N/m^2 day\f$) for individual PFTs + bareground
    real, dimension(nlat,nmos,iccp1) :: n2o_nitveg             !< N2O loss through nitrification(\f$g N/m^2 day\f$) for individual PFTs + bareground
    real, dimension(nlat,nmos,iccp1) :: n2o_denitveg           !< N2O loss through denitrification(\f$g N/m^2 day\f$) for individual PFTs + bareground
    real, dimension(nlat,nmos,iccp1) :: n2o_nitdenitveg        !< total N2O loss from denitrification and nitrification(\f$g N/m^2 day\f$) for individual PFTs + bareground
    real, dimension(nlat,nmos,iccp1) :: n2_denitveg            !< N2 loss through denitrification(\f$g N/m^2 day\f$) for individual PFTs + bareground
    real, dimension(nlat,nmos,iccp1) :: nvolveg                !< nitrogen volatilization(\f$g N/m^2 day\f$) for individual PFTs + bareground
    real, dimension(nlat,nmos,iccp1) :: nleachveg              !< nitrogen leaching(\f$g N/m^2 day\f$) for individual PFTs + bareground
    real, dimension(nlat,nmos,iccp1) :: appl_fert              !< applied nitrogen fertilizer \f$(g N m^{-2} cropland day^{-1}\f$
    real, dimension(nlat,nmos,iccp1) :: ndep_nh4               !< deposition influx into the Ammonium pool for individual PFTs + bareground \f$(g N m^{-2} day^{-1}\f$
    real, dimension(nlat,nmos,iccp1) :: ndep_no3               !< deposition influx into the Nitrate pool for individual PFTs + bareground \f$(g N m^{-2} day^{-1}\f$
    real, dimension(nlat,nmos,icc) :: ndemandveg_wp_npp      !< whole plant npp-based nitrogen demand for individual PFTs(\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,icc) :: nuptakeveg_p_nh4       !< passive nh4+ uptake for individual PFTs(\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,icc) :: nuptakeveg_p_no3       !< passive no3- uptake for individual PFTs(\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,icc) :: nuptakeveg_a_actl_nh4  !< actual active nh4+ uptake for individual PFTs(\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,icc) :: nuptakeveg_a_actl_no3  !< actual active no3- uptake for individual PFTs(\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,icc) :: nuptakeveg             !< total N uptake (active+passive, NH4+NO3) for individual PFTs (\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,icc) :: nleafns2sveg           !< nitrogen flux from non-structural to structural leaf pool(\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,icc) :: nstemns2sveg           !< nitrogen flux from non-structural to structural stem pool(\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,icc) :: nrootns2sveg           !< nitrogen flux from non-structural to structural root pool(\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,icc) :: nallocveg_l            !< nitrogen allocation to leaves for individual PFTs(\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,icc) :: nallocveg_s            !< nitrogen allocation to stem for individual PFTs(\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,icc) :: nallocveg_r            !< nitrogen allocation to root for individual PFTs(\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,icc) :: nresorpedveg_s         !< resorped N from leaves to be allocated to stem(\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,icc) :: nresorpedveg_r         !< resorped N from leaves to be allocated to root(\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,icc) :: nre_allocveg_s2l       !< reallocated N from S to L during leaf out period(\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,icc) :: nre_allocveg_r2l       !< reallocated N from R to L during leaf out period(\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,icc) :: nlitrveg_l             !< leaf N litterfall(\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,icc) :: nlitrveg_s             !< stem N litterfall(\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,icc) :: nlitrveg_r             !< root N litterfall(\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,icc) :: nlitrveg               !< total N litterfall (\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,icc) :: gl2bl_grass_nflux      !< N flux from ngleafmas to nbleafmas(\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,icc) :: c2nveg_l               !< simulated C:N ratio for leaves(\f$g C/g N\f$)
    real, dimension(nlat,nmos,icc) :: c2nveg_s               !< simulated C:N ratio for stem(\f$g C/g N\f$)
    real, dimension(nlat,nmos,icc) :: c2nveg_r               !< simulated C:N ratio for roots(\f$g C/g N\f$)
    real, dimension(nlat,nmos,icc) :: c2nveg_wp              !< simulated C:N ratio for the whole plant(\f$g C/g N\f$)
    real, dimension(nlat,nmos,iccp1) :: c2nveg_litr            !< simulated C:N ratio for litter mass(\f$g C/g N\f$)
    real, dimension(nlat,nmos,iccp1) :: c2nveg_humus           !< simulated C:N ratio for humus mass(\f$g C/g N\f$)
    real, dimension(nlat,nmos,iccp1) :: nhumtrsveg             !< N humification for individual PFTs + bare(\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,iccp1) :: nmineralveg_litr       !< N mineralization from litter pool(\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,iccp1) :: nmineralveg_humus      !< N mineralization from organic soil pool(\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,iccp1) :: nimmobilveg_nh4        !< N immobilization from NH4+ pool to soilnmas(\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,iccp1) :: nimmobilveg_no3        !< N immobilization from NO3- pool to soilnmas(\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,iccp1) :: netnmineralveg         !< N mineralization (\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,icc) :: nvgbiomas_veg          !< N vegetation biomass for individual PFTs(\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,iccp1) :: fNnetlandveg           !< net terrestrial N flux(\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos,icc) :: redcoeff_vcmax         !< Reduction coeff. passed to the Photosynthesis subroutine

    ! allocated with nlat,nmos:
    real, dimension(nlat,nmos) :: redcoeff_vcmaxMoss    !< Reduction coeff. for moss passed to the Photosynthesis subroutine
    real, dimension(nlat,nmos) :: fNleach               !< grid averaged total N leaching \f$(kg N m^{-2} s^{-1}\f$
    real, dimension(nlat,nmos) :: fNvol                 !< grid averaged total N volatilization from ammonium pool (\f$g N/m^2 day\f$)
    real, dimension(nlat,nmos) :: gavglai               !< grid averaged green leaf area index
    real, dimension(nlat,nmos) :: co2conc               !< ATMOS. CO2 CONC. IN PPM
    real, dimension(nlat,nmos) :: ch4conc               !< atmospheric methane concentration \f$kg c/m^2\f$
    real, dimension(nlat,nmos) :: canres                !<
    real, dimension(nlat,nmos) :: vgbiomas              !< grid averaged vegetation biomass, \f$kg c/m^2\f$
    real, dimension(nlat,nmos) :: gavgltms              !< grid averaged litter mass, \f$kg c/m^2\f$
    real, dimension(nlat,nmos) :: gavgscms              !< grid averaged soil c mass, \f$kg c/m^2\f$
    real, dimension(nlat,nmos) :: burnfrac              !< areal :: fraction burned due to fire for every grid cell(%)
    real, dimension(nlat,nmos) :: popdin                !< population density \f$(people / km^2)\f$
    real, dimension(nlat,nmos) :: soilpH                !< soil PH
    real, dimension(nlat,nmos) :: nfertil               !< nitrogen fertilizer \f$(g N m^{-2} cropland yr^{-1}\f$
    real, dimension(nlat,nmos) :: ndeposit              !< nitrogen deposition \f$(g N m^{-2} yr^{-1}\f$
    real, dimension(nlat) :: timharvrow                 !< the annual fractional area where timber was harvested (0-1,read in from the external forcing file)
    real, dimension(nlat,nmos) :: timharvarearow        !< the fractional area of the cell where timber will be harvested (0-1,this is the actual area harvested per tile by the model, which is determined from timharvrow by the harvest subroutines)
    real, dimension(nlat) :: prsfirerow                 !< the annual fractional area where fire burned (0-1, read in from the external forcing file)
    real, dimension(nlat,nmos) :: prsfirearearow        !< the fractional area of the cell where fire burned (0-1,this is the actual area burned per tile by the model, which is determined from prsfirerow by the prescribed fire subroutines)
    real, dimension(nlat,nmos) :: tileAgerow            !< the age of the tile since the start of the run in months this is reset by harvest and fire and incremented at the CTEM timestep by indexTileAge
    real, dimension(nlat,nmos) :: lterm                 !< lightning term for fire probabilty calc
    real, dimension(nlat,nmos) :: extnprob              !< fire extingusinging probability
    real, dimension(nlat,nmos) :: prbfrhuc              !< probability of fire due to human causes
    real, dimension(nlat,nmos) :: rml                   !< leaf maintenance respiration(\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, dimension(nlat,nmos) :: rms                   !< stem maintenance respiration(\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, dimension(nlat,nmos) :: rmr                   !< root maintenance respiration(\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, dimension(nlat,nmos) :: ch4WetSpec            !< methane flux from wetlands calculated using hetrores in umol ch4/m2.s
    real, dimension(nlat,nmos) :: wetfdyn               !< dynamic wetland fraction
    real, dimension(nlat,nmos) :: wetfrac_pres          !< Prescribed wetland fraction read in from OBSWETFFile
    real, dimension(nlat,nmos) :: ch4WetDyn                !< methane flux from wetlands calculated using hetrores and wetfdyn, in umol ch4/m2.s
    real, dimension(nlat,nmos) :: ch4_soills            !< Methane uptake into the soil column(\f$mg CH_4 m^{-2} s^{-1}\f$)
    real, dimension(nlat,nmos) :: lucemcom              !< land use change(luc) related combustion emission losses(\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, dimension(nlat,nmos) :: lucltrin              !< luc related inputs to litter pool(\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, dimension(nlat,nmos) :: lucsocin              !< luc related inputs to soil c pool(\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, dimension(nlat,nmos) :: cproduct              !< luc related product pool, sum of paper and furniture, (\f$\ Kg C m^{-2} \f$)
    real, dimension(nlat,nmos) :: fproductdecomp        !< luc related decompsition from its product pools, sum of paper and furniture decomposition, (\f$\ mu mol CO2 m^{-2} s^{-1} \f$)
    real, dimension(nlat,nmos) :: cleaf               !< Grid averaged leaf C pool , \f$kg C/m^2\f$
    real, dimension(nlat,nmos) :: cstem               !< Grid averaged stem C pool , \f$kg C/m^2\f$
    real, dimension(nlat,nmos) :: croot               !< Grid averaged root C pool , \f$kg C/m^2\f$
    real, dimension(nlat,nmos) :: nppleaf             !< Grid averaged leaf NPP, \f$\mu mol CO2 m^{-2} s^{-1}\f$
    real, dimension(nlat,nmos) :: nppstem             !< Grid averaged stem NPP, \f$\mu mol CO2 m^{-2} s^{-1}\f$
    real, dimension(nlat,nmos) :: npproot             !< Grid averaged root NPP, \f$\mu mol CO2 m^{-2} s^{-1}\f$
    real, dimension(nlat,nmos) :: nMineralNH4         !< Grid average NH4 mass (for CanESM), \f$kg N/m^2\f$
    real, dimension(nlat,nmos) :: nMineralNO3         !< Grid average NO3 mass (for CanESM), \f$kg N/m^2\f$
    real, dimension(nlat,nmos) :: nVeg                !< Grid average vegetation N mass (for CanESM), \f$kg N/m^2\f$
    real, dimension(nlat,nmos) :: nLitter             !< Grid average litter N mass (for CanESM), \f$kg N/m^2\f$
    real, dimension(nlat,nmos) :: nSoil               !< Grid average soil N mass mass (for CanESM), \f$kg N/m^2\f$
    real, dimension(nlat,nmos) :: nLeaf               !< Grid average leaf N mass mass (for CanESM), \f$kg N/m^2\f$
    real, dimension(nlat,nmos) :: nStem               !< Grid average stem N mass mass (for CanESM), \f$kg N/m^2\f$
    real, dimension(nlat,nmos) :: nRoot               !< Grid average root N mass mass (for CanESM), \f$kg N/m^2\f$
    real, dimension(nlat,nmos) :: fBNF                !< grid average total biological nitrogen fixation \f$(kg N m^{-2} s^{-1}\f$
    real, dimension(nlat,nmos) :: fNdep               !< grid average total nitrogen deposition \f$(kg N m^{-2} s^{-1}\f$
    real, dimension(nlat,nmos) :: fNfert              !< grid average total nitrogen fertilization \f$(kg N m^{-2} s^{-1}\f$
    real, dimension(nlat,nmos) :: fNgasNonFire        !< grid average total nitrogen flux excluding fire \f$(kg N m^{-2} s^{-1}\f$
    real, dimension(nlat,nmos) :: fNgasFire           !< grid average total nitrogen flux from fire only \f$(kg N m^{-2} s^{-1}\f$
    real, dimension(nlat,nmos) :: fNgas               !< grid average total nitrogen flux from all sources including fire \f$(kg N m^{-2} s^{-1}\f$
    real, dimension(nlat,nmos) :: fNnetmin            !< grid average total net mineralization \f$(kg N m^{-2} s^{-1}\f$
    real, dimension(nlat,nmos) :: fNOx                !< grid average total NOx flux from all sources including fire \f$(kg N m^{-2} s^{-1}\f$
    real, dimension(nlat,nmos) :: fNup                !< grid average pasive N updatake from the NH4 and NO3 pools \f$(kg N m^{-2} s^{-1}\f$
    real, dimension(nlat,nmos) :: fNvegSoil           !< grid average total nitrogen flux from vegetation to soil [0 at present] \f$(kg N m^{-2} s^{-1}\f$
    real, dimension(nlat,nmos) :: fN2o                !< grid average total N20 flux from all sources \f$(kg N m^{-2} s^{-1}\f$
    real, dimension(nlat,nmos) :: lucemcomn   !< land use change (luc) related combustion emission losses of N (g N m^-2 day^-1)
    real, dimension(nlat,nmos) :: lucltrinn   !< luc related inputs to litter pool of N (g N m^-2 day^-1)
    real, dimension(nlat,nmos) :: lucsocinn   !< luc related inputs to soil N pool of N (g N m^-2 day^-1)
    real, dimension(nlat,nmos) :: npp                   !< net primary productivity
    real, dimension(nlat,nmos) :: nep                   !< net ecosystem productivity
    real, dimension(nlat,nmos) :: nepCMIP               !< net ecosystem productivity 
    real, dimension(nlat,nmos) :: nbp                   !< net biome productivity
    real, dimension(nlat,nmos) :: gpp                   !< gross primary productivity
    real, dimension(nlat,nmos) :: hetrores              !< heterotrophic respiration
    real, dimension(nlat,nmos) :: autores               !< autotrophic respiration
    real, dimension(nlat,nmos) :: soilcresp             !<
    real, dimension(nlat,nmos) :: rm                    !< maintenance respiration
    real, dimension(nlat,nmos) :: rg                    !< growth respiration
    real, dimension(nlat,nmos) :: litres                !< litter respiration
    real, dimension(nlat,nmos) :: socres                !< soil carbon respiration
    real, dimension(nlat,nmos) :: dstcemls              !< carbon emission losses due to disturbance, mainly fire
    real, dimension(nlat,nmos) :: litrfall              !< total litter fall(from leaves, stem, and root) due to
    !< all causes(mortality,turnover,and disturbance)(\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(nlat,nmos) :: humiftrs              !< transfer of humidified litter from litter to soil c pool
    real, dimension(nlat,nmos) :: cfluxcg               !<
    real, dimension(nlat,nmos) :: cfluxcs               !<
    real, dimension(nlat,nmos) :: ROFB                  !< Base flow from bottom of soil column \f$[kg m^{-2} s^{-1} ]\f$
    real, dimension(nlat,nmos) :: dstcemls3             !< carbon emission losses due to disturbance(fire at present) from litter pool
    real, dimension(nlat,nmos) :: uvaccrow_m            !<
    real, dimension(nlat,nmos) :: vvaccrow_m            !<
    real, dimension(nlat,nmos) :: qevpacc_m_save        !<
    real, dimension(nlat)      :: twarmm                !< temperature of the warmest month(c)
    real, dimension(nlat)      :: tcoldm                !< temperature of the coldest month(c)
    real, dimension(nlat)      :: gdd5                  !< growing degree days above 5 c
    real, dimension(nlat)      :: aridity               !< aridity index, ratio of potential evaporation to precipitation
    real, dimension(nlat)      :: srplsmon              !< number of months in a year with surplus water i.e. precipitation more than potential evaporation
    real, dimension(nlat)      :: defctmon              !< number of months in a year with water deficit i.e. precipitation less than potential evaporation
    real, dimension(nlat)      :: anndefct              !< annual water deficit(mm)
    real, dimension(nlat)      :: annsrpls              !< annual water surplus(mm)
    real, dimension(nlat)      :: annpcp                !< annual precipitation(mm)
    real, dimension(nlat)      :: dry_season_length     !< length of dry season(months)

    character(len=8), dimension(nlat,nmos) :: peatlandType !< Peatland type (else 'None')
    real, dimension(nlat,nmos) :: peatSoilC         !< peat soil C mass, \f$kg C/m^2\f$
    real, dimension(nlat,nmos) :: peatdep

    character(len=8), dimension(nlat,nmos) :: mossPresent !< Type of moss present (else 'None') 
    real, dimension(nlat,nmos,ignd) :: litrmsmoss      !< Moss litter C mass, \f$kg C/m^2\f$ 
    real, dimension(nlat,nmos,ignd) :: upMossSoilC     !< Moss soil C mass (non-peat), \f$kg C/m^2\f$
    real, dimension(nlat,nmos) :: Cmossmas

    real, dimension(nlat,nmos) :: litrmsmossN          !< Moss litter N mass, \f$g N/m^2\f$
    real, dimension(nlat,nmos) :: upMossSoilN          !< Moss soil N mass (non-peat), \f$g N/m^2\f$
    real, dimension(nlat,nmos) :: Nmossmas             !< Moss N mass, \f$g N/m^2\f$
    
    real, dimension(nlat,nmos) :: dmoss
    real, dimension(nlat,nmos) :: nppmoss
    real, dimension(nlat,nmos) :: armoss

    integer, dimension(nlat,nmos) :: colddays_leaffall          !< cold days counter for tracking days below a certain
    !< temperature threshold for ndl dcd tree.
    integer, dimension(nlat,nmos) :: colddays_harvest          !< cold days counter for tracking days below a certain
    !< temperature threshold for crops.

    ! allocated with nlat,nmos,ican:
    real, dimension(nlat,nmos,ican) :: zolnc            !< lumped log of roughness length for class' 4 pfts
    real, dimension(nlat,nmos,ican) :: ailc             !< lumped lai for class' 4 pfts
    real, dimension(nlat,nmos,ican) :: cmasvegc         !< total canopy mass for each of the 4 class pfts. recall that
    !< class requires canopy mass as an input,and this is now provided by ctem. \f$kg/m^2\f$.
    real, dimension(nlat,nmos,ican) :: alvsctm          !<
    real, dimension(nlat,nmos,ican) :: paic             !< plant area index for class' 4 pfts. this is the sum of leaf
    !< area index and stem area index.
    real, dimension(nlat,nmos,ican) :: slaic            !< storage lai. this will be used as min. lai that class sees
    !< so that it doesn't blow up in its stomatal conductance calculations.
    real, dimension(nlat,nmos,ican) :: alirctm          !<

    ! allocated with nlat,nmos,ican,ignd:
    real, dimension(nlat,nmos,ican,ignd) :: rmatc       !< fraction of roots for each of class' 4 pfts in each soil layer

    ! allocated with nlat,nmos,icc,ignd:
    real, dimension(nlat,nmos,icc,ignd) :: rmatctem     !< fraction of roots for each of ctem's 9 pfts in each soil layer

    ! allocated with nlat,nmos,iccp1:
    real, dimension(nlat,nmos,iccp1) :: nepveg      !< net ecosystem productity for each pft
    real, dimension(nlat,nmos,iccp1) :: nbpveg      !< net biome productity for bare fraction OR net biome productity for each pft
    real, dimension(nlat,nmos,iccp1) :: hetroresveg !<

    ! allocated with nlat,nmos,iccp2,ignd:
    real, dimension(nlat,nmos,iccp2,ignd) :: litrmass    !< litter mass for each of the 9 ctem pfts + bare, \f$kg c/m^2\f$
    real, dimension(nlat,nmos,iccp2,ignd) :: soilcmas    !< soil carbon mass for each of the 9 ctem pfts + bare, \f$kg c/m^2\f$
    real, dimension(nlat,nmos,iccp2,ignd) :: litresveg   !<
    real, dimension(nlat,nmos,iccp2,ignd) :: soilcresveg !<
    real, dimension(nlat,nmos,iccp2,ignd) :: humiftrsveg !< Transfer of humidified litter from litter to soil C pool 

    real, dimension(nlat,nmos,iccp1) :: nh4_mass    !< ammonium mass for individual PFTs + bare, \f$g N/m^2\f$
    real, dimension(nlat,nmos,iccp1) :: no3_mass    !< nitrate mass for individual PFTs + bare, \f$g N/m^2\f$
    real, dimension(nlat,nmos,iccp2) :: nlitrmass   !< litter nitrogen mass for individual PFTs + bare, \f$g N/m^2\f$
    real, dimension(nlat,nmos,iccp2) :: soilnmas    !< soil organic nitrogen mass for individual PFTs + bare, \f$g N/m^2\f$

    ! allocated with nlat,nmos,{some number}:
    real, dimension(nlat,nmos,8)  :: slopefrac          !< prescribed fraction of wetlands based on slope
    !< only(0.025,0.05,0.1,0.15,0.20,0.25,0.3 and 0.35 percent slope thresholds)

    ! allocated with nlat:
    real, dimension(nlat)    :: dayl_max        !< maximum daylength for that location(hours)
    real, dimension(nlat)    :: dayl            !< daylength for that location(hours)
    real, dimension(nlat)    :: grclarearow     !< area of the grid cell, \f$km^2\f$

    ! allocated with nlat,nmos,icc:
    real, dimension(nlat,nmos,icc) :: lygleafmasmax        !< last year maximum of the gleafmas
    real, dimension(nlat,nmos,icc) :: lystemmassmax        !< last year maximum of the stemmass
    real, dimension(nlat,nmos,icc) :: lyrootmassmax        !< last year maximum of the rootmass
    real, dimension(nlat,nmos,icc) :: lmaxt, smaxt, rmaxt  !< temp vars to find previous year C pool max

    !input variables that control the dynamic tiling subtoutine
    integer, dimension(nlat,nmos) :: controlVector(nlat,nmos) !< this is the control vector which is specified in the initializtion file
    !integer, dimension(nlat,nmos) :: inputVector !< this is the input vector which is generated by the script (1 = operate on tile, 0 = do not operate on tile) 
    !character(len=5), dimension(nlat) :: mode !< mode is a character string that specifies the operation (i.e. move, split or copy) 
    !real, dimension(nlat) :: FAREAROUT !< a real value greater than zero and less than one which specifies the desired size of the output tile
    logical :: DynTilInitializeFlag !< this flag is used to re-initialize the model after dynamic tiling runs

  end type veg_rot

  !=================================================================================
  !> CTEM's 'gat' vars
  type veg_gat(ilg, ignd)

    integer, len :: ilg
    integer, len :: ignd
    
    ! This is the basic data structure that contains the state variables
    ! for the Plant Functional Type (PFT). The dimensions are ilg,{icc,iccp1,iccp2}

    real, dimension(ilg,icc) :: gleafmas   !< Green leaf mass for each of the CTEM PFTs, \f$kg c/m^2\f$
    real, dimension(ilg,icc) :: gleafmas_ns  !< non-structural green leaf mass for each of the CTEM PFTs, \f$kg c/m^2\f$

    
    real, dimension(ilg,icc) :: gleafmas_s   !< structural green leaf mass for each of the CTEM PFTs, \f$kg c/m^2\f$
    real, dimension(ilg,icc) :: leafns2s     !< carbon flux from non-structural to structural leaf pool, \f$kg c/m^2.day\f$
    real, dimension(ilg,icc) :: stemns2s     !< carbon flux from non-structural to structural stem pool, \f$kg c/m^2.day\f$
    real, dimension(ilg,icc) :: rootns2s     !< carbon flux from non-structural to structural root pool, \f$kg c/m^2.day\f$
    real, dimension(ilg,icc) :: re_alloc_s2l !< amount of nsc reallocated from stem to leaves during leaf out, \f$g c/m^2\f$.
    real, dimension(ilg,icc) :: re_alloc_r2l !< amount of nsc reallocated from root to leaves during leaf out, \f$g c/m^2\f$.
    real, dimension(ilg,icc) :: re_alloc_sr2l!< amount of nsc reallocated from stem and root to leaves during leaf out, \f$g c/m^2\f$.
    real, dimension(ilg,icc) :: bleafmas   !< Brown leaf mass for each of the CTEM PFTs, \f$kg c/m^2\f$
    real, dimension(ilg,icc) :: stemmass   !< Stem mass for each of the CTEM PFTs, \f$kg c/m^2\f$
    real, dimension(ilg,icc) :: stemmass_ns  !< non-structural stem mass for each of the CTEM PFTs, \f$kg c/m^2\f$
    real, dimension(ilg,icc) :: stemmass_s   !< structural stem mass for each of the CTEM PFTs, \f$kg c/m^2\f$
    real, dimension(ilg,icc) :: rootmass   !< Root mass for each of the CTEM PFTs, \f$kg c/m^2\f$
    real, dimension(ilg,icc) :: rootmass_ns  !< non-structural root mass for each of the CTEM PFTs, \f$kg c/m^2\f$
    real, dimension(ilg,icc) :: rootmass_s   !< structural root mass for each of the CTEM PFTs, \f$kg c/m^2\f$
    real, dimension(ilg,icc) :: pstemmass  !< Stem mass from previous timestep, is value before fire. used by burntobare subroutine
    real, dimension(ilg,icc) :: pgleafmass !< Root mass from previous timestep, is value before fire. used by burntobare subroutine
    real, dimension(ilg,icc) :: fcancmx    !< max. fractional coverage of ctem's 9 pfts, but this can be
    !< modified by land-use change,and competition between pfts
    real, dimension(ilg,icc) :: ngleafmas    !< green leaf nitrogen mass for each of the CTEM PFTs, \f$g N/m^2\f$
    real, dimension(ilg,icc) :: ngleafmas_ns !< non-structural green leaf nitrogen mass for each of the CTEM PFTs, \f$g N/m^2\f$
    real, dimension(ilg,icc) :: ngleafmas_s  !< structural green leaf nitrogen mass for each of the CTEM PFTs, \f$g N/m^2\f$
    real, dimension(ilg,icc) :: nbleafmas    !< brown leaf nitrogen mass for each of the CTEM PFTs, \f$g N/m^2\f$
    real, dimension(ilg,icc) :: nstemmass    !< stem nitrogen mass for each of the CTEM PFTs, \f$g N/m^2\f$
    real, dimension(ilg,icc) :: nstemmass_ns !< non-structural stem nitrogen mass for each of the CTEM PFTs, \f$g N/m^2\f$
    real, dimension(ilg,icc) :: nstemmass_s  !< structural stem nitrogen mass for each of the CTEM PFTs, \f$g N/m^2\f$
    real, dimension(ilg,icc) :: nrootmass    !< root nitrogen mass for each of the CTEM PFTs, \f$g N/m^2\f$
    real, dimension(ilg,icc) :: nrootmass_ns !< non-structural nitrogen root mass for each of the CTEM PFTs, \f$g N/m^2\f$
    real, dimension(ilg,icc) :: nrootmass_s  !< structural root nitrogen mass for each of the CTEM PFTs, \f$g N/m^2\f$
    real, dimension(ilg,icc) :: nvgbiomas_veg!< N vegetation biomass for each of the CTEM PFTs, \f$g N/m^2\f$
    real, dimension(ilg) :: gavglai        !< grid averaged green leaf area index

    real, dimension(ilg) :: lightng        !< total lightning frequency, flashes/km2.year

    real, dimension(ilg,ican) :: zolnc     !< lumped log of roughness length for class' 4 pfts
    real, dimension(ilg,ican) :: ailc      !< lumped lai for class' 4 pfts

    real, dimension(ilg,icc) :: ailcg      !< green lai for ctem's 9 pfts
    real, dimension(ilg,icc) :: ailcgs     !< GREEN LAI FOR CANOPY OVER SNOW SUB-AREA
    real, dimension(ilg,icc) :: fcancs     !< FRACTION OF CANOPY OVER SNOW FOR CTEM's 9 PFTs
    real, dimension(ilg,icc) :: fcanc      !< FRACTIONAL COVERAGE OF 8 CARBON PFTs, CANOPY OVER GROUND

    real, dimension(ilg)   :: co2conc    !< ATMOS. CO2 CONC. IN PPM
    real, dimension(ilg)   :: ch4conc    !<

    real, dimension(ilg,icc) :: co2i1cg    !< Intercellular CO2 concentration for CTEM PFTs for canopy over ground subarea (Pa) - for single/sunlit leaf
    real, dimension(ilg,icc) :: co2i1cs    !< Intercellular CO2 concentration for CTEM PFTs for canopy over snow subarea (Pa) - for single/sunlit leaf
    real, dimension(ilg,icc) :: co2i2cg    !< Intercellular CO2 concentration for CTEM PFTs for canopy over ground subarea (Pa) - for shaded leaf
    real, dimension(ilg,icc) :: co2i2cs    !< Intercellular CO2 concentration for CTEM PFTs for canopy over snow subarea (Pa) - for shaded leaf
    real, dimension(ilg,icc) :: ancsveg    !< net photosynthetic rate for ctems 9 pfts for canopy over snow subarea
    real, dimension(ilg,icc) :: ancgveg    !< net photosynthetic rate for ctems 9 pfts for canopy over ground subarea
    real, dimension(ilg,icc) :: rmlcsveg   !< leaf respiration rate for ctems 9 pfts forcanopy over snow subarea
    real, dimension(ilg,icc) :: rmlcgveg   !< leaf respiration rate for ctems 9 pfts forcanopy over ground subarea
    real, dimension(ilg,icc) :: slai       !< storage/imaginary lai for phenology purposes
    real, dimension(ilg,icc) :: ailcb      !< brown lai for ctem's 9 pfts. for now we assume only grasses can have brown lai
    real, dimension(ilg)   :: canres     !<
    real, dimension(ilg,icc) :: flhrloss   !< fall or harvest loss for deciduous trees and crops, respectively, \f$kg c/m^2\f$il1
    real, dimension(ilg,icc) :: flhrloss_ns     !< fall or harvest loss for deciduous trees and crops, respectively, \f$kg c/m^2\f$il1
    real, dimension(ilg,icc) :: flhrloss_s     !< fall or harvest loss for deciduous trees and crops, respectively, \f$kg c/m^2\f$il1
    real, dimension(ilg,icc) :: grwtheff   !< growth efficiency. change in biomass per year per unit max.
    !< lai (\f$kg c/m^2\f$)/(m2/m2),for use in mortality subroutine
    real, dimension(ilg,icc) :: lystmmas   !< stem mass at the end of last year
    real, dimension(ilg,icc) :: lyrotmas   !< root mass at the end of last year
    real, dimension(ilg,icc) :: tymaxlai   !< this year's maximum lai
    real, dimension(ilg)   :: vgbiomas   !< grid averaged vegetation biomass, \f$kg c/m^2\f$
    real, dimension(ilg)   :: gavgltms   !< grid averaged litter mass, \f$kg c/m^2\f$
    real, dimension(ilg)   :: gavgscms   !< grid averaged soil c mass, \f$kg c/m^2\f$
    real, dimension(ilg,icc) :: stmhrlos   !< stem harvest loss for crops, \f$kg c/m^2\f$
    real, dimension(ilg,ican,ignd) :: rmatc    !< fraction of roots for each of class' 4 pfts in each soil layer
    real, dimension(ilg,icc,ignd) :: rmatctem !< fraction of roots for each of ctem's 9 pfts in each soil layer
    real, dimension(ilg,iccp2,ignd) :: litrmass   !< Litter mass for each of the CTEM PFTs + bare + LUC product pools, \f$kg c/m^2\f$
    real, dimension(ilg,iccp2,ignd) :: soilcmas   !< Soil carbon mass for each of the CTEM PFTs + bare + LUC product pools, \f$kg c/m^2\f$
    real, dimension(ilg,icc) :: vgbiomas_veg !< vegetation biomass for each pft
    real, dimension(ilg,iccp1) :: nh4_mass   !< ammonium mass for each of the 9 ctem pfts + bare, \f$g N/m^2\f$
    real, dimension(ilg,iccp1) :: no3_mass   !< nitrate mass for each of the 9 ctem pfts + bare, \f$g N/m^2\f$
    real, dimension(ilg,iccp2) :: nlitrmass  !< litter nitrogen mass for each of the 9 ctem pfts + bare, \f$g N/m^2\f$
    real, dimension(ilg,iccp2) :: soilnmas   !< soil nitrogen mass for each of the 9 ctem pfts + bare, \f$g N/m^2\f$

    real, dimension(ilg,icc) :: emit_co2   !< carbon dioxide (kg <species> $m^{-2}$$s^{-1}$)
    real, dimension(ilg,icc) :: emit_co    !< carbon monoxide (kg <species> $m^{-2}$$s^{-1}$)
    real, dimension(ilg,icc) :: emit_ch4   !< methane (kg <species> $m^{-2}$$s^{-1}$)
    real, dimension(ilg,icc) :: emit_nmhc  !< non-methane hydrocarbons (kg <species> $m^{-2}$$s^{-1}$)
    real, dimension(ilg,icc) :: emit_h2    !< hydrogen gas (kg <species> $m^{-2}$$s^{-1}$)
    real, dimension(ilg,icc) :: emit_nox   !< nitrogen oxides (kg <species> $m^{-2}$$s^{-1}$)
    real, dimension(ilg,icc) :: emit_n2o   !< nitrous oxide (kg <species> $m^{-2}$$s^{-1}$)
    real, dimension(ilg,icc) :: emit_nh3   !< ammonia (kg <species> $m^{-2}$$s^{-1}$)
    real, dimension(ilg,icc) :: emit_pm25  !< particulate matter less than 2.5 um in diameter (kg <species> $m^{-2}$$s^{-1}$)
    real, dimension(ilg,icc) :: emit_tpm   !< total particulate matter (kg <species> $m^{-2}$$s^{-1}$)
    real, dimension(ilg,icc) :: emit_tc    !< total carbon (kg <species> $m^{-2}$$s^{-1}$)
    real, dimension(ilg,icc) :: emit_bc    !< black carbon (kg <species> $m^{-2}$$s^{-1}$)
    real, dimension(ilg,icc) :: emit_oc    !< organic carbon (kg <species> $m^{-2}$$s^{-1}$)
    real, dimension(ilg,icc) :: emit_so2   !< sulfur dioxide (kg <species> $m^{-2}$$s^{-1}$)
    real, dimension(ilg,icc) :: fiengat    !< fire energy per PFT \f$J\f$
    real, dimension(ilg)     :: fabgat     !< fire area burned \f$km^2\f$
    real, dimension(ilg)     :: elcggat    !< prognostic cloud-to-ground lightning \f$flashes/km^2/year\f$
    real, dimension(ilg,icc) :: nmfigat    !< number of fires per grid cell
    real, dimension(ilg)     :: burnfrac   !< areal :: fraction burned due to fire for every grid cell (%)
    real, dimension(ilg,icc) :: burnvegf   !< per PFT fraction burned of that PFT's area
    real, dimension(ilg,icc) :: smfuncveg  !<
        real, dimension(ilg,icc) :: drgtstrs !< soil dryness factor for pfts
    real, dimension(ilg,ignd) :: betadrgt !< dryness term for soil layers

    real :: SRH(ilg), & 
                         isi(ilg), & 
                         bui(ilg), &
                         fwi(ilg), &
                         ffmc(ilg), &
                         dmc(ilg), &
                         dc(ilg)
    real, dimension(ilg)   :: popdin     !< population density (people / \f$km^2\f$)
    real, dimension(ilg)   :: soilpH     !< soil PH
    real, dimension(ilg)   :: nfertil    !< nitrogen fertilizer \f$(g N m^{-2} cropland yr^{-1}\f$
    real, dimension(ilg)   :: ndeposit   !< nitrogen deposition \f$(g N m^{-2} yr^{-1}\f$
    real, dimension(ilg) :: timharvareagat   !< the same as timharvarea, but in CTEM's 'gat' format 
    real, dimension(ilg) :: tileAgegat       !< the age of the tile since the start of the run in months this is reset by harvest and fire and incremented at the CTEM timestep by indexTileAge
    real, dimension(ilg) :: prsfireareagat   !< the same as prsfirearea, but in CTEM's 'gat' format
    real, dimension(ilg,icc) :: bterm      !< biomass term for fire probabilty calc
    real, dimension(ilg)   :: lterm      !< lightning term for fire probabilty calc
    real, dimension(ilg,icc) :: mterm      !< moisture term for fire probabilty calc
    real, dimension(ilg,icc) :: glcaemls  !< green leaf carbon emission disturbance losses, \f$kg c/m^2\f$
    real, dimension(ilg,icc) :: blcaemls  !< brown leaf carbon emission disturbance losses, \f$kg c/m^2\f$
    real, dimension(ilg,icc) :: rtcaemls  !< root carbon emission disturbance losses, \f$kg c/m^2\f$
    real, dimension(ilg,icc) :: stcaemls  !< stem carbon emission disturbance losses, \f$kg c/m^2\f$
    real, dimension(ilg,icc) :: ltrcemls  !< litter carbon emission disturbance losses, \f$kg c/m^2\f$
    real, dimension(ilg,icc) :: ntchlveg  !< fluxes for each pft: Net change in leaf biomass, u-mol CO2/m2.sec
    real, dimension(ilg,icc) :: ntchsveg  !< fluxes for each pft: Net change in stem biomass, u-mol CO2/m2.sec
    real, dimension(ilg,icc) :: ntchrveg  !< fluxes for each pft: Net change in root biomass,
    !! the net change is the difference between allocation and
    !! autotrophic respiratory fluxes, u-mol CO2/m2.sec

    real, dimension(ilg)   :: extnprob   !< fire extingusinging probability
    real, dimension(ilg)   :: prbfrhuc   !< probability of fire due to human causes
    real, dimension(ilg)   :: dayl_max   !< maximum daylength for that location (hours)
    real, dimension(ilg)   :: dayl       !< daylength for that location (hours)

    real, dimension(ilg,icc) :: bmasveg    !< total (gleaf + stem + root) biomass for each ctem pft, \f$kg c/m^2\f$
    real, dimension(ilg,ican) :: cmasvegc   !< total canopy mass for each of the 4 class pfts. recall that
    !< class requires canopy mass as an input,and this is now provided by ctem. \f$kg/m^2\f$.
    real, dimension(ilg,icc) :: veghght    !< vegetation height (meters)
    real, dimension(ilg,icc) :: rootdpth   !< 99% soil rooting depth (meters)
    !< both veghght & rootdpth can be used as diagnostics to see
    !< how vegetation grows above and below ground, respectively
    real, dimension(ilg)   :: rml        !< leaf maintenance respiration (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, dimension(ilg)   :: rms        !< stem maintenance respiration (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, dimension(ilg,icc) :: tltrleaf   !< total leaf litter fall rate (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, dimension(ilg,icc) :: blfltrdt !< brown leaf litter generated due to disturbance \f$(kg c/m^2)\f$
    real, dimension(ilg,icc) :: glfltrdt !< brown leaf litter generated due to disturbance \f$(kg c/m^2)\f$
    real, dimension(ilg,icc) :: tltrstem   !< total stem litter fall rate (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, dimension(ilg,icc) :: tltrroot   !< total root litter fall rate (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, dimension(ilg,icc) :: leaflitr   !< leaf litter fall rate (\f$\mu mol CO2 m^{-2} s^{-1}\f$). this leaf litter
    !< does not include litter generated due to mortality/fire
    real, dimension(ilg,icc) :: roottemp   !< root temperature, k
    real, dimension(ilg,icc) :: afrleaf    !< allocation fraction for leaves
    real, dimension(ilg,icc) :: afrstem    !< allocation fraction for stem
    real, dimension(ilg,icc) :: afrroot    !< allocation fraction for root
    real, dimension(ilg,icc) :: wtstatus   !< soil water status used for calculating allocation fractions
    real, dimension(ilg,icc) :: ltstatus   !< light status used for calculating allocation fractions
    real, dimension(ilg)   :: rmr        !< root maintenance respiration (\f$\mu mol CO2 m^{-2} s^{-1}\f$)

    real, dimension(ilg,8) :: slopefrac      !< prescribed fraction of wetlands based on slope
    !< only(0.025,0.05,0.1,0.15,0.20,0.25,0.3 and 0.35 percent slope thresholds)
    real, dimension(ilg)   :: wetfrac_pres  !< Prescribed fraction of wetlands in a grid cell
    real, dimension(ilg)   :: ch4WetSpec       !< methane flux from wetlands calculated using hetrores (\f$\mu mol CH_4 m^{-2} s^{-1}\f$)
    real, dimension(ilg)   :: wetfdyn       !< dynamic wetland fraction
    real, dimension(ilg)   :: ch4WetDyn       !< methane flux from wetlands calculated using hetrores
    !< and wetfdyn, (\f$\mu mol CH_4 m^{-2} s^{-1}\f$)
    real, dimension(ilg)   :: ch4_soills    !< Methane uptake into the soil column (\f$mg CH_4 m^{-2} s^{-1}\f$)

    real, dimension(ilg) :: fNleach               !< grid averaged total N leaching \f$(kg N m^{-2} s^{-1}\f$
    real, dimension(ilg) :: fNvol                 !< grid averaged total N volatilization from ammonium pool (\f$g N/m^2 day\f$)
    real, dimension(ilg)   :: lucemcom   !< land use change (luc) related combustion emission losses, (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, dimension(ilg)   :: lucltrin   !< luc related inputs to litter pool, (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, dimension(ilg)   :: lucsocin   !< luc related inputs to soil c pool, (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, dimension(ilg)   :: cproduct   !< luc related product pool, sum of paper and furniture, (\f$\ Kg C m^{-2} \f$)
    real, dimension(ilg)   :: fproductdecomp  !< luc related decompsition from its product pools, sum of paper and furniture decomposition, (\f$\ mu mol CO2 m^{-2} s^{-1} \f$)
    real, dimension(ilg)   :: cleaf               !< Grid averaged leaf C pool , \f$kg C/m^2\f$
    real, dimension(ilg)   :: cstem               !< Grid averaged stem C pool , \f$kg C/m^2\f$
    real, dimension(ilg)   :: croot               !< Grid averaged root C pool , \f$kg C/m^2\f$
    real, dimension(ilg)   :: nppleaf             !< Grid averaged leaf NPP, \f$\mu mol CO2 m^{-2} s^{-1}\f$
    real, dimension(ilg)   :: nppstem             !< Grid averaged stem NPP, \f$\mu mol CO2 m^{-2} s^{-1}\f$
    real, dimension(ilg)   :: npproot             !< Grid averaged root NPP, \f$\mu mol CO2 m^{-2} s^{-1}\f$
    real, dimension(ilg)   :: nMineralNH4         !< Grid average NH4 mass (for CanESM), \f$kg N/m^2\f$
    real, dimension(ilg)   :: nMineralNO3         !< Grid average NO3 mass (for CanESM), \f$kg N/m^2\f$
    real, dimension(ilg)   :: nVeg                !< Grid average vegetation N mass (for CanESM), \f$kg N/m^2\f$
    real, dimension(ilg)   :: nLitter             !< Grid average litter N mass (for CanESM), \f$kg N/m^2\f$
    real, dimension(ilg)   :: nSoil               !< Grid average soil N mass mass (for CanESM),\f$kg N/m^2\f$
    real, dimension(ilg)   :: nLeaf               !< Grid average leaf N mass mass (for CanESM), \f$kg N/m^2\f$
    real, dimension(ilg)   :: nStem               !< Grid average stem N mass mass (for CanESM), \f$kg N/m^2\f$
    real, dimension(ilg)   :: nRoot               !< Grid average root N mass mass (for CanESM), \f$kg N/m^2\f$
    real, dimension(ilg)   :: fBNF                !< grid average total biological nitrogen fixation \f$(kg N m^{-2} s^{-1}\f$
    real, dimension(ilg)   :: fNdep               !< grid average total nitrogen deposition \f$(kg N m^{-2} s^{-1}\f$
    real, dimension(ilg)   :: fNfert              !< grid average total nitrogen fertilization \f$(kg N m^{-2} s^{-1}\f$
    real, dimension(ilg)   :: fNgasNonFire        !< grid average total nitrogen flux excluding fire \f$(kg N m^{-2} s^{-1}\f$
    real, dimension(ilg)   :: fNgasFire           !< grid average total nitrogen flux from fire only \f$(kg N m^{-2} s^{-1}\f$
    real, dimension(ilg)   :: fNgas               !< grid average total nitrogen flux from all sources including fire \f$(kg N m^{-2} s^{-1}\f$
    real, dimension(ilg)   :: fNnetmin            !< grid average total net mineralization \f$(kg N m^{-2} s^{-1}\f$
    real, dimension(ilg)   :: fNOx                !< grid average total NOx flux from all sources including fire \f$(kg N m^{-2} s^{-1}\f$
    real, dimension(ilg)   :: fNup                !< grid average pasive N updatake from the NH4 and NO3 pools \f$(kg N m^{-2} s^{-1}\f$
    real, dimension(ilg)   :: fNvegSoil           !< grid average total nitrogen flux from vegetation to soil [0 at present] \f$(kg N m^{-2} s^{-1}\f$
    real, dimension(ilg)   :: fN2o                !< grid average total N20 flux from all sources \f$(kg N m^{-2} s^{-1}\f$
    real, dimension(ilg)   :: lucemcomn   !< land use change (luc) related combustion emission losses of N (g N m^-2 day^-1)
    real, dimension(ilg)   :: lucltrinn   !< luc related inputs to litter pool of N (g N m^-2 day^-1)
    real, dimension(ilg)   :: lucsocinn   !< luc related inputs to soil N pool of N (g N m^-2 day^-1)

    real, dimension(ilg)   :: npp        !< net primary productivity
    real, dimension(ilg)   :: nep        !< net ecosystem productivity
    real, dimension(ilg)   :: nepCMIP        !< net ecosystem productivity
    real, dimension(ilg)   :: nbp        !< net biome productivity
    real, dimension(ilg)   :: gpp        !< gross primary productivity
    real, dimension(ilg)   :: hetrores   !< heterotrophic respiration
    real, dimension(ilg)   :: autores    !< autotrophic respiration
    real, dimension(ilg)   :: soilcresp  !<
    real, dimension(ilg)   :: rm         !< maintenance respiration
    real, dimension(ilg)   :: rg         !< growth respiration
    real, dimension(ilg)   :: litres     !< litter respiration
    real, dimension(ilg)   :: socres     !< soil carbon respiration
    real, dimension(ilg)   :: dstcemls   !< carbon emission losses due to disturbance, mainly fire
    real, dimension(ilg)   :: litrfall   !< total litter fall (from leaves, stem, and root) due to
    !< all causes (mortality,turnover,and disturbance)(\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg)   :: humiftrs   !< transfer of humidified litter from litter to soil c pool

    real, dimension(ilg,icc) :: gppveg     !< gross primary productity for each pft
    real, dimension(ilg,icc) :: vcmax0     !< max. photosynthetic rate at the top of canopy (\f$(mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg,icc) :: avesmfunc    !< soil mositure scalar on photosnthesis (uniless)
    real, dimension(ilg,iccp1) :: nepveg     !< net ecosystem productity for bare fraction expnbaln(i)=0.0 amount
    !< of c related to spatial expansion Not used JM Jun 2014
    !< OR net ecosystem productity for each pft

    real, dimension(ilg,iccp1) :: bnf_free         !< free-living biological nitrogen fixation (\f$g N/m^2day\f$)
    real, dimension(ilg,icc) :: bnf_ant          !< anthropogenic biological nitrogen fixation (\f$g N/m^2day\f$)
    real, dimension(ilg,icc) :: bnf_nat          !< natural biological nitrogen fixation (\f$g N/m^2day\f$)
    real, dimension(ilg,iccp1) :: bnf_tot          !< total biological nitrogen fixation (\f$g N/m^2day\f$)
    real, dimension(ilg)       :: bnf_moss       !< moss-related biological nitrogen fixation (\f$g N/m^2day\f$)
    real, dimension(ilg,icc) :: nstress          !< N stress
    real, dimension(ilg,iccp1) :: nitrifveg             !< nitrification for individual PFTs + bareground (\f$g N/m^2 day\f$)
    real, dimension(ilg,iccp1) :: no_nitveg             !< NO loss through nitrification (\f$g N/m^2 day\f$) for individual PFTs + bareground
    real, dimension(ilg,iccp1) :: no_denitveg           !< NO loss through denitrification (\f$g N/m^2 day\f$) for individual PFTs + bareground
    real, dimension(ilg,iccp1) :: no_nitdenitveg        !< total NO loss from denitrification and nitrification (\f$g N/m^2 day\f$) for individual PFTs + bareground
    real, dimension(ilg,iccp1) :: n2o_nitveg            !< N2O loss through nitrification (\f$g N/m^2 day\f$) for individual PFTs + bareground
    real, dimension(ilg,iccp1) :: n2o_denitveg          !< N2O loss through denitrification (\f$g N/m^2 day\f$) for individual PFTs + bareground
    real, dimension(ilg,iccp1) :: n2o_nitdenitveg       !< total N2O loss from denitrification and nitrification (\f$g N/m^2 day\f$) for individual PFTs + bareground
    real, dimension(ilg,iccp1) :: n2_denitveg           !< N2 loss through denitrification (\f$g N/m^2 day\f$) for individual PFTs + bareground
    real, dimension(ilg,iccp1) :: nvolveg               !< volatilization loss from ammonium pool (\f$g N/m^2 day\f$) for individual PFTs + bareground
    real, dimension(ilg,iccp1) :: nleachveg             !< leaching from nitrate pool (\f$g N/m^2 day\f$) for individual PFTs + bareground
    real, dimension(ilg,iccp1) :: appl_fert             !< applied nitrogen fertilizer \f$(g N m^{-2} cropland day^{-1}\f$
    real, dimension(ilg,iccp1) :: ndep_nh4              !< deposition influx into the Ammonium pool for individual PFTs + bareground \f$(g N m^{-2} day^{-1}\f$
    real, dimension(ilg,iccp1) :: ndep_no3              !< deposition influx into the Nitrate pool for individual PFTs + bareground \f$(g N m^{-2} day^{-1}\f$
    real, dimension(ilg,icc) :: ndemandveg_wp_npp     !< whole plant npp-based nitrogen demand for each of the 9 ctem pfts (\f$g N/m^2 day\f$)
    real, dimension(ilg,icc) :: nuptakeveg_p_nh4      !< passive nh4+ uptake for individual PFTs (\f$g N/m^2 day\f$)
    real, dimension(ilg,icc) :: nuptakeveg_p_no3      !< passive no3- uptake for individual PFTs (\f$g N/m^2 day\f$)
    real, dimension(ilg,icc) :: nuptakeveg_a_actl_nh4 !< actual active nh4+ uptake for individual PFTs (\f$g N/m^2 day\f$)
    real, dimension(ilg,icc) :: nuptakeveg_a_actl_no3 !< actual active no3- uptake for individual PFTs (\f$g N/m^2 day\f$)
    real, dimension(ilg,icc) :: nuptakeveg            !< total N uptake (active+passive, NH4+NO3) for individual PFTs (\f$g N/m^2 day\f$)
    real, dimension(ilg,icc) :: nleafns2sveg          !< nitrogen flux from non-structural to structural leaf pool (\f$g N/m^2 day\f$)
    real, dimension(ilg,icc) :: nstemns2sveg          !< nitrogen flux from non-structural to structural stem pool (\f$g N/m^2 day\f$)
    real, dimension(ilg,icc) :: nrootns2sveg          !< nitrogen flux from non-structural to structural root pool (\f$g N/m^2 day\f$)
    real, dimension(ilg,icc) :: nallocveg_l           !< nitrogen allocation to leaves for individual PFTs (\f$g N/m^2 day\f$)
    real, dimension(ilg,icc) :: nallocveg_s           !< nitrogen allocation to stem for individual PFTs (\f$g N/m^2 day\f$)
    real, dimension(ilg,icc) :: nallocveg_r           !< nitrogen allocation to root for individual PFTs (\f$g N/m^2 day\f$)
    real, dimension(ilg,icc) :: nresorpedveg_s        !< resorped N from leaves to be allocated to stem (\f$g N/m^2 day\f$)
    real, dimension(ilg,icc) :: nresorpedveg_r        !< resorped N from leaves to be allocated to root (\f$g N/m^2 day\f$)
    real, dimension(ilg,icc) :: nre_allocveg_s2l      !< reallocated N from S to L during leaf out period (\f$g N/m^2 day\f$)
    real, dimension(ilg,icc) :: nre_allocveg_r2l      !< reallocated N from R to L during leaf out period (\f$g N/m^2 day\f$)
    real, dimension(ilg,icc) :: nlitrveg_l            !< leaf N litterfall (\f$g N/m^2 day\f$)
    real, dimension(ilg,icc) :: nlitrveg_s            !< stem N litterfall (\f$g N/m^2 day\f$)
    real, dimension(ilg,icc) :: nlitrveg_r            !< root N litterfall (\f$g N/m^2 day\f$)
    real, dimension(ilg,icc) :: nlitrveg              !< total N litterfall (\f$g N/m^2 day\f$)
    real, dimension(ilg,icc) :: gl2bl_grass_nflux     !< N flux from ngleafmas to nbleafmas (\f$g N/m^2 day\f$)
    real, dimension(ilg,icc) :: c2nveg_l              !< simulated C:N ratio for leaves (\f$g C/g N\f$)
    real, dimension(ilg,icc) :: c2nveg_s              !< simulated C:N ratio for stem (\f$g C/g N\f$)
    real, dimension(ilg,icc) :: c2nveg_r              !< simulated C:N ratio for roots (\f$g C/g N\f$)
    real, dimension(ilg,icc) :: c2nveg_wp             !< simulated C:N ratio for the whole plant (\f$g C/g N\f$)
    real, dimension(ilg,iccp1) :: c2nveg_litr           !< simulated C:N ratio for litter mass (\f$g C/g N\f$)
    real, dimension(ilg,iccp1) :: c2nveg_humus          !< simulated C:N ratio for humus mass (\f$g C/g N\f$)
    real, dimension(ilg,iccp1) :: nhumtrsveg            !< N humification for individual PFTs + bare (\f$g N/m^2 day\f$)
    real, dimension(ilg,iccp1) :: nmineralveg_litr      !< N mineralization from litter pool (\f$g N/m^2 day\f$)
    real, dimension(ilg,iccp1) :: nmineralveg_humus     !< N mineralization from organic soil pool (\f$g N/m^2 day\f$)
    real, dimension(ilg,iccp1) :: netnmineralveg        !< N mineralization (\f$g N/m^2 day\f$)
    real, dimension(ilg,iccp1) :: nimmobilveg_nh4       !< N immobilization from NH4+ pool to soilnmas (\f$g N/m^2 day\f$)
    real, dimension(ilg,iccp1) :: nimmobilveg_no3       !< N immobilization from NO3- pool to soilnmas (\f$g N/m^2 day\f$)
    real, dimension(ilg,iccp1) :: fNnetlandveg          !< net terrestrial N flux (\f$g N/m^2 day\f$)
    real, dimension(ilg,icc) :: redcoeff_vcmax        !< Reduction coeff. passed to the Photosynthesis subroutine

    character(8), dimension(ilg) :: peatlandType !< Peatland type (else 'None')
    real, dimension(ilg) :: peatdep      !< Depth of peat column (m)
    real, dimension(ilg) :: peatSoilC  !< peat soil C mass, \f$kg C/m^2\f$
    real, dimension(ilg) :: redcoeff_vcmaxMoss      !< Reduction coeff. for moss passed to the Photosynthesis subroutine
    character(len=8), dimension(ilg) :: mossPresent !< Type of moss present (else 'None') 
    real, dimension(ilg) :: nppmoss    !< net primary production of moss (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, dimension(ilg) :: armoss     !< autotrophic respiration of moss (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, dimension(ilg,ignd) :: litrmsmoss  !< moss litter mass, \f$kg C/m^2\f$
    real, dimension(ilg,ignd) :: upMossSoilC !< Moss soil C mass (non-peat), \f$kg C/m^2\f$
    real, dimension(ilg) :: Cmossmas   !< C in moss biomass, \f$kg C/m^2\f$

    real, dimension(ilg) :: litrmsmossN  !< Moss litter N mass, \f$g N/m^2\f$
    real, dimension(ilg) :: upMossSoilN  !< Moss soil C mass (non-peat), \f$g N/m^2\f$
    real, dimension(ilg) :: Nmossmas     !< N in moss biomass, \f$g N/m^2\f$

    real, dimension(ilg) :: dmoss      !< depth of living moss (m)
    real, dimension(ilg) :: ancsmoss   !< moss net photosynthesis in canopy snow subarea (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, dimension(ilg) :: angsmoss   !< moss net photosynthesis in snow ground subarea (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, dimension(ilg) :: ancmoss    !< moss net photosynthesis in canopy ground subarea (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, dimension(ilg) :: angmoss    !< moss net photosynthesis in bare ground subarea (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, dimension(ilg) :: rmlcsmoss  !< moss maintenance respiration in canopy snow subarea (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, dimension(ilg) :: rmlgsmoss  !< moss maintenance respiration in ground snow subarea (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, dimension(ilg) :: rmlcmoss   !< moss maintenance respiration in canopy ground subarea (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, dimension(ilg) :: rmlgmoss   !< moss maintenance respiration in bare ground subarea (\f$\mu mol CO2 m^{-2} s^{-1}\f$)

    real, dimension(ilg,iccp1) :: nbpveg     !< net biome productity for bare fraction OR net biome productity for each pft
    real, dimension(ilg,icc) :: nppveg     !< npp for individual pfts, (\f$\mu mol CO2 m^{-2} s^{-1}\f$)
    real, dimension(ilg,iccp1) :: hetroresveg !<
    real, dimension(ilg,icc) :: autoresveg !<
    real, dimension(ilg,iccp2,ignd) :: litresveg  !<
    real, dimension(ilg,iccp2,ignd) :: soilcresveg !<
    real, dimension(ilg,iccp2,ignd) :: humiftrsveg !<
    real, dimension(ilg,icc) :: rmlvegacc  !<
    real, dimension(ilg,icc) :: rmsveg     !< stem maintenance resp. rate for each pft
    real, dimension(ilg,icc) :: rmrveg     !< root maintenance resp. rate for each pft
    real, dimension(ilg,icc) :: rgveg      !< growth resp. rate for each pft
    real, dimension(ilg,icc) :: litrfallveg !< litter fall for each pft (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg,icc) :: reprocost   !< Cost of making reproductive tissues, only non-zero when NPP is positive (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)

    real, dimension(ilg,icc) :: rothrlos !< root death as crops are harvested, \f$kg c/m^2/d\f$
    real, dimension(ilg,icc) :: pfcancmx !< previous year's fractional coverages of pfts
    real, dimension(ilg,icc) :: nfcancmx !< next year's fractional coverages of pfts
    real, dimension(ilg,ican) :: alvsctm  !<
    real, dimension(ilg,ican) :: paic     !< plant area index for class' 4 pfts. this is the sum of leaf
    !< area index and stem area index.
    real, dimension(ilg,ican) :: slaic    !< storage lai. this will be used as min. lai that class sees
    !< so that it doesn't blow up in its stomatal conductance calculations.
    real, dimension(ilg,ican) :: alirctm  !<
    real, dimension(ilg)   :: cfluxcg  !<
    real, dimension(ilg)   :: cfluxcs  !<
    real, dimension(ilg)   :: CFLUX_GA !< Product of surface drag coefficient and wind speed \f$[m s^{-1} ]\f$
    real, dimension(ilg)   :: USTARBS_GA !< Friction velocity to be used in nitrogen volatilization  \f$[m s^{-1} ]\f$
    real, dimension(ilg)   :: ROFB     !< Base flow from bottom of soil column \f$[kg m^{-2} s^{-1} ]\f$
    real, dimension(ilg)   :: dstcemls3 !< carbon emission losses due to disturbance (fire at present) from litter pool
    real, dimension(ilg,icc) :: anveg    !< net photosynthesis rate for each pft
    real, dimension(ilg,icc) :: rmlveg   !< leaf maintenance resp. rate for each pft

    real, dimension(ilg) :: twarmm            !< temperature of the warmest month (c)
    real, dimension(ilg) :: tcoldm            !< temperature of the coldest month (c)
    real, dimension(ilg) :: gdd5              !< growing degree days above 5 c
    real, dimension(ilg) :: aridity           !< aridity index, ratio of potential evaporation to precipitation
    real, dimension(ilg) :: srplsmon          !< number of months in a year with surplus water i.e. precipitation more than potential evaporation
    real, dimension(ilg) :: defctmon          !< number of months in a year with water deficit i.e. precipitation less than potential evaporation
    real, dimension(ilg) :: anndefct          !< annual water deficit (mm)
    real, dimension(ilg) :: annsrpls          !< annual water surplus (mm)
    real, dimension(ilg) :: annpcp            !< annual precipitation (mm)
    real, dimension(ilg) :: dry_season_length !< length of dry season (months)
    integer, dimension(ilg) :: colddays_leaffall !< cold days counter for tracking days below a certain
    !< temperature threshold for ndl dcd.
    integer, dimension(ilg) :: colddays_harvest !< cold days counter for tracking days below a certain
    !< temperature threshold for crop.

    ! These go into CTEM and are used to keep track of the bioclim limits.
    real, dimension(ilg) :: tcurm     !< temperature of the current month (c)
    real, dimension(ilg) :: srpcuryr  !< water surplus for the current year
    real, dimension(ilg) :: dftcuryr  !< water deficit for the current year
    real, dimension(12,ilg) :: tmonth  !< monthly temperatures
    real, dimension(ilg) :: anpcpcur  !< annual precipitation for current year (mm)
    real, dimension(ilg) :: anpecur   !< annual potential evaporation for current year (mm)
    real, dimension(ilg) :: gdd5cur   !< growing degree days above 5 c for current year
    real, dimension(ilg) :: surmncur  !< number of months with surplus water for current year
    real, dimension(ilg) :: defmncur  !< number of months with water deficit for current year
    real, dimension(ilg) :: srplscur  !< water surplus for the current month
    real, dimension(ilg) :: defctcur  !< water deficit for the current month

    real, dimension(ilg,icc) :: geremort !< growth efficiency related mortality (1/day)
    real, dimension(ilg,icc) :: intrmort !< intrinsic (age related) mortality (1/day)
    real, dimension(ilg,icc) :: cc       !< colonization rate & mortality rate
    real, dimension(ilg,icc) :: mm       !< colonization rate & mortality rate

    logical, dimension(ilg,icc) :: pftexist !< logical array indicating pfts exist (t) or not (f)
    integer, dimension(ilg,icc) :: lfstatus !< leaf phenology status
    integer, dimension(ilg,icc) :: pandays  !< days with positive net photosynthesis (an) for use in
    !< the phenology subroutine
    real, dimension(ilg) :: grclarea      !< area of the grid cell, \f$km^2\f$

    integer, dimension(ilg) :: altotcount_ctm ! nlat  !< Counter used for calculating total albedo
    real, dimension(ilg,icc)  :: todfrac  !(ilg,icc)   !< Max. fractional coverage of ctem's 9 pfts by the end of the day, for use by land use subroutine
    real, dimension(ilg)    :: fsinacc_gat !(ilg)    !<
    real, dimension(ilg)    :: flutacc_gat !(ilg)    !<
    real, dimension(ilg)    :: flinacc_gat !(ilg)    !<
    ! real, dimension(ilg)    :: pregacc_gat !(ilg)    !<
    real, dimension(ilg)    :: altotacc_gat !(ilg)   !<
    real, dimension(ilg)    :: netrad_gat !(ilg)     !<
    real, dimension(ilg)    :: preacc_gat !(ilg)     !<
    real, dimension(ilg)    :: sdepgat !(ilg)        !<
    !     real, dimension(ilg,ignd)  :: rgmgat !(ilg,ignd)    !<
    real, dimension(ilg,ignd)  :: sandgat !(ilg,ignd)   !<
    real, dimension(ilg)    :: xdiffusgat !(ilg)
    real, dimension(ilg)    :: faregat !(ilg)

    real, dimension(ilg,icc) :: lygleafmasmax       !< last year maximum of the gleafmas
    real, dimension(ilg,icc) :: lystemmassmax       !< last year maximum of the stemmass
    real, dimension(ilg,icc) :: lyrootmassmax       !< last year maximum of the rootmass
    real, dimension(ilg,icc) :: lmaxt, smaxt, rmaxt !< temp vars to find previous year C pool

  end type veg_gat

  !=================================================================================
  type tracersType(nlat, nmos, ignd, ilg)
    !   Simple tracer variables. Only written to if useTracer > 0.

    ! NOTE: Units may vary depending on the tracer used, see convertTracerUnits in tracer.f90

    integer, len :: nlat
    integer, len :: nmos
    integer, len :: ignd
    integer, len :: ilg
    
    ! Pools:
    ! allocated with nlat, nmos, ...:
    real, dimension(nlat,nmos) :: mossCMassrot      !< Tracer mass in moss biomass, 14C: \f$ng ^{14}C/m^2\f$
    real, dimension(nlat,nmos) :: mossLitrMassrot   !< Tracer mass in moss litter, 14C: \f$ng ^{14}C/m^2\f$
    real, dimension(nlat,nmos) :: tracerCO2rot     !< Atmopspheric tracer CO2 concentration (units vary)

    ! allocated with nlat, nmos, icc:
    real, dimension(nlat,nmos,icc) :: gLeafMassrot      !< Tracer mass in the green leaf pool for each of the CTEM pfts, 14C: \f$ng ^{14}C/m^2\f$
    real, dimension(nlat,nmos,icc) :: bLeafMassrot      !< Tracer mass in the brown leaf pool for each of the CTEM pfts, 14C: \f$ng ^{14}C/m^2\f$
    real, dimension(nlat,nmos,icc) :: stemMassrot       !< Tracer mass in the stem for each of the CTEM pfts, 14C: \f$ng ^{14}C/m^2\f$
    real, dimension(nlat,nmos,icc) :: rootMassrot       !< Tracer mass in the roots for each of the CTEM pfts, 14C: \f$ng ^{14}C/m^2\f$
    real, dimension(nlat,nmos,icc) :: rothrlosrot       !< Tracer root death as crops are harvested, 14C: \f$ng ^{14}C/m^2\f$
    real, dimension(nlat,nmos,icc) :: stmhrlosrot       !< Tracer stem harvest loss as crops are harvested, 14C: \f$ng ^{14}C/m^2\f$
    real, dimension(nlat,nmos,icc) :: FlHrLossrot       !< Tracer fall & harvest loss for bdl dcd plants and crops, respectively,14C: \f$ng ^{14}C/m^2\f$
    
    ! allocated with nlat, nmos, iccp2, ignd:
    real, dimension(nlat,nmos,iccp2,ignd) :: litrMassrot       !< Tracer mass in the litter pool for each of the CTEM pfts + bareground and LUC products, 14C: \f$ng ^{14}C/m^2\f$
    real, dimension(nlat,nmos,iccp2,ignd) :: soilCMassrot      !< Tracer mass in the soil carbon pool for each of the CTEM pfts + bareground and LUC products, 14C: \f$ng ^{14}C/m^2\f$

    ! allocated with ilg, ...:
    real, dimension(ilg) :: mossCMassgat      !< Tracer mass in moss biomass, 14C: \f$ng ^{14}C/m^2\f$
    real, dimension(ilg) :: mossLitrMassgat   !< Tracer mass in moss litter, 14C: \f$ng ^{14}C/m^2\f$
    real, dimension(ilg) :: tracerCO2gat     !< Atmopspheric tracer CO2 concentration (14C: 1E12 * 14C/C)

    ! allocated with nlat, nmos, icc:
    real, dimension(ilg,icc) :: gLeafMassgat      !< Tracer mass in the green leaf pool for each of the CTEM pfts, 14C: \f$ng ^{14}C/m^2\f$
    real, dimension(ilg,icc) :: bLeafMassgat      !< Tracer mass in the brown leaf pool for each of the CTEM pfts, 14C: \f$ng ^{14}C/m^2\f$
    real, dimension(ilg,icc) :: stemMassgat       !< Tracer mass in the stem for each of the CTEM pfts, 14C: \f$ng ^{14}C/m^2\f$
    real, dimension(ilg,icc) :: rootMassgat       !< Tracer mass in the roots for each of the CTEM pfts, 14C: \f$ng ^{14}C/m^2\f$
    real, dimension(ilg,icc) :: rothrlosgat       !< Tracer root death as crops are harvested, 14C: \f$ng ^{14}C/m^2\f$
    real, dimension(ilg,icc) :: stmhrlosgat       !< Tracer stem harvest loss as crops are harvested, 14C: \f$ng ^{14}C/m^2\f$
    real, dimension(ilg,icc) :: FlHrLossgat       !< Tracer fall & harvest loss for bdl dcd plants and crops, respectively, 14C: \f$ng ^{14}C/m^2\f$

    ! allocated with nlat, nmos, iccp2, ignd:
    real, dimension(ilg,iccp2,ignd) :: litrMassgat       !< Tracer mass in the litter pool for each of the CTEM pfts + bareground and LUC products, 14C: \f$ng ^{14}C/m^2\f$
    real, dimension(ilg,iccp2,ignd) :: soilCMassgat      !< Tracer mass in the soil carbon pool for each of the CTEM pfts + bareground and LUC products, 14C: \f$ng ^{14}C/m^2\f$

  end type tracersType

  !=================================================================================
  !> CTEM's variables per tile
  type ctem_tile_level(ilg, ignd)

    integer, len :: ilg
    integer, len :: ignd
    
    !   Tile-level variables (denoted by an ending of "_t")

    real, dimension(ilg) :: fsnowacc_t       !<
    real, dimension(ilg) :: CFLUX_GAacc_t    !<daily accu. boundary layer aerodynamic conductance \f$[m day^{-1} ]\f$
    real, dimension(ilg) :: USTARBS_GAacc_t  !<daily accu. friction velocity to be used in nitrogen volatilization  \f$[m day^{-1} ]\f$
    real, dimension(ilg) :: ROFBacc_t        !<daily accu. base flow from bottom of soil column \f$[kg m^{-2} day^{-1} ]\f$
    real, dimension(ilg) :: tcansacc_t       !<
    real, dimension(ilg) :: tcanoaccgat_t    !<
    real, dimension(ilg) :: taaccgat_t       !<
    real, dimension(ilg) :: uvaccgat_t       !<
    real, dimension(ilg) :: vvaccgat_t       !<
    real, dimension(ilg) :: anmossac_t       !< daily averaged moss net photosynthesis accumulated (\f$\mu mol /m^2 /s\f$)
    real, dimension(ilg) :: rmlmossac_t      !< daily averaged moss maintainence respiration (\f$\mu mol /m^2 /s\f$)
    real, dimension(ilg) :: gppmossac_t      !< daily averaged gross primary production (\f$\mu mol /m^2 /s\f$)

    ! allocated with ilg, ignd:
    real, dimension(ilg,ignd) :: QFCacc_t    !< daily accu. water removed from soil layers by transpiration \f$[kg m^{-2} day^{-1}]\f$
    real, dimension(ilg,ignd) :: tbaraccgat_t !<
    real, dimension(ilg,ignd) :: thliqacc_t  !<
    real, dimension(ilg,ignd) :: thiceacc_t  !< Added in place of YW's thicaccgat_m. EC Dec 23 2016.

    ! allocated with ilg, icc:
    real, dimension(ilg,icc) :: ancgvgac_t  !<
    real, dimension(ilg,icc) :: rmlcgvga_t  !<

    real, dimension(ilg) :: NoonTempc           !< area of the grid cell, \f$km^2\f$
    real, dimension(ilg) :: NoonPrecip            !< centre latitude of grid cells in degrees
    real, dimension(ilg) :: NoonPreciprun            !< centre latitude of grid cells in degrees
    real, dimension(ilg) :: NoonWindtot               !< u wind speed, m/s
    real, dimension(ilg) :: NoonRhum               !< v wind speed, m/s
    real, dimension(ilg) :: NoonTempcrun

  end type ctem_tile_level

contains

  ! -----------------------------------------------------

  !> \ingroup ctemstatevars_initRowVarsBioGeoChem
  
  !> Initializes 'row' variables
  subroutine initRowVarsBioGeoChem(vrot, nlat, nmos, ignd)

    implicit none
    
    integer, intent(in) :: nmos
    integer, intent(in) :: nlat
    integer, intent(in) :: ignd
    class(veg_rot(nlat, nmos, ignd)), intent(inout) :: vrot

    vrot%fNleach = 0.0
    vrot%fNvol = 0.0
    vrot%co2conc = 0.0
    vrot%npp = 0.0
    vrot%nep = 0.0
    vrot%nepCMIP = 0.0
    vrot%hetrores = 0.0
    vrot%autores = 0.0
    vrot%soilcresp = 0.0
    vrot%rm = 0.0
    vrot%rg = 0.0
    vrot%nbp = 0.0
    vrot%litres = 0.0
    vrot%socres = 0.0
    vrot%gpp = 0.0
    vrot%dstcemls = 0.0
    vrot%dstcemls3 = 0.0
    vrot%litrfall = 0.0
    vrot%humiftrs = 0.0
    vrot%canres = 0.0
    vrot%rml = 0.0
    vrot%rms = 0.0
    vrot%rmr = 0.0
    vrot%lucemcom = 0.0
    vrot%lucltrin = 0.0
    vrot%lucsocin = 0.0
    vrot%cproduct = 0.0
    vrot%fproductdecomp = 0.0
    vrot%cleaf = 0.0
    vrot%cstem = 0.0
    vrot%croot = 0.0
    vrot%nppleaf = 0.0
    vrot%nppstem = 0.0
    vrot%npproot = 0.0
    vrot%nMineralNH4 = 0.0
    vrot%nMineralNO3 = 0.0
    vrot%nVeg = 0.0
    vrot%nLitter = 0.0
    vrot%nSoil = 0.0
    vrot%nLeaf = 0.0
    vrot%nStem = 0.0
    vrot%nRoot = 0.0
    vrot%fBNF = 0.0
    vrot%fNdep = 0.0
    vrot%fNfert = 0.0
    vrot%fNgasNonFire = 0.0
    vrot%fNgasFire = 0.0
    vrot%fNgas = 0.0
    vrot%fNnetmin = 0.0
    vrot%fNOx = 0.0
    vrot%fNup = 0.0
    vrot%fNvegSoil = 0.0
    vrot%fN2o = 0.0
    vrot%lucemcomn = 0.0
    vrot%lucltrinn = 0.0
    vrot%lucsocinn = 0.0
    vrot%burnfrac = 0.0
    vrot%lterm = 0.0
    !vrot%cfluxcg = 0.0
    !vrot%cfluxcs = 0.0
    vrot%ROFB = 0.0
    vrot%ch4WetSpec = 0.0
    vrot%wetfdyn = 0.0
    vrot%ch4WetDyn = 0.0
    vrot%ch4_soills = 0.0
    vrot%nppmoss = 0.0
    vrot%armoss = 0.0
    vrot%peatdep = 0.0
    vrot%peatSoilC = 0.0
    vrot%ZOLNC = 0.0
    vrot%AILC = 0.0
    vrot%CMASVEGC = 0.0
    vrot%ALVSCTM = 0.0
    vrot%ALIRCTM = 0.0
    vrot%PAIC = 0.0
    vrot%SLAIC = 0.0
    vrot%RMATC = 0.0
    vrot%smfuncveg = 0.0
    vrot%drgtstrs = 0.0
    vrot%betadrgt = 0.0

    vrot%SRHrow = 0.0
    vrot%isirow = 0.0
    vrot%buirow = 0.0
    vrot%fwirow = 0.0
    vrot%ffmcrow = 0.0
    vrot%dmcrow = 0.0
    vrot%dcrow = 0.0
    vrot%gleafmas = 0.0
    vrot%gleafmas_ns = 0.0
    vrot%gleafmas_s = 0.0
    vrot%bleafmas = 0.0
    vrot%stemmass = 0.0
    vrot%stemmass_ns = 0.0
    vrot%stemmass_s = 0.0
    vrot%rootmass = 0.0
    vrot%rootmass_ns = 0.0
    vrot%rootmass_s = 0.0
    vrot%pstemmass = 0.0
    vrot%pgleafmass = 0.0
    vrot%litrfallveg = 0.0
    vrot%bterm = 0.0
    vrot%mterm = 0.0
    vrot%ailcg = 0.0
    vrot%ailcgs = 0.0
    vrot%fcancs = 0.0
    vrot%fcanc = 0.0
    vrot%fcancmx = 0.0
    !vrot%co2i1cg = 0.0
    !vrot%co2i1cs = 0.0
    !vrot%co2i2cg = 0.0
    !vrot%co2i2cs = 0.0
    vrot%ancsveg = 0.0
    vrot%ancgveg = 0.0
    vrot%rmlcsveg = 0.0
    vrot%rmlcgveg = 0.0
    vrot%ailcb = 0.0
    vrot%grwtheff = 0.0
    vrot%bmasveg = 0.0
    vrot%tltrleaf = 0.0
    vrot%tltrstem = 0.0
    vrot%tltrroot = 0.0
    vrot%leaflitr = 0.0
    vrot%roottemp = 0.0
    vrot%afrleaf = 0.0
    vrot%afrstem = 0.0
    vrot%afrroot = 0.0
    vrot%wtstatus = 0.0
    vrot%ltstatus = 0.0
    vrot%pfcancmx = 0.0
    vrot%nfcancmx = 0.0
    vrot%nppveg = 0.0
    vrot%veghght = 0.0
    vrot%rootdpth = 0.0
    vrot%anveg = 0.0
    vrot%rmlveg = 0.0
    vrot%rmlvegacc = 0.0
    vrot%rmsveg = 0.0
    vrot%rmrveg = 0.0
    vrot%rgveg = 0.0
    vrot%vgbiomas_veg = 0.0
    vrot%gppveg = 0.0
    vrot%vcmax0 = 0.0
    vrot%avesmfunc = 0.0
    vrot%autoresveg = 0.0
    vrot%emit_co2 = 0.0
    vrot%emit_co = 0.0
    vrot%emit_ch4 = 0.0
    vrot%emit_nmhc = 0.0
    vrot%emit_h2 = 0.0
    vrot%emit_nox = 0.0
    vrot%emit_n2o = 0.0
    vrot%emit_nh3 = 0.0
    vrot%emit_pm25 = 0.0
    vrot%emit_tpm = 0.0
    vrot%emit_tc = 0.0
    vrot%emit_bc = 0.0
    vrot%emit_oc = 0.0
    vrot%emit_so2 = 0.0
    vrot%burnvegf = 0.0
    vrot%fFireCVeg   = 0.0
    vrot%fFireLitter = 0.0
    vrot%fFireCSoil  = 0.0
    vrot%leafns2s = 0.0
    vrot%stemns2s = 0.0
    vrot%rootns2s = 0.0
    vrot%re_alloc_s2l = 0.0
    vrot%re_alloc_r2l = 0.0
    vrot%re_alloc_sr2l = 0.0
    vrot%ngleafmas = 0.0
    vrot%ngleafmas_ns = 0.0
    vrot%ngleafmas_s = 0.0
    vrot%nbleafmas = 0.0
    vrot%nstemmass = 0.0
    vrot%nstemmass_ns = 0.0
    vrot%nstemmass_s = 0.0
    vrot%nrootmass = 0.0
    vrot%nrootmass_ns = 0.0
    vrot%nrootmass_s = 0.0
    vrot%nvgbiomas_veg = 0.0
    vrot%ndemandveg_wp_npp = 0.0
    vrot%nuptakeveg_p_nh4 = 0.0
    vrot%nuptakeveg_p_no3 = 0.0
    vrot%nuptakeveg_a_actl_nh4 = 0.0
    vrot%nuptakeveg_a_actl_no3 = 0.0
    vrot%nuptakeveg = 0.0
    vrot%nallocveg_l = 0.0
    vrot%nallocveg_s = 0.0
    vrot%nallocveg_r = 0.0
    vrot%nresorpedveg_s = 0.0
    vrot%nresorpedveg_r = 0.0
    vrot%nre_allocveg_s2l = 0.0
    vrot%nre_allocveg_r2l = 0.0
    vrot%nleafns2sveg = 0.0
    vrot%nstemns2sveg = 0.0
    vrot%nrootns2sveg = 0.0
    vrot%gl2bl_grass_nflux = 0.0
    vrot%nlitrveg_l = 0.0
    vrot%nlitrveg_s = 0.0
    vrot%nlitrveg_r = 0.0
    vrot%nlitrveg   = 0.0
    vrot%c2nveg_l = 0.0
    vrot%c2nveg_s = 0.0
    vrot%c2nveg_r = 0.0
    vrot%c2nveg_wp = 0.0
    vrot%rmatctem = 0.0
    vrot%hetroresveg = 0.0
    vrot%nepveg = 0.0
    vrot%nbpveg = 0.0
    vrot%litrmass = 0.0
    vrot%soilcmas = 0.0
    vrot%litresveg = 0.0
    vrot%soilcresveg = 0.0
    vrot%humiftrsveg = 0.0
    vrot%bnf_free = 0.0
    vrot%bnf_ant = 0.0
    vrot%bnf_nat = 0.0
    vrot%bnf_tot = 0.0
    vrot%bnf_moss = 0.0
    vrot%nstress = 0.0
    vrot%nitrifveg = 0.0
    vrot%no_nitveg = 0.0
    vrot%no_denitveg = 0.0
    vrot%no_nitdenitveg = 0.0
    vrot%n2o_nitveg = 0.0
    vrot%n2o_denitveg = 0.0
    vrot%n2o_nitdenitveg = 0.0
    vrot%n2_denitveg = 0.0
    vrot%nvolveg = 0.0
    vrot%nleachveg = 0.0
    vrot%nh4_mass = 0.0
    vrot%no3_mass = 0.0
    vrot%nlitrmass = 0.0
    vrot%soilnmas = 0.0
    vrot%appl_fert = 0.0
    vrot%ndep_nh4 = 0.0
    vrot%ndep_no3 = 0.0
    vrot%c2nveg_litr = 0.0
    vrot%c2nveg_humus = 0.0
    vrot%nhumtrsveg = 0.0
    vrot%nmineralveg_litr = 0.0
    vrot%nmineralveg_humus = 0.0
    vrot%nimmobilveg_nh4 = 0.0
    vrot%nimmobilveg_no3 = 0.0
    vrot%netnmineralveg  = 0.0
    vrot%fNnetlandveg = 0.0
    vrot%redcoeff_vcmax = 0.0
    vrot%redcoeff_vcmaxMoss = 0.0
    vrot%twarmm = 0.0
    vrot%tcoldm = 0.0
    vrot%gdd5 = 0.0
    vrot%aridity = 0.0
    vrot%srplsmon = 0.0
    vrot%defctmon = 0.0
    vrot%anndefct = 0.0
    vrot%annsrpls = 0.0
    vrot%annpcp = 0.0
    vrot%dry_season_length = 0.0
    vrot%tileAgerow = 0.0

  end subroutine initRowVarsBioGeoChem
  

  
  !==================================================
  !> \ingroup ctemstatevars_resetMosaicAccum
  
  !> Resets physics accumulator variables (used as input to CTEM) after CTEM has been called
  subroutine resetMosaicAccum(vgat, ctem_tile, ilg, ignd)

    implicit none
    
    integer, intent(in) :: ilg
    integer, intent(in) :: ignd
    type(veg_gat(ilg, ignd)),         intent(inout) :: vgat
    type(ctem_tile_level(ilg, ignd)), intent(inout) :: ctem_tile

    vgat%fsinacc_gat(:) = 0.
    vgat%flinacc_gat(:) = 0.
    vgat%flutacc_gat(:) = 0.
    vgat%preacc_gat(:) = 0.
    ctem_tile%fsnowacc_t(:) = 0.0      ! daily accu. fraction of snow
    ctem_tile%taaccgat_t(:) = 0.0
    ctem_tile%vvaccgat_t(:) = 0.0
    ctem_tile%uvaccgat_t(:) = 0.0
    ctem_tile%CFLUX_GAacc_t(:) = 0.0   ! daily accu. boundary layer aerodynamic conductance \f$[m day^{-1} ]\f$
    ctem_tile%USTARBS_GAacc_t(:) = 0.0 ! daily accu. friction velocity to be used in nitrogen volatilization  \f$[m day^{-1} ]\f$
    ctem_tile%ROFBacc_t(:) = 0.0       ! daily accu. base flow from bottom of soil column \f$[kg m^{-2} day^{-1} ]\f$
    vgat%altotacc_gat(:) = 0.0
    vgat%altotcount_ctm(:) = 0
    ctem_tile%QFCacc_t(:,:) = 0.0
    ctem_tile%tbaraccgat_t(:,:) = 0.0
    ctem_tile%thliqacc_t(:,:) = 0.0
    ctem_tile%thiceacc_t(:,:) = 0.0

    ctem_tile%ancgvgac_t(:,:) = 0.0    ! daily accu. net photosyn.
    ctem_tile%rmlcgvga_t(:,:) = 0.0    ! daily accu. leaf respiration

    !-reset the accumulators for moss daily C fluxes-------------------------------
    ctem_tile%anmossac_t(:)  = 0.0
    ctem_tile%rmlmossac_t(:) = 0.0
    ctem_tile%gppmossac_t(:) = 0.0

  end subroutine resetMosaicAccum
  
  !=================================================================================

  !==================================================
  !> \ingroup ctemstatevars_ctemdump
  
  !> dumps all the ctemdump statevars for diagnostic purposes.
  subroutine ctemdump(vrot, vgat, tracer, ctem_tile)

    use classicParams,  only : c_switch, nlat, nmos, ignd, ilg
  
    implicit none

    type(veg_rot(nlat, nmos, ignd)),          intent(inout) :: vrot
    type(veg_gat(ilg, ignd)),                 intent(inout) :: vgat
    type(tracersType(nlat, nmos, ignd, ilg)), intent(inout) :: tracer
    type(ctem_tile_level(ilg, ignd)),         intent(inout) :: ctem_tile
    
    print *, 'begin ctem dump'
    print *, 'projectedGrid ', c_switch%projectedGrid
    print *, 'ctem_on ', c_switch%ctem_on
    print *, 'Ncycle_on ', c_switch%Ncycle_on
    print *, 'metLoop ', c_switch%metLoop
    print *, 'leap ', c_switch%leap
    print *, 'spinfast ', c_switch%spinfast
    print *, 'useTracer ', c_switch%useTracer
    print *, 'tracerCO2file ', c_switch%tracerCO2file
    print *, 'tracertimeAdjust', c_switch%tracertimeAdjust
    print *, 'actualMetStartYear ', c_switch%actualMetStartYear
    print *, 'actualMetEndYear ', c_switch%actualMetEndYear
    print *, 'runStartYear ', c_switch%runStartYear
    print *, 'runEndYear ', c_switch%runEndYear
    print *, 'metOrder ', c_switch%metOrder
    print *, 'transientCO2 ', c_switch%transientCO2
    print *, 'CO2File ', c_switch%CO2File
    print *, 'co2timeAdjust', c_switch%co2timeAdjust
    print *, 'fixedYearCO2 ', c_switch%fixedYearCO2
    print *, 'doMethane ', c_switch%doMethane
    print *, 'transientCH4 ', c_switch%transientCH4
    print *, 'CH4File ', c_switch%CH4File
    print *, 'fixedYearCH4 ', c_switch%fixedYearCH4
    print *, 'dofire ', c_switch%dofire
    print *, 'transientPOPD ', c_switch%transientPOPD
    print *, 'POPDFile ', c_switch%POPDFile
    print *, 'fixedYearPOPD ', c_switch%fixedYearPOPD
    print *, 'transientLGHT ', c_switch%transientLGHT
    print *, 'LGHTFile ', c_switch%LGHTFile
    print *, 'fixedYearLGHT ', c_switch%fixedYearLGHT
    print *, 'PFTCompetition ', c_switch%PFTCompetition
    print *, 'start_bare ', c_switch%start_bare
    print *, 'inibioclim ', c_switch%inibioclim
    print *, 'lnduseon ', c_switch%lnduseon
    print *, 'LUCFile ', c_switch%LUCFile
    print *, 'fixedYearLUC ', c_switch%fixedYearLUC
    print *, 'fertilizeron ', c_switch%fertilizeron
    print *, 'FERFile ', c_switch%FERFile
    print *, 'fixedYearFER ', c_switch%fixedYearFER
    print *, 'transientFER ', c_switch%transientFER
    print *, 'depositionon ', c_switch%depositionon
    print *, 'DEPFile ', c_switch%DEPFile
    print *, 'fixedYearDEP ', c_switch%fixedYearDEP
    print *, 'transientDEP ', c_switch%transientDEP
    print *, 'transientOBSWETF ', c_switch%transientOBSWETF
    print *, 'OBSWETFFile ', c_switch%OBSWETFFile
    print *, 'fixedYearOBSWETF ', c_switch%fixedYearOBSWETF
    print *, 'allLocalTime ', c_switch%allLocalTime
    print *, 'metFileFss ', c_switch%metFileFss
    print *, 'metFilefracFsf ', c_switch%metFilefracFsf
    print *, 'metFileFdl ', c_switch%metFileFdl
    print *, 'metFileSnow ', c_switch%metFileSnow
    print *, 'metFilePre ', c_switch%metFilePre
    print *, 'metFileTa ', c_switch%metFileTa
    print *, 'metFileQa ', c_switch%metFileQa
    print *, 'metFileUv ', c_switch%metFileUv
    print *, 'metFilePres ', c_switch%metFilePres
    print *, 'init_file ', c_switch%init_file
    print *, 'rs_file_to_overwrite ', c_switch%rs_file_to_overwrite
    print *, 'runparams_file ', c_switch%runparams_file
    print *, 'Comment ', c_switch%Comment
    print *, 'output_directory ', c_switch%output_directory
    print *, 'xmlFile ', c_switch%xmlFile
    print *, 'doperpftoutput ', c_switch%doperpftoutput
    print *, 'dopertileoutput ', c_switch%dopertileoutput
    print *, 'doChecksums ', c_switch%doChecksums
    print *, 'doAnnualOutput ', c_switch%doAnnualOutput
    print *, 'doMonthOutput ', c_switch%doMonthOutput
    print *, 'jmosty ', c_switch%jmosty
    print *, 'doDayOutput ', c_switch%doDayOutput
    print *, 'jdstd ', c_switch%jdstd
    print *, 'jdendd ', c_switch%jdendd
    print *, 'jdsty ', c_switch%jdsty
    print *, 'jdendy ', c_switch%jdendy
    print *, 'doHhOutput ', c_switch%doHhOutput
    print *, 'jhhstd ', c_switch%jhhstd
    print *, 'jhhendd ', c_switch%jhhendd
    print *, 'jhhsty ', c_switch%jhhsty
    print *, 'jhhendy ', c_switch%jhhendy
    print *, 'idisp ', c_switch%idisp
    print *, 'izref ', c_switch%izref
    print *, 'islfd ', c_switch%islfd
    print *, 'ipcp ', c_switch%ipcp
    print *, 'iwf ', c_switch%iwf
    print *, 'ITC ', c_switch%ITC
    print *, 'ITCG ', c_switch%ITCG
    print *, 'ITG ', c_switch%ITG
    print *, 'IPAI ', c_switch%IPAI
    print *, 'IHGT ', c_switch%IHGT
    print *, 'IALC ', c_switch%IALC
    print *, 'IALS ', c_switch%IALS
    print *, 'IALG ', c_switch%IALG
    print *, 'fracSnowParam ', c_switch%fracSnowParam
    print *, 'snoAlbedoParam ', c_switch%snoAlbedoParam
    print *, 'alb4BandParamsFile ', c_switch%alb4BandParamsFile
    print *, 'KsatScalaron ', c_switch%KsatScalaron
    print *, 'specifiedSoilPermDepth', c_switch%specifiedSoilPermDepth
    print *, 'EcTResistancePart ', c_switch%EcTResistancePart

    print *, 'pftexist ', vrot%pftexist
    print *, 'lfstatus ', vrot%lfstatus
    print *, 'pandays ', vrot%pandays
    print *, 'gleafmas ', vrot%gleafmas
    print *, 'gleafmas_ns ', vrot%gleafmas_ns
    print *, 'gleafmas_s ', vrot%gleafmas_s
    print *, 'leafns2s ', vrot%leafns2s
    print *, 'stemns2s ', vrot%stemns2s
    print *, 'rootns2s ', vrot%rootns2s
    print *, 're_alloc_s2l ', vrot%re_alloc_s2l
    print *, 're_alloc_r2l ', vrot%re_alloc_r2l
    print *, 're_alloc_sr2l ', vrot%re_alloc_sr2l
    print *, 'bleafmas ', vrot%bleafmas
    print *, 'stemmass ', vrot%stemmass
    print *, 'stemmass_ns ', vrot%stemmass_ns
    print *, 'stemmass_s ', vrot%stemmass_s
    print *, 'rootmass ', vrot%rootmass
    print *, 'rootmass_ns ', vrot%rootmass_ns
    print *, 'rootmass_s ', vrot%rootmass_s
    print *, 'pstemmass ', vrot%pstemmass
    print *, 'pgleafmass ', vrot%pgleafmass
    print *, 'fcancmx ', vrot%fcancmx
    print *, 'ngleafmas ', vrot%ngleafmas
    print *, 'ngleafmas_ns ', vrot%ngleafmas_ns
    print *, 'ngleafmas_s ', vrot%ngleafmas_s
    print *, 'nbleafmas ', vrot%nbleafmas
    print *, 'nstemmass ', vrot%nstemmass
    print *, 'nstemmass_ns ', vrot%nstemmass_ns
    print *, 'nstemmass_s ', vrot%nstemmass_s
    print *, 'nrootmass ', vrot%nrootmass
    print *, 'nrootmass_ns ', vrot%nrootmass_ns
    print *, 'nrootmass_s ', vrot%nrootmass_s
    print *, 'ailcg ', vrot%ailcg
    print *, 'ailcgs ', vrot%ailcgs
    print *, 'fcancs ', vrot%fcancs
    print *, 'fcanc ', vrot%fcanc
    print *, 'co2i1cg ', vrot%co2i1cg
    print *, 'co2i1cs ', vrot%co2i1cs
    print *, 'co2i2cg ', vrot%co2i2cg
    print *, 'co2i2cs ', vrot%co2i2cs
    print *, 'ancsveg ', vrot%ancsveg
    print *, 'ancgveg ', vrot%ancgveg
    print *, 'rmlcsveg ', vrot%rmlcsveg
    print *, 'rmlcgveg ', vrot%rmlcgveg
    print *, 'slai ', vrot%slai
    print *, 'ailcb ', vrot%ailcb
    print *, 'flhrloss ', vrot%flhrloss
    print *, 'flhrloss_ns ', vrot%flhrloss_ns
    print *, 'flhrloss_s ', vrot%flhrloss_s
    print *, 'grwtheff ', vrot%grwtheff
    print *, 'lystmmas ', vrot%lystmmas
    print *, 'lyrotmas ', vrot%lyrotmas
    print *, 'tymaxlai ', vrot%tymaxlai
    print *, 'stmhrlos ', vrot%stmhrlos
    print *, 'vgbiomas_veg ', vrot%vgbiomas_veg
    print *, 'emit_co2 ', vrot%emit_co2
    print *, 'emit_co ', vrot%emit_co
    print *, 'emit_ch4 ', vrot%emit_ch4
    print *, 'emit_nmhc ', vrot%emit_nmhc
    print *, 'emit_h2 ', vrot%emit_h2
    print *, 'emit_nox ', vrot%emit_nox
    print *, 'emit_n2o ', vrot%emit_n2o
    print *, 'emit_pm25 ', vrot%emit_pm25
    print *, 'emit_tpm ', vrot%emit_tpm
    print *, 'emit_tc ', vrot%emit_tc
    print *, 'emit_bc ', vrot%emit_bc
    print *, 'emit_oc ', vrot%emit_oc
    print *, 'emit_so2 ', vrot%emit_so2
    print *, 'elcg ', vrot%elcgrow
    print *, 'burnvegf ', vrot%burnvegf
    print *, 'smfuncveg ', vrot%smfuncveg
    print *, 'bterm ', vrot%bterm
    print *, 'mterm ', vrot%mterm
    print *, 'fFireCveg', vrot%fFireCveg
    print *, 'fFireLitter', vrot%fFireLitter
    print *, 'fFireCsoil', vrot%fFireCsoil  
    print *, 'bmasveg ', vrot%bmasveg
    print *, 'veghght ', vrot%veghght
    print *, 'rootdpth ', vrot%rootdpth
    print *, 'tltrleaf ', vrot%tltrleaf
    print *, 'tltrstem ', vrot%tltrstem
    print *, 'tltrroot ', vrot%tltrroot
    print *, 'leaflitr ', vrot%leaflitr
    print *, 'roottemp ', vrot%roottemp
    print *, 'afrleaf ', vrot%afrleaf
    print *, 'afrstem ', vrot%afrstem
    print *, 'afrroot ', vrot%afrroot
    print *, 'wtstatus ', vrot%wtstatus
    print *, 'ltstatus ', vrot%ltstatus
    print *, 'gppveg ', vrot%gppveg
    print *, 'vcmax0 ', vrot%vcmax0
    print *, 'nppveg ', vrot%nppveg
    print *, 'autoresveg ', vrot%autoresveg
    print *, 'rmlvegacc ', vrot%rmlvegacc
    print *, 'rmsveg ', vrot%rmsveg
    print *, 'rmrveg ', vrot%rmrveg
    print *, 'rgveg ', vrot%rgveg
    print *, 'litrfallveg ', vrot%litrfallveg
    print *, 'rothrlos ', vrot%rothrlos
    print *, 'pfcancmx ', vrot%pfcancmx
    print *, 'nfcancmx ', vrot%nfcancmx
    print *, 'anveg ', vrot%anveg
    print *, 'rmlveg ', vrot%rmlveg
    print *, 'bnf_free ', vrot%bnf_free
    print *, 'bnf_ant ', vrot%bnf_ant
    print *, 'bnf_nat ', vrot%bnf_nat
    print *, 'bnf_tot ', vrot%bnf_tot
    print *, 'bnf_moss ', vrot%bnf_moss
    print *, 'nstress ', vrot%nstress
    print *, 'nitrifveg ', vrot%nitrifveg
    print *, 'no_nitveg ', vrot%no_nitveg
    print *, 'no_denitveg ', vrot%no_denitveg
    print *, 'no_nitdenitveg ', vrot%no_nitdenitveg
    print *, 'n2o_nitveg ', vrot%n2o_nitveg
    print *, 'n2o_denitveg ', vrot%n2o_denitveg
    print *, 'n2o_nitdenitveg ', vrot%n2o_nitdenitveg
    print *, 'n2_denitveg ', vrot%n2_denitveg
    print *, 'nvolveg ', vrot%nvolveg
    print *, 'nleachveg ', vrot%nleachveg
    print *, 'appl_fert ', vrot%appl_fert
    print *, 'ndep_nh4 ', vrot%ndep_nh4
    print *, 'ndep_no3 ', vrot%ndep_no3
    print *, 'ndemandveg_wp_npp ', vrot%ndemandveg_wp_npp
    print *, 'nuptakeveg_p_nh4 ', vrot%nuptakeveg_p_nh4
    print *, 'nuptakeveg_p_no3 ', vrot%nuptakeveg_p_no3
    print *, 'nuptakeveg_a_actl_nh4 ', vrot%nuptakeveg_a_actl_nh4
    print *, 'nuptakeveg_a_actl_no3 ', vrot%nuptakeveg_a_actl_no3
    print *, 'nuptakeveg ', vrot%nuptakeveg
    print *, 'nleafns2sveg ', vrot%nleafns2sveg
    print *, 'nstemns2sveg ', vrot%nstemns2sveg
    print *, 'nrootns2sveg ', vrot%nrootns2sveg
    print *, 'nallocveg_l ', vrot%nallocveg_l
    print *, 'nallocveg_s ', vrot%nallocveg_s
    print *, 'nallocveg_r ', vrot%nallocveg_r
    print *, 'nresorpedveg_s ', vrot%nresorpedveg_s
    print *, 'nresorpedveg_r ', vrot%nresorpedveg_r
    print *, 'nre_allocveg_s2l ', vrot%nre_allocveg_s2l
    print *, 'nre_allocveg_r2l ', vrot%nre_allocveg_r2l
    print *, 'nlitrveg_l ', vrot%nlitrveg_l
    print *, 'nlitrveg_s ', vrot%nlitrveg_s
    print *, 'nlitrveg_r ', vrot%nlitrveg_r
    print *, 'nlitrveg ', vrot%nlitrveg
    print *, 'gl2bl_grass_nflux ', vrot%gl2bl_grass_nflux
    print *, 'c2nveg_l ', vrot%c2nveg_l
    print *, 'c2nveg_s ', vrot%c2nveg_s
    print *, 'c2nveg_r ', vrot%c2nveg_r
    print *, 'c2nveg_wp ', vrot%c2nveg_wp
    print *, 'c2nveg_litr ', vrot%c2nveg_litr
    print *, 'c2nveg_humus ', vrot%c2nveg_humus
    print *, 'nhumtrsveg ', vrot%nhumtrsveg
    print *, 'nmineralveg_litr ', vrot%nmineralveg_litr
    print *, 'nmineralveg_humus ', vrot%nmineralveg_humus
    print *, 'netnmineralveg ', vrot%netnmineralveg
    print *, 'nimmobilveg_nh4 ', vrot%nimmobilveg_nh4
    print *, 'nimmobilveg_no3 ', vrot%nimmobilveg_no3
    print *, 'nvgbiomas_veg ', vrot%nvgbiomas_veg
    print *, 'fNnetlandveg ', vrot%fNnetlandveg
    print *, 'redcoeff_vcmax ', vrot%redcoeff_vcmax
    print *, 'redcoeff_vcmaxMoss ', vrot%redcoeff_vcmaxMoss    
    print *, 'gavglai ', vrot%gavglai
    print *, 'co2conc ', vrot%co2conc
    print *, 'ch4conc ', vrot%ch4conc
    print *, 'canres ', vrot%canres
    print *, 'vgbiomas ', vrot%vgbiomas
    print *, 'gavgltms ', vrot%gavgltms
    print *, 'gavgscms ', vrot%gavgscms
    print *, 'burnfrac ', vrot%burnfrac
    print *, 'timharvarearot ', vrot%timharvarearow
    print *, 'tileAgerow ', vrot%tileAgerow
    print *, 'popdin ', vrot%popdin
    print *, 'soilpH ', vrot%soilpH
    print *, 'nfertil ', vrot%nfertil
    print *, 'ndeposit ', vrot%ndeposit
    print *, 'lterm ', vrot%lterm
    print *, 'extnprob ', vrot%extnprob
    print *, 'prbfrhuc ', vrot%prbfrhuc
    print *, 'rml ', vrot%rml
    print *, 'rms ', vrot%rms
    print *, 'rmr ', vrot%rmr
    print *, 'ch4WetSpec ', vrot%ch4WetSpec
    print *, 'wetfdyn ', vrot%wetfdyn
    print *, 'wetfrac_pres ', vrot%wetfrac_pres
    print *, 'ch4WetDyn ', vrot%ch4WetDyn
    print *, 'ch4_soills ', vrot%ch4_soills
    print *, 'lucemcom ', vrot%lucemcom
    print *, 'lucltrin ', vrot%lucltrin
    print *, 'lucsocin ', vrot%lucsocin
    print *, 'cproduct ', vrot%cproduct
    print *, 'fproductdecomp ', vrot%fproductdecomp
    print *, 'cleaf ', vrot%cleaf
    print *, 'cstem ', vrot%cstem
    print *, 'croot ', vrot%croot
    print *, 'nppleaf ', vrot%nppleaf
    print *, 'nppstem ', vrot%nppstem
    print *, 'npproot ', vrot%npproot
    print *, 'npp ', vrot%npp
    print *, 'nep ', vrot%nep
    print *, 'nbp ', vrot%nbp
    print *, 'gpp ', vrot%gpp
    print *, 'hetrores ', vrot%hetrores
    print *, 'autores ', vrot%autores
    print *, 'soilcresp ', vrot%soilcresp
    print *, 'rm ', vrot%rm
    print *, 'rg ', vrot%rg
    print *, 'litres ', vrot%litres
    print *, 'socres ', vrot%socres
    print *, 'dstcemls ', vrot%dstcemls
    print *, 'litrfall ', vrot%litrfall
    print *, 'humiftrs ', vrot%humiftrs
    print *, 'cfluxcg ', vrot%cfluxcg
    print *, 'cfluxcs ', vrot%cfluxcs
    print *, 'ROFB ', vrot%ROFB
    print *, 'dstcemls3 ', vrot%dstcemls3
    print *, 'uvaccrow_m ', vrot%uvaccrow_m
    print *, 'vvaccrow_m ', vrot%vvaccrow_m
    print *, 'qevpacc_m_save ', vrot%qevpacc_m_save
    print *, 'twarmm ', vrot%twarmm
    print *, 'tcoldm ', vrot%tcoldm
    print *, 'gdd5 ', vrot%gdd5
    print *, 'aridity ', vrot%aridity
    print *, 'srplsmon ', vrot%srplsmon
    print *, 'defctmon ', vrot%defctmon
    print *, 'anndefct ', vrot%anndefct
    print *, 'annsrpls ', vrot%annsrpls
    print *, 'annpcp ', vrot%annpcp
    print *, 'dry_season_length ', vrot%dry_season_length
    print *, 'peatlandType ', vrot%peatlandType
    print *, 'litrmsmoss ', vrot%litrmsmoss
    print *, 'upMossSoilC ', vrot%upMossSoilC
    print *, 'Cmossmas ', vrot%Cmossmas
    print *, 'litrmsmossN ', vrot%litrmsmossN
    print *, 'upMossSoilN ', vrot%upMossSoilN
    print *, 'Nmossmas ', vrot%Nmossmas
    print *, 'dmoss ', vrot%dmoss
    print *, 'peatSoilC ', vrot%peatSoilC
    print *, 'nppmoss ', vrot%nppmoss
    print *, 'armoss ', vrot%armoss
    print *, 'peatdep ', vrot%peatdep
    print *, 'colddays_leaffall ', vrot%colddays_leaffall
    print *, 'colddays_harvest ', vrot%colddays_harvest
    print *, 'zolnc ', vrot%zolnc
    print *, 'ailc ', vrot%ailc
    print *, 'cmasvegc ', vrot%cmasvegc
    print *, 'alvsctm ', vrot%alvsctm
    print *, 'paic ', vrot%paic
    print *, 'slaic ', vrot%slaic
    print *, 'alirctm ', vrot%alirctm
    print *, 'rmatc ', vrot%rmatc
    print *, 'rmatctem ', vrot%rmatctem
    print *, 'nepveg ', vrot%nepveg
    print *, 'nbpveg ', vrot%nbpveg
    print *, 'hetroresveg ', vrot%hetroresveg
    print *, 'litrmass ', vrot%litrmass
    print *, 'soilcmas ', vrot%soilcmas
    print *, 'litresveg ', vrot%litresveg
    print *, 'soilcresveg ', vrot%soilcresveg
    print *, 'humiftrsveg ', vrot%humiftrsveg
    print *, 'nh4_mass ', vrot%nh4_mass
    print *, 'no3_mass ', vrot%no3_mass
    print *, 'nlitrmass ', vrot%nlitrmass
    print *, 'soilnmas ', vrot%soilnmas
    print *, 'slopefrac ', vrot%slopefrac
    print *, 'dayl_max ', vrot%dayl_max
    print *, 'dayl ', vrot%dayl
    print *, 'grclarearow ', vrot%grclarearow
    print *, 'lygleafmasmax ', vrot%lygleafmasmax
    print *, 'lystemmassmax ', vrot%lystemmassmax
    print *, 'lyrootmassmax ', vrot%lyrootmassmax
    print *, 'lmaxt ', vrot%lmaxt
    print *, 'smaxt ', vrot%smaxt
    print *, 'rmaxt ', vrot%rmaxt

    print *, 'gleafmas ', vgat%gleafmas
    print *, 'gleafmas_ns ', vgat%gleafmas_ns
    print *, 'gleafmas_s ', vgat%gleafmas_s
    print *, 'leafns2s ', vgat%leafns2s
    print *, 'stemns2s ', vgat%stemns2s
    print *, 'rootns2s ', vgat%rootns2s
    print *, 're_alloc_s2l ', vgat%re_alloc_s2l
    print *, 're_alloc_r2l ', vgat%re_alloc_r2l
    print *, 're_alloc_sr2l ', vgat%re_alloc_sr2l
    print *, 'bleafmas ', vgat%bleafmas
    print *, 'stemmass ', vgat%stemmass
    print *, 'stemmass_ns ', vgat%stemmass_ns
    print *, 'stemmass_s ', vgat%stemmass_s
    print *, 'rootmass ', vgat%rootmass
    print *, 'rootmass_ns ', vgat%rootmass_ns
    print *, 'rootmass_s ', vgat%rootmass_s
    print *, 'pstemmass ', vgat%pstemmass
    print *, 'pgleafmass ', vgat%pgleafmass
    print *, 'fcancmx ', vgat%fcancmx
    print *, 'ngleafmas ', vgat%ngleafmas
    print *, 'ngleafmas_ns ', vgat%ngleafmas_ns
    print *, 'ngleafmas_s ', vgat%ngleafmas_s
    print *, 'nbleafmas ', vgat%nbleafmas
    print *, 'nstemmass ', vgat%nstemmass
    print *, 'nstemmass_ns ', vgat%nstemmass_ns
    print *, 'nstemmass_s ', vgat%nstemmass_s
    print *, 'nrootmass ', vgat%nrootmass
    print *, 'nrootmass_ns ', vgat%nrootmass_ns
    print *, 'nrootmass_s ', vgat%nrootmass_s
    print *, 'nvgbiomas_veg ', vgat%nvgbiomas_veg
    print *, 'gavglai ', vgat%gavglai
    print *, 'lightng ', vgat%lightng
    print *, 'zolnc ', vgat%zolnc
    print *, 'ailc ', vgat%ailc
    print *, 'ailcg ', vgat%ailcg
    print *, 'ailcgs ', vgat%ailcgs
    print *, 'fcancs ', vgat%fcancs
    print *, 'fcanc ', vgat%fcanc
    print *, 'co2conc ', vgat%co2conc
    print *, 'ch4conc ', vgat%ch4conc
    print *, 'co2i1cg ', vgat%co2i1cg
    print *, 'co2i1cs ', vgat%co2i1cs
    print *, 'co2i2cg ', vgat%co2i2cg
    print *, 'co2i2cs ', vgat%co2i2cs
    print *, 'ancsveg ', vgat%ancsveg
    print *, 'ancgveg ', vgat%ancgveg
    print *, 'rmlcsveg ', vgat%rmlcsveg
    print *, 'rmlcgveg ', vgat%rmlcgveg
    print *, 'slai ', vgat%slai
    print *, 'ailcb ', vgat%ailcb
    print *, 'canres ', vgat%canres
    print *, 'flhrloss ', vgat%flhrloss
    print *, 'flhrloss_ns ', vgat%flhrloss_ns
    print *, 'flhrloss_s ', vgat%flhrloss_s
    print *, 'grwtheff ', vgat%grwtheff
    print *, 'lystmmas ', vgat%lystmmas
    print *, 'lyrotmas ', vgat%lyrotmas
    print *, 'tymaxlai ', vgat%tymaxlai
    print *, 'vgbiomas ', vgat%vgbiomas
    print *, 'gavgltms ', vgat%gavgltms
    print *, 'gavgscms ', vgat%gavgscms
    print *, 'stmhrlos ', vgat%stmhrlos
    print *, 'rmatc ', vgat%rmatc
    print *, 'rmatctem ', vgat%rmatctem
    print *, 'litrmass ', vgat%litrmass
    print *, 'soilcmas ', vgat%soilcmas
    print *, 'vgbiomas_veg ', vgat%vgbiomas_veg
    print *, 'nh4_mass ', vgat%nh4_mass
    print *, 'no3_mass ', vgat%no3_mass
    print *, 'nlitrmass ', vgat%nlitrmass
    print *, 'soilnmas ', vgat%soilnmas
    print *, 'emit_co2 ', vgat%emit_co2
    print *, 'emit_co ', vgat%emit_co
    print *, 'emit_ch4 ', vgat%emit_ch4
    print *, 'emit_nmhc ', vgat%emit_nmhc
    print *, 'emit_h2 ', vgat%emit_h2
    print *, 'emit_nox ', vgat%emit_nox
    print *, 'emit_n2o ', vgat%emit_n2o
    print *, 'emit_pm25 ', vgat%emit_pm25
    print *, 'emit_tpm ', vgat%emit_tpm
    print *, 'emit_tc ', vgat%emit_tc
    print *, 'emit_bc ', vgat%emit_bc
    print *, 'emit_oc ', vgat%emit_oc
    print *, 'emit_so2 ', vgat%emit_so2
    print *, 'fien ', vgat%fiengat
    print *, 'fab ', vgat%fabgat
    print *, 'elcg ', vgat%elcggat
    print *, 'nmfi ', vgat%nmfigat
    print *, 'burnfrac ', vgat%burnfrac
    print *, 'timharvareagat ', vgat%timharvareagat
    print *, 'tileAgegat ', vgat%tileAgegat
    print *, 'burnvegf ', vgat%burnvegf
    print *, 'smfuncveg ', vgat%smfuncveg
    print *, 'popdin ', vgat%popdin
    print *, 'soilpH ', vgat%soilpH
    print *, 'nfertil ', vgat%nfertil
    print *, 'ndeposit ', vgat%ndeposit
    print *, 'bterm ', vgat%bterm
    print *, 'lterm ', vgat%lterm
    print *, 'mterm ', vgat%mterm
    print *, 'glcaemls ', vgat%glcaemls
    print *, 'blcaemls ', vgat%blcaemls
    print *, 'rtcaemls ', vgat%rtcaemls
    print *, 'stcaemls ', vgat%stcaemls
    print *, 'ltrcemls ', vgat%ltrcemls
    print *, 'ntchlveg ', vgat%ntchlveg
    print *, 'ntchsveg ', vgat%ntchsveg
    print *, 'extnprob ', vgat%extnprob
    print *, 'prbfrhuc ', vgat%prbfrhuc
    print *, 'dayl_max ', vgat%dayl_max
    print *, 'dayl ', vgat%dayl
    print *, 'bmasveg ', vgat%bmasveg
    print *, 'cmasvegc ', vgat%cmasvegc
    print *, 'veghght ', vgat%veghght
    print *, 'rootdpth ', vgat%rootdpth
    print *, 'rml ', vgat%rml
    print *, 'rms ', vgat%rms
    print *, 'tltrleaf ', vgat%tltrleaf
    print *, 'blfltrdt ', vgat%blfltrdt
    print *, 'glfltrdt ', vgat%glfltrdt
    print *, 'tltrstem ', vgat%tltrstem
    print *, 'tltrroot ', vgat%tltrroot
    print *, 'leaflitr ', vgat%leaflitr
    print *, 'roottemp ', vgat%roottemp
    print *, 'afrleaf ', vgat%afrleaf
    print *, 'afrstem ', vgat%afrstem
    print *, 'afrroot ', vgat%afrroot
    print *, 'wtstatus ', vgat%wtstatus
    print *, 'ltstatus ', vgat%ltstatus
    print *, 'rmr ', vgat%rmr
    print *, 'slopefrac ', vgat%slopefrac
    print *, 'wetfrac_pres ', vgat%wetfrac_pres
    print *, 'ch4WetSpec ', vgat%ch4WetSpec
    print *, 'wetfdyn ', vgat%wetfdyn
    print *, 'ch4WetDyn ', vgat%ch4WetDyn
    print *, 'ch4_soills ', vgat%ch4_soills
    print *, 'lucemcom ', vgat%lucemcom
    print *, 'lucltrin ', vgat%lucltrin
    print *, 'lucsocin ', vgat%lucsocin
    print *, 'cproduct ', vgat%cproduct
    print *, 'fproductdecomp ', vgat%fproductdecomp
    print *, 'cleaf ', vgat%cleaf
    print *, 'cstem ', vgat%cstem
    print *, 'croot ', vgat%croot
    print *, 'nppleaf ', vgat%nppleaf
    print *, 'nppstem ', vgat%nppstem
    print *, 'npproot ', vgat%npproot
    print *, 'npp ', vgat%npp
    print *, 'nep ', vgat%nep
    print *, 'nbp ', vgat%nbp
    print *, 'gpp ', vgat%gpp
    print *, 'hetrores ', vgat%hetrores
    print *, 'autores ', vgat%autores
    print *, 'soilcresp ', vgat%soilcresp
    print *, 'rm ', vgat%rm
    print *, 'rg ', vgat%rg
    print *, 'litres ', vgat%litres
    print *, 'socres ', vgat%socres
    print *, 'dstcemls ', vgat%dstcemls
    print *, 'litrfall ', vgat%litrfall
    print *, 'humiftrs ', vgat%humiftrs
    print *, 'gppveg ', vgat%gppveg
    print *, 'vcmax0 ', vgat%vcmax0
    print *, 'nepveg ', vgat%nepveg
    print *, 'bnf_free ', vgat%bnf_free
    print *, 'bnf_ant ', vgat%bnf_ant
    print *, 'bnf_nat ', vgat%bnf_nat
    print *, 'bnf_tot ', vgat%bnf_tot
    print *, 'bnf_moss ', vgat%bnf_moss
    print *, 'nstress ', vgat%nstress
    print *, 'nitrifveg ', vgat%nitrifveg
    print *, 'no_nitveg ', vgat%no_nitveg
    print *, 'no_denitveg ', vgat%no_denitveg
    print *, 'no_nitdenitveg ', vgat%no_nitdenitveg
    print *, 'n2o_nitveg ', vgat%n2o_nitveg
    print *, 'n2o_denitveg ', vgat%n2o_denitveg
    print *, 'n2o_nitdenitveg ', vgat%n2o_nitdenitveg
    print *, 'n2_denitveg ', vgat%n2_denitveg
    print *, 'nvolveg ', vgat%nvolveg
    print *, 'nleachveg ', vgat%nleachveg
    print *, 'appl_fert ', vgat%appl_fert
    print *, 'ndep_nh4 ', vgat%ndep_nh4
    print *, 'ndep_no3 ', vgat%ndep_no3
    print *, 'ndemandveg_wp_npp ', vgat%ndemandveg_wp_npp
    print *, 'nuptakeveg_p_nh4 ', vgat%nuptakeveg_p_nh4
    print *, 'nuptakeveg_p_no3 ', vgat%nuptakeveg_p_no3
    print *, 'nuptakeveg_a_actl_nh4 ', vgat%nuptakeveg_a_actl_nh4
    print *, 'nuptakeveg_a_actl_no3 ', vgat%nuptakeveg_a_actl_no3
    print *, 'nuptakeveg ', vgat%nuptakeveg
    print *, 'nleafns2sveg ', vgat%nleafns2sveg
    print *, 'nstemns2sveg ', vgat%nstemns2sveg
    print *, 'nrootns2sveg ', vgat%nrootns2sveg
    print *, 'nallocveg_l ', vgat%nallocveg_l
    print *, 'nallocveg_s ', vgat%nallocveg_s
    print *, 'nallocveg_r ', vgat%nallocveg_r
    print *, 'nresorpedveg_s ', vgat%nresorpedveg_s
    print *, 'nresorpedveg_r ', vgat%nresorpedveg_r
    print *, 'nre_allocveg_s2l ', vgat%nre_allocveg_s2l
    print *, 'nre_allocveg_r2l ', vgat%nre_allocveg_r2l
    print *, 'nlitrveg_l ', vgat%nlitrveg_l
    print *, 'nlitrveg_s ', vgat%nlitrveg_s
    print *, 'nlitrveg_r ', vgat%nlitrveg_r
    print *, 'nlitrveg ', vgat%nlitrveg
    print *, 'gl2bl_grass_nflux ', vgat%gl2bl_grass_nflux
    print *, 'c2nveg_l ', vgat%c2nveg_l
    print *, 'c2nveg_s ', vgat%c2nveg_s
    print *, 'c2nveg_r ', vgat%c2nveg_r
    print *, 'c2nveg_wp ', vgat%c2nveg_wp
    print *, 'c2nveg_litr ', vgat%c2nveg_litr
    print *, 'c2nveg_humus ', vgat%c2nveg_humus
    print *, 'nhumtrsveg ', vgat%nhumtrsveg
    print *, 'nmineralveg_litr ', vgat%nmineralveg_litr
    print *, 'nmineralveg_humus ', vgat%nmineralveg_humus
    print *, 'netnmineralveg ', vgat%netnmineralveg
    print *, 'nimmobilveg_nh4 ', vgat%nimmobilveg_nh4
    print *, 'nimmobilveg_no3 ', vgat%nimmobilveg_no3
    print *, 'fNnetlandveg ', vgat%fNnetlandveg
    print *, 'redcoeff_vcmax ', vgat%redcoeff_vcmax
    print *, 'redcoeff_vcmaxMoss ', vgat%redcoeff_vcmaxMoss
    print *, 'peatlandType ', vgat%peatlandType
    print *, 'peatdep ', vgat%peatdep
    print *, 'nppmoss ', vgat%nppmoss
    print *, 'armoss ', vgat%armoss
    print *, 'litrmsmoss ', vgat%litrmsmoss
    print *, 'upMossSoilC ', vgat%upMossSoilC
    print *, 'Cmossmas ', vgat%Cmossmas
    print *, 'litrmsmossN ', vgat%litrmsmossN
    print *, 'upMossSoilN ', vgat%upMossSoilN
    print *, 'Nmossmas ', vgat%Nmossmas    
    print *, 'dmoss ', vgat%dmoss
    print *, 'peatSoilC ', vgat%peatSoilC
    print *, 'ancsmoss ', vgat%ancsmoss
    print *, 'angsmoss ', vgat%angsmoss
    print *, 'ancmoss ', vgat%ancmoss
    print *, 'angmoss ', vgat%angmoss
    print *, 'rmlcsmoss ', vgat%rmlcsmoss
    print *, 'rmlgsmoss ', vgat%rmlgsmoss
    print *, 'rmlcmoss ', vgat%rmlcmoss
    print *, 'rmlgmoss ', vgat%rmlgmoss
    print *, 'nbpveg ', vgat%nbpveg
    print *, 'nppveg ', vgat%nppveg
    print *, 'hetroresveg ', vgat%hetroresveg
    print *, 'autoresveg ', vgat%autoresveg
    print *, 'litresveg ', vgat%litresveg
    print *, 'soilcresveg ', vgat%soilcresveg
    print *, 'humiftrsveg ', vgat%humiftrsveg
    print *, 'rmlvegacc ', vgat%rmlvegacc
    print *, 'rmsveg ', vgat%rmsveg
    print *, 'rmrveg ', vgat%rmrveg
    print *, 'rgveg ', vgat%rgveg
    print *, 'litrfallveg ', vgat%litrfallveg
    print *, 'reprocost ', vgat%reprocost
    print *, 'rothrlos ', vgat%rothrlos
    print *, 'pfcancmx ', vgat%pfcancmx
    print *, 'nfcancmx ', vgat%nfcancmx
    print *, 'alvsctm ', vgat%alvsctm
    print *, 'paic ', vgat%paic
    print *, 'slaic ', vgat%slaic
    print *, 'alirctm ', vgat%alirctm
    print *, 'cfluxcg ', vgat%cfluxcg
    print *, 'cfluxcs ', vgat%cfluxcs
    print *, 'CFLUX_GA ', vgat%CFLUX_GA
    print *, 'USTARBS_GA ', vgat%USTARBS_GA
    print *, 'ROFB ', vgat%ROFB
    print *, 'dstcemls3 ', vgat%dstcemls3
    print *, 'anveg ', vgat%anveg
    print *, 'rmlveg ', vgat%rmlveg
    print *, 'twarmm ', vgat%twarmm
    print *, 'tcoldm ', vgat%tcoldm
    print *, 'gdd5 ', vgat%gdd5
    print *, 'aridity ', vgat%aridity
    print *, 'srplsmon ', vgat%srplsmon
    print *, 'defctmon ', vgat%defctmon
    print *, 'anndefct ', vgat%anndefct
    print *, 'annsrpls ', vgat%annsrpls
    print *, 'annpcp ', vgat%annpcp
    print *, 'dry_season_length ', vgat%dry_season_length
    print *, 'colddays_leaffall ', vgat%colddays_leaffall
    print *, 'colddays_harvest ', vgat%colddays_harvest
    print *, 'tcurm ', vgat%tcurm
    print *, 'srpcuryr ', vgat%srpcuryr
    print *, 'dftcuryr ', vgat%dftcuryr
    print *, 'tmonth ', vgat%tmonth
    print *, 'anpcpcur ', vgat%anpcpcur
    print *, 'anpecur ', vgat%anpecur
    print *, 'gdd5cur ', vgat%gdd5cur
    print *, 'surmncur ', vgat%surmncur
    print *, 'defmncur ', vgat%defmncur
    print *, 'srplscur ', vgat%srplscur
    print *, 'defctcur ', vgat%defctcur
    print *, 'geremort ', vgat%geremort
    print *, 'intrmort ', vgat%intrmort
    print *, 'cc ', vgat%cc
    print *, 'mm ', vgat%mm
    print *, 'pftexist ', vgat%pftexist
    print *, 'lfstatus ', vgat%lfstatus
    print *, 'pandays ', vgat%pandays
    print *, 'grclarea ', vgat%grclarea
    print *, 'altotcount_ctm ', vgat%altotcount_ctm
    print *, 'todfrac ', vgat%todfrac
    print *, 'fsinacc_gat ', vgat%fsinacc_gat
    print *, 'flutacc_gat ', vgat%flutacc_gat
    print *, 'flinacc_gat ', vgat%flinacc_gat
    print *, 'altotacc_gat ', vgat%altotacc_gat
    print *, 'netrad_gat ', vgat%netrad_gat
    print *, 'preacc_gat ', vgat%preacc_gat
    print *, 'sdepgat ', vgat%sdepgat
    print *, 'sandgat ', vgat%sandgat
    print *, 'xdiffusgat ', vgat%xdiffusgat
    print *, 'faregat ', vgat%faregat
    print *, 'lygleafmasmax ', vgat%lygleafmasmax
    print *, 'lystemmassmax ', vgat%lystemmassmax
    print *, 'lyrootmassmax ', vgat%lyrootmassmax
    print *, 'lmaxt ', vgat%lmaxt
    print *, 'smaxt ', vgat%smaxt
    print *, 'rmaxt ', vgat%rmaxt

    print *, 'tracer mossCMassrot ', tracer%mossCMassrot
    print *, 'tracer mossLitrMassrot ', tracer%mossLitrMassrot
    print *, 'tracer tracerCO2rot ', tracer%tracerCO2rot
    print *, 'tracer gLeafMassrot ', tracer%gLeafMassrot
    print *, 'tracer bLeafMassrot ', tracer%bLeafMassrot
    print *, 'tracer stemMassrot ', tracer%stemMassrot
    print *, 'tracer rootMassrot ', tracer%rootMassrot
    print *, 'tracer litrMassrot ', tracer%litrMassrot
    print *, 'tracer soilCMassrot ', tracer%soilCMassrot
    print *, 'tracer rothrlosrot ', tracer%rothrlosrot
    print *, 'tracer stmhrlosrot ', tracer%stmhrlosrot
    print *, 'tracer mossCMassgat ', tracer%mossCMassgat
    print *, 'tracer mossLitrMassgat ', tracer%mossLitrMassgat
    print *, 'tracerCO2gat ', tracer%tracerCO2gat
    print *, 'tracer gLeafMassgat ', tracer%gLeafMassgat
    print *, 'tracer bLeafMassgat ', tracer%bLeafMassgat
    print *, 'tracer stemMassgat ', tracer%stemMassgat
    print *, 'tracer rootMassgat ', tracer%rootMassgat
    print *, 'tracer litrMassgat ', tracer%litrMassgat
    print *, 'tracer soilCMassgat ', tracer%soilCMassgat
    print *, 'tracer rothrlosgat ', tracer%rothrlosgat
    print *, 'tracer stmhrlosgat ', tracer%stmhrlosgat

    print *, 'taaccgat_t ', ctem_tile%taaccgat_t
    print *, 'uvaccgat_t ', ctem_tile%uvaccgat_t
    print *, 'vvaccgat_t ', ctem_tile%vvaccgat_t
    print *, 'anmossac_t ', ctem_tile%anmossac_t
    print *, 'rmlmossac_t ', ctem_tile%rmlmossac_t
    print *, 'gppmossac_t ', ctem_tile%gppmossac_t
    print *, 'QFCacc_t ', ctem_tile%QFCacc_t
    print *, 'tbaraccgat_t ', ctem_tile%tbaraccgat_t
    print *, 'thliqacc_t ', ctem_tile%thliqacc_t
    print *, 'thiceacc_t ', ctem_tile%thiceacc_t
    print *, 'ancgvgac_t ', ctem_tile%ancgvgac_t
    print *, 'rmlcgvga_t ', ctem_tile%rmlcgvga_t

    print *, 'end ctem dump'

  end subroutine ctemdump
  
  !=================================================================================

  !> \namespace ctemstatevars
  !> Contains the biogeochemistry-related variable type structures.
  !! 1. vrot - CTEM's 'rot' vars
  !! 2. vgat - CTEM's 'gat' vars
  !! 3. ctem_tile - CTEM's variables per tile
  !
end module ctemStateVars
