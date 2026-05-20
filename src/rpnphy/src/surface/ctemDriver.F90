!> \file
!> Principle driver for CLASSIC's biogeochemistry (the CTEM sub-model)
!! @author V. Arora, J. Melton, Y. Peng, R. Shrestha, A. Asaadi
!!
module ctemDriver

  implicit none

  ! subroutines contained in this module:

  public  :: ctem
  public  :: calcNPP
  public  :: calcNEP
  public  :: calcNBP

contains
  !> \ingroup ctemdriver_ctem
  
  subroutine ctem (fsnow, sand, ilmos, & ! In
                   ilg, il1, il2, iday, radj, & ! In
                   ta, ancgveg, rmlcgveg, & ! In
                   zbotw, doMethane, & ! In
                   uwind, vwind, lightng, & ! In 
                   elcggat, tbar, & ! In
                   spinfast, todfrac, & ! In
                   netrad, precip, psisat, & ! In
                   grclarea, popdin, isand, & ! In
                   wetfrac, slopefrac, bi, & ! In
                   thpor, currlat, ch4conc, & ! In
                   THFC, THLW, thliq, thice, & ! In
                   peatlandType, mossPresent, anmoss, gppmoss, & ! In
                   wtable, maxAnnualActLyr, & ! In
                   PFTCompetition, dofire, fireSubmodule, lnduseon, inibioclim, & ! In
                   leapnow, useTracer, tracerCO2, useStaticPeatDep, & ! In
                   pfcancmx, nfcancmx, & ! In/Out
                   Ncycle_on, CFLUX_GA, soilpH, & ! In
                   USTARBS_GA, nfertil, ndeposit, & !In
                   ROFB, QFC, co2conc, dayl, dayl_max, zbot, & ! In
                   prsfireareagat, prescribedFire, & !In
                   delzw, soildpth, & ! In
                   stemmass, stemmass_ns, stemmass_s, & ! In/ Out
                   rootmass, rootmass_ns, rootmass_s, litrmass, & ! In/ Out
                   gleafmas, gleafmas_ns, gleafmas_s, & ! In/ Out
                   bleafmas, soilcmas, ailcg, ailc, & ! In/ Out
                   zolnc, rmatctem, rmatc, ailcb, & ! In/ Out
                   flhrloss, flhrloss_ns, flhrloss_s, tracerFlHrLoss, pandays, lfstatus, grwtheff, & ! In/ Out
                   lystmmas, lyrotmas, lmaxt, smaxt, rmaxt, & ! In/ Out
                   lygleafmasmax, lystemmassmax, lyrootmassmax, & ! In/ Out
                   tymaxlai, vgbiomas, & ! In/ Out
                   gavgltms, gavgscms, stmhrlos, tracerStmhrlos, slai, & ! In/ Out
                   bmasveg, cmasvegc, colddays_leaffall, & ! In/Out
                   colddays_harvest, rothrlos, tracerRothrlos, & ! In/ Out
                   fcanmx, alvisc, alnirc, gavglai, & ! In/ Out
                   Cmossmas, litrmsmoss, upMossSoilC, peatdep, peatSoilC, fcancmx, &! In/ Out
                   Nmossmas, litrmsmossN, upMossSoilN, rmlmoss, &! In/ Out
                   geremort, intrmort, pstemmass, pgleafmass, &! In/ Out
                   tcurm, srpcuryr, dftcuryr, &! In/ Out
                   tmonth, anpcpcur, anpecur, gdd5cur, &! In/ Out
                   surmncur, defmncur, srplscur, defctcur, &! In/ Out
                   aridity, srplsmon, defctmon, anndefct, &! In/ Out
                   annsrpls, annpcp, dry_season_length, &! In/ Out
                   pftexist, twarmm, tcoldm, gdd5, nppveg, &! In/ Out
                   tracerStemMass, tracerRootMass, tracerGLeafMass, tracerBLeafMass, & ! In/Out
                   tracerSoilCMass, tracerLitrMass, tracerMossCMass, tracerMossLitrMass, & ! In/Out
                   leafns2s, stemns2s, rootns2s, &! In/Out
                   re_alloc_s2l, re_alloc_r2l, re_alloc_sr2l, &! In/Out
                   bnf_nat, bnf_ant, nstress, bnf_moss, & ! In/Out
                   trackTileAge, dynamicTilingOn, tileAgegat, timberharvest, timharvareagat, timharvarearow, & ! In/Out
                   npp, nep, hetrores, autores, & ! Out (Primary)
                   soilresp, rm, rg, nbp, & ! Out (Primary)
                   litres, socres, gpp, dstcemls1, & ! Out (Primary)
                   litrfall, humiftrs, veghght, rootdpth, & ! Out (Primary)
                   rml, rms, rmr, tltrleaf, & ! Out (Primary)
                   tltrstem, tltrroot, leaflitr, roottemp, & ! Out (Primary)
                   burnfrac, lucemcom, lucltrin, & ! Out (Primary)
                   cproduct, fproductdecomp, nMineralNH4, nMineralNO3, nVeg, nLitter, nSoil, nLeaf, nStem, nRoot, & ! Out (Primary)
                   fBNF, fNdep, fNfert, fNgasNonFire, fNgasFire, fNgas, fNnetmin, fNOx, fNup, fNvegSoil, fN2o, &  ! Out (Primary)
                   cleaf, cstem, croot, nppleaf, nppstem, npproot, & ! Out (Secondary)
                   lucsocin, lucemcomn, lucltrinn, lucsocinn, dstcemls3, & ! Out (Primary)
                   ch4WetSpec, ch4WetDyn, wetfdyn, ch4soills, & ! Out (Primary)
                   paicgat, slaicgat, & ! Out (Primary)
                   emit_co2, emit_ch4, reprocost, blfltrdt, glfltrdt, &  ! Out (Primary)
                   glcaemls, blcaemls, rtcaemls, stcaemls, ltrcemls, &  ! Out (Primary)
                   fFireCveg, fFireLitter, fFireCsoil, & ! Out (Primary)
                   ntchlveg, ntchsveg, ntchrveg, &  ! Out (Primary)
                   emit_co, emit_nmhc, smfunc_veg, & ! Out (Secondary)
                   emit_h2, emit_nox, emit_n2o, emit_nh3, emit_pm25, & ! Out (Secondary)
                   emit_tpm, emit_tc, emit_bc, emit_oc, emit_so2, & ! Out (Secondary)
                   fiengat, fabgat, nmfigat, & ! Out (Secondary)
                   bterm_veg, lterm, mterm_veg, burnvegf, & ! Out (Secondary)
                   litrfallveg, humtrsvg, ltstatus, & ! Out (Secondary)
                   afrleaf, afrstem, afrroot, wtstatus, & ! Out (Secondary)
                   rmlveg, rmsveg, rmrveg, rgveg, & ! Out (Secondary)
                   vgbiomas_veg, & ! Out (Secondary)
                   gppveg, nepveg, nbpveg, & ! Out (Secondary)
                   hetrsveg, autoresveg, ltresveg, scresveg, & ! Out (Secondary)
                   nppmoss, armoss, & ! Out (Secondary)
                   colrate, mortrate, & ! Out (Secondary)
                   ngleafmas, ngleafmas_ns, ngleafmas_s, nbleafmas, & ! Out
                   nstemmass, nstemmass_ns, nstemmass_s, & ! Out
                   nrootmass, nrootmass_ns, nrootmass_s, & ! Out
                   nlitrmass, soilnmas, & ! Out
                   bnf_free, bnf_tot, & !Out
                   nh4_mass, no3_mass, nvolveg, nleachveg, fNleach, fNvol, & ! Out
                   nitrifveg, no_nitveg, no_denitveg, & ! Out
                   no_nitdenitveg, n2o_nitveg, n2o_denitveg, & ! Out
                   n2o_nitdenitveg, n2_denitveg, appl_fert, & ! Out
                   ndep_nh4, ndep_no3, ndemandveg_wp_npp, & ! Out
                   nuptakeveg_p_nh4, nuptakeveg_p_no3, & ! Out
                   nuptakeveg_a_actl_nh4, & ! Out
                   nuptakeveg_a_actl_no3, nuptakeveg, nallocveg_l, nallocveg_s, & ! Out
                   nallocveg_r, nresorpedveg_s, & ! Out
                   nresorpedveg_r, nre_allocveg_s2l, & ! Out
                   nre_allocveg_r2l, nleafns2sveg, nstemns2sveg, & ! Out
                   nrootns2sveg, nlitrveg_l, nlitrveg_s, & ! Out
                   nlitrveg_r, nlitrveg, gl2bl_grass_nflux, & ! Out
                   c2nveg_l, c2nveg_s, c2nveg_r, & ! Out
                   c2nveg_wp, c2nveg_litr, c2nveg_humus, & ! Out
                   nhumtrsveg, nmineralveg_litr, & ! Out
                   nmineralveg_humus, netnmineralveg, nimmobilveg_nh4, & ! Out
                   nimmobilveg_no3, nvgbiomas_veg, & ! Out
                   fNnetlandveg, redcoeff_vcmax, redcoeff_vcmaxMoss, faregat, nepCMIP, drgtstrs, betadrgt, & ! Out
                   QACGAT, PRESGAT, isi, bui, fwi, ffmc, dmc, dc, IMONTH, NoonTempc, &
                   NoonPrecip, NoonWindtot, NoonRhum, snwct) ! For FWI test 


    use classicParams,             only : kk, pi, zero, icp1, &
                                          iccp1, ican, nlat, &
                                          ignd, icc, nmos, l2max, grescoef, &
                                          humicfac, laimin, laimax, &
                                          crop, repro_fraction, &
                                          rmortmoss, humicfacmoss, GRAV, RHOW, RHOICE, &
                                          Fmax_leaf, min_ns2t_l, gleafmasmax, &
                                          Fmax_stem, min_ns2t_s, stemmassmax, &
                                          Fmax_root, min_ns2t_r, rootmassmax, &
                                          ccost_coeff, &
                                          classpfts, ctempfts, iccp2, humicfac_bg, &
                                          deltat, tolrance, convertg2kg, convertkgN, &
                                          convertn2o2N, convertnh32N, convertnox2N, c_switch
    use landuseChange,             only : luc
    use competitionScheme,         only : bioclim, existence, competition, expansion
    use disturbance_scheme,        only : disturbance
    use heterotrophicRespirationMod, only : heterotrophicRespiration, updatePoolsHetResp
    use peatlandsMod,              only : peatDayEnd, peatDepth, peatStorage
    use mossMod,                   only : updateMossC
    use ctemUtilities,             only : genSortIndex
    use autotrophicRespiration,    only : mainres
    use balanceCarbon,             only : balcar, prepBalanceC
    use mortality,                 only : mortalty, updatePoolsMortality
    use turnoverMod,               only : turnoverStemRoot, updatePoolsTurnover
    use applyAllometry,            only : allometry
    use tracerModule,              only : prepTracer, doTracerBalance, checkTracerBalance
    use methaneProcesses,          only : soil_ch4uptake, wetlandMethane
    use soilCProcesses,            only : turbation
    use allocateCarbon,            only : allocate, updatePoolsAllocateRepro
    use tiledDisturbance,          only : harvestTile
    use dynamicTiling,             only : indexTileAge

    implicit none

    ! inputs

    logical, intent(in) :: lnduseon                         !< logical switch to run the land use change subroutine or not.
    logical, intent(in) :: Ncycle_on                        !< logical switch to run the nitrogen cycle subroutines or not.
    logical, intent(in) :: PFTCompetition                   !< logical boolean telling if competition between pfts is on or not
    logical, intent(in) :: dofire                           !< boolean, if true allow fire, if false no fire.
    character(10), intent(in) :: fireSubmodule !< This selects the fire submodule. The options here are:
                                                            !! 'Default' - this is the present global fire model (e.g. Arora and Melton 2018 Nat. Commm.)
                                                            !! 'FWI' - Canada-specific fire weather index based fire model. Presently NOT set up for global
                                                            !! runs.
    logical, intent(in) :: leapnow                          !< true if this year is a leap year. Only used if the switch 'leap' is true.
    logical, intent(in) :: doMethane                        !< true if you wish to do the methane calculations.
    integer, intent(in) :: iday                             !< day of year
    integer, intent(in) :: spinfast                         !< spinup factor for soil carbon & nitrogen whose default value is 1. as this factor increases the
                                                            !< soil c & n pools will come into equilibrium faster. reasonable value for spinfast is
                                                            !< between 5 and 10. when spinfast/=1 then the balcar & balnitro subroutines are not run.
    integer, intent(in) :: ilg                              !< ilg=no. of grid cells in latitude circle
    integer, intent(in) :: il1                              !< il1=1
    integer, intent(in) :: il2                              !< il2=ilg (no. of grid cells in latitude circle)
    integer, intent(in), dimension(:) :: ilmos              !< Index of gridcell corresponding to current element of gathered vector of land surface variables [ ]
    integer, dimension(ilg,ignd), intent(in) :: isand       !<
    character(8), dimension(ilg), intent(in) :: peatlandType !< Peatland type (else 'None')
    character(8), dimension(ilg), intent(in) :: mossPresent !< Type of moss present (else 'None')
    real, dimension(ilg), intent(in) :: dayl_max            !< maximum daylength at each location
    real, dimension(ilg), intent(in) :: dayl                !< daylength at each location
    real, dimension(ilg), intent(in) :: fsnow               !< fraction of snow simulated by class
    real, dimension(ilg,ignd), intent(in) :: sand           !< percentage sand
    real, dimension(ilg), intent(in) :: radj                !< latitude in radians
    real, dimension(ilg,ignd), intent(in) ::  tbar          !< Soil temperature, K
    real, dimension(ilg,ignd), intent(in) :: psisat         !< Saturated soil matric potential (m)
    real, dimension(ilg,ignd), intent(in) :: bi             !< Brooks and Corey/Clapp and Hornberger b term
    real, dimension(ilg,ignd), intent(in) :: thpor          !< Soil total porosity \f$(cm^3 cm^{-3})\f$ - daily average
    real, dimension(ilg), intent(in) :: ta                  !< air temp, K
    real, dimension(ilg), intent(inout) :: soildpth            !< soil depth (m)
    real, dimension(ilg,ignd), intent(in) :: zbotw          !< bottom of soil layers
    real, dimension(ilg,ignd), intent(inout) :: delzw          !< thicknesses of the soil layers
    real, dimension(ignd), intent(in) :: zbot               !< Depth of to the bottom of soil layer [m]
    !
    real, dimension(ilg,ignd), intent(in) :: thliq          !< liquid mois. content of soil layers
    real, dimension(ilg,ignd), intent(in) :: thice          !< Frozen soil moisture content
    real, dimension(ilg), intent(in) ::  grclarea           !< area of the grid cell, \f$km^2\f$
    real, dimension(ilg), intent(in) ::  currlat            !< centre latitude of grid cells in degrees
    real, dimension(ilg), intent(in) :: uwind               !< u wind speed, m/s
    real, dimension(ilg), intent(in) :: vwind               !< v wind speed, m/s
    real, dimension(ilg), intent(in) ::  precip             !< daily precipitation (mm/day)
    real, dimension(ilg), intent(in) ::  netrad             !< daily net radiation (w/m2)
    real, dimension(ilg), intent(in) :: lightng             !< total lightning frequency, flashes/km2.year
    real, dimension(ilg), intent(in) :: elcggat             !< cloud-to-ground lightning \f$[flashes/cm2/sec]\f$

    real, dimension(ilg,icc), intent(in) :: todfrac         !< max. fractional coverage of ctem's 9 pfts by the end of the day, for use by land use subroutine
    real, dimension(ilg), intent(in) :: ch4conc             !< Atmospheric \f$CH_4\f$ concentration at the soil surface (ppmv)
    real, dimension(ilg), intent(in) :: wetfrac             !< Prescribed fraction of wetlands in a grid cell
    real, dimension(ilg,8), intent(in) :: slopefrac         !<
    real, dimension(ilg), intent(in) :: anmoss              !< moss net photoysnthesis -daily averaged C fluxes rates (umol/m2/s)
    real, dimension(ilg), intent(in) :: gppmoss             !< moss GPP -daily averaged C fluxes rates (umol/m2/s)
    real, dimension(ilg), intent(in) :: wtable              !< water table (m)
    real, dimension(ilg,icc), intent(in) :: ancgveg         !< net photosynthetic rate for CTEM's pfts
    real, dimension(ilg,icc), intent(in) :: rmlcgveg        !< leaf respiration rate for CTEM's pfts
    real, dimension(ilg,ignd), intent(in) :: THFC           !<
    real, dimension(ilg,ignd), intent(in) :: THLW           !<
    real, dimension(ilg), intent(in) :: maxAnnualActLyr     !< Active layer depth maximum over the e-folding period specified by parameter eftime (m).
    real, intent(in) :: tracerCO2(:)       !< Tracer CO2 value read in from tracerCO2File, units vary (simple: ppm, 14C \f$\Delta ^{14}C\f$)
    character(6), intent(in) :: useTracer !< character: Switch for use of a model tracer. If useTracer is 'None' then the
                                          !! tracer code is not used. useTracer = 'Simple' turns on a simple tracer that tracks
                                          !! pools and fluxes. The simple tracer then requires that the tracer values in
                                          !! the init_file and the tracerCO2file are set to meaningful values for the experiment being run.
                                          !! useTracer = '14C' means the tracer is 14C and will then call a 14C decay scheme.
                                          !! useTracer = '13C' means the tracer is 13C and will then call a 13C fractionation scheme.
    logical, intent(in) :: useStaticPeatDep  !< True will keep the peat depth at the sdep value. This is used for spinup to ensure
                                              !! that the disequilibrium between the peatland and the driving climate, time period 
                                              !! since the peatland initiation, or sub-grid heterogeneity impacts don't adversely impact the ability 
                                              !! of the model to spinup.

    real, dimension(ilg), intent(in) :: CFLUX_GA       !< aerodynamic conductance (\f$[m s^{-1} ]\f$), inverse of boundary layer aerodynamic resistance (ra)
    real, dimension(ilg), intent(in) :: USTARBS_GA     !< friction velocity to be used in nitrogen volatilization  \f$[m s^{-1} ]\f$
    real, dimension(ilg), intent(in) :: ROFB           !< base flow from bottom of soil column \f$[kg m^{-2} s^{-1} ]\f$
    real, dimension(ilg,ignd), intent(in) :: QFC       !< water removed from soil layers by transpiration \f$[kg m^{-2} s^{-1}]\f$
    real, dimension(ilg), intent(in) :: co2conc        !< ATMOS. CO2 CONC. IN PPM
    real, intent(in) :: faregat(ilg)
    real, intent(in) :: prsfireareagat(ilg)             !< FLAG rsc  please add info for these
    logical, intent(in) :: prescribedFire
    logical, intent(in) :: trackTileAge

    integer, dimension(ilg), intent(inout) :: snwct

    !     updates
    !
    logical, intent(inout) :: pftexist(ilg,icc)             !< True if PFT is present in tile.
    logical, intent(inout) :: inibioclim                    !< switch telling if bioclimatic parameters are being initialized from scratch (false)
    !< or being initialized from some spun up values(true).
    integer, dimension(ilg,icc), intent(inout) :: pandays   !< days with positive net photosynthesis (an) for use in the phenology subroutine
    integer, dimension(ilg), intent(inout) :: colddays_leaffall !< cold days counter for tracking days below a certain temperature threshold for ndl dcd 
    integer, dimension(ilg), intent(inout) :: colddays_harvest  !< cold days counter for tracking days below a certain temperature threshold for crops
    integer, dimension(ilg,icc), intent(inout) :: lfstatus  !< leaf phenology status
    real, dimension(ilg,icc), intent(inout) :: fcancmx      !< max. fractional coverage of CTEM's pfts, but this can be
    !< modified by land-use change,and competition between pfts
    real, dimension(ilg,icc), intent(inout) :: pfcancmx        !< previous year's fractional coverages of pfts
    real, dimension(ilg,icc), intent(inout) :: nfcancmx        !< next year's fractional coverages of pfts
    real, dimension(ilg,ican,ignd), intent(inout) :: rmatc  !< fraction of roots for each of class' 4 pfts in each soil layer
    real, dimension(ilg), intent(inout) :: surmncur         !< number of months with surplus water for current year
    real, dimension(ilg), intent(inout) :: defmncur         !< number of months with water deficit for current year
    real, dimension(ilg), intent(inout) :: tcurm            !< temperature of the current month (c)
    real, dimension(ilg), intent(inout) :: annpcp           !< annual precipitation (mm)
    real, dimension(ilg), intent(inout) :: dry_season_length !< length of the dry season (months)
    real, dimension(ilg), intent(inout) :: twarmm           !< temperature of the warmest month (c)
    real, dimension(ilg), intent(inout) :: tcoldm           !< temperature of the coldest month (c)
    real, dimension(ilg), intent(inout) :: gdd5             !< growing degree days above 5 c
    real, dimension(ilg), intent(inout) :: aridity          !< aridity index, ratio of potential evaporation to precipitation
    real, dimension(ilg), intent(inout) :: srplsmon         !< number of months in a year with surplus water i.e. precipitation more than potential evaporation
    real, dimension(ilg), intent(inout) :: defctmon         !< number of months in a year with water deficit i.e. precipitation less than potential evaporation
    real, dimension(12,ilg), intent(inout) :: tmonth        !< monthly temperatures
    real, dimension(ilg), intent(inout) :: anpcpcur         !< annual precipitation for current year (mm)
    real, dimension(ilg), intent(inout) :: anpecur          !< annual potential evaporation for current year (mm)
    real, dimension(ilg), intent(inout) :: gdd5cur          !< growing degree days above 5 c for current year
    real, dimension(ilg), intent(inout) :: srplscur         !< water surplus for the current month
    real, dimension(ilg), intent(inout) :: defctcur         !< water deficit for the current month
    real, dimension(ilg), intent(inout) :: srpcuryr         !< water surplus for the current year
    real, dimension(ilg), intent(inout) :: dftcuryr         !< water deficit for the current year
    real, dimension(ilg), intent(inout) :: anndefct         !< annual water deficit (mm)
    real, dimension(ilg), intent(inout) :: annsrpls         !< annual water surplus (mm)
    real, dimension(ilg,icc), intent(inout) :: stemmass     !< stem mass for each of the ctem pfts, \f$(kg C/m^2)\f$
    real, dimension(ilg,icc), intent(inout) :: stemmass_ns  !< non-structural stem mass for each of the ctem pfts, \f$(kg C/m^2)\f$
    real, dimension(ilg,icc), intent(inout) :: stemmass_s   !< structural stem mass for each the 9 ctem pfts, \f$(kg C/m^2)\f$
    real, dimension(ilg,icc), intent(inout) :: rootmass     !< root mass for each of the ctem pfts, \f$(kg C/m^2)\f$
    real, dimension(ilg,icc), intent(inout) :: rootmass_ns  !< non-structural root mass for each of the ctem pfts, \f$(kg C/m^2)\f$
    real, dimension(ilg,icc), intent(inout) :: rootmass_s   !< structural root mass for each of the ctem pfts, \f$(kg C/m^2)\f$
    real, dimension(ilg,iccp2,ignd), intent(inout) :: litrmass   !< litter mass for each of the ctem pfts + bare + LUC product pools, \f$(kg C/m^2)\f$
    real, dimension(ilg,icc), intent(inout) :: gleafmas     !< green leaf mass for each of the ctem pfts, \f$(kg C/m^2)\f$
    real, dimension(ilg,icc), intent(inout) :: gleafmas_ns  !< non-structural green leaf mass for each of the ctem pfts, \f$(kg C/m^2)\f$
    real, dimension(ilg,icc), intent(inout) :: gleafmas_s   !< structural green leaf mass for each of the ctem pfts, \f$(kg C/m^2)\f$
    real, dimension(ilg,icc) :: re_alloc_sr2l(ilg,icc)      !< reallocated nsc from stem and root to leaves during leaf out
    real, dimension(ilg,icc) :: re_alloc_s2l(ilg,icc)       !< reallocated nsc from stem to leaves during leaf out
    real, dimension(ilg,icc) :: re_alloc_r2l(ilg,icc)       !< reallocated nsc from root to leaves during leaf out
    real, dimension(ilg,icc), intent(inout) :: bleafmas     !< brown leaf mass for each of the ctem pfts, \f$(kg C/m^2)\f$
    real, dimension(ilg,iccp2,ignd), intent(inout) :: soilcmas   !< soil carbon mass for each of the ctem pfts + bare + LUC product pools, \f$(kg C/m^2)\f$

    real, dimension(ilg,icc), intent(inout) :: ailcg        !< Green LAI for ctem's pfts \f$(m^2 leaf/m^2 ground)\f$
    real, dimension(ilg,ican), intent(inout) :: ailc        !< lumped lai for class' 4 pfts
    real, dimension(ilg,icc,ignd), intent(inout) :: rmatctem !< fraction of roots for each of ctem's 9 pfts in each soil layer
    real, dimension(ilg,ican), intent(inout) :: zolnc       !< lumped log of roughness length for class' 4 pfts
    real, dimension(ilg,icc), intent(inout) :: ailcb        !< brown lai for ctem's 9 pfts. for now we assume only grasses can have brown lai
    real, dimension(ilg), intent(inout) :: vgbiomas         !< grid averaged vegetation biomass, \f$(kg C/m^2)\f$
    real, dimension(ilg), intent(inout) :: gavgltms         !< grid averaged litter mass, \f$(kg C/m^2)\f$
    real, dimension(ilg), intent(inout) :: gavgscms         !< grid averaged soil c mass, \f$(kg C/m^2)\f$
    real, dimension(ilg), intent(inout) :: gavglai          !< grid averaged green leaf area index
    real, dimension(ilg,icc), intent(inout) :: bmasveg      !< total (gleaf + stem + root) biomass for each ctem pft, \f$(kg C/m^2)\f$
    real, dimension(ilg,ican), intent(inout) :: cmasvegc    !< total canopy mass for each of the 4 class pfts. recall that class requires canopy
    !< mass as an input,and this is now provided by ctem. \f$kg/m^2\f$.
    real, dimension(ilg,icp1), intent(inout) :: fcanmx      !< fractional coverage of class' 4 pfts
    real, dimension(ilg,ican), intent(inout) :: alvisc      !< visible albedo for class' 4 pfts
    real, dimension(ilg,ican), intent(inout) :: alnirc      !< near ir albedo for class' 4 pfts
    real, dimension(ilg,icc), intent(inout) :: pstemmass    !< stem mass from previous timestep, is value before fire. used by burntobare subroutine
    real, dimension(ilg,icc), intent(inout) :: pgleafmass   !< root mass from previous timestep, is value before fire. used by burntobare subroutine
    real, dimension(ilg,icc), intent(inout) :: flhrloss     !< fall or harvest loss for deciduous trees and crops, respectively \f$(kg C/m^2)\f$
    real, dimension(ilg,icc), intent(inout) :: flhrloss_ns  !< fall or harvest loss for deciduous trees and crops, respectively, non-structural \f$(kg C/m^2)\f$
    real, dimension(ilg,icc), intent(inout) :: flhrloss_s   !< fall or harvest loss for deciduous trees and crops, respectively, structural \f$(kg C/m^2)\f$
    real, dimension(ilg,icc), intent(inout) :: tracerFlHrLoss !< Tracer fall & harvest loss for bdl dcd plants and crops, respectively, 14C: \f$ng ^{14}C/m^2/d\f$.
    real, dimension(ilg,icc), intent(inout) :: stmhrlos     !< stem harvest loss for crops, \f$(kg C/m^2/d)\f$
    real, dimension(ilg,icc), intent(inout) :: rothrlos     !< root death as crops are harvested, \f$(kg C/m^2/d)\f$
    real, dimension(ilg,icc), intent(inout) :: tracerStmhrlos !< tracer stem harvest loss for crops, 14C: \f$ng ^{14}C/m^2/d\f$
    real, dimension(ilg,icc), intent(inout) :: tracerRothrlos !< tracer root death as crops are harvested, 14C: \f$ng ^{14}C/m^2/d\f$
    real, dimension(ilg,icc), intent(inout) :: grwtheff     !< growth efficiency. change in biomass per year per unit max. lai (\f$(kg C/m^2)\f$)/(m2/m2),
    !< for use in mortality subroutine
    real, dimension(ilg,icc), intent(inout) :: lystmmas     !< stem mass at the end of last year
    real, dimension(ilg,icc), intent(inout) :: lyrotmas     !< root mass at the end of last year
    real, dimension(ilg,icc), intent(inout) :: tymaxlai     !< this year's maximum lai
    real, dimension(ilg,icc), intent(inout) ::  geremort    !<
    real, dimension(ilg,icc), intent(inout) :: intrmort     !<
    real, dimension(ilg,icc), intent(inout) ::  burnvegf    !< per PFT fraction burned of that PFT's area
    real, dimension(ilg), intent(inout) :: popdin           !< population density \f$(people / km^2)\f$
    real, dimension(ilg), intent(inout) :: Cmossmas         !< moss biomass C pool (kgC/m2)
    real, dimension(ilg,ignd), intent(inout) :: litrmsmoss   !< moss litter C pool (kgC/m2)
    real, dimension(ilg,ignd), intent(inout) :: upMossSoilC  !< Moss soil C mass (non-peat), \f$kg C/m^2\f$
    ! ECFLAG: Intent should be in, and move it "in" section? real, dimension(ilg), intent(in) :: rmlmoss           !< moss maintainance respiration -daily averaged C fluxes rates (umol/m2/s)
    real, dimension(ilg), intent(inout) :: rmlmoss           !< moss maintainance respiration -daily averaged C fluxes rates (umol/m2/s)
    real, dimension(ilg), intent(inout) :: peatdep           !< peat depth (m)
    real, intent(inout) :: peatSoilC(ilg)                    !< Peat soil C pool (kgC/m2)
    real, dimension(ilg), intent(inout) :: Nmossmas          !< moss biomass N pool, \f$g N/m^2\f$
    real, dimension(ilg), intent(inout) :: litrmsmossN       !< moss litter N pool, \f$g N/m^2\f$
    real, dimension(ilg), intent(inout) :: upMossSoilN       !< Moss soil N mass (non-peat), \f$g N/m^2\f$
    real, dimension(ilg),intent(inout) :: bnf_moss           !< moss-associated biological nitrogen fixation, \f$g N m^{-2} day^{-1}\f$

    real, dimension(ilg, icc), intent(inout) :: lmaxt, smaxt, rmaxt  !< vars to find previous year C pool max

    real, intent(inout) :: tracerGLeafMass(:,:)      !< Tracer mass in the green leaf pool for each of the CTEM pfts, 14C: \f$ng ^{14}C/m^2\f$
    real, intent(inout) :: tracerBLeafMass(:,:)      !< Tracer mass in the brown leaf pool for each of the CTEM pfts, 14C: \f$ng ^{14}C/m^2\f$
    real, intent(inout) :: tracerStemMass(:,:)       !< Tracer mass in the stem for each of the CTEM pfts, 14C: \f$ng ^{14}C/m^2\f$
    real, intent(inout) :: tracerRootMass(:,:)       !< Tracer mass in the roots for each of the CTEM pfts, 14C: \f$ng ^{14}C/m^2\f$
    real, intent(inout) :: tracerLitrMass(:,:,:)     !< Tracer mass in the litter pool for each of the CTEM pfts + bareground and LUC products, 14C: \f$ng ^{14}C/m^2\f$
    real, intent(inout) :: tracerSoilCMass(:,:,:)    !< Tracer mass in the soil carbon pool for each of the CTEM pfts + bareground and LUC products, 14C: \f$ng ^{14}C/m^2\f$
    real, intent(inout) :: tracerMossCMass(:)        !< Tracer mass in moss biomass, 14C: \f$ng ^{14}C/m^2\f$
    real, intent(inout) :: tracerMossLitrMass(:)     !< Tracer mass in moss litter, 14C: \f$ng ^{14}C/m^2\f$

    real, dimension(ilg,icc), intent(inout) :: slai           !< storage/imaginary lai for phenology purposes
    real, dimension(ilg,icc), intent(inout) :: nppveg         !< NPP for individual pfts, (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(inout), dimension(ilg) :: timharvareagat      !< FLAG rsc  please add info for these
    real, intent(inout), dimension(nlat,nmos) :: timharvarearow
    real, intent(inout), dimension(ilg) :: tileAgegat            !< the age of the tile since the start of the run in months
    logical, intent(in) :: timberHarvest                       !< FLAG rsc  please add info for these
    logical, intent(in) :: dynamicTilingOn
    !
    !   outputs
    !
    real, dimension(ilg), intent(out) :: ch4soills          !< Methane uptake into the soil column \f$(mg CH_4 m^{-2} s^{-1})\f$
    real, dimension(ilg), intent(out) :: rml                !< Tile level leaf maintenance respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg), intent(out) :: gpp                !< Tile level gross primary productivity (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)

    real, dimension(ilg,icc), intent(out) :: veghght        !< vegetation height (meters)
    real, dimension(ilg,icc), intent(out) :: rootdpth       !< 99% soil rooting depth (meters) both veghght & rootdpth can be used as diagnostics
    !< to see how vegetation grows above and below ground, respectively
    real, dimension(ilg), intent(out) :: npp                !< Tile-level net primary productivity (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg), intent(out) :: nep                !< Tile-level net ecosystem productivity (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg), intent(out) :: hetrores           !< Tile-level heterotrophic respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg), intent(out) :: autores             !< Tile level autotrophic respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg), intent(out) :: soilresp           !< Soil respiration. This includes root respiration and respiration from
    !! litter and soil carbon pools. Note that soilresp is different from
    !! socres,which is respiration from the soil C pool.(\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg), intent(out) :: rm                 !< Tile level maintenance respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg), intent(out) :: rg                 !< Tile level growth respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg), intent(out) :: nbp                !< Tile level net biome productivity (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg), intent(out) :: dstcemls1          !< carbon emission losses due to disturbance (fire at present) from vegetation
    real, dimension(ilg), intent(out) :: litrfall           !< total litter fall (from leaves, stem, and root) due to all causes (mortality, turnover, and disturbance)
    real, dimension(ilg), intent(out) :: humiftrs           !< Transfer of humidified litter from litter to soil C pool (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg), intent(out) :: lucemcom           !< land use change (luc) related combustion emission losses, u-mol co2/m2.sec
    real, dimension(ilg), intent(out) :: lucltrin           !< luc related inputs to litter pool, u-mol co2/m2.sec
    real, dimension(ilg), intent(out) :: lucsocin           !< luc related inputs to soil c pool, u-mol co2/m2.sec
    real, dimension(ilg), intent(out) :: cproduct           !< LUC products in paper and furniture pools, Kg C/m2
    real, dimension(ilg), intent(out) :: fproductdecomp     !< Decomposition of LUC product pools, u-mol CO2-C/m2.s
    real, dimension(ilg), intent(out) :: cleaf              !< Grid averaged leaf biomass kg C/m2
    real, dimension(ilg), intent(out) :: cstem              !< Grid averaged stem biomass kg C/m2
    real, dimension(ilg), intent(out) :: croot              !< Grid averaged root biomass kg C/m2
    real, dimension(ilg), intent(out) :: nppleaf            !< Grid averaged leaf NPP, u-mol CO2-C/m2.s
    real, dimension(ilg), intent(out) :: nppstem            !< Grid averaged stem NPP, u-mol CO2-C/m2.s
    real, dimension(ilg), intent(out) :: npproot            !< Grid averaged root NPP, u-mol CO2-C/m2.s
    real, dimension(ilg), intent(out) :: nMineralNH4        !< Grid average NH4 mass (for CanESM), Kg N/m2
    real, dimension(ilg), intent(out) :: nMineralNO3        !< Grid average NO3 mass (for CanESM), Kg N/m2
    real, dimension(ilg), intent(out) :: nVeg               !< Grid average vegetation N mass (for CanESM), Kg N/m2
    real, dimension(ilg), intent(out) :: nLitter            !< Grid average litter N mass (for CanESM), Kg N/m2
    real, dimension(ilg), intent(out) :: nSoil              !< Grid average soil N mass mass (for CanESM), Kg N/m2
    real, dimension(ilg), intent(out) :: fBNF                !< grid average total biological nitrogen fixation \f$(kg N m^{-2} s^{-1}\f$
    real, dimension(ilg), intent(out) :: fNdep               !< grid average total nitrogen deposition \f$(kg N m^{-2} s^{-1}\f$
    real, dimension(ilg), intent(out) :: fNfert              !< grid average total nitrogen fertilization \f$(kg N m^{-2} s^{-1}\f$
    real, dimension(ilg), intent(out) :: fNgasNonFire        !< grid average total nitrogen flux excluding fire \f$(kg N m^{-2} s^{-1}\f$
    real, dimension(ilg), intent(out) :: nLeaf              !< Grid average vegetation N mass (for CanESM), Kg N/m2
    real, dimension(ilg), intent(out) :: nStem              !< Grid average vegetation N mass (for CanESM), Kg N/m2
    real, dimension(ilg), intent(out) :: nRoot              !< Grid average vegetation N mass (for CanESM), Kg N/m2
    real, dimension(ilg), intent(out) :: fNgasFire           !< grid average total nitrogen flux from fire only \f$(kg N m^{-2} s^{-1}\f$
    real, dimension(ilg), intent(out) :: fNgas               !< grid average total nitrogen flux from all sources including fire \f$(kg N m^{-2} s^{-1}\f$
    real, dimension(ilg), intent(out) :: fNnetmin            !< grid average total net mineralization \f$(kg N m^{-2} s^{-1}\f$
    real, dimension(ilg), intent(out) :: fNOx                !< grid average total NOx flux from all sources including fire \f$(kg N m^{-2} s^{-1}\f$
    real, dimension(ilg), intent(out) :: fNup                !< grid average pasive N updatake from the NH4 and NO3 pools \f$(kg N m^{-2} s^{-1}\f$
    real, dimension(ilg), intent(out) :: fNvegSoil           !< grid average total nitrogen flux from vegetation to soil [0 at present] \f$(kg N m^{-2} s^{-1}\f$
    real, dimension(ilg), intent(out) :: fN2o                !< grid average total N20 flux from all sources \f$(kg N m^{-2} s^{-1}\f$
    real, dimension(ilg), intent(out) :: lucemcomn           !< land use change (luc) related combustion emission losses of N (g N m^-2 day^-1)
    real, dimension(ilg), intent(out) :: lucltrinn           !< luc related inputs to litter pool of N (g N m^-2 day^-1)
    real, dimension(ilg), intent(out) :: lucsocinn           !< luc related inputs to soil N pool of N (g N m^-2 day^-1)
    real, dimension(ilg), intent(out) :: dstcemls3          !< carbon emission losses due to disturbance (fire at present) from litter pool
    real, dimension(ilg), intent(out) :: rms                !< Tile level stem maintenance respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg), intent(out) :: rmr                !< Tile level root maintenance respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg), intent(out) :: litres             !< Litter respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg), intent(out) :: socres             !< Soil carbon respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg,icc), intent(out) :: rmsveg         !< Maintenance respiration for stem for the CTEM pfts in u mol co2/m2. sec
    real, dimension(ilg,icc), intent(out) :: rmrveg         !< Maintenance respiration for root for the CTEM pfts in u mol co2/m2. sec
    real, dimension(ilg,icc), intent(out) :: rmlveg         !< Leaf maintenance respiration per PFT (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg,icc), intent(out) :: gppveg         !< Gross primary productivity per PFT (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg,icc), intent(out) :: rgveg          !< PFT level growth respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg,iccp1), intent(out) :: nepveg       !<
    real, dimension(ilg,iccp1), intent(out) :: nbpveg       !<
    real, dimension(ilg,iccp2,ignd), intent(out) :: ltresveg     !< fluxes for each pft: litter respiration for each pft + bare fraction
    real, dimension(ilg,iccp2,ignd), intent(out) :: scresveg     !< soil carbon respiration for the given sub-area in umol co2/m2.s, for ctem's pfts
    real, dimension(ilg,iccp1), intent(out) :: hetrsveg     !< Vegetation averaged litter and soil C respiration rates (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg,iccp2,ignd), intent(out) :: humtrsvg     !< transfer of humidified litter from litter to soil c pool per PFT. (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg,icc), intent(out) :: autoresveg     !< PFT level autotrophic respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg,icc), intent(out) :: litrfallveg    !<
    real, dimension(ilg,icc), intent(out) :: roottemp       !< root temperature, k
    real, dimension(ilg,icc), intent(out) :: emit_co2       !< carbon dioxide emitted from biomass burning in g of compound
    real, dimension(ilg,icc), intent(out) :: emit_co        !< carbon monoxide emitted from biomass burning in g of compound
    real, dimension(ilg,icc), intent(out) :: emit_ch4       !< methane emitted from biomass burning in g of compound
    real, dimension(ilg,icc), intent(out) :: emit_nmhc      !< non-methane hydrocarbons emitted from biomass burning in g of compound
    real, dimension(ilg,icc), intent(out) :: emit_h2        !< hydrogen gas emitted from biomass burning in g of compound
    real, dimension(ilg,icc), intent(out) :: emit_nox       !< nitrogen oxides emitted from biomass burning in g of compound
    real, dimension(ilg,icc), intent(out) :: emit_n2o       !< nitrous oxide emitted from biomass burning in g of compound
    real, dimension(ilg,icc), intent(out) :: emit_nh3       !< ammonia emitted from biomass burning in g of compound
    real, dimension(ilg,icc), intent(out) :: emit_pm25      !< particulate matter less than 2.5 um in diameter emitted from biomass burning in g of compound
    real, dimension(ilg,icc), intent(out) :: emit_tpm       !< total particulate matter emitted from biomass burning in g of compound
    real, dimension(ilg,icc), intent(out) :: emit_tc        !< total carbon emitted from biomass burning in g of compound
    real, dimension(ilg,icc), intent(out) :: emit_bc        !< black carbon emitted from biomass burning in g of compound
    real, dimension(ilg,icc), intent(out) :: emit_oc        !< organic carbon emitted from biomass burning in g of compound
    real, dimension(ilg,icc), intent(out) :: emit_so2       !< sulfur dioxide emitted from biomass burning in g of compound
    real, dimension(ilg,icc), intent(out) :: bterm_veg      !< biomass term for fire probabilty calc
    real, dimension(ilg), intent(out) :: lterm              !< lightning term for fire probabilty calc
    real, dimension(ilg,icc), intent(out) :: mterm_veg      !< moisture term for fire probabilty calc
    real, dimension(ilg), intent(out) :: ch4WetSpec            !<
    real, dimension(ilg), intent(out) :: wetfdyn            !<
    real, dimension(ilg), intent(out) :: ch4WetDyn            !<
    real, dimension(ilg,icc), intent(out) :: colrate        !< colonization rate (1/day)
    real, dimension(ilg,icc), intent(out) :: mortrate       !< mortality rate
    real, dimension(ilg,icc), intent(out) :: afrleaf        !< allocation fraction for leaves
    real, dimension(ilg,icc), intent(out) :: afrstem        !< allocation fraction for stem
    real, dimension(ilg,icc), intent(out) :: afrroot        !< allocation fraction for root
    real, dimension(ilg,icc), intent(out) :: wtstatus       !< soil water status used for calculating allocation fractions
    real, dimension(ilg,icc), intent(out) :: ltstatus       !< light status used for calculating allocation fractions
    real, dimension(ilg), intent(out) :: burnfrac           !< areal :: fraction burned due to fire for every grid cell (%)
    real, dimension(ilg,icc), intent(out) :: leaflitr       !< leaf litter fall rate (\f$\mu mol CO_2 m^{-2} s^{-1}\f$). this leaf litter does not
    !< include litter generated due to mortality/fire
    real, dimension(ilg,icc) :: leaflitr_ns    !<non-structural leaf litter \f$kg c/m^2\f$
    real, dimension(ilg,icc) :: leaflitr_s     !<structural leaf litter \f$kg c/m^2\f$
    real, dimension(ilg,icc), intent(out) :: smfunc_veg     !< soil moisture dependence on fire spread rate
    real, dimension(ilg,icc), intent(out) :: tltrleaf       !< total leaf litter fall rate (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg,icc), intent(out) :: tltrstem       !< total stem litter fall rate (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg,icc), intent(out) :: tltrroot       !< total root litter fall rate (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg,ican), intent(out) :: paicgat       !<
    real, dimension(ilg,ican), intent(out) :: slaicgat      !<
    real, dimension(ilg,icc), intent(out) :: vgbiomas_veg   !<
    real, dimension(ilg), intent(out) :: armoss             !< autotrophic respiration of moss (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, dimension(ilg), intent(out) :: nppmoss            !< net primary production of moss (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: reprocost(ilg,icc) !< Cost of making reproductive tissues, only non-zero when NPP is positive (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: blfltrdt(ilg,icc)  !< brown leaf litter generated due to disturbance \f$(kg c/m^2)\f$
    real, intent(out) :: glfltrdt(ilg,icc)  !< green leaf litter generated due to disturbance \f$(kg c/m^2)\f$
    real, intent(out) :: glcaemls(ilg,icc)  !< green leaf carbon emission disturbance losses, \f$kg c/m^2\f$
    real, intent(out) :: blcaemls(ilg,icc)  !< brown leaf carbon emission disturbance losses, \f$kg c/m^2\f$
    real, intent(out) :: rtcaemls(ilg,icc)  !< root carbon emission disturbance losses, \f$kg c/m^2\f$
    real, intent(out) :: stcaemls(ilg,icc)  !< stem carbon emission disturbance losses, \f$kg c/m^2\f$
    real, intent(out) :: ltrcemls(ilg,icc)  !< litter carbon emission disturbance losses, \f$kg c/m^2\f$
    real, intent(out) :: fFireCVeg(ilg,icc)    !< Live biomass combusted by fire, \f$kg c/m^2\f$
    real, intent(out) :: fFireLitter(ilg,icc)  !< above ground necromass combusted by fire (considers layer 1 litter as surface), \f$kg c/m^2\f$
    real, intent(out) :: fFireCSoil(ilg,icc)   !< soil carbon combusted by fire (considers soil c + litter(2:ignd)), \f$kg c/m^2\f$
    real, intent(out) :: ntchlveg(ilg,icc)  !< fluxes for each pft: Net change in leaf biomass, u-mol CO2/m2.sec
    real, intent(out) :: ntchsveg(ilg,icc)  !< fluxes for each pft: Net change in stem biomass, u-mol CO2/m2.sec
    real, intent(out) :: ntchrveg(ilg,icc)  !< fluxes for each pft: Net change in root biomass,
    !! the net change is the difference between allocation and
    !! autotrophic respiratory fluxes, u-mol CO2/m2.sec

    real, intent(out) :: leafns2s(ilg,icc)  !< Carbon flux from non-structural to structural leaf pool, \f$(kg C/m^2.day)\f$
    real, intent(out) :: stemns2s(ilg,icc)  !< Carbon flux from non-structural to structural stem pool, \f$(kg C/m^2.day)\f$
    real, intent(out) :: rootns2s(ilg,icc)  !< Carbon flux from non-structural to structural root pool, \f$(kg C/m^2.day)\f$
    
    !  Variables related to the Nitrogen cycle
    real, intent(in) :: nfertil(ilg)                !< N fertilizer \f$(g N m^{-2} cropland yr^{-1})\f$
    real, intent(in) :: ndeposit(ilg)               !< N deposition \f$(g N m^{-2} yr^{-1})\f$
    real, intent(inout) :: soilpH(ilg)                 !< Soil PH
    real, intent(inout) :: nh4_mass(ilg,iccp1)      !< ammonium mass for individual PFTs + bare (\f$g N/m^2\f$)
    real, intent(inout) :: no3_mass(ilg,iccp1)      !< nitrate mass for individual PFTs + bare (\f$g N/m^2\f$)
    real, intent(inout) :: ngleafmas(ilg,icc)       !< green leaf nitrogen mass for individual PFTs (\f$g N/m^2\f$)
    real, intent(inout) :: ngleafmas_ns(ilg,icc)    !< non-structural green leaf nitrogen mass for individual PFTs (\f$g N/m^2\f$)
    real, intent(inout) :: ngleafmas_s(ilg,icc)     !< structural green leaf nitrogen mass for individual PFTs (\f$g N/m^2\f$)
    real, intent(inout) :: nbleafmas(ilg,icc)       !< brown leaf nitrogen mass for individual PFTs (\f$g N/m^2\f$)
    real, intent(inout) :: nstemmass(ilg,icc)       !< stem nitrogen mass for individual PFTs (\f$g N/m^2\f$)
    real, intent(inout) :: nstemmass_ns(ilg,icc)    !< non-structural stem nitrogen mass for individual PFTs (\f$g N/m^2\f$)
    real, intent(inout) :: nstemmass_s(ilg,icc)     !< structural stem nitrogen mass for individual PFTs (\f$g N/m^2\f$)
    real, intent(inout) :: nrootmass(ilg,icc)       !< root nitrogen mass for individual PFTs (\f$g N/m^2\f$)
    real, intent(inout) :: nrootmass_ns(ilg,icc)    !< non-structural root nitrogen mass for individual PFTs (\f$g N/m^2\f$)
    real, intent(inout) :: nrootmass_s(ilg,icc)     !< structural root nitrogen mass for individual PFTs (\f$g N/m^2\f)
    real, intent(inout) :: nlitrmass(ilg,iccp2)     !< litter nitrogen mass for individual PFTs + bare (\f$g N/m^2\f$)
    real, intent(inout) :: soilnmas(ilg,iccp2)      !< soil organic nitrogen mass for individual PFTs + bare (\f$g N/m^2\f$)
    real, intent(inout) :: nvgbiomas_veg(ilg,icc)   !< vegetation nitrogen for each PFT, \f$g N m^{-2}\f$
    real,intent(out), dimension(ilg,iccp1) :: bnf_tot          !< total biological nitrogen fixation for PFTs + bare soil, \f$g N m^{-2} day^{-1}\f$
    real,intent(out), dimension(ilg,iccp1) :: bnf_free         !< free-living biological nitrogen fixation for PFTs + bare soil, \f$g N m^{-2} day^{-1}\f$
    real,intent(inout), dimension(ilg,icc) :: bnf_nat          !< natural symbiotic biological nitrogen fixation for PFTs + bare soil, \f$g N m^{-2} day^{-1}\f$
    real,intent(inout), dimension(ilg,icc) :: bnf_ant          !< anthropogenic symbiotic biological nitrogen fixation for PFTs + bare soil, \f$g N m^{-2} day^{-1}\f$
    real,intent(inout), dimension(ilg,icc) :: nstress          !< plant nitrogen stress
    real,intent(out), dimension(ilg,iccp1) :: nitrifveg      !< nitrification flux from NH4+ to NO3- for individual PFTs + bare soil, \f$g N m^{-2} day^{-1}\f$
    real,intent(out), dimension(ilg,iccp1)  :: no_nitveg     !< nitrification NO loss for individual PFTs + bare soil, \f$g N m^{-2} day^{-1}\f$
    real,intent(out), dimension(ilg,iccp1) :: no_denitveg    !< denitrification NO loss for individual PFTs + bare soil, \f$g N m^{-2} day^{-1}\f$
    real,intent(out), dimension(ilg,iccp1) :: no_nitdenitveg !< total NO loss from denitrification and nitrification for individual PFTs + bare soil, \f$g N m^{-2} day^{-1}\f$
    real,intent(out), dimension(ilg,iccp1) :: n2o_nitveg     !< nitrification N2O loss for individual PFTs + bare soil, \f$g N m^{-2} day^{-1}\f$
    real,intent(out), dimension(ilg,iccp1) :: n2o_denitveg   !< denitrification N2O loss for individual PFTs + bare soil, \f$g N m^{-2} day^{-1}\f$
    real,intent(out), dimension(ilg,iccp1) :: n2o_nitdenitveg!< total N2O loss from denitrification and nitrification for individual PFTs + bare, \f$g N m^{-2} day^{-1}\f$
    real,intent(out), dimension(ilg,iccp1) :: n2_denitveg    !< N2 loss from denitrification for individual PFTs + bare soil, \f$g N m^{-2} day^{-1}\f$
    real,intent(out), dimension(ilg,iccp1) :: nvolveg        !< NH3 volatilization loss from NH4+, \f$g N m^{-2} day^{-1}\f$
    real,intent(out), dimension(ilg,iccp1) :: nleachveg      !< nitrogen leaching from NO3-, \f$g N m^{-2} day^{-1}\f$
    real,intent(out), dimension(ilg) :: fNleach      !< grid averaged total N leaching \f$(kg N m^{-2} s^{-1}\f$
    real,intent(out), dimension(ilg) :: fNvol      !< grid averaged total N volatilization from ammonium pool (\f$(kg N m^{-2} s^{-1}\f$)
    real,intent(out), dimension(ilg,iccp1) :: appl_fert      !< daily nitrogen fertilizer for individual PFTs + bare soil, \f$g N m^{-2} day^{-1}\f$
    real,intent(out), dimension(ilg,iccp1) :: ndep_nh4       !< daily nitrogen deposition into NH4+ for individual PFTs + bare soil, \f$g N m^{-2} day^{-1}\f$
    real,intent(out), dimension(ilg,iccp1) :: ndep_no3       !< daily nitrogen deposition into NO3- for individual PFTs + bare soil, \f$g N m^{-2} day^{-1}\f$
    real,intent(out), dimension(ilg,icc) :: ndemandveg_wp_npp!< whole plant nitrogen demand based on NPP for each PFT, \f$g N m^{-2} day^{-1}\f$
    real,intent(out), dimension(ilg,icc) :: nuptakeveg_p_nh4 !< passive nh4+ uptake for each PFT, \f$g N m^{-2} day^{-1}\f$
    real,intent(out), dimension(ilg,icc) :: nuptakeveg_p_no3 !< passive no3- uptake for each PFT, \f$g N m^{-2} day^{-1}\f$
    real,intent(out), dimension(ilg,icc) :: nuptakeveg_a_actl_nh4 !< actual active nh4+ uptake for each PFT, \f$g N m^{-2} day^{-1}\f$
    real,intent(out), dimension(ilg,icc) :: nuptakeveg_a_actl_no3 !< actual active no3- uptake for each PFT, \f$g N m^{-2} day^{-1}\f$
    real,intent(out), dimension(ilg,icc) :: nuptakeveg       !< total N uptake (active+passive, NH4+NO3) for each PFT, \f$g N m^{-2} day^{-1}\f$
    real,intent(out), dimension(ilg,icc) :: nleafns2sveg     !< nitrogen flux from non-structural to structural leaf nitrogen, \f$g N m^{-2} day^{-1}\f$
    real,intent(out), dimension(ilg,icc) :: nstemns2sveg     !< nitrogen flux from non-structural to structural stem nitrogen, \f$g N m^{-2} day^{-1}\f$
    real,intent(out), dimension(ilg,icc) :: nrootns2sveg     !< nitrogen flux from non-structural to structural root nitrogen, \f$g N m^{-2} day^{-1}\f$
    real,intent(out), dimension(ilg,icc) :: nallocveg_l      !< nitrogen allocation to leaf for each PFT, \f$g N m^{-2} day^{-1}\f$
    real,intent(out), dimension(ilg,icc) :: nallocveg_s      !< nitrogen allocation to stem for each PFT, \f$g N m^{-2} day^{-1}\f$
    real,intent(out), dimension(ilg,icc) :: nallocveg_r      !< nitrogen allocation to root for each PFT, \f$g N m^{-2} day^{-1}\f$
    real,intent(out), dimension(ilg,icc) :: nresorpedveg_s   !< resorbed nitrogen from leaf allocated to stem, \f$g N m^{-2} day^{-1}\f$
    real,intent(out), dimension(ilg,icc) :: nresorpedveg_r   !< resorbed nitrogen from leaf allocated to root, \f$g N m^{-2} day^{-1}\f$
    real,intent(out), dimension(ilg,icc) :: nre_allocveg_s2l !< non-structural nitrogen reallocation from stem to leaf during leaf out, \f$g N m^{-2} day^{-1}\f$
    real,intent(out), dimension(ilg,icc) :: nre_allocveg_r2l !< non-structural nitrogen reallocation from root to leaf during leaf out, \f$g N m^{-2} day^{-1}\f$
    real,intent(out), dimension(ilg,icc) :: nlitrveg_l       !< total leaf nitrogen litter fall rate, \f$g N m^{-2} day^{-1}\f$
    real,intent(out), dimension(ilg,icc) :: nlitrveg_s       !< total stem nitrogen litter fall rate, \f$g N m^{-2} day^{-1}\f$
    real,intent(out), dimension(ilg,icc) :: nlitrveg_r       !< total root nitrogen litter fall rate, \f$g N m^{-2} day^{-1}\f$
    real,intent(out), dimension(ilg,icc) :: nlitrveg         !< total nitrogen litter fall rate, \f$g N m^{-2} day^{-1}\f$
    real,intent(out), dimension(ilg,icc) :: gl2bl_grass_nflux!< nitrogen flux from green to brown leaf nitrogen for each PFT, \f$g N m^{-2} day^{-1}\f$
    real,intent(out), dimension(ilg,icc) :: c2nveg_l         !< simulated C:N ratio for leaf, \f$g C g N^{-1}\f$
    real,intent(out), dimension(ilg,icc) :: c2nveg_s         !< simulated C:N ratio for stem, \f$g C g N^{-1}\f$
    real,intent(out), dimension(ilg,icc) :: c2nveg_r         !< simulated C:N ratio for root, \f$g C g N^{-1}\f$
    real,intent(out), dimension(ilg,icc) :: c2nveg_wp        !< simulated C:N ratio for the whole plant, \f$g C g N^{-1}\f$
    real,intent(out), dimension(ilg,iccp1) :: c2nveg_litr    !< simulated C:N ratio for detritus, \f$g C g N^{-1}\f$
    real,intent(out), dimension(ilg,iccp1) :: c2nveg_humus   !< simulated C:N ratio for SOM, \f$g C g N^{-1}\f$
    real,intent(out), dimension(ilg,iccp1) :: nhumtrsveg     !< nitrogen flux from detritus to SOM for individual PFTs + bare soil, \f$g N m^{-2} day^{-1}\f$
    real,intent(out), dimension(ilg,iccp1) :: nmineralveg_litr  !< nitrogen mineralization from detritus, \f$g N m^{-2} day^{-1}\f$
    real,intent(out), dimension(ilg,iccp1) :: nmineralveg_humus !< nitrogen mineralization from SOM, \f$g N m^{-2} day^{-1}\f$
    real,intent(out), dimension(ilg,iccp1) :: netnmineralveg       !< nitrogen mineralization, \f$g N m^{-2} day^{-1}\f$
    real,intent(out), dimension(ilg,iccp1) :: nimmobilveg_nh4   !< nitrogen immobilization from NH4+ to soil, \f$g N m^{-2} day^{-1}\f$
    real,intent(out), dimension(ilg,iccp1) :: nimmobilveg_no3   !< nitrogen immobilization from NO3- to soil, \f$g N m^{-2} day^{-1}\f$
    real,intent(out), dimension(ilg,iccp1) :: fNnetlandveg   !< net terrestrial nitrogen flux, \f$g N m^{-2} day^{-1}\f$
    real,intent(out), dimension(ilg,icc) :: redcoeff_vcmax   !< reduction coefficient of V_c,max passed to the photosynthesis subroutine
    real, intent(out), dimension(ilg,icc) :: fiengat          !< fire plume energy per PFT for plume height calculation \f$[J]\f$
    real, intent(out), dimension(ilg)     :: fabgat           !< total area burned due to fire\f$[km^2]\f$
    real, intent(out), dimension(ilg,icc) :: nmfigat          !< number of fires per PFT in each grid cell
    real,intent(out), dimension(ilg)     :: redcoeff_vcmaxMoss   !< reduction coefficient of V_c,max passed to the photosynthesis subroutine for moss

    real, intent(out) :: nepCMIP(ilg)
    real, intent(out) :: drgtstrs(ilg,icc) !< soil dryness factor for pfts
    real, intent(out) :: betadrgt(ilg,ignd)!< dryness term for soil layers

    real, intent(in) :: QACGAT(ilg)     !< Specific humidity
    real, intent(in) :: PRESGAT(ilg)    !<atmospheric pressure
    real, intent(out) :: isi(ilg)     !< Initial Spread Index (FWI-based fire scheme)
    real, intent(out) :: bui(ilg)     !< Build-Up Index (FWI-based fire scheme) 
    real, intent(out) :: fwi(ilg)     !< Fire Weather Index (FWI-based fire scheme)
    real, intent(inout) :: ffmc(ilg)     !<  Fine Fuel Moisture Code (FWI-based fire scheme)
    real, intent(inout) :: dmc(ilg)     !<  Duff Moisture Code (FWI-based fire scheme)
    real, intent(inout) :: dc(ilg)     !< Drought Code (FWI-based fire scheme)
    integer, intent(in) :: IMONTH        !< Month of the year simulation is in.

    real, dimension(ilg), intent(in) :: NoonTempc   !< FLAG rsc  please add info for these
    real, dimension(ilg), intent(in) :: NoonPrecip
    real, dimension(ilg), intent(in) :: NoonWindtot
    real, dimension(ilg), intent(in) :: NoonRhum

    ! ---------------------------------------------
    ! Local variables:

    integer :: i, j, k, n, m
    integer :: sort(icc)
    !character(8) :: pftkind

    real :: yesfrac_comp(ilg,icc) !<
    real :: fc(ilg)  !< Fraction of grid cell that is covered by canopy
    real :: fg(ilg)  !< Fraction of grid cell that is bare ground.
    real :: pglfmass(ilg,icc)  !< Prior timestep green leaf mass for each of the CTEM pfts, \f$(kg C/m^2)\f$
    real :: pblfmass(ilg,icc)  !< Prior timestep brown leaf mass for each of the CTEM pfts, \f$(kg C/m^2)\f$
    real :: pstemass(ilg,icc)  !< Prior timestep stem mass for each of the CTEM pfts, \f$(kg C/m^2)\f$
    real :: protmass(ilg,icc)  !< Prior timestep root mass for each of the CTEM pfts, \f$(kg C/m^2)\f$
    real :: plitmass(ilg,iccp2)!< Prior timestep litter mass for each of the CTEM pfts + bare + LUC product pools, \f$(kg C/m^2)\f$
    real :: psocmass(ilg,iccp2)!< Prior timestep soil carbon mass for each of the CTEM pfts + bare + LUC product pools, \f$(kg C/m^2)\f$
    real :: pvgbioms(ilg)      !< Prior timestep
    real :: pgavltms(ilg)      !< Prior timestep
    real :: pgavscms(ilg)      !< Prior timestep
    real :: pheanveg(ilg,icc) !<
    real :: rootlitr(ilg,icc) !< root litter \f$(kg C/m^2)\f$
    real :: rootlitr_ns(ilg,icc) !<
    real :: rootlitr_s(ilg,icc)  !<
    real :: stemlitr(ilg,icc) !< stem litter \f$(kg C/m^2)\f$
    real :: stemlitr_ns(ilg,icc) !< non-structural stem litter \f$(kg C/m^2)\f$
    real :: stemlitr_s(ilg,icc)  !< structural stem litter \f$(kg C/m^2)\f$
    real :: stemltrm(ilg,icc) !<
    real :: stemltrm_ns(ilg,icc) !<
    real :: stemltrm_s(ilg,icc)  !<
    real :: rootltrm(ilg,icc) !<
    real :: rootltrm_ns(ilg,icc) !<
    real :: rootltrm_s(ilg,icc)  !<
    real :: glealtrm(ilg,icc) !<
    real :: glealtrm_ns(ilg,icc) !<
    real :: glealtrm_s(ilg,icc)  !<
    real :: stemltrm_n(ilg,icc) !<
    real :: stemltrm_ns_n(ilg,icc) !<
    real :: stemltrm_s_n(ilg,icc)  !<
    real :: rootltrm_n(ilg,icc) !<
    real :: rootltrm_ns_n(ilg,icc) !<
    real :: rootltrm_s_n(ilg,icc)  !<
    real :: glealtrm_n(ilg,icc) !<
    real :: glealtrm_ns_n(ilg,icc) !<
    real :: glealtrm_s_n(ilg,icc)  !<
    real :: stemltdt(ilg,icc)  !<
    real :: nstemltdt(ilg,icc) !<
    real :: rootltdt(ilg,icc)  !<
    real :: nrootltdt(ilg,icc) !<
    real :: dscemlv1(ilg,icc)  !<
    real :: dscemlv2(ilg,icc)  !<
    real :: add2allo(ilg,icc)  !<
    real :: repro_cost_g(ilg)  !< Tile-level cost of making reproductive tissues, only non-zero when NPP is positive (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real :: litresmoss(ilg,ignd)  !< moss litter respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real :: socres_peat(ilg)   !< heterotrophic repsiration from peat soil (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real :: socres_moss(ilg,ignd) !< heterotrophic repsiration from moss-derived soil C (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real :: resoxic(ilg)       !< oxic respiration from peat soil (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real :: resanoxic(ilg)     !< anoxic respiration from peat soil (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real :: litrfallmoss(ilg)  !< moss litter fall (kgC/m2/timestep)
    real :: ltrestepmoss(ilg,ignd)  !< litter respiration from moss (kgC/m2/timestep)
    real :: humstepmoss(ilg,ignd)   !< moss humification (transfer from litter to SOM) (kgC/m2/timestep)
    real :: pCmossmas(ilg)     !< moss biomass C at the previous time step (kgC/m2)
    real :: plitrmsmoss(ilg)   !< moss litter C at the previous time step (kgC/m2)
    real :: nppmosstep(ilg)    !< moss npp (kgC/m2/timestep)
    real :: socrestep(ilg)     !< heterotrophic respiration from soil (kgC/m2/timestep)
    real :: nblfltrdt(ilg,icc) !< N brown leaf litter generated due to disturbance \f$(g N/m^2)\f$
    real :: nglfltrdt(ilg,icc) !< N green leaf litter generated due to disturbance \f$(g N/m^2)\f$
    real :: afrleaf_redu_nh(ilg)  !< coeff. for reducing afrleaf after summer solstice in northern hemisphere
    real :: afrleaf_redu_sh(ilg)  !< coeff. for reducing afrleaf after summer solstice in southern hemisphere
    real :: deficit
    real :: lygleafmasmax(ilg,icc)!< last year maximum of the gleafmas
    real :: lystemmassmax(ilg,icc)!< last year maximum of the stemmass
    real :: lyrootmassmax(ilg,icc)!< last year maximum of the rootmass
    real :: ccost_nuptake(ilg,icc)!< carbon cost associated to active nitrogen uptake (\f$kg C/m^2\f$)
    real :: lfthrs(ilg,icc)       !< threshold lai for finding leaf status

    real :: rmsTracer(ilg,icc)                !< Tracer maintenance respiration for stem (14C: \f$amol ^{14}CO_2 m^{-2} s^{-1}\f$)
    real :: rmrTracer(ilg,icc)                !< Tracer maintenance respiration for root (14C: \f$amol ^{14}CO_2 m^{-2} s^{-1}\f$)
    real :: ltResTracer(ilg,iccp2,ignd)       !< Tracer litter respiration for each pft + bare fraction (14C: \f$amol ^{14}CO_2 m^{-2} s^{-1}\f$)
    real :: sCResTracer(ilg,iccp2,ignd)       !< Tracer soil carbon respiration for the given sub-area, 14C: \f$amol ^{14}CO_2 m^{-2} s^{-1}\f$
    real :: litResMossTracer(ilg)             !< Tracer moss litter respiration (14C: \f$amol ^{14}CO_2 m^{-2} s^{-1}\f$)
    real :: soCResPeatTracer(ilg)             !< Tracer heterotrophic repsiration from peat soil (14C: \f$amol ^{14}CO_2 m^{-2} s^{-1}\f$)
    real :: tracerNPP(ilg,icc)                !< tracer NPP for individual pfts, (14C: \f$amol ^{14}CO_2 m^{-2} s^{-1}\f$)
    real :: tracerLeafLitr(ilg,icc)           !< Tracer leaf litter generated by normal turnover, cold
                                              !! and drought stress,and leaf fall/harvest,14C:  \f$ng ^{14}C/m^2\f$
    real :: tracerStemLitr(ilg,icc)           !< Tracer stem litter14C:  \f$ng ^{14}C/m^2\f$
    real :: tracerRootLitr(ilg,icc)           !< Tracer root litter14C:  \f$ng ^{14}C/m^2\f$
    real :: tracerReproCost(ilg,icc)          !< Tracer cost of making reproductive tissues, only non-zero when
                                              !! NPP is positive (14C: \f$amol ^{14}CO_2 m^{-2} s^{-1}\f$)
    real :: tracerStemMort(ilg,icc)           !< Tracer stem litter from mortality, 14C: \f$ng ^{14}C/m^2\f$
    real :: tracerRootMort(ilg,icc)           !< Tracer stem litter from mortality, 14C: \f$ng ^{14}C/m^2\f$
    real :: tracerGLeafMort(ilg,icc)          !< Tracer stem litter from mortality, 14C: \f$ng ^{14}C/m^2\f$
    real :: tracerValue(ilg)                  !< Tracer CO2 value updated from the prepTracer subroutine, (14C: 1E12* 14C/C ratio)
    real :: tracerRML(ilg,icc)                !< Tracer leaf maintenance respiration (14C: \f$amol ^{14}CO_2 m^{-2} s^{-1}\f$)
    real :: tracerGPP(ilg,icc)                !< Tracer GPP (14C: \f$amol ^{14}CO_2 m^{-2} s^{-1}\f$)
    real :: glniemls(ilg,icc)                 !< green leaf nitrogen emission disturbance losses, \f$g N/m^2\f$
    real :: rtniemls(ilg,icc)                 !< root nitrogen emission disturbance losses, \f$g N/m^2\f$
    real :: stniemls(ilg,icc)                 !< stem nitrogen emission disturbance losses, \f$g N/m^2\f$
    real :: ltrnemls(ilg,icc)                 !< litter nitrogen emission disturbance losses, \f$g N/m^2\f$
    real :: blniemls(ilg,icc)                 !< brown leaf nitrogen emission disturbance losses, \f$g N/m^2\f$
    real :: nh4_mass_g(ilg)                   !< grid avg. NH4+ mass (\f$g N/m^2\f$)
    real :: no3_mass_g(ilg)                   !< grid avg. NO3- mass (\f$g N/m^2\f$)
    real :: nlitrmass_g(ilg)                  !< grid avg. litter nitrogen mass (\f$g N/m^2\f$)
    real :: soilnmas_g(ilg)                   !< grid avg. soil organic nitrogen mass (\f$g N/m^2\f$)
    real :: nvgbiomas(ilg)

    real :: gl2bl_grass_cflux(ilg,icc)        !< carbon flux from gleafmas to bleafmass, to be used in nitrogen cycle (\f$kg C/m^2 day\f$)
    real :: bl2ltr_grass_cflux(ilg,icc)       !< litter carbon flux from bleafmas to litrmass, to be used in nitrogen cycle (\f$kg C/m^2 day\f$)
    real :: barefrac(ilg)                     !< bare fraction of each grid cell
    real :: dstcemls2(ilg)                    !< grid ave. carbon emission losses due to disturbance, total,  \f$\mu mol CO_2 m^{-2} s^{-1}\f$
    real :: pcproduct_before_luc(ilg)         !< grid averaged previous cProduct pool, Kg C/m2
    real :: pool_based_nbp(ilg)               !< grid averaged pool based NBP -u mol CO2-C/m2.s
    real :: pool_based_npp(ilg)               !< grid averaged pool based NPP -u mol CO2-C/m2.s
    real :: pvgbioms_before_luc(ilg) 
    real :: pgavltms_before_luc(ilg) 
    real :: pgavscms_before_luc(ilg) 
    real :: scaled_nppleaf, scaled_npproot, scaled_nppstem

#if !defined without_agcm_
    real :: diff1, t1, t2, t3, t4
#endif

    associate( &
    turbationON => c_switch%turbationON                 & !< logical: true if turbation is on. 
    )

    !> initialize required arrays to zero
    do i = il1,il2
       barefrac(i) = 1.0
       nvgbiomas(i) = 0.0
       nlitrmass_g(i) = 0.0
       soilnmas_g(i) = 0.0
       nh4_mass_g(i) = 0.0
       no3_mass_g(i) = 0.0
       plitrmsmoss(i) = 0.0
       fNleach(i) = 0.0
       fNvol(i) = 0.0
       do j = 1,icc
           vgbiomas_veg(i,j) = 0.0
           nvgbiomas_veg(i,j) = 0.0
       end do
    end do

    !> Calculate grid averaged nitrogen pool sizes for competition:
    do i = il1,il2
      do j = 1,icc
        nvgbiomas(i) = nvgbiomas(i) + fcancmx(i,j) * (ngleafmas(i,j) + nbleafmas(i,j) + &
          nrootmass(i,j) + nstemmass(i,j))
        nlitrmass_g(i) = nlitrmass_g(i) + nlitrmass(i,j) * fcancmx(i,j)
        soilnmas_g(i) = soilnmas_g(i) + soilnmas(i,j) * fcancmx(i,j)
        nh4_mass_g(i) = nh4_mass_g(i) + nh4_mass(i,j) * fcancmx(i,j)
        no3_mass_g(i) = no3_mass_g(i) + no3_mass(i,j) * fcancmx(i,j)
        barefrac(i) = barefrac(i) - fcancmx(i,j)
      end do
      nlitrmass_g(i) = nlitrmass_g(i) + nlitrmass(i,iccp1) * barefrac(i) 
      soilnmas_g(i) = soilnmas_g(i) + soilnmas(i,iccp1) * barefrac(i) 
      nh4_mass_g(i) = nh4_mass_g(i) + nh4_mass(i,iccp1) * barefrac(i)
      no3_mass_g(i) = no3_mass_g(i) + no3_mass(i,iccp1) * barefrac(i)
    end do

    

    ! ---------------------------------------------------------------------------
    ! Begin calculations

    !> Generate the sort index for correspondence between CTEM pfts and the
    !>  values in the parameter vectors
    sort = genSortIndex()

    if (useTracer /= 'None') then
      !> Set up the tracer for today
      call prepTracer(il1, il2, ilg, tracerCO2, & ! In
                      stemmass_s, stemmass_ns, rootmass_s, rootmass_ns, litrmass, & ! In  
                      gleafmas_s, gleafmas_ns, bleafmas, soilcmas, peatlandType, & ! In
                      tracerValue, & ! Out
                      tracerStemMass, tracerRootMass, tracerLitrMass, & ! InOut
                      tracerGLeafMass, tracerBLeafMass, tracerSoilCMass) ! InOut
    end if

    if (PFTCompetition) then

      !> Calculate bioclimatic parameters for estimating pfts existence
      call  bioclim(iday, ta, precip, netrad, &
                    1, il2, ilg, leapnow, &
                    tcurm, srpcuryr, dftcuryr, inibioclim, &
                    tmonth, anpcpcur, anpecur, gdd5cur, &
                    surmncur, defmncur, srplscur, defctcur, &
                    twarmm, tcoldm, gdd5, aridity, &
                    srplsmon, defctmon, anndefct, annsrpls, &
                    annpcp, dry_season_length)

      if (inibioclim) then

        !> If first day of year then based on updated bioclimatic parameters
        !! find if pfts can exist or not.
        !! If .not. inibioclim then it is the first year of a run that you do not have the
        !! climatological means already in the CTM file. After one
        !! year inibioclim is set to true and the climatological means
        !! are used from the first year.
        !!
        call existence(iday, 1, il2, ilg, &
                       sort, twarmm, tcoldm, &
                       gdd5, aridity, srplsmon, defctmon, &
                       anndefct, annsrpls, annpcp, pftexist, &
                       dry_season_length)

        !> Call competition subroutine which on the basis of previous day's
        !! NPP estimates changes in fractional coverage of pfts
        call competition(iday, 1, il2, ilg, nppveg, dofire, leapnow, useTracer, & ! In
                         pftexist, geremort, intrmort, pgleafmass, rmatctem, & ! In
                         grclarea, ailcg, lfstatus, burnvegf, sort, pstemmass, & ! In
                         gleafmas, gleafmas_ns, gleafmas_s, bleafmas, & ! In/Out
                         stemmass, stemmass_ns, stemmass_s, & ! In/Out
                         rootmass, rootmass_ns, rootmass_s, & ! In/Out
                         litrmass, soilcmas, fcancmx, fcanmx, & ! In/Out
                         tracerGLeafMass, tracerBLeafMass, tracerStemMass, tracerRootMass, & ! In/Out
                         tracerLitrMass, tracerSoilCMass, & ! In/Out
                         vgbiomas, gavgltms, gavgscms, bmasveg, & ! In/Out
                         add2allo, colrate, mortrate, & ! Out
                         ngleafmas, ngleafmas_ns, ngleafmas_s, nbleafmas, & ! In/Out
                         nstemmass, nstemmass_ns, nstemmass_s, nrootmass, & ! In/Out
                         nrootmass_ns, nrootmass_s, nlitrmass, soilnmas, & ! In/Out
                         nh4_mass, no3_mass, nvgbiomas, nlitrmass_g, & ! In/Out
                         soilnmas_g, nh4_mass_g, no3_mass_g) ! In/Out

      end if ! inibioclim
    end if  ! if (PFTCompetition)
   
    do i = il1,il2
       pcproduct_before_luc(i) = litrmass(i,iccp2,1) + soilcmas(i,iccp2,1) !< previous cProduct pool, litrmass(i,iccp2,1) + soilcmas(i,iccp2,1)
                                                                ! have not been updated by LUC yet 
                                                                ! These are grid-averaged because LUC product pools are assumed to
                                                                ! occur over the entire grid cell
       pvgbioms_before_luc(i)=vgbiomas(i) ! from previous day coming i ! from previous day coming in
       pgavltms_before_luc(i)=gavgltms(i)
       pgavscms_before_luc(i)=gavgscms(i)
    enddo

    !>
    !! If landuse is on, then implelement luc, change fractional coverages,
    !! move biomasses around, and estimate luc related combustion emission losses.
    !!
    if (lnduseon) then

      do j = 1,icc
        do i = il1,il2
          yesfrac_comp(i,j) = fcancmx(i,j)
        end do
      end do

      call luc(il1, il2, ilg, PFTCompetition, leapnow, useTracer, & ! In
               grclarea, iday, todfrac, yesfrac_comp, .true., & ! In
               pfcancmx, nfcancmx, &! In/Out
               gleafmas, gleafmas_ns, gleafmas_s, & ! In/Out
               bleafmas, stemmass, stemmass_ns, stemmass_s, & ! In/Out
               rootmass, rootmass_ns, rootmass_s, litrmass, & ! In/Out
               soilcmas, vgbiomas, nvgbiomas, gavgltms, nlitrmass_g, gavgscms, soilnmas_g, nh4_mass_g, no3_mass_g, & ! In/Out
               fcancmx, fcanmx, ngleafmas, ngleafmas_ns, & ! In/Out
               ngleafmas_s, nbleafmas, nstemmass, nstemmass_ns, & ! In/Out
               nstemmass_s, nrootmass, nrootmass_ns, nrootmass_s, & ! In/Out
               nlitrmass, soilnmas, nh4_mass, no3_mass, & ! In / Out
               tracerLitrMass, tracerSoilCMass, & ! In/Out
               tracerGLeafMass, tracerBLeafMass, tracerStemMass, tracerRootMass, & ! In / Out
               lucemcom, lucltrin, lucsocin, lucemcomn, lucltrinn, lucsocinn) ! Out
    else
      lucemcom = 0.
      lucltrin = 0.
      lucsocin = 0.
      lucemcomn = 0.
      lucltrinn = 0.
      lucsocinn = 0.
    end if ! lnduseon

    ! If dynamic tiling is active or the trackTileAge flag is set (i.e. during dynamic tiling spin up)
    ! increment tile age at the CTEM time step
    if(dynamicTilingOn .or. trackTileAge) then
      call indexTileAge(FAREGAT,tileAgegat, il1, il2, leapnow)
    end if

    ! If timber harvesting is used (with or without tiling)
    if (timberHarvest) then
      call harvestTile(Ncycle_on, il1, il2, timharvareagat, timharvarearow, tileAgegat, &
                    gleafmas, gleafmas_ns, gleafmas_s, &
                    bleafmas, stemmass, stemmass_ns, stemmass_s, &
                    rootmass, rootmass_ns, rootmass_s, &
                    litrmass, soilcmas, vgbiomas, gavgltms, &
                    gavgscms, fcancmx, ngleafmas, ngleafmas_ns, &
                    ngleafmas_s, nbleafmas, nstemmass, nstemmass_ns, &
                    nstemmass_s, nrootmass, nrootmass_ns, nrootmass_s, &
                    nlitrmass, soilnmas, &
                    lucemcom, lucltrin, lucsocin, lucemcomn, lucltrinn, lucsocinn, faregat, rmatctem, &
                    trackTileAge, dynamicTilingOn, lfstatus, pandays, flhrloss, flhrloss_ns, flhrloss_s) 
    end if
    
    !> Store green and brown leaf, stem, and root biomass, and litter and
    !! soil c pool mass in arrays. knowing initial sizes of all pools and
    !! final sizes at the end of this subroutine, we check for conservation of mass.
    plitmass = 0.0
    psocmass = 0.0
    do j = 1,icc
      do i = il1,il2
        pglfmass(i,j) = gleafmas(i,j)    !< green leaf mass from last time step
        pblfmass(i,j) = bleafmas(i,j)    !< brown leaf mass from last time step
        pstemass(i,j) = stemmass(i,j)    !< stem mass from last time step
        protmass(i,j) = rootmass(i,j)    !< root mass from last time step
        do k = 1,ignd ! FLAG at this stage keep as per pft and per tile. JM Feb8 2016.
          plitmass(i,j)=plitmass(i,j) + litrmass(i,j,k)    ! litter mass from last time step
          psocmass(i,j)=psocmass(i,j) + soilcmas(i,j,k)    ! soil c mass from last time step
          !plitmasspl(i,j,k)= litrmass(i,j,k)
          !psocmasspl(i,j,k)= soilcmas(i,j,k)
        end do
      end do ! loop 140
    end do ! loop 130
    !
    do i = il1,il2
      pvgbioms(i) = vgbiomas(i)          !< vegetation biomass from last time step
      vgbiomas(i) = 0.0

      cproduct(i) = 0.0                  !< initialize to zero, update later in loop 198
      fproductdecomp(i) = 0.0            !< initialize to zero, update later in loop 198
      cleaf(i) = 0.0            !< initialize to zero, update later in loop 199
      cstem(i) = 0.0            !< initialize to zero, update later in loop 199
      croot(i) = 0.0            !< initialize to zero, update later in loop 199
      nppleaf(i) = 0.0            !< initialize to zero, update later in loop 199
      npproot(i) = 0.0            !< initialize to zero, update later in loop 199
      nppstem(i) = 0.0            !< initialize to zero, update later in loop 199

      nMineralNH4(i) = 0.0            !< initialize all to zero, update later right after N cycle is called
      nMineralNO3(i) = 0.0
      nVeg(i) = 0.0
      nLitter(i) = 0.0
      nSoil(i) = 0.0
      nLeaf(i) = 0.0
      nStem(i) = 0.0
      nRoot(i) = 0.0
      fBNF(i) = 0.0
      fNdep(i) = 0.0
      fNfert(i) = 0.0
      fNgasNonFire(i) = 0.0
      fNgasFire(i) = 0.0
      fNgas(i) = 0.0
      fNnetmin(i) = 0.0
      fNOx(i) = 0.0
      fNup(i) = 0.0
      fNvegSoil(i) = 0.0
      fN2o(i) = 0.0             !< end initialize all to zero, update later right after N cycle is called

      pgavltms(i) = gavgltms(i)          !< litter mass from last time step
      pgavscms(i) = gavgscms(i)          !< soil c mass from last time step
      do j = iccp1,iccp2 ! do over the bare fraction and the LUC pool
        do k = 1,ignd ! FLAG at this stage keep as per pft and per tile. JM Feb8 2016.
          plitmass(i,j)=plitmass(i,j) + litrmass(i,j,k)  ! litter mass over bare fraction
          psocmass(i,j)=psocmass(i,j) + soilcmas(i,j,k)  ! soil c mass over bare fraction
        end do
      end do

      pCmossmas(i)  = Cmossmas(i)
      do k = 1,ignd
        plitrmsmoss(i) = plitrmsmoss(i) + litrmsmoss(i,k)
      end do 
    end do ! loop 145

    ! Find the canopy covered fraction and the bare fraction of the tiles:
    do i = il1,il2
      fc(i) = sum(fcancmx(i,:))
      fg(i) = 1.0 - fc(i)
    end do

    !     ------------------------------------------------------------------
    !> Initialization ends

    !> Autotrophic respiration
    !!
    !! Leaf respiration is calculated in phtsyn subroutine, while stem
    !! and root maintenance respiration are calculated here. We use air
    !! temperature as a surrogate for  stem temperature
    !!
    !! Find stem and root maintenance respiration in umol co2/m2/sec
    !! If the tracer is being calculated then we also determine the
    !! tracer flux (14C units are amol 14CO2/m2/s).

    call mainres(fcancmx, fc, stemmass, rootmass, & ! In
                 il1, il2, ilg, leapnow, & ! In
                 ta, tbar, rmatctem, sort, isand, & ! In
                 useTracer, tracerStemMass, tracerRootMass, & ! In
                 rmsveg, rmrveg, roottemp, & ! Out
                 rmsTracer, rmrTracer) ! Out

    !> Calculate NPP (net primary productivity), difference between GPP and
    !! autotrophic respriation, for each pft and tile. Also sets the net
    !! photosynthesis to be used by phenology and determines the total
    !! autotrophic respiration fluxes at the PFT and tile-levels. Calculates
    !! values for both upland and peatland sites.
    call calcNPP(il1, il2, ilg, ancgveg, lfstatus, rmlcgveg, & ! In
                 slai, ailcg, sort, rmsveg, rmrveg, mossPresent, & ! In
                 anmoss, fcancmx, rmlmoss, gppmoss, useTracer, & ! In
                 gleafmas, tracerGLeafMass, tracerValue, & ! In
                 rmsTracer, rmrTracer, & ! In/Out
                 bnf_nat, bnf_ant, rootmass_ns, Ncycle_on, bnf_moss, & ! In
                 pheanveg, rmlveg, rml, rms, rmr, rm, rg, npp, gpp, & ! Out
                 autores, autoresveg, rgveg, nppmoss, armoss, & ! Out
                 gppveg, nppmosstep, nppveg, tracerNPP, tracerRML, tracerGPP) ! Out

    !! Find heterotrophic respiration rates (umol co2/m2/sec)
    ! If tracer is being used then calculate the tracer fluxes.
    call heterotrophicRespiration(il1, il2, ilg, peatlandType, mossPresent, fcancmx, & ! In
                                  litrmass, soilcmas, delzw, thpor, tbar, & ! In
                                  psisat, thliq, sort, bi, isand, thice, & ! In
                                  fg, litrmsmoss, upMossSoilC, peatdep, wtable, zbotw, & ! In
                                  useTracer, tracerLitrMass, tracerSoilCMass, tracerMossLitrMass, & ! In
                                  ltresveg, scresveg, litresmoss, socres_peat, socres_moss, & ! Out
                                  resoxic, resanoxic , & ! Out
                                  ltResTracer, sCResTracer, litResMossTracer, soCResPeatTracer) ! Out
    
    !> Find vegetation and tile averaged litter and soil C respiration rates
    !! using values from canopy over ground and canopy over snow subareas.
    !! Also adds the moss and peat soil respiration to the tile level quantities.
    !! Next the litter and soil C pools are updated based on litter and soil C respiration rates.
    !! The humidified litter is then transferred to the soil C pool.
    !! Soil respiration is estimated as the sum of heterotrophic respiration and root maintenance respiration.
    !! For peatlands, we additionally add moss values to the grid (litter respiration
    !! and moss root respiration).

    call updatePoolsHetResp(il1, il2, ilg, fcancmx, ltresveg, scresveg, & ! In
                            peatlandType, mossPresent, fg, litresmoss, socres_peat, socres_moss, & ! In
                            sort, spinfast, rmrveg, rmr, leapnow, & ! In
                            useTracer, ltResTracer, sCResTracer, litResMossTracer, soCResPeatTracer, & ! In
                            litrmass, soilcmas, Cmossmas, litrmsmoss, peatSoilC, upMossSoilC, & ! In / Out
                            tracerLitrMass, tracerSoilCMass, tracerMossCMass, & ! In / Out
                            hetrsveg, litres, socres, hetrores, humtrsvg, soilresp, & ! Out
                            humiftrs, litrfallmoss, ltrestepmoss, &
                            humstepmoss, socrestep) ! Out

    !> Calculate NEP (net ecosystem productivity), difference between NPP and
    !! heterotrophic respiration, for each pft and tile
    call calcNEP(il1, il2, ilg, nppveg, hetrsveg, fg, npp, hetrores, & ! In
                 nep, nepveg) ! Out

    if (doMethane) then 
      
      !> Find CH4 wetland area (if not prescribed) and emissions:
      call  wetlandMethane(hetrores, il1, il2, ilg, & ! In
                           wetfrac, thliq, currlat, sand, &   ! In
                           slopefrac, ta, & ! In
                           ch4WetSpec, wetfdyn, ch4WetDyn) ! Out

      !> Calculate the methane that is oxidized by the soil sink
      call soil_ch4uptake(il1, il2, ilmos, ilg, tbar, & ! In
                        bi, thliq, thice, psisat, & ! In
                        fcanmx, wetfdyn, wetfrac, & ! In
                        isand, ch4conc, thpor, & ! In
                        ch4soills) ! Out
    end if 

    !> Estimate allocation fractions for leaf, stem, and root components.
    call allocate(lfstatus, thliq, ailcg, & ! In
                  il1, il2, ilg, iday, & ! In
                  rmatctem, gleafmas, gleafmas_ns, gleafmas_s, & ! In
                  stemmass, stemmass_ns, stemmass_s, & ! In
                  rootmass, rootmass_ns,rootmass_s, & ! In
                  radj, sort, fcancmx, & ! In
                  isand, THFC, THLW, peatlandType, dayl, dayl_max, & ! In
                  afrleaf, afrstem, afrroot, wtstatus, ltstatus) ! In/Out
 
    call updatePoolsAllocateRepro(il1, il2, ilg, sort, ailcg, lfstatus, nppveg, & ! In
                                  pftexist, gppveg, rmsveg, rmrveg, rmlveg, fcancmx, iday, & ! In
                                  useTracer, tracerNPP, rmsTracer, rmrTracer, tracerRML, tracerGPP, leapnow, & ! In
                                  afrleaf, afrstem, afrroot, gleafmas, gleafmas_ns, & ! In/out
                                  gleafmas_s, bleafmas, stemmass, stemmass_ns, stemmass_s, & ! In /Out
                                  rootmass, rootmass_ns, rootmass_s, lygleafmasmax, & ! In/Out
                                  lystemmassmax, lyrootmassmax, lmaxt, smaxt, rmaxt, leafns2s, & ! In/Out
                                  stemns2s, rootns2s, tracerGLeafMass, tracerStemMass, & ! In/Out
                                  tracerRootMass, tracerBLeafMass, & ! In/Out
                                  reprocost, ntchlveg, ntchsveg, ntchrveg, & ! Out
                                  repro_cost_g, tracerReproCost) ! Out
 
    !! The phenology subroutine determines leaf status for each pft and calculates leaf litter.
    !! the phenology subroutine uses soil temperature (tbar) and root temperature. however,
    !! since CTEM doesn't make the distinction between canopy over ground, and canopy over
    !! snow sub-areas for phenology purposes (for example, leaf onset is not assumed to occur
    !! at different times over these sub-areas) we use average soil and root temperature in
    !! the phenology subroutine.

    call phenology(il1, il2, ilg, leapnow, tbar, thice, &! In
                   thliq, THLW, THFC, ta, peatlandType, &! In
                   pheanveg, iday, radj, roottemp, useTracer, &! In
                   rmatctem, sort, fcancmx, isand, dayl_max, &! In
                   gleafmas, gleafmas_ns, gleafmas_s, & ! In/Out
                   bleafmas, stemmass, stemmass_ns, stemmass_s, tracerStemMass, &! In/Out
                   rootmass, rootmass_ns, rootmass_s, tracerRootMass, &! In/Out
                   lfstatus, pandays, colddays_leaffall, colddays_harvest, lfthrs, & ! In/Out
                   tracerGLeafMass, tracerBLeafMass, & ! In/Out
                   re_alloc_s2l, re_alloc_r2l, re_alloc_sr2l, & !Out
                   gl2bl_grass_cflux, bl2ltr_grass_cflux, & !Out
                   flhrloss, flhrloss_ns, flhrloss_s, tracerFlHrLoss,leaflitr, leaflitr_ns, leaflitr_s, & !In/Out
                   tracerLeafLitr) ! Out
    
    !> While leaf litter is calculated in the phenology subroutine, stem
    !! and root turnover is calculated in the turnoverStemRoot subroutine.
    call turnoverStemRoot(stemmass, stemmass_ns, stemmass_s, & ! In
                          rootmass, rootmass_ns, rootmass_s, & ! In
                          lfstatus, ailcg, & ! In
                          il1, il2, ilg, leapnow, useTracer, &! In
                          sort, fcancmx, tracerStemMass, tracerRootMass, &! In
                          stmhrlos, rothrlos, tracerRothrlos, tracerStmhrlos, & ! In/Out
                          stemlitr, stemlitr_ns, stemlitr_s, & ! In/Out
                          rootlitr, rootlitr_ns, rootlitr_s, & ! In/Out
                          tracerStemLitr, tracerRootLitr) ! Out

    !> Update green leaf biomass for trees and crops, brown leaf biomass for grasses,
    !! stem and root biomass for litter deductions, and update litter pool with leaf
    !! litter calculated in the phenology subroutine and stem and root litter
    !! calculated in the turnoverStemRoot subroutine. Also add the reproduction
    !!  carbon directly to the litter pool. We only add to non-perennially frozen soil
    !! layers so first check which layers are unfrozen and then do the allotment
    !! appropriately. For defining which layers are frozen, we use the active layer depth.
    call updatePoolsTurnover(il1, il2, ilg, reprocost, rmatctem, useTracer, tracerReproCost, & ! In
                             stemmass, stemmass_ns, stemmass_s, & ! In
                             rootmass, rootmass_ns, rootmass_s, & ! In
                             litrmass, rootlitr, rootlitr_ns, rootlitr_s, & ! In/Out
                             gleafmas, gleafmas_ns, gleafmas_s, bleafmas, & ! In/Out
                             leaflitr, leaflitr_ns, leaflitr_s, stemlitr, & ! In/Out
                             stemlitr_ns, stemlitr_s, & ! In/Out
                             tracerGLeafMass, tracerBLeafMass, tracerLitrMass, & ! In/Out
                             tracerRootMass, tracerStemMass, & ! In/Out
                             tracerLeafLitr, tracerStemLitr, tracerRootLitr) ! In/Out 
    
    !> Call the mortality subroutine which calculates mortality due to reduced growth and aging.
    !! Exogenous mortality due to fire and other disturbances and the subsequent litter
    !! that is generated is calculated in the disturb subroutine.
    !!
    !! Set maxage >0 in classicParams.f90 to switch on mortality due to age and
    !! reduced growth. Mortality is linked to the competition parameterization and generates bare fraction.
    call mortalty(stemmass, stemmass_ns, stemmass_s, rootmass, & ! In
                  rootmass_ns, rootmass_s, ailcg, gleafmas, & ! In
                  gleafmas_ns, gleafmas_s, & ! In
                  nstemmass, nstemmass_ns, nstemmass_s, & ! In
                  nrootmass, nrootmass_ns, nrootmass_s, & ! In
                  ngleafmas, ngleafmas_ns, ngleafmas_s, & ! In
                  il1, il2, ilg, & ! In
                  leapnow, iday, sort, fcancmx, peatlandType, & ! In
                  useTracer, tracerStemMass, tracerRootMass, tracerGLeafMass, & ! In
                  lystmmas, lyrotmas, tymaxlai, grwtheff, & ! In/Out
                  stemltrm, stemltrm_ns, stemltrm_s, rootltrm, & ! Out
                  rootltrm_ns, rootltrm_s, glealtrm, glealtrm_ns, & ! Out
                  glealtrm_s, geremort, & ! Out
                  intrmort, tracerStemMort, tracerRootMort, tracerGLeafMort, & ! Out
                  stemltrm_n, stemltrm_ns_n, stemltrm_s_n, & ! Out
                  rootltrm_n, rootltrm_ns_n, rootltrm_s_n, & ! Out
                  glealtrm_n, glealtrm_ns_n, glealtrm_s_n) ! Out

    !> Update leaf, stem, and root biomass pools to take into account losses due to mortality, and put the
    !! litter into the litter pool. Mortality of green grasses doesn't generate litter, instead they turn brown.
    if(.not. PFTCompetition) then
      call updatePoolsMortality(il1, il2, ilg, stemltrm, stemltrm_ns, stemltrm_s, & ! In
                              rootltrm, rootltrm_ns, rootltrm_s, & ! In
                              stemltrm_n, stemltrm_ns_n, & ! In
                              stemltrm_s_n, rootltrm_n, rootltrm_ns_n, rootltrm_s_n, & ! In
                              useTracer, & ! In
                              rmatctem, tracerStemMort, tracerRootMort, tracerGLeafMort, & ! In
                              stemmass, stemmass_ns, stemmass_s, rootmass, rootmass_ns, & ! In/Out
                              rootmass_s, litrmass, glealtrm, glealtrm_ns, glealtrm_s, & ! In/Out
                              gleafmas, gleafmas_ns, gleafmas_s, bleafmas, & ! In/Out
                              nstemmass, nstemmass_ns, nstemmass_s, & ! In/Out
                              nrootmass, nrootmass_ns, nrootmass_s, nlitrmass, & ! In/Out
                              glealtrm_n, glealtrm_ns_n, glealtrm_s_n, & ! In/Out
                              ngleafmas, ngleafmas_ns, ngleafmas_s, nbleafmas, & ! In/Out
                              tracerLitrMass, & ! In/Out
                              tracerStemMass, tracerRootMass, tracerGLeafMass, tracerBLeafMass) ! In/Out

    else
      do j = 1,icc
        do i = il1,il2
          stemltrm(i,j) = 0.0
          stemltrm_ns(i,j) = 0.0
          stemltrm_s(i,j) = 0.0
          rootltrm(i,j) = 0.0
          rootltrm_ns(i,j) = 0.0
          rootltrm_s(i,j) = 0.0
          glealtrm(i,j) = 0.0
          glealtrm_ns(i,j) = 0.0
          glealtrm_s(i,j) = 0.0
          stemltrm_n(i,j) = 0.0
          stemltrm_ns_n(i,j) = 0.0
          stemltrm_s_n(i,j) = 0.0
          rootltrm_n(i,j) = 0.0
          rootltrm_ns_n(i,j) = 0.0
          rootltrm_s_n(i,j) = 0.0
          glealtrm_n(i,j) = 0.0
          glealtrm_ns_n(i,j) = 0.0
          glealtrm_s_n(i,j) = 0.0
        end do
      end do
    end if
    !> Call the disturbance subroutine, which calculates mortality due to fire and other disturbances.
    !! The primary outputs from the disturbance subroutine are litter generated, C emissions due to fire
    !! and area burned, which may be used to estimate change in fractional coverages.
    !!
    !! Disturbance is spatial and requires the area of the gcm grid cell and areas of different pfts present in
    !! a given grid cell. However, when ctem is operated at a point scale then it is assumed that the
    !! spatial scale is 1 hectare = 10,000 m2. the disturbance subroutine may be stopped from simulating
    !! any fire by specifying fire extinguishing probability equal to 1.
    call disturbance (thliq, THLW, THFC, THPOR, uwind, useTracer, peatlandType, & ! In
                      vwind, lightng, elcggat, fcancmx, isand, & ! In
                      rmatctem, ilg, il1, il2, sort, & ! In
                      grclarea, thice, popdin, lucemcom, lucemcomn, & ! In
                      dofire, fireSubmodule, currlat, iday, fsnow, prsfireareagat, prescribedFire, & ! In
                      stemmass, stemmass_ns, stemmass_s, & ! In/Out
                      rootmass, rootmass_ns, rootmass_s, & ! In/Out
                      gleafmas, gleafmas_ns, gleafmas_s, & ! In/Out
                      bleafmas, litrmass, nstemmass, nstemmass_ns, & ! In/Out
                      nstemmass_s, nrootmass, nrootmass_ns, nrootmass_s, & ! In/Out
                      ngleafmas, ngleafmas_ns, ngleafmas_s, nbleafmas, & ! In/Out
                      nlitrmass, & ! In/Out
                      tracerStemMass, tracerRootMass, tracerGLeafMass, tracerBLeafMass, tracerLitrMass, & ! In/Ou
                      stemltdt, rootltdt, glfltrdt, blfltrdt, & ! Out (Primary)
                      glcaemls, rtcaemls, stcaemls, & ! Out (Primary)
                      blcaemls, ltrcemls, burnfrac, & ! Out (Primary)
                      fFireCveg, fFireLitter, fFireCsoil, & ! Out (Primary)
                      nstemltdt, nrootltdt, nglfltrdt, nblfltrdt, & ! Out (Primary)
                      glniemls, rtniemls, stniemls, & ! Out (Primary)
                      blniemls, ltrnemls, & ! Out (Primary)
                      pstemmass, pgleafmass, emit_co2, emit_ch4, & ! Out (Primary)
                      emit_co, emit_nmhc, emit_h2, emit_nox, & ! Out (Secondary)
                      emit_n2o, emit_nh3, emit_pm25, emit_tpm, emit_tc, & ! Out (Secondary)
                      emit_bc, emit_oc, emit_so2, burnvegf, bterm_veg, & ! Out (Secondary)
                      fiengat, fabgat, nmfigat, & ! Out
                      mterm_veg, lterm, smfunc_veg, & ! Out (Secondary)
                      trackTileAge, dynamicTilingOn, tileAgegat, & ! In/Out (logical flags and tile age)
                      drgtstrs, betadrgt, precip, ta, QACGAT, PRESGAT, &
                      isi, bui, fwi, ffmc, dmc, dc, IMONTH, NoonTempc, & 
                      NoonPrecip, NoonWindtot, NoonRhum, snwct) !FWI vars


    !> Calculate NBP (net biome production) for each pft by taking into account
    !! C emission losses. The disturbance routine produces emissions due to fire
    !! while the land use change subroutine calculates emissions due to LUC.
    !! The LUC related combustion flux is assumed to be spread uniformly over the
    !! tile as it is no longer associated with any one PFT. To calculate NBP,
    !! we do not subtract LUC emissions from the PFT-level NBP, but we do subtract
    !! it from the per tile NBP.
    call calcNBP(il1, il2, ilg, nepveg, fcancmx, & ! In
                 lucemcom, ltresveg, scresveg, nep, & ! In
                 glcaemls, blcaemls, stcaemls, rtcaemls, ltrcemls, & ! In/Out
                 nbpveg, dstcemls1, dstcemls3, nbp, dstcemls2,  nepCMIP) ! Out

    if (turbationON) then
      !> Allow cryoturbation and bioturbation to move the soil C between
      !! layers. Since this is neither consuming nor adding C, this does not
      !! affect our C balance in balcar. There is also an internal C balance check.
      call turbation(il1, il2, zbotw, isand, &
                    maxAnnualActLyr, spinfast, peatlandType, &! In
                    litrmass, soilcmas)! In/Out
    end if
    
    if (useTracer /= "None") call turbation(il1, il2, zbotw, isand, &
                                      maxAnnualActLyr, spinfast, peatlandType, &! In
                                      tracerLitrMass, tracerSoilCMass) ! In/Out

    ! if (any(mossPresent /= 'None')) call turbation(il1, il2, zbotw, isand, &
    !                                     maxAnnualActLyr, spinfast, peatlandType, &! In
    !                                     litrmsmoss, upMossSoilC) ! In/Out

    !> If there is a moss coverage, update the moss carbon pool for npp and litterfall.
    if (any((mossPresent) /= 'None')) call updateMossC(Cmossmas,nppmosstep,litrfallmoss)   

    !> Prepare for the carbon balance check. Calculate total litter fall from each
    !! component (leaves, stem, and root) from all causes (normal turnover, drought
    !! and cold stress for leaves, mortality, and disturbance), calculate grid-average
    !! vegetation biomass, litter mass, and soil carbon mass, and litter fall rate.
    !! Also add the bare ground values to the grid-average. If a peatland or moss 
    !! covered, we assume no bareground and add the moss values instead. 
    !! Note: peatland soil C is not aggregated from plants but updated
    !! by humification and respiration from the previous stored value
    call prepBalanceC(il1, il2, ilg, fcancmx, glealtrm, glfltrdt, &  ! In
                      blfltrdt, stemltrm, stemltdt, rootltrm, rootltdt, & ! In
                      peatlandType, mossPresent, pgavscms, peatSoilC, & ! In
                      stemlitr, rootlitr, rootmass, & ! In
                      litrmass, soilCmas, stemmass, bleafmas, & ! In
                      gleafmas, socrestep, fg, litrfallmoss, & ! In
                      Cmossmas, litrmsmoss, upMossSoilC, & ! In
                      leaflitr, & ! In/Out
                      tltrleaf, tltrstem, tltrroot, & ! Out
                      vgbiomas, litrfall, gavgltms, litrfallveg, &! Out
                      gavgscms, vgbiomas_veg) ! Out


    !! Prepare additional CMIP6 output variables related to land use change
    do i = il1,il2 ! loop 198
       cproduct(i) = litrmass(i,iccp2,1) + soilcmas(i,iccp2,1) !< LUC products in paper and furniture pools
       fproductdecomp(i) = ltresveg(i,iccp2,1) + scresveg(i,iccp2,1) !< Decomposition of LUC product pools
    enddo ! loop 198

    !! Prepare additional CMIP6 output variables, cleaf, cstem, croot, nppleaf, nppstem, npproot 
    do j = 1, icc ! loop 199
      do i = il1, il2
        !
        cleaf(i)=cleaf(i)+fcancmx(i,j)*gleafmas(i,j)
        cstem(i)=cstem(i)+fcancmx(i,j)*stemmass(i,j)
        croot(i)=croot(i)+fcancmx(i,j)*rootmass(i,j)
        !
        nppleaf(i) = nppleaf(i) + fcancmx(i,j)*ntchlveg(i,j)
        nppstem(i) = nppstem(i) + fcancmx(i,j)*ntchsveg(i,j)
        npproot(i) = npproot(i) + fcancmx(i,j)*ntchrveg(i,j)
      end do 
    end do ! loop 199

#if !defined without_agcm_
    !----------------------------------------------------- \
    !
    ! kludge 2.0
    !
    ! When running coupled to CanESM (32bit), 
    ! adjust gpp and autores (forced C conservation) and other fluxes to
    ! make them consistent with forced C conservation changes
    ! adjusting only GPP and autores because that's where the non C conservation
    ! bug appears to be. We are pretty sure about this.

    do i = il1,il2 ! loop 200

      ! get NBP from pools

      pool_based_nbp(i) = ( pvgbioms_before_luc(i) + pgavltms_before_luc(i) + pgavscms_before_luc(i) + pcproduct_before_luc(i) ) - &
                           ( vgbiomas(i) + gavgltms(i) + gavgscms(i) + cproduct(i) )

      ! convert pool_based_nbp from Kg C/m2.day to u-mol CO2-C/m2.s

      pool_based_nbp(i)=pool_based_nbp(i)*(963.62/deltat)*-1.0 ! (I should have subtracted old from new pools)

      ! now work backward from this pool_based_nbp to calculate NPP
      ! been subtracted from vgbiomas when it came of the LUC subroutine, so the pool_based_nbp
      ! at this stage doesn't include lucemcom(i) so we do not subtract it

      pool_based_npp(i) = pool_based_nbp(i)  &
                          + hetrores(i) + fProductDecomp(i) + dstcemls2(i) & 
                          + lucemcom(i)

      !! we are assuming based on observing the nature of non conservation that the bug
      !! lies in GPP and Ra.

      !! Find the difference between pool-based and flux based npp

      diff1 = pool_based_npp(i) - npp(i)

      !! this diff1 needs to proportionally adjusted between GPP and Ra

      if ( (gpp(i)+autores(i)).gt.0.0 ) then
        t1=gpp(i)/(gpp(i)+autores(i))
        t2=autores(i)/(gpp(i)+autores(i))
      else
        t1=0.0
        t2=0.0
      endif

      t4=diff1*t1
      gpp(i)=gpp(i)+t4
      t4=diff1*t2
      autores(i)=autores(i)-t4

      ! also adjust nppleaf, nppstem, and npproot
      ! this is more tricky. If all the 3 components of NPP are +ve or -ve it's easy to 
      ! scale them but if some of them are +ve and some -ve then it's very tricky.
      ! read comments for subroutine scale_numbers
      ! Regardless, this scaling seems to work 99.9993 % of the time.

      call scale_numbers(nppleaf(i), npproot(i), nppstem(i), npp(i), pool_based_npp(i), & 
                         scaled_nppleaf, scaled_npproot, scaled_nppstem)
      
      nppleaf(i)=scaled_nppleaf
      npproot(i)=scaled_npproot
      nppstem(i)=scaled_nppstem


      !! now find the grid averaged fluxes back again which should, in principle, be now
      !! consistent with pools (well cLand to be exact). 

      npp(i)=gpp(i)-autores(i) 
      nep(i)=npp(i)-hetrores(i) 
      nbp(i)=pool_based_nbp(i)

      if ( (rg(i)+rm(i)).gt.0 ) then
        t1=rg(i)/(rg(i)+rm(i))
        t2=rm(i)/(rg(i)+rm(i))
        rg(i)=t1*autores(i) 
        rm(i)=t2*autores(i) 
      endif
      if ( rm(i).gt.0 ) then
        t1=rml(i)/(rml(i)+rms(i)+rmr(i))
        t2=rms(i)/(rml(i)+rms(i)+rmr(i))
        t3=rmr(i)/(rml(i)+rms(i)+rmr(i))
        rml(i)=t1*rm(i)
        rms(i)=t2*rm(i)
        rmr(i)=t3*rm(i)
      endif


    enddo  ! loop 200   
    !----------------------------------------------------- /
#endif
 
    !> At this stage we have all required fluxes in u-mol co2/m2.sec and initial (loop 140 and 145)
    !! and updated sizes of all pools (in \f$(kg C/m^2)\f$). Now we call the balcar subroutine and make sure
    !! that C in leaves, stem, root, litter and soil C pool balances within a certain tolerance.
    if (spinfast == 1) then
      call  balcar(gleafmas, gleafmas_ns, gleafmas_s, stemmass, &
                   stemmass_ns, stemmass_s, rootmass, rootmass_ns, &
                   rootmass_s, bleafmas, &
                   litrmass, soilcmas, ntchlveg, ntchsveg, &
                   ntchrveg, tltrleaf, tltrstem, tltrroot, &
                   glcaemls, blcaemls, stcaemls, rtcaemls, &
                   ltrcemls, ltresveg, scresveg, humtrsvg, &
                   pglfmass, pblfmass, pstemass, protmass, &
                   plitmass, psocmass, vgbiomas, reprocost, &
                   pvgbioms, gavgltms, pgavltms, gavgscms, &
                   pgavscms, dstcemls3, repro_cost_g, &
                   autores, hetrores, gpp, &
                   litres, socres, dstcemls1, &
                   litrfall, humiftrs, &
                   il1, il2, ilg, &
                   peatlandType, mossPresent, Cmossmas, pCmossmas, &
                   nppmosstep, litrfallmoss, litrmsmoss, &
                   plitrmsmoss, ltrestepmoss, humstepmoss, &
                   re_alloc_s2l, re_alloc_r2l, re_alloc_sr2l)

      !> Check for mass balance for the tracer if the simple tracer is being used and
      !! doTracerBalance is true. If you are going to use the tracer, it is recommended
      !! that you perform this check prior to any simulations with other tracers (14C, 13C)

      if (useTracer == 'Simple' .and. doTracerBalance) then
        call checkTracerBalance(il1, il2, useTracer, &
                                stemmass_s, stemmass_ns, rootmass_s, rootmass_ns, litrmass, &
                                gleafmas_s, gleafmas_ns, bleafmas, soilcmas, &
                                tracerGLeafMass, tracerStemMass, tracerRootMass, &
                                tracerLitrMass, tracerBLeafMass, tracerSoilCMass)
      end if

    end if

    !> Finally find vegetation structural attributes which can be passed
    !! to the land surface scheme using leaf, stem, and root biomass.
    !>
    call allometry(gleafmas, gleafmas_ns, gleafmas_s, bleafmas, & ! In
                   stemmass, stemmass_ns, stemmass_s, rootmass, & ! In
                   rootmass_ns, rootmass_s, il1, il2, ilg, zbotw, & ! In
                   soildpth, fcancmx, & ! In
                   peatlandType, maxAnnualActLyr, & ! In
                   ailcg, ailcb, ailc, zolnc, & ! Out
                   rmatc, rmatctem, slai, bmasveg, & ! Out
                   cmasvegc, veghght, rootdpth, alvisc, & ! Out
                   alnirc, paicgat, slaicgat) ! Out

    !> Calculation of gavglai is moved from loop 1100 to here since ailcg is updated by allometry
    gavglai (:) = 0.0
    do j = 1,icc
      do i = il1,il2
        gavglai(i) = gavglai(i) + fcancmx(i,j) * ailcg(i,j)
      end do
    end do

    !> At the end of the day, find the depth of the peat (can be dynamic or static), 
    do i = il1,il2
      if (peatlandType(i) /= 'None') then     
        if (.not. useStaticPeatDep) then
          ! Dynamically determine the peat depth based upon the amount of peat soil C
          peatdep(i) = peatDepth(peatSoilC(i)) ! only bogs and fen
        else 

          ! First remove the peatSoilC value (since it has been changed during the timestep) from
          ! the grid average soil C mass (used for carbon balance):
          gavgscms(i) = gavgscms(i) - peatSoilC(i)

          ! Don't allow the peat depth to dynamically change during the run. This is used for the purpose
          ! of spinup to ensure that the disequilibrium between the peatland and the driving climate, time period 
          ! since the peatland initiation, or sub-grid heterogeneity impacts don't adversely impact the ability 
          ! of the model to spinup. 
          peatdep(i) = soildpth(i)

          ! Also ensure the peatSoilC remains based upon this peat depth.
          peatSoilC(i) = peatStorage(peatdep(i))

          ! Now redefine the grid average soil C mass using the defined peatSoilC. The changed value is
          ! used as the previous value for the next day so must be corrected to the constant peatSoilC value.
          gavgscms(i) = gavgscms(i) + peatSoilC(i)
          
        end if 
      end if 
    end do

    !> and the peat bottom layer depth
    if (any(peatlandType /= 'None')) call peatDayEnd(il1, il2, delzw, peatdep, soildpth, zbot)

    !> Now that the model has passed all the caculations, processes, fluxes, and pools associated with the
    !! Carbon cycle, the "n_processes" subroutine is called to start off the entire Nitrogen cycle.
    if (Ncycle_on) then

      call n_processes(il1, il2, spinfast, thliq, thice, thpor, zbotw, tbar, isand, & ! In
                       iday, radj, lfstatus, THFC, THLW, fcancmx, & ! In
                       delzw, sort, peatlandType, ROFB, CFLUX_GA, USTARBS_GA, leapnow, & ! In
                       QFC, rootdpth, soilpH, nfertil, ndeposit, ntchlveg, ntchsveg, & ! In
                       ntchrveg, gleafmas, gleafmas_ns, gleafmas_s, bleafmas, stemmass, & ! In
                       stemmass_ns, stemmass_s, rootmass, rootmass_ns, rootmass_s, & ! In
                       litrmass, soilcmas, leafns2s, stemns2s, rootns2s, & ! In
                       tltrleaf, tltrstem, tltrroot, gl2bl_grass_cflux, bl2ltr_grass_cflux, & ! In
                       re_alloc_sr2l, ailcg, lfthrs, co2conc, & ! In
                       humtrsvg, ltresveg, scresveg, psisat, bi, & ! In
                       humstepmoss, ltrestepmoss, socres_moss, & ! In
                       mossPresent, Cmossmas, litrmsmoss, & ! In
                       upMossSoilC, litrfallmoss, nppmosstep, & ! In
                       ngleafmas, ngleafmas_ns, ngleafmas_s, nbleafmas, & ! In/Out
                       nstemmass, nstemmass_ns, nstemmass_s, & ! In/Out
                       nrootmass, nrootmass_ns, nrootmass_s, & ! In/Out
                       nlitrmass, soilnmas, Nmossmas, & ! In/Out
                       litrmsmossN, upMossSoilN, &! In/ Out
                       bnf_free, bnf_nat, bnf_ant, bnf_tot, nstress, & !Out
                       nh4_mass, no3_mass, nvolveg, nleachveg, & ! Out
                       nitrifveg, no_nitveg, no_denitveg, & ! Out
                       no_nitdenitveg, n2o_nitveg, n2o_denitveg, & ! Out
                       n2o_nitdenitveg, n2_denitveg, appl_fert, & ! Out
                       ndep_nh4, ndep_no3, ndemandveg_wp_npp, & ! Out
                       nuptakeveg_p_nh4, nuptakeveg_p_no3, & ! Out
                       nuptakeveg_a_actl_nh4, & ! Out
                       nuptakeveg_a_actl_no3, nuptakeveg, nallocveg_l, nallocveg_s, & ! Out
                       nallocveg_r, nresorpedveg_s, & ! Out
                       nresorpedveg_r, nre_allocveg_s2l, & ! Out
                       nre_allocveg_r2l, nleafns2sveg, nstemns2sveg, & ! Out
                       nrootns2sveg, nlitrveg_l, nlitrveg_s, & ! Out
                       nlitrveg_r, nlitrveg, gl2bl_grass_nflux, & ! Out
                       c2nveg_l, c2nveg_s, c2nveg_r, & ! Out
                       c2nveg_wp, c2nveg_litr, c2nveg_humus, & ! Out
                       nhumtrsveg, nmineralveg_litr, & ! Out
                       nmineralveg_humus, netnmineralveg, nimmobilveg_nh4, & ! Out
                       nimmobilveg_no3, nvgbiomas_veg, & ! Out
                       fNnetlandveg, redcoeff_vcmax, redcoeff_vcmaxMoss, bnf_moss) ! Out

      ! Calculate grid average values for CanESM
      do i = il1,il2
        do j = 1,icc
          nMineralNH4(i) = nMineralNH4(i) + nh4_mass(i,j) * fcancmx(i,j) * convertg2kg
          nMineralNO3(i) = nMineralNO3(i) + no3_mass(i,j) * fcancmx(i,j) * convertg2kg
          nVeg(i) = nVeg(i) + nvgbiomas_veg(i,j) * fcancmx(i,j) * convertg2kg
          nLitter(i) = nLitter(i) + nlitrmass(i,j) * fcancmx(i,j) * convertg2kg
          nSoil(i) = nSoil(i) + soilnmas(i,j) * fcancmx(i,j) * convertg2kg
          nLeaf(i) = nLeaf(i) + ngleafmas(i,j) * fcancmx(i,j) * convertg2kg
          nStem(i) = nStem(i) + nstemmass(i,j) * fcancmx(i,j) * convertg2kg
          nRoot(i) = nRoot(i) + nrootmass(i,j) * fcancmx(i,j) * convertg2kg
          fBNF(i) = fBNF(i) + bnf_tot(i,j) * fcancmx(i,j) * convertkgN  
          fNdep(i) = fNdep(i) + (ndep_nh4(i,j) + ndep_no3(i,j)) * fcancmx(i,j) * convertkgN
          fNfert(i) = fNfert(i) + appl_fert(i,j) * fcancmx(i,j) * convertkgN
          fNgasNonFire(i) = fNgasNonFire(i) + (nvolveg(i,j)+n2o_nitdenitveg(i,j)+n2_denitveg(i,j)+no_denitveg(i,j)) * fcancmx(i,j) * convertkgN
          fNgasFire(i) = fNgasFire(i) + (emit_n2o(i,j)*convertn2o2N) + (emit_nh3(i,j)*convertnh32N) + (emit_nox(i,j)*convertnox2N)  * fcancmx(i,j) * convertkgN
          !fNgas(i) only done below in iccp1 section
          fNnetmin(i) = fNnetmin(i) + netnmineralveg(i,j) * fcancmx(i,j) * convertkgN
          fNOx(i) = fNOx(i) + ((emit_nox(i,j)*convertnox2N) + n2o_nitdenitveg(i,j)) * fcancmx(i,j) * convertkgN
          fNup(i) = fNup(i) + nuptakeveg(i,j) * fcancmx(i,j) * convertkgN
          fNvegSoil(i) = 0.0 !this is zero in the model at present (i.e. no fluxes from veg to soil)
          fN2o(i) = fN2o(i) + (n2o_nitdenitveg(i,j) + (emit_n2o(i,j)*convertn2o2N)) * fcancmx(i,j) * convertkgN
          fNleach(i) = fNleach(i) + nleachveg(i,j) * fcancmx(i,j) * convertkgN
          fNVol(i) = fNVol(i) + nvolveg(i,j) * fcancmx(i,j) * convertkgN

        end do
          nMineralNH4(i) = nMineralNH4(i) + nh4_mass(i,iccp1) * barefrac(i) * convertg2kg
          nMineralNO3(i) = nMineralNO3(i) + no3_mass(i,iccp1) * barefrac(i) * convertg2kg
          nLitter(i) = nLitter(i) + nlitrmass(i,iccp1) * barefrac(i) * convertg2kg
          nSoil(i) = nSoil(i) + soilnmas(i,iccp1) * barefrac(i) * convertg2kg
          fBNF(i) = fBNF(i) + bnf_tot(i,iccp1) * barefrac(i) * convertkgN  
          fNdep(i) = fNdep(i) + (ndep_nh4(i,iccp1) + ndep_no3(i,iccp1)) * barefrac(i) * convertkgN
          fNfert(i) = fNfert(i) + appl_fert(i,iccp1) * barefrac(i) * convertkgN
          fNgasNonFire(i) = fNgasNonFire(i) + (nvolveg(i,iccp1)+n2o_nitdenitveg(i,iccp1)+n2_denitveg(i,iccp1)+no_denitveg(i,iccp1)) * barefrac(i) * convertkgN
          !fNgasFire(i) calculation omitted because emit_* variables don't include bare ground
          fNgas(i) = fNgasNonFire(i) + fNgasFire(i) !this is only done here after bare ground is considered
          fNnetmin(i) = fNnetmin(i) + netnmineralveg(i,iccp1) * barefrac(i) * convertkgN
          fNOx(i) = fNOx(i) + n2o_nitdenitveg(i,iccp1) * barefrac(i) * convertkgN !emit_nox omitted because emit_* variables don't include bare ground
          fN2o(i) = fN2o(i) + n2o_nitdenitveg(i,iccp1) * barefrac(i) * convertkgN !emit_n2o omitted because emit_* variables don't include bare ground
          fNleach(i) = fNleach(i) + nleachveg(i,iccp1) * barefrac(i) * convertkgN
          fNVol(i) = fNVol(i) + nvolveg(i,iccp1) * barefrac(i) * convertkgN

      enddo

    end if ! N cycle

    return

    end associate
    return
  end subroutine ctem
  
  ! ----------------------------------------------------------------------------
  !> \ingroup ctemdriver_calcnpp
  
  !> Calculate NPP (net primary productivity), difference between GPP and
  !! autotrophic respriation, for each pft and tile. Also sets the net
  !! photosynthesis to be used by phenology and determines the total
  !! autotrophic respiration fluxes at the PFT and tile-levels. Calculates
  !! values for both upland and peatland sites.
  !! @author V. Arora, J. Melton
  subroutine calcNPP (il1, il2, ilg, ancgveg, lfstatus, rmlcgveg, & ! In
                      slai, ailcg, sort, rmsveg, rmrveg, mossPresent, & ! In
                      anmoss, fcancmx, rmlmoss, gppmoss, useTracer, & ! In
                      gleafmas, tracerGLeafMass, tracerValue, & ! In
                      rmsTracer, rmrTracer, & ! In/Out
                      bnf_nat, bnf_ant, rootmass_ns, Ncycle_on, bnf_moss, & ! In
                      pheanveg, rmlveg, rml, rms, rmr, rm, rg, npp, gpp, & ! Out
                      autores, autoresveg, rgveg, nppmoss, armoss, & ! Out
                      gppveg, nppmosstep, nppveg, tracerNPP, tracerRML, tracerGPP) ! Out

    use classicParams, only : icc, iccp1, kn, zero, grescoefmoss, deltat, &
                              grescoef, C_cost_bnf_s, C_cost_mossBNF
    use autotrophicRespiration, only : growthRespiration

    implicit none

    ! arguments
    logical, intent(in) :: Ncycle_on
    integer, intent(in) :: il1             !< il1=1
    integer, intent(in) :: il2             !< il2=ilg (no. of grid cells in latitude circle)
    integer, intent(in) :: ilg
    character(6), intent(in) :: useTracer !< character: Switch for use of a model tracer. If useTracer is 'None' then the
                                          !! tracer code is not used. useTracer = 'Simple' turns on a simple tracer that tracks
                                          !! pools and fluxes. The simple tracer then requires that the tracer values in
                                          !! the init_file and the tracerCO2file are set to meaningful values for the experiment being run.
                                          !! useTracer = '14C' means the tracer is 14C and will then call a 14C decay scheme.
                                          !! useTracer = '13C' means the tracer is 13C and will then call a 13C fractionation scheme.
    real, intent(in) :: ancgveg(:,:)      !< net photosynthetic rate for CTEM's pfts
    integer, intent(in) :: lfstatus(:,:)  !< leaf phenology status
    real, intent(in) :: rmlcgveg(:,:)     !< leaf respiration rate for CTEM's pfts
    real, intent(in) :: slai(:,:)         !< storage/imaginary lai for phenology purposes
    real, intent(in) :: ailcg(:,:)        !< green lai for ctem's 9 pfts
    integer, intent(in) :: sort(:)        !< index for correspondence between biogeochem pfts and the number of values in parameters vectors in run params file.
    real, intent(inout) :: rmrveg(:,:)       !< Maintenance respiration for leaf for the CTEM pfts in u mol co2/m2. sec
    real, intent(in) :: rmsveg(:,:)       !< Maintenance respiration for stem for the CTEM pfts in u mol co2/m2. sec
    character(8), intent(in) :: mossPresent(ilg)   !< Type of moss present (else 'None')
    real, intent(in) :: anmoss(:)         !< moss net photoysnthesis -daily averaged C fluxes rates (umol/m2/s)
    real, intent(in) :: fcancmx(:,:)      !< max. fractional coverage of CTEM's pfts, but this can be
                                          !! modified by land-use change,and competition between pfts
    real, intent(inout) :: rmlmoss(:)        !< moss maintainance respiration -daily averaged C fluxes rates (umol/m2/s)
    real, intent(in) :: gppmoss(:)        !< moss GPP -daily averaged C fluxes rates (umol/m2/s)
    real, intent(in) :: gleafmas(:,:)     !< green leaf mass for each of the ctem pfts, \f$(kg C/m^2)\f$
    real, intent(in) :: tracerGLeafMass(:,:) !< Tracer mass in the green leaf pool, 14C: \f$ng ^{14}C/m^2\f$
    real, intent(in) :: tracerValue(:)       !< Tracer CO2 value read in from tracerCO2File, units vary (14C: 1E12* 14C/C ratio)
    real, intent(in) :: bnf_nat(ilg,icc) !< natural symbiotic biological nitrogen fixation for PFTs + bare soil, \f$g N m^{-2} day^{-1}\f$
    real, intent(in) :: bnf_ant(ilg,icc) !< anthropogenic symbiotic biological nitrogen fixation for PFTs + bare soil, \f$g N m^{-2} day^{-1}\f$
    real, intent(in) :: rootmass_ns(ilg,icc)  !< non-structural root mass for each of the ctem pfts, \f$(kg C/m^2)\f$
    real, intent(in) :: bnf_moss(ilg)    !< moss-associated biological nitrogen fixation, \f$g N m^{-2} day^{-1}\f$

    real, intent(inout) :: rmsTracer(:,:)   !< Tracer maintenance respiration for stem for the CTEM pfts (14C: \f$amol ^{14}CO_2 m^{-2} s^{-1}\f$)
    real, intent(inout) :: rmrTracer(:,:)   !< Tracer maintenance respiration for root for the CTEM pfts both (14C: \f$amol ^{14}CO_2 m^{-2} s^{-1}\f$)

    real, intent(out) :: pheanveg(ilg,icc) !< net photosynthesis passed to phenology (adopts value of ancgveg)
    real, intent(out) :: rmlveg(ilg,icc)   !< Leaf maintenance respiration per PFT (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: gppveg(ilg,icc)   !< Gross primary productivity per PFT (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: rml(ilg)          !< Tile level leaf maintenance respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: rms(ilg)         !< Tile level stem maintenance respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: rmr(ilg)         !< Tile level Root maintenance respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: rm(ilg)          !< Tile level maintenance respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: rg(ilg)          !< Tile level growth respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: npp(ilg)         !< Tile-level net primary productivity (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: gpp(ilg)         !< Tile level gross primary productivity (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: autores(ilg)     !< Tile level autotrophic respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: autoresveg(ilg,icc)  !< PFT level autotrophic respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: rgveg(ilg,icc)   !< PFT level growth respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: nppmoss(ilg)     !< Net primary production of moss (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: armoss(ilg)      !< autotrophic respiration of moss (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: nppmosstep(ilg)  !< moss NPP (kgC/m2/timestep)
    real, intent(out) :: nppveg(ilg,icc)  !< NPP for individual pfts, (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: tracerNPP(ilg,icc) !< tracer NPP for individual pfts, (14C: \f$amol ^{14}CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: tracerRML(ilg,icc) !< Tracer leaf maintenance respiration (14C: \f$amol ^{14}CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: tracerGPP(ilg,icc) !< Tracer GPP (14C: \f$amol ^{14}CO_2 m^{-2} s^{-1}\f$)

    ! Local
    integer :: j,i
    real :: anveg(ilg,icc)       !<
    real :: term                 !<
    real :: rgmoss(ilg)          !< moss growth respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real :: rmveg(ilg,icc)       !<
    real :: frac                 !< temp var.
    real :: tracerRG(ilg,icc)    !< Tracer growth respiration (14C: \f$amol ^{14}CO_2 m^{-2} s^{-1}\f$)
    real :: tracerRM(ilg,icc)    !<  Tracer total maintenance respiration (14C: \f$amol ^{14}CO_2 m^{-2} s^{-1}\f$)
    real :: extratracer          !< Temp var used in tracer calcs. (14C: \f$amol ^{14}CO_2 m^{-2} s^{-1}\f$)
    real :: tracerRespRootStem   !< Temp var used in tracer calcs. maintenance resp for root and stems (14C: \f$amol ^{14}CO_2 m^{-2} s^{-1}\f$)

    !--------

    ! NOTE: This next bit is a little tricky. Remember ancgveg is the net photosynthesis
    ! so it is ancgveg = gpp - rmlcgveg. Here we assign the daily mean net photosynthesis
    ! (ancgveg) and mean rml (rmlveg) to the anveg and rmlveg variables
    ! and their sum to gpp veg since gpp = anveg + rml if we both have some coverage of the
    ! PFT (fcancmx > 1) and the leaves are not imaginary (lfstatus /= 4).
    ! If no real leaves are in existence we leave rml, anveg and gpp set to 0.
    ! However, if we have real leaves, but they are quite small we are going to
    ! reduce rml, but we need to use the original rml to find the gppveg otherwise
    ! the gpp will not be corrected properly. We store the original ancgveg for
    ! phenology to test if the leaves should be coming out.

    pheanveg = 0.
    rmlveg = 0.
    anveg = 0.
    gppveg = 0.
    nppveg = 0. 
    tracerRM = 0.
    do j = 1,icc
      do i = il1,il2

        if (fcancmx(i,j) > zero) then

          pheanveg(i,j) = ancgveg(i,j) ! to be used for phenology purposes

          if (lfstatus(i,j) /= 4) then ! real leaves so use values

            anveg(i,j) = ancgveg(i,j)
            rmlveg(i,j) = rmlcgveg(i,j)
            gppveg(i,j) = anveg(i,j) + rmlveg(i,j)
            
            if (slai(i,j) > ailcg(i,j)) then
              term = ((1.0 / kn(sort(j))) * (1.0 - exp( - kn(sort(j)) * ailcg(i,j))) &
                     / (1.0 / kn(sort(j))) * (1.0 - exp( - kn(sort(j)) * slai(i,j))))
              rmlveg(i,j) = rmlveg(i,j) * term
            end if
            ! else
            ! the leaves were imaginary so leave variables set to the initialized 0 value.
          end if
        end if
      end do ! loop 190
    end do ! loop 180

    !! Find total maintenance respiration values and net primary productivity.
    rml(:) = 0.
    rms(:) = 0.
    rmr(:) = 0.
    rm(:) = 0.
    rg(:) = 0.
    npp(:) = 0.
    gpp(:) = 0.
    tracerNPP(:,:) = 0.

    do j = 1,icc
      do i = il1,il2

        ! Remove C cost of BNF:
        if(rootmass_ns(i,j) > 0.0 .and. Ncycle_on) then
          rmrveg(i,j) = rmrveg(i,j) + (bnf_nat(i,j) + bnf_ant(i,j))*(C_cost_bnf_s/1000.*963.62)
        end if

        rmveg(i,j)  = rmlveg(i,j) + rmrveg(i,j) + rmsveg(i,j)
        nppveg(i,j) = gppveg(i,j) - rmveg(i,j)

        if (useTracer /= 'None') then
          !> Determine the NPP for the tracer. To find the value of rml for the tracer
          !! we need to scale the 'normal' rml value by the proportion of tracer to
          !! 12C, if 14C it is in units of 1E12*14C/C ratio which gives the RML in
          !! the amount of 14C. 
          tracerRML(i,j) = rmlveg(i,j) * tracerValue(i) 

          ! Remove C cost of BNF:
          if (rootmass_ns(i,j) > 0.0 .and. Ncycle_on) then
            rmrTracer(i,j) = rmrTracer(i,j) + (bnf_nat(i,j) + bnf_ant(i,j))*(C_cost_bnf_s/1000.*963.62)
          end if

          tracerRM(i,j) = tracerRML(i,j) + rmrTracer(i,j) + rmsTracer(i,j)
        
          !> Now we can find the tracerGPP by multiplying the gppveg value by the
          !! tracerValue.
          tracerGPP(i,j) = gppveg(i,j) * tracerValue(i)
          
          !> And from that calculate the tracerNPP, which is what we are needing.
          tracerNPP(i,j) = tracerGPP(i,j) - tracerRM(i,j)

          !> In rare circumstances the NPP of the tracer can be more than is possible considering the 
          !! proportion of the tracer in the atmosphere (assuming that the NPP drawn from a short-lived pool). 
          !! To prevent these excursions, check the tracer NPP to make sure it is not greater than could 
          !! be coming from the atmosphere. Any extra is taken to increase the respiratory costs in proportion to
          !! their amount of the total tracerRM. We only do this to the root and stem costs as the leaf is 
          !! already pinned to that from the total C pools.
          if (abs(nppveg(i,j)) > 0.) then 
            if (tracerNPP(i,j) < 0 .and. nppveg(i,j) > 0. .or. abs(tracerNPP(i,j) / nppveg(i,j)) >  tracerValue(i)) then   
              extratracer = max(0., (tracerNPP(i,j) - nppveg(i,j) * tracerValue(i)))
              tracerNPP(i,j) = nppveg(i,j) * tracerValue(i)
              tracerRespRootStem = rmrTracer(i,j) + rmsTracer(i,j)
              if (tracerRespRootStem > 0.) then 
                rmrTracer(i,j) = rmrTracer(i,j) + extratracer * (rmrTracer(i,j) / tracerRespRootStem)
                rmsTracer(i,j) = rmsTracer(i,j) + extratracer * (rmsTracer(i,j) / tracerRespRootStem)
              end if 
            end if 
          end if 
          
          if (tracerNPP(i,j) < 0 .and. nppveg(i,j) > 0.) print*,'Mismatch in tracer/total NPP!',j,tracerNPP(i,j),nppveg(i,j)

        end if

      end do ! loop 280
    end do ! loop 270

        !> Now that we know maintenance respiration from leaf, stem, and root
        !! and gpp, we can find growth respiration for each vegetation type
        if (useTracer /= 'None') then 
          call growthRespiration(il1,il2,sort,useTracer,nppveg,rgveg,tracerNPP,tracerRG)
        else 
          call growthRespiration(il1,il2,sort,useTracer,nppveg,rgveg)
        end if 

    do j = 1,icc
      do i = il1,il2
        ! Now update the NPP for the costs of growth respiration
        nppveg(i,j) = nppveg(i,j) - rgveg(i,j)
        if (useTracer /= 'None') tracerNPP(i,j) =  tracerNPP(i,j) - tracerRG(i,j)  
        
        !> Calculate grid/tile-averaged rates of rm, rg, npp, and gpp
        rml(i) = rml(i) + fcancmx(i,j) * rmlveg(i,j)
        rms(i) = rms(i) + fcancmx(i,j) * rmsveg(i,j)
        rmr(i) = rmr(i) + fcancmx(i,j) * rmrveg(i,j)
        rm(i) = rm(i) + fcancmx(i,j) * rmveg(i,j)
        rg(i) = rg(i) + fcancmx(i,j) * rgveg(i,j)
        npp(i) = npp(i) + fcancmx(i,j) * nppveg(i,j)
        gpp(i) = gpp(i) + fcancmx(i,j) * gppveg(i,j)
        autores(i) = rg(i) + rm(i)
        autoresveg(i,j) = rmveg(i,j) + rgveg(i,j)

      end do 
    end do
    !
    !>    Add moss GPP and rml to the grid/tile average C fluxes
    !>    for grid cells which have moss
    !
    do i = il1,il2
      if (mossPresent(i) /= 'None') then          
        if (mossPresent(i) == 'Sphagnum') then
          rgmoss(i) = anmoss(i) * grescoefmoss(1)
        else if (mossPresent(i) == 'Feather') then
          rgmoss(i) = anmoss(i) * grescoefmoss(2)
        else 
          print*,'Unknown moss type: ',mossPresent(i)
          call errorHandler('calcNPP', -1)
        end if 

        ! When N cycling is on, add the C cost of BNF to the growth respiration for moss.
        if (Ncycle_on) then 
          rgmoss(i) = rgmoss(i) + bnf_moss(i) * (C_cost_mossBNF/1000*963.62)
        end if 

        rml(i) = rml(i) + rmlmoss(i)
        rm(i) = rm(i)  + rmlmoss(i)
        rg(i) = rg(i)  + rgmoss (i)
        armoss(i) = rmlmoss(i) + rgmoss(i)
        nppmoss(i) = anmoss(i) - rgmoss(i)
        npp(i) = npp(i) + nppmoss (i)
        gpp(i) = gpp(i) + gppmoss(i)
        autores(i) = autores(i) + armoss(i)

        nppmosstep(i) = nppmoss(i) * (1.0/963.62) * deltat    ! kgC/m2/dt

      else
        nppmosstep(i) = 0.
      end if
    end do ! loop 335

  end subroutine calcNPP
  
  ! ----------------------------------------------------------------------------

  !> \ingroup ctemdriver_calcnep
  
  !> Calculate NEP (net ecosystem productivity), difference between NPP and
  !! heterotrophic respriation, for each pft and tile
  !!
  !! @author V. Arora, J. Melton
  subroutine calcNEP (il1, il2, ilg, nppveg, hetrsveg, fg, npp, hetrores, & ! In
                      nep, nepveg) ! Out

    use classicParams,   only : icc, iccp1, zero

    implicit none

    ! arguments
    integer, intent(in) :: il1             !< il1=1
    integer, intent(in) :: il2             !< il2=ilg (no. of grid cells in latitude circle)
    integer, intent(in) :: ilg
    real, intent(in) :: nppveg(:,:)    !< NPP for individual pfts,  (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(in) :: hetrsveg(:,:) !< Vegetation averaged litter and soil C respiration rates (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(in) :: fg(:)              !< Fraction of tile that is bare ground.
    real, intent(in) :: npp(:)            !< Tile-level net primary productivity (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(in) :: hetrores(:)       !< Tile-level heterotrophic respiration (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)

    real, intent(out) :: nep(ilg) !< Tile-level net ecosystem productivity (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)
    real, intent(out) :: nepveg(ilg,iccp1) !< PFT-level net ecosystem productivity (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)

    integer :: i, j

    nep(:) = 0.
    do i = il1,il2
      do j = 1,icc
        nepveg(i,j) = nppveg(i,j) - hetrsveg(i,j)
      end do
      if (fg(i) > zero) then
        nepveg(i,iccp1) = 0. - hetrsveg(i,iccp1)
      end if
      nep(i) = npp(i) - hetrores(i)
    end do

  end subroutine calcNEP
  
  ! ----------------------------------------------------------------------------
  !> \ingroup ctemdriver_calcnbp
  
  !> Calculate NBP (net biome production) for each pft by taking into account
  !! C emission losses. The disturbance routine produces emissions due to fire
  !! and while the land use change subroutine calculates emissions due to LUC.
  !! The LUC related combustion flux is assumed to be spread uniformly over the
  !! tile as it is no longer associated with any one PFT. To calculate the NBP
  !! we do not subtract LUC emissions from the PFT-level NBP but we do subtract
  !! it from the per tile NBP.
  !! @author V. Arora, J. Melton
  subroutine calcNBP (il1, il2, ilg, nepveg, fcancmx, & ! In
                      lucemcom, ltresveg, scresveg, nep, & ! In
                      glcaemls, blcaemls, stcaemls, rtcaemls, ltrcemls, & ! In/Out
                      nbpveg, dstcemls1, dstcemls3, nbp, dstcemls2, nepCMIP) ! Out

    use classicParams, only : icc, iccp1, iccp2, deltat

    implicit none

    ! arguments
    integer, intent(in) :: il1             !< il1=1
    integer, intent(in) :: il2             !< il2=ilg (no. of grid cells in latitude circle)
    integer, intent(in) :: ilg
    real, intent(in)    :: fcancmx(:,:)    !< max. fractional coverage of ctem's 9 pfts, but this can be
    !! modified by land-use change,and competition between pfts
    !real, intent(in) :: deltat             !< CTEM (biogeochemical) time step (days)
    real, intent(in) :: nepveg(:,:)        !< Net ecosystem productivity,  \f$\mu mol CO_2 m^{-2} s^{-1}\f$
    real, intent(in) :: lucemcom(:)        !< Land use change (LUC) related combustion emission losses,  \f$\mu mol CO_2 m^{-2} s^{-1}\f$
    real, intent(in) :: ltresveg(:,:,:)    !< Litter respiration for each pft, bare fraction, and LUC product pool, \f$\mu mol CO_2 m^{-2} s^{-1}\f$
    real, intent(in) :: scresveg(:,:,:)    !< Soil carbon respiration for each pft, bare fraction, and LUC product pool, \f$\mu mol CO_2 m^{-2} s^{-1}\f$
    real, intent(in) :: nep(:)             !< Net ecosystem productivity, tile average,  \f$\mu mol CO_2 m^{-2} s^{-1}\f$

    real, intent(inout) ::  glcaemls(:,:)  !< Green leaf carbon emission losses, \f$(kg C/m^2)\f$
    real, intent(inout) ::  blcaemls(:,:)  !< Brown leaf carbon emission losses, \f$(kg C/m^2)\f$
    real, intent(inout) ::  rtcaemls(:,:)  !< Root carbon emission losses, \f$(kg C/m^2)\f$
    real, intent(inout) ::  stcaemls(:,:)  !< Stem carbon emission losses, \f$(kg C/m^2)\f$
    real, intent(inout) ::  ltrcemls(:,:)  !< Litter carbon emission losses, \f$(kg C/m^2)\f$

    real, intent(out) :: nbpveg(ilg,iccp1) !< Net biome productivity,  \f$\mu mol CO_2 m^{-2} s^{-1}\f$
    real, intent(out) :: dstcemls1(ilg)    !< grid ave. carbon emission losses due to disturbance, vegetation,  \f$\mu mol CO_2 m^{-2} s^{-1}\f$
    real, intent(out) :: dstcemls3(ilg)    !< grid ave. carbon emission losses due to disturbance, litter,  \f$\mu mol CO_2 m^{-2} s^{-1}\f$
    real, intent(out) :: nbp(ilg)          !< Net biome productivity, tile average, \f$\mu mol CO_2 m^{-2} s^{-1}\f$
    real, intent(out) :: dstcemls2(ilg)     !yy< grid ave. carbon emission losses due to disturbance, total,  \f$\mu mol CO_2 m^{-2} s^{-1}\f$
    real, intent(out) :: nepCMIP(ilg)      !< FLAG rsc please add in description

    ! Local vars
    integer :: i, j
    real :: dscemlv1(ilg,icc)  !< Disturbance emission losses from plants, \f$(kg C/m^2)\f$
    real :: dscemlv2(ilg,icc)  !< Disturbance emission losses from plants and litter, \f$(kg C/m^2)\f$

    !--

    do i = il1,il2
      do j = 1,icc
        dscemlv1(i,j) = glcaemls(i,j) + blcaemls(i,j) + stcaemls(i,j) + rtcaemls(i,j)
        dscemlv2(i,j) = dscemlv1(i,j) + ltrcemls(i,j)

        ! Convert \f$(kg C/m^2)\f$ emitted in one day into u mol co2/m2.sec before
        ! subtracting emission losses from nep.
        nbpveg(i,j) = nepveg(i,j) - dscemlv2(i,j) * (963.62/deltat)

      end do ! loop 101

      ! For accounting purposes, we also need to account for the bare fraction
      ! NBP. Since there is no fire on the bare, we use 0.
      nbpveg(i,iccp1) = nepveg(i,iccp1)   - 0.

    end do ! loop 100
    !>
    !! Calculate grid. averaged rate of carbon emissions due to fire in u-mol co2/m2.sec.
    !! Convert all emission losses from \f$(kg C/m^2)\f$ emitted in 1 day to u-mol co2/m2.sec.
    !! Calculate grid averaged carbon emission losses from litter.
    !!
    dstcemls1 = 0.0
    dstcemls2 = 0.0
    do j = 1,icc
      do i = il1,il2
        dstcemls1(i) = dstcemls1(i) + fcancmx(i,j) * dscemlv1(i,j) * (963.62 / deltat)
        dstcemls2(i) = dstcemls2(i) + fcancmx(i,j) * dscemlv2(i,j) * (963.62 / deltat)
        glcaemls(i,j) = glcaemls(i,j) * (963.62/deltat)
        blcaemls(i,j) = blcaemls(i,j) * (963.62/deltat)
        stcaemls(i,j) = stcaemls(i,j) * (963.62/deltat)
        rtcaemls(i,j) = rtcaemls(i,j) * (963.62/deltat)
        ltrcemls(i,j) = ltrcemls(i,j) * (963.62/deltat)
      end do ! loop 104
    end do ! loop 103

    ! For the tile-level NBP, we include the disturbance emissions as well as
    ! respiration from the paper (litter) and furniture (soil carbon) pools (LUC
    ! product pools). Also include here the instantaneous emissions due to LUC.
    nbp(:) = 0.
    dstcemls3 = 0.0
    do i = il1,il2
      nbp(i) = nep(i) - dstcemls2(i) - (ltresveg(i,iccp2,1) + scresveg(i,iccp2,1)) - lucemcom(i)
      nepCMIP(i) = nep(i) - dstcemls2(i)
      dstcemls3(i) = dstcemls2(i) - dstcemls1(i)  ! litter is total - vegetation.
    end do ! loop 105

  end subroutine calcNBP
  

  !! This subroutine is meant to scale the sum of nppleaf(a), npproot(b), and
  !! nppstem (c) from d to e, post kludge 2. The scaled values of a, b, and c are then outputted.
  !! If the model conserved C this wouldn't be needed. But it doesn't at present. 
  !! The logic doesn't work when d and e are of different signs in addition to
  !! a, b, and c of being different sign. It can happen. In this case the subroutine prints
  !! a warning message. We will have to see how often does this happen.
  !! In any case nppleaf/root/stem or Tier 2 variables.

  subroutine scale_numbers(a, b, c, d, e, scaled_a, scaled_b, scaled_c)
    implicit none
    real, intent(in) :: a, b, c, d, e
    real, intent(out) :: scaled_a, scaled_b, scaled_c
    ! logical :: different_signs

    ! ! Initialize the warning flag
    ! different_signs = .false.

    ! ! Check if d and e have different signs
    ! if (d * e < 0.0) then
    !   ! Check if a, b, and c are not of the same sign
    !   if ((a > 0.0 .and. b < 0.0 .and. c > 0.0) .or. &
    !       (a < 0.0 .and. b > 0.0 .and. c < 0.0) .or. &
    !       (a > 0.0 .and. b < 0.0 .and. c < 0.0) .or. &
    !       (a < 0.0 .and. b > 0.0 .and. c > 0.0)) then
    !     different_signs = .true.
    !   end if
    ! end if

    ! ! Print warning if both conditions are met
    ! if (different_signs) then
    !   print *, "Warning: in scale_numbers in side ctemDriver, d and e have different signs, and a, b, c are not of the same sign!"
    ! end if

    ! Handle the case where all numbers are zero
    if (d == 0.0) then
      scaled_a = e / 3.0
      scaled_b = e / 3.0
      scaled_c = e / 3.0
    else
      ! Scale each number by e/d
      scaled_a = a * e / d
      scaled_b = b * e / d
      scaled_c = c * e / d
    end if

  end subroutine scale_numbers

  !> \namespace ctemdriver
  !! Central module that contains the ctem driver and associated subroutines.
  !!
  !! The basic model structure of CTEM includes three live vegetation components
  !! (leaf (L), stem (S) and root (R)) and two dead carbon pools (litter or
  !! detritus (D) and soil carbon (H)). The amount of carbon in these pools
  !! (\f$C_\mathrm{L}\f$, \f$C_\mathrm{S}\f$, \f$C_\mathrm{R}\f$, \f$C_\mathrm{D}\f$,
  !! \f$C_\mathrm{H}\f$, \f$kgC m^{-2}\f$) is tracked prognostically through the
  !! fluxes in and out of them. The rate change equations for carbon in these
  !! pools are summarized in Sect. \ref{rate_change_eqns} after the processes
  !! leading to the calculation of fluxes in and out of these pools are introduced
  !! in the following sections.
  !!
end module ctemDriver
