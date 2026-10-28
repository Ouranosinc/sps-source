!> \file                                                                               
!> Peatland specific parameterizations 
module peatlandsMod

  ! J. Melton. Sep 26, 2016

  ! J. Melton. Nov 6 2023 - pull out moss subroutines and make this only peatland-specific

  implicit none

  ! Subroutines contained in this module:
  public  :: peatDayEnd
  public  :: peatDepth
  public  :: peatStorage

contains

  ! ------------------------------------------------------------------
  !> \ingroup peatlandsmod_peatDayEnd
  
  !> At the end of the day update the degree days for moss photosynthesis and the peat bottom layer depth
  !> @author Yuanqiao Wu
  subroutine peatDayEnd (il1, il2, delzw, peatdep, soildpth, zbot)

    use classicParams,  only : ignd, ilg

    implicit none

    integer, intent(in) :: il1
    integer, intent(in) :: il2 

    real, dimension(ilg),      intent(in) :: peatdep        !< Depth of peat column [m]
    real, dimension(ignd),     intent(in) :: zbot           !< Bottom of soil layers (m)
    real, dimension(ilg,ignd), intent(inout) :: delzw       !< Permeable thickness of soil layer [m]
    real, dimension(ilg),      intent(inout) :: soildpth    !< soil depth (m)

    integer :: i, k
    integer :: botlyr

    do i = il1,il2 

      !> Update peatland bottom layer depth for fens and bogs       
        botlyr = 1
        do k = 1,ignd
          if (peatdep(i) < zbot(k)) exit
          botlyr = k
        end do
        
        ! Small error checking, don't let it set delzw to less than the second layer (generally 20 cm).
        ! if 1 is kept it will cause cause an indexing issue below.
        if (botlyr == 1) botlyr = 2
        
        ! The permeable thickness of the soil layer just above the peat is set to the 
        ! difference between the peat depth and the bottom of the layer above.
        delzw(i,botlyr) = peatdep(i) - zbot(botlyr-1)
        
        ! NOTE: The soil permeable depth is just the peat depth as we assume that the 
        ! total soil column in peatlands is peat. If useStaticPeatDep is true then this
        ! will just set soildpth to itself, since peatdep was set to soildpth up in ctemDriver.
        soildpth(i) = peatdep(i)

    end do

    return
  end subroutine peatDayEnd

  ! ---------------------------------------------------------------------------------------------------
  !> \ingroup peatlandsmod_peatDepth
  
  !> Calculate the peat depth based on equation 18 in Wu, Verseghy, Melton 2016 GMD.
  !! @author Y. Wu, J. Melton
  !!
  real function peatDepth (peatSoilC)

    implicit none

    real, intent(in) :: peatSoilC  !< Peat Soil C mass, \f$kg C/m^2\f$

    ! Calculate the peat depth based on equation 18 in Wu, Verseghy, Melton 2016 GMD.
    ! Since peatlands are run over tiles this uses the tile average soil mass (so only over
    ! the peat tile). 
    !NOTE peatSoilC is, itself, based on peat depth (see peatStorage below), so this is just 
    ! used to update the peat depth as the peat soil C pool changes

    peatDepth = ( -72067.0 + sqrt((72067.0 **2) - (4.0 * 4056.6 &
                * ( -peatSoilC * 1000 / 0.487)))) / (2. * 4056.6)
  end function peatDepth
  

  ! ---------------------------------------------------------------------------------------------------
  !> \ingroup peatlandsmod_peatStorage
  !! Finds the carbon storage in the peat based on depth (or oxic and anoxic compartments)
  !! The water table depth delineates the oxic and anoxic compartments.
  !! functions (R**2 = 0.9999) determines the carbon content of each
  !! compartment from a peat bulk density profile based on unpulished
  !! data from P.J.H. Richard (described in fig. 1, Frokling et al.(2001)
  !! conversion of peat into carbon with 48.7% (Mer Bleue unpublished data,
  !! Moore)
  !>
  !! @author Y. Wu, J. Melton
  !!
  real function peatStorage (depth)

    implicit none

    real, intent(in) :: depth  !< Peat compartment depth (m). Either total column or oxic/anoxic.

    peatStorage = (4056.6 * depth**2 + 72067.0 * depth) * 0.487 / 1000.0

  end function peatStorage
  
  ! ---------------------------------------------------------------------------------------------------

  !> \namespace peatlandsmod
  !! Peatland specific processes 
  !! @author Y. Wu, J. Melton
  !!
  !! The peatland module is published in Geoscientific Model Development (Wu et al. 2016) \cite Wu2016-zt.
  !!
  !! To account for the eco-hydrological and biogeochemical interactions among
  !! vegetation, atmosphere and soil in peatlands, the following modifications
  !! were made to the coupled CLASS3.6--CTEM2.0 modelling framework:
  !!
  !! 1. The top soil layer was characterized as a moss layer with a higher heat
  !! and hydraulic capacity than a mineral soil layer. The moss layer buffers the
  !! exchange of energy and water at the soil surface and regulates the soil
  !! temperature and moisture (Turetsky et al. 2012) \cite Turetsky2012-qh.
  !!
  !! 2. Three peatland vascular PFTs (evergreen shrubs, deciduous shrubs and
  !! sedges) as well as mosses were added to the existing nine CTEM PFTs. These
  !! peatland-specific PFTs are adapted to cold climate and inundated soil with
  !! optimized plant structure (shoot/root ratio, rooting depth), growth strategy
  !! and metabolic acclimations to light, water and temperature.
  !!
  !! 3. We considered the soil inundation stress on microbial respiration in the
  !! litter C pool. The original CTEM assumed that litter respiration was not
  !! affected by oxygen deficit as a result of flooding, since litter was always
  !! assumed to have access to air. This assumption does not hold for peatlands
  !! where high water table positions occur routinely.
  !!
  !! 4. To provide the framework for future runs coupled to the global earth system model, we separated the soil C balance and heterotrophic respiration
  !! (HR) calculations for peatland and non-peatland fractions for each grid cell
  !! in the global model. Over the non-peatland fraction, we use the original CTEM
  !! approach that aggregates the HR from each PFT weighted by the fractional
  !! cover. Over the peatland fraction the soil C pool and decomposition are
  !! controlled by the water table position, following the two-compartment
  !! approach used in the MWM (St. Hilaire et al. 2010) \cite St-Hilaire2010-5e9.
  !!
  !! The standard configuration of soil layers in CLASS consists of three layers with
  !! thickness of 0.10, 0.25, and 3.75, m. Organic soil in CLASS was parameterized
  !! by Letts et al. (2000) \cite Letts2000-pg as fibric, hemic and sapric peat in the three soil
  !! layers respectively, representing fresh, moderately decomposed and highly
  !! decomposed organic matter. Tests of CLASS on peatlands revealed improved
  !! performance in the energy simulations for fens and bogs with this organic
  !! soil parameterization. However, the model overestimated energy and water
  !! fluxes at bog surfaces during dry periods due to the neglect of the moss
  !! cover (Comer et al., 2000) \cite Comer2000-mz.
  !!
  !! To take into account the interaction amongst the moss and the soil layers and
  !! the overlying atmosphere for energy and water transfer, we added a new soil
  !! layer 0.10 m thick above the fibric organic soil to represent living and dead
  !! peatland bryophytes, such as Sphagnum mosses and true mosses
  !! (Bryopsida). The physical characteristics of mosses differ from those of
  !! either the shoots or the roots of vascular plants (Rice et al., 2008). In
  !! particular, mosses can hold more than 30 g of water per gram of biomass
  !! (Robroek et al., 2009). More than 90 \% of the moss leaf volume is
  !! occupied by the water-holding hyaline cells (Rice et al., 2008), which retain
  !! water even when the water table depth declines to 1--10 m below the surface
  !! (Hayward and Clymo, 1982).
  !!
  !! The parameter values of the moss layer for water and energy properties were
  !! derived from a number of recent experiments measuring the hydraulic
  !! properties of mosses (Price et al., 2008 \cite Price2008-fr; Price and Whittington, 2010
  !! McCarter and Price, 2012). Living mosses range from 2--3 to over
  !! 5 cm in height (Rice et al., 2008) and have lower values of dry bulk density
  !! and field capacity than fibric peat (Price et al., 2008) \cite Price2008-fr. Compared to fibric
  !! peat, the saturated hydraulic conductivity of living moss is higher by orders
  !! of magnitude (Price et al., 2008) \cite Price2008-fr and the thermal conductivity is more
  !! affected by the water content (O'Donnell et al., 2009). To fully account for
  !! the effect of mosses, we set the depth of the living moss (\f$z_{m}\f$)
  !! within the top soil (i.e. moss) layer to 3 cm for fens and 4 cm for bogs,
  !! and interpolated its water content \f$w_m\f$ (kg water /(kg dry mass)) from the water content of the overall layer
  !! \f$\theta_{l, 1}\f$ (m\f$^3\f$ water /(m soil)\f$^3\f$) and the depth of the
  !! living moss:
  !!
  !! \f$w_m =\frac{z_m\theta_{l, 1}\rho_w}{B_m} \f$
  !!
  !! where the dry moss biomass (\f$B_m) \f$ is converted from moss C
  !! (C\f$_m \f$) using the standard conversion factor of 0.46 kg C per kg dry
  !! biomass, \f$\theta _{l, 1} \f$ (m\f$^3\f$\, m\f$^{-3}\f$) is the liquid water
  !! content of the top soil layer, and \f$\rho_w \f$ is the density of water
  !! (1000 kg m\f$^{-3} \f$). The maximum and minimum moss water contents were
  !! estimated from a number of observed moss water contents (e.g. Williams and Flanagan, 1998; Robroek et al., 2009). In CLASS, evaporation at the soil
  !! surface is controlled by a soil evaporation efficiency coefficient \f$\beta \f$
  !! (Verseghy, 2012). This parameter is calculated from the liquid water content
  !! and the field capacity of the first soil layer following Lee and
  !! Pielke (1992). For peatlands, \f$\beta \f$ was assumed to be regulated by the
  !! relative moisture of the living moss rather than the ratio of relative liquid
  !! water content of the first soil layer:
  !!
  !! \f$ \beta = 0.25 [ 1- \cos \left( \frac{w_m
  !! -w_{m, min}}{w_m-w_{m, max}} \right)]^{2} \f$
  !!
  !! where \f$w_m \f$, \f$w_{m, max} \f$, and \f$w_{m, min} \f$ are the water content
  !! and the maximum and minimum water contents of the living moss in kg water / (kg dry mass).
  !!
  !!
  !> \file
end module peatlandsMod
