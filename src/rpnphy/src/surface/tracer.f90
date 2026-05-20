!> \file
!> Module containing all relevant subroutines for the model
!! tracer.
!! @author Joe Melton
!!

module tracerModule

  implicit none

  public :: prepTracer
  public :: decay14C
  public :: fractionate13C
  public :: convertTracerUnits
  public :: convertTracerOutput_Scalar
  public :: convertTracerOutput
  public :: checkTracerBalance

  ! Set to true if you want to check for tracer balance.
  ! NOTE: See the comments in the prepTracer header for how to set
  ! up run for this check.
  logical :: doTracerBalance  = .false. !< Logical to determine if the tracer pool balance check is performed.

contains

  !> \ingroup tracermodule
  
  !> Tracks C flow through the system.
  !! No fractionation effects are applied in this subroutine.
  !! The tracer's value depends on how
  !! the model is initialized and the input file used.
  !!
  !! The tracer trackes the C movement through the green leaf,
  !! brown leaf, root, stem, litter and soil C. Carbon that is
  !! incorporated into the plants are given a tracer value of
  !! tracerValue, which corresponds to that read in from the
  !! tracerCO2 file for the year simulated. As the simulation
  !! runs and C is transferred from the living pools to the detrital
  !! pools, the tracer also is transfered.
  !!
  !! The subroutine contains a tracer pools balance check that should be
  !! run before using the tracer subroutines, in case any other model
  !! developments have not be added to the fluxes for the tracer pools but
  !! have been applied elsewhere. To do the balance check, set useTracer
  !! to 1 (for simple tracer). This ensures no fractionation or decay is
  !! performed on the tracer. Ensure that the initFile tracer pool sizes
  !! are the same as the model C pools otherwise it will fail immediately.
  !! Next set checkTracerBalance to True in the prepTracer code. This
  !! replaces the read in tracerCO2 value with a value of 1, which means
  !! the tracer is given the same amount of C as the normal C pools. Then
  !! at the end of prepTracer a small section of code is called
  !! that checks that the tracer C pools are the same size as the model C pools.
  !! If they differ it indicates a missing flux term in prepTracer or
  !! that the initFile was not set up with identical tracer and normal C pools
  !! at the start of the balance check.
  !!
  !> @author Joe Melton
  subroutine prepTracer (il1, il2, ilg, tracerCO2, & ! In
                         stemmass_s, stemmass_ns, rootmass_s, rootmass_ns, litrmass, & ! In  
                         gleafmas_s, gleafmas_ns, bleafmas, soilcmas, peatlandType, & ! In
                         tracerValue, & ! Out
                         tracerStemMass, tracerRootMass, tracerLitrMass, & ! InOut
                         tracerGLeafMass, tracerBLeafMass, tracerSoilCMass) ! InOut

    use classicParams, only : icc, iccp2, ignd, zero

    implicit none

    integer, intent(in) :: il1 !< other variables: il1=1
    integer, intent(in) :: il2 !< other variables: il2=ilg
    integer, intent(in) :: ilg !< number of grid cells in this latitude band.

    real, intent(out) :: tracerValue(ilg) !< Temporary variable containing the tracer CO2 value to be used. (14C: $f1E12 ^{14}C/C\$f)
    real, intent(in) :: tracerCO2(:)      !< Tracer CO2 value read in from tracerCO2File, (simple: ppm, 14C: $f1E12 ^{14}C/C\$f)

    real, dimension(ilg,icc), intent(in) :: stemmass_s      !< Structural stem mass for each of the CTEM PFTs, \f$(kg c/m^2)\f$
    real, dimension(ilg,icc), intent(in) :: stemmass_ns     !< Non-structural stem mass for each of the CTEM PFTs, \f$(kg C/m^2)\f$
    real, dimension(ilg,icc), intent(in) :: rootmass_s      !< Structural root mass for each of the CTEM PFTs, \f$(kg C/m^2)\f$
    real, dimension(ilg,icc), intent(in) :: rootmass_ns     !< Non-structural root mass for each of the CTEM PFTs, \f$(kg C/m^2)\f$
    real, dimension(ilg,iccp2,ignd), intent(in) :: litrmass !< Litter mass for each of the CTEM PFTs + bare + LUC product pools, \f$(kg C/m^2)\f$
    real, dimension(ilg,icc), intent(in) :: gleafmas_s      !< Structural green leaf mass for each of the CTEM PFTs, \f$(kg C/m^2)\f$
    real, dimension(ilg,icc), intent(in) :: gleafmas_ns     !< Non-structural green leaf mass for each of the CTEM PFTs, \f$(kg C/m^2)\f$
    real, dimension(ilg,icc), intent(in) :: bleafmas        !< Brown leaf mass for each of the CTEM PFTs, \f$(kg C/m^2)\f$
    real, dimension(ilg,iccp2,ignd), intent(in) :: soilcmas !< Soil carbon mass for each of the CTEM PFTs + bare + LUC product pools, \f$(kg C/m^2)\f$

    real, intent(inout) :: tracerGLeafMass(:,:)      !< Tracer mass in the green leaf pool for each of the CTEM PFTs, 14C: \f$ng ^{14}C/m^2\f$ 
    real, intent(inout) :: tracerBLeafMass(:,:)      !< Tracer mass in the brown leaf pool for each of the CTEM PFTs, 14C: \f$ng ^{14}C/m^2\f$ 
    real, intent(inout) :: tracerStemMass(:,:)       !< Tracer mass in the stem for each of the CTEM PFTs, 14C: \f$ng ^{14}C/m^2\f$ 
    real, intent(inout) :: tracerRootMass(:,:)       !< Tracer mass in the roots for each of the CTEM PFTs, 14C: \f$ng ^{14}C/m^2\f$ 
    real, intent(inout) :: tracerLitrMass(:,:,:)     !< Tracer mass in the litter pool for each of the CTEM PFTs + bareground and LUC products, 14C: \f$ng ^{14}C/m^2\f$ 
    real, intent(inout) :: tracerSoilCMass(:,:,:)    !< Tracer mass in the soil carbon pool for each of the CTEM PFTs + bareground and LUC products, 14C: \f$ng ^{14}C/m^2\f$ 

    character(8), dimension(ilg), intent(in) :: peatlandType !< Peatland type (else 'None')

    ! Local
    integer :: i, j, k

    ! ---------

    !> If doTracerBalance is true, Check for mass balance for the tracer.
    !!  First set the tracer CO2 value to 1 so it gets the same inputs as the
    !! model CO2 pools. Otherwise just use the read in tracerValue.
    if (doTracerBalance) then
      tracerValue = 1.
    else
      tracerValue = tracerCO2
    end if

    !> We don't need to do any conversion of the carbon that is uptaked
    !! in this timestep. We assume that the tracerValue should just be applied as is.

    !> If the normal pools are empty, make the tracer pools the same.
    do i = il1,il2
      do j = 1,iccp2
        if (j <= icc) then ! these are just icc sized arrays.
          if (rootmass_ns(i,j) + rootmass_s(i,j) < zero) tracerRootMass(i,j) = 0.
          if (stemmass_ns(i,j) + stemmass_s(i,j) < zero) tracerStemMass(i,j) = 0.
          if (bleafmas(i,j) < zero) tracerBLeafMass(i,j) = 0.
          if (gleafmas_ns(i,j) + gleafmas_s(i,j) < zero) tracerGLeafMass(i,j) = 0.
        end if
        if (sum(litrmass(i,j,:)) < zero) tracerLitrMass(i,j,:) = 0.
        if (sum(soilcmas(i,j,:)) < zero) tracerSoilCMass(i,j,:) = 0.
      end do ! j

      if (peatlandType(i) /= 'None') print * ,'Tracer not set up yet for peatlands.'

    end do ! i

    return
  end subroutine prepTracer
  
  ! -------------------------------------------------------
  !> \ingroup tracermodule
  
  !> Calculates the decay of \f$^{14}C\f$ in the tracer pools.
  !!
  !! Once a year we calculate the decay of \f$^{14}C\f$ in the tracer pools.
  !! This calculation is only called if useTracer == '14C'.
  !! If using spinfast to equilibrate the model we need to
  !! adjust the decay of 14C in the model soil C pool. Since
  !! spinfast increases the turnover time of soil C by the factor
  !! 1/spinfast, the 14C in the soils generated with a spinfast > 1
  !! will be too young by the same factor so we increase the decay proportionally
  !> @author Joe Melton
  subroutine decay14C (il1, il2, tracer, vgat)

    use ctemStateVars, only : tracersType, veg_gat
    use classicParams, only : iccp2, icc, nlat, nmos, ignd, ilg, c_switch, lambda14C, zero

    implicit none

    integer, intent(in) :: il1 !< other variables: il1=1
    integer, intent(in) :: il2 !< other variables: il2=ilg
    type(tracersType(nlat, nmos, ignd, ilg)), intent(inout) :: tracer
    type(veg_gat(ilg, ignd)),                 intent(inout) :: vgat

    integer  :: i, j, k, l    ! counters
    real :: dfac    !< decay constant. \f$y^{-1}\f$

    ! Associate names with variables defined in derived types.
    associate( &
    tracerGLeafMass   => tracer%gLeafMassgat,           & !< real(:,:): Tracer mass in the green leaf pool for each of the CTEM pfts, 14C: \f$ng ^{14}C/m^2\f$  
    tracerBLeafMass   => tracer%bLeafMassgat,           & !< real(:,:): Tracer mass in the brown leaf pool for each of the CTEM pfts, 14C: \f$ng ^{14}C/m^2\f$  
    tracerStemMass    => tracer%stemMassgat,            & !< real(:,:): Tracer mass in the stem for each of the CTEM pfts, 14C: \f$ng ^{14}C/m^2\f$  
    tracerRootMass    => tracer%rootMassgat,            & !< real(:,:): Tracer mass in the roots for each of the CTEM pfts, 14C: \f$ng ^{14}C/m^2\f$  
    tracerLitrMass    => tracer%litrMassgat,            & !< real(:,:,:): Tracer mass in the litter pool for each of the CTEM pfts + bareground and LUC products, 14C: \f$ng ^{14}C/m^2\f$  
    tracerSoilCMass   => tracer%soilCMassgat,           & !< real(:,:,:): Tracer mass in the soil carbon pool for each of the CTEM pfts + bareground and LUC products, 14C: \f$ng ^{14}C/m^2\f$  
    tracerMossCMass   => tracer%mossCMassgat,           & !< real(:): Tracer mass in moss biomass, 14C: \f$ng ^{14}C/m^2\f$  
    tracerMossLitrMass => tracer%mossLitrMassgat,       & !< real(:): Tracer mass in moss litter, 14C: \f$ng ^{14}C/m^2\f$  
    spinfast          => c_switch%spinfast,              & !< integer: Set this to a higher number up to 10 to spin up 
                                                          !! soil carbon pool faster
    gleafmas          => vgat%gleafmas,                 & !< real(:,:): Green leaf mass for each of the CTEM PFTs, \f$kg c/m^2\f$
    bleafmas          => vgat%bleafmas,                 & !< real(:,:): Brown leaf mass for each of the CTEM PFTs, \f$kg c/m^2\f$
    stemmass          => vgat%stemmass,                 & !< real(:,:): Stem mass for each of the CTEM PFTs, \f$kg c/m^2\f$
    rootmass          => vgat%rootmass,                 & !< real(:,:): Root mass for each of the CTEM PFTs, \f$kg c/m^2\f$
    litrmass          => vgat%litrmass,                 & !< real(:,:,:): Litter for each of the CTEM PFTs, \f$kg c/m^2\f$
    soilcmas          => vgat%soilcmas,                 & !< real(:,:,:): Soil C mass for each of the CTEM PFTs, \f$kg c/m^2\f$
    Cmossmas          => vgat%Cmossmas,                 & !< real(:): Moss live C mass, \f$kg c/m^2\f$
    litrmsmoss        => vgat%litrmsmoss                & !< real(:): Green leaf mass for each of the CTEM PFTs, \f$kg c/m^2\f$
    )
    !----------
    ! begin calculations

    dfac = exp( - 1. / lambda14C)

    do i = il1,il2
      if (tracerMossCMass(i) > zero) then
        tracerMossCMass(i) = ((tracerMossCMass(i) + 1.) * dfac) - 1.
      else 
        tracerMossCMass(i) = 0.
      end if 
    
      ! If there is too little 14C, compared to the amount of total C (i.e. Delta 14C <= -1000), then
      ! set the 14C pool to zero.    
      if (convertTracerOutput_Scalar(tracerMossCMass(i),Cmossmas(i)) <= -1000.) tracerMossCMass(i) = 0.  

      if (tracerMossLitrMass(i) > zero) then 
        tracerMossLitrMass(i) = ((tracerMossLitrMass(i) + 1.) * dfac) - 1.
      else 
        tracerMossLitrMass(i) = 0.
      end if 
      do k = 1,ignd
        if (convertTracerOutput_Scalar(tracerMossLitrMass(i),litrmsmoss(i,k)) <= -1000.) tracerMossLitrMass(i) = 0.
      end do 
      do j = 1,iccp2
        if (j <= icc) then

          if (tracerGLeafMass(i,j) > zero) then
            tracerGLeafMass(i,j) = ((tracerGLeafMass(i,j) + 1.) * dfac) - 1.
          else 
            tracerGLeafMass(i,j) = 0.
          end if 

          if (convertTracerOutput_Scalar(tracerGLeafMass(i,j),gleafmas(i,j)) <= -1000.) tracerGLeafMass(i,j) = 0.

          if (tracerBLeafMass(i,j) > zero) then 
            tracerBLeafMass(i,j) = ((tracerBLeafMass(i,j) + 1.) * dfac) - 1.
          else 
            tracerBLeafMass(i,j) = 0.
          end if 

          if (convertTracerOutput_Scalar(tracerBLeafMass(i,j),bleafmas(i,j)) <= -1000.) tracerBLeafMass(i,j) = 0.

          if (tracerStemMass(i,j) > zero) then
            tracerStemMass(i,j) = ((tracerStemMass(i,j) + 1.) * dfac) - 1.
          else 
            tracerStemMass(i,j) = 0.
          end if 

          if (convertTracerOutput_Scalar(tracerStemMass(i,j),stemmass(i,j)) <= -1000.) tracerStemMass(i,j) = 0.

          if (tracerRootMass(i,j) > zero) then
            tracerRootMass(i,j) = ((tracerRootMass(i,j) + 1.) * dfac) - 1.
          else 
            tracerRootMass(i,j) = 0.
          end if 

          if (convertTracerOutput_Scalar(tracerRootMass(i,j),rootmass(i,j)) <= -1000.) tracerRootMass(i,j) = 0.

        end if
        do k = 1,ignd
          if (tracerLitrMass(i,j,k) > zero) then
            tracerLitrMass(i,j,k) = ((tracerLitrMass(i,j,k) + 1.) * dfac) - 1.
          else 
            tracerLitrMass(i,j,k) = 0.
          end if 

          if (convertTracerOutput_Scalar(tracerLitrMass(i,j,k),litrmass(i,j,k)) <= -1000.) tracerLitrMass(i,j,k) = 0.

          if (tracerSoilCMass(i,j,k) > zero) then
            do l = 1,spinfast 
              tracerSoilCMass(i,j,k) = ((tracerSoilCMass(i,j,k) + 1.) * dfac) - 1.
            end do
          else 
            tracerSoilCMass(i,j,k) = 0.
          end if 

          if (convertTracerOutput_Scalar(tracerSoilCMass(i,j,k),soilcmas(i,j,k)) <= -1000.) tracerSoilCMass(i,j,k) = 0.

        end do
      end do
    end do
 
    end associate
  end subroutine decay14C
  
  ! -------------------------------------------------------
  !> \ingroup tracermodule
  !! 13C tracer -- NOT IMPLEMENTED YET.
  subroutine fractionate13C

    use classicParams, only : c_switch

    implicit none

    print * ,'fractionate13C: Not implemented yet. Usetracer cannot == 13C !'
    call errorHandler('tracer', - 1)
  end subroutine fractionate13C
  
  ! -------------------------------------------------------
  !> \ingroup tracermodule
  
  !> Converts the units of the tracers, depending on the tracer
  !! being simulated.
  !!
  !! If the tracer is a simple tracer, no conversion of units is needed.
  !! If the tracer is \f$^{14}C\f$ then we apply a conversion as follows.
  !! Incoming units expected are \f$\Delta ^{14}C\f$ reported
  !! relative to the Modern standard, including corrections for
  !! age and fractionation following Stuiver and Polach 1977 \cite Stuiver1977-yj
  !!
  !! \f$\Delta ^{14}C = 1000 \left( \left[ 1 + \frac{\delta ^{14}C}{1000} \right]
  !! \frac{0.975^2}{1 + \frac{\delta ^{13}C}{1000}} - 1 \right)\f$
  !!
  !! Here \f$\delta ^{14}C\f$ is  the  measured value and \f$\Delta ^{14}C\f$ corrects
  !! for isotopic fractionation by mass-dependent processes. The 0.975 term is the
  !! fractionation of \f$^{13}C\f$ by photosynthesis. Since we have only one tracer
  !! we make a few simplifying assumptions. First we assume that no fractionation of
  !! \f$^{14}C\f$ occurs, thus 0.975 becomes 1 and the \f$\delta ^{13}C\f$ becomes 0.
  !! \f$\Delta ^{14}C\f$ then simplifies to \f$\delta ^{14}C\f$ and becomes:
  !!
  !! \f$\Delta ^{14}C = \delta ^{14}C = 1000 \left(\frac{A_s}{A_{abs}} - 1 \right)\f$
  !!
  !! where the \f$A_s\f$ is the \f$^{14}C/C\f$ ratio in a given sample. We follow Koven
  !! et al. (2013) \cite Koven2013-dd in assuming a background preindustrial atmospheric
  !! \f$^{14}C/C\f$ ratio (\f$A_{abs}\f$) of \f$10^{-12}\f$. Technically, the standard value of this is
  !! 0.2260 \f$\pm\f$ 0.0012 Bq/gC where a Bq is 433.2 x \f$10^{-15}\f$ mole (\f$^{14}C\f$).
  !! The \f$^{14}C/C\f$ ratio, or \f$A_s\f$ can then be solved to be:
  !!
  !! \f$ A_s = 10^{-12}\left( \frac{\Delta ^{14}C}{1000} + 1 \right)\f$
  !! 
  !! However, to avoid issues due to numerical precision, we remove the \f$10^{-12}\f$ and
  !! thus use \f$10^{12}^{14}C/C\f$ for all model calculations. This results in most quantities
  !! being tracked in units of \f$ng ^{14}C m^{-2}\f$ and \f$amol ^{14}CO_{2} m^{-2} s^{-1}\f$.
  !!
  !! @author Joe Melton
  function convertTracerUnits (tracerco2conc)

    use classicParams, only : nlat, nmos, c_switch

    implicit none

    real, dimension(nlat,nmos)       :: convertTracerUnits
    real, dimension(:,:), intent(in) :: tracerco2conc

    associate( &
    useTracer         => c_switch%useTracer             & !< character: Switch for use of a model tracer. If useTracer is 'None' then the
                                                          !! tracer code is not used. useTracer = 'Simple' turns on a simple tracer that tracks
                                                          !! pools and fluxes. The simple tracer then requires that the tracer values in
                                                          !! the init_file and the tracerCO2file are set to meaningful values for the experiment being run.
                                                          !! useTracer = '14C' means the tracer is 14C and will then call a 14C decay scheme.
                                                          !! useTracer = '13C' means the tracer is 13C and will then call a 13C fractionation scheme.
    )

    ! ------------------------

    select case (useTracer)

    case ('Simple') ! Simple tracer

      ! Simple tracer, no conversion of units needed.
      convertTracerUnits = tracerco2conc  

    case ('14C') ! 14C

      ! Incoming units expected are \Delta ^{14}C so need to convert to 14C/C ratio
      ! We exclude the multiplication by 1E-12 for numerical precision and thus produce
      ! 1E12 * 14C/C as our final units.
      convertTracerUnits = (tracerco2conc / 1000. + 1.) !* 1E-12 <- this factor is purposefully left out.

    case ('13C') ! 13C

      print * ,'convertTracerUnits: Error, 13C not implemented (useTracer == 3 used)'
      call errorHandler('tracer', 2)

    case default

      print * ,'Error,entered convertTracerUnits with a non-valid useTracer value', useTracer
      call errorHandler('tracer',1)

    end select

    end associate
  end function convertTracerUnits
  
  ! -------------------------------------------------------
  !> \ingroup tracermodule
  
  !> Converts the units of the tracers, depending on the tracer
  !! being simulated, for writing to output files.
  !!
  !! If the tracer is a simple tracer, no conversion of units is needed.
  !! If the tracer is \f$^{14}C\f$ then we undo the conversion from convertTracerUnits.
  !!
  !! @author Joe Melton
  function convertTracerOutput (tracerVal,cPoolVal)

    use classicParams, only : c_switch

    implicit none

    real, allocatable       :: convertTracerOutput(:)
    real, dimension(:), intent(in) :: tracerVal   !< Tracer pool
    real, dimension(:), intent(in) :: cPoolVal    !< Total C pool
    integer :: counter, k

    associate( &
    useTracer         => c_switch%useTracer             & !< character: Switch for use of a model tracer. If useTracer is 'None' then the
                                                          !! tracer code is not used. useTracer = 'Simple' turns on a simple tracer that tracks
                                                          !! pools and fluxes. The simple tracer then requires that the tracer values in
                                                          !! the init_file and the tracerCO2file are set to meaningful values for the experiment being run.
                                                          !! useTracer = '14C' means the tracer is 14C and will then call a 14C decay scheme.
                                                          !! useTracer = '13C' means the tracer is 13C and will then call a 13C fractionation scheme.
    )

    ! ------------------------

    allocate(convertTracerOutput(size(tracerVal)))
    counter = size(tracerVal)

    select case (useTracer)

    case ('Simple') ! Simple tracer

      ! Simple tracer, no conversion of units needed.
      convertTracerOutput = tracerVal  

    case ('14C') ! 14C

      ! Incoming units are 1E12 * 14C/C ratio so convert back to \Delta ^{14}C
      ! using the total C pool. A value of -1000 is used when there are no 14C present.
      do k = 1, counter
        if (cPoolVal(k) > 0.) then 
          convertTracerOutput(k) =((tracerVal(k) / cPoolVal(k)) - 1.) * 1000.
        else 
          convertTracerOutput(k) = -1000.
        end if 
      end do

    case ('13C') ! 13C

      print * ,'convertTracerSoilOutput: Error, 13C not implemented (useTracer == 13C used)'
      call errorHandler('tracer', 2)

    case default

      print * ,'Error,entered convertTracerSoilOutput with a non-valid useTracer value', useTracer
      call errorHandler('tracer',1)

    end select

    end associate
  end function convertTracerOutput
  

    ! -------------------------------------------------------
  !> \ingroup tracermodule
  
  !> Converts the units of the tracers, depending on the tracer
  !! being simulated, for scalars
  !!
  !! If the tracer is a simple tracer, no conversion of units is needed.
  !! If the tracer is \f$^{14}C\f$ then we undo the conversion from convertTracerUnits.
  !!
  !! @author Joe Melton
  function convertTracerOutput_Scalar (tracerVal,cPoolVal)

    use classicParams, only : c_switch

    implicit none

    real      :: convertTracerOutput_Scalar
    real, intent(in) :: tracerVal   !< Tracer pool
    real, intent(in) :: cPoolVal    !< Total C pool

    associate( &
    useTracer         => c_switch%useTracer             & !< character: Switch for use of a model tracer. If useTracer is 'None' then the
                                                          !! tracer code is not used. useTracer = 'Simple' turns on a simple tracer that tracks
                                                          !! pools and fluxes. The simple tracer then requires that the tracer values in
                                                          !! the init_file and the tracerCO2file are set to meaningful values for the experiment being run.
                                                          !! useTracer = '14C' means the tracer is 14C and will then call a 14C decay scheme.
                                                          !! useTracer = '13C' means the tracer is 13C and will then call a 13C fractionation scheme.
    )

    ! ------------------------


    select case (useTracer)

    case ('Simple') ! Simple tracer

      ! Simple tracer, no conversion of units needed.
      convertTracerOutput_Scalar = tracerVal  

    case ('14C') ! 14C

      ! Incoming units are 1E12 * 14C/C ratio so convert back to \Delta ^{14}C
      ! using the total C pool. A value of -1000 is used when there are no 14C present.
      if (cPoolVal > 0.) then 
        convertTracerOutput_Scalar =((tracerVal / cPoolVal) - 1.) * 1000.
      else 
        convertTracerOutput_Scalar = -1000.
      end if 

    case ('13C') ! 13C

      print * ,'convertTracerOutput_Scalar: Error, 13C not implemented (useTracer == 13C used)'
      call errorHandler('tracer', 2)

    case default

      print * ,'Error,entered convertTracerOutput_Scalar with a non-valid useTracer value', useTracer
      call errorHandler('tracer',1)

    end select

    end associate
  end function convertTracerOutput_Scalar
  
  ! -------------------------------------------------------
  !> \ingroup tracermodule
  
  !> Checks for balance between the tracer pools and the
  !! model normal C pools. The subroutine is called when
  !! the doTracerBalance logical is set to true in updateSimpleTracer
  !! and the model initFile is set up as described in the
  !! notes of updateSimpleTracer. checkTracerBalance uses the
  !! same tolerance for comparisions as balcar.
  !! @author Joe Melton
  subroutine checkTracerBalance (il1, il2, useTracer, &
                                 stemmass_s, stemmass_ns, rootmass_s, rootmass_ns, litrmass, &
                                 gleafmas_s, gleafmas_ns, bleafmas, soilcmas, &
                                 tracerGLeafMass, tracerStemMass, tracerRootMass, &
                                 tracerLitrMass, tracerBLeafMass, tracerSoilCMass)

    use classicParams, only : icc, iccp2, ilg, ignd, tolrance, c_switch

    implicit none

    integer, intent(in) :: il1 !< other variables: il1=1
    integer, intent(in) :: il2 !< other variables: il2=ilg
    
    real, dimension(ilg,icc), intent(in) :: stemmass_s      !< Structural stem mass for each of the CTEM PFTs, \f$kg c/m^2\f$ 
    real, dimension(ilg,icc), intent(in) :: stemmass_ns     !< Non-structural stem mass for each of the CTEM PFTs, \f$kg c/m^2\f$ 
    real, dimension(ilg,icc), intent(in) :: rootmass_s      !< Structrual root mass for each of the CTEM PFTs, \f$kg c/m^2\f$ 
    real, dimension(ilg,icc), intent(in) :: rootmass_ns     !< Non-structural root mass for each of the CTEM PFTs, \f$kg c/m^2\f$ 
    real, dimension(ilg,iccp2,ignd), intent(in) :: litrmass !< Litter mass for each of the CTEM PFTs + bare + LUC product pools, \f$kg c/m^2\f$ 
    real, dimension(ilg,icc), intent(in) :: gleafmas_s      !< Structural green leaf mass for each of the CTEM PFTs, \f$kg c/m^2\f$ 
    real, dimension(ilg,icc), intent(in) :: gleafmas_ns     !< Non-structural green leaf mass for each of the CTEM PFTs, \f$kg c/m^2\f$ 
    real, dimension(ilg,icc), intent(in) :: bleafmas        !< Brown leaf mass for each of the CTEM PFTs, \f$kg c/m^2\f$ 
    real, dimension(ilg,iccp2,ignd), intent(in) :: soilcmas !< Soil carbon mass for each of the CTEM PFTs + bare + LUC product pools, \f$kg c/m^2\f$ 
    real, intent(in) :: tracerGLeafMass(:,:)      !< Tracer mass in the green leaf pool for each of the CTEM pfts, 14C: \f$ng ^{14}C/m^2\f$  
    real, intent(in) :: tracerBLeafMass(:,:)      !< Tracer mass in the brown leaf pool for each of the CTEM pfts, 14C: \f$ng ^{14}C/m^2\f$  
    real, intent(in) :: tracerStemMass(:,:)       !< Tracer mass in the stem for each of the CTEM pfts, 14C: \f$ng ^{14}C/m^2\f$
    real, intent(in) :: tracerRootMass(:,:)       !< Tracer mass in the roots for each of the CTEM pfts, 14C: \f$ng ^{14}C/m^2\f$  
    real, intent(in) :: tracerLitrMass(:,:,:)     !< Tracer mass in the litter pool for each of the CTEM pfts + bareground and LUC products, 14C: \f$ng ^{14}C/m^2\f$  
    real, intent(in) :: tracerSoilCMass(:,:,:)    !< Tracer mass in the soil carbon pool for each of the CTEM pfts + bareground and LUC products, 14C: \f$ng ^{14}C/m^2\f$  

    integer :: i, j, k
    real :: limit

    character(6), intent(in) :: useTracer !< character: Switch for use of a model tracer. If useTracer is 'None' then the
                                          !! tracer code is not used. useTracer = 'Simple' turns on a simple tracer that tracks
                                          !! pools and fluxes. The simple tracer then requires that the tracer values in
                                          !! the init_file and the tracerCO2file are set to meaningful values for the experiment being run.
                                          !! useTracer = '14C' means the tracer is 14C and will then call a 14C decay scheme.
                                          !! useTracer = '13C' means the tracer is 13C and will then call a 13C fractionation scheme.

    ! ---------

    if (useTracer /= 'Simple') then
      print * ,'ERROR! checkTracerBalance: useTracer must == Simple for this check, you have:',useTracer
      print * ,'Either turn off checkTracerBalance (doTracerBalance == .false.) or set useTracer == Simple'
      print * ,'Ending run.'
      print * ,' *** ^^^ *** ^^^ *** ^^^ *** ^^^ *** ^^^ *** ^^^ '
      call errorHandler('checkTracerBalance', 1)
    end if
    limit = tolrance * 1E-9
    do i = il1,il2
      do j = 1,iccp2
        if (j <= icc) then

          ! green leaf mass
          if (abs(tracerGLeafMass(i,j) - gleafmas_s(i,j) - gleafmas_ns(i,j)) > limit) then
            print * ,'checkTracerBalance: Tracer balance fail for green leaf mass'
            print * ,i,'PFT',j,'tracerGLeafMass',tracerGLeafMass(i,j), &
                    'gleafmas',gleafmas_s(i,j) + gleafmas_ns(i,j)
            call errorHandler('checkTracerBalance', - 1)
          end if

          ! brown leaf mass
          if (abs(tracerBLeafMass(i,j) - bleafmas(i,j)) > limit) then
            print * ,'checkTracerBalance: Tracer balance fail for brown leaf mass'
            print * ,i,'PFT',j,'tracerBLeafMass',tracerBLeafMass(i,j), &
                    'bleafmas',bleafmas(i,j)
            call errorHandler('checkTracerBalance', - 2)
          end if

          ! stem mass
          if (abs(tracerStemMass(i,j) - stemmass_s(i,j) - stemmass_ns(i,j)) > limit) then
            print * ,'checkTracerBalance: Tracer balance fail for stem mass'
            print * ,i,'PFT',j,'tracerStemMass',tracerStemMass(i,j), &
                    'stemmass',stemmass_s(i,j) + stemmass_ns(i,j)
            call errorHandler('checkTracerBalance', - 1)
          end if

          ! root mass
          if (abs(tracerRootMass(i,j) - rootmass_s(i,j) - rootmass_ns(i,j)) > limit) then
            print * ,'checkTracerBalance: Tracer balance fail for root mass'
            print * ,i,'PFT',j,'tracerRootMass',tracerRootMass(i,j), &
                    'rootmass',rootmass_s(i,j) + rootmass_ns(i,j)
            call errorHandler('checkTracerBalance', - 1)
          end if
        end if

        do k = 1,ignd
          ! litter mass
          if (abs(tracerLitrMass(i,j,k) - litrmass(i,j,k)) > limit) then
            print*,'checkTracerBalance: Tracer balance fail for litter mass'
            print*,i,'PFT',j,'layer',k,'tracerLitrMass',tracerLitrMass(i,j,k), &
                      'litrmass',litrmass(i,j,k)
            call errorHandler('checkTracerBalance',-1)
          end if
          ! soil C mass
          if (abs(tracerSoilCMass(i,j,k) - soilcmas(i,j,k)) > limit) then
            print*,'checkTracerBalance: Tracer balance fail for soil C mass'
            print*,i,'PFT',j,'layer',k,'tracerSoilCMass',tracerSoilCMass(i,j,k), &
                      'soilcmas',soilcmas(i,j,k)
            call errorHandler('checkTracerBalance',-1)
          end if
        end do
      end do
    end do

    return
  end subroutine checkTracerBalance
  
  ! -------------------------------------------------------
  !> \namespace tracermodule
  !!
  !! Contains all relevant subroutines for the model tracer.
  !!
  !> \file
end module tracerModule
