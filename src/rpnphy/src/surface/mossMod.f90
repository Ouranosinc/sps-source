!> \file                                                                              
!> Moss specific parameterizations
module mossMod

  ! J. Melton. Nov 11, 2023 - created by splitting moss subroutines from peatlandsMod

  implicit none

  ! Subroutines contained in this module:
  public  :: mossPht
  public  :: updateMossC
  public  :: mossLitrFall

contains

  ! ------------------------------------------------------------------

  !> \ingroup peatlandsmod_mossPht
  
  !> Moss photosynthesis subroutine (equations are in module-level description)
  !> @author Yuanqiao Wu, J. Melton
  subroutine mossPht (il1, il2, iday, qswnv, thliq, co2conc, tsurfk, zsnow, & ! In 
                      pres, Cmossmas, dmoss, mossPresent, Nmossmas,  & ! In 
                      daylength, daylength_max, redcoeff_vcmaxMoss, Ncycle_on, & ! In 
                      ievapmoss, anmoss, rmlmoss, cevapmoss) ! Out 

    use classicParams,  only : rmlmoss25, tau25m, ektau, gasc, kc25, ko25, ec, ej, eo, evc, sj, &
                               hj, alpha_moss, thpmoss, thmmoss, ilg, ignd, TFREZ, RHOW, mossMaxSnow, &
                               mossMinTemp, pi, vcmax_moss, alpha_vcmax_moss, beta_vcmax_moss

    implicit none

    ! arguments:

    integer, intent(in) ::  il1                     !<
    integer, intent(in) ::  il2                     !<
    integer, intent(in) ::  iday                    !< model day of year
    character(8), dimension(ilg), intent(in) :: mossPresent  !< Type of moss present (else 'None')
    real, dimension(ilg), intent(in) :: qswnv       !< visible short wave radiation = qswnv in energBalNoVegSolve
                                                    !! and qswnvg in energBalVegSolve (W/m2)
    real, dimension(ilg,ignd), intent(in) :: thliq  !<
    real, dimension(ilg), intent(in) :: zsnow       !< Snow depth (m)
    real, dimension(ilg), intent(in) :: co2conc     !<
    real, dimension(ilg), intent(in) :: daylength   !< daylength for this location (hours)
    real, dimension(ilg), intent(in) :: daylength_max !< maximum daylength for this location (hours)
    real, dimension(ilg), intent(in) :: tsurfk      !< grid average ground surface temprature in K
    real, dimension(ilg), intent(in) :: pres        !<
    real, dimension(ilg), intent(in) :: Cmossmas    !< unit kg moss C updated in ctem
    real, dimension(ilg), intent(in)  :: Nmossmas   !< moss biomass N pool, \f$g N/m^2\f$
    real, dimension(ilg), intent(in) :: dmoss       !< unit m, depth of living moss. assume = 2 cm
                                                    !! can be related to mmoss as a variable
    real, dimension(ilg), intent(in) :: redcoeff_vcmaxMoss !< Reduction coefficient of the VCMAX0 for moss
    logical, intent(in) :: Ncycle_on                !< True if simulation includes Nitrogen Cycle processes
    
    integer, dimension(ilg), intent(out) ::  ievapmoss   !< Value is 0 is no evaporation from moss, 1 otherwise
    real, dimension(ilg), intent(out) :: anmoss     !< net photosynthesis (umol CO2/m2/s)
    real, dimension(ilg), intent(out) :: rmlmoss    !< moss autotrophic respiration (umol CO2/m2/s)
    real, dimension(ilg), intent(out) :: cevapmoss  !< evaporation coefficent for moss surface

    ! temporary variables depending on moss type
    real :: rmlmoss25Temp        !> Temp var for rmlmoss25
    real :: alpha_mossTemp       !> Temp var for alpha_moss
    real :: thpmossTemp          !> Temp var for thpmoss 
    real :: thmmossTemp          !> Temp var for thmmoss
    real :: alpha_vcmax_mossTemp !> Temp var for alpha_vcmax_moss
    real :: beta_vcmax_mossTemp  !> Temp var for beta_vcmax_moss
    real :: vcmax_mossTemp       !> Temp var for vcmax_moss

    ! Local variables

    integer :: i               !>
    integer :: j               !>
    real :: mmoss(ilg)       !<  dry moss biomasss (kg dry moss biomass)
    real :: parm(ilg)        !< par at the ground (moss layer) surface umol/m2/s
    real :: tsurf(ilg)       !< grid average ground surface temperature in C
    real :: wmoss(ilg)       !< water content extraporated from the surface
                             !! humidity qg and thliq of the first soil layer
                             !! unit kg water/ kg dw
    real :: wmosmin(ilg)     !< residual water content kg water /kg moss
    real :: wmosmax(ilg)     !< maximum water content kg water /kg moss
    real :: fwmoss(ilg)      !< relative water content of mosses in g fw /g dw
    real :: dsmoss(ilg)      !< degree of moss saturation = relative water
                             !! content/maximum relative water content
    real :: g_moss(ilg)      !< moss conductance umol CO2/m2/s (based on
                             !! Williams and Flanagan, 1998 for Sphagnum)
    real :: mwce (ilg)       !< moisture function of dark respiration of moss
    real :: tmoss(ilg)       !< moss temperature extraporated from the tbar 1
                             !! and grid averaged ground surface temperature, tsurf
    real :: tmossk(ilg)      !< moss temperature in K
    real :: q10rmlmos(ilg)   !< temperature function of the moss dark respiration
    real :: gamma(ilg)       !< compensation point for gross photosynthesis (Pa)
    real :: o2(ilg)          !< partial presure of oxygen (Pa)
    real :: co2a(ilg)        !< partical pressure of co2 (pa) same as in PHTSYN
    real :: tau(ilg)         !< arrhenius funciton of temperature
    real :: kc(ilg)          !< kinetic coeffficient of CO2 for photosynthesis
    real :: ko(ilg)          !< kinetic coeffficient of O2 for photosynthesis
    real :: bc(ilg)          !< coefficient used for wc
    real :: vcmax25(ilg)     !< seasonal varied maximum carboylation at 25
                             !! sphagnum (fig. 6, Williams and Flanagan, 1998)
    real :: vcmax(ilg)       !< max carboxylation rate (umol/m2/s)
    real :: jmax25(ilg)      !< maximum electorn transport rate at 25 degrees (umol/m2/s)
    real :: jmax(ilg)        !< maximum electorn transport rate (umol/m2/s)
    real :: wj(ilg)          !< net co2 assimilation rate limited by
                             !! electron transport (umol/m2/s)=jE in PHTSYN
    real :: wc(ilg)          !< net co2 assimilation rate limited by
    real :: ws(ilg)          !< net co2 assimilation rate limited by
                             !! sucrose availability (umol/m2/s)=JE IN PHTSYN
    real :: photon(ilg)      !< electron transport rate (umol/m2/s)
    real :: term1(ilg)       !< temporary terms for photosynthesis calculations
    real :: term2(ilg)       !< temporary terms for photosynthesis calculations
    real :: term3(ilg)       !< temporary terms for photosynthesis calculations
    real :: psna(ilg)        !< coefficients for quadratric solution of net photosynthesis
    real :: psnb(ilg)        !< coefficients for quadratric solution of net photosynthesis
    real :: psne(ilg)        !< coefficients for quadratric solution of net photosynthesis
    real :: mI(ilg)          !< coefficients of the solutions for net psn
    real :: mII(ilg)         !< coefficients of the solutions for net psn

    ! Local parameters:
    real, parameter :: tref = 298.16    !< unit K

    !     PHOTOSYNTHESIS COUPLING OR CURVATURE COEFFICIENTS
    ! real, parameter :: BETA1 = 0.950
    ! real, parameter :: BETA2 = 0.990

    ! ...........................
    ! Begin calculations:

    do i = il1,il2 ! loop 100
      wj(i)     = 0.0
      ws(i)     = 0.0
      wc(i)     = 0.0
      anmoss(i) = 0.0
      rmlmoss(i) = 0.0
      cevapmoss(i) = 0.0
    end do ! loop 100

    do i =  il1,il2 ! loop 200

      if (mossPresent(i) /= 'None') then ! Only do for regions with moss. 
        
        ! Depending on the moss type, we pick different parameters
        if (mossPresent(i) == 'Sphagnum') then 
          rmlmoss25Temp = rmlmoss25(1)
          alpha_mossTemp = alpha_moss(1)
          thpmossTemp = thpmoss(1)
          thmmossTemp = thmmoss(1)
          vcmax_mossTemp = vcmax_moss(1)
          alpha_vcmax_mossTemp = alpha_vcmax_moss(1)
          beta_vcmax_mossTemp  = beta_vcmax_moss(1)
        else if (mossPresent(i) == 'Feather') then ! Feather
          rmlmoss25Temp = rmlmoss25(2)
          alpha_mossTemp = alpha_moss(2)
          thpmossTemp = thpmoss(2)
          thmmossTemp = thmmoss(2)
          vcmax_mossTemp = vcmax_moss(2)
          alpha_vcmax_mossTemp = alpha_vcmax_moss(2)
          beta_vcmax_mossTemp  = beta_vcmax_moss(2)
        else
          print*,'Unknown moss type in mossPht', mossPresent(i)
          call errorHandler('mossMod', - 1)
        end if

        tsurf(i) = tsurfk(i) - tfrez
        tmossk(i) = tsurfk(i)
        tmoss(i) = tsurf(i)

        !>    phenology water factor on mosses, grow when temperature > mossMinTemp
        !!    and snowpack < mossMaxSnow
        
        if (zsnow(i) <= mossMaxSnow(1) .and. tsurf (i) >= mossMinTemp(1)) then   

          !> Find the light level (parm) at the ground surface for moss photosynthesis,  
          !! and a scaling factor degree of saturation of the moss layer phenology       
          !! parm is in umol/m2/s and converted from qswnv in W/m2
          parm(i) = qswnv(i) * 4.6
          o2(i)  = 20.9/100.0 * pres(i)
          co2a(i) = co2conc(i)/1000000.0 * pres(i)

          !> Water content used for the living moss (depth dmoss)
          !! dmoss is an input and site specific. Preferably make dmoss a function
          !! of Cmoss and mossPresent
          !! observed range of wmoss: 5 to 40 in Robrek (2007, 2009), 5 to 25          
          !! (Flanagen and Williams 1998)
          !! dmoss is between 2.5 to 5cm based on the species (Lamberty et al. 2006)    

          ! Convert from moss C to moss dry biomass using a standard conversion of 0.46 kg C 
          ! per kg dry biomass
          mmoss(i) = Cmossmas(i)/0.46
          
          ! Find the water content of the moss layer based on the total layer content and 
          ! the depth of the living moss. It is constrained to be within the minumum and 
          ! maximum water holding capacity of the moss.
          wmoss(i) = thliq(i,1) * rhow/(mmoss(i)/dmoss(i))
          wmosmax(i) = min(45.0,thpmossTemp * dmoss(i) * rhow/mmoss(i))    
          wmosmin(i) = max(5.0,thmmossTemp * dmoss(i) * rhow/mmoss(i))
          wmoss(i) = min(wmosmax(i),max(wmosmin(i),wmoss(i)))
          
          ! Lastly convert to the relative water content.
          fwmoss(i) = wmoss(i) + 1.     ! g fresh weight /g dry weight

          !> Find moss conductance (g_moss) in umol CO2/m2/s
          !! (Williams and Flanagan, 1998 for Sphagnum). follow MWM, fwmoss is
          !! the mosswat_fd in MWM. Empirical equation is only valid up to
          !! fwmoss=13, above 13 apply a linear extension to the equation.

          if (fwmoss(i) <= 13.0) then                                            
            g_moss(i) = - 0.195 + 0.134 * fwmoss(i) - 0.0256 * (fwmoss(i)) &
                        ** 2 + 0.00228 * (fwmoss(i)) ** 3 - 0.0000984 * &
                        (fwmoss(i)) ** 4 + 0.00000168 * (fwmoss(i)) ** 5
          else
            g_moss (i) = - 0.000447 * fwmoss(i) + 0.0489
          end if

          ! Rose: same structure as above for sphagnum, but with feather moss equation. also from Williams & Flanagan, 1998.
          ! make if/else statement with if ipeatland == 1 or 2 then sphagnum, else feather? (currently in an "if 1 or 2" statement)
          
          ! if (fwmoss(i) <= 12.0) then     ! Rose: 12.0? sphagnum: W&F98 say 17 but here it's 13, and linear instead of cst after...
            !g_moss = - 0.0194 + 0.021 * fwmoss(i) - 0.00446 * (fwmoss(i)) &
                      !** 2 + 0.000455 * (fwmoss(i)) ** 3 - 0.0000248 * &
                      !(fwmoss(i)) ** 4 + 0.000000592 * (fwmoss(i)) ** 5
          !else
            ! g_moss (i) = (linear equation?)
          !end if

          g_moss (i) = g_moss(i) * 1000000.0
          g_moss(i) = max(0.0,g_moss(i))

          !> Find moss surface evaporation coefficient
          !! controled by the degree of saturation in moss, pass to energBalVegSolve and energBalNoVegSolve
          !! apply a similar equation of soil surface cevap in energyBudgetPrep
          !! CEVAP = 0.25*[1 – cos(THLIQ*pi/THFC)]^2

          ! Fist get degree of moss saturation                        
          dsmoss(i) = (wmoss(i) - wmosmin(i))/(wmosmax(i) - wmosmin(i))
          
          ! then dsmoss is used to set cevap to 0, 1, or somewhere in between.
          if (dsmoss(i) <   0.001) then
            ievapmoss(i) = 0
            cevapmoss(i) = 0.
          else if (dsmoss(i) >= 1.0) then
            ievapmoss(i) = 1
            cevapmoss(i) = 1.0
          else
            ievapmoss(i) = 1
            cevapmoss(i) = 0.25 * (1.0 - cos(pi * dsmoss(i))) ** 2
          end if


          !> Calculate the moss water content effect on dark respiration
          !! in MWM and PDM an optimal wmoss is at 5.8 gw/gdw(fig. 2e, Frolking et al., 1996)
          !! Recent studies show weak but significant increases of sphagnum dark respiration     
          !! with moss water content above 5.8 gw/gdw (Adkinson and Humphreys, 2011 and ref.)
          !! this change has improved the ER simulation greatly
          if (wmoss(i) < 0.4) then
            mwce (i) = 0.0
          else if (wmoss(i) < 5.8 .and. wmoss(i) > 0.4) then
            mwce(i) = 0.35 * wmoss(i) ** (2.0/3.0) - 0.14
          else
            mwce(i) = 0.01 * wmoss(i) + 0.942
          end if

          !> Moss dark respiration
          !! observed range of rmlmoss 0.60 to 1.60 umol/m2/s (e.g. Adkinson 2006)       

          q10rmlmos(i) = (3.22 - (0.046 * tmoss(i))) ** ((tmoss(i) - 25.0)/10.0)
          rmlmoss(i) = rmlmoss25Temp * mwce(i) * q10rmlmos(i)

          !> Moss photosynthesis
          !! calculate bc (coefficient used for Wc, limited by Rubisco)

          tau(i) = tau25m * exp((tmossk(i) - tref) * ektau/(tref * gasc * tmossk(i)))     
          gamma(i) = 0.5 * o2(i)/ tau(i)
          kc(i) = kc25 * exp((tmossk(i) - tref) * ec/(tref * gasc * tmossk(i)))
          ko(i) = ko25 * exp((tmossk(i) - tref) * eo/(tref * gasc * tmossk(i)))
          bc(i)  = kc(i) * (1.0 + (o2(i)/ ko(i)))

          if (Ncycle_on) then 

            !> Allow the N status of the moss to influence the Vcmax.
            if (Cmossmas(i) > 0. .and. Nmossmas(i) > 0.) then 

              vcmax25(i) = alpha_vcmax_mossTemp * Nmossmas(i) + beta_vcmax_mossTemp
              vcmax25(i) = max(vcmax25(i), 0.5 * vcmax_mossTemp)
              vcmax25(i) = min(vcmax25(i), 1.5 * vcmax_mossTemp)

              if (redcoeff_vcmaxMoss(i) > 0.) vcmax25(i) = vcmax25(i) * redcoeff_vcmaxMoss(i)

            else ! Following the convention with the vascular PFTs, we just set to the max.
              vcmax25(i) = vcmax_mossTemp
            end if 

          else ! Nitrogen cycle off so just set to prescribed Vcmax25 value.

              vcmax25(i) = vcmax_mossTemp

          end if 

          !! Account for the seasonal change of Vcmax (e.g. fig. 6, Williams and Flanagan, 1998)
          !! Here we use the same adjustments as done for the vascular PFTs. It appears to give 
          !! a good approximation of the seasonal changes in the figure mentioned above.
          select case (mossPresent(i))
          case ('Sphagnum') 
            vcmax25(i) = vcmax25(i) * (daylength(i) / daylength_max(i)) ** 2
          case ('Feather')
            ! Following Fig 6 in Williams and Flanagan, 1998, we assume no seasonal 
            ! variation in Vcmax of feather. So do nothing.
          case default 
            print*,'Unknown moss type in mossPht', mossPresent(i)
            call errorHandler('mossMod', - 1)
          end select  

          ! Adjust the Vcmax to the temperature conditions
          vcmax(i) = vcmax25(i) * exp((tmossk(i) - tref) * evc/(tref * gasc * tmossk(i)))

          !> Calculate ws (phototysnthesis rate limited by transport capacity)
          !! = js in photosynCanopyConduct
          
          if (qswnv(i) > 0.) then
            ws(i) = 0.5 * vcmax(i)
          end if


          !> Calculate the maximum electron transport rate Jmax (umol/m2/s)
          !! 1.67 = vcmax25m/jmax25m ratio                                     

          jmax25(i) = 1.67 * vcmax25(i)
          term1(i) = exp(((tmossk(i)/tref) - 1.) * ej/(gasc * tmossk(i)))
          term2(i) = 1. + exp(((tref * sj) - hj)/(tref * gasc))
          term3(i) = 1. + exp(((sj * tmossk(i)) - hj)/(gasc * tmossk(i)))
          jmax(i) = jmax25(i) * term1(i) * term2(i) * term3(i)
           
          ! find the  electron trasport rate in mosses
          if (jmax(i) > 0.0) then
            photon(i) = alpha_mossTemp * parm(i)/sqrt(1.0 + (alpha_mossTemp ** 2 * parm(i) ** 2/(jmax(i) ** 2)))
          else
            photon(i) = 0.0    
          end if

          !> Calculate Wj, Wc (Farquhar and Caemmerer 1982)
          !! wj = light limited, = je in photosynCanopyConduct

          wj(i) = photon(i) * (co2a(i) - gamma(i))/(4. * co2a(i) + (8. * gamma(i)))

          !> Carboxylase(rubisco) limitation = jc in photosynCanopyConduct
          wc(i) = vcmax(i) * (co2a(i) - gamma(i))/(co2a(i) + bc(i))

          !> Choose the minimum of Wj and Wj both having the form:
          !! W = (a Ci - ad) / (e Ci + b)
          !! Then set psna, psnb, psnd and psne for the quadratic solution for net photosynthesis.

          if (wj(i) < wc(i)) then           
            psnb(i) = 8. * gamma(i)
            psna(i) = photon(i)
            psne(i) = 4.0
          else if (wc(i) < wj(i)) then
            psnb(i) = bc(i)
            psna(i) = vcmax(i)
            psne(i) = 1.0
          end if

          !> Calculate net and gross photosynthesis by solve the quadratic equation
          !! first root of solution is net photosynthesis An= min(Wj,Wc) - Rd
          !! gross photosynthesis GPP = min(Wc,Wj) = An + Rd

          if (psna(i) > 0.0) then                                                
            mI(i) = rmlmoss(i) - &
                    (psnb(i) * g_moss(i)/pres(i)/psne(i)) - &
                    (co2a(i) * g_moss(i)/pres(i)) - (psna(i)/psne(i))

            mII(i) = (psna(i) * co2a(i) * g_moss(i)/pres(i)/psne(i)) - &
                    (rmlmoss(i) * co2a(i) * g_moss(i)/pres(i)) - &
                     (rmlmoss(i) * psnb(i) * g_moss(i)/pres(i)/psne(i)) - &
                      (psna(i) * gamma(i) * &
                     g_moss(i)/pres(i)/psne(i))
          else
            mI(i) = 0.0
            mII(i) = 0.0
          end if
          anmoss(i) = ( - mI(i) - (mI(i) * mI(i) - 4 * mII(i)) ** 0.5)/2
          anmoss(i) = min(anmoss(i),ws(i))
        end if !check if the phenology permits photosynthesis
      end if ! mossPresent
    end do ! loop 200
    return
  end subroutine mossPht
  
  ! ---------------------------------------------------------------------------------------------------
  !> \ingroup mossMod_updateMossC
  
  !> At the end of the day update the moss carbon pool
  !> @author Joe Melton 
  subroutine updateMossC(Cmossmas,nppmosstep,litrfallmoss)
  
    real, intent(in) :: nppmosstep(:) !< moss npp \f$(kg C/m^2/timestep)\f$
    real, intent(in) :: litrfallmoss(:) !< moss litter fall \f$(kg C/m^2/timestep)\f$
    real, intent(inout) :: Cmossmas(:) !< C in moss biomass, \f$kg C/m^2\f$
    
    Cmossmas = Cmossmas + nppmosstep - litrfallmoss
  
  end subroutine updateMossC
  

  ! ---------------------------------------------------------------------------------------------------
  !> \ingroup mossmod_mossLitrFall
  
  !> Calculate the moss litter fall depending on which type of moss is present.
  !! @author Y. Wu, J. Melton
  real function mossLitrFall (Cmossmas,mossPresent,ldoy)

    use classicParams, only : rmortmoss, deltat

    implicit none

    real, intent(in) :: Cmossmas  !< Moss biomass C pool (kgC/m2)
    character(8), intent(in) :: mossPresent  !< Type of moss present (else 'None')
    real, intent(in) :: ldoy   !< Last day of the year

    if (mossPresent == 'Sphagnum') then
      mossLitrFall = Cmossmas * rmortmoss(1) / ldoy * deltat ! kgC/m2/day(dt)
      ! if (useTracer > 0) FLAG,not connected up !
    else if (mossPresent == 'Feather') then
      mossLitrFall = Cmossmas * rmortmoss(2) / ldoy * deltat ! kgC/m2/day(dt)
      ! if (useTracer > 0) FLAG,not connected up !
    else  
      print*,'Moss type unknown',mossPresent
      call errorHandler('mossLitrFall (mossMod)', -1)
    end if

  end function mossLitrFall
  



! ---------------------------------------------------------------------------------------------------
  !> \namespace mossMod
  !! Moss specific processes 
  !! @author Y. Wu, J. Melton, R. Lefebvre
  !!
  !! Moss photosynthesis subroutine
  !!
  !! Mosses are an important contributor to the primary production and the C
  !! sequestration in peatlands, owing to the low decomposability of the moss
  !! tissue. Sphagnum in peatlands grows at
  !! 20--1600 g biomass m\f$^{-2} \f$/yr and accounts for about 50 \% of
  !! the total peat volume (Turetsky, 2003). We have modified CTEM to include a
  !! moss C pool and moss litter pool along with the related C fluxes, i.e.
  !! photosynthesis, autotrophic respiration, heterotrophic respiration, and
  !! humification. The net photosynthesis of moss (\f$G_m)\f$ is calculated
  !! from the gross photosynthesis (\f$G_{0, m})\f$ and dark respiration
  !! (\f$R_{d, m}\f$):
  !!
  !! \f$ G_m=G_{0, m}-R_{d, m} \f$.
  !!
  !! The moss photosynthesis and dark respiration are calculated using the
  !! Farquhar~(1989) biochemical approach following the MWM (St-Hilaire et al.,
  !! 2010) \cite St-Hilaire2010-5e9 and CTEM (Melton and Arora, 2016) \cite Melton2016-zx, with modifications for integration
  !! with CLASS--CTEM and moss phenology. The leaf-level gross photosynthesis rate
  !! \f$G_{0, m} \f$ (\f$\mu \f$ mol CO\f$_2 \f$ /m \f$^2 \f$/s) is obtained as
  !! the minimum of the transportation limited photosynthesis rates (\f$J_s)\f$ and
  !! the first root of the quadratic solution of the light-limited rate (\f$J_e)\f$
  !! and the Rubisco limited rate (\f$J_c)\f$. A logistic factor (\f$\varsigma\f$) is
  !! added with values 0 or 1 to introduce a seasonal control of moss
  !! photosynthesis. In the MWM, spring photosynthesis starts when the snow depth
  !! is below 0.05 m and the soil temperature at 5 cm depth goes above
  !! 0.5 C (Moore et al., 2006). Since in our case CLASS sets the
  !! minimum depth for melting, discontinuous snow to 0.10 m, this limits the
  !! spring photosynthesis to starting only once the snow is completely melted.
  !!
  !! \f$ G_{0, m} = \quad \varsigma \min \left(J_{s}, \frac{(J_{c}+J_{e})\pm
  !! \sqrt{{(J_c +J_ e)}^2 - 4(J_c +J_e)} }{2}\right) \f$
  !!
  !! The dark respiration in mosses (\f$R_{d, m})\f$ is calculated as a function
  !! of the base dark respiration rate (\f$R_{d, m, 0})\f$, which has a value of
  !! 1.1 (\f$\mu \f$ mol CO\f$_2 \f$ /m \f$^2 \f$/s) (Adkinson and Humphreys, 2011) scaled
  !! by the moss moisture (\f$f_{m, rd}\f$) and soil temperature functions
  !! (\f$f_{T, rd})\f$. The moss moisture function is based on the volumetric
  !! water content of the moss, \f$\theta_m\f$ (kg water /(kg dry
  !! mass)). The MWM models the relation between water content in mosses and dark
  !! respiration with optimal water content at 5.8 g water per g dry weight,
  !! following the approach in Frolking et al. (1996). We modified the relation
  !! for water content above the optimal water content, based on a recent
  !! discovery of a weak linear positive relation between the dark respiration
  !! rate and the water content above the optimal water content during the late
  !! summer and fall (Adkinson and Humphreys, 2011):
  !!
  !! \f$R_{d, m} = R_{d, m, 0} f_{m, rd} f_{T, {rd}}\f$
  !!
  !! \f$ f_{T, {rd}} = (3.22-(0.046 \cdot
  !! T_{moss})^{(T_{moss}-25/10)} \f$
  !!
  !! \f$ f_{m, rd} = 0\f$ for\f$ \theta_{m} <0.4 \f$
  !!
  !! \f$ f_{m, rd} = 0.35 \theta_{m}^{2/3}-0.14\f$ for\f$ 0.4\le \theta_m<5.8 \f$
  !!
  !! \f$ f_{m, rd} = 0.01 \theta_{m} +0.942\f$ for \f$5.8<\theta _{m} \f$
  !!
  !! Photosynthetic photon flux density (PPFD) is measured by the
  !! photosynthetically active radiation (PAR), which is defined as the solar
  !! radiation between 0.4 to 0.7 \f$\mu \f$ mol that can be used by plants via
  !! photosynthesis. In the coupled CLASS--CTEM system, the PAR received by the
  !! moss (PAR\f$_m\f$, unit \f$\mu \f$mol photons m\f$^-2\f$/s is
  !! converted from the visible short-wave radiation reaching the ground (\f$K_{\ast g}\f$,
  !! unit W/m\f$^2\f$) in CLASS by a factor of
  !! 4.6 \f$\mu\f$ mol /m\f$^2\f$/s per W/m\f$^2\f$ (McCree, 1972).
  !! \f$K_{\ast g}\f$ is a function of the incoming short-wave radiation
  !! (\f$K\downarrow \f$, unit: W/m\f$^2\f$), the surface albedo (\f$\alpha_g\f$),
  !! and the canopy transmissivity (\f$\tau _c\f$):
  !!
  !! \f$K_{\ast g} = K \downarrow \tau_{c}\left( 1-\alpha_{g}\right)\f$
  !!
  !! The energy uptake by the moss layer is thus a function of the total incoming
  !! short-wave radiation, the aggregated LAI of the PFTs present, the snow depth,
  !! the fractional vegetation cover, and the soil water content (Verseghy, 2012).
  !! In peatland C models that do not consider vegetation dynamics, the
  !! transmissivity of the vegetation canopy is usually assumed to be constant
  !! (e.g. St-Hilaire et al., 2010) \cite St-Hilaire2010-5e9. Compared with such models, CLASS enables a
  !! more detailed representation of light incident on the moss surface since it
  !! includes partitioning of direct/diffuse and visible/near-IR radiation,
  !! PFT-specific transmissivities, and time-varying LAI and fractional PFT
  !! coverages (Verseghy, 2012) \cite Verseghy2012-c0e.
  !!
  !> \file
end module mossMod
