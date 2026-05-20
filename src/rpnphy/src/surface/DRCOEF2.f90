!> \file
!> Calculates drag coefficients and related variables
!! This work has been described in \cite Abdella1996-em.
!! @author K. Abdella, N. Macfarlane, M. Lazare, D. Verseghy, E. Chan
!
subroutine DRCOEF2(cdm,cdh,rib,cflux,qg,qa,zomin,zohin, &
                        crib,tvirtg,tvirta,va,fi,iter, &
                        ilg,il1,il2)
  !
  !     * NOV 04/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
  !     * SEP 10/02 - K.ABDELLA.  BUGFIX IN CALCULATION OF "OLS" (2 PLACES).
  !     * MAY 22/02 - N.MCFARLANE.USE THE ASYMPTOTIC VALUE FOR CDH
  !     *                         GENERALLY TO LIMIT FLUXES AND PREVENT
  !     *                         OVERSHOOT. (COMMMENTED OUT FOR
  !     *                         OFF-LINE RUNS.)
  !     * APR 11/01 - M.LAZARE.   SHORTENED "CLASS2" COMMON BLOCK.
  !     * OCT 26/99 - E. CHAN.    COMMENT OUT ARTIFICIAL DAMPING OF
  !     *                         TURBULENT FLUXES FOR STAND-ALONE TESTING.
  !     * JUL 18/97 - M. LAZARE, CLASS 2.7. PASS IN ADDITIONAL WORK ARRAYS
  !     *             D.VERSEGHY. ZOMIN AND ZOHIN TO USE INTERNALLY,
  !     *                         SO THAT INPUT ZOMIN AND ZOHIN DO NOT
  !     *                         CHANGE WHEN PASSED BACK TO THE
  !     *                         ITERATION (BUG FROM PREVIOUS
  !     *                         VERSIONS). PREVIOUS ZOM AND ZOH
  !     *                         BECOME WORK ARRAYS.
  !     * MAR 10/97 - M. LAZARE.  PASS IN QG AND QA AND USE TO ONLY
  !     *                         DAMP TURBULENT FLUXES USING CDHMOD
  !     *                         UNDER STABLE CONDITIONS IF MOISTURE
  !     *                         FLUX IS UPWARDS. ALSO, BETTER
  !     *                         DEFINITION OF SURFACE LAYER TOP
  !     *                         FROM K.ABDELLA.
  !     * MAY 21/96 - K. ABDELLA. MODIFICATION FOR FREE-CONVECTIVE
  !     *                         LIMIT ON UNSTABLE SIDE ADDED.
  !     * JAN 10/96 - K. ABDELLA. CORRECT ERROR IN AU1 (UNSTABLE
  !     *                         SIDE) AND PUT IN PRANDTL NUMBER
  !     *                         RANGE (0.74->1)/
  !     *                         "CDHMOD" USED ON BOTH STABLE AND
  !     *                         UNSTABLE SIDE, TO LIMIT FLUXES
  !     *                         OVER LONG TIMESTEP.
  !     * M. LAZARE - FEB 14/95.  USE VARIABLES "ZOLN" AND "ZMLN"
  !     *                         ON UNSTABLE SIDE, FOR OPTIMIZATION.
  !     *                         THIS IS PREVIOUS VERSION "DRCOEFX".
  !     * K. ABDELLA/M. LAZARE. - NOV 30/94.
  !
  use classicParams, only : GRAV, VKC, rmos1, rmos2, use_tke

  implicit none

  !     * calculates drag coefficients and related variables for class.

  !     * output fields are:
  !     *    cdm    : stability-dependent drag coefficient for momentum.
  !     *    cdh    : stability-dependent drag coefficient for heat.
  !     *    rib    : bulk richardson number.
  !     *    cflux  : cd * mod(v), bounded by free-convective limit.

  !     * input fields are:
  !     *    zomin/: roughness heights for momentum/heat normalized by
  !     *    zohin   reference height.
  !     *    crib   : -rgas*slthkef/(va**2), where
  !     *             slthkef=-log(max(sgj(ilev),shj(ilev)))
  !     *    tvirtg : "SURFACE" virtual temperature.
  !     *    tvirta : lowest level virtual temperature.
  !     *    va     : amplitude of lowest level wind.
  !     *    fi     : fraction of surface type being studied.
  !     *    qg     : saturation specific humidity at ground temperature.
  !     *    qa     : lowest level specific humidity.
  !     *    iter   : index array indicating if point is undergoing
  !     *             further iteration or not.

  !     *    zom/  : work arrays used for scaling zomin/zohin
  !     *    zoh     on stable side, as part of calculation.
  !     * integer :: constants.

  integer, intent(in) :: ilg   !< Total number of atmospheric columns \f$[unitless]\f$
  integer, intent(in) :: il1   !< Index of first atmospheric column for calculations \f$[unitless]\f$
  integer, intent(in) :: il2   !< Index of last atmospheric column for calculations \f$[unitless]\f$

  !     * output arrays.

  real, intent(out) :: cdm    (ilg) !< stability-dependent drag coefficient for momentum
 ! real, intent(inout) :: z0_momentum    (ilg) !< surface roughness length for momentum
 ! real, intent(inout) :: z0_heat    (ilg) !< surface roughness length for heat
  real, intent(out) :: cdh    (ilg) !< stability-dependent drag coefficient for heat
  real, intent(inout) :: rib    (ilg) !< Bulk Richardson number
  real, intent(out) :: cflux  (ilg) !< cd * mod(v), bounded by free-convective limit: term in momentum flux (surface stress) calculations

  !     * input arrays.

  real, intent(in) :: zomin  (ilg) !< roughness height for momentum normalized by reference height
  real, intent(in) :: zohin  (ilg) !< roughness height for heat normalized by reference height
  real, intent(in) :: crib   (ilg) !< -RGAS*SLTHKEF/(VA**2), WHERE SLTHKEF=-LOG(MAX(SGJ(ILEV), SHJ(ILEV)))
  real, intent(in) :: tvirtg (ilg) !< "SURFACE" VIRTUAL TEMPERATURE
  real, intent(in) :: tvirta (ilg) !< LOWEST LEVEL VIRTUAL TEMPERATURE
  real, intent(in) :: va     (ilg) !< AMPLITUDE OF LOWEST LEVEL WIND
  real, intent(in) :: fi     (ilg) !< FRACTION OF SURFACE TYPE BEING STUDIED
  real, intent(in) :: qg     (ilg) !< SATURATION SPECIFIC HUMIDITY AT GROUND TEMPERATURE
  real, intent(in) :: qa     (ilg) !< LOWEST LEVEL SPECIFIC HUMIDITY

  integer, intent(in) :: iter(ilg) !< INDEX ARRAY INDICATING IF POINT IS UNDERGOING FURTHER ITERATION OR NOT

  !     * work arrays.
  !> ZOM/ZOH: WORK ARRAYS USED FOR SCALING ZOMIN/ZOHIN ON STABLE SIDE, AS PART OF CALCULATION
  integer :: i !<
  real :: zom    (ilg) !<
  real :: zoh    (ilg) !<

  !     * temporary variables.

  real :: aa !<
  real :: aa1 !<
  real :: beta !<
  real :: pr !<
  real :: zlev !<
  real :: zs !<
  real :: zoln !<
  real :: zmln !<
  real :: cpr !<
  real :: zi !<
  real :: olsf !<
  real :: olfact !<
  real :: zl !<
  real :: zmol !<
  real :: zhol !<
  real :: xm !<
  real :: xh !<
  real :: bh1 !<
  real :: bh2 !<
  real :: bh !<
  real :: wb !<
  real :: wstar !<
  real :: rib0 !<
  real :: wspeed !<
  real :: au1 !<
  real :: ols !<
  real :: psim1 !<
  real :: psim0 !<
  real :: psih1 !<
  real :: psih0 !<
  real :: ustar !<
  real :: tstar !<
  real :: wts !<
  real :: as1 !<
  real :: as2 !<
  real :: as3 !<
  real :: climit !<

  !     * additional temporary variables for tke.
  real :: etae,chie,phime,phifac1,phifac2,phihe,ffe
  real :: a1,a2,dffe,chif

  !-------------------------------------------------------------
  aa = 9.5285714
  aa1 = 14.285714
  beta = 1.2
  pr = 1.
  !

  if (use_tke) then

   aa=2.381
   do i=il1,il2
    if(fi(i).gt.0. .and. iter(i).eq.1) then
      rib(i)=crib(i)*(tvirtg(i)-tvirta(i))
      zlev=-crib(i)*tvirta(i)*(va(i)**2)/grav
      if(rib(i).ge.0.0) then
        zs=max(10.,5.*max(zomin(i)*zlev, zohin(i)*zlev))
        zs=zlev*(1.+rib(i))/(1.+(zlev/zs)*rib(i))
        zom(i)=zomin(i)*zlev/zs
        zoh(i)=zohin(i)*zlev/zs
        rib(i)=rib(i)*zs/zlev
      else
        zom(i)=zomin(i)
        zoh(i)=zohin(i)
      endif
      zoln=log(zoh(i))
      zmln=log(zom(i))
      if(rib(i).lt.0.0) then
        cpr=min(max(zoln/zmln,0.74),1.0)
        zi=1000.0
        olsf=beta**3*zi*vkc**2/zmln**3
        olfact=1.7*(log(1.+zom(i)/zoh(i)))**0.5+0.9
        olsf=olsf*olfact
        zl = -crib(i)*tvirta(i)*(va(i)**2)/grav
        zmol=zom(i)*zl/olsf
        zhol=zoh(i)*zl/olsf
        xm=(1.00-rmos1*zmol)**(0.250)
        xh=(1.00-rmos2*zhol)**0.25
        bh1=-log(-2.41*zmol)+log(((1.+xm)/2.)**2*(1.+xm**2)/2.)
        bh1=bh1-2.*atan(xm)+atan(1.)*2.
        bh1=bh1**1.5
        bh2=-log(-0.25*zhol)+2.*log(((1.00+xh**2)/2.00))
        bh=vkc**3.*beta**1.5/(bh1*(bh2)**1.5)
        wb=sqrt(grav*(tvirtg(i)-tvirta(i))*zi/tvirtg(i))
        wstar=bh**(0.333333)*wb
        rib0=rib(i)

        wspeed=sqrt(va(i)**2+(beta*wstar)**2)
        rib(i)=rib0*va(i)**2/wspeed**2
        au1=1.+5.0*(zoln-zmln)*rib(i)*(zoh(i)/zom(i))**0.25
        ols=-rib(i)*zmln**2/(cpr*zoln)*(1.0+au1/(1.0-rib(i)/(zom(i)*zoh(i))**0.25))
        if (ols > 0.066) then
          print * ,'Warning: OLS in DRCOEF2.f90 is > 1/15, which should cause a crash in the psim1 calculation. It might not crash, when using the intel compiler, &
          but results might still be unreasonable. If this happens, check for example whether the ZRFM and ZRFH values are reasonable.'
        end if  
        psim1=log(((1.00+(1.00-rmos1*ols)**0.250)/2.00)**2* &
             (1.0+(1.00-rmos1*ols)**0.5)/2.0)-2.0*atan(    &
             (1.00-rmos1*ols)**0.250)+atan(1.00)*2.00
        psim0=log(((1.00+(1.00-rmos1*ols*zom(i))**0.250)/2.00)**2 &
             *(1.0+(1.00-rmos1*ols*zom(i))**0.5)/2.0)-2.0* &
             atan((1.00-rmos1*ols*zom(i))**0.250)+atan(1.00)*2.0
        psih1=log(((1.00+(1.00-rmos2*ols)**0.50)/2.00)**2)
        psih0=log(((1.00+(1.00-rmos2*ols*zoh(i))**0.50)/2.00)**2)

        ustar=vkc/(-zmln-psim1+psim0)
        tstar=vkc/(-zoln-psih1+psih0)
        cdh(i)=ustar*tstar/pr
        wts=cdh(i)*wspeed*(tvirtg(i)-tvirta(i))
        wstar=(grav*zi/tvirtg(i)*wts)**(0.333333)

        wspeed=sqrt(va(i)**2+(beta*wstar)**2)
        rib(i)=rib0*va(i)**2/wspeed**2
        au1=1.+5.0*(zoln-zmln)*rib(i)*(zoh(i)/zom(i))**0.25
        ols=-rib(i)*zmln**2/(cpr*zoln)*(1.0+au1/(1.0-rib(i)/(zom(i)*zoh(i))**0.25))
        if (ols > 0.066) then
          print * ,'Warning: OLS in DRCOEF2.f90 is > 1/15, which should cause a crash in the psim1 calculation. It might not crash, when using the intel compiler, &
          but results might still be unreasonable. If this happens, check for example whether the ZRFM and ZRFH values are reasonable.'
        end if  
        psim1=log(((1.00+(1.00-rmos1*ols)**0.250)/2.00)**2* &
             (1.0+(1.00-rmos1*ols)**0.5)/2.0)-2.0*atan( &
             (1.00-rmos1*ols)**0.250)+atan(1.00)*2.00
        psim0=log(((1.00+(1.00-rmos1*ols*zom(i))**0.250)/2.00)**2 &
             *(1.0+(1.00-rmos1*ols*zom(i))**0.5)/2.0)-2.0* &
             atan((1.00-rmos1*ols*zom(i))**0.250)+atan(1.00)*2.0
        psih1=log(((1.00+(1.00-rmos2*ols)**0.50)/2.00)**2)
        psih0=log(((1.00+(1.00-rmos2*ols*zoh(i))**0.50)/2.00)**2)

      else

        wspeed=va(i)
        as1=10.0*zmln*(zom(i)-1.0)
        as2=5.00/(2.0-8.53*rib(i)*exp(-3.35*rib(i))+0.05*rib(i)**2)
        as2=as2*pr*sqrt(-zmln)/2.
        as3=27./(8.*pr*pr)
        etae=4.*rib(i)
        chie=etae*((zmln**2)/(-pr*zoln)+as3*etae)
        phime=-zmln+chie*(1.-zom(i))
        phifac1=(1.+2.*chie/3.)**0.5
        phifac2=(1.+2.*chie*zoh(i)/3.)**0.5
        phihe=-zoln+phifac1**3-phifac2**3
        ffe=chie*phihe/phime**2-etae
        a1=phifac1-zoh(i)*phifac2
        a2=(a1/phihe-2.*(1-zom(i))/phime)
        dffe=(phihe/phime**2)*(1.+chie*a2)
        chif=chie-ffe/dffe
        ! ** note that in the formulae below chif = 4z(ilev)/l where l is the monin-obukov length. **
        ! ** that is, ols=chif/4., where ols is the the ratio of the lowest model level height (z1) to the monin-obukov length.
        psim1=-chif
        psim0=-chif*zom(i)
        psih1=-(1.+2.*chif/3.)**1.5 +1.
        psih0=-(1.+3.*chif*zoh(i)/3.)**1.5 +1.
      endif

      ustar=vkc/(-zmln-psim1+psim0)
      tstar=vkc/(-zoln-psih1+psih0)

      cdm(i)=ustar**2.0
      cdh(i)=ustar*tstar/pr
      !
      !z0_momentum(i)=exp(zmln)*zlev
      !z0_heat(i)=exp(zoln)*zlev
      !* calculate cd*mod(v) under free-convective limit.
      !
      if(tvirtg(i).gt.tvirta(i))    then
        climit=1.9e-3*(tvirtg(i)-tvirta(i))**0.333333
      else
        climit=0.
      endif
      cflux(i)=max(cdh(i)*wspeed,climit)
    endif
    enddo
  endif
!
  if (.not. use_tke) then
    aa=9.5285714
    do i=il1,il2
      if (fi(i)>0. .and. iter(i)==1) then
        rib(i)=crib(i)*(tvirtg(i)-tvirta(i))
        zlev=-crib(i)*tvirta(i)*(va(i)**2)/grav
        if (rib(i)>=0.0) then
          zs=max(10.,5.*max(zomin(i)*zlev, zohin(i)*zlev))
          zs=zlev*(1.+rib(i))/(1.+(zlev/zs)*rib(i))
          zom(i)=zomin(i)*zlev/zs
          zoh(i)=zohin(i)*zlev/zs
          rib(i)=rib(i)*zs/zlev
        else
          zom(i)=zomin(i)
          zoh(i)=zohin(i)
      end if
        zoln=log(zoh(i))
        zmln=log(zom(i))
        if (rib(i)<0.0) then
          cpr=max(zoln/zmln,0.74)
          cpr=min(cpr,1.0)
          zi=1000.0
          olsf=beta**3*zi*vkc**2/zmln**3
          olfact=1.7*(log(1.+zom(i)/zoh(i)))**0.5+0.9
          olsf=olsf*olfact
          zl = -crib(i)*tvirta(i)*(va(i)**2)/grav
          zmol=zom(i)*zl/olsf
          zhol=zoh(i)*zl/olsf
          xm=(1.00-15.0*zmol)**(0.250)
          xh=(1.00-9.0*zhol)**0.25
          bh1=-log(-2.41*zmol)+log(((1.+xm)/2.)**2*(1.+xm**2)/2.)
          bh1=bh1-2.*atan(xm)+atan(1.)*2.
          bh1=bh1**1.5
          bh2=-log(-0.25*zhol)+2.*log(((1.00+xh**2)/2.00))
          bh=vkc**3.*beta**1.5/(bh1*(bh2)**1.5)
          wb=sqrt(grav*(tvirtg(i)-tvirta(i))*zi/tvirtg(i))
          wstar=bh**(0.333333)*wb
          rib0=rib(i)

          wspeed=sqrt(va(i)**2+(beta*wstar)**2)
          rib(i)=rib0*va(i)**2/wspeed**2
          au1=1.+5.0*(zoln-zmln)*rib(i)*(zoh(i)/zom(i))**0.25
          ols=-rib(i)*zmln**2/(cpr*zoln)*(1.0+au1/ &
              (1.0-rib(i)/(zom(i)*zoh(i))**0.25))
          if (ols > 0.066) then
            print * ,'Warning: OLS in DRCOEF2.f90 is > 1/15, which should cause a crash in the psim1 calculation. It might not crash, when using the intel compiler, &
          but results might still be unreasonable. If this happens, check for example whether the ZRFM and ZRFH values are reasonable.'
          end if  
          psim1=log(((1.00+(1.00-15.0*ols)**0.250)/2.00)**2* &
               (1.0+(1.00-15.0*ols)**0.5)/2.0)-2.0*atan( &
               (1.00-15.0*ols)**0.250)+atan(1.00)*2.00
          psim0=log(((1.00+(1.00-15.0*ols*zom(i))**0.250)/2.00)**2 &
               *(1.0+(1.00-15.0*ols*zom(i))**0.5)/2.0)-2.0* &
               atan((1.00-15.0*ols*zom(i))**0.250)+atan(1.00)*2.0
          psih1=log(((1.00+(1.00-9.0*ols)**0.50)/2.00)**2)
          psih0=log(((1.00+(1.00-9.0*ols*zoh(i))**0.50)/2.00)**2)

          ustar=vkc/(-zmln-psim1+psim0)
          tstar=vkc/(-zoln-psih1+psih0)
          cdh(i)=ustar*tstar/pr
          wts=cdh(i)*wspeed*(tvirtg(i)-tvirta(i))
          wstar=(grav*zi/tvirtg(i)*wts)**(0.333333)

          wspeed=sqrt(va(i)**2+(beta*wstar)**2)
          rib(i)=rib0*va(i)**2/wspeed**2
          au1=1.+5.0*(zoln-zmln)*rib(i)*(zoh(i)/zom(i))**0.25
          ols=-rib(i)*zmln**2/(cpr*zoln)*(1.0+au1/ &
              (1.0-rib(i)/(zom(i)*zoh(i))**0.25))
          if (ols > 0.066) then
            print * ,'Warning: OLS in DRCOEF2.f90 is > 1/15, which should cause a crash in the psim1 calculation. It might not crash, when using the intel compiler, &
          but results might still be unreasonable. If this happens, check for example whether the ZRFM and ZRFH values are reasonable.'
          end if      
          psim1=log(((1.00+(1.00-15.0*ols)**0.250)/2.00)**2* &
               (1.0+(1.00-15.0*ols)**0.5)/2.0)-2.0*atan( &
               (1.00-15.0*ols)**0.250)+atan(1.00)*2.00
          psim0=log(((1.00+(1.00-15.0*ols*zom(i))**0.250)/2.00)**2 &
               *(1.0+(1.00-15.0*ols*zom(i))**0.5)/2.0)-2.0* &
               atan((1.00-15.0*ols*zom(i))**0.250)+atan(1.00)*2.0
          psih1=log(((1.00+(1.00-9.0*ols)**0.50)/2.00)**2)
          psih0=log(((1.00+(1.00-9.0*ols*zoh(i))**0.50)/2.00)**2)

      else

          wspeed=va(i)
          as1=10.0*zmln*(zom(i)-1.0)
          as2=5.00/(2.0-8.53*rib(i)*exp(-3.35*rib(i))+0.05*rib(i)**2)
          !< <<
          as2=as2*pr*sqrt(-zmln)/2.
          as3=27./(8.*pr*pr)
          !> >>
          ols=rib(i)*(zmln**2+as3*as1*(rib(i)**2+as2*rib(i))) &
             /(as1*rib(i)-pr*zoln)
          psim1=-0.667*(ols-aa1)*exp(-0.35*ols)-aa-ols
          psim0=-0.667*(ols*zom(i)-aa1)*exp(-0.35*ols*zom(i)) &
               -aa-ols*zom(i)
          psih1=-(1.0+2.0*ols/3.0)**1.5-0.667*(ols-aa1) &
               *exp(-0.35*ols)-aa+1.0
          psih0=-(1.0+2.0*ols*zoh(i)/3.0)**1.5-0.667*(ols*zoh(i)-aa1) &
               *exp(-0.35*ols*zoh(i))-aa+1.0

      end if

      ustar=vkc/(-zmln-psim1+psim0)
        tstar=vkc/(-zoln-psih1+psih0)

        cdm(i)=ustar**2.0
        cdh(i)=ustar*tstar/pr

        !z0_momentum(i)=exp(zmln)*zlev
        !z0_heat(i)=exp(zoln)*zlev
        !
        !         * CALCULATE CD*MOD(V) UNDER FREE-CONVECTIVE LIMIT.
        !
        if (tvirtg(i)>tvirta(i)) then
          climit=1.9e-3*(tvirtg(i)-tvirta(i))**0.333333
        else
            climit=0.
        end if
        cflux(i)=max(cdh(i)*wspeed,climit)
      end if
    end do ! loop 100
  end if
  return
end subroutine DRCOEF2
