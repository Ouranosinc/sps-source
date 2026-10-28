!> \fireenergycl.F90
!> \This subroutine uses burned area and fuel consumed per area to calculate fire energy, 
!>  which will be used in CanAM/atmosphere for fire plume rise. Same as fireenergy.F90 in ctem/ 
!>  but modified for use with CLASSIC.
!!
!! @author Cynthia Whaley
!
! This subroutine is called by disturbance.f90 when cffeps_plume=.true.

subroutine fireenergycl (tot_emit_dom,burnarea,fiengat,icc,num_fires)  ! ,burnveg
  
  implicit none

  integer, intent(in) :: icc   !< Index of 9 plant functional types (1-9) \f$[unitless]\f$
  real, intent(in)    :: burnarea	!< area burned per PFT \f$[km^2]\f$
  real, intent(in)    :: tot_emit_dom	!< fuel (dry organic matter, dom) consumed per area per PFT \f$[kg/m^2]\f$
  real, intent(in)    :: num_fires      !< number of fires per PFT
  real, intent(out)   :: fiengat	!< fire energy that goes into the plume per PFT \f$[J]\f$
  
  ! internal work fields
  integer :: j
 ! real, dimension(ilg,icc) :: pftscaler !< scale factor for fire energy depending on PFT in the grid cell \f$[unitless]\f$
 ! real, dimension(ilg,icc) :: kgperpft !< PFT consumed per area \f$[kg/m^2]\f$

  real :: toten !< total fire energy \f$[joules]\f$
  real :: A     !< burned area in \f$[m^2]\f$

  ! Coefficients for the energy equation:
  real, parameter :: H = 18000000. !< heat of combustion = 1.8x10^7 J/kg for dry wood\f$[J/kg]\f$
 ! real, parameter, dimension(9) :: hocs = (/1.0, 1.0, 1.0, 1.0, 1.0, 0.01, 0.01, 1.0, 1.0/) !< scale factor for heat of combustion of the 9 PFTs
  ! E.g., if most of the grid cell is crops, fire energy will be scaled down considerably.
  !
  ! initialize
  toten = 0.0
  A = 0.0

! removed PFT-dependence on fire energy for now.
!  do j = 1, icc 
!     kgperpft(j) = 0.0
!     pftscaler(j) = 0.0
!  end do
    
!  ! Determine appropriate scale factor based on amount of PFT burned in the grid cell and its heat of combustion
!  do j = 1, icc 
!      if(burnarea>0.0) then
!        kgperpft(j) = tot_emit_dom * burnveg(j)/burnarea
!        pftscaler(j) = kgperpft(j) * hocs(j)
!      end if
!  end do
  

  ! Fire energy calculation
   if(burnarea>0.0 .and. tot_emit_dom>0.0) then
     A = 1000.0 * 1000.0 * burnarea     ! convert from km^2 to m^2 
      ! test to try out differences in PFTs:
      ! do j = 1, icc  
      !	  toten = toten + (H * pftscaler(j) * A)
      ! end do
     toten = H * tot_emit_dom * A
     if(num_fires > 1.) then
        toten = toten / num_fires   ! avg fire energy per fire
     end if
     fiengat = 0.2 * toten       ! [CW] see note below.
   else
     fiengat = 0.0
   end if

  
  ! Note: not all fire energy goes into the plume. There are 5 other 'compartments' for the fire energy. See
  ! mach_cffeps_energy.F90 on how to calculate the plume energy specifically. For now, just to get the 'plumbing'
  ! working, a simplification of 0.2*toten is used, as it's dry fuel consumed, but energy calculation expects wet. 
  ! This is noted in the CFFEPS report to be a reasonable value (they say 10-20%). 
  ! A scale factor based on PFT and the PFT area burned is also included, since dry wood "H" is used for all. But this
  ! is commented out for now.  
  
  
  return
end subroutine fireenergycl
!> \file
!> Energy Balance of a Fire required for Plume Rise
!!     Anderson, K.R.; Pankratz, A; Mooney, C. 2011.   
!!     A thermodynamic approach to estimating smoke plume heights.
!!     In 9th Symp. on Fire and Forest Meteorology, 
!!     Oct 18-20, 2011.  Palm Springs, CA.
!!     Am. Meteorol. Soc., Boston, MS.
!! Note: 9 PFTs in CTEM are: needleleaf evergreen trees, needleleaf deciduous trees, 
!!       broadleaf evergreen trees, broadleaf cold deciduous trees, broadleaf drought/dry 
!!       deciduous trees, C3 crop, C4 crop, C3 grass, and C4 grass.
