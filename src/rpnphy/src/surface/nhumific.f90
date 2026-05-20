!> \file
!! Calculates N flux from litter pool to soil N pool.
!> @author A. Asaadi, S. Kou-Giesbrecht, J. Melton
!!

subroutine nhumific(il1, il2, sort, fcancmx, spinfast, humtrsvg,& ! In
                    c2nveg_litr, mossPresent, humstepmoss_bulk, c2nMoss_litr, & ! In
                    nlitrmass, soilnmas, litrmsmossN, upMossSoilN, & ! In/Out
                    nhumtrsveg, nhumtrsmoss) ! Out


  use classicParams, only : ilg, ignd, icc, iccp1, iccp2, ctempfts, &
                            deltat, ntolerance

  implicit none

  integer,intent(in) :: il1            !< il1=1
  integer,intent(in) :: il2            !< il2=ilg (no. of grid cells in latitude circle)
  integer,intent(in) :: spinfast       !< spinup factor for soil carbon & nitrogen pools
  integer,intent(in),dimension(icc) :: sort            !< index for correspondence between PFTs and the 12 values in parameters vectors
  real,intent(in),dimension(ilg,icc) :: fcancmx        !< max. fractional coverage of PFTs
  real,intent(in),dimension(ilg,iccp1) :: humtrsvg     !< carbon flux from litter to soil for individual PFTs + bare soil, \f$umol CO2 m^{-2} sec^{-1 }\f$
  real,intent(in),dimension(ilg,iccp1) :: c2nveg_litr  !< simulated C:N ratio for litter, \f$g C g N^{-1}\f$
  character(8), intent(in) :: mossPresent(ilg)         !< Type of moss present (else 'None')
  real,intent(in), dimension(ilg) :: c2nMoss_litr      !< Moss simulated C:N ratio for detritus, \f$g C g N^{-1}\f$
  real,intent(in), dimension(ilg) :: humstepmoss_bulk  !< moss humification (transfer from litter to SOM) (kgC/m2/timestep)

  real,intent(out),dimension(ilg,iccp2) :: nhumtrsveg  !< nitrogen flux from litter to soil for individual PFTs + bare soil + LUC product, \f$g N m^{-2} day^{-1}\f$ 
  real,intent(out), dimension(ilg)      :: nhumtrsmoss !< nitrogen flux from litter to soil for moss, \f$g N m^{-2} day^{-1}\f$ 
  
  real,intent(inout),dimension(ilg,iccp2) :: nlitrmass   !< litter nitrogen for individual PFTs + bare soil, \f$g N m^{-2}\f$
  real,intent(inout),dimension(ilg,iccp2) :: soilnmas    !< soil nitrogen for individual PFTs + bare soil, \f$g N m^{-2}\f$
  real,intent(inout), dimension(ilg)      :: litrmsmossN !< moss litter N pool, \f$g N/m^2\f$
  real,intent(inout), dimension(ilg)      :: upMossSoilN !< Moss soil N mass (non-peat), \f$g N/m^2\f$

  ! local variables
  integer :: i, j
  real :: temp
  character(8) :: pftkind
  real,dimension(ilg,iccp1) :: nlitrmass_in !< previous litter nitrogen for individual PFTs + bare soil, \f$g N m^{-2}\f$
  real,dimension(ilg,iccp1) :: soilnmas_in  !< previous soil nitrogen for individual PFTs + bare soil, \f$g N m^{-2}\f$
  real, dimension(ilg)      :: litrmsmossN_in !< previous moss litter N pool, \f$g N/m^2\f$
  real, dimension(ilg)      :: upMossSoilN_in !< previous Moss soil N mass (non-peat), \f$g N/m^2\f$

 temp = 0.0

  do i = il1,il2
     do j = 1,iccp1
      nhumtrsveg(i,j) = 0.0
      nlitrmass_in(i,j) = 0.0
      soilnmas_in(i,j) = 0.0
     end do ! loop 110
     nhumtrsmoss(i) = 0.
  end do ! loop 100

  !! Calculate N flux from litter pool to soil N pool
  !! (see equation 18 in Asaadi and Arora (2021); https://doi.org/10.5194/bg-2020-147)
  !! \f$ H_{N,D2H} = \frac{H_{C,D2H}}{C:N_D} \f$
  do j = 1,icc
     do i = il1,il2

      if (nlitrmass(i,j) > 0.0 .and. fcancmx(i,j) > 0.0 .and. &
                                   c2nveg_litr(i,j) > 0.0) then
         nhumtrsveg(i,j) = humtrsvg(i,j) * (deltat/963.62) &
                           * 1000. / c2nveg_litr(i,j)
         nhumtrsveg(i,j) = max(0.0, min(nhumtrsveg(i,j) , nlitrmass(i,j)))
      endif

     end do ! loop 130
  end do ! loop 120

  !! bare ground
  do i = il1,il2
   if (nlitrmass(i,iccp1) > 0.0 .and. c2nveg_litr(i,iccp1) > 0.0) then
      nhumtrsveg(i,iccp1) = humtrsvg(i,iccp1)  * (deltat/963.62) &
                            * 1000. / c2nveg_litr(i,iccp1)
      nhumtrsveg(i,iccp1) = max(0.0, min(nhumtrsveg(i,iccp1) , nlitrmass(i,iccp1)))
   endif
  end do ! loop 135

  !! Moss, if present
  do i = il1,il2
   if (mossPresent(i) /= 'None') then 
      if (litrmsmossN(i) > 0.0 .and. c2nMoss_litr(i) > 0.0) then
         nhumtrsmoss(i) = humstepmoss_bulk(i) * 1000. / c2nMoss_litr(i)
         nhumtrsmoss(i) = max(0.0, min(nhumtrsmoss(i) , litrmsmossN(i)))
      end if 
   end if
  end do ! loop 135

  !> Updating litter and soil nitrogen and checking conservation 
  do j = 1,iccp1
     do i = il1,il2

      nlitrmass_in(i,j) = nlitrmass(i,j)
      soilnmas_in(i,j) = soilnmas(i,j)

      !! 1) nlitrmass

      nlitrmass(i,j)= nlitrmass(i,j) - nhumtrsveg(i,j)

      !! Adjusting to prevent negative nlitrmass(i,j)
      if (nlitrmass(i,j) < 0.0) then
         write(*,*)'nlitrmass lt zero in nhumific.f90 at i=',i,' for pft=',j
         ! write(*,*)'nlitrmass(i,j)  = ',nlitrmass(i,j)
         ! write(*,*)'nhumtrsveg(i,j) = ',nhumtrsveg(i,j)
         write(*,*)'adjusting to prevent negative nlitrmass(i,j)'
         temp = nhumtrsveg(i,j)
         nhumtrsveg(i,j) = nhumtrsveg(i,j) - nhumtrsveg(i,j) &
                           * abs(nlitrmass(i,j) / temp)
         nlitrmass(i,j) = 0.0
         temp = 0.0
      end if

      !! Checking conservation for nlitrmass(i,j)
      if ((nlitrmass_in(i,j) - nlitrmass(i,j) - nhumtrsveg(i,j)) &
                            > ntolerance) then
         write(*,*)'imbalance in nlitrmass pool at (i)=',i,'pft=',j
         call errorHandler('nhumific', - 1)
      end if


      !! 2) soilnmas

      soilnmas(i,j) = soilnmas(i,j) + real(spinfast) * nhumtrsveg(i,j)

      if (soilnmas(i,j) < 0.0) then
         write(*,*)'soilnmas lt zero in nhumific.f90 at i=',i,' for pft=',j
         ! write(*,*)'soilnmas(i,j)  = ',soilnmas(i,j)
         ! write(*,*)'nhumtrsveg(i,j) = ',nhumtrsveg(i,j)
         call errorHandler('nhumific', - 2)
      end if

      ! Checking conservation for soilnmas(i,j)
      if ((soilnmas_in(i,j) - soilnmas(i,j) + nhumtrsveg(i,j)) &
                         > ntolerance .and. spinfast == 1) then
         write(*,*)'imbalance in soilnmas pool at (i)=',i,'pft=',j
         call errorHandler('nhumific', - 3)
      end if

    end do ! loop 150
  end do ! loop 140

  !> Update the moss litter and soil N pools and check conservation
  do i = il1, il2
    if (mossPresent(i) /= 'None') then

      litrmsmossN_in(i) = litrmsmossN(i)
      upMossSoilN_in(i) = upMossSoilN(i)

      ! Update moss litter 
      litrmsmossN(i) = litrmsmossN(i) - nhumtrsmoss(i)

      ! make sure moss litter doesn't go negative
      if (litrmsmossN(i) < 0.) then 
         write(*,*)'litrmsmossN negative in nhumific.f90 at i=',i
         write(*,*)'adjusting to prevent negative litrmsmossN(i)'
         temp = nhumtrsmoss(i)
         nhumtrsmoss(i) = nhumtrsmoss(i) - nhumtrsmoss(i) &
                           * abs(litrmsmossN(i) / temp)
         litrmsmossN(i) = 0.0
      end if 

      !! Checking conservation for litrmsmossN(i)
      if ((litrmsmossN_in(i) - litrmsmossN(i) - nhumtrsmoss(i)) &
                            > ntolerance) then
         write(*,*)'imbalance in litrmsmossN pool at (i)=',i
         call errorHandler('nhumific', -4)
      end if

      ! Update moss soil N mass
      upMossSoilN(i) = upMossSoilN(i) + real(spinfast) * nhumtrsmoss(i)

      ! make sure moss soil N doesn't go negative
      if (upMossSoilN(i) < 0.) then 
         write(*,*)'upMossSoilN negative in nhumific.f90 at i=',i
         call errorHandler('nhumific', - 5)
      end if 

      ! Check conservation of upMossSoilN
      if ((upMossSoilN_in(i) - upMossSoilN(i) + nhumtrsmoss(i)) &
                         > ntolerance .and. spinfast == 1) then
         write(*,*)'imbalance in upMossSoilN pool at (i)=',i
         call errorHandler('nhumific', - 6)
      end if

    end if 
  end do 

  return
end subroutine nhumific
