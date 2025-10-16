subroutine iniclass

  !@Object Initialize CLASS fields
  !@Author K. Winger
  !@Revisions
  ! 001 K. Winger (UQAM/ESCER) Jun 2020 - Initial version
  ! 002 C. Gauthier (Ouranos)  Oct 2020 - Adapt paths for SPS file structure

  use sfc_options  , only : delt
  use sfclayer , only : set_class_const
  use classicParams, only : runParamsFile, prepareGlobalParams, GROWYR

  implicit none
!  character(len=*), intent(in) :: F_path
!  integer, intent(in) :: F_myproc
!  integer :: F_istat
  character(len=255) :: cwd_path, base_path, table_p
  integer :: ierror, pos
  logical, save :: first = .true.

if (first) then
  first = .false.
  print *,'In iniclass.F90'
endif

! Get path to CLASS_input_table following SPS file structure
  call GETCWD(cwd_path, ierror)

  if (ierror == 0) then
    cwd_path = trim(cwd_path)

    pos = index(cwd_path, '/RUNMOD/')

    if (pos > 0) then
      base_path = cwd_path(1:pos + 7)
      table_p = trim(base_path) // 'input/cfg_0000/CLASS_input_table'
      print *, "Base path: ", trim(base_path)
      print *, "New path (table_p): ", trim(table_p)
    end if
  else
    print *, "Error getting current working directory. Error code:", ierror
  end if
  runParamsFile = table_p

  ! Pass some constants from GEM to CLASS
  call set_class_const(delt)

  ! Set pure CLASS/CTEM constants
  call prepareGlobalParams

print *,'GROWYR:',GROWYR(1,1,:)

end subroutine iniclass
