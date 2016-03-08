!*******************************************************************************
! This program calculates the phonon mode energies and eigenfunction of CNTs
! based on graphene phonon modes and using cutting line concept.

! Amirhossein Davoody
! Last modified: 7/24/2014
!*******************************************************************************

program cnt_electron_phonon_interaction
  use cnt_electron_band_structure_mod, only: cnt_electron_band_structure
  use cnt_geometry_mod, only: cnt_geometry
  use cnt_phonon_mod, only: cnt_phonon_dispersion
  use comparams_mod, only: currcnt
  use parse_input_file_mod, only: parse_input_file
  use write_log_mod, only: writeLog
  implicit none

  real :: starttime,endtime !time and date of simulation
	character(len=200) :: logInput

  call CPU_time(starttime)

  call parse_input_file()

  call cnt_geometry(currcnt)

  call cnt_electron_band_structure(currcnt)

  call cnt_phonon_dispersion(currcnt)

  ! call fnInput
  !
  ! call fnOpenFiles
  !
  ! call fnPhysConst
  !
  ! call fnGeomProp
  !
  ! call fnCNTphonon
  !
  ! call fnCNTelectron
  !
  ! call fnScattering
  !
  ! call fnCloseFiles

  call CPU_time(endtime)
  write(logInput,'("Run time = ",f10.3," seconds.")'),endtime-starttime
  call writeLog(trim(logInput))

end program cnt_electron_phonon_interaction
