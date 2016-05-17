!*******************************************************************************
! This program calculates the phonon mode energies and eigenfunction of CNTs
! based on graphene phonon modes and using cutting line concept.

! Amirhossein Davoody
! Last modified: 7/24/2014
!*******************************************************************************

program CNT_Exciton
  use comparams
  implicit none
  
  call fnInput
  
  call fnOpenFiles
  
  call fnPhysConst
  
  call fnGeomProp
  
  call fnCNTphonon
  
  call fnCNTelectron
  
  call fnScattering
  
  call fnCloseFiles
  
  print *,'Finish!!!!'
  pause
endprogram

