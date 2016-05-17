!*******************************************************************************
! This subroutines sets the physical constants of the simulation
!*******************************************************************************
subroutine fnPhysConst()
  use comparams
  implicit none
  
  eV=1.6d-19 ![Joules]
  hb=6.5d-16*eV ![eV.s]
  a_cc=1.42d-10 ![meters]
  a_l=sqrt(3.d0)*a_cc
  
  ! these are the numbers used in Zlatan code. I don't know the units of carbon atom mass, but these numbers produce correct phonon energies
  mc=1.d0/sqrt(2.d0) ![?]
  conv_coeff=1.d4 ![dyn/cm]
  
  g0= 5.3*eV/(1.d-10) ![J/meters]
  mass=1.201d1*1.66d-27 ![kg]
  
  kbt=2.5d-2*eV
  
  ! these are the SI units which produce large energies for phonon modes 
  !mc=1.201d1*1.66d-27 ![kg]
  !conv_coeff=1.d1 ![N/m]
  
  return
end