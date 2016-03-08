!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Declaration of physical variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module physical_constant_mod
  implicit none
  private

  !Mathematical constants
  real*8, public, parameter :: pi=3.141592d0
  complex*16, public, parameter :: i1=(0.d0,1.d0)

  !Physical constants
  real*8, public, parameter :: eV=1.6d-19 ![Joules]
  real*8, public, parameter :: hb=6.5d-16*eV ![eV.s]
  real*8, public, parameter :: kb=1.3865d-23 ![J/K]

  real*8, public, parameter :: a_cc=1.42d-10 !carbon-carbon distance [meters]
  real*8, public, parameter :: a_l=dsqrt(3.d0)*a_cc !lattice constants

  real*8, public, parameter :: e2p = 0.d0 !tight binding constants
  real*8, public, parameter :: t0 = 2.7d0*eV  !tight binding constants
  real*8, public, parameter :: s0 = 0.d0 !tight binding constants

  real*8, public, parameter :: Upp = 11.3d0*eV !constant used in the Ohno potential
  real*8, public, parameter :: eps0 = 8.85d-12 !permittivity of free space
  real*8, public, parameter :: q0 = 1.6d-19 !charge of electron

  ! these are the numbers used in Zlatan code for calculating phonon dispersion. I don't know the units of carbon atom mass, but these numbers produce correct phonon energies
  real*8, public, parameter :: m_carbon_dispersion = 1.d0/sqrt(2.d0) ![?] !carbon atom mass used in phonon dispersion calculation
  real*8, public, parameter :: spring_conv_coeff = 1.d4 ![dyn/cm] !the normalization factor for spring constants between carbon atoms

  ! these are the SI units which produce unrealistically large energies for phonon dispersion curves
  !m_carbon_dispersion=1.201d1*1.66d-27 ![kg]
  !conv_coeff=1.d1 ![N/m]

  ! the following "mass" parameter is the carbon atom masses used for calculating electron-phonon scattering rates.
  real*8, public, parameter :: m_carbon = 1.201d1*1.66d-27 ![kg] !carbon atom mass used in scattering rate calculation
  real*8, public, parameter :: g0 = 5.3*eV/(1.d-10) ![J/meters] !electron-phonon coupling constant

end module physical_constant_mod
