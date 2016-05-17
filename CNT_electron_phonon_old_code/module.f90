!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Declaration of global variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module comparams
  INCLUDE 'link_fnl_static.h'
  implicit none
  
  real*8, parameter :: pi=3.141592d0
  complex*16, parameter :: i1=(0.d0,1.d0)
  
  !Input parameters
  integer :: n_ch,m_ch !chiral vector parameters
  integer :: nkg,nr !reciprocal and real space mesh size in graphene
  integer :: i_sub !the subband number used in exciton energy calculation
  real*8 :: E_th !threshold energy
  real*8 :: Kcm_max
  
  !Physical constants
  real*8 :: eV,hb
  real*8 :: a_cc,a_l !lattice constants
  real*8 :: t0 !tight binding constants
  real*8 :: mc !carbon atom mass used in phonon dispersion calculation
  real*8 :: conv_coeff !the normalization factor for spring constants between carbon atoms
  real*8 :: mass !carbon atom mass used in scattering rate calculation
  real*8 :: g0 !electron-phonon coupling constant
  real*8 :: kbt !thermal energy used in Bose-Einestein distribution
  
  !Geometrical properties
  real*8, dimension(2) :: a1,a2,b1,b2,ch_vec,t_vec,aCC_vec
  real*8 :: len_ch,radius
  integer :: dR,Nu
  integer :: t1,t2
  real*8, dimension(:,:), allocatable :: posA,posB,posAA,posBB,posAB,posBA
  
  !Reciprocal lattice properties
  real*8 :: dk
  real*8, dimension(2) :: K1, K2
  integer :: ikc_max, ikc_min
  
  !CNT band structure properties
  real*8, dimension(:,:,:), allocatable:: Ek
  
  !CNT phonon dispersion
  real*8, dimension (:,:), allocatable:: omega
  
  !These are the unit numbers for the input/output files in the program
  integer :: fh1,fh2,fh3,fh4,fh5,fh6,fh7,fh8,fh9,fh10,fh11,fh12,fh13,fh14,fh15,fh16,fh17,fh18
  
  !mu of the lowest band
  integer :: mu0
  
end module comparams  