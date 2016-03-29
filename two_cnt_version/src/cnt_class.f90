module cnt_class
	implicit none
	private

	public  :: cnt, free_cnt_memory

	type cnt
		integer :: n_ch,m_ch !chiral vector parameters
		integer :: i_sub !subband index used in exciton energy calculation

		!Geometrical properties
		real*8, dimension(2) :: a1,a2,b1,b2,ch_vec,t_vec,aCC_vec
		real*8 :: len_ch,radius
		integer :: Nu !number of graphene unit cells in cnt unit cell.
		integer :: nr !length of cnt in terms of its unit cell.
		real*8, dimension(:,:), allocatable :: posA,posB,posAA,posBB,posAB,posBA
		real*8, dimension(:,:), allocatable :: posA3, posB3
		real*8, dimension(:,:,:), allocatable :: pos2d, pos3d
		real*8, dimension(:,:), allocatable :: r_posA3, ur_posA3 ! this is the rotated and unrotated position of carbon atoms in 3D.
		real*8, dimension(:), allocatable :: az_angle ! this is azimuthal angle of carbon atoms in roled CNT

		!Length and location of cnt for calculating the resonance energy transfer rate
		real*8 :: Length
		real*8 :: center_position

		!Environment properties
		real*8 :: kappa !this is the dielectric factor that accounts for core electron and environment screening
		real*8 :: Ckappa !this is the factor that is multiplied to a kappa_coeff to yield kappa. This is varient under different environments.
		real*8 :: kappa_coeff !this is the scaling factor for calculating kappa and is different for each different cnt chirality.

		!Reciprocal lattice properties
		integer :: nkg
		real*8 :: dk !reciprocal lattice mesh size for calculating self-energy, Fourier transform of coulomb interaction v_FT, and dielectric function.
		real*8 :: dkx !reciprocal lattice mesh size for calculating exciton dispersion. This is used to determine the spacing between calculated exciton spacing: K_cm = iKcm * dkx
		integer :: dk_dkx_ratio ! this is the ratio of dk and dkx: dk = dkx * dk_dkx_ratio
		real*8, dimension(2) :: K1, K2

		!CNT band structure properties
		integer, dimension(:), allocatable :: min_sub
		integer :: ikc_max, ikc_min !these are index limits for the wave vector inside the carbon nanotube brillouine zone in the direction of the carbon nanotube axis.
		integer :: ik_max, ik_min, iKcm_max, iKcm_min, ik_high, ik_low, ikr_high, ikr_low, iq_max, iq_min
		integer :: iKcm_min_fine, iKcm_max_fine
		integer :: mu_cm

		!Dielectric function
		! quantities that have _fine at the end are those that are calculated via interpolation of the quantities without _fine at the end.
		real*8, dimension(:,:), allocatable :: eps_q, eps_q_fine
		complex*16, dimension(:,:,:,:), allocatable :: v_FT, v_FT_fine ! v_FT(mu,q,n,m) stores the Fourier transform of the Coulomb potential at the wavevector determined by band index "mu" and wavenumber "q" for atoms of type A (n or m = 1) or type B (n or m = 2)

		!CNT self energy and tight binding coefficients
		! quantities that have _fine at the end are those that are calculated via interpolation of the quantities without _fine at the end.
		real*8, dimension(:,:,:), allocatable :: Ek, Ek_fine!Tight-binding energy , Ek(mu,k,n) stores the tight-binding energy of the band "mu" with wavevector "k" in conduction band (n=1) or the valence band (n=2).
		real*8, dimension(:,:,:), allocatable :: Sk, Sk_fine!Self-energy
		complex*16, dimension(:,:,:), allocatable :: Cc, Cv, Cc_fine, Cv_fine !Cc(mu,k,b) is the conduction band tight-binding wavefunction coefficients where "mu" is the band index (1 is +mu and 2 is -mu), "k" is the wave vector along the CNT axis, "b" is the atom index in graphene unit cell (1 is A type atom) and (2 is B type atom)

		!A-type exciton wavefunction and energies
		real*8, dimension(:,:), allocatable :: Ex_A1, Ex0_A2, Ex1_A2 !the first index is subband, the second index is iKcm
		complex*16, dimension(:,:,:), allocatable :: Psi_A1, Psi0_A2, Psi1_A2 !the first index is ikr, the scond index is the subband, the third index is iKcm

		!E-type exciton wavefunction and energies
		real*8, dimension(:,:), allocatable :: Ex0_Ep, Ex0_Em, Ex1_Ep, Ex1_Em !the first index is subband, the second index is iKcm
		complex*16, dimension(:,:,:), allocatable :: Psi0_Ep, Psi0_Em, Psi1_Ep, Psi1_Em !the first index is ikr, the scond index is the subband, the third index is iKcm

		!Target exciton wavefunction and energies
		real*8, dimension(:,:), allocatable :: Ex_t !the first index is subband, the second index is iKcm
		complex*16, dimension(:,:,:), allocatable :: Psi_t !the first index is ikr, the scond index is the subband, the third index is iKcm
		character(len=20) :: targetExcitonType !this is the type of target exciton which should be one this options: Ex_A1, Ex0_A2, Ex1_A2
		real*8 :: ex_symmetry

		!number of exciton bands below free-electron free-hole energy level
		integer :: nX_a, nX_e, nX_t
		real*8 :: E_th
		real*8 :: Kcm_max

		!cnt phonon dispersion
		real*8, dimension (:,:,:), allocatable :: omega_phonon

		!directory that the CNT information is stored
		character(len=1000) :: directory

	end type cnt

contains

	!**************************************************************************************************************************
	! subroutine to free all allocatable quantities in cnt_class
	!**************************************************************************************************************************

	subroutine free_cnt_memory(currcnt)

		type(cnt), intent(inout) :: currcnt

		if (allocated(currcnt%posA)) deallocate(currcnt%posA)
		if (allocated(currcnt%posB)) deallocate(currcnt%posB)
		if (allocated(currcnt%posAA)) deallocate(currcnt%posAA)
		if (allocated(currcnt%posBB)) deallocate(currcnt%posBB)
		if (allocated(currcnt%posAB)) deallocate(currcnt%posAB)
		if (allocated(currcnt%posBA)) deallocate(currcnt%posBA)
		if (allocated(currcnt%posA3)) deallocate(currcnt%posA3)
		if (allocated(currcnt%posB3)) deallocate(currcnt%posB3)
		if (allocated(currcnt%pos2d)) deallocate(currcnt%pos2d)
		if (allocated(currcnt%pos3d)) deallocate(currcnt%pos3d)
		if (allocated(currcnt%r_posA3)) deallocate(currcnt%r_posA3)
		if (allocated(currcnt%ur_posA3)) deallocate(currcnt%ur_posA3)
		if (allocated(currcnt%az_angle)) deallocate(currcnt%az_angle)
		if (allocated(currcnt%min_sub)) deallocate(currcnt%min_sub)
		if (allocated(currcnt%eps_q)) deallocate(currcnt%eps_q)
		if (allocated(currcnt%eps_q_fine)) deallocate(currcnt%eps_q_fine)
		if (allocated(currcnt%v_FT)) deallocate(currcnt%v_FT)
		if (allocated(currcnt%v_FT_fine)) deallocate(currcnt%v_FT_fine)
		if (allocated(currcnt%Ek)) deallocate(currcnt%Ek)
		if (allocated(currcnt%Ek_fine)) deallocate(currcnt%Ek_fine)
		if (allocated(currcnt%Sk)) deallocate(currcnt%Sk)
		if (allocated(currcnt%Sk_fine)) deallocate(currcnt%Sk_fine)
		if (allocated(currcnt%Cc)) deallocate(currcnt%Cc)
		if (allocated(currcnt%Cc_fine)) deallocate(currcnt%Cc_fine)
		if (allocated(currcnt%Cv)) deallocate(currcnt%Cv)
		if (allocated(currcnt%Cv_fine)) deallocate(currcnt%Cv_fine)
		if (allocated(currcnt%Ex_A1)) deallocate(currcnt%Ex_A1)
		if (allocated(currcnt%Ex0_A2)) deallocate(currcnt%Ex0_A2)
		if (allocated(currcnt%Ex1_A2)) deallocate(currcnt%Ex1_A2)
		if (allocated(currcnt%Psi_A1)) deallocate(currcnt%Psi_A1)
		if (allocated(currcnt%Psi0_A2)) deallocate(currcnt%Psi0_A2)
		if (allocated(currcnt%Psi1_A2)) deallocate(currcnt%Psi1_A2)
		if (allocated(currcnt%Ex0_Em)) deallocate(currcnt%Ex0_Em)
		if (allocated(currcnt%Ex0_Ep)) deallocate(currcnt%Ex0_Ep)
		if (allocated(currcnt%Ex1_Em)) deallocate(currcnt%Ex1_Em)
		if (allocated(currcnt%Ex1_Ep)) deallocate(currcnt%Ex1_Ep)
		if (allocated(currcnt%Psi0_Em)) deallocate(currcnt%Psi0_Em)
		if (allocated(currcnt%Psi0_Ep)) deallocate(currcnt%Psi0_Ep)
		if (allocated(currcnt%Psi1_Em)) deallocate(currcnt%Psi1_Em)
		if (allocated(currcnt%Psi1_Ep)) deallocate(currcnt%Psi1_Ep)
		if (allocated(currcnt%Ex_t)) deallocate(currcnt%Ex_t)
		if (allocated(currcnt%Psi_t)) deallocate(currcnt%Psi_t)
		if (allocated(currcnt%omega_phonon)) deallocate(currcnt%omega_phonon)

	end subroutine free_cnt_memory

end module cnt_class
