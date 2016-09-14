module cnt_class
	implicit none
	private

	public  :: cnt, exciton, free_cnt_memory, free_exciton_memory

	!***************************************************************************
	!	-definition of exciton class (derived type) that is used in the cnt
	!	 derived type.
	!***************************************************************************
	type exciton
		character(len=1000) :: name
		integer :: i_sub
		integer :: spin
		integer :: mu_cm
		integer :: iKcm_max, iKcm_min, ikr_high, ikr_low
		integer :: n_mu_r ! this will hold the size of mu_r
		integer :: nx ! this is the number of exciton state that is calculated at each iKcm
		integer :: dkr_dKcm_ratio
		real*8 :: dKcm, dkr
		integer, dimension(:), allocatable :: mu_r
		real*8, dimension(:,:), allocatable :: ex !exciton energy: the first index is subband, the second index is iKcm
		complex*16, dimension(:,:,:,:), allocatable :: psi !exciton wavefunction in k-space: the first index is ikr, the scond index is the subband, the third index is iKcm, the fourth index is mu_r
	end type exciton

	!***************************************************************************
	!	-definition of carbon nanotube class (derived type)
	!***************************************************************************
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

		! ik_max			:the higher limit of k-vector that is below E_th
		! ik_min			:the lower limit of k-vector that is below E_th
		! iKcm_max			:the higher limit of center of mass wave vector that we calculate
		! iKcm_min			:the lower limit of center of mass wave vector that we calculate
		! ikr_high			:the maximum index that the relative wavenumber in the entire simulation.
		! ikr_low			:the minimum index that the relative wavenumber in the entire simulation.
		! ik_high			:the maximum index that the wavenumber can get in the entire simulation.
		! ik_low			:the minimum index that the wavenumber can get in the entire simulation.
		! iq_max			:the higher limit of the index in v_FT and esp_q
		! iq_min			:the lower limit of the index in v_FT and esp_q
		! iKcm_max_fine		:the upper limit of center of mass wave vector that we calculate when using a finer mesh size for exciton center of mass momentum
		! iKcm_min_fine		:the lower limit of center of mass wave vector that we calculate when using a finer mesh size for exciton center of mass momentum

		!Dielectric function
		! quantities that have _fine at the end are those that are calculated via interpolation of the quantities without _fine at the end.
		real*8, dimension(:,:), allocatable :: eps_q, eps_q_fine
		complex*16, dimension(:,:,:,:), allocatable :: v_FT, v_FT_fine ! v_FT(mu,q,n,m) stores the Fourier transform of the Coulomb potential at the wavevector determined by band index "mu" and wavenumber "q" for atoms of type A (n or m = 1) or type B (n or m = 2)

		!CNT self energy and tight binding coefficients
		! quantities that have _fine at the end are those that are calculated via interpolation of the quantities without _fine at the end.
		real*8, dimension(:,:,:), allocatable :: Ek, Ek_fine!Tight-binding energy , Ek(mu,k,n) stores the tight-binding energy of the band "mu" with wavevector "k" in conduction band (n=1) or the valence band (n=2).
		real*8, dimension(:,:,:), allocatable :: Sk, Sk_fine!Self-energy
		complex*16, dimension(:,:,:), allocatable :: Cc, Cv, Cc_fine, Cv_fine !Cc(mu,k,b) is the conduction band tight-binding wavefunction coefficients where "mu" is the band index (1 is +mu and 2 is -mu), "k" is the wave vector along the CNT axis, "b" is the atom index in graphene unit cell (1 is A type atom) and (2 is B type atom)

		! -	excitons(ex_type,alpha) is the array of all excitons. the first
		!	index (ex_type)	represents symmetry or center of mass. The second
		!	index (alpha) is the singlet or triplet. specifically we take:
		!	ex_type=1 : A1 exciton.
		!	ex_type=2 : A2 exciton.
		!	ex_type=3 : Ep exciton.
		!	ex_type=4 : Em exciton.
		!	alpha=0 : singlet.
		!	alpha=1 : triplet.
		type(exciton), dimension(1:4,0:1) :: excitons

		type(exciton), pointer :: selected_exciton

		character(len=1000) :: selected_exciton_name !this is the type of target exciton which should be one this options: Ex_A1, Ex0_A2, Ex1_A2

		real*8 :: E_th
		real*8 :: Kcm_max

		!cnt phonon dispersion
		real*8, dimension (:,:,:), allocatable :: omega_phonon

		!directory that the CNT information is stored
		character(len=1000) :: directory

		!name of the CNT for writing the simulation resutls in the output directory
		character(len=1000) :: name

	end type cnt


	!***************************************************************************
	!	- carbon nanotube objects that are used in the simulation
	!***************************************************************************
	type (cnt), target, public :: cnt1, cnt2

contains

	!***************************************************************************
	! subroutine to free all allocatable quantities in cnt_class
	!***************************************************************************

	subroutine free_cnt_memory(currcnt)

		type(cnt), intent(inout) :: currcnt
		integer :: exciton_type, alpha

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
		if (allocated(currcnt%omega_phonon)) deallocate(currcnt%omega_phonon)

		do alpha = 0,1
			do exciton_type = 1,4
				if(allocated(currcnt%excitons(exciton_type, alpha)%mu_r)) deallocate(currcnt%excitons(exciton_type, alpha)%mu_r)
				if(allocated(currcnt%excitons(exciton_type, alpha)%ex)) deallocate(currcnt%excitons(exciton_type, alpha)%ex)
				if(allocated(currcnt%excitons(exciton_type, alpha)%psi)) deallocate(currcnt%excitons(exciton_type, alpha)%psi)
			enddo
		enddo

	end subroutine free_cnt_memory

	subroutine free_exciton_memory(my_exciton)
		type(exciton), intent(inout) :: my_exciton

		if(allocated(my_exciton%mu_r)) deallocate(my_exciton%mu_r)
		if(allocated(my_exciton%ex)) deallocate(my_exciton%ex)
		if(allocated(my_exciton%psi)) deallocate(my_exciton%psi)

	end subroutine

end module cnt_class
