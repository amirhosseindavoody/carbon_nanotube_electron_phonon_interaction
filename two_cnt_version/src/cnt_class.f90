module cnt_class
	implicit none
	private

	public  :: cnt, cnt_geometry, cnt_band, free_cnt_memory

	type cnt
		integer, public :: n_ch,m_ch !chiral vector parameters
		integer, public :: i_sub !subband index used in exciton energy calculation

		!Geometrical properties
		real*8, dimension(2) :: a1,a2,b1,b2,ch_vec,t_vec
		real*8, dimension(2), public :: aCC_vec
		real*8 :: len_ch,radius
		integer :: Nu
		integer, public :: nr
		real*8, dimension(:,:), allocatable, public :: posA,posB,posAA,posBB,posAB,posBA
		real*8, dimension(:,:), allocatable, public :: posA3, posB3
		real*8, dimension(:,:,:), allocatable, public :: pos2d, pos3d
		real*8, dimension(:,:), allocatable, public :: r_posA3, ur_posA3 ! this is the rotated and unrotated position of carbon atoms in 3D.
		real*8, dimension(:), allocatable, public :: az_angle ! this is azimuthal angle of carbon atoms in roled CNT

		!Length and location of cnt for calculating the resonance energy transfer rate
		real*8, public :: Length
		real*8, public :: center_position

		!Environment properties
		real*8, public :: kappa, Ckappa

		!Reciprocal lattice properties
		integer, public :: nkg
		integer, public :: dk_dkx_ratio
		real*8, public :: dk, dkx
		real*8, dimension(2) :: K1, K2

		!CNT band structure properties
		integer, dimension(:), allocatable :: min_sub
		integer, public :: ikc_max, ikc_min, ik_max, ik_min, iKcm_max, iKcm_min, ik_high, ik_low, ikr_high, ikr_low, iq_max, iq_min
		integer, public :: iKcm_min_fine, iKcm_max_fine
		integer, public :: mu_cm

		!CNT self energy and tight binding coefficients
		real*8, dimension(:,:,:), allocatable, public :: Ek  !Tight-binding energy
		real*8, dimension(:,:,:), allocatable:: Sk  !Self-energy
		complex*16, dimension(:,:,:), allocatable:: Cc,Cv !Tight-binding wavefunction coefficients

		!Exciton wavefunction and energies
		real*8, dimension(:,:), allocatable, public :: Ex_A1, Ex0_A2, Ex1_A2 !the first index is subband, the second index is iKcm
		complex*16, dimension(:,:,:), allocatable, public :: Psi_A1, Psi0_A2, Psi1_A2 !the first index is ikr, the scond index is the subband, the third index is iKcm
		real*8, dimension(:,:), allocatable, public :: Ex0_Em, Ex0_Ep, Ex1_Em, Ex1_Ep !the first index is subband, the second index is iKcm
		complex*16, dimension(:,:,:), allocatable, public :: Psi0_Em, Psi0_Ep, Psi1_Em, Psi1_Ep !the first index is ikr, the scond index is the subband, the third index is iKcm

		!Target exciton wavefunction and energies
		real*8, dimension(:,:), allocatable, public :: Ex_t !the first index is subband, the second index is iKcm
		complex*16, dimension(:,:,:), allocatable, public :: Psi_t !the first index is ikr, the scond index is the subband, the third index is iKcm
		character(len=20), public :: targetExcitonType !this is the type of target exciton which should be one this options: Ex_A1, Ex0_A2, Ex1_A2
		real*8, public :: ex_symmetry

		!number of exciton bands below free-electron free-hole energy level
		integer, public :: nX_a, nX_e, nX_t
		real*8, public :: E_th
		real*8, public :: Kcm_max

		!directory that the CNT information is stored
		character(len=200), public :: directory

	end type cnt

contains

		!**************************************************************************************************************************
		! initialize CNT by calculating its geometrical properties
		!**************************************************************************************************************************

		subroutine cnt_geometry(currcnt)
			use math_Functions_mod, only: gcd, my_norm2
			use physicalConstants, only: a_l, pi

			type(cnt), intent(inout) :: currcnt

			integer :: dR = 1
			integer :: t1,t2
			real*8 :: cosTh, sinTh
			real*8, dimension(2,2) :: Rot
			integer :: i,j,k

			! unit vectors and reciprocal lattice vectors.
			currcnt%a1=(/dsqrt(3.d0)/2.d0*a_l, 1.d0/2.d0*a_l/)
			currcnt%a2=(/dsqrt(3.d0)/2.d0*a_l, -1.d0/2.d0*a_l/)
			currcnt%b1=(/1.d0/dsqrt(3.d0)*2.d0*pi/a_l, +1.d0*2.d0*pi/a_l/)
			currcnt%b2=(/1.d0/dsqrt(3.d0)*2.d0*pi/a_l, -1.d0*2.d0*pi/a_l/)
			currcnt%aCC_vec=1.d0/3.d0*(currcnt%a1+currcnt%a2)

			! calculate chirality and translational vectors of CNT unit cell.
			currcnt%ch_vec=dble(currcnt%n_ch)*currcnt%a1+dble(currcnt%m_ch)*currcnt%a2
			currcnt%len_ch=a_l*dsqrt(dble(currcnt%n_ch)**2+dble(currcnt%m_ch)**2+dble(currcnt%n_ch)*dble(currcnt%m_ch))
			currcnt%radius=currcnt%len_ch/2.d0/pi

			call gcd(dR,2*currcnt%n_ch+currcnt%m_ch,2*currcnt%m_ch+currcnt%n_ch)

			t1=(2*currcnt%m_ch+currcnt%n_ch)/dR
			t2=-(2*currcnt%n_ch+currcnt%m_ch)/dR

			currcnt%t_vec=dble(t1)*currcnt%a1+dble(t2)*currcnt%a2

			currcnt%Nu=2*(currcnt%n_ch**2+currcnt%m_ch**2+currcnt%n_ch*currcnt%m_ch)/dR


			! rotate basis vectors so that ch_vec is along x-axis.
			cosTh=currcnt%ch_vec(1)/my_norm2(currcnt%ch_vec)
			sinTh=currcnt%ch_vec(2)/my_norm2(currcnt%ch_vec)
			Rot=reshape((/ cosTh, -sinTh , sinTh, cosTh /), (/2,2/))
			currcnt%ch_vec=matmul(Rot,currcnt%ch_vec)
			currcnt%t_vec=matmul(Rot,currcnt%t_vec)
			currcnt%a1=matmul(Rot,currcnt%a1)
			currcnt%a2=matmul(Rot,currcnt%a2)
			currcnt%b1=matmul(Rot,currcnt%b1)
			currcnt%b2=matmul(Rot,currcnt%b2)
			currcnt%aCC_vec=matmul(Rot,currcnt%aCC_vec)

			! calculate reciprocal lattice of CNT.
			currcnt%dk=my_norm2(currcnt%b1)/(dble(currcnt%nkg)-1.d0)
			currcnt%dkx=currcnt%dk / currcnt%dk_dkx_ratio
			currcnt%K1=(-t2*currcnt%b1+t1*currcnt%b2)/(dble(currcnt%Nu))
			currcnt%K2=(dble(currcnt%m_ch)*currcnt%b1-dble(currcnt%n_ch)*currcnt%b2)/dble(currcnt%Nu)
			currcnt%K2=currcnt%K2/my_norm2(currcnt%K2)
			currcnt%ikc_max=floor(pi/my_norm2(currcnt%t_vec)/currcnt%dk)
			currcnt%ikc_min=-currcnt%ikc_max

			! calculate coordinates of atoms in the unwarped CNT unit cell.
			allocate(currcnt%posA(currcnt%Nu,2))
			allocate(currcnt%posB(currcnt%Nu,2))

			k=0
			do i=0,t1+currcnt%n_ch
				do j=t2,currcnt%m_ch
				if ((dble(t2)/dble(t1)*dble(i) .le. dble(j)) .and. (dble(currcnt%m_ch)/dble(currcnt%n_ch)*dble(i) .ge. dble(j)) .and. (dble(t2)/dble(t1)*dble(i-currcnt%n_ch) .gt. dble(j-currcnt%m_ch)) .and. (dble(currcnt%m_ch)/dble(currcnt%n_ch)*dble(i-t1) .lt. dble(j-t2))) then
					k=k+1
					currcnt%posA(k,1)=dble(i)*currcnt%a1(1)+dble(j)*currcnt%a2(1)
					currcnt%posA(k,2)=dble(i)*currcnt%a1(2)+dble(j)*currcnt%a2(2)
					currcnt%posB(k,1)=currcnt%posA(k,1)+currcnt%aCC_vec(1)
					currcnt%posB(k,2)=currcnt%posA(k,2)+currcnt%aCC_vec(2)

					if (currcnt%posA(k,1) .gt. currcnt%ch_vec(1)) currcnt%posA(k,1)=currcnt%posA(k,1)-currcnt%ch_vec(1)
					if (currcnt%posA(k,1) .lt. 0) currcnt%posA(k,1)=currcnt%posA(k,1)+currcnt%ch_vec(1)
					if (currcnt%posA(k,2) .gt. currcnt%t_vec(2)) currcnt%posA(k,2)=currcnt%posA(k,2)-currcnt%t_vec(2)
					if (currcnt%posA(k,2) .lt. 0) currcnt%posA(k,2)=currcnt%posA(k,2)+currcnt%t_vec(2)

					if (currcnt%posB(k,1) .gt. currcnt%ch_vec(1)) currcnt%posB(k,1)=currcnt%posB(k,1)-currcnt%ch_vec(1)
					if (currcnt%posB(k,1) .lt. 0) currcnt%posB(k,1)=currcnt%posB(k,1)+currcnt%ch_vec(1)
					if (currcnt%posB(k,2) .gt. currcnt%t_vec(2)) currcnt%posB(k,2)=currcnt%posB(k,2)-currcnt%t_vec(2)
					if (currcnt%posB(k,2) .lt. 0) currcnt%posB(k,2)=currcnt%posB(k,2)+currcnt%t_vec(2)

				endif
				enddo
			enddo

			allocate(currcnt%pos2d(2,currcnt%Nu,2))
			currcnt%pos2d(1,:,:) = currcnt%posA(:,:)
			currcnt%pos2d(2,:,:) = currcnt%posB(:,:)

			if (k .ne. currcnt%Nu) then
				write(*,*) "*** Error in calculating atom positions ***"
				call exit()
			endif

			! calculate distances between atoms in a warped CNT unit cell.
			allocate(currcnt%posAA(currcnt%Nu,2))
			allocate(currcnt%posAB(currcnt%Nu,2))
			allocate(currcnt%posBA(currcnt%Nu,2))
			allocate(currcnt%posBB(currcnt%Nu,2))

			do i=1,currcnt%Nu
				currcnt%posAA(i,:)=currcnt%posA(i,:)-currcnt%posA(1,:)
				currcnt%posAB(i,:)=currcnt%posA(i,:)-currcnt%posB(1,:)
				currcnt%posBA(i,:)=currcnt%posB(i,:)-currcnt%posA(1,:)
				currcnt%posBB(i,:)=currcnt%posB(i,:)-currcnt%posB(1,:)
				if (currcnt%posAA(i,1) .gt. currcnt%ch_vec(1)/2.d0) currcnt%posAA(i,1)=currcnt%posAA(i,1)-currcnt%ch_vec(1)
				if (currcnt%posAB(i,1) .gt. currcnt%ch_vec(1)/2.d0) currcnt%posAB(i,1)=currcnt%posAB(i,1)-currcnt%ch_vec(1)
				if (currcnt%posBA(i,1) .gt. currcnt%ch_vec(1)/2.d0) currcnt%posBA(i,1)=currcnt%posBA(i,1)-currcnt%ch_vec(1)
				if (currcnt%posBB(i,1) .gt. currcnt%ch_vec(1)/2.d0) currcnt%posBB(i,1)=currcnt%posBB(i,1)-currcnt%ch_vec(1)
			end do

			! calculate coordinates of atoms in 3D unit cell
			allocate(currcnt%posA3(currcnt%Nu,3))
			allocate(currcnt%posB3(currcnt%Nu,3))

			do i=1,currcnt%Nu
				currcnt%posA3(i,1) = currcnt%radius*sin(2*pi*currcnt%posA(i,1)/currcnt%len_ch)
				currcnt%posA3(i,2) = currcnt%posA(i,2)
				currcnt%posA3(i,3) = -currcnt%radius*cos(2*pi*currcnt%posA(i,1)/currcnt%len_ch)
			end do

			do i=1,currcnt%Nu
				currcnt%posB3(i,1) = currcnt%radius*sin(2*pi*currcnt%posB(i,1)/currcnt%len_ch)
				currcnt%posB3(i,2) = currcnt%posB(i,2)
				currcnt%posB3(i,3) = -currcnt%radius*cos(2*pi*currcnt%posB(i,1)/currcnt%len_ch)
			end do

			allocate(currcnt%pos3d(2,currcnt%Nu,3))
			currcnt%pos3d(1,:,:) = currcnt%posA3(:,:)
			currcnt%pos3d(2,:,:) = currcnt%posB3(:,:)

		end subroutine cnt_geometry


		!**************************************************************************************************************************
		! calculate band structure of the CNT
		!**************************************************************************************************************************

		subroutine cnt_band(currcnt)
			use physicalConstants
			use math_functions_mod, only: my_norm2
			type(cnt), intent(inout) :: currcnt

			integer :: nkc, imin_sub
			integer :: i,j,mu,ik,tmpi
			integer, dimension(:), allocatable :: min_loc
			real*8 :: tmpr
			real*8, dimension(2) :: k, E1_tmp, E2_tmp, e_tmp
			real*8, dimension(:), allocatable :: k_vec,min_energy
			real*8, dimension(:,:,:), allocatable :: E_k
			complex*16, dimension(:,:,:), allocatable :: Cc_k,Cv_k
			complex*16, dimension(2) :: Cc_tmp, Cv_tmp


			! calculate CNT energy dispersion.
			nkc=2*currcnt%ikc_max+1

			allocate(k_vec(currcnt%ikc_min:currcnt%ikc_max))
			allocate(E_k(1-currcnt%Nu/2:currcnt%Nu/2,currcnt%ikc_min:currcnt%ikc_max,2))
			allocate(Cc_k(1-currcnt%Nu/2:currcnt%Nu/2,currcnt%ikc_min:currcnt%ikc_max,2))
			allocate(Cv_k(1-currcnt%Nu/2:currcnt%Nu/2,currcnt%ikc_min:currcnt%ikc_max,2))
			allocate(min_loc(0:currcnt%Nu/2))

			do ik=currcnt%ikc_min,currcnt%ikc_max
				k_vec(ik)=dble(ik)*currcnt%dk
			end do

			do mu=1-currcnt%Nu/2,currcnt%Nu/2
				do ik=currcnt%ikc_min,currcnt%ikc_max
					k = dble(mu) * currcnt%K1 + dble(ik) * currcnt%dk * currcnt%K2
					call graphene_electron(currcnt, e_tmp, Cc_tmp, Cv_tmp, k)
					E_k(mu,ik,:) = e_tmp
					Cc_k(mu,ik,:) = Cc_tmp
					Cv_k(mu,ik,:) = Cv_tmp
				enddo
			enddo

			! find the subbands with a minimum energy.
			min_loc=minloc(E_k(0:currcnt%Nu/2,:,1),2)
			imin_sub=count((min_loc .lt. nkc) .and. (min_loc .gt. 1))
			allocate(currcnt%min_sub(imin_sub))
			allocate(min_energy(imin_sub))

			! store the value of mu for subbands with minimums in the variable min_sub
			i=1
			do mu=0,currcnt%Nu/2
				if ((min_loc(mu) .gt. 1) .and. (min_loc(mu) .lt. nkc)) then
				currcnt%min_sub(i)=mu
				min_energy(i)=minval(E_k(mu,:,1))
				i=i+1
				end if
			end do

			! sort the subbands
			do i=imin_sub,2,-1
				do j=i-1,1,-1
				if (min_energy(i) .lt. min_energy(j)) then
					tmpr=min_energy(i)
					tmpi=currcnt%min_sub(i)
					min_energy(i)=min_energy(j)
					currcnt%min_sub(i)=currcnt%min_sub(j)
					min_energy(j)=tmpr
					currcnt%min_sub(j)=tmpi
				end if
				end do
			end do
			! find the max k-index that energy is below threshold energy (E_th).
			ik=0
			E1_tmp=(/ min_energy(currcnt%i_sub),0.d0 /)
			E2_tmp=(/ min_energy(currcnt%i_sub),0.d0 /)
			do while ((min(E1_tmp(1),E2_tmp(1))-min_energy(currcnt%i_sub)) .le. currcnt%E_th )
				k=dble(currcnt%min_sub(currcnt%i_sub))*currcnt%K1+dble(ik)*currcnt%dk*currcnt%K2
				call graphene_electron(currcnt,E1_tmp,Cc_tmp,Cv_tmp,k)
				k=dble(currcnt%min_sub(currcnt%i_sub))*currcnt%K1-dble(ik)*currcnt%dk*currcnt%K2
				call graphene_electron(currcnt,E2_tmp(:),Cc_tmp,Cv_tmp,k)
				ik=ik+1
			end do

			! set the index boundaries for some arrays and kernels.
			currcnt%ik_max=ik                              		!the higher limit of k-vector that is below E_th
			currcnt%ik_min=-ik                             		!the lower limit of k-vector that is below E_th
			currcnt%iKcm_max=floor(currcnt%Kcm_max/currcnt%dk)  !the higher limit of center of mass wave vector that we calculate
			currcnt%iKcm_min = - currcnt%iKcm_max             	!the lower limit of center of mass wave vector that we calculate
			currcnt%ikr_high=currcnt%iKcm_max-currcnt%ik_min    !the maximum index that the relative wavenumber in the entire simulation.
			currcnt%ikr_low=-currcnt%ikr_high                 	!the minimum index that the relative wavenumber in the entire simulation.
			currcnt%ik_high=currcnt%ikr_high+currcnt%iKcm_max   !the maximum index that the wavenumber in the entire simulation.
			currcnt%ik_low=-currcnt%ik_high                   	!the minimum index that the wavenumber in the entire simulation.
			currcnt%iq_max=2*currcnt%ikr_high                 	!the higher limit of the index in v_FT and esp_q
			currcnt%iq_min=-currcnt%iq_max                    	!the lower limit of the index in v_FT and esp_q

			currcnt%iKcm_max_fine = currcnt%iKcm_max * currcnt%dk_dkx_ratio
			currcnt%iKcm_min_fine = currcnt%iKcm_min * currcnt%dk_dkx_ratio

			select case (trim(currcnt%targetExcitonType))
			case('Ex_A1', 'Ex0_A2', 'Ex1_A2')
				currcnt%mu_cm = 0
			case('Ex0_Ep', 'Ex1_Ep')
				currcnt%mu_cm = +1 * currcnt%min_sub(currcnt%i_sub)
	! 			currcnt%mu_cm = 0
			case('Ex0_Em', 'Ex1_Em')
				currcnt%mu_cm = -1 * currcnt%min_sub(currcnt%i_sub)
	! 			currcnt%mu_cm = 0
			case default
				write(*,*) "ERROR: undetermined target exciton type!!!!"
				call exit()
			end select

			! calculate the tight-binding energies and coefficients.
			allocate(currcnt%Ek(2,currcnt%ik_low*currcnt%dk_dkx_ratio:currcnt%ik_high*currcnt%dk_dkx_ratio,2))
			allocate(currcnt%Cc(2,currcnt%ik_low*currcnt%dk_dkx_ratio:currcnt%ik_high*currcnt%dk_dkx_ratio,2))
			allocate(currcnt%Cv(2,currcnt%ik_low*currcnt%dk_dkx_ratio:currcnt%ik_high*currcnt%dk_dkx_ratio,2))

			do ik=currcnt%ik_low*currcnt%dk_dkx_ratio,currcnt%ik_high*currcnt%dk_dkx_ratio
				mu=currcnt%min_sub(currcnt%i_sub) !first band
				k=dble(mu)*currcnt%K1+dble(ik)*currcnt%dkx*currcnt%K2
				call graphene_electron(currcnt,e_tmp,Cc_tmp,Cv_tmp,k)
				currcnt%Ek(1,ik,:) = e_tmp
				currcnt%Cc(1,ik,:) = Cc_tmp
				currcnt%Cv(1,ik,:) = Cv_tmp


				mu=-currcnt%min_sub(currcnt%i_sub) !second band
				k=dble(mu)*currcnt%K1+dble(ik)*currcnt%dkx*currcnt%K2
				call graphene_electron(currcnt,e_tmp,Cc_tmp,Cv_tmp,k)
				currcnt%Ek(2,ik,:) = e_tmp
				currcnt%Cc(2,ik,:) = Cc_tmp
				currcnt%Cv(2,ik,:) = Cv_tmp
			enddo

		end subroutine cnt_band

		!**************************************************************************************************************************
		! private subroutine to calculate Bloch functions and energy in graphene
		!**************************************************************************************************************************

		subroutine graphene_electron(currCNT,E,Cc,Cv,k)
			use physicalConstants, only: i1, t0
			type(cnt), intent(in) :: currCNT
			complex*16 :: f_k
			real*8, dimension(2), intent(in) :: k
			real*8, dimension(2), intent(out) :: E
			complex*16, dimension(2), intent(out) :: Cv
			complex*16, dimension(2), intent(out) :: Cc

			f_k=exp(i1*dcmplx(dot_product(k,(currCNT%a1+currCNT%a2)/3.d0)))+exp(i1*dcmplx(dot_product(k,(currCNT%a1-2.d0*currCNT%a2)/3.d0)))+exp(i1*dcmplx(dot_product(k,(currCNT%a2-2.d0*currCNT%a1)/3.d0)))

			E(1)=+t0*abs(f_k)
			E(2)=-t0*abs(f_k)

			Cc(1)=dcmplx(+1.d0/sqrt(2.d0))
			Cc(2)=dcmplx(+1.d0/sqrt(2.d0)/abs(f_k))*conjg(f_k)
			Cv(1)=dcmplx(+1.d0/sqrt(2.d0))
			Cv(2)=dcmplx(-1.d0/sqrt(2.d0)/abs(f_k))*conjg(f_k)
		end subroutine graphene_electron

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
			if (allocated(currcnt%Ek)) deallocate(currcnt%Ek)
			if (allocated(currcnt%Sk)) deallocate(currcnt%Sk)
			if (allocated(currcnt%Cc)) deallocate(currcnt%Cc)
			if (allocated(currcnt%Cv)) deallocate(currcnt%Cv)
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

		end subroutine free_cnt_memory

end module cnt_class
