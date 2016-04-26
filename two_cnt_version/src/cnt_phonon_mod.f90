module cnt_phonon_mod
	implicit none
	private

	public :: cnt_phonon_dispersion

contains

	subroutine cnt_phonon_dispersion(currcnt)
		use cnt_class, only : cnt
		use graphene_mod, only: graphene_phonon
		use write_log_mod, only: write_log, log_input

		type(cnt) :: currcnt
		integer :: ik
		integer :: mu, ib
		real*8, dimension(2) :: k
		real*8, dimension(:), allocatable :: k_vec
		real*8, dimension(6) :: omega_tmp
		complex*16, dimension(6,6) :: u_ph
		real*8, dimension (:,:,:), allocatable :: omega_2nd_brillouin_zone
		character(len=1000) :: filename


		! ! calculate the phonon dispersion ONLY FOR THE FIRST brillouin zones
		! allocate(k_vec(currcnt%ikc_min:currcnt%ikc_max))
		! allocate(currcnt%omega_phonon(1-currcnt%Nu/2:currcnt%Nu/2,currcnt%ikc_min:currcnt%ikc_max,6))
		!
		! do ik=currcnt%ikc_min,currcnt%ikc_max
		! 	k_vec(ik)=dble(ik)*currcnt%dk
		! end do
		!
		! do mu=1-currcnt%Nu/2,currcnt%Nu/2
		! 	do ik=currcnt%ikc_min,currcnt%ikc_max
		! 		k=dble(mu)*currcnt%K1+dble(ik)*currcnt%dk*currcnt%K2
		! 		call graphene_phonon(omega_tmp,u_ph,k,currcnt%aCC_vec)
		! 		currcnt%omega_phonon(mu,ik,:) = omega_tmp(:)
		! 	enddo
		! enddo

		! calculate the phonon dispersion for two brillouin zones so that it can conserve momentum for the entire range of electron scattering
		allocate(k_vec(2*currcnt%ikc_min:2*currcnt%ikc_max))
		allocate(currcnt%omega_phonon(1-currcnt%Nu:currcnt%Nu-1,2*currcnt%ikc_min:2*currcnt%ikc_max,6))

		do ik=2*currcnt%ikc_min,2*currcnt%ikc_max
			k_vec(ik)=dble(ik)*currcnt%dk
		end do

		do mu=1-currcnt%Nu,currcnt%Nu-1
			do ik=2*currcnt%ikc_min,2*currcnt%ikc_max
				k=dble(mu)*currcnt%K1+dble(ik)*currcnt%dk*currcnt%K2
				call graphene_phonon(omega_tmp,u_ph,k,currcnt%aCC_vec)
				currcnt%omega_phonon(mu,ik,:) = omega_tmp(:)
			enddo
		enddo

		! save the CNT phonon energy dispersion*************************************************************************************
		write(filename,'(A)') trim(currcnt%name)//".phonon_k_vector.dat"
		open(unit=100,file=trim(filename),status="unknown")

		do ik=lbound(k_vec,dim=1),ubound(k_vec,dim=1)
			write(100,'(E16.8)', advance='no') k_vec(ik)
		end do

		close(100)

		write(filename,'(A)') trim(currcnt%name)//".phonon_energy.dat"
		open(unit=100,file=trim(filename),status="unknown")

		do ib=1,6
			do mu=lbound(currcnt%omega_phonon,dim=1),ubound(currcnt%omega_phonon,dim=1)
				do ik=lbound(currcnt%omega_phonon,dim=2),ubound(currcnt%omega_phonon,dim=2)
					write(100,'(E16.8)', advance='no') currcnt%omega_phonon(mu,ik,ib)
				end do
				write(100,*)
			enddo
		enddo

		close(100)

		write(log_input,'(A,I2.2,A,I2.2,A)') new_line('A')//"phonon dispersion calculated for carbon nanotube with chirality: (", currcnt%n_ch, ",", currcnt%m_ch, ")"
		call write_log(trim(log_input))

	end subroutine cnt_phonon_dispersion

end module cnt_phonon_mod
