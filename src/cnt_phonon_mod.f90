module cnt_phonon_mod
	implicit none
	private

	public :: cnt_phonon_dispersion

contains

	!***************************************************************************
	! - this subroutine calculates the phonon energy dispersion
	!***************************************************************************
	subroutine cnt_phonon_dispersion(currcnt, dq, iq_max, iq_min, mu_min, mu_max, save_dispersion)
		use cnt_class, only : cnt
		use graphene_mod, only: graphene_phonon
		use write_log_mod, only: write_log, log_input

		type(cnt), intent(inout) :: currcnt
		real*8, optional, intent(in) :: dq
		integer, optional, intent(in) :: iq_max, iq_min, mu_min, mu_max
		logical, optional, intent(in) :: save_dispersion

		integer :: iq
		integer :: mu, ib
		real*8, dimension(2) :: q
		real*8, dimension(:), allocatable :: q_vec
		real*8, dimension(6) :: omega_tmp
		complex*16, dimension(6,6) :: u_ph
		character(len=1000) :: filename

		! these are the input variables that are used in the subroutines
		real*8 :: dq_used
		integer :: iq_max_used, iq_min_used, mu_min_used, mu_max_used
		logical :: save_dispersion_used

		!***********************************************************************
		! set default values for optional arguments
		dq_used = currcnt%dk
		iq_max_used = currcnt%ikc_max
		iq_min_used = currcnt%ikc_min
		mu_max_used = currcnt%Nu/2
		mu_min_used = 1-currcnt%Nu/2
		save_dispersion_used = .false.

		if(present(dq)) dq_used = currcnt%dk
		if(present(iq_max)) iq_max_used = iq_max
		if(present(iq_min)) iq_min_used = iq_min
		if(present(mu_max)) mu_max_used = mu_max
		if(present(mu_min)) mu_min_used = mu_min
		if(present(save_dispersion)) save_dispersion_used = save_dispersion

		!***********************************************************************
		! calculate the phonon dispersion
		allocate(q_vec(iq_min_used:iq_max_used))
		if(allocated(currcnt%omega_phonon)) deallocate(currcnt%omega_phonon)
		allocate(currcnt%omega_phonon(mu_min_used:mu_max_used,iq_min_used:iq_max_used,6))

		do iq=iq_min_used,iq_max_used
			q_vec(iq)=dble(iq)*dq_used
		end do

		do mu=mu_min_used,mu_max_used
			do iq=iq_min_used,iq_max_used
				q=dble(mu)*currcnt%K1+dble(iq)*dq_used*currcnt%K2
				call graphene_phonon(omega_tmp,u_ph,q,currcnt%aCC_vec)
				currcnt%omega_phonon(mu,iq,:) = omega_tmp(:)
			enddo
		enddo

		!***********************************************************************
		! save the CNT phonon energy dispersion
		if(save_dispersion_used) then
			write(filename,'(A)') trim(currcnt%name)//".phonon_k_vector.dat"
			open(unit=100,file=trim(filename),status="unknown")

			do iq=lbound(q_vec,dim=1),ubound(q_vec,dim=1)
				write(100,'(E16.8)', advance='no') q_vec(iq)
			end do

			close(100)

			write(filename,'(A)') trim(currcnt%name)//".phonon_energy.dat"
			open(unit=100,file=trim(filename),status="unknown")

			do ib=1,6
				do mu=lbound(currcnt%omega_phonon,dim=1),ubound(currcnt%omega_phonon,dim=1)
					do iq=lbound(currcnt%omega_phonon,dim=2),ubound(currcnt%omega_phonon,dim=2)
						write(100,'(E16.8)', advance='no') currcnt%omega_phonon(mu,iq,ib)
					end do
					write(100,*)
				enddo
			enddo

			close(100)

			write(log_input,'(A,I0,A,I0,A)') new_line('A')//"phonon dispersion calculated for carbon nanotube with chirality: (", currcnt%n_ch, ",", currcnt%m_ch, ")"
			call write_log(trim(log_input))
		endif

	end subroutine cnt_phonon_dispersion

end module cnt_phonon_mod
