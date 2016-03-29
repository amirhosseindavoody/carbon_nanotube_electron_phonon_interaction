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

		allocate(k_vec(currcnt%ikc_min:currcnt%ikc_max))
		allocate(currcnt%omega_phonon(1-currcnt%Nu/2:currcnt%Nu/2,currcnt%ikc_min:currcnt%ikc_max,6))

		do ik=currcnt%ikc_min,currcnt%ikc_max
			k_vec(ik)=dble(ik)*currcnt%dk
		end do

		! write(*,*) !this write command is necessary to prevent a crash in the eigen value solver inside graphene_phonon_dispersion function.

		do mu=1-currcnt%Nu/2,currcnt%Nu/2
			do ik=currcnt%ikc_min,currcnt%ikc_max
				k=dble(mu)*currcnt%K1+dble(ik)*currcnt%dk*currcnt%K2
				call graphene_phonon(omega_tmp,u_ph,k,currcnt%aCC_vec)
				currcnt%omega_phonon(mu,ik,:) = omega_tmp(:)
			enddo
		enddo

		! save the CNT phonon energy dispersion*************************************************************************************
		open(unit=100,file='phonon_dispersion.dat',status="unknown")

		do ik=currcnt%ikc_min,currcnt%ikc_max
			write(100,'(E16.8)', advance='no') k_vec(ik)
		end do

		write(100,*)
		write(101,*)

		do ib=1,6
			do mu=1-currcnt%Nu/2,currcnt%Nu/2
				do ik=currcnt%ikc_min,currcnt%ikc_max
					write(100,'(E16.8)', advance='no') currcnt%omega_phonon(mu,ik,ib)
				end do
				write(100,*)
			enddo
		enddo
		
		close(100)

	end subroutine cnt_phonon_dispersion

end module cnt_phonon_mod
