module cnt_scattering_electron_phonon_mod
	implicit none
	private
	public :: cnt_scattering_electron_phonon

contains

	!**************************************************************************************************************************
	! calculate band structure of the CNT
	!**************************************************************************************************************************

	subroutine cnt_scattering_electron_phonon(currcnt)
		use cnt_class, only: cnt
		use constants_mod
		use graphene_mod, only: graphene_electron, graphene_electron_phonon
		! use write_log_mod, only: write_log, log_input

		type(cnt), intent(inout) :: currcnt

		integer :: nkc
		integer :: mu,ik
		integer :: mu_q, iq
		integer :: ib
		real*8, dimension(2) :: k, q
		real*8, dimension(2) :: e_tmp
		real*8, dimension(:), allocatable :: k_vec
		real*8, dimension(:,:,:), allocatable :: E_k
		complex*16, dimension(2) :: Cc_tmp, Cv_tmp
		character(len=1000) :: filename
		complex*16, dimension(6) :: f_tilde_1_tmp, f_tilde_2_tmp
		complex*16, dimension(:,:,:), allocatable :: f_tilde_1, f_tilde_2


		! calculate CNT energy dispersion.
		nkc=2*currcnt%ikc_max+1

		allocate(k_vec(currcnt%ikc_min:currcnt%ikc_max))
		allocate(E_k(1-currcnt%Nu/2:currcnt%Nu/2,currcnt%ikc_min:currcnt%ikc_max,2))
		! allocate(Cc_k(1-currcnt%Nu/2:currcnt%Nu/2,currcnt%ikc_min:currcnt%ikc_max,2))
		! allocate(Cv_k(1-currcnt%Nu/2:currcnt%Nu/2,currcnt%ikc_min:currcnt%ikc_max,2))
		! allocate(min_loc(0:currcnt%Nu/2))

		do ik=currcnt%ikc_min,currcnt%ikc_max
			k_vec(ik)=dble(ik)*currcnt%dk
		end do

		do mu=1-currcnt%Nu/2,currcnt%Nu/2
			do ik=currcnt%ikc_min,currcnt%ikc_max
				k = dble(mu) * currcnt%K1 + dble(ik) * currcnt%dk * currcnt%K2
				call graphene_electron(e_tmp,Cc_tmp,Cv_tmp,k,currcnt%a1,currcnt%a2)
				E_k(mu,ik,:) = e_tmp
			enddo
		enddo

		! calculate electron-phonon interaction matrix element
		allocate(f_tilde_1(1-currcnt%Nu/2:currcnt%Nu/2,currcnt%ikc_min:currcnt%ikc_max,6))
		allocate(f_tilde_2(1-currcnt%Nu/2:currcnt%Nu/2,currcnt%ikc_min:currcnt%ikc_max,6))
		k = (/0.d0, 0.d0/)
		do mu_q=1-currcnt%Nu/2,currcnt%Nu/2
			do iq=currcnt%ikc_min,currcnt%ikc_max
				q = dble(mu_q) * currcnt%K1 + dble(iq) * currcnt%dk * currcnt%K2
				call graphene_electron_phonon(f_tilde_1_tmp,f_tilde_2_tmp,k,q,currcnt%a1,currcnt%a2)
				f_tilde_1(mu_q,iq,:) = f_tilde_1_tmp
				f_tilde_2(mu_q,iq,:) = f_tilde_2_tmp
			enddo
		enddo

		! save the CNT electron energy dispersion*************************************************************************************
		write(filename,'(A)') trim(currcnt%name)//".electron_k_vector.dat"
		open(unit=100,file=trim(filename),status="unknown")

		do ik=currcnt%ikc_min,currcnt%ikc_max
			write(100,'(E16.8)', advance='no') k_vec(ik)
		end do

		close(100)


		write(filename,'(A)') trim(currcnt%name)//".electron_conduction_band.dat"
		open(unit=100,file=trim(filename),status="unknown")
		write(filename,'(A)') trim(currcnt%name)//".electron_valence_band.dat"
		open(unit=101,file=trim(filename),status="unknown")

		do mu=1-currcnt%Nu/2,currcnt%Nu/2
			do ik=currcnt%ikc_min,currcnt%ikc_max
				write(100,'(E16.8)', advance='no') E_k(mu,ik,1)
				write(101,'(E16.8)', advance='no') E_k(mu,ik,2)
			enddo
			write(100,*)
			write(101,*)
		enddo

		close(100)
		close(101)

		! save the CNT electron-phonon matrix element*************************************************************************************
		write(filename,'(A)') trim(currcnt%name)//".electron_phonon_matrix_element_branch_1.dat"
		open(unit=100,file=trim(filename),status="unknown")
		write(filename,'(A)') trim(currcnt%name)//".electron_phonon_matrix_element_branch_2.dat"
		open(unit=101,file=trim(filename),status="unknown")
		write(filename,'(A)') trim(currcnt%name)//".electron_phonon_matrix_element_branch_3.dat"
		open(unit=102,file=trim(filename),status="unknown")
		write(filename,'(A)') trim(currcnt%name)//".electron_phonon_matrix_element_branch_4.dat"
		open(unit=103,file=trim(filename),status="unknown")
		write(filename,'(A)') trim(currcnt%name)//".electron_phonon_matrix_element_branch_5.dat"
		open(unit=104,file=trim(filename),status="unknown")
		write(filename,'(A)') trim(currcnt%name)//".electron_phonon_matrix_element_branch_6.dat"
		open(unit=105,file=trim(filename),status="unknown")

		! do mu_q=1-currcnt%Nu/2,currcnt%Nu/2
			mu_q=0
			do iq=currcnt%ikc_min,currcnt%ikc_max
				write(100,'(E16.8)', advance='no') (abs(f_tilde_1(mu_q,iq,1)+f_tilde_2(mu_q,iq,1)))**2
				write(101,'(E16.8)', advance='no') (abs(f_tilde_1(mu_q,iq,2)+f_tilde_2(mu_q,iq,2)))**2
				write(102,'(E16.8)', advance='no') (abs(f_tilde_1(mu_q,iq,3)+f_tilde_2(mu_q,iq,3)))**2
				write(103,'(E16.8)', advance='no') (abs(f_tilde_1(mu_q,iq,4)+f_tilde_2(mu_q,iq,4)))**2
				write(104,'(E16.8)', advance='no') (abs(f_tilde_1(mu_q,iq,5)+f_tilde_2(mu_q,iq,5)))**2
				write(105,'(E16.8)', advance='no') (abs(f_tilde_1(mu_q,iq,6)+f_tilde_2(mu_q,iq,6)))**2
			enddo
			write(100,*)
			write(101,*)
			write(102,*)
			write(103,*)
			write(104,*)
			write(105,*)
		! enddo

		close(100)
		close(101)
		close(102)
		close(103)
		close(104)
		close(105)

	end subroutine cnt_scattering_electron_phonon

end module cnt_scattering_electron_phonon_mod
