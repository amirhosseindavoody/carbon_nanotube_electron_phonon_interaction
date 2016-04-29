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
		use math_functions_mod, only: find_all_roots
		use write_log_mod, only: write_log, log_input

		type(cnt), intent(inout) :: currcnt

		integer :: nkc
		integer :: mu_e,ik_e, mu_e_2
		integer :: mu_ph, iq_ph
		integer :: ib
		integer :: number_of_scattered_states
		integer :: n_root
		integer, dimension(:), allocatable :: root_idx
		integer :: tmp_i
		real*8, dimension(2) :: k_e, q_ph
		real*8, dimension(2) :: e_tmp
		real*8, dimension(:), allocatable :: k_vec
		real*8, dimension(:), allocatable :: tmp_real_array
		real*8, dimension(:,:,:), allocatable :: E_k
		complex*16, dimension(2) :: Cc_tmp, Cv_tmp
		character(len=1000) :: filename
		complex*16, dimension(6) :: f_tilde_1_tmp, f_tilde_2_tmp
		complex*16, dimension(:,:,:), allocatable :: f_tilde_1, f_tilde_2
		integer, dimension(:,:), allocatable :: scattering_state_list

		! calculate CNT electronic energy dispersion.
		nkc=2*currcnt%ikc_max+1

		allocate(k_vec(currcnt%ikc_min:currcnt%ikc_max))
		allocate(E_k(1-currcnt%Nu/2:currcnt%Nu/2,currcnt%ikc_min:currcnt%ikc_max,2))

		do ik_e=currcnt%ikc_min,currcnt%ikc_max
			k_vec(ik_e)=dble(ik_e)*currcnt%dk
		end do

		do mu_e=1-currcnt%Nu/2,currcnt%Nu/2
			do ik_e=currcnt%ikc_min,currcnt%ikc_max
				k_e = dble(mu_e) * currcnt%K1 + dble(ik_e) * currcnt%dk * currcnt%K2
				call graphene_electron(e_tmp,Cc_tmp,Cv_tmp,k_e,currcnt%a1,currcnt%a2)
				E_k(mu_e,ik_e,:) = e_tmp
			enddo
		enddo

		! save the CNT electron energy dispersion*************************************************************************************
		write(filename,'(A)') trim(currcnt%name)//".electron_k_vector.dat"
		open(unit=100,file=trim(filename),status="unknown")

		do ik_e=currcnt%ikc_min,currcnt%ikc_max
			write(100,'(E16.8)', advance='no') k_vec(ik_e)
		end do

		close(100)


		write(filename,'(A)') trim(currcnt%name)//".electron_conduction_band.dat"
		open(unit=100,file=trim(filename),status="unknown")
		write(filename,'(A)') trim(currcnt%name)//".electron_valence_band.dat"
		open(unit=101,file=trim(filename),status="unknown")

		do mu_e=1-currcnt%Nu/2,currcnt%Nu/2
			do ik_e=currcnt%ikc_min,currcnt%ikc_max
				write(100,'(E16.8)', advance='no') E_k(mu_e,ik_e,1)
				write(101,'(E16.8)', advance='no') E_k(mu_e,ik_e,2)

			enddo
			write(100,*)
			write(101,*)

		enddo

		close(100)
		close(101)

		! calculate electron-phonon interaction matrix element*************************************************************************************
		allocate(f_tilde_1(1-currcnt%Nu/2:currcnt%Nu/2,currcnt%ikc_min:currcnt%ikc_max,6))
		allocate(f_tilde_2(1-currcnt%Nu/2:currcnt%Nu/2,currcnt%ikc_min:currcnt%ikc_max,6))
		k_e = (/0.d0, 0.d0/)
		do mu_ph=1-currcnt%Nu/2,currcnt%Nu/2
			do iq_ph=currcnt%ikc_min,currcnt%ikc_max
				q_ph = dble(mu_ph) * currcnt%K1 + dble(iq_ph) * currcnt%dk * currcnt%K2
				call graphene_electron_phonon(f_tilde_1_tmp,f_tilde_2_tmp,k_e,q_ph,currcnt%a1,currcnt%a2)
				f_tilde_1(mu_ph,iq_ph,:) = f_tilde_1_tmp
				f_tilde_2(mu_ph,iq_ph,:) = f_tilde_2_tmp
			enddo
		enddo

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

		do mu_ph=1-currcnt%Nu/2,currcnt%Nu/2
			do iq_ph=currcnt%ikc_min,currcnt%ikc_max
				write(100,'(E16.8)', advance='no') (abs(f_tilde_1(mu_ph,iq_ph,1)+f_tilde_2(mu_ph,iq_ph,1)))**2
				write(101,'(E16.8)', advance='no') (abs(f_tilde_1(mu_ph,iq_ph,2)+f_tilde_2(mu_ph,iq_ph,2)))**2
				write(102,'(E16.8)', advance='no') (abs(f_tilde_1(mu_ph,iq_ph,3)+f_tilde_2(mu_ph,iq_ph,3)))**2
				write(103,'(E16.8)', advance='no') (abs(f_tilde_1(mu_ph,iq_ph,4)+f_tilde_2(mu_ph,iq_ph,4)))**2
				write(104,'(E16.8)', advance='no') (abs(f_tilde_1(mu_ph,iq_ph,5)+f_tilde_2(mu_ph,iq_ph,5)))**2
				write(105,'(E16.8)', advance='no') (abs(f_tilde_1(mu_ph,iq_ph,6)+f_tilde_2(mu_ph,iq_ph,6)))**2
			enddo
			write(100,*)
			write(101,*)
			write(102,*)
			write(103,*)
			write(104,*)
			write(105,*)
		enddo

		close(100)
		close(101)
		close(102)
		close(103)
		close(104)
		close(105)

		!find phonon momentum that conserves both momentum and energy*************************************************************************************

		allocate(tmp_real_array(currcnt%ikc_min:currcnt%ikc_max))
		allocate(root_idx(2*currcnt%ikc_max+1))

		!first count the number of scattered states to be able to allocate the need memory for scattered states.
		number_of_scattered_states = 0
		do ib=1,6
			do mu_e=1-currcnt%Nu/2,currcnt%Nu/2
				do mu_e_2=1-currcnt%Nu/2,currcnt%Nu/2
					mu_ph = mu_e-mu_e_2
					do ik_e = currcnt%ikc_min,currcnt%ikc_max

						tmp_real_array = E_k(mu_e_2,:,1)-E_k(mu_e,ik_e,1)+currcnt%omega_phonon(mu_ph,ik_e+currcnt%ikc_max:ik_e+currcnt%ikc_min:-1,ib)

						call find_all_roots (tmp_real_array, currcnt%ikc_min, currcnt%ikc_max, 0.d0, n_root, root_idx)
						number_of_scattered_states = number_of_scattered_states + n_root

					enddo
				enddo
			enddo
		enddo

		!now allocate the scattering_state_list and save the information of the scattering states.
		allocate(scattering_state_list(number_of_scattered_states,8))
		scattering_state_list = 0
		number_of_scattered_states = 0

		do ib=1,6
			do mu_e=1-currcnt%Nu/2,currcnt%Nu/2
				do mu_e_2=1-currcnt%Nu/2,currcnt%Nu/2
					mu_ph = mu_e-mu_e_2
					do ik_e = currcnt%ikc_min,currcnt%ikc_max

						tmp_real_array = E_k(mu_e_2,:,1)-E_k(mu_e,ik_e,1)+currcnt%omega_phonon(mu_ph,ik_e+currcnt%ikc_max:ik_e+currcnt%ikc_min:-1,ib)

						call find_all_roots (tmp_real_array, currcnt%ikc_min, currcnt%ikc_max, 0.d0, n_root, root_idx)

						do tmp_i= 1, n_root
							iq_ph = ik_e - root_idx(tmp_i)
							scattering_state_list(number_of_scattered_states + tmp_i,:) = (/ number_of_scattered_states+tmp_i, ib, mu_e, mu_e_2, mu_ph, ik_e, root_idx(tmp_i), iq_ph   /)
						enddo

						number_of_scattered_states = number_of_scattered_states + n_root

					enddo
				enddo
			enddo
		enddo


		!save the calculated list of scattering states that conserve both momentum and energy *************************************************************************************
		write(filename,'(A)') trim(currcnt%name)//".electron_phonon_scattering_states.dat"
		open(unit=100,file=trim(filename),status="unknown")

		do tmp_i = 1, number_of_scattered_states
			write(100,'(SP,I10.9, A, I10.9, A, I10.9, A, I10.9, A, I10.9, A, I10.9, A, I10.9, A, I10.9)') scattering_state_list(tmp_i,1), '   ', scattering_state_list(tmp_i,2), '   ', scattering_state_list(tmp_i,3), '   ', scattering_state_list(tmp_i,4), '   ', scattering_state_list(tmp_i,5), '   ', scattering_state_list(tmp_i,6), '   ', scattering_state_list(tmp_i,7), '   ', scattering_state_list(tmp_i,8)
		enddo

		close(100)

		write(log_input,'(A,I10.10)') "Number of scattered states is : ", number_of_scattered_states
		call write_log(log_input)


	end subroutine cnt_scattering_electron_phonon

end module cnt_scattering_electron_phonon_mod
