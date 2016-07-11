module cnt_scattering_electron_phonon_mod
	implicit none
	private
	public :: cnt_electron_phonon_scattering_rate_emission, cnt_electron_phonon_scattering_rate_absorption, cnt_electron_phonon_matrix_element, cnt_electron_phonon_scattering_states

	real*8, private :: energy_mesh_min, energy_mesh_max
	real*8, private :: energy_mesh_length = 2.0d0 ! this is the energy distance between energy_mesh_max and energy_mesh_min

	real*8, dimension(:,:,:), allocatable :: E_k

contains

	!***************************************************************************
	! calculate the electron-phonon scattering rates due to phonon emission
	!***************************************************************************
	subroutine cnt_electron_phonon_scattering_rate_emission(currcnt)
		use cnt_class, only: cnt
		use cnt_phonon_mod, only: cnt_phonon_dispersion
		use constants_mod
		use graphene_mod, only: graphene_electron, graphene_electron_phonon_matrix_element
		use math_functions_mod, only: find_all_roots, first_derivative
		use sim_properties_mod, only: temperature, energy_mesh_size
		use write_log_mod, only: write_log, log_input

		type(cnt), intent(inout) :: currcnt

		integer :: i, j, k
		real*8, dimension(:), allocatable :: tmp_real_array_1
		real*8, dimension(:), allocatable :: tmp_real_array_2

		integer :: mu_e,ik_e
		integer :: mu_e_2, ik_e_2
		integer :: mu_ph, iq_ph
		integer :: ib
		integer :: n_root
		integer, dimension(:), allocatable :: root_idx
		real*8, dimension(2) :: k_e, q_ph
		character(len=1000) :: filename

		real*8, dimension(:), allocatable :: energy_mesh
		integer :: n_scattering_state
		integer, dimension(:), allocatable :: scattering_state_idx
		real*8, dimension(:,:), allocatable :: electron_phonon_scattering_rate
		real*8 :: Ec_derivative, omega_phonon_derivative
		real*8 :: omega_tmp
		complex*16 :: matrix_element

		!***********************************************************************
		! calculate CNT electronic energy dispersion.
		call calculate_cnt_electron_energy_dispersion(currcnt)

		!***********************************************************************
		! calculate CNT phonon energy dispersion for two brillouine zones.
		call cnt_phonon_dispersion(currcnt, iq_max=2*currcnt%ikc_max, iq_min=2*currcnt%ikc_min, mu_max=currcnt%Nu-1, mu_min=1-currcnt%Nu )

		!***********************************************************************
		!calculate electron-phonon scattering rate as a function of energy

		!create the electron energy mesh
		allocate(energy_mesh(energy_mesh_size))

		energy_mesh_min = minval(E_k(:,:,1))
		energy_mesh_max = energy_mesh_min + energy_mesh_length*eV

		do i=1,energy_mesh_size
			energy_mesh(i) = energy_mesh_min + dble(i-1)*(energy_mesh_max-energy_mesh_min)/dble(energy_mesh_size-1)
		enddo

		allocate(electron_phonon_scattering_rate(energy_mesh_size,6))
		electron_phonon_scattering_rate = 0.d0

		!now solve find electronic states with energy equal to each one of the energy mesh points
		if (allocated(tmp_real_array_1)) deallocate(tmp_real_array_1)
		allocate(tmp_real_array_1(currcnt%ikc_min:currcnt%ikc_max))

		if (allocated(tmp_real_array_2)) deallocate(tmp_real_array_2)
		allocate(tmp_real_array_2(2*currcnt%ikc_min:2*currcnt%ikc_max))

		if(allocated(root_idx)) deallocate(root_idx)
		allocate(root_idx(2*currcnt%ikc_max+1))

		if(allocated(scattering_state_idx)) deallocate(scattering_state_idx)
		allocate(scattering_state_idx(2*currcnt%ikc_max+1))

		do i=1,energy_mesh_size
			do mu_e=1-currcnt%Nu/2,currcnt%Nu/2
				tmp_real_array_1 = E_k(mu_e,:,1)
				call find_all_roots (tmp_real_array_1, currcnt%ikc_min, currcnt%ikc_max, energy_mesh(i), n_root, root_idx)
				!find all the scattering states for the roots that has been found before
				do j=1,n_root
					ik_e = root_idx(j)
					k_e = dble(mu_e) * currcnt%K1 + dble(ik_e) * currcnt%dk * currcnt%K2
					do ib=1,6
						do mu_e_2=1-currcnt%Nu/2,currcnt%Nu/2
							mu_ph = mu_e-mu_e_2
							tmp_real_array_1 = E_k(mu_e_2,:,1)-E_k(mu_e,ik_e,1)+currcnt%omega_phonon(mu_ph,ik_e+currcnt%ikc_max:ik_e+currcnt%ikc_min:-1,ib)
							call find_all_roots (tmp_real_array_1, currcnt%ikc_min, currcnt%ikc_max, 0.d0, n_scattering_state, scattering_state_idx)

							do k = 1, n_scattering_state
								ik_e_2 = scattering_state_idx(k)
								iq_ph = ik_e - ik_e_2

								q_ph = dble(mu_ph) * currcnt%K1 + dble(iq_ph) * currcnt%dk * currcnt%K2
								call graphene_electron_phonon_matrix_element(matrix_element, k_e, q_ph, ib, currcnt%a1, currcnt%a2)

								tmp_real_array_1 = E_k(mu_e_2,:,1)
								call first_derivative(tmp_real_array_1, currcnt%ikc_min, currcnt%ikc_max, ik_e_2, currcnt%dk, Ec_derivative)
								tmp_real_array_2 = currcnt%omega_phonon(mu_ph, :, ib)
								call first_derivative(tmp_real_array_2, 2*currcnt%ikc_min, 2*currcnt%ikc_max, iq_ph, currcnt%dk, omega_phonon_derivative)

								omega_tmp = currcnt%omega_phonon(mu_ph, iq_ph, ib)

								electron_phonon_scattering_rate(i,ib) = electron_phonon_scattering_rate(i,ib) + (1.d0+1.d0/(exp(omega_tmp/temperature/kb)-1.d0)) * ((abs(matrix_element))**2)/(abs(-Ec_derivative+omega_phonon_derivative)*omega_tmp/hb)


							enddo

						enddo
					enddo

				enddo
			enddo
		enddo

		electron_phonon_scattering_rate = ((g0**2)*A_u/(16.d0*pi*m_carbon*currcnt%radius))*electron_phonon_scattering_rate

		!***********************************************************************
		!save the calculated electron-phonon scattering rate

		write(filename,'(A)') trim(currcnt%name)//".electron_phonon_scattering_rate_emission.dat"
		open(unit=100,file=trim(filename),status="unknown")

		do i = 1, energy_mesh_size
			write(100,'(SP, E20.8)', advance='no') energy_mesh(i)
		enddo

		write(100,*)

		do ib=1,6
			do i = 1, energy_mesh_size
				write(100,'(SP, E20.8)', advance='no') electron_phonon_scattering_rate(i,ib)
			enddo
			write(100,*)
		enddo

		close(100)

		write(log_input,'(A)') "Electron-phonon scattering rates due to phonon emission calculated and saved!";		call write_log(log_input)
	end subroutine cnt_electron_phonon_scattering_rate_emission


	!***************************************************************************
	! calculate the electron-phonon scattering rates due to phonon absorption
	!***************************************************************************
	subroutine cnt_electron_phonon_scattering_rate_absorption(currcnt)
		use cnt_class, only: cnt
		use cnt_phonon_mod, only: cnt_phonon_dispersion
		use constants_mod
		use graphene_mod, only: graphene_electron, graphene_electron_phonon_matrix_element
		use math_functions_mod, only: find_all_roots, first_derivative
		use sim_properties_mod, only: temperature, energy_mesh_size
		use write_log_mod, only: write_log, log_input

		type(cnt), intent(inout) :: currcnt

		integer :: i, j, k
		real*8, dimension(:), allocatable :: tmp_real_array_1
		real*8, dimension(:), allocatable :: tmp_real_array_2

		integer :: mu_e,ik_e
		integer :: mu_e_2, ik_e_2
		integer :: mu_ph, iq_ph
		integer :: ib
		integer :: n_root
		integer, dimension(:), allocatable :: root_idx
		real*8, dimension(2) :: k_e, q_ph
		character(len=1000) :: filename

		real*8, dimension(:), allocatable :: energy_mesh
		integer :: n_scattering_state
		integer, dimension(:), allocatable :: scattering_state_idx
		real*8, dimension(:,:), allocatable :: electron_phonon_scattering_rate
		real*8 :: Ec_derivative, omega_phonon_derivative
		real*8 :: omega_tmp
		complex*16 :: matrix_element


		!***********************************************************************
		! calculate CNT electronic energy dispersion.
		call calculate_cnt_electron_energy_dispersion(currcnt)

		!***********************************************************************
		! calculate CNT phonon energy dispersion for two brillouine zones.
		call cnt_phonon_dispersion(currcnt, iq_max=2*currcnt%ikc_max, iq_min=2*currcnt%ikc_min, mu_max=currcnt%Nu-1, mu_min=1-currcnt%Nu )

		!***********************************************************************
		!calculate electron-phonon scattering rate as a function of energy

		!create the electron energy mesh
		allocate(energy_mesh(energy_mesh_size))

		energy_mesh_min = minval(E_k(:,:,1))
		energy_mesh_max = energy_mesh_min + energy_mesh_length*eV

		do i=1,energy_mesh_size
			energy_mesh(i) = energy_mesh_min + dble(i-1)*(energy_mesh_max-energy_mesh_min)/dble(energy_mesh_size-1)
		enddo

		allocate(electron_phonon_scattering_rate(energy_mesh_size,6))
		electron_phonon_scattering_rate = 0.d0

		!now solve find electronic states with energy equal to each one of the energy mesh points
		if (allocated(tmp_real_array_1)) deallocate(tmp_real_array_1)
		allocate(tmp_real_array_1(currcnt%ikc_min:currcnt%ikc_max))

		if (allocated(tmp_real_array_2)) deallocate(tmp_real_array_2)
		allocate(tmp_real_array_2(2*currcnt%ikc_min:2*currcnt%ikc_max))

		if(allocated(root_idx)) deallocate(root_idx)
		allocate(root_idx(2*currcnt%ikc_max+1))

		if(allocated(scattering_state_idx)) deallocate(scattering_state_idx)
		allocate(scattering_state_idx(2*currcnt%ikc_max+1))

		do i=1,energy_mesh_size
			do mu_e=1-currcnt%Nu/2,currcnt%Nu/2
				tmp_real_array_1 = E_k(mu_e,:,1)
				call find_all_roots (tmp_real_array_1, currcnt%ikc_min, currcnt%ikc_max, energy_mesh(i), n_root, root_idx)
				!find all the scattering states for the roots that has been found before
				do j=1,n_root
					ik_e = root_idx(j)
					k_e = dble(mu_e) * currcnt%K1 + dble(ik_e) * currcnt%dk * currcnt%K2
					do ib=1,6
						do mu_e_2=1-currcnt%Nu/2,currcnt%Nu/2
							mu_ph = mu_e-mu_e_2
							tmp_real_array_1 = E_k(mu_e_2,:,1)-E_k(mu_e,ik_e,1)-currcnt%omega_phonon(-mu_ph,currcnt%ikc_min-ik_e:currcnt%ikc_max-ik_e,ib)
							call find_all_roots (tmp_real_array_1, currcnt%ikc_min, currcnt%ikc_max, 0.d0, n_scattering_state, scattering_state_idx)

							do k = 1, n_scattering_state
								ik_e_2 = scattering_state_idx(k)
								iq_ph = ik_e - ik_e_2

								q_ph = dble(mu_ph) * currcnt%K1 + dble(iq_ph) * currcnt%dk * currcnt%K2
								call graphene_electron_phonon_matrix_element(matrix_element, k_e, q_ph, ib, currcnt%a1, currcnt%a2)

								tmp_real_array_1 = E_k(mu_e_2,:,1)
								call first_derivative(tmp_real_array_1, currcnt%ikc_min, currcnt%ikc_max, ik_e_2, currcnt%dk, Ec_derivative)
								tmp_real_array_2 = currcnt%omega_phonon(-mu_ph, :, ib)
								call first_derivative(tmp_real_array_2, 2*currcnt%ikc_min, 2*currcnt%ikc_max, -iq_ph, currcnt%dk, omega_phonon_derivative)

								omega_tmp = currcnt%omega_phonon(-mu_ph, -iq_ph, ib)

								electron_phonon_scattering_rate(i,ib) = electron_phonon_scattering_rate(i,ib) + (1.d0/(exp(omega_tmp/temperature/kb)-1.d0)) * ((abs(matrix_element))**2)/(abs(-Ec_derivative+omega_phonon_derivative)*omega_tmp/hb)

							enddo

						enddo
					enddo

				enddo
			enddo
		enddo

		electron_phonon_scattering_rate = ((g0**2)*A_u/(16.d0*pi*m_carbon*currcnt%radius))*electron_phonon_scattering_rate

		!***********************************************************************
		!save the calculated electron-phonon scattering rate

		write(filename,'(A)') trim(currcnt%name)//".electron_phonon_scattering_rate_absorption.dat"
		open(unit=100,file=trim(filename),status="unknown")

		do i = 1, energy_mesh_size
			write(100,'(SP, E20.8)', advance='no') energy_mesh(i)
		enddo

		write(100,*)

		do ib = 1,6
			do i = 1, energy_mesh_size
				write(100,'(SP, E20.8)', advance='no') electron_phonon_scattering_rate(i,ib)
			enddo
			write(100,*)
		enddo

		close(100)

		write(log_input,'(A)') "Electron-phonon scattering rates due to phonon absorption calculated and saved!"
		call write_log(log_input)
	end subroutine cnt_electron_phonon_scattering_rate_absorption

	!***************************************************************************
	! - this subroutine calculates the electron energy dispersion that is used
	!   in this calculating the electron scattering rates.
	!***************************************************************************
	subroutine calculate_cnt_electron_energy_dispersion(currcnt)
		use cnt_class, only: cnt
		use constants_mod
		use graphene_mod, only: graphene_electron
		! use write_log_mod, only: write_log, log_input

		type(cnt), intent(in) :: currcnt

		integer :: nkc
		integer :: mu_e,ik_e
		real*8, dimension(2) :: k_e
		real*8, dimension(2) :: e_tmp
		real*8, dimension(:), allocatable :: k_vec
		complex*16, dimension(2) :: Cc_tmp, Cv_tmp
		character(len=1000) :: filename


		!***********************************************************************
		! calculate CNT electronic energy dispersion.
		nkc=2*currcnt%ikc_max+1

		allocate(k_vec(currcnt%ikc_min:currcnt%ikc_max))

		if (allocated(E_k)) deallocate(E_k)
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

		!***********************************************************************
		! save the CNT electron energy dispersion

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
	end subroutine calculate_cnt_electron_energy_dispersion

	!***************************************************************************
	! calculate and save the electron-phonon matrix elements
	!***************************************************************************
	subroutine cnt_electron_phonon_matrix_element(currcnt)
		use cnt_class, only: cnt
		use constants_mod
		use graphene_mod, only: graphene_electron, graphene_electron_phonon_matrix_element
		use math_functions_mod, only: find_all_roots, first_derivative
		! use write_log_mod, only: write_log, log_input

		type(cnt), intent(inout) :: currcnt

		integer :: mu_e,ik_e
		integer :: mu_ph, iq_ph, ib
		real*8, dimension(2) :: k_e, q_ph
		character(len=1000) :: filename
		complex*16, dimension(:,:,:), allocatable :: matrix_element

		allocate(matrix_element(1-currcnt%Nu/2:currcnt%Nu/2,currcnt%ikc_min:currcnt%ikc_max,6))


		mu_e = 0
		ik_e = 0
		k_e = dble(mu_e) * currcnt%K1 + dble(ik_e) * currcnt%dk * currcnt%K2

		do mu_ph=1-currcnt%Nu/2,currcnt%Nu/2
			do iq_ph=currcnt%ikc_min,currcnt%ikc_max
				do ib = 1,6
					q_ph = dble(mu_ph) * currcnt%K1 + dble(iq_ph) * currcnt%dk * currcnt%K2
					call graphene_electron_phonon_matrix_element(matrix_element(mu_ph, iq_ph, ib),k_e,q_ph, ib, currcnt%a1, currcnt%a2)
				enddo
			enddo
		enddo

		!***********************************************************************
		! save the CNT electron-phonon matrix element

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
				write(100,'(E16.8)', advance='no') (abs(matrix_element(mu_ph,iq_ph,1)))**2
				write(101,'(E16.8)', advance='no') (abs(matrix_element(mu_ph,iq_ph,2)))**2
				write(102,'(E16.8)', advance='no') (abs(matrix_element(mu_ph,iq_ph,3)))**2
				write(103,'(E16.8)', advance='no') (abs(matrix_element(mu_ph,iq_ph,4)))**2
				write(104,'(E16.8)', advance='no') (abs(matrix_element(mu_ph,iq_ph,5)))**2
				write(105,'(E16.8)', advance='no') (abs(matrix_element(mu_ph,iq_ph,6)))**2
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
	end subroutine cnt_electron_phonon_matrix_element

	!***************************************************************************
	! -	this subroutine calculates and saves a list of the electronic and
	!	phononic states that conserve both momentum and energy in an
	!	electron-phonon scattering process.
	!***************************************************************************
	subroutine cnt_electron_phonon_scattering_states(currcnt)
		use cnt_class, only: cnt
		use constants_mod
		use graphene_mod, only: graphene_electron
		use math_functions_mod, only: find_all_roots, first_derivative
		use write_log_mod, only: write_log, log_input

		type(cnt), intent(inout) :: currcnt

		integer :: i
		real*8, dimension(:), allocatable :: tmp_real_array_1

		integer :: mu_e,ik_e
		integer :: mu_e_2
		integer :: mu_ph, iq_ph
		integer :: ib
		integer :: number_of_scattered_states
		integer :: n_root
		integer, dimension(:), allocatable :: root_idx
		character(len=1000) :: filename
		integer, dimension(:,:), allocatable :: scattering_state_list

		!***********************************************************************
		! calculate CNT electronic energy dispersion.
		call calculate_cnt_electron_energy_dispersion(currcnt)

		!***********************************************************************
		!find phonon momentum that conserves both momentum and energy

		allocate(tmp_real_array_1(currcnt%ikc_min:currcnt%ikc_max))
		allocate(root_idx(2*currcnt%ikc_max+1))

		!first count the number of scattered states to be able to allocate the need memory for scattered states.
		number_of_scattered_states = 0
		do ib=1,6
			do mu_e=1-currcnt%Nu/2,currcnt%Nu/2
				do mu_e_2=1-currcnt%Nu/2,currcnt%Nu/2
					mu_ph = mu_e-mu_e_2
					do ik_e = currcnt%ikc_min,currcnt%ikc_max

						tmp_real_array_1 = E_k(mu_e_2,:,1)-E_k(mu_e,ik_e,1)+currcnt%omega_phonon(mu_ph,ik_e+currcnt%ikc_max:ik_e+currcnt%ikc_min:-1,ib)

						call find_all_roots (tmp_real_array_1, currcnt%ikc_min, currcnt%ikc_max, 0.d0, n_root, root_idx)
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

						tmp_real_array_1 = E_k(mu_e_2,:,1)-E_k(mu_e,ik_e,1)+currcnt%omega_phonon(mu_ph,ik_e+currcnt%ikc_max:ik_e+currcnt%ikc_min:-1,ib)

						call find_all_roots (tmp_real_array_1, currcnt%ikc_min, currcnt%ikc_max, 0.d0, n_root, root_idx)

						do i= 1, n_root
							iq_ph = ik_e - root_idx(i)
							scattering_state_list(number_of_scattered_states + i,:) = (/ number_of_scattered_states+i, ib, mu_e, mu_e_2, mu_ph, ik_e, root_idx(i), iq_ph   /)
						enddo

						number_of_scattered_states = number_of_scattered_states + n_root

					enddo
				enddo
			enddo
		enddo

		!***********************************************************************
		!save the calculated list of scattering states that conserve both momentum and energy

		write(filename,'(A)') trim(currcnt%name)//".electron_phonon_scattering_states.dat"
		open(unit=100,file=trim(filename),status="unknown")

		do i = 1, number_of_scattered_states
			write(100,'(SP,I10.9, A, I10.9, A, I10.9, A, I10.9, A, I10.9, A, I10.9, A, I10.9, A, I10.9)') scattering_state_list(i,1), '   ', scattering_state_list(i,2), '   ', scattering_state_list(i,3), '   ', scattering_state_list(i,4), '   ', scattering_state_list(i,5), '   ', scattering_state_list(i,6), '   ', scattering_state_list(i,7), '   ', scattering_state_list(i,8)
		enddo

		close(100)

		write(log_input,'(A,I10.10)') "Number of scattered states is : ", number_of_scattered_states
		call write_log(log_input)
	end subroutine cnt_electron_phonon_scattering_states

end module cnt_scattering_electron_phonon_mod
