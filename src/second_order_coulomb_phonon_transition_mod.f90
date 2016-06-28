!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate first order Coulomb mediated transition rates
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module second_order_coulomb_phonon_transition_mod
	implicit none
	private
	public :: calculate_second_order_transition_rates

	real*8, private :: energy_mesh_min, energy_mesh_max
	real*8, private :: energy_mesh_length = 0.7d0 ! this is the energy distance between energy_mesh_max and energy_mesh_min
	integer, private :: energy_mesh_size = 100

contains
	!**************************************************************************************************************************
	! calculate transition table
	!**************************************************************************************************************************

	subroutine calculate_second_order_transition_rates (cnt_1,cnt_2)
		use constants_mod
		use cnt_class, only: cnt, exciton
		use cnt_phonon_mod, only: cnt_phonon_dispersion
		use cnt_scattering_exciton_phonon_mod, only: cnt_exciton_phonon_matrix_element
		use geometric_matrix_element_mod, only: calculate_infinite_geometric_matrix_element
		use kspace_matrix_element_mod, only: calculate_Q_tilde
		use math_functions_mod, only: find_all_roots, first_derivative
		use sim_properties_mod, only: temperature, c2c_min
		use write_log_mod, only: write_log, log_input

		type(cnt), target, intent(inout) :: cnt_1,cnt_2

		integer :: ix1, ix2, iKcm1, iKcm2
		integer :: ix_m, iKcm_m
		complex*16 :: total_matrix_element, geometric_matrix_element, Q_tilde, coulomb_matrix_element, phonon_matrix_element

		integer :: t_ex_type
		type(exciton), pointer :: i_exciton, m_exciton, f_exciton

		integer :: mu_ph, ib, iq_ph

		real*8 :: min_energy, max_energy

		real*8, dimension(:), allocatable :: energy_mesh
		real*8, dimension(:), allocatable :: scattering_rate

		integer :: i, j, k

		real*8, dimension(:), allocatable :: tmp_array_1, tmp_array_2
		integer, dimension(:), allocatable :: initial_state_idx, mid_state_idx
		integer :: n_initial_state, n_mid_state

		real*8 :: c2c_distance
		real*8 :: omega_phonon_derivative, omega_tmp
		character(len=1000) :: filename

		c2c_distance = c2c_min

		i_exciton => cnt_1%selected_exciton
		f_exciton => cnt_2%selected_exciton

		! create the electron energy mesh
		allocate(energy_mesh(energy_mesh_size))

		energy_mesh_min = minval(i_exciton%ex(:,:))
		energy_mesh_length = (-1.d0) * log(1.d-3) * kb * temperature
		energy_mesh_max = energy_mesh_min + energy_mesh_length

		do i=1,energy_mesh_size
			energy_mesh(i) = energy_mesh_min + dble(i-1)*(energy_mesh_max-energy_mesh_min)/dble(energy_mesh_size-1)
		enddo

		write(log_input, '(A, E10.3)') "energy_mesh_length = ", energy_mesh_length/eV
		call write_log(log_input)

		allocate(scattering_rate(energy_mesh_size))
		scattering_rate = 0.d0

		!allocate some temporary arrays that are used for finding scattering states
		if (allocated(tmp_array_1)) deallocate(tmp_array_1)
		allocate(tmp_array_1(i_exciton%iKcm_min:i_exciton%iKcm_max))

		if (allocated(tmp_array_2)) deallocate(tmp_array_2)
		allocate(tmp_array_2(2*i_exciton%iKcm_min:2*i_exciton%iKcm_max))

		if(allocated(initial_state_idx)) deallocate(initial_state_idx)
		allocate(initial_state_idx(i_exciton%iKcm_max-i_exciton%iKcm_min+1))

		if (allocated(mid_state_idx)) deallocate(mid_state_idx)
		allocate(mid_state_idx(2*i_exciton%iKcm_max-2*i_exciton%iKcm_min+1))

		do t_ex_type = 2, 2
			m_exciton => cnt_1%excitons(t_ex_type,0)

			mu_ph = 2*(i_exciton%mu_cm-m_exciton%mu_cm)
			call cnt_phonon_dispersion(cnt_1, dq=2.d0*i_exciton%dKcm, iq_max=2*i_exciton%iKcm_max, iq_min=2*i_exciton%iKcm_min, mu_max=mu_ph, mu_min=mu_ph )

			min_energy = minval(cnt_1%omega_phonon)
			max_energy = maxval(cnt_1%omega_phonon)

			do i = 1, energy_mesh_size

				do ix1 = 1, i_exciton%nx
					tmp_array_1 = i_exciton%ex(ix1,:)
					call find_all_roots (tmp_array_1, lbound(tmp_array_1,dim=1), ubound(tmp_array_1, dim=1), energy_mesh(i), n_initial_state, initial_state_idx)

					do j=1,n_initial_state
						iKcm1 = initial_state_idx(j)

						do ix2 = 1,f_exciton%nx
							do iKcm2 = f_exciton%iKcm_min, f_exciton%iKcm_max
								if (((i_exciton%ex(ix1,iKcm1)-f_exciton%ex(ix2,iKcm2)) .le. max_energy) .and. ((i_exciton%ex(ix1,iKcm1)-f_exciton%ex(ix2,iKcm2)) .ge. min_energy)) then

									do ib = 1,6
										tmp_array_2 = cnt_1%omega_phonon(mu_ph,:,ib)
										call find_all_roots(tmp_array_2, lbound(tmp_array_2, dim=1), ubound(tmp_array_2, dim=1), i_exciton%ex(ix1,iKcm1)-f_exciton%ex(ix2,iKcm2), n_mid_state, mid_state_idx)
										do k=1,n_mid_state
											iq_ph = 2*mid_state_idx(k)
											iKcm_m = iKcm1-iq_ph/2

											if ((iKcm_m .ge. m_exciton%iKcm_min) .and. (iKcm_m .le. m_exciton%iKcm_max)) then

												total_matrix_element = (0.d0, 0.d0)
												do ix_m = 1,20 !m_exciton%nx
													call cnt_exciton_phonon_matrix_element(phonon_matrix_element, cnt_1, i_exciton, m_exciton, ix1, iKcm1, ix_m, iKcm_m, ib )
													call calculate_Q_tilde(Q_tilde, cnt_1, cnt_2, m_exciton, f_exciton, ix_m, iKcm_m, ix2, iKcm2)
													call calculate_infinite_geometric_matrix_element(geometric_matrix_element, cnt_1, cnt_2, m_exciton%mu_cm, iKcm_m, f_exciton%mu_cm, iKcm2, pi/2.d0, c2c_distance)
													coulomb_matrix_element = Q_tilde * geometric_matrix_element * (A_u**2) /(4.d0*(pi**2)*cnt_1%radius*cnt_2%radius)
													total_matrix_element = total_matrix_element + phonon_matrix_element * coulomb_matrix_element / dcmplx(i_exciton%ex(ix1,iKcm1)-m_exciton%ex(ix_m,iKcm_m)-cnt_1%omega_phonon(mu_ph,mid_state_idx(k),ib))
												enddo

												omega_tmp = cnt_1%omega_phonon(mu_ph, mid_state_idx(k), ib)

												tmp_array_2 = cnt_1%omega_phonon(mu_ph, :, ib)
												call first_derivative(tmp_array_2, lbound(tmp_array_2,dim=1), ubound(tmp_array_2, dim=1), mid_state_idx(k), 2.d0*i_exciton%dKcm, omega_phonon_derivative)
												if (omega_phonon_derivative .eq. 0.d0) then
													call first_derivative(tmp_array_2, lbound(tmp_array_2,dim=1), ubound(tmp_array_2, dim=1), mid_state_idx(k)+1, 2.d0*i_exciton%dKcm, omega_phonon_derivative)
												endif
												scattering_rate(i) = scattering_rate(i) + (1.d0+1.d0/(exp(omega_tmp/temperature/kb)-1.d0)) * ((abs(total_matrix_element))**2)/(abs(omega_phonon_derivative))

											endif

										enddo
									enddo

								endif
							enddo
						enddo

					enddo

				enddo

				scattering_rate(i) = scattering_rate(i) * f_exciton%dKcm / (2*pi*hb)

				write(log_input, '(A,I0,A,I0,A,E8.3)') "second-order exciton scattering rate(", i, "/", energy_mesh_size, ")=", scattering_rate(i)
				call write_log(log_input)
			enddo

		enddo

		!***********************************************************************
		!save the calculated exciton-phonon scattering rate

		write(filename,'(A, A, A, A, A)') "phonon_assisted_scattering_rate_emission.", trim(i_exciton%name), "_to_", trim(f_exciton%name), ".dat"
		open(unit=100,file=trim(filename),status="unknown")

		do i = 1, energy_mesh_size
			write(100,'(SP, E20.8)', advance='no') energy_mesh(i)
		enddo

		write(100,*)
		do i = 1, energy_mesh_size
			write(100,'(SP, E20.8)', advance='no') scattering_rate(i)
		enddo

		close(100)

		write(log_input,'(A)') "phonon-assisted scattering rates due to phonon emission calculated!"
		call write_log(log_input)

	end subroutine calculate_second_order_transition_rates

end module second_order_coulomb_phonon_transition_mod
