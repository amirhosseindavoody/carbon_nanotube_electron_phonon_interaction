!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate first order Coulomb mediated transition rates
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module second_order_coulomb_phonon_transition_mod
	implicit none
	private
	public :: calculate_phonon_emission_Coulomb_coupling_perpendicular, calculate_Coulomb_coupling_phonon_emission_perpendicular, calculate_phonon_emission_Coulomb_coupling_parallel, calculate_Coulomb_coupling_phonon_emission_parallel, calculate_phonon_absorption_Coulomb_coupling_parallel

	real*8, private :: energy_mesh_min, energy_mesh_max

contains
	

	!**************************************************************************************************************************
	! calculate transition table using integration over phononic states.
	! the geometry of the tubes is perpendicular (unparallel) in this case.
	! the length of the tubes are infinite.
	! the physical process is phonon emission in the donor CNT then Coulomb coupling to acceptor CNT.
	!**************************************************************************************************************************

	subroutine calculate_phonon_emission_Coulomb_coupling_perpendicular (cnt_1,cnt_2)
		use constants_mod
		use cnt_class, only: cnt, exciton
		use cnt_phonon_mod, only: cnt_phonon_dispersion
		use cnt_scattering_exciton_phonon_mod, only: cnt_exciton_phonon_matrix_element, i_exciton_phonon_scattering_rate
		use geometric_matrix_element_mod, only: calculate_infinite_geometric_matrix_element
		use input_cnt_mod, only: input_exciton, smooth_exciton_dispersion
		use kspace_matrix_element_mod, only: calculate_Q_tilde
		use math_functions_mod, only: find_all_roots, first_derivative
		use sim_properties_mod, only: temperature, c2c_min, energy_mesh_size, energy_mesh_length
		use write_log_mod, only: write_log, log_input

		type(cnt), target, intent(inout) :: cnt_1,cnt_2

		integer :: ix1, ix2, iKcm1, iKcm2
		integer :: ix_m, iKcm_m
		complex*16 :: total_matrix_element, geometric_matrix_element, Q_tilde, coulomb_matrix_element, phonon_matrix_element

		integer :: t_ex_type
		type(exciton), pointer :: i_exciton, m_exciton, f_exciton

		integer :: mu_ph, ib, iq_ph, iq_ph_half

		real*8 :: min_energy, max_energy

		real*8, dimension(:), allocatable :: energy_mesh
		real*8, dimension(:), allocatable :: scattering_rate

		integer :: i, j, k

		real*8, dimension(:), allocatable :: tmp_array_1, tmp_array_2
		integer, dimension(:), allocatable :: initial_state_idx, final_state_idx
		integer :: n_initial_state, n_final_state

		real*8 :: c2c_distance
		character(len=1000) :: filename

		real*8 :: omega_tmp
		real*8 :: f_exciton_energy, f_exciton_energy_derivative
		integer :: counter
		complex*16 :: matrix_element

		real*8 :: tmp_scattering_rate

		real*8 :: partition_function
		real*8 :: dE, i_exciton_energy_derivative
		real*8, dimension(:), allocatable :: occupation

		c2c_distance = c2c_min

		i_exciton => cnt_1%selected_exciton
		f_exciton => cnt_2%selected_exciton
		m_exciton => cnt_1%excitons(2,0)
		if (.not. allocated(m_exciton%ex)) then
			call input_exciton(ex_type=2, alpha=0, currcnt=cnt_1, exciton_energy_filename='Ex0_A2.dat', exciton_wavefunction_filename='Psi0_A2.dat')
		endif

		call smooth_exciton_dispersion(i_exciton, "i_exciton.dispersion.dat")
		call smooth_exciton_dispersion(m_exciton, "m_exciton.dispersion.dat")
		call smooth_exciton_dispersion(f_exciton, "f_exciton.dispersion.dat")

		write(log_input, '(A, A)') "exciton dispersions are smoothed!!!!", new_line('A')
		call write_log(log_input)

		! create the electron energy mesh
		allocate(energy_mesh(energy_mesh_size))

		energy_mesh_min = minval(i_exciton%ex(:,:))
		energy_mesh_max = energy_mesh_min + energy_mesh_length

		do i=1,energy_mesh_size
			energy_mesh(i) = energy_mesh_min + dble(i-1)*(energy_mesh_max-energy_mesh_min)/dble(energy_mesh_size-1)
		enddo

		allocate(scattering_rate(energy_mesh_size))
		scattering_rate = 0.d0

		!allocate some temporary arrays that are used for finding scattering states
		if (allocated(tmp_array_1)) deallocate(tmp_array_1)
		allocate(tmp_array_1(i_exciton%iKcm_min:i_exciton%iKcm_max))

		if (allocated(tmp_array_2)) deallocate(tmp_array_2)
		allocate(tmp_array_2(f_exciton%iKcm_min:f_exciton%iKcm_max))

		if(allocated(initial_state_idx)) deallocate(initial_state_idx)
		allocate(initial_state_idx(i_exciton%iKcm_max-i_exciton%iKcm_min+1))

		if (allocated(final_state_idx)) deallocate(final_state_idx)
		allocate(final_state_idx(f_exciton%iKcm_max-f_exciton%iKcm_min+1))

		mu_ph = 2*(i_exciton%mu_cm-m_exciton%mu_cm)
		call cnt_phonon_dispersion(cnt_1, dq=2.d0*i_exciton%dKcm, iq_max=i_exciton%iKcm_max, iq_min=i_exciton%iKcm_min, mu_max=mu_ph, mu_min=mu_ph, save_dispersion=.true. )

		! calculate occupation number
		partition_function = 0.d0
		dE = energy_mesh(2) - energy_mesh(1)
		allocate(occupation(energy_mesh_size))
		occupation = 0.d0
		do i = 1, energy_mesh_size
			do ix1 = 1, 1 !i_exciton%nx
				tmp_array_1 = i_exciton%ex(ix1,:)
				call find_all_roots (tmp_array_1, lbound(tmp_array_1,dim=1), ubound(tmp_array_1, dim=1), energy_mesh(i), n_initial_state, initial_state_idx)

				do j=1, n_initial_state
					iKcm1 = initial_state_idx(j)
					if (iKcm1 .eq. 0) then
						call first_derivative(tmp_array_1, lbound(tmp_array_1,dim=1), ubound(tmp_array_1, dim=1), iKcm1+2, i_exciton%dKcm, i_exciton_energy_derivative, single=.true.)
					else
						call first_derivative(tmp_array_1, lbound(tmp_array_1,dim=1), ubound(tmp_array_1, dim=1), iKcm1, i_exciton%dKcm, i_exciton_energy_derivative, single=.true.)
					endif
					partition_function = partition_function + dE / abs(i_exciton_energy_derivative) * exp(-i_exciton%ex(ix1,iKcm1)/kb/temperature)
					! partition_function = partition_function + dE * exp(-i_exciton%ex(ix1,iKcm1)/kb/temperature)
				enddo
			enddo
		enddo

		do i = 1, energy_mesh_size
			do ix1 = 1, 1 !i_exciton%nx
				tmp_array_1 = i_exciton%ex(ix1,:)
				call find_all_roots (tmp_array_1, lbound(tmp_array_1,dim=1), ubound(tmp_array_1, dim=1), energy_mesh(i), n_initial_state, initial_state_idx)
				do j=1, n_initial_state
					iKcm1 = initial_state_idx(j)
					
					if (iKcm1 .eq. 0) then
						call first_derivative(tmp_array_1, lbound(tmp_array_1,dim=1), ubound(tmp_array_1, dim=1), iKcm1+2, i_exciton%dKcm, i_exciton_energy_derivative, single=.true.)
					else
						call first_derivative(tmp_array_1, lbound(tmp_array_1,dim=1), ubound(tmp_array_1, dim=1), iKcm1, i_exciton%dKcm, i_exciton_energy_derivative, single=.true.)
					endif
					! occupation(i) = occupation(i) + 1.d0 / abs(i_exciton_energy_derivative) * exp(-i_exciton%ex(ix1,iKcm1)/kb/temperature) / partition_function
					! occupation(i) = occupation(i) + dE * exp(-i_exciton%ex(ix1,iKcm1)/kb/temperature) / partition_function
					occupation(i) = occupation(i) + dE / abs(i_exciton_energy_derivative) * exp(-i_exciton%ex(ix1,iKcm1)/kb/temperature) / partition_function
					! write(*,*) iKcm1, i_exciton_energy_derivative
					! occupation(i) = occupation(i) + 1 / i_exciton_energy_derivative
				enddo

			enddo

			write(log_input, '(A, A, I0, A, I0, A, E10.4, A, E10.4)') trim(i_exciton%name), "   ", i, "/", energy_mesh_size, "  energy=", energy_mesh(i)/eV, "   occupation = ",occupation(i)
			call write_log(log_input)

		enddo

		write(*,*) "sum of occupations = ", sum(occupation)

		call exit()


		do i = 1, energy_mesh_size

			do ix1 = 1, 2 !i_exciton%nx
				do ix2= 1,2 !f_exciton%nx
			
					do ib= 1,6

				
						tmp_array_1 = i_exciton%ex(ix1,:)
						call find_all_roots (tmp_array_1, lbound(tmp_array_1,dim=1), ubound(tmp_array_1, dim=1), energy_mesh(i), n_initial_state, initial_state_idx)

						do j=1, n_initial_state
							iKcm1 = initial_state_idx(j)



							do iq_ph_half = lbound(cnt_1%omega_phonon, dim=2), ubound(cnt_1%omega_phonon, dim=2)
								omega_tmp = cnt_1%omega_phonon(mu_ph,iq_ph_half,ib)
								iq_ph = 2*iq_ph_half
								iKcm_m = iKcm1-iq_ph_half
								if ((iKcm_m .ge. m_exciton%iKcm_min) .and. (iKcm_m .le. m_exciton%iKcm_max)) then
									f_exciton_energy = i_exciton%ex(ix1,iKcm1)-cnt_1%omega_phonon(mu_ph,iq_ph_half,ib)
									tmp_array_2 = f_exciton%ex(ix2,:)
									call find_all_roots (tmp_array_2, lbound(tmp_array_2,dim=1), ubound(tmp_array_2, dim=1), f_exciton_energy, n_final_state, final_state_idx)

									do k=1, n_final_state
										iKcm2 = final_state_idx(k)

										total_matrix_element = (0.d0, 0.d0)
										do ix_m = 1,2 !m_exciton%nx

											call cnt_exciton_phonon_matrix_element(phonon_matrix_element, cnt_1, i_exciton, m_exciton, ix1, iKcm1, ix_m, iKcm_m, ib )
											call calculate_Q_tilde(Q_tilde, cnt_1, cnt_2, m_exciton, f_exciton, ix_m, iKcm_m, ix2, iKcm2)
											call calculate_infinite_geometric_matrix_element(geometric_matrix_element, cnt_1, cnt_2, m_exciton%mu_cm, iKcm_m, f_exciton%mu_cm, iKcm2, pi/2.d0, c2c_distance)
											coulomb_matrix_element = Q_tilde * geometric_matrix_element * (A_u**2) /(4.d0*(pi**2)*cnt_1%radius*cnt_2%radius)
											matrix_element = phonon_matrix_element * coulomb_matrix_element / (dcmplx(f_exciton_energy-m_exciton%ex(ix_m,iKcm_m))+ i1*dcmplx(0.03d0*eV))
											! matrix_element = phonon_matrix_element * coulomb_matrix_element / (dcmplx(f_exciton_energy-m_exciton%ex(ix_m,iKcm_m))+ i1*dcmplx(i_exciton_phonon_scattering_rate(i)*hb*2.d0*pi/2.0d0)) 
											total_matrix_element = total_matrix_element + matrix_element

										enddo

										tmp_array_2 = f_exciton%ex(ix2,:)
										call first_derivative(tmp_array_2, lbound(tmp_array_2,dim=1), ubound(tmp_array_2, dim=1), iKcm2, f_exciton%dKcm, f_exciton_energy_derivative)
										if (f_exciton_energy_derivative .eq. 0.d0) then
											write(log_input, '(A)') "final exciton energy derivative is zero!!"
											call write_log(log_input)
											call first_derivative(tmp_array_2, lbound(tmp_array_2,dim=1), ubound(tmp_array_2, dim=1), iKcm2+1, f_exciton%dKcm, f_exciton_energy_derivative)
										endif

										tmp_scattering_rate = 2.d0*i_exciton%dKcm / (2*pi*hb) *(1.d0+1.d0/(exp(omega_tmp/temperature/kb)-1.d0)) * ((abs(total_matrix_element))**2)/(abs(f_exciton_energy_derivative))

										! if (((i.eq. 11) .and. (tmp_scattering_rate .ge. 1.d+0)) .or. ((i .eq. 20) .and. (tmp_scattering_rate .ge. 1.d+1))) then

											! write(log_input, '(A, E8.2, A, E8.2, A, E8.2, A, I0, A, I0, A, I0, A, I0, A, I0, A, I0, A, I0, A, I0, A, I0, A, I0, A, I0)') &
											!  		"scattering rate=", tmp_scattering_rate, &
											!  		"   matrix_element=", abs(total_matrix_element), &
											!  		"   f_exciton_energy_derivative=", abs(f_exciton_energy_derivative), &
											!  		"   k=", k, "/", n_final_state, &
											!  		"   ib=", ib, &
											!  		"   iq_ph_half=", iq_ph_half, "/", ubound(cnt_1%omega_phonon, dim=2), &
											!  		"   iKcm_m=", iKcm_m, "/", ubound(m_exciton%ex, dim=2), &
											!  		"   iKcm1=", iKcm1, "/", ubound(i_exciton%ex, dim=2), &
											!  		"   iKcm2=", iKcm2, "/", ubound(f_exciton%ex, dim=2)
											! call write_log(log_input)

											scattering_rate(i) = scattering_rate(i) + tmp_scattering_rate

										! endif

									enddo

								endif

							enddo

						enddo
					enddo
				enddo

			enddo

			write(log_input, '(A, A, A, A, I0, A, I0, A, E10.4, A, E10.4)') trim(i_exciton%name), " to ", trim(f_exciton%name), " phonon_emission_coulomb_coupling_perpendicular (", i, "/", energy_mesh_size, ")  energy=", energy_mesh(i)/eV, " rate=",scattering_rate(i)
			call write_log(log_input)

		enddo

		!***********************************************************************
		!save the calculated exciton-phonon scattering rate

		write(filename,'(A, A, A, A, A)') "phonon_emission_coulomb_coupling_perpendicular.", trim(i_exciton%name), "_to_", trim(f_exciton%name), ".dat"
		open(unit=100,file=trim(filename),status="unknown")

		do i = 1, energy_mesh_size
			write(100,'(SP, E20.8)', advance='no') energy_mesh(i)
		enddo

		write(100,*)
		do i = 1, energy_mesh_size
			write(100,'(SP, E20.8)', advance='no') scattering_rate(i)
		enddo

		close(100)

	end subroutine calculate_phonon_emission_Coulomb_coupling_perpendicular



	!**************************************************************************************************************************
	! calculate transition table by integration over phononic states
	! the geometry of the tubes is perpendicular (unparallel) in this case.
	! the length of the tubes are infinite.
	! the physical process is Coulomb coupling then phonon emission in the acceptor CNT.
	!**************************************************************************************************************************

	subroutine calculate_Coulomb_coupling_phonon_emission_perpendicular (cnt_1,cnt_2)
		use constants_mod
		use cnt_class, only: cnt, exciton
		use cnt_phonon_mod, only: cnt_phonon_dispersion
		use cnt_scattering_exciton_phonon_mod, only: cnt_exciton_phonon_matrix_element, i_exciton_phonon_scattering_rate
		use geometric_matrix_element_mod, only: calculate_infinite_geometric_matrix_element
		use input_cnt_mod, only: input_exciton, smooth_exciton_dispersion
		use kspace_matrix_element_mod, only: calculate_Q_tilde
		use math_functions_mod, only: find_all_roots, first_derivative
		use sim_properties_mod, only: temperature, c2c_min, energy_mesh_size, energy_mesh_length
		use write_log_mod, only: write_log, log_input

		type(cnt), target, intent(inout) :: cnt_1,cnt_2

		integer :: ix1, ix2, iKcm1, iKcm2
		integer :: ix_m, iKcm_m
		complex*16 :: total_matrix_element, geometric_matrix_element, Q_tilde, coulomb_matrix_element, phonon_matrix_element

		integer :: t_ex_type
		type(exciton), pointer :: i_exciton, m_exciton, f_exciton

		integer :: mu_ph, ib, iq_ph, iq_ph_half

		real*8 :: min_energy, max_energy

		real*8, dimension(:), allocatable :: energy_mesh
		real*8, dimension(:), allocatable :: scattering_rate

		integer :: i, j, k

		real*8, dimension(:), allocatable :: tmp_array_1, tmp_array_2
		integer, dimension(:), allocatable :: initial_state_idx, final_state_idx
		integer :: n_initial_state, n_final_state

		real*8 :: c2c_distance
		character(len=1000) :: filename

		real*8 :: omega_tmp
		real*8 :: f_exciton_energy, f_exciton_energy_derivative
		integer :: counter
		complex*16 :: matrix_element

		real*8 :: tmp_scattering_rate

		c2c_distance = c2c_min

		i_exciton => cnt_1%selected_exciton
		f_exciton => cnt_2%selected_exciton
		m_exciton => cnt_2%excitons(2,0)
		if (.not. allocated(m_exciton%ex)) then
			call input_exciton(ex_type=2, alpha=0, currcnt=cnt_2, exciton_energy_filename='Ex0_A2.dat', exciton_wavefunction_filename='Psi0_A2.dat')
		endif

		call smooth_exciton_dispersion(i_exciton, "i_exciton.dispersion.dat")
		call smooth_exciton_dispersion(m_exciton, "m_exciton.dispersion.dat")
		call smooth_exciton_dispersion(f_exciton, "f_exciton.dispersion.dat")

		write(log_input, '(A, A)') "exciton dispersions are smoothed!!!!", new_line('A')
		call write_log(log_input)

		! create the electron energy mesh
		allocate(energy_mesh(energy_mesh_size))

		energy_mesh_min = minval(i_exciton%ex(:,:))
		energy_mesh_max = energy_mesh_min + energy_mesh_length

		do i=1,energy_mesh_size
			energy_mesh(i) = energy_mesh_min + dble(i-1)*(energy_mesh_max-energy_mesh_min)/dble(energy_mesh_size-1)
		enddo

		allocate(scattering_rate(energy_mesh_size))
		scattering_rate = 0.d0

		!allocate some temporary arrays that are used for finding scattering states
		if (allocated(tmp_array_1)) deallocate(tmp_array_1)
		allocate(tmp_array_1(i_exciton%iKcm_min:i_exciton%iKcm_max))

		if (allocated(tmp_array_2)) deallocate(tmp_array_2)
		allocate(tmp_array_2(f_exciton%iKcm_min:f_exciton%iKcm_max))

		if(allocated(initial_state_idx)) deallocate(initial_state_idx)
		allocate(initial_state_idx(i_exciton%iKcm_max-i_exciton%iKcm_min+1))

		if (allocated(final_state_idx)) deallocate(final_state_idx)
		allocate(final_state_idx(f_exciton%iKcm_max-f_exciton%iKcm_min+1))


		mu_ph = 2*(m_exciton%mu_cm-f_exciton%mu_cm)
		call cnt_phonon_dispersion(cnt_2, dq=2.d0*f_exciton%dKcm, iq_max=f_exciton%iKcm_max, iq_min=f_exciton%iKcm_min, mu_max=mu_ph, mu_min=mu_ph, save_dispersion=.true. )


		do i = 1, energy_mesh_size

			do ix1 = 1, 2 !i_exciton%nx
				do ix2= 1,2 !f_exciton%nx
			
					do ib= 1,6

				
						tmp_array_1 = i_exciton%ex(ix1,:)
						call find_all_roots (tmp_array_1, lbound(tmp_array_1,dim=1), ubound(tmp_array_1, dim=1), energy_mesh(i), n_initial_state, initial_state_idx)

						do j=1, n_initial_state
							iKcm1 = initial_state_idx(j)



							do iq_ph_half = lbound(cnt_2%omega_phonon, dim=2), ubound(cnt_2%omega_phonon, dim=2)
								omega_tmp = cnt_2%omega_phonon(mu_ph,iq_ph_half,ib)
								iq_ph = 2*iq_ph_half

									
								f_exciton_energy = i_exciton%ex(ix1,iKcm1)-cnt_2%omega_phonon(mu_ph,iq_ph_half,ib)
								tmp_array_2 = f_exciton%ex(ix2,:)
								call find_all_roots (tmp_array_2, lbound(tmp_array_2,dim=1), ubound(tmp_array_2, dim=1), f_exciton_energy, n_final_state, final_state_idx)

								do k=1, n_final_state
									iKcm2 = final_state_idx(k)
									iKcm_m = iKcm2+iq_ph_half


									if ((iKcm_m .ge. m_exciton%iKcm_min) .and. (iKcm_m .le. m_exciton%iKcm_max)) then
										total_matrix_element = (0.d0, 0.d0)
										do ix_m = 1,2 !m_exciton%nx

											call cnt_exciton_phonon_matrix_element(phonon_matrix_element, cnt_2, m_exciton, f_exciton, ix_m, iKcm_m, ix2, iKcm2, ib )
											call calculate_Q_tilde(Q_tilde, cnt_1, cnt_2, i_exciton, m_exciton, ix1, iKcm1, ix_m, iKcm_m)
											call calculate_infinite_geometric_matrix_element(geometric_matrix_element, cnt_1, cnt_2, i_exciton%mu_cm, iKcm1, m_exciton%mu_cm, iKcm_m, pi/2.d0, c2c_distance)
											coulomb_matrix_element = Q_tilde * geometric_matrix_element * (A_u**2) /(4.d0*(pi**2)*cnt_1%radius*cnt_2%radius)
											matrix_element = phonon_matrix_element * coulomb_matrix_element / (dcmplx(i_exciton%ex(ix1,iKcm1)-m_exciton%ex(ix_m,iKcm_m))+ i1*dcmplx(0.03d0*eV))
											! matrix_element = phonon_matrix_element * coulomb_matrix_element / (dcmplx(i_exciton%ex(ix1,iKcm1)-m_exciton%ex(ix_m,iKcm_m))+ i1*dcmplx(i_exciton_phonon_scattering_rate(i)*hb*2.d0*pi/2.0d0))
											total_matrix_element = total_matrix_element + matrix_element

										enddo

										tmp_array_2 = f_exciton%ex(ix2,:)
										
										if ((tmp_array_2(iKcm2+1)-tmp_array_2(iKcm2))*(tmp_array_2(iKcm2)-tmp_array_2(iKcm2-1)) .lt. 0) then
											iKcm2 = iKcm2 + 1
										endif

										call first_derivative(tmp_array_2, lbound(tmp_array_2,dim=1), ubound(tmp_array_2, dim=1), iKcm2, f_exciton%dKcm, f_exciton_energy_derivative)

										tmp_scattering_rate = 2.d0*i_exciton%dKcm / (2*pi*hb) *(1.d0+1.d0/(exp(omega_tmp/temperature/kb)-1.d0)) * ((abs(total_matrix_element))**2)/(abs(f_exciton_energy_derivative))

										! if (((i.eq. 11) .and. (tmp_scattering_rate .ge. 1.d0)) .or. ((i .eq. 12) .and. (tmp_scattering_rate .ge. 1.d0))) then

										! 	write(log_input, '(A, E8.2, A, E8.2, A, E8.2, A, I0, A, I0, A, I0, A, I0, A, I0, A, I0, A, I0, A, I0, A, I0, A, I0, A, I0)') &
										! 	 		"scattering rate=", tmp_scattering_rate, &
										! 	 		"   matrix_element=", abs(total_matrix_element), &
										! 	 		! "   phonon_matrix_element=", abs(phonon_matrix_element), &
										! 	 		! "   Q_tilde=", abs(Q_tilde), &
										! 	 		! "   geometric_matrix_element=", abs(geometric_matrix_element), &
										! 	 		"   f_exciton_energy_derivative=", abs(f_exciton_energy_derivative), &
										! 	 		"   k=", k, "/", n_final_state, &
										! 	 		"   ib=", ib, &
										! 	 		"   iq_ph_half=", iq_ph_half, "/", ubound(cnt_2%omega_phonon, dim=2), &
										! 	 		"   iKcm_m=", iKcm_m, "/", ubound(m_exciton%ex, dim=2), &
										! 	 		"   iKcm1=", iKcm1, "/", ubound(i_exciton%ex, dim=2), &
										! 	 		"   iKcm2=", iKcm2, "/", ubound(f_exciton%ex, dim=2)
										! 	call write_log(log_input)

											scattering_rate(i) = scattering_rate(i) + tmp_scattering_rate

										! endif

									endif

								enddo

							enddo

						enddo
					enddo
				enddo

			enddo

			write(log_input, '(A, A, A, A, I0, A, I0, A, E10.4, A, E10.4)') trim(i_exciton%name), " to ", trim(f_exciton%name), " coulomb_coupling_phonon_emission_perpendicular (", i, "/", energy_mesh_size, ")  energy=", energy_mesh(i)/eV, " rate=",scattering_rate(i)
			call write_log(log_input)

		enddo

		!***********************************************************************
		!save the calculated exciton-phonon scattering rate

		write(filename,'(A, A, A, A, A)') "coulomb_coupling_phonon_emission_perpendicular.", trim(i_exciton%name), "_to_", trim(f_exciton%name), ".dat"
		open(unit=100,file=trim(filename),status="unknown")

		do i = 1, energy_mesh_size
			write(100,'(SP, E20.8)', advance='no') energy_mesh(i)
		enddo

		write(100,*)
		do i = 1, energy_mesh_size
			write(100,'(SP, E20.8)', advance='no') scattering_rate(i)
		enddo

		close(100)

	end subroutine calculate_Coulomb_coupling_phonon_emission_perpendicular


	!**************************************************************************************************************************
	! calculate transition table using integration over phononic states.
	! the geometry of the tubes is parallel in this case.
	! the length of the tubes are infinite.
	! the physical process is phonon emission in the donor CNT then Coulomb coupling to acceptor CNT.
	!**************************************************************************************************************************

	subroutine calculate_phonon_emission_Coulomb_coupling_parallel (cnt_1,cnt_2)
		use constants_mod
		use cnt_class, only: cnt, exciton
		use cnt_phonon_mod, only: cnt_phonon_dispersion
		use cnt_scattering_exciton_phonon_mod, only: cnt_exciton_phonon_matrix_element, i_exciton_phonon_scattering_rate
		use geometric_matrix_element_mod, only: calculate_infinite_parallel_geometric_matrix_element
		use input_cnt_mod, only: input_exciton, smooth_exciton_dispersion
		use kspace_matrix_element_mod, only: calculate_Q_tilde
		use math_functions_mod, only: find_all_roots, first_derivative
		use sim_properties_mod, only: temperature, c2c_min, energy_mesh_size, energy_mesh_length
		use write_log_mod, only: write_log, log_input

		type(cnt), target, intent(inout) :: cnt_1,cnt_2

		integer :: ix1, ix2, iKcm1, iKcm2
		integer :: ix_m, iKcm_m
		complex*16 :: total_matrix_element, geometric_matrix_element, Q_tilde, coulomb_matrix_element, phonon_matrix_element

		integer :: t_ex_type
		type(exciton), pointer :: i_exciton, m_exciton, f_exciton

		integer :: mu_ph, ib, iq_ph, iq_ph_half
		integer :: iq_ph_min, iq_ph_max

		real*8 :: min_energy, max_energy

		real*8, dimension(:), allocatable :: energy_mesh
		real*8, dimension(:), allocatable :: scattering_rate

		integer :: i, j, k

		real*8, dimension(:), allocatable :: tmp_array_1, tmp_array_2, tmp_array_3
		integer, dimension(:), allocatable :: initial_state_idx, final_state_idx
		integer :: n_initial_state, n_final_state

		real*8 :: c2c_distance
		character(len=1000) :: filename

		real*8 :: omega_tmp
		real*8 :: f_exciton_energy, f_exciton_energy_derivative
		integer :: counter
		complex*16 :: matrix_element

		real*8 :: tmp_scattering_rate

		c2c_distance = c2c_min

		i_exciton => cnt_1%selected_exciton
		f_exciton => cnt_2%selected_exciton
		m_exciton => cnt_1%excitons(2,0)
		if (.not. allocated(m_exciton%ex)) then
			call input_exciton(ex_type=2, alpha=0, currcnt=cnt_1, exciton_energy_filename='Ex0_A2.dat', exciton_wavefunction_filename='Psi0_A2.dat')
		endif

		call smooth_exciton_dispersion(i_exciton, "i_exciton.dispersion.dat")
		call smooth_exciton_dispersion(m_exciton, "m_exciton.dispersion.dat")
		call smooth_exciton_dispersion(f_exciton, "f_exciton.dispersion.dat")

		write(log_input, '(A, A)') "exciton dispersions are smoothed!!!!", new_line('A')
		call write_log(log_input)

		! create the electron energy mesh
		allocate(energy_mesh(energy_mesh_size))

		energy_mesh_min = minval(i_exciton%ex(:,:))
		energy_mesh_max = energy_mesh_min + energy_mesh_length

		do i=1,energy_mesh_size
			energy_mesh(i) = energy_mesh_min + dble(i-1)*(energy_mesh_max-energy_mesh_min)/dble(energy_mesh_size-1)
		enddo

		allocate(scattering_rate(energy_mesh_size))
		scattering_rate = 0.d0

		!allocate some temporary arrays that are used for finding scattering states
		if (allocated(tmp_array_1)) deallocate(tmp_array_1)
		allocate(tmp_array_1(i_exciton%iKcm_min:i_exciton%iKcm_max))

		if (allocated(tmp_array_3)) deallocate(tmp_array_3)
		allocate(tmp_array_3(f_exciton%iKcm_min:f_exciton%iKcm_max))

		if(allocated(initial_state_idx)) deallocate(initial_state_idx)
		allocate(initial_state_idx(i_exciton%iKcm_max-i_exciton%iKcm_min+1))

		mu_ph = 2*(i_exciton%mu_cm-m_exciton%mu_cm)
		call cnt_phonon_dispersion(cnt_1, dq=2.d0*i_exciton%dKcm, iq_max=i_exciton%iKcm_max, iq_min=i_exciton%iKcm_min, mu_max=mu_ph, mu_min=mu_ph, save_dispersion=.true. )

		do i = 1, energy_mesh_size

			do ix1 = 1, 10 !i_exciton%nx
				do ix2= 1,10 !f_exciton%nx
			
					do ib= 1,6

				
						tmp_array_1 = i_exciton%ex(ix1,:)
						call find_all_roots (tmp_array_1, lbound(tmp_array_1,dim=1), ubound(tmp_array_1, dim=1), energy_mesh(i), n_initial_state, initial_state_idx)

						do j=1, n_initial_state
							iKcm1 = initial_state_idx(j)

							iq_ph_min = max((iKcm1-f_exciton%iKcm_max),lbound(cnt_1%omega_phonon, dim=2))
							iq_ph_max = min((iKcm1-f_exciton%iKcm_min),ubound(cnt_1%omega_phonon, dim=2))

							! tmp_array_2 = 0.d0
							if (allocated(tmp_array_2)) deallocate(tmp_array_2)
							allocate(tmp_array_2(iq_ph_min:iq_ph_max))
							if (allocated(final_state_idx)) deallocate(final_state_idx)
							allocate(final_state_idx(iq_ph_max-iq_ph_min+1))
							tmp_array_2 = i_exciton%ex(ix1,iKcm1) - cnt_1%omega_phonon(mu_ph,iq_ph_min:iq_ph_max,ib) - f_exciton%ex(ix2,iKcm1-iq_ph_min:iKcm1-iq_ph_max:-1)
							call find_all_roots (tmp_array_2, lbound(tmp_array_2,dim=1), ubound(tmp_array_2, dim=1), 0.d0, n_final_state, final_state_idx)

							do k = 1, n_final_state


								iq_ph_half = final_state_idx(k)
								omega_tmp = cnt_1%omega_phonon(mu_ph,iq_ph_half,ib)
								iq_ph = 2*iq_ph_half
								iKcm_m = iKcm1-iq_ph_half
								iKcm2 = iKcm_m
								if ( ((iKcm_m .ge. m_exciton%iKcm_min) .and. (iKcm_m .le. m_exciton%iKcm_max)) .and. ((iKcm2 .ge. f_exciton%iKcm_min) .and. (iKcm2 .le. f_exciton%iKcm_max)) )then
									f_exciton_energy = i_exciton%ex(ix1,iKcm1)-cnt_1%omega_phonon(mu_ph,iq_ph_half,ib)

									total_matrix_element = (0.d0, 0.d0)
									do ix_m = 1,10 !m_exciton%nx

										call cnt_exciton_phonon_matrix_element(phonon_matrix_element, cnt_1, i_exciton, m_exciton, ix1, iKcm1, ix_m, iKcm_m, ib )
										call calculate_Q_tilde(Q_tilde, cnt_1, cnt_2, m_exciton, f_exciton, ix_m, iKcm_m, ix2, iKcm2)
										call calculate_infinite_parallel_geometric_matrix_element(geometric_matrix_element, cnt_1, cnt_2, m_exciton%mu_cm, iKcm_m, f_exciton%mu_cm, iKcm2, c2c_distance)
										coulomb_matrix_element = Q_tilde * geometric_matrix_element * (A_u**2) /(4.d0*(pi**2)*cnt_1%radius*cnt_2%radius)
										matrix_element = phonon_matrix_element * coulomb_matrix_element / (dcmplx(f_exciton_energy-m_exciton%ex(ix_m,iKcm_m))+ i1*dcmplx(0.03d0*eV))
										! matrix_element = phonon_matrix_element * coulomb_matrix_element / (dcmplx(f_exciton_energy-m_exciton%ex(ix_m,iKcm_m))+ i1*dcmplx(i_exciton_phonon_scattering_rate(i)*hb*2.d0*pi/2.0d0))
										total_matrix_element = total_matrix_element + matrix_element

									enddo

									tmp_array_3 = f_exciton%ex(ix2,:)
									call first_derivative(tmp_array_3, lbound(tmp_array_3,dim=1), ubound(tmp_array_3, dim=1), iKcm2, f_exciton%dKcm, f_exciton_energy_derivative)
									if (f_exciton_energy_derivative .eq. 0.d0) then
										write(log_input, '(A)') "final exciton energy derivative is zero!!"
										call write_log(log_input)
										call first_derivative(tmp_array_3, lbound(tmp_array_3,dim=1), ubound(tmp_array_3, dim=1), iKcm2+1, f_exciton%dKcm, f_exciton_energy_derivative)
									endif

									tmp_scattering_rate = (1/hb) *(1.d0+1.d0/(exp(omega_tmp/temperature/kb)-1.d0)) * ((abs(total_matrix_element))**2)/(abs(f_exciton_energy_derivative)/2.d0)

									! if (((i.eq. 19) .and. (tmp_scattering_rate .ge. 1.d+0)) .or. ((i .eq. 20) .and. (tmp_scattering_rate .ge. 1.d+1))) then

										! write(log_input, '(A, E8.2, A, E8.2, A, E8.2, A, I0, A, I0, A, I0, A, I0, A, I0, A, I0, A, I0, A, I0, A, I0, A, I0, A, I0)') &
										!  		"scattering rate=", tmp_scattering_rate, &
										!  		"   matrix_element=", abs(total_matrix_element), &
										!  		"   f_exciton_energy_derivative=", abs(f_exciton_energy_derivative), &
										!  		"   k=", k, "/", n_final_state, &
										!  		"   ib=", ib, &
										!  		"   iq_ph_half=", iq_ph_half, "/", ubound(cnt_1%omega_phonon, dim=2), &
										!  		"   iKcm_m=", iKcm_m, "/", ubound(m_exciton%ex, dim=2), &
										!  		"   iKcm1=", iKcm1, "/", ubound(i_exciton%ex, dim=2), &
										!  		"   iKcm2=", iKcm2, "/", ubound(f_exciton%ex, dim=2)
										! call write_log(log_input)

										scattering_rate(i) = scattering_rate(i) + tmp_scattering_rate

									! endif

								endif

							enddo

						enddo
					enddo
				enddo

			enddo

			write(log_input, '(A, A, A, A, I0, A, I0, A, E10.4, A, E10.4)') trim(i_exciton%name), " to ", trim(f_exciton%name), " phonon_emission_coulomb_coupling_parallel (", i, "/", energy_mesh_size, ")  energy=", energy_mesh(i)/eV, " rate=",scattering_rate(i)
			call write_log(log_input)

		enddo

		!***********************************************************************
		!save the calculated exciton-phonon scattering rate

		write(filename,'(A, A, A, A, A)') "phonon_emission_coulomb_coupling_parallel.", trim(i_exciton%name), "_to_", trim(f_exciton%name), ".dat"
		open(unit=100,file=trim(filename),status="unknown")

		do i = 1, energy_mesh_size
			write(100,'(SP, E20.8)', advance='no') energy_mesh(i)
		enddo

		write(100,*)
		do i = 1, energy_mesh_size
			write(100,'(SP, E20.8)', advance='no') scattering_rate(i)
		enddo

		close(100)

	end subroutine calculate_phonon_emission_Coulomb_coupling_parallel


	!**************************************************************************************************************************
	! calculate transition table using integration over phononic states.
	! the geometry of the tubes is parallel in this case.
	! the length of the tubes are infinite.
	! the physical process is Coulomb coupling then phonon emission in the acceptor CNT.
	!**************************************************************************************************************************

	subroutine calculate_Coulomb_coupling_phonon_emission_parallel (cnt_1,cnt_2)
		use constants_mod
		use cnt_class, only: cnt, exciton
		use cnt_phonon_mod, only: cnt_phonon_dispersion
		use cnt_scattering_exciton_phonon_mod, only: cnt_exciton_phonon_matrix_element, i_exciton_phonon_scattering_rate
		use geometric_matrix_element_mod, only: calculate_infinite_parallel_geometric_matrix_element
		use input_cnt_mod, only: input_exciton, smooth_exciton_dispersion
		use kspace_matrix_element_mod, only: calculate_Q_tilde
		use math_functions_mod, only: find_all_roots, first_derivative
		use sim_properties_mod, only: temperature, c2c_min, energy_mesh_size, energy_mesh_length
		use write_log_mod, only: write_log, log_input

		type(cnt), target, intent(inout) :: cnt_1,cnt_2

		integer :: ix1, ix2, iKcm1, iKcm2
		integer :: ix_m, iKcm_m
		complex*16 :: total_matrix_element, geometric_matrix_element, Q_tilde, coulomb_matrix_element, phonon_matrix_element

		integer :: t_ex_type
		type(exciton), pointer :: i_exciton, m_exciton, f_exciton

		integer :: mu_ph, ib, iq_ph, iq_ph_half
		integer :: iq_ph_min, iq_ph_max

		real*8 :: min_energy, max_energy

		real*8, dimension(:), allocatable :: energy_mesh
		real*8, dimension(:), allocatable :: scattering_rate

		integer :: i, j, k

		real*8, dimension(:), allocatable :: tmp_array_1, tmp_array_2, tmp_array_3
		integer, dimension(:), allocatable :: initial_state_idx, final_state_idx
		integer :: n_initial_state, n_final_state

		real*8 :: c2c_distance
		character(len=1000) :: filename

		real*8 :: omega_tmp
		real*8 :: f_exciton_energy, f_exciton_energy_derivative
		integer :: counter
		complex*16 :: matrix_element

		real*8 :: tmp_scattering_rate

		c2c_distance = c2c_min

		i_exciton => cnt_1%selected_exciton
		f_exciton => cnt_2%selected_exciton
		m_exciton => cnt_2%excitons(2,0)
		if (.not. allocated(m_exciton%ex)) then
			call input_exciton(ex_type=2, alpha=0, currcnt=cnt_2, exciton_energy_filename='Ex0_A2.dat', exciton_wavefunction_filename='Psi0_A2.dat')
		endif

		call smooth_exciton_dispersion(i_exciton, "i_exciton.dispersion.dat")
		call smooth_exciton_dispersion(m_exciton, "m_exciton.dispersion.dat")
		call smooth_exciton_dispersion(f_exciton, "f_exciton.dispersion.dat")

		write(log_input, '(A, A)') "exciton dispersions are smoothed!!!!", new_line('A')
		call write_log(log_input)

		! create the electron energy mesh
		allocate(energy_mesh(energy_mesh_size))

		energy_mesh_min = minval(i_exciton%ex(:,:))
		energy_mesh_max = energy_mesh_min + energy_mesh_length

		do i=1,energy_mesh_size
			energy_mesh(i) = energy_mesh_min + dble(i-1)*(energy_mesh_max-energy_mesh_min)/dble(energy_mesh_size-1)
		enddo

		allocate(scattering_rate(energy_mesh_size))
		scattering_rate = 0.d0

		!allocate some temporary arrays that are used for finding scattering states
		if (allocated(tmp_array_1)) deallocate(tmp_array_1)
		allocate(tmp_array_1(i_exciton%iKcm_min:i_exciton%iKcm_max))

		if (allocated(tmp_array_3)) deallocate(tmp_array_3)
		allocate(tmp_array_3(f_exciton%iKcm_min:f_exciton%iKcm_max))

		if(allocated(initial_state_idx)) deallocate(initial_state_idx)
		allocate(initial_state_idx(i_exciton%iKcm_max-i_exciton%iKcm_min+1))

		mu_ph = 2*(m_exciton%mu_cm-f_exciton%mu_cm)
		call cnt_phonon_dispersion(cnt_2, dq=2.d0*i_exciton%dKcm, iq_max=i_exciton%iKcm_max, iq_min=i_exciton%iKcm_min, mu_max=mu_ph, mu_min=mu_ph, save_dispersion=.true. )

		do i = 1, energy_mesh_size

			do ix1 = 1, 2 !i_exciton%nx
				do ix2= 1,2 !f_exciton%nx
			
					do ib= 1,6

				
						tmp_array_1 = i_exciton%ex(ix1,:)
						call find_all_roots (tmp_array_1, lbound(tmp_array_1,dim=1), ubound(tmp_array_1, dim=1), energy_mesh(i), n_initial_state, initial_state_idx)

						do j=1, n_initial_state
							iKcm1 = initial_state_idx(j)

							iq_ph_min = max((iKcm1-f_exciton%iKcm_max),lbound(cnt_2%omega_phonon, dim=2))
							iq_ph_max = min((iKcm1-f_exciton%iKcm_min),ubound(cnt_2%omega_phonon, dim=2))

							if (allocated(tmp_array_2)) deallocate(tmp_array_2)
							allocate(tmp_array_2(iq_ph_min:iq_ph_max))
							if (allocated(final_state_idx)) deallocate(final_state_idx)
							allocate(final_state_idx(iq_ph_max-iq_ph_min+1))
							tmp_array_2 = i_exciton%ex(ix1,iKcm1) - cnt_2%omega_phonon(mu_ph,iq_ph_min:iq_ph_max,ib) - f_exciton%ex(ix2,iKcm1-iq_ph_min:iKcm1-iq_ph_max:-1)
							call find_all_roots (tmp_array_2, lbound(tmp_array_2,dim=1), ubound(tmp_array_2, dim=1), 0.d0, n_final_state, final_state_idx)

							do k = 1, n_final_state


								iq_ph_half = final_state_idx(k)
								omega_tmp = cnt_2%omega_phonon(mu_ph,iq_ph_half,ib)
								iq_ph = 2*iq_ph_half
								iKcm_m = iKcm1
								iKcm2 = iKcm1-iq_ph_half
								if ( ((iKcm_m .ge. m_exciton%iKcm_min) .and. (iKcm_m .le. m_exciton%iKcm_max)) .and. ((iKcm2 .ge. f_exciton%iKcm_min) .and. (iKcm2 .le. f_exciton%iKcm_max)) )then
									f_exciton_energy = i_exciton%ex(ix1,iKcm1)-cnt_2%omega_phonon(mu_ph,iq_ph_half,ib)

									total_matrix_element = (0.d0, 0.d0)
									do ix_m = 1,2 !m_exciton%nx

										call cnt_exciton_phonon_matrix_element(phonon_matrix_element, cnt_2, m_exciton, f_exciton, ix_m, iKcm_m, ix2, iKcm2, ib )
										call calculate_Q_tilde(Q_tilde, cnt_1, cnt_2, i_exciton, m_exciton, ix1, iKcm1, ix_m, iKcm_m)
										call calculate_infinite_parallel_geometric_matrix_element(geometric_matrix_element, cnt_1, cnt_2, i_exciton%mu_cm, iKcm1, m_exciton%mu_cm, iKcm_m, c2c_distance)
										coulomb_matrix_element = Q_tilde * geometric_matrix_element * (A_u**2) /(4.d0*(pi**2)*cnt_1%radius*cnt_2%radius)
										matrix_element = phonon_matrix_element * coulomb_matrix_element / (dcmplx(i_exciton%ex(ix1,iKcm1)-m_exciton%ex(ix_m,iKcm_m))+ i1*dcmplx(0.03d0*eV))
										! matrix_element = phonon_matrix_element * coulomb_matrix_element / (dcmplx(i_exciton%ex(ix1,iKcm1)-m_exciton%ex(ix_m,iKcm_m))+ i1*dcmplx(i_exciton_phonon_scattering_rate(i)*hb*2.d0*pi/2.0d0))
										total_matrix_element = total_matrix_element + matrix_element

									enddo

									tmp_array_3 = f_exciton%ex(ix2,:)
									call first_derivative(tmp_array_3, lbound(tmp_array_3,dim=1), ubound(tmp_array_3, dim=1), iKcm2, f_exciton%dKcm, f_exciton_energy_derivative)
									if (f_exciton_energy_derivative .eq. 0.d0) then
										write(log_input, '(A)') "final exciton energy derivative is zero!!"
										call write_log(log_input)
										call first_derivative(tmp_array_3, lbound(tmp_array_3,dim=1), ubound(tmp_array_3, dim=1), iKcm2+1, f_exciton%dKcm, f_exciton_energy_derivative)
									endif

									tmp_scattering_rate = (1/hb) *(1.d0+1.d0/(exp(omega_tmp/temperature/kb)-1.d0)) * ((abs(total_matrix_element))**2)/(abs(f_exciton_energy_derivative)/2.d0)

									! if (((i.eq. 19) .and. (tmp_scattering_rate .ge. 1.d+0)) .or. ((i .eq. 20) .and. (tmp_scattering_rate .ge. 1.d+1))) then

										! write(log_input, '(A, E8.2, A, E8.2, A, E8.2, A, I0, A, I0, A, I0, A, I0, A, I0, A, I0, A, I0, A, I0, A, I0, A, I0, A, I0)') &
										!  		"scattering rate=", tmp_scattering_rate, &
										!  		"   matrix_element=", abs(total_matrix_element), &
										!  		"   f_exciton_energy_derivative=", abs(f_exciton_energy_derivative), &
										!  		"   k=", k, "/", n_final_state, &
										!  		"   ib=", ib, &
										!  		"   iq_ph_half=", iq_ph_half, "/", ubound(cnt_2%omega_phonon, dim=2), &
										!  		"   iKcm_m=", iKcm_m, "/", ubound(m_exciton%ex, dim=2), &
										!  		"   iKcm1=", iKcm1, "/", ubound(i_exciton%ex, dim=2), &
										!  		"   iKcm2=", iKcm2, "/", ubound(f_exciton%ex, dim=2)
										! call write_log(log_input)

										scattering_rate(i) = scattering_rate(i) + tmp_scattering_rate

									! endif

								endif

							enddo

						enddo
					enddo
				enddo

			enddo

			write(log_input, '(A, A, A, A, I0, A, I0, A, E10.4, A, E10.4)') trim(i_exciton%name), " to ", trim(f_exciton%name), " coulomb_coupling_phonon_emission_parallel (", i, "/", energy_mesh_size, ")  energy=", energy_mesh(i)/eV, " rate=",scattering_rate(i)
			call write_log(log_input)

		enddo

		!***********************************************************************
		!save the calculated exciton-phonon scattering rate

		write(filename,'(A, A, A, A, A)') "coulomb_coupling_phonon_emission_parallel.", trim(i_exciton%name), "_to_", trim(f_exciton%name), ".dat"
		open(unit=100,file=trim(filename),status="unknown")

		do i = 1, energy_mesh_size
			write(100,'(SP, E20.8)', advance='no') energy_mesh(i)
		enddo

		write(100,*)
		do i = 1, energy_mesh_size
			write(100,'(SP, E20.8)', advance='no') scattering_rate(i)
		enddo

		close(100)

	end subroutine calculate_Coulomb_coupling_phonon_emission_parallel


	!**************************************************************************************************************************
	! calculate transition table using integration over phononic states.
	! the geometry of the tubes is parallel in this case.
	! the length of the tubes are infinite.
	! the physical process is phonon absorption in the donor CNT then Coulomb coupling to acceptor CNT.
	!**************************************************************************************************************************

	subroutine calculate_phonon_absorption_Coulomb_coupling_parallel (cnt_1,cnt_2)
		use constants_mod
		use cnt_class, only: cnt, exciton
		use cnt_phonon_mod, only: cnt_phonon_dispersion
		use cnt_scattering_exciton_phonon_mod, only: cnt_exciton_phonon_matrix_element, i_exciton_phonon_scattering_rate
		use geometric_matrix_element_mod, only: calculate_infinite_parallel_geometric_matrix_element
		use input_cnt_mod, only: input_exciton, smooth_exciton_dispersion
		use kspace_matrix_element_mod, only: calculate_Q_tilde
		use math_functions_mod, only: find_all_roots, first_derivative
		use sim_properties_mod, only: temperature, c2c_min, energy_mesh_size, energy_mesh_length
		use write_log_mod, only: write_log, log_input

		type(cnt), target, intent(inout) :: cnt_1,cnt_2

		integer :: ix1, ix2, iKcm1, iKcm2
		integer :: ix_m, iKcm_m
		complex*16 :: total_matrix_element, geometric_matrix_element, Q_tilde, coulomb_matrix_element, phonon_matrix_element

		integer :: t_ex_type
		type(exciton), pointer :: i_exciton, m_exciton, f_exciton

		integer :: mu_ph, ib, iq_ph, iq_ph_half
		integer :: iq_ph_min, iq_ph_max

		real*8 :: min_energy, max_energy

		real*8, dimension(:), allocatable :: energy_mesh
		real*8, dimension(:), allocatable :: scattering_rate

		integer :: i, j, k

		real*8, dimension(:), allocatable :: tmp_array_1, tmp_array_2, tmp_array_3
		integer, dimension(:), allocatable :: initial_state_idx, final_state_idx
		integer :: n_initial_state, n_final_state

		real*8 :: c2c_distance
		character(len=1000) :: filename

		real*8 :: omega_tmp
		real*8 :: f_exciton_energy, f_exciton_energy_derivative
		integer :: counter
		complex*16 :: matrix_element

		real*8 :: tmp_scattering_rate

		c2c_distance = c2c_min

		i_exciton => cnt_1%selected_exciton
		f_exciton => cnt_2%selected_exciton
		m_exciton => cnt_1%excitons(2,0)
		if (.not. allocated(m_exciton%ex)) then
			call input_exciton(ex_type=2, alpha=0, currcnt=cnt_1, exciton_energy_filename='Ex0_A2.dat', exciton_wavefunction_filename='Psi0_A2.dat')
		endif

		call smooth_exciton_dispersion(i_exciton, "i_exciton.dispersion.dat")
		call smooth_exciton_dispersion(m_exciton, "m_exciton.dispersion.dat")
		call smooth_exciton_dispersion(f_exciton, "f_exciton.dispersion.dat")

		write(log_input, '(A, A)') "exciton dispersions are smoothed!!!!", new_line('A')
		call write_log(log_input)

		! create the electron energy mesh
		allocate(energy_mesh(energy_mesh_size))

		energy_mesh_min = minval(i_exciton%ex(:,:))
		energy_mesh_max = energy_mesh_min + energy_mesh_length

		do i=1,energy_mesh_size
			energy_mesh(i) = energy_mesh_min + dble(i-1)*(energy_mesh_max-energy_mesh_min)/dble(energy_mesh_size-1)
		enddo

		allocate(scattering_rate(energy_mesh_size))
		scattering_rate = 0.d0

		!allocate some temporary arrays that are used for finding scattering states
		if (allocated(tmp_array_1)) deallocate(tmp_array_1)
		allocate(tmp_array_1(i_exciton%iKcm_min:i_exciton%iKcm_max))

		if (allocated(tmp_array_3)) deallocate(tmp_array_3)
		allocate(tmp_array_3(f_exciton%iKcm_min:f_exciton%iKcm_max))

		if(allocated(initial_state_idx)) deallocate(initial_state_idx)
		allocate(initial_state_idx(i_exciton%iKcm_max-i_exciton%iKcm_min+1))

		mu_ph = 2*(i_exciton%mu_cm-m_exciton%mu_cm)
		call cnt_phonon_dispersion(cnt_1, dq=2.d0*i_exciton%dKcm, iq_max=i_exciton%iKcm_max, iq_min=i_exciton%iKcm_min, mu_max=mu_ph, mu_min=mu_ph, save_dispersion=.true. )

		do i = 1, energy_mesh_size

			do ix1 = 1, 2 !i_exciton%nx
				do ix2= 1,2 !f_exciton%nx
			
					do ib= 1,6

				
						tmp_array_1 = i_exciton%ex(ix1,:)
						call find_all_roots (tmp_array_1, lbound(tmp_array_1,dim=1), ubound(tmp_array_1, dim=1), energy_mesh(i), n_initial_state, initial_state_idx)

						do j=1, n_initial_state
							iKcm1 = initial_state_idx(j)

							iq_ph_min = max((iKcm1-f_exciton%iKcm_max),lbound(cnt_1%omega_phonon, dim=2))
							iq_ph_max = min((iKcm1-f_exciton%iKcm_min),ubound(cnt_1%omega_phonon, dim=2))

							! tmp_array_2 = 0.d0
							if (allocated(tmp_array_2)) deallocate(tmp_array_2)
							allocate(tmp_array_2(iq_ph_min:iq_ph_max))
							if (allocated(final_state_idx)) deallocate(final_state_idx)
							allocate(final_state_idx(iq_ph_max-iq_ph_min+1))
							tmp_array_2 = i_exciton%ex(ix1,iKcm1) + cnt_1%omega_phonon(mu_ph,iq_ph_min:iq_ph_max,ib) - f_exciton%ex(ix2,iKcm1-iq_ph_min:iKcm1-iq_ph_max:-1)
							call find_all_roots (tmp_array_2, lbound(tmp_array_2,dim=1), ubound(tmp_array_2, dim=1), 0.d0, n_final_state, final_state_idx)

							do k = 1, n_final_state


								iq_ph_half = final_state_idx(k)
								omega_tmp = cnt_1%omega_phonon(mu_ph,iq_ph_half,ib)
								iq_ph = 2*iq_ph_half
								iKcm_m = iKcm1-iq_ph_half
								iKcm2 = iKcm_m
								if ( ((iKcm_m .ge. m_exciton%iKcm_min) .and. (iKcm_m .le. m_exciton%iKcm_max)) .and. ((iKcm2 .ge. f_exciton%iKcm_min) .and. (iKcm2 .le. f_exciton%iKcm_max)) )then
									f_exciton_energy = i_exciton%ex(ix1,iKcm1)+cnt_1%omega_phonon(mu_ph,iq_ph_half,ib)

									total_matrix_element = (0.d0, 0.d0)
									do ix_m = 1,2 !m_exciton%nx

										call cnt_exciton_phonon_matrix_element(phonon_matrix_element, cnt_1, i_exciton, m_exciton, ix1, iKcm1, ix_m, iKcm_m, ib )
										call calculate_Q_tilde(Q_tilde, cnt_1, cnt_2, m_exciton, f_exciton, ix_m, iKcm_m, ix2, iKcm2)
										call calculate_infinite_parallel_geometric_matrix_element(geometric_matrix_element, cnt_1, cnt_2, m_exciton%mu_cm, iKcm_m, f_exciton%mu_cm, iKcm2, c2c_distance)
										coulomb_matrix_element = Q_tilde * geometric_matrix_element * (A_u**2) /(4.d0*(pi**2)*cnt_1%radius*cnt_2%radius)
										matrix_element = phonon_matrix_element * coulomb_matrix_element / (dcmplx(f_exciton_energy-m_exciton%ex(ix_m,iKcm_m))+ i1*dcmplx(0.03d0*eV))
										! matrix_element = phonon_matrix_element * coulomb_matrix_element / (dcmplx(f_exciton_energy-m_exciton%ex(ix_m,iKcm_m))+ i1*dcmplx(i_exciton_phonon_scattering_rate(i)*hb*2.d0*pi/2.0d0))
										total_matrix_element = total_matrix_element + matrix_element

									enddo

									tmp_array_3 = f_exciton%ex(ix2,:)
									call first_derivative(tmp_array_3, lbound(tmp_array_3,dim=1), ubound(tmp_array_3, dim=1), iKcm2, f_exciton%dKcm, f_exciton_energy_derivative)
									if (f_exciton_energy_derivative .eq. 0.d0) then
										write(log_input, '(A)') "final exciton energy derivative is zero!!"
										call write_log(log_input)
										call first_derivative(tmp_array_3, lbound(tmp_array_3,dim=1), ubound(tmp_array_3, dim=1), iKcm2+1, f_exciton%dKcm, f_exciton_energy_derivative)
									endif

									tmp_scattering_rate = (1/hb) *(1.d0/(exp(omega_tmp/temperature/kb)-1.d0)) * ((abs(total_matrix_element))**2)/(abs(f_exciton_energy_derivative)/2.d0)

									! if (((i.eq. 19) .and. (tmp_scattering_rate .ge. 1.d+0)) .or. ((i .eq. 20) .and. (tmp_scattering_rate .ge. 1.d+1))) then

										! write(log_input, '(A, E8.2, A, E8.2, A, E8.2, A, I0, A, I0, A, I0, A, I0, A, I0, A, I0, A, I0, A, I0, A, I0, A, I0, A, I0)') &
										!  		"scattering rate=", tmp_scattering_rate, &
										!  		"   matrix_element=", abs(total_matrix_element), &
										!  		"   f_exciton_energy_derivative=", abs(f_exciton_energy_derivative), &
										!  		"   k=", k, "/", n_final_state, &
										!  		"   ib=", ib, &
										!  		"   iq_ph_half=", iq_ph_half, "/", ubound(cnt_1%omega_phonon, dim=2), &
										!  		"   iKcm_m=", iKcm_m, "/", ubound(m_exciton%ex, dim=2), &
										!  		"   iKcm1=", iKcm1, "/", ubound(i_exciton%ex, dim=2), &
										!  		"   iKcm2=", iKcm2, "/", ubound(f_exciton%ex, dim=2)
										! call write_log(log_input)

										scattering_rate(i) = scattering_rate(i) + tmp_scattering_rate

									! endif

								endif

							enddo

						enddo
					enddo
				enddo

			enddo

			write(log_input, '(A, A, A, A, I0, A, I0, A, E10.4, A, E10.4)') trim(i_exciton%name), " to ", trim(f_exciton%name), " phonon_absorption_coulomb_coupling_parallel (", i, "/", energy_mesh_size, ")  energy=", energy_mesh(i)/eV, " rate=",scattering_rate(i)
			call write_log(log_input)

		enddo

		!***********************************************************************
		!save the calculated exciton-phonon scattering rate

		write(filename,'(A, A, A, A, A)') "phonon_absorption_coulomb_coupling_parallel.", trim(i_exciton%name), "_to_", trim(f_exciton%name), ".dat"
		open(unit=100,file=trim(filename),status="unknown")

		do i = 1, energy_mesh_size
			write(100,'(SP, E20.8)', advance='no') energy_mesh(i)
		enddo

		write(100,*)
		do i = 1, energy_mesh_size
			write(100,'(SP, E20.8)', advance='no') scattering_rate(i)
		enddo

		close(100)

	end subroutine calculate_phonon_absorption_Coulomb_coupling_parallel
	

end module second_order_coulomb_phonon_transition_mod
