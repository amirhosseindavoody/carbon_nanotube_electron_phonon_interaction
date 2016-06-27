!*******************************************************************************
! This program calculates the excitonic band structure of single wall carbon nanotubes through simple tight-binding method.
! Amirhossein Davoody
! Last modified: 3/17/2014
!*******************************************************************************

program cnt_phonon_assisted_energy_transfer
	use cnt_class, only: cnt, cnt1, cnt2, free_cnt_memory
	use cnt_electron_mod, only: cnt_electron_band_structure
	use cnt_geometry_mod, only: cnt_geometry
	use cnt_phonon_mod, only: cnt_phonon_dispersion
	use cnt_scattering_electron_phonon_mod, only: cnt_electron_phonon_scattering_rate_emission, cnt_electron_phonon_scattering_rate_absorption, cnt_electron_phonon_matrix_element, cnt_electron_phonon_scattering_states
	use cnt_scattering_exciton_phonon_mod, only: cnt_exction_phonon_scattering_rate_emission, cnt_exction_phonon_scattering_rate_absorption
	use first_order_coulomb_transition_mod, only: calculate_first_order_transition_rates
	use input_cnt_mod, only: input_cnt_parameters, input_exciton
	use partition_function_mod, only: calculate_partition_function
	use second_order_coulomb_phonon_transition_mod, only: calculate_second_order_transition_rates
	use sim_properties_mod, only: input_sim_properties, finalize_output_directory_name
	use write_log_mod, only: write_log, log_input

	implicit none

	character(len=1000) :: filename
	real :: start_time, end_time

	call CPU_time(start_time)

	! check the input format for the correct number of input files
	if (command_argument_count() .ne. 3) then
		write(*,'(A)') "Input format ERROR!"
		write(*,'(A)') "Correct input format is: main.exe cnt1_prop.in cnt2_prop.in simulation_prop.in"
		call exit()
	end if

	call get_command_argument(3,filename)
	call input_sim_properties(trim(filename))

	write(cnt1%name,'(A)') "cnt1"
	call get_command_argument(1,filename)
	call input_cnt_parameters(cnt1,trim(filename))
	call cnt_geometry(cnt1)
	call cnt_electron_band_structure(cnt1)
	call cnt_phonon_dispersion(cnt1, save_dispersion=.true.)

	call input_exciton(ex_type=2, alpha=0, currcnt=cnt1, exciton_energy_filename='Ex0_A2.dat', exciton_wavefunction_filename='Psi0_A2.dat')

	write(log_input,'(A)') new_line('A')//"cnt1 data loaded successfuly!!!"
	call write_log(trim(log_input))

	write(cnt2%name,'(A)') "cnt2"
	call get_command_argument(2,filename)
	call input_cnt_parameters(cnt2,trim(filename))
	call cnt_geometry(cnt2)
	call cnt_electron_band_structure(cnt2)
	call cnt_phonon_dispersion(cnt2, save_dispersion=.true.)

	call input_exciton(ex_type=2, alpha=0, currcnt=cnt2, exciton_energy_filename='Ex0_A2.dat', exciton_wavefunction_filename='Psi0_A2.dat')

	write(log_input,'(A)') new_line('A')//"cnt2 data loaded successfuly!!!"
	call write_log(trim(log_input))

	cnt1%selected_exciton => cnt1%excitons(2,0)
	cnt2%selected_exciton => cnt2%excitons(2,0)

	! call calculate_first_order_transition_rates(cnt1, cnt2)
	call calculate_second_order_transition_rates(cnt1, cnt2)

	! save information about simulation runtime
	call CPU_time(end_time)
	write(log_input,'(A,f10.3,A)') new_line('A')//"Run time = ",end_time-start_time," seconds."
	call write_log(log_input)

	! rename the output directory from a temporary name to a final name
	call finalize_output_directory_name(cnt1, cnt2)

	! deallocate all allocatable components in cnt_class
	call free_cnt_memory(cnt1)
	call free_cnt_memory(cnt2)

end program cnt_phonon_assisted_energy_transfer
