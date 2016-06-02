!*******************************************************************************
! This program calculates the excitonic band structure of single wall carbon nanotubes through simple tight-binding method.
! Amirhossein Davoody
! Last modified: 3/17/2014
!*******************************************************************************

program cnt_phonon_assisted_energy_transfer
	use cnt_class, only: cnt, free_cnt_memory
	use cnt_electron_mod, only: cnt_electron_band_structure
	use cnt_geometry_mod, only: cnt_geometry
	use cnt_phonon_mod, only: cnt_phonon_dispersion
	use cnt_scattering_electron_phonon_mod, only: cnt_electron_phonon_scattering_rate_emission, cnt_electron_phonon_scattering_rate_absorption, cnt_electron_phonon_matrix_element, cnt_electron_phonon_scattering_states
	use cnt_scattering_exciton_phonon_mod, only: cnt_exction_phonon_scattering_rate_emission, cnt_exction_phonon_scattering_rate_absorption
	use comparams, only: cnt1, cnt2
	use input_cnt_mod, only: input_cnt_parameters, input_exciton
	! use occupation_mod, only: calculate_occupation_table
	! use parse_input_file_mod, only: parse_input_file
	! use prepareForster_module, only: saveDOS
	use sim_properties_mod, only: input_sim_properties, finalize_output_directory_name
	! use transition_table_mod, only: calculate_transition_table
	use write_log_mod, only: write_log, log_input

	implicit none

	integer :: tmp_i, tmp_j

	character(len=1000) :: filename
	real :: start_time, end_time

	call CPU_time(start_time)

	! check the input format for the correct number of input files
	if (command_argument_count() .ne. 3) then
		write(*,*) "Input format ERROR!"
		write(*,*) "Correct input format is: main.exe cnt1_prop.in cnt2_prop.in simulation_prop.in"
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

	call input_exciton(ex_type=1, alpha=0, currcnt=cnt1, exciton_energy_filename='Ex_A1.dat', exciton_wavefunction_filename='Psi_A1.dat')
	call input_exciton(ex_type=2, alpha=0, currcnt=cnt1, exciton_energy_filename='Ex0_A2.dat', exciton_wavefunction_filename='Psi0_A2.dat')
	call input_exciton(ex_type=3, alpha=0, currcnt=cnt1, exciton_energy_filename='Ex0_Ep.dat', exciton_wavefunction_filename='Psi0_Ep.dat')
	call input_exciton(ex_type=4, alpha=0, currcnt=cnt1, exciton_energy_filename='Ex0_Em.dat', exciton_wavefunction_filename='Psi0_Em.dat')

	do tmp_i = 1, 4
		do tmp_j = 1, 4
			call cnt_exction_phonon_scattering_rate_emission(currcnt=cnt1, i_exciton=cnt1%excitons(tmp_i,0), f_exciton=cnt1%excitons(tmp_j,0))
			call cnt_exction_phonon_scattering_rate_absorption(currcnt=cnt1, i_exciton=cnt1%excitons(tmp_i,0), f_exciton=cnt1%excitons(tmp_j,0))
		enddo
	enddo

	call cnt_electron_phonon_scattering_rate_emission(cnt1)
	call cnt_electron_phonon_scattering_rate_absorption(cnt1)
	call cnt_electron_phonon_matrix_element(cnt1)

	write(log_input,'(A)') new_line('A')//"cnt1 data loaded successfuly!!!"
	call write_log(trim(log_input))

	! call exit()
	!
	! write(cnt2%name,'(A)') "cnt2"
	! call get_command_argument(2,filename)
	! call input_cnt_parameters(cnt2,trim(filename))
	! call cnt_geometry(cnt2)
	! call cnt_electron_band_structure(cnt2)
	! call cnt_phonon_dispersion(cnt2)
	! ! call input_a_exciton(cnt2)
	! call cnt_electron_phonon_scattering_rate_emission(cnt2)
	! call cnt_electron_phonon_scattering_rate_absorption(cnt2)
	!
	! write(log_input,'(A)') new_line('A')//"cnt2 data loaded successfuly!!!"
	! call write_log(trim(log_input))



	! call parse_input_file()
	!
	! call write_log(new_line('A')//"************** Reading cnt1 ****************")
	! call input_cnt(cnt1)
	! call saveDOS(cnt1)
	! call calculate_occupation_table(cnt1)
	!
	! call exit()
	!
	! call write_log(new_line('A')//"************** Reading cnt2 ****************")
	! call input_cnt(cnt2)
	! call saveDOS(cnt2)
	! call calculate_occupation_table(cnt2)
	!
	! call exit()
	!
	! call calculate_transition_table(cnt1,cnt2)
	!
	! call calculateKappaMatrix(cnt1,cnt2)

	! save information about simulation runtime
	call CPU_time(end_time)
	write(log_input,'(A,f10.3,A)') new_line('A')//"Run time = ",end_time-start_time," seconds."
	call write_log(log_input)

	! rename the output directory from a temporary name to a final name
	call finalize_output_directory_name()

	! deallocate all allocatable components in cnt_class
	call free_cnt_memory(cnt1)
	call free_cnt_memory(cnt2)

end program cnt_phonon_assisted_energy_transfer
