!*******************************************************************************
! This program calculates the excitonic band structure of single wall carbon nanotubes through simple tight-binding method.
! Amirhossein Davoody
! Last modified: 3/17/2014
!*******************************************************************************

program cnt_resonance_energy_transfer
	use cnt_class, only: cnt, free_cnt_memory
	use cnt_electron_mod, only: cnt_electron_band_structure
	use cnt_geometry_mod, only: cnt_geometry
	use cnt_phonon_mod, only: cnt_phonon_dispersion
	use comparams, only: cnt1, cnt2
	use input_cnt_mod, only: input_cnt_parameters
	! use occupation_mod, only: calculate_occupation_table
	! use parse_input_file_mod, only: parse_input_file
	! use prepareForster_module, only: saveDOS
	use sim_properties_mod, only: input_sim_properties
	! use transition_table_mod, only: calculate_transition_table
	use write_log_mod, only: write_log, log_input

	implicit none

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

	call get_command_argument(1,filename)
	call input_cnt_parameters(cnt1,trim(filename))

	call cnt_geometry(cnt1)
	call cnt_electron_band_structure(cnt1)
	call cnt_phonon_dispersion(cnt1)

	! don't go beyond this point. program not tested!!!!
	! call exit()

	! call parse_input_file()

	! call write_log(new_line('A')//"************** Reading cnt1 ****************")
	! call input_cnt(cnt1)
! 	call saveDOS(cnt1)
! 	call calculate_occupation_table(cnt1)

! 	call exit()

	! call write_log(new_line('A')//"************** Reading cnt2 ****************")
	! call input_cnt(cnt2)
! 	call saveDOS(cnt2)
! 	call calculate_occupation_table(cnt2)

!  	call exit()

 ! 	call calculate_transition_table(cnt1,cnt2)

! 	call calculateKappaMatrix(cnt1,cnt2)

	call CPU_time(end_time)
	write(log_input,'("Run time = ",f10.3," seconds.")'),end_time-start_time
	call write_log(log_input)

	! deallocate all allocatable components in cnt_class
	call free_cnt_memory(cnt1)
	call free_cnt_memory(cnt2)

end program cnt_resonance_energy_transfer
