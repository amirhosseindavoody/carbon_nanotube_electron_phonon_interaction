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
	use cnt_scattering_electron_phonon_mod, only: cnt_scattering_electron_phonon
	use comparams, only: cnt1, cnt2
	use input_cnt_mod, only: input_cnt_parameters, input_a_exciton
	! use occupation_mod, only: calculate_occupation_table
	! use parse_input_file_mod, only: parse_input_file
	! use prepareForster_module, only: saveDOS
	use sim_properties_mod, only: input_sim_properties
	! use transition_table_mod, only: calculate_transition_table
	use write_log_mod, only: write_log, log_input

	implicit none

	character(len=1000) :: filename
	real :: start_time, end_time

	! complex*16, dimension(3) :: a_vec, b_vec, c_vec
	!
	!
	!
	! a_vec = (/(1.d0,-0.d0), (2.d0,-0.d0), (3.d0,-0.d0) /)
	! ! a_vec = (/(1.d0,-1.d0), (2.d0,-2.d0), (3.d0,-3.d0) /)
	! b_vec = (/(1.d1,-1.d1), (1.d2,-1.d2), (1.d3,-1.d3) /)
	! c_vec = (/(0.d0,-0.d0), (0.d0,-0.d0), (0.d0,-0.d0) /)
	!
	! write(*,'(A,E11.3,E11.3,E11.3,E11.3,E11.3,E11.3)') new_line('A')//'a_vec = ', a_vec
	! write(*,'(A,E11.3,E11.3,E11.3,E11.3,E11.3,E11.3)') new_line('A')//'b_vec = ', b_vec
	! ! write(*,'(A,E11.3,E11.3,E11.3,E11.3,E11.3,E11.3)') new_line('A')//'c_vec = ', c_vec
	!
	! c_vec = a_vec * b_vec
	! ! write(*,'(A,E11.3,E11.3,E11.3,E11.3,E11.3,E11.3)') new_line('A')//'a_vec = ', a_vec
	! ! write(*,'(A,E11.3,E11.3,E11.3,E11.3,E11.3,E11.3)') new_line('A')//'b_vec = ', b_vec
	! write(*,'(A,E11.3,E11.3,E11.3,E11.3,E11.3,E11.3)') new_line('A')//'c_vec = ', c_vec
	!
	! call exit()

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
	call cnt_phonon_dispersion(cnt1)
	! call input_a_exciton(cnt1)
	call cnt_scattering_electron_phonon(cnt1)

	write(log_input,'(A)') new_line('A')//"cnt1 data loaded successfuly!!!"
	call write_log(trim(log_input))

	call exit()

	call get_command_argument(2,filename)
	call input_cnt_parameters(cnt2,trim(filename))
	call cnt_geometry(cnt2)
	call cnt_electron_band_structure(cnt2)
	call cnt_phonon_dispersion(cnt2)
	call input_a_exciton(cnt2)
	write(cnt2%name,'(A)') "cnt2"

	write(log_input,'(A)') new_line('A')//"cnt2 data loaded successfuly!!!"
	call write_log(trim(log_input))



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
	write(log_input,'(A,f10.3,A)') new_line('A')//"Run time = ",end_time-start_time," seconds."
	call write_log(log_input)

	! deallocate all allocatable components in cnt_class
	call free_cnt_memory(cnt1)
	call free_cnt_memory(cnt2)

end program cnt_resonance_energy_transfer
