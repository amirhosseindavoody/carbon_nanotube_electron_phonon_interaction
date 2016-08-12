module sim_properties_mod
	implicit none

	character(len=1000) :: outdir_tmp
	character(len=1000) :: outdir_final

	real :: starttime,endtime !duration of simulation
	real*8 :: temperature = 300.d0 !temperature of the system in Kelvin units, this is the default value
	real*8 :: ppLen = 10.d-9 !length per perpendicular tubes, this is the default value

	real*8 :: c2c_min, c2c_max
	integer :: n_c2c

	real*8 :: theta_min, theta_max
	integer :: n_theta

	integer :: energy_mesh_size

contains

	!***************************************************************************
	! -	this subroutine reads the input file that contains the properties of
	!	the simulation for calculating scattering rate.
	!***************************************************************************
	subroutine input_sim_properties(filename)
		use constants_mod, only: pi
		use write_log_mod, only : write_log, log_input

		character(len=*) :: filename
		integer :: istat=0
		integer :: ios=0
		integer :: pos_comma=0, pos_equal=0
		character(len=1000) :: buffer, command, label, value
		integer :: i_tmp=0
		logical :: folder_exists=.true.
		integer, dimension(3) :: date, time


		! set some default values for simulation properties. These values are most likely going to be over written later.
		temperature = 300.d0
		ppLen = 10.d-9
		c2c_min = 1.2d-9
		c2c_max = 1.2d-9
		n_c2c = 1
		theta_min = 0.d0
		theta_max = pi/2.d0
		n_theta = 2


		! read and set the simulation properties variables
		open(unit=100,file=filename,status="old", action="read", iostat=istat)
		if (istat .ne. 0) then
			write(*,*) ""
			write(*,*) "Unable to read simulation properties input file:", filename
			call exit()
		end if

		do while (ios == 0)
			read (100,'(A)',iostat=ios) buffer
			if (ios == 0) then
				if (buffer(1:1) .ne. '#') then
					pos_comma = scan(buffer,',')
					pos_equal = scan(buffer,'=')
					command = adjustl(buffer(1:pos_comma-1))
					label = adjustl(buffer(pos_comma+1:pos_equal-1))
					value = adjustl(buffer(pos_equal+1:))

					select case (trim(command))
					case('directory')
						select case (trim(label))
						case ('output')
							outdir_final = trim(value)
						end select
					case ('geometry')
						select case (trim(label))
						case('c2c_min[nm]')
							read(value, *) c2c_min
							c2c_min = c2c_min * 1.d-9
						case('c2c_max[nm]')
							read(value, *) c2c_max
							c2c_max = c2c_max * 1.d-9
						case('n_c2c')
							read(value, *) n_c2c
						case('theta_min[degree]')
							read(value, *) theta_min
							theta_min = theta_min * pi/180
						case('theta_max[degree]')
							read(value, *) theta_max
							theta_max = theta_max * pi/180
						case ('n_theta')
							read(value, *) n_theta
						case ('len_per_perpendicular_tube[nm]')
							read(value, *) ppLen
							ppLen = ppLen * 1d-9
						end select
					case('simulation')
						select case (trim(label))
						case ('temperature[Kelvin]')
							read(value,*) temperature
						case ('energy_mesh_size')
							read(value,*) energy_mesh_size
						end select
					end select
				end if
			else if (ios .gt. 0) then
				write (*,'(A)') "Error in reading input file!"
				call exit()
			end if
		end do
		close(100)

		! create the temporary name of the output directory
		folder_exists = .true.
		i_tmp = 0
		do while (folder_exists)
			i_tmp = i_tmp+1
			write(outdir_tmp,"(A, A, I0)") trim(outdir_final),"tmp_",i_tmp
			inquire(file=trim(outdir_tmp)//'/.', exist=folder_exists)

			if (.not. folder_exists) then
				!create the output directory
				write(command,'(A ,A, A)') "mkdir '", trim(outdir_tmp), "'"
				call system(trim(command), status=istat)
				if (istat .ne. 0) then
					folder_exists = .true.
				endif
			endif

		end do

		! !create the output directory
		! write(command,'(A ,A, A)') "mkdir '", trim(outdir_tmp), "'"
		! call system(trim(command))

		! copy the input files to the output directory
		call get_command_argument(1,buffer)
		write(command,'(A, A, A, A)') "cp ", trim(buffer), " ", trim(outdir_tmp)
		call system(trim(command))
		call get_command_argument(2,buffer)
		write(command,'(A, A, A, A)') "cp ", trim(buffer), " ", trim(outdir_tmp)
		call system(trim(command))
		call get_command_argument(3,buffer)
		write(command,'(A, A, A, A)') "cp ", trim(buffer), " ", trim(outdir_tmp)
		call system(trim(command))

		!change the working directory to the output directory
		istat=chdir(trim(outdir_tmp))
		if (istat .ne. 0) then
			write(*,'(A)') "Directory did not changed!!!"
			write(*,'(A)') "Simulation stopped!!!"
			call exit()
		end if

		! get time and date of start of simulation and write it to the log file
		call idate(date)
		call itime(time)
		write(log_input,'(A, I0, A, I0, A, I0, A, I0, A, I0, A, I0)') "Simulation started at:"//new_line('A')//"date = ",date(1),"/",date(2),"/",date(3),new_line('A')//"time = ", time(1),":",time(2),":",time(3)
		call write_log(log_input)


	end subroutine input_sim_properties

	!***************************************************************************
	! -	this subroutine renames the output directory from a temporary name to a
	!	final name that indicates the simulation has run successfully.
	!***************************************************************************
	subroutine finalize_output_directory_name(cnt_1, cnt_2)
		use cnt_class, only: cnt
		use constants_mod, only: pi

		type(cnt), intent(in) :: cnt_1, cnt_2
		logical :: folder_exists
		character(len=1000) :: command
		integer :: istat

		! write(outdir_final,"(A, A)") trim(outdir_final), "final_result"
		write(outdir_final,'(A, I0, A, I0, A, A, A, I0, A, I0, A, I0, A, F0.1, A, I0, A, I0, A, I0, A, A, A, I0, A, I0, A, I0, A, F0.1, A, I0, A, F0.1, A, F0.1, A, I0, A, I0 )') trim(outdir_final)//"transfer_", cnt_1%n_ch, "_", cnt_1%m_ch, "_", trim(cnt_1%selected_exciton%name), "_iSub_", cnt_1%i_sub, "_length_", nint(cnt_1%length*1.d9), "nm_center_", nint(cnt_1%center_position*1.d9), "nm_Ckappa_", cnt_1%Ckappa,"_dk_ratio_", cnt_1%dk_dkx_ratio, "_to_", cnt_2%n_ch, "_", cnt_2%m_ch, "_", trim(cnt_2%selected_exciton%name), "_iSub_", cnt_2%i_sub, "_length_", nint(cnt_2%length*1.d9), "nm_center_", nint(cnt_2%center_position*1.d9), "nm_Ckappa_", cnt_2%Ckappa, "_dk_ratio_", cnt_2%dk_dkx_ratio, "_C2C_", c2c_min*1.d9, "nm_", c2c_max*1.d9, "nm_theta_", nint(theta_min*180/pi), "_", nint(theta_max*180/pi)

		! remove the final output directory if it already exists
		folder_exists = .true.
		inquire(file=trim(outdir_final)//'/.', exist=folder_exists)

		if (folder_exists) then
			write(command, "(A, A)") "rm -r ", trim(outdir_final)
			call system(trim(command))
		end if

		!rename the temporary output directory to the final output directory
		write(command, '(A, A, A, A)') "mv ", trim(outdir_tmp), " ", trim(outdir_final)
		call system(trim(command))

		!change the working directory to the final output directory
		istat=chdir(trim(outdir_final))
		if (istat .ne. 0) then
			write(*,'(A)') "Directory did not changed!!!"
			write(*,'(A)') "Simulation stopped!!!"
			call exit()
		end if

	end subroutine finalize_output_directory_name

end module sim_properties_mod
