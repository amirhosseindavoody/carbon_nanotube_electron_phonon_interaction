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

	integer :: partition_function_type


contains
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

		open(unit=100,file=filename,status="old", action="read", iostat=istat)
		if (istat .ne. 0) then
			write(*,*) ""
			write(*,*) "Unable to read simulation properties input file:", filename
			call exit()
		end if

		! read and set the simulation properties variables
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
						end select
					case('simulation')
						select case (trim(label))
						case ('partition_function')
							partition_function_type = -1
							if (trim(value) .eq. 'total') partition_function_type = 0
							if (trim(value) .eq. 'target') partition_function_type = 1
						case ('temperature[Kelvin]')
							read(value,*) temperature
						end select
					end select
				end if
			else if (ios .gt. 0) then
				write (*,*) "Error in reading input file!"
				call exit()
			end if
		end do
		close(100)

		! create the temporary name of the output directory
		folder_exists = .true.
		i_tmp = 0
		do while (folder_exists)
			i_tmp = i_tmp+1
			write(outdir_tmp,"(A,A,I3.3)") trim(outdir_final),"tmp_",i_tmp
			inquire(file=trim(outdir_tmp)//'/.', exist=folder_exists)
		end do

		!create the output directory
		write(command,'("mkdir ''",A,"''")') trim(outdir_tmp)
		call system(trim(command))

		!change the working directory to the output directory
		istat=chdir(trim(outdir_tmp))
		if (istat .ne. 0) then
			write(*,'(A)') "Directory did not changed!!!"
			write(*,'(A)') "Simulation stopped!!!"
			call exit()
		end if

		! copy the input files to the output directory
		call get_command_argument(1,buffer)
		write(command,'("cp ",A," ",A)')trim(buffer), trim(outdir_tmp)
		call system(trim(command))
		call get_command_argument(2,buffer)
		write(command,'("cp ",A," ",A)')trim(buffer), trim(outdir_tmp)
		call system(trim(command))
		call get_command_argument(3,buffer)
		write(command,'("cp ",A," ",A)')trim(buffer), trim(outdir_tmp)
		call system(trim(command))

		! get time and date of start of simulation and write it to the log file
		call idate(date)
		call itime(time)
		write(log_input,'(A,I2.2,A,I2.2,A,I4.4,A,I2.2,A,I2.2,A,I2.2)') "Simulation started at:"//new_line('A')//"date = ",date(1),"/",date(2),"/",date(3),new_line('A')//"time = ", time(1),":",time(2),":",time(3)
		call write_log(log_input)


	end subroutine input_sim_properties

end module sim_properties_mod
