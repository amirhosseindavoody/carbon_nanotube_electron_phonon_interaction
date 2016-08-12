!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Declaration of input parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module input_cnt_mod
	implicit none
	private
	public :: input_cnt_parameters, input_exciton, input_selected_exciton

contains
	!**************************************************************************************************************************
	! read the CNT information from the log file of exciton calculation which is done previously
	!**************************************************************************************************************************

	subroutine input_cnt_parameters(currcnt,filename)
		use cnt_class, only: cnt
		use constants_mod, only: eV
		use write_log_mod, only: write_log, log_input

		type(cnt) :: currcnt
		character(len=*),intent(in) :: filename

		character(len=1000) :: buffer, label, command, value
		integer :: ios
		integer :: istat
		integer :: pos_comma, pos_equal
		logical :: flg_dielectric, folder_exists

		ios=0
		istat=0
		pos_comma=0
		pos_equal=0

		open(unit=101,file=filename,status="old", action="read", iostat=istat)
		if (istat .ne. 0) then
			write(*,'(A,I2.2)') "istat = ", istat
			write(*,'(A,A)') "Unable to read CNT input file: ", filename
			call exit()
		end if

		! set the simulation variables to default values
		currcnt%n_ch=20
		currcnt%m_ch=0
		currcnt%nkg=0501
		currcnt%nr=200
		currcnt%E_th=1.5d0
		currcnt%Kcm_max=1.5d9
		currcnt%i_sub=1
		currcnt%kappa=2.d0
		currcnt%Ckappa=0.d0
		currcnt%kappa_coeff=0.d0
		flg_dielectric=.true.  !when .true. dielectric function is calculated, when .false. dielectric function is read from file.

		do while (ios == 0)
			read (101,'(A)',iostat=ios) buffer
			if (ios == 0) then
				if (buffer(1:1) .ne. '#') then
					pos_comma = scan(buffer,',')
					pos_equal = scan(buffer,'=')
					command = adjustl(buffer(1:pos_comma-1))
					label = adjustl(buffer(pos_comma+1:pos_equal-1))
					value = adjustl(buffer(pos_equal+1:))

					! set the target cnt
					select case (trim(command))
					case('directory')
						select case (trim(label))
						case ('output')
							currcnt%directory = trim(value)
						case default
							write(*,*) "ERROR in 'directory' input arguments!!!"
							write(*,*) "simulation STOPPED!!!"
							call exit()
						end select
					case ('cnt')
						select case (trim(label))
						case ('n_ch')
							read(value, *) currcnt%n_ch
						case ('m_ch')
							read(value, *) currcnt%m_ch
						case ('nkg')
							read(value, *) currcnt%nkg
						case ('dk/dkx')
							read(value, *) currcnt%dk_dkx_ratio
						case ('nr')
							read(value, *) currcnt%nr
						case ('E_th[eV]')
							read(value, *) currcnt%E_th
							currcnt%E_th = currcnt%E_th * eV
						case ('Kcm_max[1/nm]')
							read(value, *) currcnt%Kcm_max
							currcnt%Kcm_max = currcnt%Kcm_max * 1.d9
						case ('i_sub')
							read(value, *) currcnt%i_sub
						case ('Ckappa')
							read(value, *) currcnt%Ckappa
						case ('kappa_coeff')
							read(value, *) currcnt%kappa_coeff
						case ('selected_exciton')
							read(value, *) currcnt%selected_exciton_name
						case ('length[nm]')
							read(value, *) currcnt%length
							currcnt%length = currcnt%length*1.d-9
						case ('center_position[nm]')
							read(value, *) currcnt%center_position
						case default
							write(*,'(A, A, A, A)') "ERROR in interpreting this line:", new_line('A'), trim(buffer), new_line('A')
						end select
					case ('flg')
						select case (trim(label))
						case ('flg_dielectric')
							read(value, *) flg_dielectric
						case default
							write(*,'(A, A, A, A)') "ERROR in interpreting this line:", new_line('A'), trim(buffer), new_line('A')
						end select
					end select
				end if
			else if (ios .gt. 0) then
				write (*,'(A)') "Error in reading input file!"
				call exit()
			end if
		end do
		close(101)

		! calculate kappa based on input parameters
		if ((currcnt%Ckappa .gt. 0.d0 ) .and. (currcnt%kappa_coeff .gt. 0.d0)) then
			currcnt%kappa = currcnt%Ckappa*currcnt%kappa_coeff
		end if

		! create the name of directory in which the cnt information is saved
		! write(currcnt%directory,"( A, 'CNT(', I2.2, ',', I2.2, ')-nkg(', I4.4, ')-nr(', I4.4, ')-E_th(', F3.1, ')-Kcm_max(', F3.1, ')-i_sub(', I1.1, ')-Ckappa(', F3.1, ')/' )") trim(currcnt%directory), currcnt%n_ch, currcnt%m_ch, currcnt%nkg, currcnt%nr, currcnt%E_th/eV, currcnt%Kcm_max*1.d-9, currcnt%i_sub, currcnt%Ckappa
		write(currcnt%directory,'(A, A, I0, A, I0, A, I0, A, I0, A, I0, A, F0.1, A, F0.1, A, I0, A, F0.1, A)') trim(currcnt%directory), "exciton_", currcnt%n_ch, "_", currcnt%m_ch, "_nkg_", currcnt%nkg, "_dk_ratio_", currcnt%dk_dkx_ratio, "_nr_", currcnt%nr, "_Eth_", currcnt%E_th/eV, "_Kcm_max_", currcnt%Kcm_max*1.d-9, "_sub_", currcnt%i_sub, "_Ckappa_", currcnt%Ckappa, "/"

		! check if the directory exists
		folder_exists = .true.
		inquire(file=trim(currcnt%directory)//'/.', exist=folder_exists)
		if (.not. folder_exists) then
			write(log_input,'(A)') new_line('A')//"input folder for cnt exciton dispersion not found:"//new_line('A')//trim(currcnt%directory)//new_line('A');			call write_log(log_input)
			call exit()
		end if

	end subroutine input_cnt_parameters

	!***************************************************************************
	! load the exciton wavefunction and energy from the ExcitonEnergy calculation
	!***************************************************************************

	subroutine input_exciton(ex_type, alpha, currcnt, exciton_energy_filename, exciton_wavefunction_filename)
		use cnt_class, only: cnt, exciton
		use write_log_mod, only: write_log

		integer, intent(in) :: ex_type
		integer, intent(in) :: alpha
		type(cnt), target, intent(inout) :: currcnt
		character(len=*), intent(in) :: exciton_energy_filename
		character(len=*), intent(in) :: exciton_wavefunction_filename
		integer :: iX, iKcm, ikr
		real*8 :: tmpr
		character(len=1000) :: tmp_txt
		type(exciton), pointer :: my_exciton

		my_exciton => currcnt%excitons(ex_type,alpha)

		my_exciton%i_sub = currcnt%i_sub

		my_exciton%iKcm_max = currcnt%iKcm_max_fine
		my_exciton%iKcm_min = currcnt%iKcm_min_fine
		my_exciton%ikr_high = currcnt%ikr_high
		my_exciton%ikr_low = currcnt%ikr_low

		my_exciton%nx = my_exciton%ikr_high-my_exciton%ikr_low+1

		my_exciton%dKcm = currcnt%dkx
		my_exciton%dkr = currcnt%dk
		my_exciton%dkr_dKcm_ratio = currcnt%dk_dkx_ratio

		my_exciton%spin = alpha

		! set information mu_cm and size of mu_r array
		if (ex_type .eq. 1) then
			my_exciton%n_mu_r = 2
			my_exciton%mu_cm = 0
		elseif (ex_type .eq. 2) then
			my_exciton%n_mu_r = 2
			my_exciton%mu_cm = 0
		elseif(ex_type .eq. 3) then
			my_exciton%n_mu_r = 1
			my_exciton%mu_cm = +1 * currcnt%min_sub(currcnt%i_sub)
		elseif(ex_type .eq. 4) then
			my_exciton%n_mu_r = 1
			my_exciton%mu_cm = -1 * currcnt%min_sub(currcnt%i_sub)
		endif

		! set value of mu_r
		if (.not. allocated(my_exciton%mu_r)) allocate(my_exciton%mu_r(my_exciton%n_mu_r))
		if (ex_type .eq. 1) then
			my_exciton%mu_r(1) = +1*currcnt%min_sub(currcnt%i_sub)
			my_exciton%mu_r(2) = -1*currcnt%min_sub(currcnt%i_sub)
		elseif(ex_type .eq. 2) then
			my_exciton%mu_r(1) = +1*currcnt%min_sub(currcnt%i_sub)
			my_exciton%mu_r(2) = -1*currcnt%min_sub(currcnt%i_sub)
		else
			my_exciton%mu_r(1) = 0
		endif

! 		! read exciton energy and wavefunction information
! 		if ((.not. allocated(my_exciton%ex)) .and. (.not. allocated(my_exciton%psi))) then
! 			allocate(my_exciton%ex(1:my_exciton%ikr_high-my_exciton%ikr_low+1,my_exciton%iKcm_min:my_exciton%iKcm_max))
! 			allocate(my_exciton%psi(my_exciton%ikr_low:my_exciton%ikr_high, my_exciton%nx, my_exciton%iKcm_min:my_exciton%iKcm_max, my_exciton%n_mu_r))

! 			write(tmp_txt,'(A, A)') trim(currcnt%directory), trim(exciton_energy_filename)
! 			open(unit=100,file=trim(tmp_txt),status="old")
! 			write(tmp_txt,'(A, A)') trim(currcnt%directory), trim(exciton_wavefunction_filename)
! 			open(unit=101,file=trim(tmp_txt),status="old")

! 			do iKcm=my_exciton%iKcm_min, my_exciton%iKcm_max
! 				do iX=1,my_exciton%nx
! 					read(100,'(E16.8)', advance='no') my_exciton%ex(iX,iKcm)
! 					do ikr=my_exciton%ikr_low,my_exciton%ikr_high
! 						read(101,'(E16.8,E16.8)', advance='no') my_exciton%psi(ikr,iX,iKcm,1)
! 					enddo
! 				enddo

! 				read(100,'(E16.8)')
! 				read(101,'(E16.8)')
! 			enddo
! 			close(100)
! 			close(101)

! 			!make sure the exciton wavefunctions are normalized
! 			do iX=1,my_exciton%ikr_high-my_exciton%ikr_low+1
! 				do iKcm=my_exciton%iKcm_min, my_exciton%iKcm_max
! 					tmpr = 0.d0
! 					do ikr=my_exciton%ikr_low,my_exciton%ikr_high
! 						tmpr = tmpr + (abs(my_exciton%psi(ikr,iX,iKcm,1)))**2
! 					enddo
! 					my_exciton%psi(:,iX,iKcm,1) = my_exciton%psi(:,iX,iKcm,1) / dcmplx(sqrt(tmpr))
! 				enddo
! 			enddo

! 			!make the coefficients of electronic states for the cutting lines close to K' point in A-type excitons
! 			if (ex_type .eq. 1) then
! 				my_exciton%psi(my_exciton%ikr_low:my_exciton%ikr_high,:,:,2) = dcmplx(-1.d0)*my_exciton%psi(my_exciton%ikr_high:my_exciton%ikr_low:-1,:,:,1)
! 				my_exciton%psi = my_exciton%psi/dcmplx(sqrt(2.d0))
! 			elseif(ex_type .eq. 2) then
! 				my_exciton%psi(my_exciton%ikr_low:my_exciton%ikr_high,:,:,2) = dcmplx(+1.d0)*my_exciton%psi(my_exciton%ikr_high:my_exciton%ikr_low:-1,:,:,1)
! 				my_exciton%psi = my_exciton%psi/dcmplx(sqrt(2.d0))
! 			endif
! 		endif

		! read exciton energy and wavefunction information
		if ((.not. allocated(my_exciton%ex)) .and. (.not. allocated(my_exciton%psi))) then
			allocate(my_exciton%ex(1:my_exciton%ikr_high-my_exciton%ikr_low+1,my_exciton%iKcm_min:my_exciton%iKcm_max))
			allocate(my_exciton%psi(my_exciton%ikr_low:my_exciton%ikr_high, my_exciton%nx, my_exciton%iKcm_min:my_exciton%iKcm_max, my_exciton%n_mu_r))

			write(tmp_txt,'(A, A)') trim(currcnt%directory), trim(exciton_energy_filename)
			open(unit=100,file=trim(tmp_txt),status="old", form='unformatted')
			write(tmp_txt,'(A, A)') trim(currcnt%directory), trim(exciton_wavefunction_filename)
			open(unit=101,file=trim(tmp_txt),status="old", form='unformatted')

			read(100) my_exciton%ex(:,:)
			read(101) my_exciton%psi(:,:,:,1)

			close(100)
			close(101)

			!make sure the exciton wavefunctions are normalized
			do iX=1,my_exciton%ikr_high-my_exciton%ikr_low+1
				do iKcm=my_exciton%iKcm_min, my_exciton%iKcm_max
					tmpr = 0.d0
					do ikr=my_exciton%ikr_low,my_exciton%ikr_high
						tmpr = tmpr + (abs(my_exciton%psi(ikr,iX,iKcm,1)))**2
					enddo
					my_exciton%psi(:,iX,iKcm,1) = my_exciton%psi(:,iX,iKcm,1) / dcmplx(sqrt(tmpr))
				enddo
			enddo

			!make the coefficients of electronic states for the cutting lines close to K' point in A-type excitons
			if (ex_type .eq. 1) then
				my_exciton%psi(my_exciton%ikr_low:my_exciton%ikr_high,:,:,2) = dcmplx(-1.d0)*my_exciton%psi(my_exciton%ikr_high:my_exciton%ikr_low:-1,:,:,1)
				my_exciton%psi = my_exciton%psi/dcmplx(sqrt(2.d0))
			elseif(ex_type .eq. 2) then
				my_exciton%psi(my_exciton%ikr_low:my_exciton%ikr_high,:,:,2) = dcmplx(+1.d0)*my_exciton%psi(my_exciton%ikr_high:my_exciton%ikr_low:-1,:,:,1)
				my_exciton%psi = my_exciton%psi/dcmplx(sqrt(2.d0))
			endif
		endif

		!set exciton_name
		select case (ex_type)
		case(1)
			write(my_exciton%name, '(A)') "A1"
		case(2)
			write(my_exciton%name, '(A)') "A2"
		case(3)
			write(my_exciton%name, '(A)') "Ep"
		case(4)
			write(my_exciton%name, '(A)') "Em"
		case default
			call exit()
		end select

		select case (alpha)
		case(0)
			write(my_exciton%name, '(A, A)') trim(my_exciton%name), "_singlet"
		case(1)
			write(my_exciton%name, '(A, A)') trim(my_exciton%name), "_triplet"
		case default
			call exit()
		end select

	end subroutine input_exciton


	!***************************************************************************
	! input the exciton information for the selected exciton type.
	!***************************************************************************

	subroutine input_selected_exciton(my_cnt)
		use cnt_class, only: cnt, exciton
		use write_log_mod, only: write_log, log_input

		type(cnt), target, intent(inout) :: my_cnt

		select case (trim(my_cnt%selected_exciton_name))
		case("A1_singlet")
			call input_exciton(ex_type=1, alpha=0, currcnt=my_cnt, exciton_energy_filename='Ex_A1.dat', exciton_wavefunction_filename='Psi_A1.dat')
			my_cnt%selected_exciton => my_cnt%excitons(1,0)
		case("A2_singlet")
			call input_exciton(ex_type=2, alpha=0, currcnt=my_cnt, exciton_energy_filename='Ex0_A2.dat', exciton_wavefunction_filename='Psi0_A2.dat')
			my_cnt%selected_exciton => my_cnt%excitons(2,0)
		case("Ep_singlet")
			call input_exciton(ex_type=3, alpha=0, currcnt=my_cnt, exciton_energy_filename='Ex0_Ep.dat', exciton_wavefunction_filename='Psi0_Ep.dat')
			my_cnt%selected_exciton => my_cnt%excitons(3,0)
		case("Em_singlet")
			call input_exciton(ex_type=4, alpha=0, currcnt=my_cnt, exciton_energy_filename='Ex0_Em.dat', exciton_wavefunction_filename='Psi0_Em.dat')
			my_cnt%selected_exciton => my_cnt%excitons(4,0)
		case("A1_triplet")
			call input_exciton(ex_type=1, alpha=1, currcnt=my_cnt, exciton_energy_filename='Ex_A1.dat', exciton_wavefunction_filename='Psi_A1.dat')
			my_cnt%selected_exciton => my_cnt%excitons(1,1)
		case("A2_triplet")
			call input_exciton(ex_type=2, alpha=1, currcnt=my_cnt, exciton_energy_filename='Ex1_A2.dat', exciton_wavefunction_filename='Psi1_A2.dat')
			my_cnt%selected_exciton => my_cnt%excitons(2,1)
		case("Ep_triplet")
			call input_exciton(ex_type=3, alpha=1, currcnt=my_cnt, exciton_energy_filename='Ex1_Ep.dat', exciton_wavefunction_filename='Psi1_Ep.dat')
			my_cnt%selected_exciton => my_cnt%excitons(3,1)
		case("Em_triplet")
			call input_exciton(ex_type=4, alpha=1, currcnt=my_cnt, exciton_energy_filename='Ex1_Em.dat', exciton_wavefunction_filename='Psi1_Em.dat')
			my_cnt%selected_exciton => my_cnt%excitons(4,1)
		case default
			write(log_input, '(A, A)') "Incorrect selected_exciton_name:", trim(my_cnt%selected_exciton_name)
			call write_log(log_input)
			call exit()
		end select

	end subroutine input_selected_exciton

end module input_cnt_mod
