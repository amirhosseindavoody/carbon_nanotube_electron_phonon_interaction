!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Declaration of input parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module input_cnt_mod
	implicit none
	private
	public :: input_cnt_parameters, input_exciton

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
						case ('target_exciton_type')
							read(value, *) currcnt%target_exciton_type
						case ('length[nm]')
							read(value, *) currcnt%length
						case ('center_position[nm]')
							read(value, *) currcnt%center_position
						case default
							write(*,'(A,A)') "label = ", trim(label)
							write(*,'(A)') "ERROR in 'cnt' input arguments!!!"
							write(*,'(A)') "simulation STOPPED!!!"
							call exit()
						end select
					case ('flg')
						select case (trim(label))
						case ('flg_dielectric')
							read(value, *) flg_dielectric
						case default
							write(*,*) "ERROR in 'flg' input arguments!!!"
							write(*,*) "simulation STOPPED!!!"
							call exit()
						end select
					end select
				end if
			else if (ios .gt. 0) then
				write (*,*) "Error in reading input file!"
				call exit()
			end if
		end do
		close(101)

		! calculate kappa based on input parameters
		if ((currcnt%Ckappa .gt. 0.d0 ) .and. (currcnt%kappa_coeff .gt. 0.d0)) then
			currcnt%kappa = currcnt%Ckappa*currcnt%kappa_coeff
		end if

		! create the name of directory in which the cnt information is saved
		write(currcnt%directory,"( A, 'CNT(', I2.2, ',', I2.2, ')-nkg(', I4.4, ')-nr(', I4.4, ')-E_th(', F3.1, ')-Kcm_max(', F3.1, ')-i_sub(', I1.1, ')-Ckappa(', F3.1, ')/' )") trim(currcnt%directory), currcnt%n_ch, currcnt%m_ch, currcnt%nkg, currcnt%nr, currcnt%E_th/eV, currcnt%Kcm_max*1.d-9, currcnt%i_sub, currcnt%Ckappa

		! check if the directory exists
		folder_exists = .true.
		inquire(file=trim(currcnt%directory)//'/.', exist=folder_exists)
		if (.not. folder_exists) then
			write(log_input,'(A)') new_line('A')//"input folder for cnt exciton dispersion not found:"//new_line('A')//trim(currcnt%directory)//new_line('A')
			call write_log(log_input)
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
		type(cnt), intent(inout) :: currcnt
		character(len=*), intent(in) :: exciton_energy_filename
		character(len=*), intent(in) :: exciton_wavefunction_filename
		integer :: iX, iKcm, ikr
		real*8 :: tmpr
		character(len=1000) :: tmp_txt

		currcnt%excitons(ex_type,alpha)%i_sub = currcnt%i_sub

		currcnt%excitons(ex_type,alpha)%iKcm_max = currcnt%iKcm_max
		currcnt%excitons(ex_type,alpha)%iKcm_min = currcnt%iKcm_min
		currcnt%excitons(ex_type,alpha)%ikr_high = currcnt%ikr_high
		currcnt%excitons(ex_type,alpha)%ikr_low = currcnt%ikr_low
		currcnt%excitons(ex_type,alpha)%iKcm_min_fine = currcnt%iKcm_min_fine
		currcnt%excitons(ex_type,alpha)%iKcm_max_fine = currcnt%iKcm_max_fine

		currcnt%excitons(ex_type,alpha)%nx = currcnt%excitons(ex_type,alpha)%ikr_high-currcnt%excitons(ex_type,alpha)%ikr_low+1

		currcnt%excitons(ex_type,alpha)%spin = alpha

		! set information mu_cm and size of mu_r array
		if (ex_type .eq. 1) then
			currcnt%excitons(ex_type,alpha)%n_mu_r = 2
			currcnt%excitons(ex_type,alpha)%mu_cm = 0
		elseif (ex_type .eq. 2) then
			currcnt%excitons(ex_type,alpha)%n_mu_r = 2
			currcnt%excitons(ex_type,alpha)%mu_cm = 0
		elseif(ex_type .eq. 3) then
			currcnt%excitons(ex_type,alpha)%n_mu_r = 1
			currcnt%excitons(ex_type,alpha)%mu_cm = +1 * currcnt%min_sub(currcnt%i_sub)
		elseif(ex_type .eq. 4) then
			currcnt%excitons(ex_type,alpha)%n_mu_r = 1
			currcnt%excitons(ex_type,alpha)%mu_cm = -1 * currcnt%min_sub(currcnt%i_sub)
		endif

		! set value of mu_r
		allocate(currcnt%excitons(ex_type,alpha)%mu_r(currcnt%excitons(ex_type,alpha)%n_mu_r))
		if (ex_type .eq. 1) then
			currcnt%excitons(ex_type,alpha)%mu_r(1) = +1*currcnt%min_sub(currcnt%i_sub)
			currcnt%excitons(ex_type,alpha)%mu_r(2) = -1*currcnt%min_sub(currcnt%i_sub)
		elseif(ex_type .eq. 2) then
			currcnt%excitons(ex_type,alpha)%mu_r(1) = +1*currcnt%min_sub(currcnt%i_sub)
			currcnt%excitons(ex_type,alpha)%mu_r(2) = -1*currcnt%min_sub(currcnt%i_sub)
		else
			currcnt%excitons(ex_type,alpha)%mu_r(1) = 0
		endif

		! read exciton energy and wavefunction information
		allocate(currcnt%excitons(ex_type,alpha)%ex(1:currcnt%excitons(ex_type,alpha)%ikr_high-currcnt%excitons(ex_type,alpha)%ikr_low+1,currcnt%excitons(ex_type,alpha)%iKcm_min_fine:currcnt%excitons(ex_type,alpha)%iKcm_max_fine))
		allocate(currcnt%excitons(ex_type,alpha)%psi(currcnt%excitons(ex_type,alpha)%ikr_low:currcnt%excitons(ex_type,alpha)%ikr_high, currcnt%excitons(ex_type,alpha)%nx, currcnt%excitons(ex_type,alpha)%iKcm_min_fine:currcnt%excitons(ex_type,alpha)%iKcm_max_fine, currcnt%excitons(ex_type,alpha)%n_mu_r))

		write(tmp_txt,'(A, A)') trim(currcnt%directory), trim(exciton_energy_filename)
		open(unit=100,file=trim(tmp_txt),status="old")
		write(tmp_txt,'(A, A)') trim(currcnt%directory), trim(exciton_wavefunction_filename)
		open(unit=101,file=trim(tmp_txt),status="old")

		do iKcm=currcnt%excitons(ex_type,alpha)%iKcm_min_fine,currcnt%excitons(ex_type,alpha)%iKcm_max_fine
			do iX=1,currcnt%excitons(ex_type,alpha)%nx
				read(100,'(E16.8)', advance='no') currcnt%excitons(ex_type,alpha)%ex(iX,iKcm)
				do ikr=currcnt%excitons(ex_type,alpha)%ikr_low,currcnt%excitons(ex_type,alpha)%ikr_high
					read(101,'(E16.8,E16.8)', advance='no') currcnt%excitons(ex_type,alpha)%psi(ikr,iX,iKcm,1)
				enddo
			enddo

			read(100,'(E16.8)')
			read(101,'(E16.8)')
		enddo
		close(100)
		close(101)

		!make sure the exciton wavefunctions are normalized
		do iX=1,currcnt%excitons(ex_type,alpha)%ikr_high-currcnt%excitons(ex_type,alpha)%ikr_low+1
			do iKcm=currcnt%excitons(ex_type,alpha)%iKcm_min_fine,currcnt%excitons(ex_type,alpha)%iKcm_max_fine
				tmpr = 0.d0
				do ikr=currcnt%excitons(ex_type,alpha)%ikr_low,currcnt%excitons(ex_type,alpha)%ikr_high
					tmpr = tmpr + (abs(currcnt%excitons(ex_type,alpha)%psi(ikr,iX,iKcm,1)))**2
				enddo
				currcnt%excitons(ex_type,alpha)%psi(:,iX,iKcm,1) = currcnt%excitons(ex_type,alpha)%psi(:,iX,iKcm,1) / dcmplx(sqrt(tmpr))
			enddo
		enddo

		!make the coefficients of electronic states for the cutting lines close to K' point in A-type excitons
		if (ex_type .eq. 1) then
			currcnt%excitons(ex_type,alpha)%psi(currcnt%excitons(ex_type,alpha)%ikr_low:currcnt%excitons(ex_type,alpha)%ikr_high,:,:,2) = dcmplx(+1.d0)*currcnt%excitons(ex_type,alpha)%psi(currcnt%excitons(ex_type,alpha)%ikr_high:currcnt%excitons(ex_type,alpha)%ikr_low:-1,:,:,1)
			currcnt%excitons(ex_type,alpha)%psi = currcnt%excitons(ex_type,alpha)%psi/dcmplx(sqrt(2.d0))
		elseif(ex_type .eq. 2) then
			currcnt%excitons(ex_type,alpha)%psi(currcnt%excitons(ex_type,alpha)%ikr_low:currcnt%excitons(ex_type,alpha)%ikr_high,:,:,2) = dcmplx(-1.d0)*currcnt%excitons(ex_type,alpha)%psi(currcnt%excitons(ex_type,alpha)%ikr_high:currcnt%excitons(ex_type,alpha)%ikr_low:-1,:,:,1)
			currcnt%excitons(ex_type,alpha)%psi = currcnt%excitons(ex_type,alpha)%psi/dcmplx(sqrt(2.d0))
		endif

		!set exciton_name
		select case (ex_type)
		case(1)
			write(currcnt%excitons(ex_type,alpha)%name, '(A)') "A1"
		case(2)
			write(currcnt%excitons(ex_type,alpha)%name, '(A)') "A2"
		case(3)
			write(currcnt%excitons(ex_type,alpha)%name, '(A)') "Ep"
		case(4)
			write(currcnt%excitons(ex_type,alpha)%name, '(A)') "Em"
		case default
			call exit()
		end select

		select case (alpha)
		case(0)
			write(currcnt%excitons(ex_type,alpha)%name, '(A, A)') trim(currcnt%excitons(ex_type,alpha)%name), "_singlet"
		case(1)
			write(currcnt%excitons(ex_type,alpha)%name, '(A, A)') trim(currcnt%excitons(ex_type,alpha)%name), "_triplet"
		case default
			call exit()
		end select

	end subroutine input_exciton

end module input_cnt_mod
