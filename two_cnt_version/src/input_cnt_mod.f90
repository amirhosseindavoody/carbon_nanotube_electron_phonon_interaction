!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Declaration of input parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module input_cnt_mod
	implicit none
	private
	public :: input_cnt_parameters, input_a_exciton

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

	!**************************************************************************************************************************
	! load the A-type exciton wavefunction and energies from the ExcitonEnergy calculation
	!**************************************************************************************************************************

	subroutine input_a_exciton(currcnt)
		use cnt_class, only: cnt
		use write_log_mod, only: write_log

		type(cnt), intent(inout) :: currcnt
		integer :: iX, iKcm, ikr
		real*8 :: tmpr1, tmpr2, tmpr3

		! ! read the information of A1 exciton
		! allocate(currcnt%Ex_A1(1:currcnt%ikr_high-currcnt%ikr_low+1,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
		! allocate(currcnt%Psi_A1(currcnt%ikr_low:currcnt%ikr_high,1:currcnt%ikr_high-currcnt%ikr_low+1,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
		!
		! open(unit=100,file=trim(currcnt%directory)//'Ex_A1.dat',status="old")
		! open(unit=101,file=trim(currcnt%directory)//'Psi_A1.dat',status="old")
		! do iKcm=currcnt%iKcm_min_fine,currcnt%iKcm_max_fine
		! 	do iX=1,currcnt%ikr_high-currcnt%ikr_low+1
		! 		read(100,'(E16.8)', advance='no') currcnt%Ex_A1(iX,iKcm)
		! 		do ikr=currcnt%ikr_low,currcnt%ikr_high
		! 			read(101,'(E16.8,E16.8)', advance='no') currcnt%Psi_A1(ikr,iX,iKcm)
		! 		enddo
		! 	enddo
		!
		! 	read(100,'(E16.8)')
		! 	read(101,'(E16.8)')
		! enddo
		! close(100)
		! close(101)

		! read the information of singlet A2 exciton
		allocate(currcnt%Ex0_A2(1:currcnt%ikr_high-currcnt%ikr_low+1,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
		allocate(currcnt%Psi0_A2(currcnt%ikr_low:currcnt%ikr_high,1:currcnt%ikr_high-currcnt%ikr_low+1,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))

		open(unit=100,file=trim(currcnt%directory)//'Ex0_A2.dat',status="old")
		open(unit=101,file=trim(currcnt%directory)//'Psi0_A2.dat',status="old")

		do iKcm=currcnt%iKcm_min_fine,currcnt%iKcm_max_fine
			do iX=1,currcnt%ikr_high-currcnt%ikr_low+1
				read(100,'(E16.8)', advance='no') currcnt%Ex0_A2(iX,iKcm)
				do ikr=currcnt%ikr_low,currcnt%ikr_high
					read(101,'(E16.8,E16.8)', advance='no') currcnt%Psi0_A2(ikr,iX,iKcm)
				enddo
			enddo

			read(100,'(E16.8)')
			read(101,'(E16.8)')
		enddo
		close(100)
		close(101)

		! ! read the information of triplet A2 exciton
		! allocate(currcnt%Ex1_A2(1:currcnt%ikr_high-currcnt%ikr_low+1,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
		! allocate(currcnt%Psi1_A2(currcnt%ikr_low:currcnt%ikr_high,1:currcnt%ikr_high-currcnt%ikr_low+1,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
		!
		! open(unit=100,file=trim(currcnt%directory)//'Ex1_A2.dat',status="old")
		! open(unit=101,file=trim(currcnt%directory)//'Psi1_A2.dat',status="old")
		! do iKcm=currcnt%iKcm_min_fine,currcnt%iKcm_max_fine
		! 	do iX=1,currcnt%ikr_high-currcnt%ikr_low+1
		! 		read(100,'(E16.8)', advance='no') currcnt%Ex1_A2(iX,iKcm)
		! 		do ikr=currcnt%ikr_low,currcnt%ikr_high
		! 			read(101,'(E16.8,E16.8)', advance='no') currcnt%Psi1_A2(ikr,iX,iKcm)
		! 		enddo
		! 	enddo
		!
		! 	read(100,'(E16.8)')
		! 	read(101,'(E16.8)')
		! enddo
		! close(100)
		! close(101)

		!make sure the exciton wavefunctions are normalized
		do iX=1,currcnt%ikr_high-currcnt%ikr_low+1
			do iKcm=currcnt%iKcm_min_fine,currcnt%iKcm_max_fine
				tmpr1 = 0.d0
				tmpr2 = 0.d0
				tmpr3 = 0.d0
				do ikr=currcnt%ikr_low,currcnt%ikr_high
					! tmpr1 = tmpr1 + abs(currcnt%Psi_A1(ikr,iX,iKcm))
					tmpr2 = tmpr2 + abs(currcnt%Psi0_A2(ikr,iX,iKcm))
					! tmpr3 = tmpr3 + abs(currcnt%Psi1_A2(ikr,iX,iKcm))
				enddo
				! currcnt%Psi_A1(:,iX,iKcm) = currcnt%Psi_A1(:,iX,iKcm) / dcmplx(sqrt(tmpr1))
				currcnt%Psi0_A2(:,iX,iKcm) = currcnt%Psi0_A2(:,iX,iKcm) / dcmplx(sqrt(tmpr2))
				! currcnt%Psi1_A2(:,iX,iKcm) = currcnt%Psi1_A2(:,iX,iKcm) / dcmplx(sqrt(tmpr3))
			enddo
		enddo


		! ! set the information of the target exciton type
		! select case (trim(currcnt%target_exciton_type))
		! case ('Ex_A1')
		! 	call write_log("Target exciton: Ex_A1")
		! 	allocate(currcnt%Ex_t(1:currcnt%ikr_high-currcnt%ikr_low+1,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
		! 	allocate(currcnt%Psi_t(currcnt%ikr_low:currcnt%ikr_high,1:currcnt%ikr_high-currcnt%ikr_low+1,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
		! 	currcnt%Ex_t = currcnt%Ex_A1
		! 	currcnt%Psi_t = currcnt%Psi_A1
		! 	currcnt%ex_symmetry = -1.d0
		! case ('Ex0_A2')
		! 	call write_log("Target exciton: Ex0_A2")
		! 	allocate(currcnt%Ex_t(1:currcnt%ikr_high-currcnt%ikr_low+1,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
		! 	allocate(currcnt%Psi_t(currcnt%ikr_low:currcnt%ikr_high,1:currcnt%ikr_high-currcnt%ikr_low+1,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
		! 	currcnt%Ex_t = currcnt%Ex0_A2
		! 	currcnt%Psi_t = currcnt%Psi0_A2
		! 	currcnt%ex_symmetry = +1.d0
		! case ('Ex1_A2')
		! 	call write_log("Target exciton: Ex1_A2")
		! 	allocate(currcnt%Ex_t(1:currcnt%ikr_high-currcnt%ikr_low+1,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
		! 	allocate(currcnt%Psi_t(currcnt%ikr_low:currcnt%ikr_high,1:currcnt%ikr_high-currcnt%ikr_low+1,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
		! 	currcnt%Ex_t = currcnt%Ex1_A2
		! 	currcnt%Psi_t = currcnt%Psi1_A2
		! 	currcnt%ex_symmetry = +1.d0
		! end select
		!
		! deallocate(currcnt%Psi_A1)
		! deallocate(currcnt%Psi0_A2)
		! deallocate(currcnt%Psi1_A2)

	end subroutine input_a_exciton

	!**************************************************************************************************************************
	! load the E-type exciton wavefunction and energies from the ExcitonEnergy calculation
	!**************************************************************************************************************************

	! subroutine input_E_exciton(currcnt)
	! 	use cnt_class, only: cnt
	! 	use write_log_mod, only: write_log
	!
	! 	type(cnt), intent(inout) :: currcnt
	! 	integer :: iX, iKcm, ikr
	! 	real*8 :: tmpr1, tmpr2, tmpr3, tmpr4
	!
	! 	! read the information of singlet E+ exciton
	! 	allocate(currcnt%Ex0_Ep(1:currcnt%ikr_high-currcnt%ikr_low+1,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
	! 	allocate(currcnt%Psi0_Ep(currcnt%ikr_low:currcnt%ikr_high,1:currcnt%ikr_high-currcnt%ikr_low+1,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
	!
	! 	open(unit=100,file=trim(currcnt%directory)//'Ex0_Ep.dat',status="old")
	! 	open(unit=101,file=trim(currcnt%directory)//'Psi0_Ep.dat',status="old")
	! 	do iKcm=currcnt%iKcm_min_fine,currcnt%iKcm_max_fine
	! 		do iX=1,currcnt%ikr_high-currcnt%ikr_low+1
	! 			read(100,'(E16.8)', advance='no') currcnt%Ex0_Ep(iX,iKcm)
	! 			do ikr=currcnt%ikr_low,currcnt%ikr_high
	! 				read(101,'(E16.8,E16.8)', advance='no') currcnt%Psi0_Ep(ikr,iX,iKcm)
	! 			enddo
	! 		enddo
	!
	! 		read(100,'(E16.8)')
	! 		read(101,'(E16.8)')
	! 	enddo
	! 	close(100)
	! 	close(101)
	!
	! 	! read the information of singlet E- exciton
	! 	allocate(currcnt%Ex0_Em(1:currcnt%ikr_high-currcnt%ikr_low+1,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
	! 	allocate(currcnt%Psi0_Em(currcnt%ikr_low:currcnt%ikr_high,1:currcnt%ikr_high-currcnt%ikr_low+1,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
	!
	! 	open(unit=100,file=trim(currcnt%directory)//'Ex0_Em.dat',status="old")
	! 	open(unit=101,file=trim(currcnt%directory)//'Psi0_Em.dat',status="old")
	! 	do iKcm=currcnt%iKcm_min_fine,currcnt%iKcm_max_fine
	! 		do iX=1,currcnt%ikr_high-currcnt%ikr_low+1
	! 			read(100,'(E16.8)', advance='no') currcnt%Ex0_Em(iX,iKcm)
	! 			do ikr=currcnt%ikr_low,currcnt%ikr_high
	! 				read(101,'(E16.8,E16.8)', advance='no') currcnt%Psi0_Em(ikr,iX,iKcm)
	! 			enddo
	! 		enddo
	!
	! 		read(100,'(E16.8)')
	! 		read(101,'(E16.8)')
	! 	enddo
	! 	close(100)
	! 	close(101)
	!
	! 	! read the information of triplet E+ exciton
	! 	allocate(currcnt%Ex1_Ep(1:currcnt%ikr_high-currcnt%ikr_low+1,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
	! 	allocate(currcnt%Psi1_Ep(currcnt%ikr_low:currcnt%ikr_high,1:currcnt%ikr_high-currcnt%ikr_low+1,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
	!
	! 	open(unit=100,file=trim(currcnt%directory)//'Ex1_Ep.dat',status="old")
	! 	open(unit=101,file=trim(currcnt%directory)//'Psi1_Ep.dat',status="old")
	! 	do iKcm=currcnt%iKcm_min_fine,currcnt%iKcm_max_fine
	! 		do iX=1,currcnt%ikr_high-currcnt%ikr_low+1
	! 			read(100,'(E16.8)', advance='no') currcnt%Ex1_Ep(iX,iKcm)
	! 			do ikr=currcnt%ikr_low,currcnt%ikr_high
	! 				read(101,'(E16.8,E16.8)', advance='no') currcnt%Psi1_Ep(ikr,iX,iKcm)
	! 			enddo
	! 		enddo
	!
	! 		read(100,'(E16.8)')
	! 		read(101,'(E16.8)')
	! 	enddo
	! 	close(100)
	! 	close(101)
	!
	! 	! read the information of triplet E- exciton
	! 	allocate(currcnt%Ex1_Em(1:currcnt%ikr_high-currcnt%ikr_low+1,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
	! 	allocate(currcnt%Psi1_Em(currcnt%ikr_low:currcnt%ikr_high,1:currcnt%ikr_high-currcnt%ikr_low+1,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
	!
	! 	open(unit=100,file=trim(currcnt%directory)//'Ex1_Em.dat',status="old")
	! 	open(unit=101,file=trim(currcnt%directory)//'Psi1_Em.dat',status="old")
	! 	do iKcm=currcnt%iKcm_min_fine,currcnt%iKcm_max_fine
	! 		do iX=1,currcnt%ikr_high-currcnt%ikr_low+1
	! 			read(100,'(E16.8)', advance='no') currcnt%Ex1_Em(iX,iKcm)
	! 			do ikr=currcnt%ikr_low,currcnt%ikr_high
	! 				read(101,'(E16.8,E16.8)', advance='no') currcnt%Psi1_Em(ikr,iX,iKcm)
	! 			enddo
	! 		enddo
	!
	! 		read(100,'(E16.8)')
	! 		read(101,'(E16.8)')
	! 	enddo
	! 	close(100)
	! 	close(101)
	!
	! 	!make sure the exciton wavefunctions are normalized
	! 	do iX=1,currcnt%ikr_high-currcnt%ikr_low+1
	! 		do iKcm=currcnt%iKcm_min_fine,currcnt%iKcm_max_fine
	! 			tmpr1 = 0.d0
	! 			tmpr2 = 0.d0
	! 			tmpr3 = 0.d0
	! 			tmpr4 = 0.d0
	! 			do ikr=currcnt%ikr_low,currcnt%ikr_high
	! 				tmpr1 = tmpr1 + abs(currcnt%Psi0_Ep(ikr,iX,iKcm))
	! 				tmpr2 = tmpr2 + abs(currcnt%Psi0_Em(ikr,iX,iKcm))
	! 				tmpr3 = tmpr3 + abs(currcnt%Psi1_Ep(ikr,iX,iKcm))
	! 				tmpr4 = tmpr4 + abs(currcnt%Psi1_Em(ikr,iX,iKcm))
	! 			enddo
	! 			currcnt%Psi0_Ep(:,iX,iKcm) = currcnt%Psi0_Ep(:,iX,iKcm) / dcmplx(sqrt(tmpr1))
	! 			currcnt%Psi0_Em(:,iX,iKcm) = currcnt%Psi0_Em(:,iX,iKcm) / dcmplx(sqrt(tmpr2))
	! 			currcnt%Psi1_Ep(:,iX,iKcm) = currcnt%Psi1_Ep(:,iX,iKcm) / dcmplx(sqrt(tmpr3))
	! 			currcnt%Psi1_Em(:,iX,iKcm) = currcnt%Psi1_Em(:,iX,iKcm) / dcmplx(sqrt(tmpr4))
	! 		enddo
	! 	enddo
	!
	!
	! 	! set the information of the target exciton type
	! 	select case (trim(currcnt%target_exciton_type))
	! 	case ('Ex0_Ep')
	! 		call write_log("Target exciton: Ex0_Ep")
	! 		allocate(currcnt%Ex_t(1:currcnt%ikr_high-currcnt%ikr_low+1,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
	! 		allocate(currcnt%Psi_t(currcnt%ikr_low:currcnt%ikr_high,1:currcnt%ikr_high-currcnt%ikr_low+1,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
	! 		currcnt%Ex_t = currcnt%Ex0_Ep
	! 		currcnt%Psi_t = currcnt%Psi0_Ep
	! 	case ('Ex0_Em')
	! 		call write_log("Target exciton: Ex0_Em")
	! 		allocate(currcnt%Ex_t(1:currcnt%ikr_high-currcnt%ikr_low+1,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
	! 		allocate(currcnt%Psi_t(currcnt%ikr_low:currcnt%ikr_high,1:currcnt%ikr_high-currcnt%ikr_low+1,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
	! 		currcnt%Ex_t = currcnt%Ex0_Em
	! 		currcnt%Psi_t = currcnt%Psi0_Em
	! 	case ('Ex1_Ep')
	! 		call write_log("Target exciton: Ex1_Ep")
	! 		allocate(currcnt%Ex_t(1:currcnt%ikr_high-currcnt%ikr_low+1,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
	! 		allocate(currcnt%Psi_t(currcnt%ikr_low:currcnt%ikr_high,1:currcnt%ikr_high-currcnt%ikr_low+1,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
	! 		currcnt%Ex_t = currcnt%Ex1_Ep
	! 		currcnt%Psi_t = currcnt%Psi1_Ep
	! 	case ('Ex1_Em')
	! 		call write_log("Target exciton: Ex1_Em")
	! 		allocate(currcnt%Ex_t(1:currcnt%ikr_high-currcnt%ikr_low+1,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
	! 		allocate(currcnt%Psi_t(currcnt%ikr_low:currcnt%ikr_high,1:currcnt%ikr_high-currcnt%ikr_low+1,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
	! 		currcnt%Ex_t = currcnt%Ex1_Em
	! 		currcnt%Psi_t = currcnt%Psi1_Em
	! 	end select
	!
	! 	deallocate(currcnt%Psi0_Ep)
	! 	deallocate(currcnt%Psi1_Ep)
	! 	deallocate(currcnt%Psi0_Em)
	! 	deallocate(currcnt%Psi1_Em)
	!
	! 	return
	! end subroutine input_E_exciton

end module input_cnt_mod
