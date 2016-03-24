!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Declaration of input parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
module input_cnt_mod
	implicit none
	private
	public :: input_cnt

contains	
	!**************************************************************************************************************************
	! read the CNT information from the log file of exciton calculation which is done previously
	!**************************************************************************************************************************

	subroutine input_cnt(currcnt)
		use cnt_class, only: cnt, cnt_geometry, cnt_band
		use physicalConstants, only: eV
		use write_log_mod, only: writeLog

		type(cnt) :: currcnt
		character(len=200) :: buffer, label
		integer :: ios=0
		integer :: pos=0
		integer, parameter :: nparam=12
		integer :: iparam=0
		character(len=100) :: logInput

		open(unit=100, file= trim(currcnt%directory)//'log.dat', status="old")

		iparam=0
		ios=0
		do while ((ios == 0) .and. (iparam .lt. nparam))
			read (100,'(A)',iostat=ios) buffer
			if (ios == 0) then
				pos = scan(buffer,'=')
				label = buffer(1:pos-1)
				buffer = buffer(pos+1:)
				label = adjustl(label)
				buffer = adjustl(buffer)

				select case (trim(label))
				case ('n_ch')
					read(buffer, *, iostat=ios) currcnt%n_ch
					write(logInput,"('n_ch = ', I2.2)") currcnt%n_ch
					call writeLog(logInput)
					iparam = iparam+1
				case ('m_ch')
					read(buffer, *, iostat=ios) currcnt%m_ch
					write(logInput,"('m_ch = ', I2.2)") currcnt%m_ch
					call writeLog(logInput)
					iparam = iparam+1
				case ('nkg')
					read(buffer, *, iostat=ios) currcnt%nkg
					write(logInput,"('nkg = ', I4.4)") currcnt%nkg
					call writeLog(logInput)
					iparam = iparam+1
				case ('dk/dkx')
					read(buffer, *, iostat=ios) currcnt%dk_dkx_ratio
					write(logInput,"('dk/dkx = ', I4.4)") currcnt%dk_dkx_ratio
					call writeLog(logInput)
					iparam = iparam+1
				case ('nr')
					read(buffer, *, iostat=ios) currcnt%nr
					write(logInput,"('nr = ', I4.4)") currcnt%nr
					call writeLog(logInput)
					iparam = iparam+1
				case ('E_th[eV]')
					read(buffer, *, iostat=ios) currcnt%E_th
					write(logInput,"('E_th[eV] = ', f3.1)") currcnt%E_th
					call writeLog(logInput)
					iparam = iparam+1
					currcnt%E_th=currcnt%E_th*eV
				case ('Kcm_max[1/nm]')
					read(buffer, *, iostat=ios) currcnt%Kcm_max
					write(logInput,"('Kcm_max[1/nm] = ', f3.1)") currcnt%Kcm_max
					call writeLog(logInput)
					currcnt%Kcm_max = currcnt%Kcm_max*1.d9
					iparam = iparam+1
				case ('i_sub')
					read(buffer, *, iostat=ios) currcnt%i_sub
					write(logInput,"('i_sub = ', I1.1)") currcnt%i_sub
					call writeLog(logInput)
					iparam = iparam+1
				case ('Ckappa')
					read(buffer, *, iostat=ios) currcnt%Ckappa
					write(logInput,"('Ckappa = ', f3.1)") currcnt%Ckappa
					call writeLog(logInput)
					iparam = iparam+1
				case ('kappa')
					read(buffer, *, iostat=ios) currcnt%kappa
					write(logInput,"('kappa = ', f3.1)") currcnt%kappa
					call writeLog(logInput)
					iparam = iparam+1
				case ('nX_a')
					read(buffer, *, iostat=ios) currcnt%nX_a
					write(logInput,"('nX_a = ', I3.3)") currcnt%nX_a
					call writeLog(logInput)
					iparam = iparam+1
				case ('nX_e')
					read(buffer, *, iostat=ios) currcnt%nX_e
					write(logInput,"('nX_e = ', I3.3)") currcnt%nX_e
					call writeLog(logInput)
					iparam = iparam+1
				end select
			end if
		end do

		close(100)

		if (iparam .lt. nparam) then
			call writeLog("Error in reading parameters!")
			call exit()
		end if

		! create the cnt object and calculate bands and load exciton wavefunction
		call cnt_geometry(currcnt)
		call cnt_band(currcnt)
		currcnt%nX_a = currcnt%ikr_high-currcnt%ikr_low+1
		currcnt%nX_e = currcnt%ikr_high-currcnt%ikr_low+1

		call input_A_exciton(currcnt)
		call input_E_exciton(currcnt)

		if (.not. allocated(currcnt%Ex_t)) then
			call writeLog("Error in setting target exciton type!!!!")
			call exit()
		endif

		currcnt%nX_t = size(currcnt%Ex_t,1)

		return
	end subroutine input_cnt

	!**************************************************************************************************************************
	! load the A-type exciton wavefunction and energies from the ExcitonEnergy calculation
	!**************************************************************************************************************************
	
	subroutine input_A_exciton(currcnt)
		use cnt_class, only: cnt
		use write_log_mod, only: writeLog

		type(cnt), intent(inout) :: currcnt
		integer :: iX, iKcm, ikr
		real*8 :: tmpr1, tmpr2, tmpr3

		! read the information of A1 exciton
		allocate(currcnt%Ex_A1(1:currcnt%nX_a,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
		allocate(currcnt%Psi_A1(currcnt%ikr_low:currcnt%ikr_high,1:currcnt%nX_a,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
		
		open(unit=100,file=trim(currcnt%directory)//'Ex_A1.dat',status="old")
		open(unit=101,file=trim(currcnt%directory)//'Psi_A1.dat',status="old")
		do iKcm=currcnt%iKcm_min_fine,currcnt%iKcm_max_fine
			do iX=1,currcnt%nX_a
				read(100,'(E16.8)', advance='no') currcnt%Ex_A1(iX,iKcm)
				do ikr=currcnt%ikr_low,currcnt%ikr_high
					read(101,'(E16.8,E16.8)', advance='no') currcnt%Psi_A1(ikr,iX,iKcm)
				enddo
			enddo
			
			read(100,'(E16.8)')
			read(101,'(E16.8)')
		enddo
		close(100)
		close(101)

		! read the information of singlet A2 exciton
		allocate(currcnt%Ex0_A2(1:currcnt%nX_a,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
		allocate(currcnt%Psi0_A2(currcnt%ikr_low:currcnt%ikr_high,1:currcnt%nX_a,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
		
		open(unit=100,file=trim(currcnt%directory)//'Ex0_A2.dat',status="old")
		open(unit=101,file=trim(currcnt%directory)//'Psi0_A2.dat',status="old")
		do iKcm=currcnt%iKcm_min_fine,currcnt%iKcm_max_fine
			do iX=1,currcnt%nX_a
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

		! read the information of triplet A2 exciton
		allocate(currcnt%Ex1_A2(1:currcnt%nX_a,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
		allocate(currcnt%Psi1_A2(currcnt%ikr_low:currcnt%ikr_high,1:currcnt%nX_a,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
		
		open(unit=100,file=trim(currcnt%directory)//'Ex1_A2.dat',status="old")
		open(unit=101,file=trim(currcnt%directory)//'Psi1_A2.dat',status="old")
		do iKcm=currcnt%iKcm_min_fine,currcnt%iKcm_max_fine
			do iX=1,currcnt%nX_a
				read(100,'(E16.8)', advance='no') currcnt%Ex1_A2(iX,iKcm)
				do ikr=currcnt%ikr_low,currcnt%ikr_high
					read(101,'(E16.8,E16.8)', advance='no') currcnt%Psi1_A2(ikr,iX,iKcm)
				enddo
			enddo
			
			read(100,'(E16.8)')
			read(101,'(E16.8)')
		enddo
		close(100)
		close(101)

		!make sure the exciton wavefunctions are normalized
		do iX=1,currcnt%nX_a
			do iKcm=currcnt%iKcm_min_fine,currcnt%iKcm_max_fine
				tmpr1 = 0.d0
				tmpr2 = 0.d0
				tmpr3 = 0.d0
				do ikr=currcnt%ikr_low,currcnt%ikr_high
					tmpr1 = tmpr1 + abs(currcnt%Psi_A1(ikr,iX,iKcm))
					tmpr2 = tmpr2 + abs(currcnt%Psi0_A2(ikr,iX,iKcm))
					tmpr3 = tmpr3 + abs(currcnt%Psi1_A2(ikr,iX,iKcm))
				enddo
				currcnt%Psi_A1(:,iX,iKcm) = currcnt%Psi_A1(:,iX,iKcm) / dcmplx(sqrt(tmpr1))
				currcnt%Psi0_A2(:,iX,iKcm) = currcnt%Psi0_A2(:,iX,iKcm) / dcmplx(sqrt(tmpr2))
				currcnt%Psi1_A2(:,iX,iKcm) = currcnt%Psi1_A2(:,iX,iKcm) / dcmplx(sqrt(tmpr3))
			enddo
		enddo


		! set the information of the target exciton type
		select case (trim(currcnt%targetExcitonType))
		case ('Ex_A1')
			call writeLog("Target exciton: Ex_A1")
			allocate(currcnt%Ex_t(1:currcnt%nX_a,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
			allocate(currcnt%Psi_t(currcnt%ikr_low:currcnt%ikr_high,1:currcnt%nX_a,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
			currcnt%Ex_t = currcnt%Ex_A1
			currcnt%Psi_t = currcnt%Psi_A1
			currcnt%ex_symmetry = -1.d0
		case ('Ex0_A2')
			call writeLog("Target exciton: Ex0_A2")
			allocate(currcnt%Ex_t(1:currcnt%nX_a,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
			allocate(currcnt%Psi_t(currcnt%ikr_low:currcnt%ikr_high,1:currcnt%nX_a,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
			currcnt%Ex_t = currcnt%Ex0_A2
			currcnt%Psi_t = currcnt%Psi0_A2
			currcnt%ex_symmetry = +1.d0
		case ('Ex1_A2')
			call writeLog("Target exciton: Ex1_A2")
			allocate(currcnt%Ex_t(1:currcnt%nX_a,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
			allocate(currcnt%Psi_t(currcnt%ikr_low:currcnt%ikr_high,1:currcnt%nX_a,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
			currcnt%Ex_t = currcnt%Ex1_A2
			currcnt%Psi_t = currcnt%Psi1_A2
			currcnt%ex_symmetry = +1.d0
		end select

		deallocate(currcnt%Psi_A1)
		deallocate(currcnt%Psi0_A2)
		deallocate(currcnt%Psi1_A2)

		return
	end subroutine input_A_exciton

	!**************************************************************************************************************************
	! load the E-type exciton wavefunction and energies from the ExcitonEnergy calculation
	!**************************************************************************************************************************
	
	subroutine input_E_exciton(currcnt)
		use cnt_class, only: cnt
		use write_log_mod, only: writeLog

		type(cnt), intent(inout) :: currcnt
		integer :: iX, iKcm, ikr
		real*8 :: tmpr1, tmpr2, tmpr3, tmpr4

		! read the information of singlet E+ exciton
		allocate(currcnt%Ex0_Ep(1:currcnt%nX_e,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
		allocate(currcnt%Psi0_Ep(currcnt%ikr_low:currcnt%ikr_high,1:currcnt%nX_e,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
		
		open(unit=100,file=trim(currcnt%directory)//'Ex0_Ep.dat',status="old")
		open(unit=101,file=trim(currcnt%directory)//'Psi0_Ep.dat',status="old")
		do iKcm=currcnt%iKcm_min_fine,currcnt%iKcm_max_fine
			do iX=1,currcnt%nX_e
				read(100,'(E16.8)', advance='no') currcnt%Ex0_Ep(iX,iKcm)
				do ikr=currcnt%ikr_low,currcnt%ikr_high
					read(101,'(E16.8,E16.8)', advance='no') currcnt%Psi0_Ep(ikr,iX,iKcm)
				enddo
			enddo
			
			read(100,'(E16.8)')
			read(101,'(E16.8)')
		enddo
		close(100)
		close(101)

		! read the information of singlet E- exciton
		allocate(currcnt%Ex0_Em(1:currcnt%nX_e,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
		allocate(currcnt%Psi0_Em(currcnt%ikr_low:currcnt%ikr_high,1:currcnt%nX_e,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
		
		open(unit=100,file=trim(currcnt%directory)//'Ex0_Em.dat',status="old")
		open(unit=101,file=trim(currcnt%directory)//'Psi0_Em.dat',status="old")
		do iKcm=currcnt%iKcm_min_fine,currcnt%iKcm_max_fine
			do iX=1,currcnt%nX_e
				read(100,'(E16.8)', advance='no') currcnt%Ex0_Em(iX,iKcm)
				do ikr=currcnt%ikr_low,currcnt%ikr_high
					read(101,'(E16.8,E16.8)', advance='no') currcnt%Psi0_Em(ikr,iX,iKcm)
				enddo
			enddo
			
			read(100,'(E16.8)')
			read(101,'(E16.8)')
		enddo
		close(100)
		close(101)

		! read the information of triplet E+ exciton
		allocate(currcnt%Ex1_Ep(1:currcnt%nX_e,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
		allocate(currcnt%Psi1_Ep(currcnt%ikr_low:currcnt%ikr_high,1:currcnt%nX_e,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
		
		open(unit=100,file=trim(currcnt%directory)//'Ex1_Ep.dat',status="old")
		open(unit=101,file=trim(currcnt%directory)//'Psi1_Ep.dat',status="old")
		do iKcm=currcnt%iKcm_min_fine,currcnt%iKcm_max_fine
			do iX=1,currcnt%nX_e
				read(100,'(E16.8)', advance='no') currcnt%Ex1_Ep(iX,iKcm)
				do ikr=currcnt%ikr_low,currcnt%ikr_high
					read(101,'(E16.8,E16.8)', advance='no') currcnt%Psi1_Ep(ikr,iX,iKcm)
				enddo
			enddo
			
			read(100,'(E16.8)')
			read(101,'(E16.8)')
		enddo
		close(100)
		close(101)

		! read the information of triplet E- exciton
		allocate(currcnt%Ex1_Em(1:currcnt%nX_e,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
		allocate(currcnt%Psi1_Em(currcnt%ikr_low:currcnt%ikr_high,1:currcnt%nX_e,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
		
		open(unit=100,file=trim(currcnt%directory)//'Ex1_Em.dat',status="old")
		open(unit=101,file=trim(currcnt%directory)//'Psi1_Em.dat',status="old")
		do iKcm=currcnt%iKcm_min_fine,currcnt%iKcm_max_fine
			do iX=1,currcnt%nX_e
				read(100,'(E16.8)', advance='no') currcnt%Ex1_Em(iX,iKcm)
				do ikr=currcnt%ikr_low,currcnt%ikr_high
					read(101,'(E16.8,E16.8)', advance='no') currcnt%Psi1_Em(ikr,iX,iKcm)
				enddo
			enddo
			
			read(100,'(E16.8)')
			read(101,'(E16.8)')
		enddo
		close(100)
		close(101)

		!make sure the exciton wavefunctions are normalized
		do iX=1,currcnt%nX_e
			do iKcm=currcnt%iKcm_min_fine,currcnt%iKcm_max_fine
				tmpr1 = 0.d0
				tmpr2 = 0.d0
				tmpr3 = 0.d0
				tmpr4 = 0.d0
				do ikr=currcnt%ikr_low,currcnt%ikr_high
					tmpr1 = tmpr1 + abs(currcnt%Psi0_Ep(ikr,iX,iKcm))
					tmpr2 = tmpr2 + abs(currcnt%Psi0_Em(ikr,iX,iKcm))
					tmpr3 = tmpr3 + abs(currcnt%Psi1_Ep(ikr,iX,iKcm))
					tmpr4 = tmpr4 + abs(currcnt%Psi1_Em(ikr,iX,iKcm))
				enddo
				currcnt%Psi0_Ep(:,iX,iKcm) = currcnt%Psi0_Ep(:,iX,iKcm) / dcmplx(sqrt(tmpr1))
				currcnt%Psi0_Em(:,iX,iKcm) = currcnt%Psi0_Em(:,iX,iKcm) / dcmplx(sqrt(tmpr2))
				currcnt%Psi1_Ep(:,iX,iKcm) = currcnt%Psi1_Ep(:,iX,iKcm) / dcmplx(sqrt(tmpr3))
				currcnt%Psi1_Em(:,iX,iKcm) = currcnt%Psi1_Em(:,iX,iKcm) / dcmplx(sqrt(tmpr4))
			enddo
		enddo


		! set the information of the target exciton type
		select case (trim(currcnt%targetExcitonType))
		case ('Ex0_Ep')
			call writeLog("Target exciton: Ex0_Ep")
			allocate(currcnt%Ex_t(1:currcnt%nX_e,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
			allocate(currcnt%Psi_t(currcnt%ikr_low:currcnt%ikr_high,1:currcnt%nX_e,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
			currcnt%Ex_t = currcnt%Ex0_Ep
			currcnt%Psi_t = currcnt%Psi0_Ep
		case ('Ex0_Em')
			call writeLog("Target exciton: Ex0_Em")
			allocate(currcnt%Ex_t(1:currcnt%nX_e,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
			allocate(currcnt%Psi_t(currcnt%ikr_low:currcnt%ikr_high,1:currcnt%nX_e,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
			currcnt%Ex_t = currcnt%Ex0_Em
			currcnt%Psi_t = currcnt%Psi0_Em
		case ('Ex1_Ep')
			call writeLog("Target exciton: Ex1_Ep")
			allocate(currcnt%Ex_t(1:currcnt%nX_e,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
			allocate(currcnt%Psi_t(currcnt%ikr_low:currcnt%ikr_high,1:currcnt%nX_e,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
			currcnt%Ex_t = currcnt%Ex1_Ep
			currcnt%Psi_t = currcnt%Psi1_Ep
		case ('Ex1_Em')
			call writeLog("Target exciton: Ex1_Em")
			allocate(currcnt%Ex_t(1:currcnt%nX_e,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
			allocate(currcnt%Psi_t(currcnt%ikr_low:currcnt%ikr_high,1:currcnt%nX_e,currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))
			currcnt%Ex_t = currcnt%Ex1_Em
			currcnt%Psi_t = currcnt%Psi1_Em
		end select

		deallocate(currcnt%Psi0_Ep)
		deallocate(currcnt%Psi1_Ep)
		deallocate(currcnt%Psi0_Em)
		deallocate(currcnt%Psi1_Em)

		return
	end subroutine input_E_exciton
	
end module input_cnt_mod