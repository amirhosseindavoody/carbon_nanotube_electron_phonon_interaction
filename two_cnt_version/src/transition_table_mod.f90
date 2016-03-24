!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate transition table
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
module transition_table_mod
	implicit none
	private
	public :: calculate_transition_table

	real*8, dimension(:,:,:), allocatable :: transitionRate	

	real*8, public :: c2cMin
	real*8, public :: c2cMax
	integer, public :: nc2c
	real*8 :: dc2c
	
	real*8, public :: thetaMax
	real*8, public :: thetaMin
	integer, public :: nTheta
	real*8 :: dTheta

	integer, public :: partition_function_type

	complex*16, dimension(:), allocatable :: kSpaceMatrixElement_sameEnergy, kSpaceMatrixElement_crossoingPoints

	abstract interface 
		subroutine calculate_kSpaceMatrixElement(nTransitionPoints, transitionPoints, kSpaceMatrixElement)
			integer, intent(in) :: nTransitionPoints
			integer, dimension(nTransitionPoints,4), intent(in) :: transitionPoints
			complex*16, dimension(:), allocatable, intent(inout) :: kSpaceMatrixElement
		end subroutine calculate_kSpaceMatrixElement
	end interface
	
contains
	!**************************************************************************************************************************
	! calculate transition table
	!**************************************************************************************************************************
	
	subroutine calculate_transition_table (cnt1,cnt2)
		use a2a_kspace_matrix_element_mod, only: calculate_a2a_kSpaceMatrixElement
		use a2ep_kspace_matrix_element_mod, only: calculate_a2ep_kSpaceMatrixElement
		use a2em_kspace_matrix_element_mod, only: calculate_a2em_kSpaceMatrixElement
		use cnt_class, only: cnt
		use comparams, only: ppLen, Temperature
		use em2a_kspace_matrix_element_mod, only: calculate_em2a_kSpaceMatrixElement
		use em2ep_kspace_matrix_element_mod, only: calculate_em2ep_kSpaceMatrixElement
		use em2em_kspace_matrix_element_mod, only: calculate_em2em_kSpaceMatrixElement
		use ep2a_kspace_matrix_element_mod, only: calculate_ep2a_kSpaceMatrixElement
		use ep2ep_kspace_matrix_element_mod, only: calculate_ep2ep_kSpaceMatrixElement
		use ep2em_kspace_matrix_element_mod, only: calculate_ep2em_kSpaceMatrixElement
		use geometric_matrix_element_mod, only: calculate_finite_geometric_matrix_element, calculate_infinite_geometric_matrix_element, calculate_infinite_parallel_geometric_matrix_element
		use physicalConstants, only: pi, kb, hb, A_u
		use prepareForster_module, only: calculatePartitionFunction, calculateDOS
		use rotate_shift_mod, only: rotate_shift_cnt
		use transition_points_mod, only: findSameEnergy, sameEnergy, findCrossings, crossingPoints
		use write_log_mod, only: writeLog

		type(cnt), intent(inout) :: cnt1,cnt2
		character(len=200) :: logInput

		integer :: nSameEnergy, iS
		integer :: nCrossing, iC
		integer :: ix1, ix2, iKcm1, iKcm2
		real*8 :: partitionFunction1, partitionFunction2
		real*8 :: dos1, dos2
		complex*16 :: matrixElement, geometricMatrixElement
		
		integer :: ic2c
		real*8 :: c2cDistance
		integer :: iTheta
		real*8 :: theta

		procedure(calculate_kSpaceMatrixElement), pointer :: k_space_melement_ptr => null()
		
		call writeLog(new_line('A')//"************** Start calculating transitionTable ****************")
		
		! set seperation properties
		if (nc2c .ne. 1) then
			dc2c = (c2cMax-c2cMin)/dble(nc2c-1)
		else
			dc2c = 0.d0
		end if

		write(logInput, '("c2cMin[nm] = ", F4.1)') c2cMin*1.d9
		call writeLog(trim(logInput))
		write(logInput, '("c2cMax[nm] = ", F4.1)') c2cMax*1.d9
		call writeLog(trim(logInput))
		write(logInput, '("nC2C = ", I3.3)') nC2C
		
		! set orientation properties
		if (nTheta .ne. 1) then
			dTheta = (thetaMax-thetaMin)/dble(nTheta-1)
		else
			dTheta = 0.d0
		end if

		call writeLog(trim(logInput))
		write(logInput, '("thetaMin = ", I3.3)') nint(thetaMin*180/pi)
		call writeLog(trim(logInput))
		write(logInput, '("thetaMax = ", I3.3)') nint(thetaMax*180/pi)
		call writeLog(trim(logInput))
		write(logInput, '("nTheta = ", I3.3)') nTheta
		call writeLog(trim(logInput))
	
		!calculate the crossing points and points with the same energy between cnt1 and cnt2
		call findSameEnergy(cnt1,cnt2)
		call findCrossings(cnt1,cnt2)
		
		select case (trim(cnt1%targetExcitonType))
		case('Ex_A1', 'Ex0_A2', 'Ex1_A2')
			select case (trim(cnt2%targetExcitonType))
			case('Ex_A1', 'Ex0_A2', 'Ex1_A2')
				k_space_melement_ptr => calculate_a2a_kSpaceMatrixElement
			case('Ex0_Ep', 'Ex1_Ep')
				k_space_melement_ptr => calculate_a2ep_kSpaceMatrixElement
			case('Ex0_Em', 'Ex1_Em')
				k_space_melement_ptr => calculate_a2em_kSpaceMatrixElement
			end select
		case('Ex0_Ep', 'Ex1_Ep')
			select case (trim(cnt2%targetExcitonType))
			case('Ex_A1', 'Ex0_A2', 'Ex1_A2')
				k_space_melement_ptr => calculate_ep2a_kSpaceMatrixElement
			case('Ex0_Ep', 'Ex1_Ep')
				k_space_melement_ptr => calculate_ep2ep_kSpaceMatrixElement
			case('Ex0_Em', 'Ex1_Em')
				k_space_melement_ptr => calculate_ep2em_kSpaceMatrixElement
			end select
		case('Ex0_Em', 'Ex1_Em')
			select case (trim(cnt2%targetExcitonType))
			case('Ex_A1', 'Ex0_A2', 'Ex1_A2')
				k_space_melement_ptr => calculate_em2a_kSpaceMatrixElement
			case('Ex0_Ep', 'Ex1_Ep')
				k_space_melement_ptr => calculate_em2ep_kSpaceMatrixElement
			case('Ex0_Em', 'Ex1_Em')
				k_space_melement_ptr => calculate_em2em_kSpaceMatrixElement
			end select
		end select

		call k_space_melement_ptr(size(sameEnergy,1), sameEnergy, kSpaceMatrixElement_sameEnergy)
		call k_space_melement_ptr(size(crossingPoints,1), crossingPoints, kSpaceMatrixElement_crossoingPoints)

		!allocate the transition rate table
		allocate(transitionRate(2,nTheta,nc2c))
		transitionRate = 0.d0

		call calculatePartitionFunction(partition_function_type, cnt1, partitionFunction1)
		call calculatePartitionFunction(partition_function_type, cnt2, partitionFunction2)

		do ic2c = 1, nc2c
			c2cDistance = c2cMin+dble(ic2c-1)*dc2c

			do iTheta = 1, nTheta				
				theta = thetaMin + dble(iTheta-1)*dTheta
				
				! calculate exciton transfer rate for finite length CNTs				
				if ((cnt1%length .lt. huge(1.d0)) .and. (cnt2%length .lt. huge(1.d0))) then
					write(logInput,*) 'Calculating finite transition rate: iTheta=', iTheta, ', nTheta=', nTheta, 'iC2C=', ic2c, ', nC2C=', nc2c
					call writeLog(logInput)

					nSameEnergy = size(sameEnergy,1)

					call rotate_shift_cnt(cnt1, 0.d0, 0.d0)
					call rotate_shift_cnt(cnt2, theta, c2cDistance)
									
					do iS = 1,nSameEnergy
													
						ix1 = sameEnergy(iS,1)
						ix2 = sameEnergy(iS,2)
						iKcm1 = sameEnergy(iS,3)
						iKcm2 = sameEnergy(iS,4)

						call calculate_finite_geometric_matrix_element(iKcm1, iKcm2, geometricMatrixElement)
						
						matrixElement = geometricMatrixElement * kSpaceMatrixElement_sameEnergy(iS)
						call calculateDOS(cnt1,iKcm1,ix1,dos1)
						call calculateDOS(cnt2,iKcm2,ix2,dos2)

						transitionRate(1,iTheta,ic2c) = transitionRate(1,iTheta,ic2c) + exp(-(cnt1%Ex_t(ix1,iKcm1))/kb/Temperature) * (abs(matrixElement)**2) * dos2 * (A_u**2/(4.d0*pi*pi*cnt1%radius*cnt2%radius))**2 / hb / cnt1%length / partitionFunction1 
						transitionRate(2,iTheta,ic2c) = transitionRate(2,iTheta,ic2c) + exp(-(cnt2%Ex_t(ix2,iKcm2))/kb/Temperature) * (abs(matrixElement)**2) * dos1 * (A_u**2/(4.d0*pi*pi*cnt1%radius*cnt2%radius))**2 / hb / cnt2%length / partitionFunction2

					end do

				! calculate exciton transfer rate for infinitely long CNTs
				else
					write(logInput,*) 'Calculating infinite transition rate: iTheta=', iTheta, ', nTheta=', nTheta, 'iC2C=', ic2c, ', nC2C=', nc2c
					call writeLog(logInput)

					if (theta .eq. 0.d0) then
						nCrossing = size(crossingPoints,1)

						do iC = 1,nCrossing
							
							ix1 = crossingPoints(iC,1)
							ix2 = crossingPoints(iC,2)
							iKcm1 = crossingPoints(iC,3)
							iKcm2 = crossingPoints(iC,4)

							call calculate_infinite_parallel_geometric_matrix_element(iKcm1,iKcm2,c2cDistance,geometricMatrixElement)

							matrixElement = geometricMatrixElement * kSpaceMatrixElement_crossoingPoints(iC)
							call calculateDOS(cnt1,iKcm1,ix1,dos1)
							call calculateDOS(cnt2,iKcm2,ix2,dos2)

							transitionRate(1,iTheta,ic2c) = transitionRate(1,iTheta,ic2c) + exp(-(cnt1%Ex_t(ix1,iKcm1))/kb/Temperature) * (abs(matrixElement)**2) * dos1 * (A_u**2/(4.d0*pi*pi*cnt1%radius*cnt2%radius))**2 * (2*pi/cnt1%dkx) / hb / partitionFunction1
							transitionRate(2,iTheta,ic2c) = transitionRate(2,iTheta,ic2c) + exp(-(cnt2%Ex_t(ix2,iKcm2))/kb/Temperature) * (abs(matrixElement)**2) * dos2 * (A_u**2/(4.d0*pi*pi*cnt1%radius*cnt2%radius))**2 * (2*pi/cnt2%dkx) / hb / partitionFunction2

						end do
					
					else
						nSameEnergy = size(sameEnergy,1)
									
						do iS = 1,nSameEnergy
						
							ix1 = sameEnergy(iS,1)
							ix2 = sameEnergy(iS,2)
							iKcm1 = sameEnergy(iS,3)
							iKcm2 = sameEnergy(iS,4)

							call calculate_infinite_geometric_matrix_element(iKcm1, iKcm2, theta, c2cDistance, geometricMatrixElement)
						
							matrixElement = geometricMatrixElement * kSpaceMatrixElement_sameEnergy(iS)
							call calculateDOS(cnt1,iKcm1,ix1,dos1)
							call calculateDOS(cnt2,iKcm2,ix2,dos2)

							transitionRate(1,iTheta,ic2c) = transitionRate(1,iTheta,ic2c) + exp(-(cnt1%Ex_t(ix1,iKcm1))/kb/Temperature) * (abs(matrixElement)**2) * dos2 * (A_u**2/(4.d0*pi*pi*cnt1%radius*cnt2%radius))**2 * (sin(theta)) / hb / ppLen/ partitionFunction1
							transitionRate(2,iTheta,ic2c) = transitionRate(2,iTheta,ic2c) + exp(-(cnt2%Ex_t(ix2,iKcm2))/kb/Temperature) * (abs(matrixElement)**2) * dos1 * (A_u**2/(4.d0*pi*pi*cnt1%radius*cnt2%radius))**2 * (sin(theta)) / hb / ppLen/ partitionFunction2

						end do
					endif

				end if

				call save_transition_rates(iTheta,ic2c)

			end do
		end do

		k_space_melement_ptr => null()
		
		call saveTransitionRates()
		return				
	end subroutine calculate_transition_table
	
	!**************************************************************************************************************************
	! save the calculated transition table
	!**************************************************************************************************************************
	
	subroutine saveTransitionRates()

		integer :: ic2c
		integer :: iTheta

		!write transition rates to the file
		open(unit=100,file='transitionRates12.dat',status="unknown")
		do ic2c = 1,nc2c
			do iTheta=1,nTheta
				write(100,10, advance='no') transitionRate(1,iTheta,ic2c)
			end do
			write(100,10)
		end do
		close(100)
		
		open(unit=100,file='transitionRates21.dat',status="unknown")
		do ic2c = 1,nc2c
			do iTheta=1,nTheta
				write(100,10, advance='no') transitionRate(2,iTheta,ic2c)
			end do
			write(100,10)
		end do
		close(100)
				
		open(unit=100,file='theta.dat',status="unknown")
		do iTheta=1,nTheta
			write(100,10, advance='no') thetaMin+dble(iTheta-1)*dTheta
		enddo
		write(100,10)
				
		open(unit=100,file='c2c.dat',status="unknown")
		do ic2c=1,nc2c
			write(100,10, advance='no') c2cMin+dble(ic2c-1)*dc2c
		enddo
		close(100)
				
10		FORMAT (E16.8)
		
		return
	end subroutine saveTransitionRates

	!**************************************************************************************************************************
	! save the calculated transition table on the fly
	!**************************************************************************************************************************
	
	subroutine save_transition_rates(iTheta, ic2c)

		integer, intent(in) :: ic2c
		integer, intent(in) :: iTheta
		logical :: flgexist

		! write 1 to 2 transition rates
		inquire(file="transition_rates_12.dat",exist=flgexist)
		if (flgexist) then
			open(unit=100, file="transition_rates_12.dat", status="old", position="append", action="write")
		else
			open(unit=100, file="transition_rates_12.dat", status="new", action="write")
		end if

		if (iTheta .ne. nTheta) then
			write(100,'(E16.8)', advance='no') transitionRate(1,iTheta,ic2c)
		else
			write(100,'(E16.8)') transitionRate(1,iTheta,ic2c)
		endif
		close(100)

		! write 2 to 1 transition rates
		inquire(file="transition_rates_21.dat",exist=flgexist)
		if (flgexist) then
			open(unit=100, file="transition_rates_21.dat", status="old", position="append", action="write")
		else
			open(unit=100, file="transition_rates_21.dat", status="new", action="write")
		end if

		if (iTheta .ne. nTheta) then
			write(100,'(E16.8)', advance='no') transitionRate(2,iTheta,ic2c)
		else
			write(100,'(E16.8)') transitionRate(2,iTheta,ic2c)
		endif
		close(100)

		! write theta values
		inquire(file="theta.dat",exist=flgexist)
		if (flgexist) then
			open(unit=100, file="theta.dat", status="old", position="append", action="write")
		else
			open(unit=100, file="theta.dat", status="new", action="write")
		end if

		write(100,'(E16.8)', advance='no') thetaMin+dble(iTheta-1)*dTheta
		close(100)
				
		! write c2c values
		inquire(file="c2c.dat",exist=flgexist)
		if (flgexist) then
			open(unit=100, file="c2c.dat", status="old", position="append", action="write")
		else
			open(unit=100, file="c2c.dat", status="new", action="write")
		end if

		write(100,'(E16.8)', advance='no') c2cMin+dble(ic2c-1)*dc2c
		close(100)
		
		return
	end subroutine save_transition_rates

end module transition_table_mod