!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate kappa matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
module kappaMatrix_module
	implicit none
	private
	public :: calculateKappaMatrix
	
	real*8, dimension(:,:), allocatable :: kappaMatrix

	real*8 :: c2cDistance !center to center distance between parallel carbon nanotubes
	
	real*8 :: theta
	real*8 :: thetaMax
	real*8 :: thetaMin
	integer :: nTheta
	real*8 :: dTheta
	
contains
	!**************************************************************************************************************************
	! calculate kappa matrix
	!**************************************************************************************************************************
	subroutine calculateKappaMatrix (cnt1,cnt2)
		use cnt_class, only: cnt
		use comparams, only: ppLen
		use parallel_geometry_mod, only: calculateParallelGeometryRate
		use physicalConstants, only: pi
		use transition_points_mod, only: findCrossings, findSameEnergy
		use write_log_mod, only: writeLog

		
		type(cnt), intent(in) :: cnt1,cnt2
		integer :: iTheta1, iTheta2
		character(len=100) :: logInput
		
		call writeLog("Calculating kappa matrix ...")
		
		! set seperation properties
		c2cDistance = 5.0d-9
		ppLen = 40.0d-9
		
		! set orientation properties
		nTheta = 10
		thetaMax = pi
		thetaMin = 0.d0
		if (nTheta .ne. 1) then
			dTheta = (thetaMax-thetaMin)/dble(nTheta)
		else
			dTheta = 0.d0
		end if
		
		ppLen  = ppLen*dble(nTheta-1)

		!! calculate kappa matrix for only one type of CNT in the system
		!allocate(kappaMatrix(nTheta, nTheta))
		!kappaMatrix = 0.d0 * kappaMatrix
		!
		!! calculate crossing points and same energy points between cnt1 and cnt2
		!call findCrossings(cnt1,cnt2)
		!call findSameEnergy(cnt1,cnt2)
		!
		!! calculate scattering between cnt1 and cnt2
		!do iTheta1 = 1, nTheta
		!	do iTheta2 = 1, nTheta
		!		theta = abs(dble(iTheta1-iTheta2))*dTheta
		!		if ((iTheta1-iTheta2) .ne. 0) then
		!			!calculate unparallel rate
		!		end if
		!		print *, 'iTheta1=', iTheta1, ', iTheta2=', iTheta2
		!	end do
		!end do
		
		
		! calculate kappa matrix for two types of CNT in the system
		allocate(kappaMatrix(2*nTheta, 2*nTheta))
		kappaMatrix = 0.d0 * kappaMatrix
		
		! calculate crossing points and same energy points between cnt1 and cnt1
		call findCrossings(cnt1,cnt1)
		call findSameEnergy(cnt1,cnt1)
		
		! calculate scattering between cnt1 and cnt1
		do iTheta1 = 1, nTheta
			do iTheta2 = 1, nTheta
				theta = abs(dble(iTheta1-iTheta2))*dTheta
				if ((iTheta1-iTheta2) .ne. 0) then
					! calculate unparallel rate
				end if
				write(logInput,*) 'iTheta1=', iTheta1, ', iTheta2=', iTheta2
				call writeLog(logInput)
			end do
		end do
		
		! calculate crossing points and same energy points between cnt2 and cnt2
		call findCrossings(cnt2,cnt2)
		call findSameEnergy(cnt2,cnt2)
		
		! calculate scattering between cnt2 and cnt2
		do iTheta1 = 1, nTheta
			do iTheta2 = 1, nTheta
				theta = abs(dble(iTheta1-iTheta2))*dTheta
				if ((iTheta1-iTheta2) .ne. 0) then
					! calculate unparallel rate
				end if
				write(logInput,*) 'iTheta1=', iTheta1, ', iTheta2=', iTheta2
				call writeLog(logInput)
			end do
		end do
		
		! calculate crossing points and same energy points between cnt1 and cnt2
		call findCrossings(cnt1,cnt2)
		call findSameEnergy(cnt1,cnt2)			
		
		! calculate scattering between cnt1 and cnt2
		do iTheta1 = 1, nTheta
			do iTheta2 = 1, nTheta
				theta = abs(dble(iTheta1-iTheta2))*dTheta
				if ((iTheta1-iTheta2) .eq. 0) then
					call calculateParallelGeometryRate(cnt1,cnt2, kappaMatrix(iTheta1,nTheta+iTheta2), kappaMatrix(nTheta+iTheta2,iTheta1), c2cDistance)
				else
					! calculate unparallel rate
				end if
				write(logInput,*) 'iTheta1=', iTheta1, ', iTheta2=', iTheta2
				call writeLog(logInput)
			end do
		end do
		
		call saveKappaMatrix()
		
		return
	end subroutine calculateKappaMatrix

	!**************************************************************************************************************************
	! save kappa matrix
	!**************************************************************************************************************************
	subroutine saveKappaMatrix()
		integer :: iTheta1, iTheta2
		integer :: nKappaMatrix
		
		nKappaMatrix = size(kappaMatrix,1)
		
		!write transition rates to the file
		open(unit=100,file='kappaMatrix.dat',status="unknown")
		do iTheta1 = 1,nKappaMatrix
			do iTheta2 = 1,nKappaMatrix
				write(100,'(E16.8)', advance='no') kappaMatrix(iTheta1,iTheta2)
			end do
			write(100,'(E16.8)')
		end do
		close(100)

		return
	end subroutine saveKappaMatrix
	
end module kappaMatrix_module