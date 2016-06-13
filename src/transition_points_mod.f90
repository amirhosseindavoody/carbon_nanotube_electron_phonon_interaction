module transition_points_mod
	implicit none
	private
    public  :: find_crossings, find_same_energy

	integer, dimension(:,:), allocatable, public :: crossing_points, same_energy

contains
	!***************************************************************************
	! find the points that the bands cross each other
	!***************************************************************************

	subroutine find_crossings(exciton_1, exciton_2)
		use cnt_class, only: exciton
		use constants_mod, only: kb
		use sim_properties_mod, only: temperature
		use write_log_mod, only: write_log, log_input

		type(exciton), pointer, intent(in) :: exciton_1,exciton_2
		integer :: ix1,ix2
		integer :: iKcm
		integer :: n_crossing
		real*8 :: rtmp1, rtmp2
		real*8 :: min_energy, deltaE

		! calculate relevant crossing points for transition from cnt1 to cnt2
		deltaE = (-1.d0) * log(1.d-3) * kb*temperature
		min_energy = minval(exciton_1%ex)

		n_crossing = 0
		do ix1 = 1,exciton_1%nx
			do ix2 = 1,exciton_2%nx
				do iKcm = 1 , min(exciton_1%iKcm_max, exciton_2%iKcm_max)
					rtmp1 = (exciton_1%ex(ix1,iKcm)-exciton_2%ex(ix2,iKcm))
					rtmp2 = (exciton_1%ex(ix1,iKcm-1)-exciton_2%ex(ix2,iKcm-1))
					if (((rtmp1 * rtmp2) .le. 0.d0) .and. ((exciton_1%ex(ix1,iKcm) - min_energy) .lt. deltaE) ) then
						n_crossing = n_crossing + 2
					end if
				end do
			end do
		end do

		write(log_input,'(A, I0)') "Number of crossing points = ", n_crossing
		call write_log(log_input)

		if(allocated(crossing_points))	deallocate(crossing_points)
		allocate(crossing_points(n_crossing ,4))

		n_crossing = 0
		do ix1 = 1,exciton_1%nx
			do ix2 = 1,exciton_2%nx
				do iKcm = 1 , min(exciton_1%iKcm_max, exciton_2%iKcm_max)
					rtmp1 = (exciton_1%ex(ix1,iKcm)-exciton_2%ex(ix2,iKcm))
					rtmp2 = (exciton_1%ex(ix1,iKcm-1)-exciton_2%ex(ix2,iKcm-1))
					if (((rtmp1 * rtmp2) .le. 0.d0) .and. ((exciton_1%ex(ix1,iKcm) - min_energy) .lt. deltaE) ) then
						n_crossing = n_crossing+1
						crossing_points(n_crossing,1) = ix1
						crossing_points(n_crossing,2) = ix2
						crossing_points(n_crossing,3) = iKcm
						crossing_points(n_crossing,4) = iKcm

						n_crossing = n_crossing+1
						crossing_points(n_crossing,1) = ix1
						crossing_points(n_crossing,2) = ix2
						crossing_points(n_crossing,3) = -iKcm
						crossing_points(n_crossing,4) = -iKcm
					end if
				end do
			end do
		end do

		call write_log(new_line('A')//"Crossing points table calculated!!!"//new_line('A'))

		! call saveTransitionPoints(cnt1,cnt2)

	end subroutine find_crossings

	!***************************************************************************
	! find the points that the bands have equal energy
	!***************************************************************************

	subroutine find_same_energy(exciton_1, exciton_2)
		use cnt_class, only: exciton
		use sim_properties_mod, only: temperature
		use math_functions_mod, only: bisect_root
		use constants_mod, only: kb
		use write_log_mod, only: write_log, log_input

		type(exciton), pointer, intent(in) :: exciton_1, exciton_2
		integer :: ix1, ix2
		integer :: iKcm1, iKcm2
		integer :: n_same_energy
		integer :: n, iKcm_raw
		real*8 :: min_energy, deltaE
		real*8, dimension(:), allocatable :: ya_tmp

		! calculate relevant same energy points for transition from exciton_1 to exciton_2
		deltaE = (-1.d0) * log(1.d-3) * kb*Temperature
		min_energy = minval(exciton_1%ex)

		n = exciton_2%iKcm_max
		allocate(ya_tmp(n))

		n_same_energy = 0
		do ix1 = 1,exciton_1%nx
			do iKcm1 = exciton_1%iKcm_min , -1
				if ((exciton_1%ex(ix1,iKcm1) - min_energy) .lt. deltaE) then
					do ix2 = 1, exciton_2%nx
						ya_tmp = exciton_2%ex(ix2,exciton_2%iKcm_min:-1)
						call bisect_root(n, ya_tmp, exciton_1%ex(ix1,iKcm1), iKcm_raw)
						if (iKcm_raw .gt. 0) then
							n_same_energy = n_same_energy + 4
						endif
					enddo
				endif
			end do
		end do

		write(log_input,'(A, I0)') "Number of same energy points = ", n_same_energy
		call write_log(log_input)

		if(allocated(same_energy))	deallocate(same_energy)
        allocate(same_energy(n_same_energy ,4))

		n_same_energy = 0
		do ix1 = 1,exciton_1%nx
			do iKcm1 = exciton_1%iKcm_min , -1
				if ((exciton_1%ex(ix1,iKcm1) - min_energy) .lt. deltaE) then
					do ix2 = 1, exciton_2%nx
						ya_tmp = exciton_2%ex(ix2,exciton_2%iKcm_min:-1)
						call bisect_root(n, ya_tmp, exciton_1%ex(ix1,iKcm1), iKcm_raw)
						if (iKcm_raw .gt. 0) then
							iKcm2 = exciton_2%iKcm_min + iKcm_raw - 1

							n_same_energy = n_same_energy + 1
							same_energy(n_same_energy, 1) = ix1
							same_energy(n_same_energy, 2) = ix2
							same_energy(n_same_energy, 3) = +iKcm1
							same_energy(n_same_energy, 4) = +iKcm2

							n_same_energy = n_same_energy + 1
							same_energy(n_same_energy, 1) = ix1
							same_energy(n_same_energy, 2) = ix2
							same_energy(n_same_energy, 3) = +iKcm1
							same_energy(n_same_energy, 4) = -iKcm2

							n_same_energy = n_same_energy + 1
							same_energy(n_same_energy, 1) = ix1
							same_energy(n_same_energy, 2) = ix2
							same_energy(n_same_energy, 3) = -iKcm1
							same_energy(n_same_energy, 4) = +iKcm2

							n_same_energy = n_same_energy + 1
							same_energy(n_same_energy, 1) = ix1
							same_energy(n_same_energy, 2) = ix2
							same_energy(n_same_energy, 3) = -iKcm1
							same_energy(n_same_energy, 4) = -iKcm2

						endif
					enddo
				endif
			end do
		end do

		call write_log(new_line('A')//"Same energy table calculated!!!"//new_line('A'))


! 		call saveTransitionPoints(cnt1,cnt2)

	end subroutine find_same_energy

	! !*************************************************************************
	! ! save CNT dispersions and the crossing points and the same energy points
	! !*************************************************************************
	!
	! subroutine saveTransitionPoints(cnt1,cnt2)
	! 	use cnt_class, only: cnt
	!
	! 	type(cnt), intent(in) :: cnt1,cnt2
	! 	integer :: iKcm, iX, i
	!
	! 	!write carbon nanotube 1 k_vector
	! 	open(unit=100,file='cnt1_kvec.dat',status="unknown")
	! 	do iKcm=cnt1%iKcm_min_fine,cnt1%iKcm_max_fine
	! 		write(100,'(E16.8)', advance='no') dble(iKcm)*cnt1%dkx
	! 	enddo
	! 	close(100)
	!
	! 	!write carbon nanotube 2 k_vector
	! 	open(unit=100,file='cnt2_kvec.dat',status="unknown")
	! 	do iKcm=cnt2%iKcm_min_fine,cnt2%iKcm_max_fine
	! 		write(100,'(E16.8)', advance='no') dble(iKcm)*cnt2%dkx
	! 	enddo
	! 	close(100)
	!
	! 	!write carbon nanotube 1 Ex_t dispersion
	! 	open(unit=100,file='cnt1_Ex_t.dat',status="unknown")
	! 	do iKcm=cnt1%iKcm_min_fine,cnt1%iKcm_max_fine
	! 		do iX=1,cnt1%nX_t
	! 			write(100,'(E16.8)', advance='no') cnt1%Ex_t(iX,iKcm)
	! 		enddo
	! 		write(100,*)
	! 	enddo
	! 	close(100)
	!
	! 	!write carbon nanotube 2 Ex_t dispersion
	! 	open(unit=100,file='cnt2_Ex_t.dat',status="unknown")
	! 	do iKcm=cnt2%iKcm_min_fine,cnt2%iKcm_max_fine
	! 		do iX=1,cnt2%nX_t
	! 			write(100,'(E16.8)', advance='no') cnt2%Ex_t(iX,iKcm)
	! 		enddo
	! 		write(100,*)
	! 	enddo
	! 	close(100)
	!
	! 	!write crossing points indexes
	! 	open(unit=100,file='crossingPoints.dat',status="unknown")
	! 	do i=lbound(crossingPoints,1),ubound(crossingPoints,1)
	! 		write(100,'(4I8, 4I8, 4I8, 4I8)') crossingPoints(i,1), crossingPoints(i,2), crossingPoints(i,3), crossingPoints(i,4)
	! 	enddo
	! 	close(100)
	!
	! 	!write same energy points indexes for transition from cnt1 to cnt2
	! 	open(unit=100,file='sameEnergy.dat',status="unknown")
	! 	do i=lbound(sameEnergy,1),ubound(sameEnergy,1)
	! 		write(100,'(4I8, 4I8, 4I8, 4I8)') sameEnergy(i,1), sameEnergy(i,2), sameEnergy(i,3), sameEnergy(i,4)
	! 	enddo
	! 	close(100)
	!
	! 	return
	! end subroutine saveTransitionPoints

end module transition_points_mod
