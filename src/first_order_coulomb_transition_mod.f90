!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate first order Coulomb mediated transition rates
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module first_order_coulomb_transition_mod
	implicit none
	private
	public :: calculate_first_order_transition_rates

	real*8, dimension(:,:), allocatable :: transition_rate

	real*8 :: dc2c

	real*8 :: dTheta

	complex*16, dimension(:), allocatable :: kspace_matrix_element_same_energy, kspace_matrix_element_crossing_points

contains
	!**************************************************************************************************************************
	! calculate transition table
	!**************************************************************************************************************************

	subroutine calculate_first_order_transition_rates (cnt_1,cnt_2)
		use constants_mod, only: pi, kb, hb, A_u
		use cnt_class, only: cnt, exciton
		use geometric_matrix_element_mod, only: calculate_finite_geometric_matrix_element, calculate_infinite_geometric_matrix_element, calculate_infinite_parallel_geometric_matrix_element
		use kspace_matrix_element_mod, only: calculate_kspace_matrix_element
		use math_functions_mod, only: first_derivative
		use partition_function_mod, only: calculate_partition_function
		use rotate_shift_mod, only: rotate_shift_cnt
		use sim_properties_mod, only: ppLen, temperature, c2c_min, c2c_max, n_c2c, theta_min, theta_max, n_theta
		use transition_points_mod, only: find_same_energy, same_energy, find_crossings, crossing_points
		use write_log_mod, only: write_log, log_input

		type(cnt), intent(inout) :: cnt_1,cnt_2

		integer :: n_same_energy, iS
		integer :: n_crossing, iC
		integer :: ix1, ix2, iKcm1, iKcm2
		real*8 :: partition_function
		real*8 :: dos
		complex*16 :: matrix_element, geometric_matrix_element
		real*8, dimension(:), allocatable :: tmp_array

		integer :: ic2c
		real*8 :: c2c_distance
		integer :: i_theta
		real*8 :: theta

		call write_log(new_line('A')//"************** Start calculating transition table ****************")

		! set seperation properties
		if (n_c2c .ne. 1) then
			dc2c = (c2c_max-c2c_min)/dble(n_c2c-1)
		else
			dc2c = 0.d0
		end if

		write(log_input, '("c2c_min[nm] = ", F0.1)') c2c_min*1.d9
		call write_log(trim(log_input))
		write(log_input, '("c2c_max[nm] = ", F0.1)') c2c_max*1.d9
		call write_log(trim(log_input))
		write(log_input, '("n_c2c = ", I0)') n_c2c
		call write_log(trim(log_input))

		! set orientation properties
		if (n_theta .ne. 1) then
			dTheta = (theta_max-theta_min)/dble(n_theta-1)
		else
			dTheta = 0.d0
		end if

		write(log_input, '("theta_min = ", I0)') nint(theta_min*180/pi)
		call write_log(trim(log_input))
		write(log_input, '("theta_max = ", I0)') nint(theta_max*180/pi)
		call write_log(trim(log_input))
		write(log_input, '("n_theta = ", I0)') n_theta
		call write_log(trim(log_input))

		!allocate the transition rate table
		allocate(transition_rate(n_theta,n_c2c))
		transition_rate = 0.d0

		!calculate the crossing points and points with the same energy between cnt_1 and cnt_2
		call find_same_energy(cnt_1%selected_exciton, cnt_2%selected_exciton)
		call find_crossings(cnt_1%selected_exciton, cnt_2%selected_exciton)

		!calculate the k-space part of matrix elements for the calculated transition points
		call calculate_kspace_matrix_element(same_energy, kspace_matrix_element_same_energy, cnt_1, cnt_2)
		call calculate_kspace_matrix_element(crossing_points, kspace_matrix_element_crossing_points, cnt_1, cnt_2)

		!calculate the partition function for the donor carbon nanotube (cnt_1)
		call calculate_partition_function(cnt_1, cnt_1%selected_exciton, partition_function)


		do ic2c = 1, n_c2c
			c2c_distance = c2c_min+dble(ic2c-1)*dc2c

			do i_theta = 1, n_theta
				theta = theta_min + dble(i_theta-1)*dTheta

				! calculate exciton transfer rate for finite length CNTs
				if ((cnt_1%length .lt. huge(1.d0)) .and. (cnt_2%length .lt. huge(1.d0))) then
					write(log_input,'(A, I0, A, I0, A, I0, A, I0)') 'calculating finite transition rate: i_theta=', i_theta, ', n_theta=', n_theta, 'ic2c=', ic2c, ', n_c2c=', n_c2c
					call write_log(log_input)

					n_same_energy = size(same_energy,1)

					call rotate_shift_cnt(cnt_1, 0.d0, 0.d0)
					call rotate_shift_cnt(cnt_2, theta, c2c_distance)

					do iS = 1,n_same_energy

						ix1 = same_energy(iS,1)
						ix2 = same_energy(iS,2)
						iKcm1 = same_energy(iS,3)
						iKcm2 = same_energy(iS,4)

						call calculate_finite_geometric_matrix_element(iKcm1, iKcm2, cnt_1, cnt_2, geometric_matrix_element)

						matrix_element = geometric_matrix_element * kspace_matrix_element_same_energy(iS)

						! calculate the density of states for exciton band.
						if (allocated(tmp_array)) deallocate(tmp_array)
						allocate(tmp_array(lbound(cnt_2%selected_exciton%ex, dim=2):ubound(cnt_2%selected_exciton%ex, dim=2)))
						tmp_array = cnt_2%selected_exciton%ex(ix2,:)
						call first_derivative(tmp_array, lbound(tmp_array,dim=1), ubound(tmp_array,dim=1), iKcm2, cnt_2%dkx, dos)
						if (dos .eq. 0.d0) then
							call first_derivative(tmp_array, lbound(tmp_array,dim=1), ubound(tmp_array,dim=1), iKcm2+1, cnt_2%dkx, dos)
						endif
						dos = abs(1/dos)

						transition_rate(i_theta,ic2c) = transition_rate(i_theta,ic2c) + exp(-(cnt_1%selected_exciton%ex(ix1,iKcm1))/kb/Temperature) * ((abs(matrix_element))**2) * dos * (A_u**2/(4.d0*pi*pi*cnt_1%radius*cnt_2%radius))**2 / hb / cnt_1%length / partition_function

					end do

				! calculate exciton transfer rate for infinitely long CNTs
				else
					write(log_input,'(A, I0, A, I0, A, I0, A, I0)') 'calculating infinite transition rate: i_theta=', i_theta, ', n_theta=', n_theta, 'ic2c=', ic2c, ', n_c2c=', n_c2c
					call write_log(log_input)

					if (theta .eq. 0.d0) then
						n_crossing = size(crossing_points,1)

						do iC = 1,n_crossing

							ix1 = crossing_points(iC,1)
							ix2 = crossing_points(iC,2)
							iKcm1 = crossing_points(iC,3)
							iKcm2 = crossing_points(iC,4)

							call calculate_infinite_parallel_geometric_matrix_element(iKcm1, iKcm2, cnt_1, cnt_2, c2c_distance, geometric_matrix_element)

							matrix_element = geometric_matrix_element * kspace_matrix_element_crossing_points(iC)

							! calculate the density of states for exciton band.
							if (allocated(tmp_array)) deallocate(tmp_array)
							allocate(tmp_array(lbound(cnt_2%selected_exciton%ex, dim=2):ubound(cnt_2%selected_exciton%ex, dim=2)))
							tmp_array = cnt_2%selected_exciton%ex(ix2,:)
							call first_derivative(tmp_array, lbound(tmp_array,dim=1), ubound(tmp_array,dim=1), iKcm2, cnt_2%dkx, dos)
							if (dos .eq. 0.d0) then
								call first_derivative(tmp_array, lbound(tmp_array,dim=1), ubound(tmp_array,dim=1), iKcm2+1, cnt_2%dkx, dos)
							endif
							dos = abs(1/dos)

							transition_rate(i_theta,ic2c) = transition_rate(i_theta,ic2c) + exp(-(cnt_1%selected_exciton%ex(ix1,iKcm1))/kb/temperature) * (abs(matrix_element)**2) * dos * (A_u**2/(4.d0*pi*pi*cnt_1%radius*cnt_2%radius))**2 * (2*pi/cnt_1%dkx) / hb / partition_function

						end do

					else
						n_same_energy = size(same_energy,1)

						do iS = 1,n_same_energy

							ix1 = same_energy(iS,1)
							ix2 = same_energy(iS,2)
							iKcm1 = same_energy(iS,3)
							iKcm2 = same_energy(iS,4)

							call calculate_infinite_geometric_matrix_element(iKcm1, iKcm2, cnt_1, cnt_2, theta, c2c_distance, geometric_matrix_element)

							matrix_element = geometric_matrix_element * kspace_matrix_element_same_energy(iS)

							! calculate the density of states for exciton band.
							if (allocated(tmp_array)) deallocate(tmp_array)
							allocate(tmp_array(lbound(cnt_2%selected_exciton%ex, dim=2):ubound(cnt_2%selected_exciton%ex, dim=2)))
							tmp_array = cnt_2%selected_exciton%ex(ix2,:)
							call first_derivative(tmp_array, lbound(tmp_array,dim=1), ubound(tmp_array,dim=1), iKcm2, cnt_2%dkx, dos)
							if (dos .eq. 0.d0) then
								call first_derivative(tmp_array, lbound(tmp_array,dim=1), ubound(tmp_array,dim=1), iKcm2+1, cnt_2%dkx, dos)
							endif
							dos = abs(1/dos)

							transition_rate(i_theta,ic2c) = transition_rate(i_theta,ic2c) + exp(-(cnt_1%selected_exciton%ex(ix1,iKcm1))/kb/temperature) * (abs(matrix_element)**2) * dos * (A_u**2/(4.d0*pi*pi*cnt_1%radius*cnt_2%radius))**2 * (sin(theta)) / hb / ppLen/ partition_function

						end do
					endif

				end if

				call save_transition_rates(transition_rate(i_theta,ic2c), theta, c2c_distance)

			end do
		end do

	end subroutine calculate_first_order_transition_rates

	!***************************************************************************
	! save the calculated transition table on the fly
	!***************************************************************************

	subroutine save_transition_rates(t_rate, theta, c2c_distance)

		real*8, intent(in) :: t_rate, theta, c2c_distance
		logical :: flgexist

		! write 1 to 2 transition rates
		inquire(file="transition_rates_coulomb.dat",exist=flgexist)
		if (flgexist) then
			open(unit=100, file="transition_rates_coulomb.dat", status="old", position="append", action="write")
		else
			open(unit=100, file="transition_rates_coulomb.dat", status="new", action="write")
		end if

		write(100,'(E16.8)', advance='no') t_rate
		close(100)

		! write theta values
		inquire(file="theta.dat",exist=flgexist)
		if (flgexist) then
			open(unit=100, file="theta.dat", status="old", position="append", action="write")
		else
			open(unit=100, file="theta.dat", status="new", action="write")
		end if

		write(100,'(E16.8)', advance='no') theta
		close(100)

		! write c2c values
		inquire(file="c2c.dat",exist=flgexist)
		if (flgexist) then
			open(unit=100, file="c2c.dat", status="old", position="append", action="write")
		else
			open(unit=100, file="c2c.dat", status="new", action="write")
		end if

		write(100,'(E16.8)', advance='no') c2c_distance
		close(100)

	end subroutine save_transition_rates

end module first_order_coulomb_transition_mod
