module partition_function_mod
	implicit none
	private
    public  :: calculate_partition_function

contains
	!**************************************************************************************************************************
	! calculate the partition function for target or all exciton types of a given carbon nanotube
	!**************************************************************************************************************************

	subroutine calculate_partition_function(currcnt, target_exciton, partition_function)
		use cnt_class, only: cnt, exciton
		use constants_mod, only : kb
		use sim_properties_mod, only : temperature

		type(cnt), intent(in) :: currcnt
		type(exciton), optional, intent(in) :: target_exciton
		real*8, intent(out) :: partition_function
		integer :: ix, iKcm
		real*8 :: min_energy, deltaE
		integer :: alpha, ex_type

		deltaE = (-1.d0) * log(1.d-3) * kb * Temperature
		min_energy = 1.d3

		partition_function = 0.d0

		if (present(target_exciton)) then
			min_energy = minval(target_exciton%ex)
			do ix = 1,target_exciton%nx
				do iKcm = currcnt%iKcm_min_fine,currcnt%iKcm_max_fine
					if((target_exciton%ex(ix,iKcm)-min_energy) .le. deltaE) then
						partition_function = partition_function + exp(-target_exciton%ex(ix,iKcm)/(kb*Temperature))
					endif
				end do
			end do
		else

			! calculate minimum energy of all exciton types.
			do alpha = 0, 1
				do ex_type = 1, 4
					min_energy = minval((/min_energy, minval(currcnt%excitons(ex_type, alpha)%ex)/))
				enddo
			enddo

			! calculate partition function for energetically relevant points
			do alpha = 0, 1
				do ex_type = 1, 4
					do ix = 1,currcnt%excitons(ex_type, alpha)%nx
						do iKcm = currcnt%iKcm_min_fine,currcnt%iKcm_max_fine
							if((currcnt%excitons(ex_type, alpha)%ex(ix,iKcm)-min_energy) .le. deltaE) then
								partition_function = partition_function + exp(-(currcnt%excitons(ex_type, alpha)%ex(ix,iKcm))/(kb*temperature))
							endif
						enddo
					enddo
				enddo
			enddo

		endif

	end subroutine calculate_partition_function

end module partition_function_mod
