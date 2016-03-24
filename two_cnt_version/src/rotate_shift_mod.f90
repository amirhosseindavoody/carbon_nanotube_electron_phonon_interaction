module rotate_shift_mod
	implicit none
	private
    public  :: rotate_shift_cnt
    
contains
	
	!**************************************************************************************************************************
	! calculate the k-space part of matrix element for the crossing point number iC
	!**************************************************************************************************************************
	
	subroutine rotate_shift_cnt(currcnt, theta, c2cDistance)
		use cnt_class, only: cnt
		use physicalConstants, only: pi
		use math_functions_mod, only: my_norm2

		type(cnt), intent(inout) :: currcnt
		real*8, intent(in) :: theta
		real*8, intent(in) :: c2cDistance

		integer :: n_cnt_unitcell
		integer :: iT, i
		real*8 :: t_vec_length
		real*8 :: tmpr
		real*8, dimension(3) :: t_vec_3D
		real*8, dimension(2,2) :: Rot

		t_vec_3D(1:2) = currcnt%t_vec
		t_vec_3D(3) = 0.d0

		t_vec_length = my_norm2(currcnt%t_vec)
		n_cnt_unitcell = nint(currcnt%length/t_vec_length)

		if(.not. allocated(currcnt%ur_posA3)) allocate(currcnt%ur_posA3(n_cnt_unitcell*currcnt%Nu,3))
		if(.not. allocated(currcnt%r_posA3)) allocate(currcnt%r_posA3(n_cnt_unitcell*currcnt%Nu,3))
		if(.not. allocated(currcnt%az_angle)) allocate(currcnt%az_angle(n_cnt_unitcell*currcnt%Nu))

		do iT = 0, n_cnt_unitcell-1
			! when posA3 is calculated we assume the CNT is along y-axis. Here, we change to axis labeling so that the CNT is along x-axis.
			currcnt%ur_posA3(iT*currcnt%Nu+1:(iT+1)*currcnt%Nu,1) = iT*t_vec_3D(2) + currcnt%posA3(:,2)
			currcnt%ur_posA3(iT*currcnt%Nu+1:(iT+1)*currcnt%Nu,2) = iT*t_vec_3D(3) + currcnt%posA3(:,3)
			currcnt%ur_posA3(iT*currcnt%Nu+1:(iT+1)*currcnt%Nu,3) = iT*t_vec_3D(1) + currcnt%posA3(:,1)
			! calculate azimuthal angle of carbon atoms
			currcnt%az_angle(iT*currcnt%Nu+1:(iT+1)*currcnt%Nu) = 2.d0*pi*currcnt%posA(:,1)/currcnt%len_ch
		enddo

		! shift the tube along z-axis
		currcnt%ur_posA3(:,3) = currcnt%ur_posA3(:,3) + c2cDistance

		! shift the tube along the x_axis so the center is at currcnt%center_position
		tmpr = currcnt%center_position - (minval(currcnt%ur_posA3(:,1))+maxval(currcnt%ur_posA3(:,1)))/2.d0
		currcnt%ur_posA3(:,1) = currcnt%ur_posA3(:,1) + tmpr

		!rotate the tube in the x-y plane by angle theta
		Rot=reshape((/ cos(theta), -sin(theta) , sin(theta), cos(theta) /), (/2,2/))
		do i = lbound(currcnt%ur_posA3,1), ubound(currcnt%ur_posA3,1)
			currcnt%r_posA3(i,1:2) = matmul(Rot,currcnt%ur_posA3(i,1:2))
			currcnt%r_posA3(i,3) = currcnt%ur_posA3(i,3)
		enddo

			
		return		
	end subroutine rotate_shift_cnt

	

				
end module rotate_shift_mod