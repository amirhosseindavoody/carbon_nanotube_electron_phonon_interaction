module cnt_geometry_mod
	implicit none
	private
	public :: cnt_geometry

contains

	!**********************************************************************************************************************
	! This subroutines calculates the position of atoms in carbon nanotube unit cell and the reciprocal lattice and related geometrical properties such as chiral vector, translational vector.
	!**********************************************************************************************************************

	subroutine cnt_geometry(currcnt)
		use cnt_class, only: cnt
		use constants_mod, only: a_l, pi
		use math_functions_mod, only: gcd, my_norm2
		use write_log_mod, only:write_log, log_input

		type(cnt), intent(inout) :: currcnt

		integer :: dR = 1
		integer :: t1,t2
		real*8 :: cosTh, sinTh
		real*8, dimension(2,2) :: Rot
		integer :: i,j,k

		! unit vectors and reciprocal lattice vectors.
		currcnt%a1=(/dsqrt(3.d0)/2.d0*a_l, +1.d0/2.d0*a_l/)
		currcnt%a2=(/dsqrt(3.d0)/2.d0*a_l, -1.d0/2.d0*a_l/)
		currcnt%b1=(/1.d0/dsqrt(3.d0)*2.d0*pi/a_l, +1.d0*2.d0*pi/a_l/)
		currcnt%b2=(/1.d0/dsqrt(3.d0)*2.d0*pi/a_l, -1.d0*2.d0*pi/a_l/)
		currcnt%aCC_vec=1.d0/3.d0*(currcnt%a1+currcnt%a2)

		! calculate chirality and translational vectors of CNT unit cell.
		currcnt%ch_vec=dble(currcnt%n_ch)*currcnt%a1+dble(currcnt%m_ch)*currcnt%a2
		currcnt%len_ch=a_l*dsqrt(dble(currcnt%n_ch)**2+dble(currcnt%m_ch)**2+dble(currcnt%n_ch)*dble(currcnt%m_ch))
		currcnt%radius=currcnt%len_ch/2.d0/pi

		call gcd(dR,2*currcnt%n_ch+currcnt%m_ch,2*currcnt%m_ch+currcnt%n_ch)

		t1=+int(dble(2*currcnt%m_ch+currcnt%n_ch)/dble(dR))
		t2=-int(dble(2*currcnt%n_ch+currcnt%m_ch)/dble(dR))

		currcnt%t_vec=dble(t1)*currcnt%a1+ dble(t2)*currcnt%a2

		currcnt%Nu=2*(currcnt%n_ch**2+currcnt%m_ch**2+currcnt%n_ch*currcnt%m_ch)/dR


		! rotate basis vectors so that ch_vec is along x-axis.
		cosTh=currcnt%ch_vec(1)/my_norm2(currcnt%ch_vec)
		sinTh=currcnt%ch_vec(2)/my_norm2(currcnt%ch_vec)
		Rot=reshape((/ cosTh, -sinTh , sinTh, cosTh /), (/2,2/))
		currcnt%ch_vec=matmul(Rot,currcnt%ch_vec)
		currcnt%t_vec=matmul(Rot,currcnt%t_vec)
		currcnt%a1=matmul(Rot,currcnt%a1)
		currcnt%a2=matmul(Rot,currcnt%a2)
		currcnt%b1=matmul(Rot,currcnt%b1)
		currcnt%b2=matmul(Rot,currcnt%b2)
		currcnt%aCC_vec=matmul(Rot,currcnt%aCC_vec)

		! calculate reciprocal lattice of CNT.
		currcnt%dk=my_norm2(currcnt%b1)/(dble(currcnt%nkg)-1.d0)
		currcnt%dkx=currcnt%dk / currcnt%dk_dkx_ratio
		currcnt%K1=(-dble(t2)*currcnt%b1+dble(t1)*currcnt%b2)/(dble(currcnt%Nu))
		currcnt%K2=(dble(currcnt%m_ch)*currcnt%b1-dble(currcnt%n_ch)*currcnt%b2)/dble(currcnt%Nu)
		currcnt%K2=currcnt%K2/my_norm2(currcnt%K2)
		currcnt%ikc_max=floor(pi/my_norm2(currcnt%t_vec)/currcnt%dk)
		currcnt%ikc_min=-currcnt%ikc_max

		! calculate coordinates of atoms in the unwarped CNT unit cell.
		allocate(currcnt%posA(currcnt%Nu,2))
		allocate(currcnt%posB(currcnt%Nu,2))

		k=0
		do i=0,t1+currcnt%n_ch
			do j=t2,currcnt%m_ch
			if ((dble(t2)/dble(t1)*dble(i) .le. dble(j)) .and. (dble(currcnt%m_ch)/dble(currcnt%n_ch)*dble(i) .ge. dble(j)) .and. (dble(t2)/dble(t1)*dble(i-currcnt%n_ch) .gt. dble(j-currcnt%m_ch)) .and. (dble(currcnt%m_ch)/dble(currcnt%n_ch)*dble(i-t1) .lt. dble(j-t2))) then
				k=k+1
				currcnt%posA(k,1)=dble(i)*currcnt%a1(1)+dble(j)*currcnt%a2(1)
				currcnt%posA(k,2)=dble(i)*currcnt%a1(2)+dble(j)*currcnt%a2(2)
				currcnt%posB(k,1)=currcnt%posA(k,1)+currcnt%aCC_vec(1)
				currcnt%posB(k,2)=currcnt%posA(k,2)+currcnt%aCC_vec(2)

				if (currcnt%posA(k,1) .gt. currcnt%ch_vec(1)) currcnt%posA(k,1)=currcnt%posA(k,1)-currcnt%ch_vec(1)
				if (currcnt%posA(k,1) .lt. 0) currcnt%posA(k,1)=currcnt%posA(k,1)+currcnt%ch_vec(1)
				if (currcnt%posA(k,2) .gt. currcnt%t_vec(2)) currcnt%posA(k,2)=currcnt%posA(k,2)-currcnt%t_vec(2)
				if (currcnt%posA(k,2) .lt. 0) currcnt%posA(k,2)=currcnt%posA(k,2)+currcnt%t_vec(2)

				if (currcnt%posB(k,1) .gt. currcnt%ch_vec(1)) currcnt%posB(k,1)=currcnt%posB(k,1)-currcnt%ch_vec(1)
				if (currcnt%posB(k,1) .lt. 0) currcnt%posB(k,1)=currcnt%posB(k,1)+currcnt%ch_vec(1)
				if (currcnt%posB(k,2) .gt. currcnt%t_vec(2)) currcnt%posB(k,2)=currcnt%posB(k,2)-currcnt%t_vec(2)
				if (currcnt%posB(k,2) .lt. 0) currcnt%posB(k,2)=currcnt%posB(k,2)+currcnt%t_vec(2)

			endif
			enddo
		enddo

		allocate(currcnt%pos2d(2,currcnt%Nu,2))
		currcnt%pos2d(1,:,:) = currcnt%posA(:,:)
		currcnt%pos2d(2,:,:) = currcnt%posB(:,:)

		if (k .ne. currcnt%Nu) then
			call write_log("*** Error in calculating atom positions ***")
			call exit()
		endif

		! calculate distances between atoms in a warped CNT unit cell.
		allocate(currcnt%posAA(currcnt%Nu,2))
		allocate(currcnt%posAB(currcnt%Nu,2))
		allocate(currcnt%posBA(currcnt%Nu,2))
		allocate(currcnt%posBB(currcnt%Nu,2))

		do i=1,currcnt%Nu
			currcnt%posAA(i,:)=currcnt%posA(i,:)-currcnt%posA(1,:)
			currcnt%posAB(i,:)=currcnt%posA(i,:)-currcnt%posB(1,:)
			currcnt%posBA(i,:)=currcnt%posB(i,:)-currcnt%posA(1,:)
			currcnt%posBB(i,:)=currcnt%posB(i,:)-currcnt%posB(1,:)
			if (currcnt%posAA(i,1) .gt. currcnt%ch_vec(1)/2.d0) currcnt%posAA(i,1)=currcnt%posAA(i,1)-currcnt%ch_vec(1)
			if (currcnt%posAB(i,1) .gt. currcnt%ch_vec(1)/2.d0) currcnt%posAB(i,1)=currcnt%posAB(i,1)-currcnt%ch_vec(1)
			if (currcnt%posBA(i,1) .gt. currcnt%ch_vec(1)/2.d0) currcnt%posBA(i,1)=currcnt%posBA(i,1)-currcnt%ch_vec(1)
			if (currcnt%posBB(i,1) .gt. currcnt%ch_vec(1)/2.d0) currcnt%posBB(i,1)=currcnt%posBB(i,1)-currcnt%ch_vec(1)
		end do

		! calculate coordinates of atoms in 3D unit cell
		allocate(currcnt%posA3(currcnt%Nu,3))
		allocate(currcnt%posB3(currcnt%Nu,3))

		do i=1,currcnt%Nu
			currcnt%posA3(i,1) = currcnt%radius*sin(2*pi*currcnt%posA(i,1)/currcnt%len_ch)
			currcnt%posA3(i,2) = currcnt%posA(i,2)
			currcnt%posA3(i,3) = -currcnt%radius*cos(2*pi*currcnt%posA(i,1)/currcnt%len_ch)
		end do

		do i=1,currcnt%Nu
			currcnt%posB3(i,1) = currcnt%radius*sin(2*pi*currcnt%posB(i,1)/currcnt%len_ch)
			currcnt%posB3(i,2) = currcnt%posB(i,2)
			currcnt%posB3(i,3) = -currcnt%radius*cos(2*pi*currcnt%posB(i,1)/currcnt%len_ch)
		end do

		allocate(currcnt%pos3d(2,currcnt%Nu,3))
		currcnt%pos3d(1,:,:) = currcnt%posA3(:,:)
		currcnt%pos3d(2,:,:) = currcnt%posB3(:,:)

		! write down important informations into the output file.************************************************************
		write(log_input,'(A, A)') new_line('A'), "geometrical properties **********************************";	call write_log(log_input)
		write(log_input,'(A, I0, A, I0,A)') "chirality = (",currcnt%n_ch," , ",currcnt%m_ch,")";		call write_log(trim(log_input))
		write(log_input,'(SP,A,ES10.3,A,ES10.3,A)') "a1 = (",currcnt%a1(1)," , ",currcnt%a1(2),")";		call write_log(trim(log_input))
		write(log_input,'(SP,A,ES10.3,A,ES10.3,A)') "a2 = (",currcnt%a2(1)," , ",currcnt%a2(2),")";		call write_log(trim(log_input))
		write(log_input,'(SP,A,ES10.3,A,ES10.3,A)') "b1 = (",currcnt%b1(1)," , ",currcnt%b1(2),")";		call write_log(trim(log_input))
		write(log_input,'(SP,A,ES10.3,A,ES10.3,A)') "b2 = (",currcnt%b2(1)," , ",currcnt%b2(2),")";		call write_log(trim(log_input))
		write(log_input,'(SP,A,ES10.3,A,ES10.3,A)') "aCC_vec = (",currcnt%aCC_vec(1)," , ",currcnt%aCC_vec(2),")";		call write_log(trim(log_input))
		write(log_input,'(SP,A,ES10.3,A,ES10.3,A)') "ch_vec = (",currcnt%ch_vec(1)," , ",currcnt%ch_vec(2),")";		call write_log(trim(log_input))
		write(log_input,'(SP,A,ES10.3,A,ES10.3,A)') "t_vec = (",currcnt%t_vec(1)," , ",currcnt%t_vec(2),")";		call write_log(trim(log_input))
		write(log_input,'(SP,A,ES10.3)') "len_ch = ",currcnt%len_ch;		call write_log(trim(log_input))
		write(log_input,'(SP,A,I0)') "Nu = ",currcnt%Nu;		call write_log(trim(log_input))
		write(log_input,'(SP,A,ES10.3)') "dk = ",currcnt%dk;		call write_log(trim(log_input))
		write(log_input,'(SP,A,ES10.3)') "dkx = ",currcnt%dkx;		call write_log(trim(log_input))

	end subroutine cnt_geometry

end module cnt_geometry_mod
