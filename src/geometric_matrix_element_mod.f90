module geometric_matrix_element_mod
	implicit none
	private
    public  :: calculate_finite_geometric_matrix_element, calculate_infinite_geometric_matrix_element, calculate_infinite_parallel_geometric_matrix_element

contains

	!***************************************************************************
	!	- calculate the geometric part of matrix element for two finite tubes
	!	  forming angle theta
	!***************************************************************************

	subroutine calculate_finite_geometric_matrix_element(iKcm1, iKcm2, cnt_1, cnt_2, geometric_matrix_element)
		use cnt_class, only: cnt
		use constants_mod, only: i1

		integer, intent(in) :: iKcm1, iKcm2
		type(cnt), intent(in) :: cnt_1, cnt_2
		complex*16, intent(out) :: geometric_matrix_element

		real*8 :: K1, K2
		integer :: iu

		K1 = dble(iKcm1)*cnt_1%dkx
		K2 = dble(iKcm2)*cnt_2%dkx

		geometric_matrix_element = (0.d0, 0.d0)
		do iu = lbound(cnt_1%r_posA3,1), ubound(cnt_1%r_posA3,1)
			geometric_matrix_element = geometric_matrix_element + sum(exp(-i1*dcmplx(2.d0*(K1*cnt_1%ur_posA3(iu,1)+dble(cnt_1%selected_exciton%mu_cm)*cnt_1%az_angle(iu))))*exp(i1*dcmplx(2.d0*(K2*cnt_2%ur_posA3(:,1)+dble(cnt_2%selected_exciton%mu_cm)*cnt_2%az_angle(:))))/dcmplx(sqrt((cnt_2%r_posA3(:,1)-cnt_1%r_posA3(iu,1))**2+(cnt_2%r_posA3(:,2)-cnt_1%r_posA3(iu,2))**2+(cnt_2%r_posA3(:,3)-cnt_1%r_posA3(iu,3))**2)))
		end do

	end subroutine calculate_finite_geometric_matrix_element


	!***************************************************************************
	!	- calculate the geometric part of matrix element for two infinite tubes
	!	  that are parallel to each other
	!***************************************************************************

	! subroutine calculate_infinite_parallel_geometric_matrix_element(iKcm1, iKcm2, cnt_1, cnt_2, c2c_distance, geometric_matrix_element)
	! 	use cnt_class, only: cnt
	! 	use constants_mod, only: i1, pi, A_u, a_l
	! 	use math_functions_mod, only: bessk0

	! 	real*8, intent(in) :: c2c_distance
	! 	integer, intent(in) :: iKcm1, iKcm2
	! 	type(cnt), intent(in) :: cnt_1, cnt_2
	! 	complex*16, intent(out) :: geometric_matrix_element

	! 	real*8 :: K1, K2
	! 	integer :: iPhi1, iPhi2
	! 	integer :: nPhi1, nPhi2
	! 	real*8 :: dPhi1, dPhi2
	! 	real*8, dimension(:), allocatable :: phi1, phi2, sinPhi1, sinPhi2, cosPhi1, cosPhi2

	! 	real*8 :: radius1, radius2

	! 	geometric_matrix_element = (0.d0, 0.d0)

	! 	if ((iKcm1 .ne. 0) .and. (iKcm2 .ne. 0)) then

	! 		radius1 = cnt_1%radius
	! 		radius2	= cnt_2%radius

	! 		nPhi1 = nint(2.d0*radius1*pi/a_l)
	! 		dPhi1 = 2.d0*pi/nPhi1

	! 		nPhi2 = nint(2.d0*radius2*pi/a_l)
	! 		dPhi2 = 2.d0*pi/nPhi2

	! 		allocate(phi1(nPhi1))
	! 		allocate(sinPhi1(nPhi1))
	! 		allocate(cosPhi1(nPhi1))
	! 		do iPhi1 = 1, nPhi1
	! 			phi1(iPhi1) = dble(iPhi1)*dPhi1
	! 			sinPhi1 = sin(phi1(iPhi1))
	! 			cosPhi1 = cos(phi1(iPhi1))
	! 		end do

	! 		allocate(phi2(nPhi2))
	! 		allocate(sinPhi2(nPhi2))
	! 		allocate(cosPhi2(nPhi2))
	! 		do iPhi2 = 1, nPhi2
	! 			phi2(iPhi2) = dble(iPhi2)*dPhi2
	! 			sinPhi2 = sin(phi2(iPhi2))
	! 			cosPhi2 = cos(phi2(iPhi2))
	! 		end do

	! 		K1 = dble(iKcm1)*cnt_1%dkx
	! 		K2 = dble(iKcm2)*cnt_2%dkx

	! 		do iPhi1 = 1,nPhi1
	! 			do iPhi2 = 1,nPhi2
	! 				geometric_matrix_element = geometric_matrix_element + (exp(i1*(2.0*(cnt_2%selected_exciton%mu_cm*phi2(iPhi2)-cnt_1%selected_exciton%mu_cm*phi1(iPhi1))))*dcmplx(bessk0(2.0*abs(K1)*sqrt((radius1*sinPhi1(iPhi1)-radius2*sinPhi2(iPhi2))**2+(c2c_distance+radius1*cosPhi1(iPhi1)-radius2*cosPhi2(iPhi2))**2))))
	! 			enddo
	! 		end do

	! 		geometric_matrix_element = geometric_matrix_element * dcmplx(dPhi1*dPhi2*2.d0*radius1*radius2/(A_u**2))

	! 	end if

	! end subroutine calculate_infinite_parallel_geometric_matrix_element

	subroutine calculate_infinite_parallel_geometric_matrix_element(geometric_matrix_element, cnt_1, cnt_2, mu_cm1, iKcm1, mu_cm2, iKcm2, c2c_distance)
		use cnt_class, only: cnt
		use constants_mod, only: i1, pi, A_u, a_l
		use math_functions_mod, only: bessk0

		type(cnt), intent(in) :: cnt_1, cnt_2
		real*8, intent(in) :: c2c_distance
		integer, intent(in) :: iKcm1, iKcm2
		integer, intent(in) :: mu_cm1, mu_cm2
		complex*16, intent(out) :: geometric_matrix_element

		real*8 :: K1
		integer :: iPhi1, iPhi2
		integer :: nPhi1, nPhi2
		real*8 :: dPhi1, dPhi2
		real*8, dimension(:), allocatable :: phi1, phi2, sinPhi1, sinPhi2, cosPhi1, cosPhi2

		real*8 :: radius1, radius2

		geometric_matrix_element = (0.d0, 0.d0)

		if (iKcm1 .ne. iKcm2) then 
			write(*,'(A)') "iKcm1 is not equal to iKcm2!!!"
			call exit()
		end if

		if (iKcm1 .ne. 0) then

			radius1 = cnt_1%radius
			radius2	= cnt_2%radius

			nPhi1 = nint(2.d0*radius1*pi/a_l)
			dPhi1 = 2.d0*pi/nPhi1

			nPhi2 = nint(2.d0*radius2*pi/a_l)
			dPhi2 = 2.d0*pi/nPhi2

			allocate(phi1(nPhi1))
			allocate(sinPhi1(nPhi1))
			allocate(cosPhi1(nPhi1))
			do iPhi1 = 1, nPhi1
				phi1(iPhi1) = dble(iPhi1)*dPhi1
				sinPhi1 = sin(phi1(iPhi1))
				cosPhi1 = cos(phi1(iPhi1))
			end do

			allocate(phi2(nPhi2))
			allocate(sinPhi2(nPhi2))
			allocate(cosPhi2(nPhi2))
			do iPhi2 = 1, nPhi2
				phi2(iPhi2) = dble(iPhi2)*dPhi2
				sinPhi2 = sin(phi2(iPhi2))
				cosPhi2 = cos(phi2(iPhi2))
			end do

			K1 = dble(iKcm1)*cnt_1%dkx

			do iPhi1 = 1,nPhi1
				do iPhi2 = 1,nPhi2
					geometric_matrix_element = geometric_matrix_element + (exp(i1*(2.0*(mu_cm2*phi2(iPhi2)-mu_cm1*phi1(iPhi1))))*dcmplx(bessk0(2.0*abs(K1)*sqrt((radius1*sinPhi1(iPhi1)-radius2*sinPhi2(iPhi2))**2+(c2c_distance+radius1*cosPhi1(iPhi1)-radius2*cosPhi2(iPhi2))**2))))
				enddo
			end do

			geometric_matrix_element = geometric_matrix_element * dcmplx(dPhi1*dPhi2*2.d0*radius1*radius2/(A_u**2))

		end if

	end subroutine calculate_infinite_parallel_geometric_matrix_element


	!***************************************************************************
	!	- calculate the geometric part of matrix element for two infinite tubes
	!	  that are unparallel and form an angle theta
	!***************************************************************************

	subroutine calculate_infinite_geometric_matrix_element(geometric_matrix_element, cnt_1, cnt_2, mu_cm1, iKcm1, mu_cm2, iKcm2, theta, c2c_distance)
		use cnt_class, only: cnt, exciton
		use constants_mod, only: i1, pi, A_u, a_l

		complex*16, intent(out) :: geometric_matrix_element
		type(cnt), intent(in) :: cnt_1, cnt_2
		integer, intent(in) :: mu_cm1, mu_cm2
		real*8, intent(in) :: theta
		real*8, intent(in) :: c2c_distance
		integer, intent(in) :: iKcm1, iKcm2

		real*8 :: K1, K2
		integer :: iPhi1, iPhi2
		integer :: nPhi1, nPhi2
		real*8 :: dPhi1, dPhi2
		real*8, dimension(:), allocatable :: phi1, phi2

		real*8 :: arg1, arg2, arg3

		real*8 :: radius1, radius2

		geometric_matrix_element = (0.d0, 0.d0)

		if ((iKcm1 .ne. 0) .or. (iKcm2 .ne. 0)) then

			radius1 = cnt_1%radius
			radius2	= cnt_2%radius

			nPhi1 = nint(2.d0*radius1*pi/a_l)
			dPhi1 = 2.d0*pi/nPhi1

			nPhi2 = nint(2.d0*radius2*pi/a_l)
			dPhi2 = 2.d0*pi/nPhi2

			allocate(phi1(nPhi1))
			do iPhi1 = 1, nPhi1
				phi1(iPhi1) = dble(iPhi1)*dPhi1
			end do

			allocate(phi2(nPhi2))
			do iPhi2 = 1, nPhi2
				phi2(iPhi2) = dble(iPhi2)*dPhi2
			end do

			K1 = dble(iKcm1)*cnt_1%dkx
			K2 = dble(iKcm2)*cnt_2%dkx

			arg1 = sqrt(K1**2+K2**2-2.d0*K1*K2*cos(theta))
			do iPhi1 = 1,nPhi1
				do iPhi2 = 1,nPhi2
					arg2 = 2.d0 * (K1*(radius2*cos(phi2(iPhi2))-radius1*cos(phi1(iPhi1))*cos(theta)) + K2*(radius1*cos(phi1(iPhi1))-radius2*cos(phi2(iPhi2))*cos(theta))) / (sin(theta))
					arg3 = 2.d0 * abs((c2c_distance+radius2*sin(phi2(iPhi2))-radius1*sin(phi1(iPhi1)))/(sin(theta)))
					geometric_matrix_element = geometric_matrix_element + exp(i1*dcmplx(-2.d0*dble(mu_cm1)*phi1(iPhi1)+2.d0*dble(mu_cm2)*phi2(iPhi2))) * exp(i1*dcmplx(arg2)) * dcmplx(exp(-arg3 * arg1))
				end do
			end do
			geometric_matrix_element = geometric_matrix_element * dcmplx(1/arg1)

			geometric_matrix_element = geometric_matrix_element * dcmplx(dPhi1*dPhi2*pi*radius1*radius2/(A_u**2))

		else
			! write(*,*) "iKcm1  and iKcm2 are ZERO!!!!"

			radius1 = cnt_1%radius
			radius2	= cnt_2%radius

			nPhi1 = nint(2.d0*radius1*pi/a_l)
			dPhi1 = 2.d0*pi/nPhi1

			nPhi2 = nint(2.d0*radius2*pi/a_l)
			dPhi2 = 2.d0*pi/nPhi2

			allocate(phi1(nPhi1))
			do iPhi1 = 1, nPhi1
				phi1(iPhi1) = dble(iPhi1)*dPhi1
			end do

			allocate(phi2(nPhi2))
			do iPhi2 = 1, nPhi2
				phi2(iPhi2) = dble(iPhi2)*dPhi2
			end do

			K1 = dble(iKcm1+1)*cnt_1%dkx
			K2 = dble(iKcm2+1)*cnt_2%dkx

			arg1 = sqrt(K1**2+K2**2-2.d0*K1*K2*cos(theta))
			do iPhi1 = 1,nPhi1
				do iPhi2 = 1,nPhi2
					arg2 = 2.d0 * (K1*(radius2*cos(phi2(iPhi2))-radius1*cos(phi1(iPhi1))*cos(theta)) + K2*(radius1*cos(phi1(iPhi1))-radius2*cos(phi2(iPhi2))*cos(theta))) / (sin(theta))
					arg3 = 2.d0 * abs((c2c_distance+radius2*sin(phi2(iPhi2))-radius1*sin(phi1(iPhi1)))/(sin(theta)))
					geometric_matrix_element = geometric_matrix_element + exp(i1*dcmplx(-2.d0*dble(mu_cm1)*phi1(iPhi1)+2.d0*dble(mu_cm2)*phi2(iPhi2))) * exp(i1*dcmplx(arg2)) * dcmplx(exp(-arg3 * arg1))
				end do
			end do
			geometric_matrix_element = geometric_matrix_element * dcmplx(1/arg1)

			geometric_matrix_element = geometric_matrix_element * dcmplx(dPhi1*dPhi2*pi*radius1*radius2/(A_u**2))

		end if

	end subroutine calculate_infinite_geometric_matrix_element


end module geometric_matrix_element_mod
