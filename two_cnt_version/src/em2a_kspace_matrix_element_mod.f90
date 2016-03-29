module em2a_kspace_matrix_element_mod
	implicit none
	private
    public  :: calculate_em2a_kSpaceMatrixElement

contains

	!**************************************************************************************************************************
	! calculate the k-space part of matrix element for the crossing point number iT
	!**************************************************************************************************************************

	subroutine calculate_em2a_kSpaceMatrixElement(nTransitionPoints, transitionPoints, kSpaceMatrixElement)
		use comparams, only: cnt1, cnt2
		use constants_mod, only: pi, eps0, q0, i1
		use write_log_mod, only: write_log, log_input

		integer, intent(in) :: nTransitionPoints
		integer, dimension(nTransitionPoints,4), intent(in) :: transitionPoints
		complex*16, dimension(:), allocatable, intent(inout) :: kSpaceMatrixElement
		complex*16 :: tmpc, tmpc1, tmpc2
		integer :: ix1,ix2
		integer :: iKcm1, iKcm2
		integer :: ikr1, ikr2
		integer :: ikc2_p, ikc2_m, ikv2_p, ikv2_m
		integer :: ikc1, ikv1
		integer :: is,isp
		integer :: iT
		real*8, dimension(2) :: Kcm1, Kcm2
		real*8 , dimension(2,2) :: ds1, ds2 ! this are relative displacement of carbon atoms in graphene unit cell

		allocate(kSpaceMatrixElement(nTransitionPoints))
		kSpaceMatrixElement = kSpaceMatrixElement * dcmplx(0.d0,0.d0)

		ds1(1,:) = 0.d0
		ds1(2,:) = cnt1%aCC_vec

		ds2(1,:) = 0.d0
		ds2(2,:) = cnt2%aCC_vec

		do iT = 1,nTransitionPoints
			if (mod(iT,100) .eq. 0) then
				write(log_input, '("Calculating k-space matrix element: iT = ", I6.6, "  nTransitionPoints = ", I6.6)') iT, nTransitionPoints
				call write_log(log_input)
			end if

			ix1 = transitionPoints(iT,1)
			ix2 = transitionPoints(iT,2)
			iKcm1 = transitionPoints(iT,3)
			iKcm2 = transitionPoints(iT,4)
			kSpaceMatrixElement(iT) = (0.d0,0.d0)

			Kcm1 = dble(cnt1%mu_cm) * cnt1%K1 + dble(iKcm1) * cnt1%dkx * cnt1%K2
			Kcm2 = dble(cnt2%mu_cm) * cnt2%K1 + dble(iKcm2) * cnt2%dkx * cnt2%K2

			tmpc = (0.d0,0.d0)
			do ikr1 = cnt1%ikr_low, cnt1%ikr_high
				ikc1 = +ikr1 * cnt1%dk_dkx_ratio + iKcm1
				ikv1 = +ikr1 * cnt1%dk_dkx_ratio - iKcm1
				do ikr2 = cnt2%ikr_low, cnt2%ikr_high
					ikc2_p = +ikr2 * cnt2%dk_dkx_ratio + iKcm2
					ikv2_p = +ikr2 * cnt2%dk_dkx_ratio - iKcm2
					ikc2_m = -ikr2 * cnt2%dk_dkx_ratio + iKcm2
					ikv2_m = -ikr2 * cnt2%dk_dkx_ratio - iKcm2
					do is = 1,2
						do isp = 1,2
							tmpc1 = conjg(cnt1%Cc(2,ikc1,is))*cnt1%Cv(1,ikv1,is)*cnt2%Cc(1,ikc2_p,isp)*conjg(cnt2%Cv(1,ikv2_p,isp))
							tmpc2 = conjg(cnt1%Cc(2,ikc1,is))*cnt1%Cv(1,ikv1,is)*cnt2%Cc(2,ikc2_m,isp)*conjg(cnt2%Cv(2,ikv2_m,isp))*dcmplx(cnt2%ex_symmetry)
							tmpc = tmpc + (tmpc1 + tmpc2)*exp(i1*dcmplx(-2.d0* dot_product(Kcm1,ds1(is,:)) + 2.d0*dot_product(Kcm2,ds2(isp,:))))
						end do
					end do
					kSpaceMatrixElement(iT) = kSpaceMatrixElement(iT) + tmpc*conjg(cnt1%Psi_t(ikr1,ix1,iKcm1))*cnt2%Psi_t(ikr2,ix2,iKcm2)/dcmplx(sqrt(2.d0),0.d0)
					tmpc = (0.d0,0.d0)
				end do
			end do
		end do

		kSpaceMatrixElement = kSpaceMatrixElement * dcmplx(q0**2/(4.d0*pi*eps0*sqrt(2.d0*pi/cnt1%dk * 2.d0*pi/cnt2%dk)))

		return
	end subroutine calculate_em2a_kSpaceMatrixElement


end module em2a_kspace_matrix_element_mod
