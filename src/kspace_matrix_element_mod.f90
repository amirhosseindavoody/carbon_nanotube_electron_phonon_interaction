module kspace_matrix_element_mod
	implicit none
	private
    public  :: calculate_kspace_matrix_element, calculate_Q_tilde

contains

	!***************************************************************************
	! calculate the k-space part of matrix element for the transition points
	!***************************************************************************

	subroutine calculate_kspace_matrix_element(transition_points, kspace_matrix_element, cnt_1, cnt_2)
		use cnt_class, only: cnt
		use constants_mod, only: pi, eps0, q0, i1
		use write_log_mod, only: write_log, log_input

		integer, dimension(:,:), intent(in) :: transition_points
		complex*16, dimension(:), allocatable, intent(inout) :: kspace_matrix_element
		type(cnt), intent(in) :: cnt_1, cnt_2
		complex*16 :: tmpc
		integer :: n_transition_points
		integer :: ix1,ix2
		integer :: iKcm1, iKcm2
		integer :: ikr1, ikr2
		integer :: mu_c1, mu_v1, mu_c2, mu_v2
		integer :: imur1, imur2
		integer :: ikc1, ikc2, ikv1, ikv2
		integer :: is,isp
		integer :: iT
		real*8, dimension(2) :: Kcm1, Kcm2
		real*8 , dimension(2,2) :: ds1, ds2 ! this are relative displacement of carbon atoms in graphene unit cell

		n_transition_points = size(transition_points,1)

		allocate(kspace_matrix_element(n_transition_points))
		kspace_matrix_element = (0.d0,0.d0)

		ds1(1,:) = 0.d0
		ds1(2,:) = cnt_1%aCC_vec

		ds2(1,:) = 0.d0
		ds2(2,:) = cnt_2%aCC_vec

		do iT = 1,n_transition_points
			if (mod(iT,100) .eq. 0) then
				write(log_input, '(A, I0, A, I0)') "calculating k-space matrix element: ", iT, " / ", n_transition_points
				call write_log(log_input)
			end if

			ix1 = transition_points(iT,1)
			ix2 = transition_points(iT,2)
			iKcm1 = transition_points(iT,3)
			iKcm2 = transition_points(iT,4)
			kspace_matrix_element(iT) = (0.d0,0.d0)

			Kcm1 = dble(cnt_1%selected_exciton%mu_cm) * cnt_1%K1 + dble(iKcm1) * cnt_1%dkx * cnt_1%K2
			Kcm2 = dble(cnt_2%selected_exciton%mu_cm) * cnt_2%K1 + dble(iKcm2) * cnt_2%dkx * cnt_2%K2

			tmpc = (0.d0,0.d0)

			do imur1 = 1,cnt_1%selected_exciton%n_mu_r
				mu_c1 = cnt_1%selected_exciton%mu_r(imur1) + cnt_1%selected_exciton%mu_cm
				mu_v1 = cnt_1%selected_exciton%mu_r(imur1) - cnt_1%selected_exciton%mu_cm

				do ikr1 = cnt_1%selected_exciton%ikr_low, cnt_1%selected_exciton%ikr_high
					ikc1 = ikr1 * cnt_1%dk_dkx_ratio + iKcm1
					ikv1 = ikr1 * cnt_1%dk_dkx_ratio - iKcm1

					do imur2 = 1,cnt_2%selected_exciton%n_mu_r
						mu_c2 = cnt_2%selected_exciton%mu_r(imur2) + cnt_2%selected_exciton%mu_cm
						mu_v2 = cnt_2%selected_exciton%mu_r(imur2) - cnt_2%selected_exciton%mu_cm

						do ikr2 = cnt_2%ikr_low, cnt_2%ikr_high
							ikc2 = ikr2 * cnt_2%dk_dkx_ratio + iKcm2
							ikv2 = ikr2 * cnt_2%dk_dkx_ratio - iKcm2

							do is = 1,2
								do isp = 1,2
									tmpc = tmpc + conjg(cnt_1%Cc(mu_c1,ikc1,is))*cnt_1%Cv(mu_v1,ikv1,is)*cnt_2%Cc(mu_c2,ikc2,isp)*conjg(cnt_2%Cv(mu_v2,ikv2,isp))*exp(i1*dcmplx(-2.d0*dot_product(Kcm1,ds1(is,:))+2.d0*dot_product(Kcm2,ds2(isp,:))))
								end do
							end do
							kspace_matrix_element(iT) = kspace_matrix_element(iT) + tmpc*conjg(cnt_1%selected_exciton%psi(ikr1,ix1,iKcm1,imur1))*cnt_2%selected_exciton%psi(ikr2,ix2,iKcm2,imur2)
							tmpc = (0.d0,0.d0)
						end do
					enddo
				enddo
			enddo

		end do

		kspace_matrix_element = kspace_matrix_element * dcmplx(q0**2/(4.d0*pi*eps0*sqrt(2.d0*pi/cnt_1%dk * 2.d0*pi/cnt_2%dk)))

	end subroutine calculate_kspace_matrix_element

	!***************************************************************************
	! calculate the k-space part of matrix element for specific exciton types
	! and center of mass momentum and band index
	!***************************************************************************

	subroutine calculate_Q_tilde(Q_tilde, cnt_1, cnt_2, i_exciton, f_exciton, ix1, iKcm1, ix2, iKcm2)
		use cnt_class, only: cnt, exciton
		use constants_mod, only: pi, eps0, q0, i1
		! use write_log_mod, only: write_log, log_input

		complex*16, intent(out) :: Q_tilde
		type(cnt), intent(in) :: cnt_1, cnt_2
		type(exciton), intent(in) :: i_exciton, f_exciton
		integer, intent(in) :: ix1, iKcm1, ix2, iKcm2

		complex*16 :: tmpc
		integer :: ikr1, ikr2
		integer :: mu_c1, mu_v1, mu_c2, mu_v2
		integer :: imur1, imur2
		integer :: ikc1, ikc2, ikv1, ikv2
		integer :: is,isp

		real*8, dimension(2) :: Kcm1, Kcm2
		real*8 , dimension(2,2) :: ds1, ds2 ! this are relative displacement of carbon atoms in graphene unit cell

		Q_tilde = (0.d0, 0.d0)

		ds1(1,:) = 0.d0
		ds1(2,:) = cnt_1%aCC_vec

		ds2(1,:) = 0.d0
		ds2(2,:) = cnt_2%aCC_vec

		Kcm1 = dble(i_exciton%mu_cm) * cnt_1%K1 + dble(iKcm1) * i_exciton%dKcm * cnt_1%K2
		Kcm2 = dble(f_exciton%mu_cm) * cnt_2%K1 + dble(iKcm2) * f_exciton%dKcm * cnt_2%K2

		do imur1 = 1,i_exciton%n_mu_r
			mu_c1 = i_exciton%mu_r(imur1) + i_exciton%mu_cm
			mu_v1 = i_exciton%mu_r(imur1) - i_exciton%mu_cm

			do ikr1 = i_exciton%ikr_low, i_exciton%ikr_high
				ikc1 = ikr1 * i_exciton%dkr_dKcm_ratio + iKcm1
				ikv1 = ikr1 * i_exciton%dkr_dKcm_ratio - iKcm1

				do imur2 = 1,f_exciton%n_mu_r
					mu_c2 = f_exciton%mu_r(imur2) + f_exciton%mu_cm
					mu_v2 = f_exciton%mu_r(imur2) - f_exciton%mu_cm

					do ikr2 = f_exciton%ikr_low, f_exciton%ikr_high
						ikc2 = ikr2 * f_exciton%dkr_dKcm_ratio + iKcm2
						ikv2 = ikr2 * f_exciton%dkr_dKcm_ratio - iKcm2

						tmpc = (0.d0, 0.d0)
						do is = 1,2
							do isp = 1,2
								tmpc = tmpc + conjg(cnt_1%Cc(mu_c1,ikc1,is))*cnt_1%Cv(mu_v1,ikv1,is)*cnt_2%Cc(mu_c2,ikc2,isp)*conjg(cnt_2%Cv(mu_v2,ikv2,isp))*exp(i1*dcmplx(-2.d0*dot_product(Kcm1,ds1(is,:))+2.d0*dot_product(Kcm2,ds2(isp,:))))
							end do
						end do
						Q_tilde = Q_tilde + tmpc*conjg(i_exciton%psi(ikr1,ix1,iKcm1,imur1))*f_exciton%psi(ikr2,ix2,iKcm2,imur2)
					end do
				enddo
			enddo
		enddo

		Q_tilde = Q_tilde * dcmplx(q0**2/(4.d0*pi*eps0*sqrt(2.d0*pi/i_exciton%dkr * 2.d0*pi/f_exciton%dkr)))

	end subroutine calculate_Q_tilde

end module kspace_matrix_element_mod
