module cnt_scattering_exciton_phonon_mod
	implicit none
	private
	public :: cnt_exction_phonon_scattering_rate_emission

	real*8, private :: energy_mesh_min, energy_mesh_max
	real*8, private :: energy_mesh_length = 0.8d0 ! this is the energy distance between energy_mesh_max and energy_mesh_min
	integer, private :: energy_mesh_size = 200

contains

	!***************************************************************************
	! calculate the exciton-phonon scattering rates due to phonon emission
	!***************************************************************************
	subroutine cnt_exction_phonon_scattering_rate_emission(currcnt)
		use cnt_class, only: cnt
		use cnt_phonon_mod, only: cnt_phonon_dispersion
		use constants_mod
		use graphene_mod, only: graphene_electron, graphene_electron_phonon_matrix_element
		use math_functions_mod, only: find_all_roots, first_derivative
		use sim_properties_mod, only: temperature
		use write_log_mod, only: write_log, log_input

		type(cnt), intent(inout) :: currcnt

		integer :: i, j, k
		real*8, dimension(:), allocatable :: tmp_real_array_1
		real*8, dimension(:), allocatable :: tmp_real_array_2

		integer :: ix, iKcm, mu_cm
		integer :: ix_2, iKcm_2, mu_cm_2
		integer :: mu_ph, iq_ph
		integer :: ib
		integer :: ikr, iq_ph_r
		integer :: n_initial_state, n_final_state
		integer, dimension(:), allocatable :: initial_state_idx, final_state_idx
		real*8, dimension(2) :: q_ph
		real*8, dimension(2) :: kc, kv

		character(len=1000) :: filename

		real*8, dimension(:), allocatable :: energy_mesh
		real*8, dimension(:,:), allocatable :: exciton_phonon_scattering_rate
		real*8 :: Ex_derivative, omega_phonon_derivative
		real*8 :: omega_tmp
		complex*16 :: matrix_element_conduction, matrix_element_valence
		complex*16 :: matrix_element_exciton


		!***********************************************************************
		! calculate CNT phonon energy dispersion for two brillouine zones.
		call cnt_phonon_dispersion(currcnt, dq=2.d0*currcnt%dkx, iq_max=2*currcnt%iKcm_max_fine, iq_min=2*currcnt%iKcm_min_fine, mu_max=0, mu_min=0 )

		! create the electron energy mesh
		allocate(energy_mesh(energy_mesh_size))

		energy_mesh_min = minval(currcnt%Ex0_A2(:,:))
		energy_mesh_max = energy_mesh_min + energy_mesh_length*eV

		do i=1,energy_mesh_size
			energy_mesh(i) = energy_mesh_min + dble(i-1)*(energy_mesh_max-energy_mesh_min)/dble(energy_mesh_size-1)
		enddo

		allocate(exciton_phonon_scattering_rate(energy_mesh_size,6))
		exciton_phonon_scattering_rate = 0.d0

		!now find excitonic states with energy equal to each one of the energy mesh points
		if (allocated(tmp_real_array_1)) deallocate(tmp_real_array_1)
		allocate(tmp_real_array_1(currcnt%iKcm_min_fine:currcnt%iKcm_max_fine))

		if (allocated(tmp_real_array_2)) deallocate(tmp_real_array_2)
		allocate(tmp_real_array_2(2*currcnt%iKcm_min_fine:2*currcnt%iKcm_max_fine))

		if(allocated(initial_state_idx)) deallocate(initial_state_idx)
		allocate(initial_state_idx(currcnt%iKcm_max_fine-currcnt%iKcm_min_fine+1))

		if(allocated(final_state_idx)) deallocate(final_state_idx)
		allocate(final_state_idx(currcnt%iKcm_max_fine-currcnt%iKcm_min_fine+1))

		mu_cm = 0
		mu_cm_2 = 0

		write(log_input,'(A)') "Exciton-phonon scattering rates due to phonon emission is being calculated!"
		call write_log(log_input)

		do i=1,energy_mesh_size

			write(log_input,'(A, I4.4, A, I4.4)') "i_energy_mesh = ", i, " , energy_mesh_size = ", energy_mesh_size
			call write_log(log_input)

			do ix=1,currcnt%ikr_high-currcnt%ikr_low+1
				tmp_real_array_1 = currcnt%Ex0_A2(ix,:)
				call find_all_roots (tmp_real_array_1, lbound(tmp_real_array_1,dim=1), ubound(tmp_real_array_1, dim=1), energy_mesh(i), n_initial_state, initial_state_idx)
				!find all the scattering states for the roots that has been found before
				do j=1,n_initial_state
					iKcm = initial_state_idx(j)
					! Kcm = dble(mu_cm) * currcnt%K1 + dble(iK_cm) * currcnt%dkx * currcnt%K2
					do ib=1,6
						do ix_2=1,currcnt%ikr_high-currcnt%ikr_low+1
							mu_ph = 2*(mu_cm-mu_cm_2)
							tmp_real_array_1 = currcnt%Ex0_A2(ix_2,:)-currcnt%Ex0_A2(ix,iKcm)+currcnt%omega_phonon(mu_ph,iKcm+currcnt%iKcm_max_fine:iKcm+currcnt%iKcm_min_fine:-1,ib)
							call find_all_roots (tmp_real_array_1, lbound(tmp_real_array_1,dim=1), ubound(tmp_real_array_1, dim=1), 0.d0, n_final_state, final_state_idx)

							do k = 1, n_final_state
								iKcm_2 = final_state_idx(k)
								iq_ph = (iKcm - iKcm_2) ! note that this is value is actually q/2=Kcm-Kcm_2. but we use q/2 because the phonon dispersion has been calculated with these indexes and finite difference elements dq=2*dkx

								q_ph = dble(mu_ph) * currcnt%K1 + dble(iq_ph) * (2.d0*currcnt%dkx) * currcnt%K2 ! note that dq here is dq=2*dkx

								tmp_real_array_1 = currcnt%Ex0_A2(ix_2,:)
								call first_derivative(tmp_real_array_1, lbound(tmp_real_array_1,dim=1), ubound(tmp_real_array_1, dim=1), iKcm_2, currcnt%dkx, Ex_derivative)
								tmp_real_array_2 = currcnt%omega_phonon(mu_ph, :, ib)
								call first_derivative(tmp_real_array_2, lbound(tmp_real_array_2,dim=1), ubound(tmp_real_array_2, dim=1), iq_ph, 2.d0*currcnt%dkx, omega_phonon_derivative)

								omega_tmp = currcnt%omega_phonon(mu_ph, iq_ph, ib)

								matrix_element_exciton = (0.d0,0.d0)

								iq_ph_r = int(iq_ph/currcnt%dk_dkx_ratio)

								do ikr = (currcnt%ikr_low+abs(iq_ph_r)),(currcnt%ikr_high-abs(iq_ph_r))
									! this is the sum over basis functions |+mu_r,+kr,K>
									! conduction band electron scattering terms
									kc = dble(currcnt%min_sub(currcnt%i_sub)+mu_cm) * currcnt%K1 + (dble(ikr*currcnt%dk_dkx_ratio+iKcm)) * currcnt%dkx * currcnt%K2
									call graphene_electron_phonon_matrix_element(matrix_element_conduction, kc, q_ph, ib, currcnt%a1, currcnt%a2)
									matrix_element_exciton = matrix_element_exciton + currcnt%Psi0_A2(ikr,ix,iKcm)*conjg(currcnt%Psi0_A2(ikr-iq_ph_r,ix_2,iKcm_2))*matrix_element_conduction/dcmplx(2.d0)
									! valence band electron scattering terms
									kv = dble(currcnt%min_sub(currcnt%i_sub)-mu_cm) * currcnt%K1 + (dble(ikr*currcnt%dk_dkx_ratio-iKcm)) * currcnt%dkx * currcnt%K2
									call graphene_electron_phonon_matrix_element(matrix_element_valence, kv, q_ph, ib, currcnt%a1, currcnt%a2)
									matrix_element_exciton = matrix_element_exciton + currcnt%Psi0_A2(ikr,ix,iKcm)*conjg(currcnt%Psi0_A2(ikr+iq_ph_r,ix_2,iKcm_2))*matrix_element_valence/dcmplx(2.d0)

									! this is the sum over basis functions |-mu_r,-kr,K>
									! conduction band electron scattering terms
									kc = dble(-currcnt%min_sub(currcnt%i_sub)+mu_cm) * currcnt%K1 + (-dble(ikr*currcnt%dk_dkx_ratio+iKcm)) * currcnt%dkx * currcnt%K2
									call graphene_electron_phonon_matrix_element(matrix_element_conduction, kc, q_ph, ib, currcnt%a1, currcnt%a2)
									matrix_element_exciton = matrix_element_exciton + currcnt%Psi0_A2(ikr,ix,iKcm)*conjg(currcnt%Psi0_A2(ikr-iq_ph_r,ix_2,iKcm_2))*matrix_element_conduction/dcmplx(2.d0)
									! valence band electron scattering terms
									kv = dble(-currcnt%min_sub(currcnt%i_sub)-mu_cm) * currcnt%K1 + (-dble(ikr*currcnt%dk_dkx_ratio-iKcm)) * currcnt%dkx * currcnt%K2
									call graphene_electron_phonon_matrix_element(matrix_element_valence, kv, q_ph, ib, currcnt%a1, currcnt%a2)
									matrix_element_exciton = matrix_element_exciton + currcnt%Psi0_A2(ikr,ix,iKcm)*conjg(currcnt%Psi0_A2(ikr+iq_ph_r,ix_2,iKcm_2))*matrix_element_valence/dcmplx(2.d0)
								enddo

								exciton_phonon_scattering_rate(i,ib) = exciton_phonon_scattering_rate(i,ib) + (1.d0+1.d0/(exp(omega_tmp/temperature/kb)-1.d0)) * ((abs(matrix_element_exciton))**2)/(abs(Ex_derivative/2.d0-omega_phonon_derivative)*omega_tmp/hb)

							enddo

						enddo
					enddo

				enddo
			enddo
		enddo

		exciton_phonon_scattering_rate = ((g0**2)*A_u/(16.d0*pi*m_carbon*currcnt%radius))*exciton_phonon_scattering_rate

		!***********************************************************************
		!save the calculated exciton-phonon scattering rate

		write(filename,'(A)') trim(currcnt%name)//".exciton_phonon_scattering_rate_emission.dat"
		open(unit=100,file=trim(filename),status="unknown")

		do i = 1, energy_mesh_size
			write(100,'(SP, E20.8)', advance='no') energy_mesh(i)
		enddo

		write(100,*)

		do ib=1,6
			do i = 1, energy_mesh_size
				write(100,'(SP, E20.8)', advance='no') exciton_phonon_scattering_rate(i,ib)
			enddo
			write(100,*)
		enddo

		do i = 1, energy_mesh_size
			write(100,'(SP, E20.8)', advance='no') sum(exciton_phonon_scattering_rate(i,:))
		enddo

		close(100)

		write(log_input,'(A)') "Exciton-phonon scattering rates due to phonon emission calculated and saved!"
		call write_log(log_input)
	end subroutine cnt_exction_phonon_scattering_rate_emission

end module cnt_scattering_exciton_phonon_mod
