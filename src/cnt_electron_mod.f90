module cnt_electron_mod
	implicit none
	private
	public :: cnt_electron_band_structure

contains

	!**************************************************************************************************************************
	! calculate band structure of the CNT
	!**************************************************************************************************************************

	subroutine cnt_electron_band_structure(currcnt)
		use cnt_class, only: cnt
		use constants_mod
		use graphene_mod, only: graphene_electron
		use math_functions_mod, only: my_norm2
		use write_log_mod, only: write_log, log_input

		type(cnt), intent(inout) :: currcnt

		integer :: nkc, imin_sub
		integer :: i,j,mu,ik,tmpi
		integer, dimension(:), allocatable :: min_loc
		real*8 :: tmpr
		real*8, dimension(2) :: k, E1_tmp, E2_tmp, e_tmp
		real*8, dimension(:), allocatable :: k_vec,min_energy
		real*8, dimension(:,:,:), allocatable :: E_k
		complex*16, dimension(:,:,:), allocatable :: Cc_k,Cv_k
		complex*16, dimension(2) :: Cc_tmp, Cv_tmp


		! calculate CNT energy dispersion.
		nkc=2*currcnt%ikc_max+1

		allocate(k_vec(currcnt%ikc_min:currcnt%ikc_max))
		allocate(E_k(1-currcnt%Nu/2:currcnt%Nu/2,currcnt%ikc_min:currcnt%ikc_max,2))
		allocate(Cc_k(1-currcnt%Nu/2:currcnt%Nu/2,currcnt%ikc_min:currcnt%ikc_max,2))
		allocate(Cv_k(1-currcnt%Nu/2:currcnt%Nu/2,currcnt%ikc_min:currcnt%ikc_max,2))
		allocate(min_loc(0:currcnt%Nu/2))

		do ik=currcnt%ikc_min,currcnt%ikc_max
			k_vec(ik)=dble(ik)*currcnt%dk
		end do

		do mu=1-currcnt%Nu/2,currcnt%Nu/2
			do ik=currcnt%ikc_min,currcnt%ikc_max
				k = dble(mu) * currcnt%K1 + dble(ik) * currcnt%dk * currcnt%K2
				call graphene_electron(e_tmp,Cc_tmp,Cv_tmp,k,currcnt%a1,currcnt%a2)
				E_k(mu,ik,:) = e_tmp
				Cc_k(mu,ik,:) = Cc_tmp
				Cv_k(mu,ik,:) = Cv_tmp
			enddo
		enddo

		! find the subbands with a minimum energy.
		min_loc=minloc(E_k(0:currcnt%Nu/2,:,1),2)
		imin_sub=count((min_loc .lt. nkc) .and. (min_loc .gt. 1))
		allocate(currcnt%min_sub(imin_sub))
		allocate(min_energy(imin_sub))

		! store the value of mu for subbands with minimums in the variable min_sub
		i=1
		do mu=0,currcnt%Nu/2
			if ((min_loc(mu) .gt. 1) .and. (min_loc(mu) .lt. nkc)) then
			currcnt%min_sub(i)=mu
			min_energy(i)=minval(E_k(mu,:,1))
			i=i+1
			end if
		end do

		! sort the subbands
		do i=imin_sub,2,-1
			do j=i-1,1,-1
			if (min_energy(i) .lt. min_energy(j)) then
				tmpr=min_energy(i)
				tmpi=currcnt%min_sub(i)
				min_energy(i)=min_energy(j)
				currcnt%min_sub(i)=currcnt%min_sub(j)
				min_energy(j)=tmpr
				currcnt%min_sub(j)=tmpi
			end if
			end do
		end do
		! find the max k-index that energy is below threshold energy (E_th).
		ik=0
		E1_tmp=(/ min_energy(currcnt%i_sub),0.d0 /)
		E2_tmp=(/ min_energy(currcnt%i_sub),0.d0 /)
		do while ((min(E1_tmp(1),E2_tmp(1))-min_energy(currcnt%i_sub)) .le. currcnt%E_th )
			k=dble(currcnt%min_sub(currcnt%i_sub))*currcnt%K1+dble(ik)*currcnt%dk*currcnt%K2
			call graphene_electron(E1_tmp,Cc_tmp,Cv_tmp,k,currcnt%a1,currcnt%a2)
			k=dble(currcnt%min_sub(currcnt%i_sub))*currcnt%K1-dble(ik)*currcnt%dk*currcnt%K2
			call graphene_electron(E2_tmp,Cc_tmp,Cv_tmp,k,currcnt%a1,currcnt%a2)
			ik=ik+1
		end do

		! set the index boundaries for some arrays and kernels.
		currcnt%ik_max=ik											!the higher limit of k-vector that is below E_th
		currcnt%ik_min=-ik											!the lower limit of k-vector that is below E_th
		currcnt%iKcm_max=floor(currcnt%Kcm_max/currcnt%dk)			!the higher limit of center of mass wave vector that we calculate
		currcnt%iKcm_min = - currcnt%iKcm_max						!the lower limit of center of mass wave vector that we calculate
		currcnt%ikr_high=currcnt%iKcm_max-currcnt%ik_min			!the maximum index that the relative wavenumber in the entire simulation.
		currcnt%ikr_low=-currcnt%ikr_high							!the minimum index that the relative wavenumber in the entire simulation.
		currcnt%ik_high=currcnt%ikr_high+currcnt%iKcm_max			!the maximum index that the wavenumber in the entire simulation.
		currcnt%ik_low=-currcnt%ik_high								!the minimum index that the wavenumber in the entire simulation.
		currcnt%iq_max=max(2*currcnt%ikr_high,currcnt%ikc_max)		!the higher limit of the index in v_FT and esp_q
		currcnt%iq_min=-currcnt%iq_max								!the lower limit of the index in v_FT and esp_q

		currcnt%iKcm_max_fine = currcnt%iKcm_max * currcnt%dk_dkx_ratio !the upper limit of center of mass wave vector that we calculate when using a finer mesh size for exciton center of mass momentum
		currcnt%iKcm_min_fine = currcnt%iKcm_min * currcnt%dk_dkx_ratio !the lower limit of center of mass wave vector that we calculate when using a finer mesh size for exciton center of mass momentum

		! calculate the tight-binding energies and coefficients.
		allocate(currcnt%Ek(1-currcnt%Nu/2:currcnt%Nu/2, currcnt%ik_low*currcnt%dk_dkx_ratio:currcnt%ik_high*currcnt%dk_dkx_ratio, 2))
		allocate(currcnt%Cc(1-currcnt%Nu/2:currcnt%Nu/2, currcnt%ik_low*currcnt%dk_dkx_ratio:currcnt%ik_high*currcnt%dk_dkx_ratio, 2))
		allocate(currcnt%Cv(1-currcnt%Nu/2:currcnt%Nu/2, currcnt%ik_low*currcnt%dk_dkx_ratio:currcnt%ik_high*currcnt%dk_dkx_ratio, 2))

		do mu=1-currcnt%Nu/2,currcnt%Nu/2
			do ik=currcnt%ik_low*currcnt%dk_dkx_ratio,currcnt%ik_high*currcnt%dk_dkx_ratio
				k=dble(mu)*currcnt%K1+dble(ik)*currcnt%dkx*currcnt%K2
				call graphene_electron(e_tmp,Cc_tmp,Cv_tmp,k,currcnt%a1,currcnt%a2)
				currcnt%Ek(mu,ik,:) = e_tmp
				currcnt%Cc(mu,ik,:) = Cc_tmp
				currcnt%Cv(mu,ik,:) = Cv_tmp
			enddo
		enddo

		! save the index boundaries and index of minimum subband to the log file. ************************************************************************
		call write_log(new_line('A')//"Index boundaries *************************************")
		write(log_input,'(A, I0)') "ikc_max = ",currcnt%ikc_max
		call write_log(trim(log_input))
		write(log_input,'(A, I0)') "ik_max = ",currcnt%ik_max
		call write_log(trim(log_input))
		write(log_input,'(A, I0)') "iKcm_max = ",currcnt%iKcm_max
		call write_log(trim(log_input))
		write(log_input,'(A, I0)') "ikr_high = ",currcnt%ikr_high
		call write_log(trim(log_input))
		write(log_input,'(A, I0)') "ik_high = ",currcnt%ik_high
		call write_log(trim(log_input))
		write(log_input,'(A, I0)') "iq_max = ",currcnt%iq_max
		call write_log(trim(log_input))
		write(log_input,'(A, I0)') "min_sub(i_sub) = ",currcnt%min_sub(currcnt%i_sub)
		call write_log(trim(log_input))

	end subroutine cnt_electron_band_structure

end module cnt_electron_mod
