module cnt_phonon_mod
	implicit none
	private

	public :: cnt_phonon_dispersion

contains

	subroutine cnt_phonon_dispersion(currcnt)
		use cnt_class, only : cnt
		use write_log_mod, only: writeLog

		type(cnt) :: currcnt
		integer :: ik
		integer :: mu, ib
		real*8, dimension(2) :: k
		real*8, dimension(:), allocatable :: k_vec
		real*8, dimension(6) :: omega_tmp
		complex*16, dimension(6,6) :: u_ph

		allocate(k_vec(currcnt%ikc_min:currcnt%ikc_max))
		allocate(currcnt%omega_phonon(1-currcnt%Nu/2:currcnt%Nu/2,currcnt%ikc_min:currcnt%ikc_max,6))

		do ik=currcnt%ikc_min,currcnt%ikc_max
			k_vec(ik)=dble(ik)*currcnt%dk
		end do

		do mu=1-currcnt%Nu/2,currcnt%Nu/2
			do ik=currcnt%ikc_min,currcnt%ikc_max
				k=dble(mu)*currcnt%K1+dble(ik)*currcnt%dk*currcnt%K2
				call graphene_phonon_dispersion(omega_tmp,u_ph,k,currcnt%aCC_vec)
				currcnt%omega_phonon(mu,ik,:) = omega_tmp(:)
			enddo
		enddo

		! save the CNT phonon energy dispersion*************************************************************************************
		open(unit=100,file='phonon_dispersion.dat',status="unknown")

		do ik=currcnt%ikc_min,currcnt%ikc_max
			write(100,'(E16.8)', advance='no') k_vec(ik)
		end do

		write(100,*)
		write(101,*)

		do ib=1,6
			do mu=1-currcnt%Nu/2,currcnt%Nu/2
				do ik=currcnt%ikc_min,currcnt%ikc_max
					write(100,'(E16.8)', advance='no') currcnt%omega_phonon(mu,ik,ib)
				end do
				write(100,*)
			enddo
		enddo

		close(100)

	end subroutine cnt_phonon_dispersion

	!**********************************************************************************************************************
	! This subroutines calculates the energy (omega) and displacement vectors (u_ph) in graphene due to a phonon with a certain wave vector (k).
	!**********************************************************************************************************************

	subroutine graphene_phonon_dispersion(omega,u_ph,k,aCC_vec)
		use math_functions_mod, only: eig
		use physical_constant_mod, only: m_carbon_dispersion, spring_conv_coeff, pi, i1, eV

		real*8, dimension(2), intent(in) :: aCC_vec
		real*8, dimension(2), intent(in) :: k
		complex*16, dimension(6,6), intent(out) :: u_ph
		real*8, dimension(6), intent(out) :: omega

		integer :: i, j
		real*8 :: theta, tmpr
		real*8, dimension(2) :: deltaR
		real*8, dimension(2,2) :: Rot
		real*8, dimension(3,3) :: K1, K2, K3, K4, K_tmp
		real*8, dimension(3,3) :: U, U_inv
		real*8, dimension(3) :: anglesA1, anglesA3, anglesB1, anglesB3
		real*8, dimension(6) :: anglesA2, anglesA4, anglesB2, anglesB4
		complex*16, dimension(3,3) :: D_AA=(0.d0,0.d0), D_AB=(0.d0,0.d0), D_BA=(0.d0,0.d0), D_BB=(0.d0,0.d0)
		complex*16, dimension(6,6) :: D_tot
		complex*16, dimension(6) :: u_tmp

		D_AA=0.d0*D_AA
		D_AB=0.d0*D_AB
		D_BA=0.d0*D_BA
		D_BB=0.d0*D_BB

		! note that these matrices are different from reciprocal lattice vectors K1 and K2 defined in comparams module
		K1=spring_conv_coeff*reshape((/  3.65d1, 0.d0, 0.d0, 0.d0,  2.45d1, 0.d0, 0.d0, 0.d0,  9.82d0 /) , (/3,3/))
		K2=spring_conv_coeff*reshape((/  8.80d0, 0.d0, 0.d0, 0.d0, -3.23d0, 0.d0, 0.d0, 0.d0, -4.0d-1 /) , (/3,3/))
		K3=spring_conv_coeff*reshape((/  3.00d0, 0.d0, 0.d0, 0.d0, -5.25d0, 0.d0, 0.d0, 0.d0,  1.5d-1 /) , (/3,3/))
		K4=spring_conv_coeff*reshape((/ -1.92d0, 0.d0, 0.d0, 0.d0,  2.29d0, 0.d0, 0.d0, 0.d0, -5.8d-1 /) , (/3,3/))

		! set the angle of 1st, 2nd, 3rd, and 4th nearest neighbors of A atoms
		anglesA1=pi/3.0d0*(/ 0.0d0, 2.0d0, 4.0d0 /)
		anglesA2=pi/6.0d0*(/ 1.0d0, 3.0d0, 5.0d0, 7.0d0, 9.0d0, 1.1d1 /)
		anglesA3=pi/3.0d0*(/ 1.0d0, 3.0d0, 5.0d0 /)
		anglesA4=pi/1.2d1*(/ 1.0d0, 7.0d0, 9.0d0, 1.5d1, 1.7d1, 2.3d1 /)

		anglesB1=pi/3.0d0*(/ 1.0d0, 3.0d0, 5.0d0 /)
		anglesB2=pi/6.0d0*(/ 1.0d0, 3.0d0, 5.0d0, 7.0d0, 9.0d0, 1.1d1 /)
		anglesB3=pi/3.0d0*(/ 0.0d0, 2.0d0, 4.0d0 /)
		anglesB4=pi/1.2d1*(/ 3.0d0, 5.0d0, 1.1d1, 1.3d1, 1.9d1, 2.1d1 /)

		! calculate D_AA, D_AB, D_BA, D_BB for 1st and 3rd nearest neighbors of atoms A and B type.
		do i=1,3
			! 1st nearest neighbor for A atoms
			theta=anglesA1(i)
			Rot=reshape((/ cos(theta), sin(theta) , -sin(theta), cos(theta) /), (/2,2/))
			U=reshape((/ cos(theta), -sin(theta), 0.d0, sin(theta), cos(theta), 0.d0, 0.d0, 0.d0, 1.d0 /) , (/3,3/))
			U_inv=reshape((/ cos(theta), sin(theta), 0.d0, -sin(theta), cos(theta), 0.d0, 0.d0, 0.d0, 1.d0 /) , (/3,3/))
			K_tmp=matmul(U_inv,matmul(K1,U))
			deltaR=matmul(Rot,aCC_vec)
			D_AA=D_AA+K_tmp
			D_AB=D_AB-exp(i1*dcmplx(dot_product(k,deltaR)))*dcmplx(K_tmp)

			! 3rd nearest neighbor for A atoms
			theta=anglesA3(i)
			Rot=reshape((/ cos(theta), sin(theta) , -sin(theta), cos(theta) /), (/2,2/))
			U=reshape((/ cos(theta), -sin(theta), 0.d0, sin(theta), cos(theta), 0.d0, 0.d0, 0.d0, 1.d0 /) , (/3,3/))
			U_inv=reshape((/ cos(theta), sin(theta), 0.d0, -sin(theta), cos(theta), 0.d0, 0.d0, 0.d0, 1.d0 /) , (/3,3/))
			K_tmp=matmul(U_inv,matmul(K3,U))
			deltaR=matmul(Rot,2.d0*aCC_vec)
			D_AA=D_AA+K_tmp
			D_AB=D_AB-exp(i1*dcmplx(dot_product(k,deltaR)))*dcmplx(K_tmp)

			! 1st nearest neighbor for B atoms
			theta=anglesB1(i)
			Rot=reshape((/ cos(theta), sin(theta) , -sin(theta), cos(theta) /), (/2,2/))
			U=reshape((/ cos(theta), -sin(theta), 0.d0, sin(theta), cos(theta), 0.d0, 0.d0, 0.d0, 1.d0 /) , (/3,3/))
			U_inv=reshape((/ cos(theta), sin(theta), 0.d0, -sin(theta), cos(theta), 0.d0, 0.d0, 0.d0, 1.d0 /) , (/3,3/))
			K_tmp=matmul(U_inv,matmul(K1,U))
			deltaR=matmul(Rot,aCC_vec)
			D_BB=D_BB+K_tmp
			D_BA=D_BA-exp(i1*dcmplx(dot_product(k,deltaR)))*dcmplx(K_tmp)

			! 3rd nearest neighbor for B atoms
			theta=anglesB3(i)
			Rot=reshape((/ cos(theta), sin(theta) , -sin(theta), cos(theta) /), (/2,2/))
			U=reshape((/ cos(theta), -sin(theta), 0.d0, sin(theta), cos(theta), 0.d0, 0.d0, 0.d0, 1.d0 /) , (/3,3/))
			U_inv=reshape((/ cos(theta), sin(theta), 0.d0, -sin(theta), cos(theta), 0.d0, 0.d0, 0.d0, 1.d0 /) , (/3,3/))
			K_tmp=matmul(U_inv,matmul(K3,U))
			deltaR=matmul(Rot,2.d0*aCC_vec)
			D_BB=D_BB+K_tmp
			D_BA=D_BA-exp(i1*dcmplx(dot_product(k,deltaR)))*dcmplx(K_tmp)
		enddo

		! calculate D_AA, D_AB, D_BA, D_BB for 2nd and 4th nearest neighbors of atoms A and B type.
		do i=1,6
			! 2st nearest neighbor for A atoms
			theta=anglesA2(i)
			Rot=reshape((/ cos(theta), sin(theta) , -sin(theta), cos(theta) /), (/2,2/))
			U=reshape((/ cos(theta), -sin(theta), 0.d0, sin(theta), cos(theta), 0.d0, 0.d0, 0.d0, 1.d0 /) , (/3,3/))
			U_inv=reshape((/ cos(theta), sin(theta), 0.d0, -sin(theta), cos(theta), 0.d0, 0.d0, 0.d0, 1.d0 /) , (/3,3/))
			K_tmp=matmul(U_inv,matmul(K2,U))
			deltaR=matmul(Rot,sqrt(3.d0)*aCC_vec)
			D_AA=D_AA+(dcmplx(1.d0)-exp(i1*dcmplx(dot_product(k,deltaR))))*dcmplx(K_tmp)

			! 4rd nearest neighbor for A atoms
			theta=anglesA4(i)
			Rot=reshape((/ cos(theta), sin(theta) , -sin(theta), cos(theta) /), (/2,2/))
			U=reshape((/ cos(theta), -sin(theta), 0.d0, sin(theta), cos(theta), 0.d0, 0.d0, 0.d0, 1.d0 /) , (/3,3/))
			U_inv=reshape((/ cos(theta), sin(theta), 0.d0, -sin(theta), cos(theta), 0.d0, 0.d0, 0.d0, 1.d0 /) , (/3,3/))
			K_tmp=matmul(U_inv,matmul(K4,U))
			deltaR=matmul(Rot,sqrt(7.d0)*aCC_vec)
			D_AA=D_AA+K_tmp
			D_AB=D_AB-exp(i1*dcmplx(dot_product(k,deltaR)))*dcmplx(K_tmp)

			! 2st nearest neighbor for B atoms
			theta=anglesB2(i)
			Rot=reshape((/ cos(theta), sin(theta) , -sin(theta), cos(theta) /), (/2,2/))
			U=reshape((/ cos(theta), -sin(theta), 0.d0, sin(theta), cos(theta), 0.d0, 0.d0, 0.d0, 1.d0 /) , (/3,3/))
			U_inv=reshape((/ cos(theta), sin(theta), 0.d0, -sin(theta), cos(theta), 0.d0, 0.d0, 0.d0, 1.d0 /) , (/3,3/))
			K_tmp=matmul(U_inv,matmul(K2,U))
			deltaR=matmul(Rot,sqrt(3.d0)*aCC_vec)
			D_BB=D_BB+(dcmplx(1.d0)-exp(i1*dcmplx(dot_product(k,deltaR))))*dcmplx(K_tmp)

			! 4rd nearest neighbor for B atoms
			theta=anglesB4(i)
			Rot=reshape((/ cos(theta), sin(theta) , -sin(theta), cos(theta) /), (/2,2/))
			U=reshape((/ cos(theta), -sin(theta), 0.d0, sin(theta), cos(theta), 0.d0, 0.d0, 0.d0, 1.d0 /) , (/3,3/))
			U_inv=reshape((/ cos(theta), sin(theta), 0.d0, -sin(theta), cos(theta), 0.d0, 0.d0, 0.d0, 1.d0 /) , (/3,3/))
			K_tmp=matmul(U_inv,matmul(K4,U))
			deltaR=matmul(Rot,sqrt(7.d0)*aCC_vec)
			D_BB=D_BB+K_tmp
			D_BA=D_BA-exp(i1*dcmplx(dot_product(k,deltaR)))*dcmplx(K_tmp)
		enddo

		D_tot(1:3,1:3)=D_AA
		D_tot(1:3,4:6)=D_AB
		D_tot(4:6,1:3)=D_BA
		D_tot(4:6,4:6)=D_BB

		D_tot=dcmplx(1.d0/2.d0)*(D_tot+conjg(transpose(D_tot)))

		call eig(6,D_tot,u_ph,omega)

		do i=1,6
			omega(i)=sqrt(abs(omega(i))/m_carbon_dispersion)
		enddo

		! sort phonon energies
		do i=6,2,-1
			do j=i-1,1,-1
				if (omega(i) .lt. omega(j)) then
					tmpr=omega(i)
					u_tmp=u_ph(:,i)
					omega(i)=omega(j)
					u_ph(:,i)=u_ph(:,j)
					omega(j)=tmpr
					u_ph(:,j)=u_tmp
				end if
			end do
		end do

		! change the units from [1/cm] to [Joules]
		omega=omega*(4.1357d-3/3.3356d1)*eV

		return
	end subroutine graphene_phonon_dispersion

end module cnt_phonon_mod
