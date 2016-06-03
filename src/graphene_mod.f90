module graphene_mod
	implicit none
	private
	public :: graphene_electron, graphene_phonon, graphene_electron_phonon_matrix_element

contains

	!***************************************************************************
	! subroutine to calculate electron Bloch functions and energy in graphene
	!***************************************************************************

	subroutine graphene_electron(E,Cc,Cv,k,a1,a2)
		use constants_mod, only: i1, t0

		real*8, dimension(2), intent(in) :: k
		real*8, dimension(2), intent(in) :: a1,a2
		real*8, dimension(2), intent(out) :: E
		complex*16, dimension(2), intent(out) :: Cv
		complex*16, dimension(2), intent(out) :: Cc
		complex*16 :: f_k

		f_k=exp(i1*dcmplx(dot_product(k,(a1+a2)/3.d0)))+exp(i1*dcmplx(dot_product(k,(a1-2.d0*a2)/3.d0)))+exp(i1*dcmplx(dot_product(k,(a2-2.d0*a1)/3.d0)))

		E(1)=+t0*abs(f_k)
		E(2)=-t0*abs(f_k)

		Cc(1)=dcmplx(+1.d0/sqrt(2.d0))
		Cc(2)=dcmplx(+1.d0/sqrt(2.d0)/abs(f_k))*conjg(f_k)
		Cv(1)=dcmplx(+1.d0/sqrt(2.d0))
		Cv(2)=dcmplx(-1.d0/sqrt(2.d0)/abs(f_k))*conjg(f_k)
	end subroutine graphene_electron

	!***************************************************************************
	!	- This subroutines calculates phonon energy (omega) and displacement
	!	  vectors (u_ph) in graphene due to a phonon with a certain
	!	  wave vector (k).
	!***************************************************************************

	subroutine graphene_phonon(omega,u_ph,k,aCC_vec)
		use constants_mod, only: m_carbon_dispersion, spring_conv_coeff, pi, i1, eV
		use math_functions_mod, only: eig

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
		real*8 :: angle_tmp
		complex*16, dimension(3,3) :: D_AA, D_AB, D_BA, D_BB
		complex*16, dimension(6,6) :: D_tot
		complex*16, dimension(6) :: u_tmp

		D_AA=(0.d0,0.d0)
		D_AB=(0.d0,0.d0)
		D_BA=(0.d0,0.d0)
		D_BB=(0.d0,0.d0)

		! note that these matrices are different from reciprocal lattice vectors K1 and K2 defined in comparams module
		K1=spring_conv_coeff*reshape((/  3.65d1, 0.d0, 0.d0, 0.d0,  2.45d1, 0.d0, 0.d0, 0.d0,  9.82d0 /) , (/3,3/))
		K2=spring_conv_coeff*reshape((/  8.80d0, 0.d0, 0.d0, 0.d0, -3.23d0, 0.d0, 0.d0, 0.d0, -4.0d-1 /) , (/3,3/))
		K3=spring_conv_coeff*reshape((/  3.00d0, 0.d0, 0.d0, 0.d0, -5.25d0, 0.d0, 0.d0, 0.d0,  1.5d-1 /) , (/3,3/))
		K4=spring_conv_coeff*reshape((/ -1.92d0, 0.d0, 0.d0, 0.d0,  2.29d0, 0.d0, 0.d0, 0.d0, -5.8d-1 /) , (/3,3/))

		! set the angle of 1st, 2nd, 3rd, and 4th nearest neighbors of A atoms
		anglesA1=pi/3.0d0*(/ 0.0d0, 2.0d0, 4.0d0 /)
		anglesA2=pi/6.0d0*(/ 1.0d0, 3.0d0, 5.0d0, 7.0d0, 9.0d0, 1.1d1 /)
		anglesA3=pi/3.0d0*(/ 1.0d0, 3.0d0, 5.0d0 /)
		angle_tmp = atan2(sqrt(3.d0),5.d0)
		anglesA4=(/ angle_tmp, 2.d0*pi/3.d0+angle_tmp, 4.d0*pi/3.d0+angle_tmp, -angle_tmp, 2.d0*pi/3.d0-angle_tmp, 4.d0*pi/3.d0-angle_tmp /)

		anglesB1=anglesA1+pi
		anglesB2=anglesA2+pi
		anglesB3=anglesA3+pi
		anglesB4=anglesA4+pi

		! calculate D_AA, D_AB, D_BA, D_BB for 1st and 3rd nearest neighbors of atoms A and B type.
		do i=1,3
			! 1st nearest neighbor for A atoms
			theta=anglesA1(i)
			Rot=reshape((/ cos(theta), sin(theta) , -sin(theta), cos(theta) /), (/2,2/))
			U=reshape((/ cos(theta), -sin(theta), 0.d0, sin(theta), cos(theta), 0.d0, 0.d0, 0.d0, 1.d0 /) , (/3,3/))
			U_inv=reshape((/ cos(theta), sin(theta), 0.d0, -sin(theta), cos(theta), 0.d0, 0.d0, 0.d0, 1.d0 /) , (/3,3/))
			K_tmp=matmul(U_inv,matmul(K1,U))
			deltaR=matmul(Rot,aCC_vec)
			D_AA=D_AA+dcmplx(K_tmp)
			D_AB=D_AB-exp(i1*dcmplx(dot_product(k,deltaR)))*dcmplx(K_tmp)

			! 3rd nearest neighbor for A atoms
			theta=anglesA3(i)
			Rot=reshape((/ cos(theta), sin(theta) , -sin(theta), cos(theta) /), (/2,2/))
			U=reshape((/ cos(theta), -sin(theta), 0.d0, sin(theta), cos(theta), 0.d0, 0.d0, 0.d0, 1.d0 /) , (/3,3/))
			U_inv=reshape((/ cos(theta), sin(theta), 0.d0, -sin(theta), cos(theta), 0.d0, 0.d0, 0.d0, 1.d0 /) , (/3,3/))
			K_tmp=matmul(U_inv,matmul(K3,U))
			deltaR=matmul(Rot,2.d0*aCC_vec)
			D_AA=D_AA+dcmplx(K_tmp)
			D_AB=D_AB-exp(i1*dcmplx(dot_product(k,deltaR)))*dcmplx(K_tmp)

			! 1st nearest neighbor for B atoms
			theta=anglesB1(i)
			Rot=reshape((/ cos(theta), sin(theta) , -sin(theta), cos(theta) /), (/2,2/))
			U=reshape((/ cos(theta), -sin(theta), 0.d0, sin(theta), cos(theta), 0.d0, 0.d0, 0.d0, 1.d0 /) , (/3,3/))
			U_inv=reshape((/ cos(theta), sin(theta), 0.d0, -sin(theta), cos(theta), 0.d0, 0.d0, 0.d0, 1.d0 /) , (/3,3/))
			K_tmp=matmul(U_inv,matmul(K1,U))
			deltaR=matmul(Rot,aCC_vec)
			D_BB=D_BB+dcmplx(K_tmp)
			D_BA=D_BA-exp(i1*dcmplx(dot_product(k,deltaR)))*dcmplx(K_tmp)

			! 3rd nearest neighbor for B atoms
			theta=anglesB3(i)
			Rot=reshape((/ cos(theta), sin(theta) , -sin(theta), cos(theta) /), (/2,2/))
			U=reshape((/ cos(theta), -sin(theta), 0.d0, sin(theta), cos(theta), 0.d0, 0.d0, 0.d0, 1.d0 /) , (/3,3/))
			U_inv=reshape((/ cos(theta), sin(theta), 0.d0, -sin(theta), cos(theta), 0.d0, 0.d0, 0.d0, 1.d0 /) , (/3,3/))
			K_tmp=matmul(U_inv,matmul(K3,U))
			deltaR=matmul(Rot,2.d0*aCC_vec)
			D_BB=D_BB+dcmplx(K_tmp)
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
			D_AA=D_AA+dcmplx(K_tmp) -exp(i1*dcmplx(dot_product(k,deltaR)))*dcmplx(K_tmp)

			! 4rd nearest neighbor for A atoms
			theta=anglesA4(i)
			Rot=reshape((/ cos(theta), sin(theta) , -sin(theta), cos(theta) /), (/2,2/))
			U=reshape((/ cos(theta), -sin(theta), 0.d0, sin(theta), cos(theta), 0.d0, 0.d0, 0.d0, 1.d0 /) , (/3,3/))
			U_inv=reshape((/ cos(theta), sin(theta), 0.d0, -sin(theta), cos(theta), 0.d0, 0.d0, 0.d0, 1.d0 /) , (/3,3/))
			K_tmp=matmul(U_inv,matmul(K4,U))
			deltaR=matmul(Rot,sqrt(7.d0)*aCC_vec)
			D_AA=D_AA+dcmplx(K_tmp)
			D_AB=D_AB-exp(i1*dcmplx(dot_product(k,deltaR)))*dcmplx(K_tmp)

			! 2st nearest neighbor for B atoms
			theta=anglesB2(i)
			Rot=reshape((/ cos(theta), sin(theta) , -sin(theta), cos(theta) /), (/2,2/))
			U=reshape((/ cos(theta), -sin(theta), 0.d0, sin(theta), cos(theta), 0.d0, 0.d0, 0.d0, 1.d0 /) , (/3,3/))
			U_inv=reshape((/ cos(theta), sin(theta), 0.d0, -sin(theta), cos(theta), 0.d0, 0.d0, 0.d0, 1.d0 /) , (/3,3/))
			K_tmp=matmul(U_inv,matmul(K2,U))
			deltaR=matmul(Rot,sqrt(3.d0)*aCC_vec)
			D_BB=D_BB+dcmplx(K_tmp) -exp(i1*dcmplx(dot_product(k,deltaR)))*dcmplx(K_tmp)

			! 4rd nearest neighbor for B atoms
			theta=anglesB4(i)
			Rot=reshape((/ cos(theta), sin(theta) , -sin(theta), cos(theta) /), (/2,2/))
			U=reshape((/ cos(theta), -sin(theta), 0.d0, sin(theta), cos(theta), 0.d0, 0.d0, 0.d0, 1.d0 /) , (/3,3/))
			U_inv=reshape((/ cos(theta), sin(theta), 0.d0, -sin(theta), cos(theta), 0.d0, 0.d0, 0.d0, 1.d0 /) , (/3,3/))
			K_tmp=matmul(U_inv,matmul(K4,U))
			deltaR=matmul(Rot,sqrt(7.d0)*aCC_vec)
			D_BB=D_BB+dcmplx(K_tmp)
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

		! fix the phase of u_ph
		do i=1,6
			u_ph(:,i) = u_ph(:,i)*exp(-i1*dcmplx(atan2(aimag(u_ph(1,i)),real(u_ph(1,i)))))
		end do

		! change the units from [1/cm] to [Joules]
		omega=omega*(4.1357d-3/3.3356d1)*eV

	end subroutine graphene_phonon

	!***************************************************************************
	! - this subroutine calculates electron-phonon matrix element in graphene:
	!	M = f_tilde_1(k-q,q,ib) + f_tilde_2(k,q,ib)
	!***************************************************************************

	subroutine graphene_electron_phonon_matrix_element(matrix_element,k,q,ib,a1,a2)
		use constants_mod, only: i1

		complex*16 :: matrix_element
		real*8, dimension(2), intent(in) :: k,q
		integer, intent(in) :: ib
		real*8, dimension(2), intent(in) :: a1,a2
		complex*16 :: f_tilde_1, f_tilde_2
		complex*16 :: f_1, f_2
		complex*16 :: f_k, f_k_minus_q
		real*8, dimension(2) :: e1,e2,e3
		complex*16, dimension(2) :: e1_normalized,e2_normalized,e3_normalized
		real*8, dimension(2) :: aCC_vec
		complex*16, dimension(3) :: eA, eB
		complex*16, dimension(6,6) :: u_ph
		real*8, dimension(6) :: omega

		f_tilde_1 = (0.d0, 0.d0)
		f_tilde_2 = (0.d0, 0.d0)
		f_1 = (0.d0, 0.d0)
		f_2 = (0.d0, 0.d0)

		e1 = (a1+a2)/3.d0
		e2 = (a1-2.d0*a2)/3.d0
		e3 = (a2-2.d0*a1)/3.d0

		e1_normalized = dcmplx(e1/sqrt(dot_product(e1,e1)))
		e2_normalized = dcmplx(e2/sqrt(dot_product(e2,e2)))
		e3_normalized = dcmplx(e3/sqrt(dot_product(e3,e3)))

		aCC_vec = e1

		call graphene_phonon(omega,u_ph,q,aCC_vec)

		eA = u_ph(1:3,ib)
		eB = u_ph(4:6,ib)

		! I am not sure at this time if we need to normalize eA and eB the following to lines do the normalization if needed.
		! eA = eA/sqrt(dot_product(eA,eA))
		! eB = eB/sqrt(dot_product(eB,eB))

		! R0 = 0
		! Ru' = 0
		! R0A = 0
		! Ru'B = e1
		f_1 = f_1 + sum(e1_normalized*(eA(1:2)-eB(1:2)))*exp(+i1*dcmplx(dot_product(k,e1)))
		f_2 = f_2 + sum(e1_normalized*(eA(1:2)-eB(1:2)))*exp(-i1*dcmplx(dot_product(k-q,e1)))

		! R0 = 0
		! Ru' = -a2
		! R0A = 0
		! Ru'B = e2
		f_1 = f_1 + sum(e2_normalized*(eA(1:2)-eB(1:2)*exp(-i1*dcmplx(dot_product(q,-a2)))))*exp(+i1*dcmplx(dot_product(k,e2)))
		f_2 = f_2 + sum(e2_normalized*(eA(1:2)-eB(1:2)*exp(-i1*dcmplx(dot_product(q,-a2)))))*exp(-i1*dcmplx(dot_product(k-q,e2)))

		! R0 = 0
		! Ru' = -a1
		! R0A = 0
		! Ru'B = e3
		f_1 = f_1 + sum(e3_normalized*(eA(1:2)-eB(1:2)*exp(-i1*dcmplx(dot_product(q,-a1)))))*exp(+i1*dcmplx(dot_product(k,e3)))
		f_2 = f_2 + sum(e3_normalized*(eA(1:2)-eB(1:2)*exp(-i1*dcmplx(dot_product(q,-a1)))))*exp(-i1*dcmplx(dot_product(k-q,e3)))



		f_k=exp(i1*dcmplx(dot_product(k,e1)))+exp(i1*dcmplx(dot_product(k,e2)))+exp(i1*dcmplx(dot_product(k,e3)))
		f_tilde_1 = f_1 * conjg(f_k)/dcmplx(abs(f_k))

		f_k_minus_q=exp(i1*dcmplx(dot_product(k-q,e1)))+exp(i1*dcmplx(dot_product(k-q,e2)))+exp(i1*dcmplx(dot_product(k-q,e3)))
		f_tilde_2 = f_2 * f_k_minus_q/dcmplx(abs(f_k_minus_q))

		matrix_element = f_tilde_1 + f_tilde_2

	end subroutine graphene_electron_phonon_matrix_element

end module graphene_mod
