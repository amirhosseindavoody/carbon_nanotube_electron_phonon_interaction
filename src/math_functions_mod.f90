module math_functions_mod
	implicit none
	private
	public :: gcd, bessk0, bisect_root, my_norm2, eig, polint, find_all_roots, first_derivative

contains
	!**********************************************************************************************************************
	! This subroutine calculates the greatest common divisor of the arguments na and nb
	!**********************************************************************************************************************

	subroutine gcd(ngcd,na,nb)
	integer, intent(in) :: na,nb
	integer, intent(out) :: ngcd
	integer :: ia,ib,itemp

		ia=na
		ib=nb
		do while (ib .ne. 0)
			itemp=ia
			ia=ib
			ib=mod(itemp,ib)
		end do
		ngcd=ia
	end subroutine gcd

	!***************************************************************************
	! This function calculates the modified bessel function of the first kind with parameter nu=0: I0(x)
	!***************************************************************************

	real*8 function bessi0(x)
		real*8 :: x
		real*8 :: ax
		real*8 :: y
		real*8, save :: p1, p2, p3, p4, p5, p6, p7, q1, q2, q3, q4, q5, q6, q7, q8, q9

		data p1, p2, p3, p4, p5, p6, p7 /1.d0, 3.5156229d0, 3.0899424d0, 1.2067492d0, 0.2659732d0, 0.360768d-1, 0.45813d-2/
		data q1, q2, q3, q4, q5, q6, q7, q8, q9 /0.39894228d0, 0.1328592d-1, 0.225319d-2, -0.157565d-2, 0.916281d-2, -0.2057706d-1, 0.2635537d-1, -0.1647633d-1, 0.392377d-2/

		if (abs(x) .lt. 3.75d0) then
			y = (x/3.75d0)**2
			bessi0 = p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))))
		else
			ax = abs(x)
			y = 3.75d0/ax
			bessi0 = (exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))))
		end if
		return
	end function bessi0

	!***************************************************************************
	! This function calculates the modified bessel function of the second kind with parameter nu=0: K0(x)
	!***************************************************************************

	real*8 function bessk0(x)
		real*8 :: x
		real*8 :: y
		real*8, save :: p1, p2, p3, p4, p5, p6, p7, q1, q2, q3, q4, q5, q6, q7

		data p1, p2, p3, p4, p5, p6, p7 /-0.57721566d0, 0.42278420d0, 0.23069756d0, 0.3488590d-1, 0.262698d-2, 0.10750d-3, 0.74d-5/
		data q1, q2, q3, q4, q5, q6, q7 /1.25331414d0, -0.7832358d-1, 0.2189568d-1, -0.1062446d-1, 0.587872d-2, -0.251540d-2, 0.53208d-3/

		if (x .le. 2.d0) then
			y = x*x/4.d0
			bessk0=(-log(x/2.d0)*bessi0(x))+(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
		else
			y=2.d0/x
			bessk0=(exp(-x)/sqrt(x))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*q7))))))
		end if
		return
	end function bessk0

	!***************************************************************************
	! This subroutine calculates crossing point of an array using bisection method
	! here we use the fact that the array is decreasing with respect to its index
	! ya: input array
	! n: size of the array ya.
	! ind: index of the element that we are looking for.
	! x0: the value that we want find to find index that ya(ind)=x0
	!***************************************************************************

	subroutine bisect_root (n, ya, x0, ind)
		integer, intent(in) :: n
		integer, intent(out) :: ind
		real*8, intent(in) :: x0
		real*8, dimension(n), intent(in) :: ya

		integer :: ix_lower, ix_upper, ix_mid
		real*8, dimension(n) :: tmpArray

		ind = 0
		ix_lower=1
		ix_upper=n
		tmpArray = ya-x0

		if (tmpArray(ix_lower) * tmpArray(ix_upper) .ge. 0.d0) then
			if (tmpArray(ix_lower) .eq. 0.d0) then
				ind = ix_lower
				return
			elseif (tmpArray(ix_upper) .eq. 0.d0) then
				ind = ix_upper
				return
			else
				ind = 0
				return
			endif
		elseif (tmpArray(ix_lower) .gt. 0.d0) then
			do while((ix_upper-ix_lower) .gt. 1)
				ix_mid = (ix_upper + ix_lower)/2
				if (tmpArray(ix_mid) .gt. 0.d0) then
					ix_lower = ix_mid
				elseif(tmpArray(ix_mid) .lt. 0.d0) then
					ix_upper = ix_mid
				else
					ind = ix_mid
					return
				endif
			enddo
			ind = ix_upper
			return
		else
			do while((ix_upper-ix_lower) .gt. 1)
				ix_mid = (ix_upper + ix_lower)/2
				if (tmpArray(ix_mid) .gt. 0.d0) then
					ix_upper = ix_mid
				else if(tmpArray(ix_mid) .lt. 0.d0) then
					ix_lower = ix_mid
				else
					ind = ix_mid
					return
				endif
			enddo
			ind = ix_lower
			return
		endif

		write(*,*) "Error in finding bisection root!!!!"
		call exit()

	end subroutine bisect_root

	!***************************************************************************
	! This subroutine finds ALL the crossing point of an array from value y0
	! ya: input array
	! y0: the value that we want find to find index that ya(ind)=y0
	!***************************************************************************

	subroutine find_all_roots (ya, l_bound, u_bound, y0, n_root, root_idx)
		integer, intent(in) :: l_bound
		integer, intent(in) :: u_bound
		real*8, intent(in) :: y0
		real*8, dimension(l_bound:u_bound), intent(in) :: ya
		integer, intent(out) :: n_root
		integer, dimension(:), intent(out) :: root_idx

		integer :: i
		real*8, dimension(l_bound:u_bound) :: tmpArray

		tmpArray = ya - y0
		root_idx = l_bound-1
		n_root = 0

		do i = l_bound,u_bound-1

			if ((tmpArray(i)*tmpArray(i+1)) .lt. -tiny(0.d0)) then

				n_root = n_root+1

				if (abs(tmpArray(i)) .le. abs(tmpArray(i+1))) then
					root_idx(n_root) = i
				else
					root_idx(n_root) = i+1
				endif

			elseif (abs(tmpArray(i)) .lt. tiny(0.d0)) then

				n_root = n_root+1
				root_idx(n_root) = i

			endif

		enddo

		if (abs(tmpArray(u_bound)) .lt. tiny(0.d0)) then

			n_root = n_root+1
			root_idx(n_root) = u_bound

		endif

	end subroutine find_all_roots


	!***************************************************************************
	! This function calculates the magnitude of a 2D real*8 vector
	!***************************************************************************

	real*8 function my_norm2(my_vec)
		real*8, dimension(2), intent(in) :: my_vec

		my_norm2 = sqrt(my_vec(1)*my_vec(1)+my_vec(2)*my_vec(2))

		return
	end function my_norm2

	!***************************************************************************
	! This subroutine calculates the eigen values and eigen vectors of matrix given the size of it (nl) and gives back eigen vectors in A.
	!***************************************************************************

	subroutine eig(nl,matrix,A,W)

		interface
			subroutine ZHEEV( JOBZ, UPLO,	N, A, LDA, W, WORK, LWORK, RWORK, INFO )
				character (len=1) :: JOBZ, UPLO
				integer :: INFO
				integer :: LDA, LWORK, N
				complex*16, dimension(N,N) :: A
				real*8, dimension(N) :: W
				complex*16, dimension(2*N-1) :: WORK
				real*8, dimension(3*N-2) :: RWORK
			end subroutine ZHEEV
		end interface

		integer, intent(in) :: nl
		complex*16, dimension(nl,nl), intent(inout) :: matrix
		complex*16, dimension(nl,nl), intent(out) :: A
		complex*16, dimension(2*nl-1) :: WORK
		character (len=1) :: JOBZ, UPLO
		integer :: INFO, LDA, LWORK, N
		real*8, dimension(nl) :: W
		real*8, dimension(3*nl-2) :: RWORK

		JOBZ = 'V'
		UPLO = 'L'
		N = size(matrix,1)
		A = matrix
		LDA = N
		LWORK = 2*size(matrix,1)-1

		call ZHEEV( JOBZ, UPLO,	N, A, LDA, W, WORK, LWORK, RWORK, INFO )

		if (INFO .ne. 0) then
			write(*,"(A50,I1.1)") "ERROR: calculation of eigen value failed , INFO = ", INFO
			call exit()
		end if

	end subroutine eig

	!***************************************************************************
	! This subroutine interpolated arrays xa and ya using polynomial interpolation technique.
	! Given arrays xa and ya, each of length n, and given a value x, this routine returns a value y, and an error estimate.
	! if P(x) is athe polynomial of degree n-1 such that P(xa) = ya, for all values of arrays xa and ya, then the returned
	! value y = P(x).
	!***************************************************************************

	subroutine polint(xa, ya, n, x, y, dy)
		integer, intent(in) :: n
		real*8, intent(in) :: x
		real*8, intent(out) :: y, dy
		real*8, dimension(n), intent(in) :: xa, ya

		integer, parameter :: nmax=10
		integer :: i, m, ns
		real*8 :: den, dif, dift, ho, hp, w
		real*8, dimension(n) :: c, d

		ns = 1
		dif = abs(x-xa(1))

		do i=1,n
			dift = abs(x-xa(i))
			if (dift .lt. dif) then
				ns = i
				dif = dift
			endif
			c(i) = ya(i)
			d(i) = ya(i)
		enddo

		y = ya(ns)
		ns = ns-1

		do m=1,n-1
			do i=1, n-m
				ho=xa(i)-x
				hp=xa(i+m)-x
				w=c(i+1)-d(i)
				den=ho-hp
				if(den .eq. 0) then
					write(*,*) "Faliure in polint"
					call exit()
				endif
				den = w/den
				d(i) = hp*den
				c(i) = ho*den
			enddo
			if (2*ns .lt. n-m) then
				dy = c(ns+1)
			else
				dy=d(ns)
				ns=ns-1
			endif
			y=y+dy
		enddo
		return
	end subroutine polint

	!***************************************************************************
	! This subroutine calculates the derivative of a function
	! ya: input array
	! l_bound: the lower limit of indices for array ya
	! u_bound: the upper limit of indices for array ya
	! point_index: index of the point in which we want to calculate the derivative
	! dx: mesh point size in the discretization process
	! derivative: the calculated first derivative dy/dx
	!***************************************************************************

	subroutine first_derivative (ya, l_bound, u_bound, point_index, dx, derivative, single)
		integer, intent(in) :: l_bound
		integer, intent(in) :: u_bound
		integer, intent(in) ::point_index
		real*8, intent(in) :: dx
		real*8, dimension(l_bound:u_bound), intent(in) :: ya
		real*8, intent(out) :: derivative
		logical, optional :: single

		logical :: single_used

! 		if (point_index .lt. l_bound+2 ) then
! 			derivative = (-ya(point_index+2)+4.d0*ya(point_index+1)-3.d0*ya(point_index))/(2.d0*dx)
! 		elseif (point_index .gt. u_bound-2) then
! 			derivative = (3.d0*ya(point_index)-4.d0*ya(point_index-1)+ya(point_index-2))/(2.d0*dx)
! 		else
! 			derivative = (-ya(point_index+2)+8.d0*ya(point_index+1)-8.d0*ya(point_index-1)+ya(point_index-2))/(12.d0*dx)
! 		endif

		

! 		if (point_index .gt. 0.d0 ) then
! 			derivative = (-ya(point_index+2)+4.d0*ya(point_index+1)-3.d0*ya(point_index))/(2.d0*dx)
! 		else
! 			derivative = (3.d0*ya(point_index)-4.d0*ya(point_index-1)+ya(point_index-2))/(2.d0*dx)
! 		endif

		single_used = .false.

		if (present(single)) single_used = single

		if ((point_index .gt. u_bound-4) .or. (point_index .lt. l_bound+4)) then
			single_used = .true.
		endif

		if (single_used) then

			if (point_index .gt. 0.d0) then
				if (point_index .le. u_bound-6) then
					derivative = (-49.d0/20.d0)*ya(point_index)+(+6.d0)*ya(point_index+1)+(-15.d0/2.d0)*ya(point_index+2)+(+20.d0/3.d0)*ya(point_index+3)+(-15.d0/4.d0)*ya(point_index+4)+(+6.d0/5.d0)*ya(point_index+5)+(-1.d0/6.d0)*ya(point_index+6)
				else
					derivative = (+49.d0/20.d0)*ya(point_index)+(-6.d0)*ya(point_index-1)+(+15.d0/2.d0)*ya(point_index-2)+(-20.d0/3.d0)*ya(point_index-3)+(+15.d0/4.d0)*ya(point_index-4)+(-6.d0/5.d0)*ya(point_index-5)+(+1.d0/6.d0)*ya(point_index-6)
				endif
			else
				if (point_index .ge. l_bound+6) then
					derivative = (+49.d0/20.d0)*ya(point_index)+(-6.d0)*ya(point_index-1)+(+15.d0/2.d0)*ya(point_index-2)+(-20.d0/3.d0)*ya(point_index-3)+(+15.d0/4.d0)*ya(point_index-4)+(-6.d0/5.d0)*ya(point_index-5)+(+1.d0/6.d0)*ya(point_index-6)
				else
					derivative = (-49.d0/20.d0)*ya(point_index)+(+6.d0)*ya(point_index+1)+(-15.d0/2.d0)*ya(point_index+2)+(+20.d0/3.d0)*ya(point_index+3)+(-15.d0/4.d0)*ya(point_index+4)+(+6.d0/5.d0)*ya(point_index+5)+(-1.d0/6.d0)*ya(point_index+6)
				endif
			endif

		else
			derivative = (1.d0/280.d0)*ya(point_index-4) + (-4.d0/105.d0)*ya(point_index-3) + (1.d0/5.d0)*ya(point_index-2) + (-4.d0/5.d0)*ya(point_index-1) + (4.d0/5.d0)*ya(point_index+1) + (-1.d0/5.d0)*ya(point_index+2) + (4.d0/105.d0)*ya(point_index+3) + (-1.d0/280.d0)*ya(point_index+4)
		endif

		derivative = derivative/dx

	end subroutine first_derivative

end module math_Functions_mod
