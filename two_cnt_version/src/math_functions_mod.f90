module math_functions_mod
	implicit none
	private
	public :: gcd, bessk0, bisect_root, my_norm2

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
	return
	end subroutine gcd
		
	!**********************************************************************************************************************
	! This function calculates the modified bessel function of the first kind with parameter nu=0: I0(x)
	!**********************************************************************************************************************
	
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

	!**********************************************************************************************************************
	! This function calculates the modified bessel function of the second kind with parameter nu=0: K0(x)
	!**********************************************************************************************************************
	
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

	!**********************************************************************************************************************
	! This subroutine calculates crossing point of an array using bisection method, here we use the fact that the array is decreasing with respect to its index
	!**********************************************************************************************************************

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

	!**********************************************************************************************************************
	! This function calculates the magnitude of a 2D real*8 vector
	!**********************************************************************************************************************
	
	real*8 function my_norm2(my_vec)
		real*8, dimension(2), intent(in) :: my_vec
		
		my_norm2 = sqrt(my_vec(1)*my_vec(1)+my_vec(2)*my_vec(2))

		return 
	end function my_norm2

end module math_Functions_mod