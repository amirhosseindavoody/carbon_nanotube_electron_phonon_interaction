module comparams
	use cnt_class, only : cnt
	implicit none


	real :: starttime,endtime !duration of simulation

	real*8 :: Temperature = 300.d0 !Temperature of the system in Kelvin units

	real*8 :: ppLen = 10.d-9 !length per perpendicular tubes

	type (cnt) :: cnt1, cnt2

end module comparams