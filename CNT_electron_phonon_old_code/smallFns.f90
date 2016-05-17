!**********************************************************************************************************************
! This subroutine calculates the greatest common divisor of the arguments na and nb
!**********************************************************************************************************************
subroutine gcd(ngcd,na,nb)
  integer :: na,nb,ngcd
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
end