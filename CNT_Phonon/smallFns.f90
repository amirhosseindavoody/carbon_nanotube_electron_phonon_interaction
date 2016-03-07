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
    
!**********************************************************************************************************************
! This subroutine reads the limits of different types of k-vectors in the simulated exciton energy code.
!**********************************************************************************************************************
subroutine fnReadMisc()
  use comparams
  implicit none

  read(fh2,10) min_sub
  read(fh2,10) ik_max
  read(fh2,10) iKcm_max
  read(fh2,10) ikr_high
  read(fh2,10) ik_high
  read(fh2,10) iq_max !iq_max is a redundant variable here which will be rewritten later.
  read(fh2,10) nX
  read(fh2,11) dk
  
  ik_min = -ik_max
  iKcm_min = -iKcm_max
  ikr_low = -ikr_high
  ik_low = -ik_high
  iq_min = -iq_max !iq_min is a redundant variable here which will be rewritten later.
  
10 FORMAT (I4.4) 
11 FORMAT (E16.8)

  return
end
    
!**********************************************************************************************************************
! This subroutine reads the dispersion of excitons.
!**********************************************************************************************************************
subroutine fnReadExDispersion()
  use comparams
  implicit none
  
  integer :: iKcm, ix, ikr

  allocate(Ex_A1(1:nX,iKcm_min:iKcm_max))
  allocate(Ex0_A2(1:nX,iKcm_min:iKcm_max))
  allocate(Ex1_A2(1:nX,iKcm_min:iKcm_max))
  
  allocate(Psi_A1(1:nX,iKcm_min:iKcm_max,ikr_low:ikr_high))
  allocate(Psi0_A2(1:nX,iKcm_min:iKcm_max,ikr_low:ikr_high))
  allocate(Psi1_A2(1:nX,iKcm_min:iKcm_max,ikr_low:ikr_high))
  
  do iKcm=iKcm_min,iKcm_max
    print *,"Exciton Disp. --> iKcm=", iKcm
    
    ! read exciton energy and wavefunction
    do iX=1,nX
      read(fh3,10, advance='no') Ex_A1(iX,iKcm)
      read(fh4,10, advance='no') Ex0_A2(iX,iKcm)
      read(fh5,10, advance='no') Ex1_A2(iX,iKcm)
      do ikr=ikr_low,ikr_high
        read(fh6,11, advance='no') Psi_A1(iX,iKcm,ikr)
        read(fh7,11, advance='no') Psi0_A2(iX,iKcm,ikr)
        read(fh8,11, advance='no') Psi1_A2(iX,iKcm,ikr)
      enddo
    enddo
    read(fh3,10)
    read(fh4,10)
    read(fh5,10)
    read(fh6,10)
    read(fh7,10)
    read(fh8,10)

  enddo
  
  10 FORMAT (E16.8)  
  11 FORMAT (E16.8,E16.8) 
  
  return
end