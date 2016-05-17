subroutine fnCNTelectron()
  use comparams
  implicit none
  
  integer :: nkc
  integer :: imin_sub
  integer :: i,j,mu,ik,tmpi
  integer, dimension(:), allocatable :: min_loc, min_sub
  real*8 :: tmpr
  real*8, dimension(2) :: k,E1_tmp,E2_tmp
  real*8, dimension(:), allocatable :: k_vec,min_energy
  complex*16, dimension(2) :: Cc_tmp,Cv_tmp
  
  
  ! calculate CNT energy dispersion.***********************************************************************************
  ikc_max=floor(pi/norm2(t_vec)/dk)
  ikc_min=-ikc_max
  nkc=2*ikc_max+1
  
  allocate(k_vec(ikc_min:ikc_max))
  allocate(Ek(1-Nu/2:Nu/2,ikc_min:ikc_max,2))
  allocate(min_loc(0:Nu/2))
  
  do ik=ikc_min,ikc_max
    k_vec(ik)=dble(ik)*dk
  end do
  
  do mu=1-Nu/2,Nu/2
    do ik=ikc_min,ikc_max
      k=dble(mu)*K1+dble(ik)*dk*K2
      call fnGrapheneElectron(Ek(mu,ik,:),Cc_tmp,Cv_tmp,k)
    enddo
  enddo
  
  ! save the CNT energy dispersion*************************************************************************************
  do ik=ikc_min,ikc_max
     write(fh2,10, advance='no') k_vec(ik) 
     write(fh3,10, advance='no') k_vec(ik)
  end do
  write(fh2,10)
  write(fh3,10)
  
  do mu=1-Nu/2,Nu/2
    do ik=ikc_min,ikc_max
      write(fh2,10, advance='no') Ek(mu,ik,1) 
      write(fh3,10, advance='no') Ek(mu,ik,2) 
    end do
    write(fh2,10)
    write(fh3,10)
  enddo
  
  ! find the subbands with a minimum energy.***************************************************************************
  min_loc=minloc(Ek(0:Nu/2,:,1),2)
  imin_sub=count((min_loc .lt. nkc) .and. (min_loc .gt. 1))
  allocate(min_sub(imin_sub))
  allocate(min_energy(imin_sub))
  
  i=1
  do mu=0,Nu/2
    if ((min_loc(mu) .gt. 1) .and. (min_loc(mu) .lt. nkc)) then
       min_sub(i)=mu
       min_energy(i)=minval(Ek(mu,:,1))
       i=i+1
    end if
  end do
  
  ! sort the subbands
  do i=imin_sub,2,-1
    do j=i-1,1,-1
      if (min_energy(i) .lt. min_energy(j)) then
        tmpr=min_energy(i)
        tmpi=min_sub(i)
        min_energy(i)=min_energy(j)
        min_sub(i)=min_sub(j)
        min_energy(j)=tmpr
        min_sub(j)=tmpi
      end if    
    end do
  end do
  
  mu0=min_sub(1)
  
10 FORMAT (E16.8)
  
  continue
  return
end