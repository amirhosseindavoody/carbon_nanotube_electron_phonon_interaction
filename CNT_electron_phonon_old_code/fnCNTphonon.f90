subroutine fnCNTphonon()
  use comparams
  implicit none
  
  integer :: ik, iq, nkc
  integer :: mu, ib
  
  real*8, dimension(2) :: k, q
  real*8, dimension(:), allocatable :: k_vec, q_vec
  real*8, dimension (:,:,:), allocatable :: omega_tmp
  complex*16, dimension(6,6) :: u_ph
  
  ! calculate CNT phonon dispersion.***********************************************************************************
  ikc_max=floor(pi/norm2(t_vec)/dk)
  ikc_min=-ikc_max
  nkc=2*ikc_max+1
  
  allocate(k_vec(ikc_min:ikc_max))
  allocate(omega_tmp(1-Nu/2:Nu/2,ikc_min:ikc_max,6))
  
  do ik=ikc_min,ikc_max
    k_vec(ik)=dble(ik)*dk
  enddo
  
  do mu=1-Nu/2,Nu/2
    do ik=ikc_min,ikc_max
      k=dble(mu)*K1+dble(ik)*dk*K2
      call fnGraphenePhonon(omega_tmp(mu,ik,:),u_ph,k)
    enddo
  enddo
  
  ! save the CNT phonon dispersion*************************************************************************************
  do ik=ikc_min,ikc_max
     write(fh1,10, advance='no') k_vec(ik) 
  enddo
  write(fh1,10)
  
  do ib=1,6
    do mu=1-Nu/2,Nu/2
      do ik=ikc_min,ikc_max
        write(fh1,10, advance='no') omega_tmp(mu,ik,ib) 
      end do
      write(fh1,10)
    enddo
  enddo
  
  10 FORMAT (E16.8)
  
  return
end