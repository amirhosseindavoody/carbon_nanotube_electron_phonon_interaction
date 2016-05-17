!**********************************************************************************************************************
! This subroutine reads the dispersion of excitons.
!**********************************************************************************************************************
subroutine fnScattering()
  use comparams
  implicit none
  
  integer :: mu1, mu2
  integer :: ik1, ik2
  integer :: i_branch
  integer :: imu1, imu2
  real*8 :: tmpr1, tmpr2
  real*8 :: dos, nq
  real*8 :: ScatRate, TotScatRate
  real*8, dimension(2) :: q, kc1, kc2
  real*8, dimension(6) :: omega_q
  complex*16 :: f1, f2, MatElm
  complex*16, dimension(6,6) :: u_ph
  
  do ik1=1,ikc_max
    Ek(:,-ik1,:)=Ek(:,ik1,:)
  enddo
  
  do mu1=1,Nu/2-1
    Ek(-mu1,:,:)=Ek(mu1,:,:)
  enddo
  
  mu1=mu0

  do ik1=0,ikc_max
    print *, 'Calculating el-ph scattering rates -->', ik1
    TotScatRate=0
    do imu1=0,1
      mu1=(2*imu1-1)*mu0
      do mu2=1-Nu/2,Nu/2 
        do ik2=ikc_min,0
          !if ((mu1 .ne. mu2) .or. (ik1 .ne. ik2)) then
            tmpr1=(Ek(mu2,ik2-1,1)-Ek(mu1,ik1,1))
            tmpr2=(Ek(mu2,ik2+1,1)-Ek(mu1,ik1,1))
            if ((tmpr1*tmpr2) .lt. 0.d0) then
              dos=2.d0*dk/abs(Ek(mu2,ik2+1,1)-Ek(mu2,ik2-1,1))
          
              kc1=dble(mu1)*K1+dble(ik1)*dk*K2
              kc2=dble(mu2)*K1+dble(ik2)*dk*K2
              q=kc1-kc2
              call fnGraphenePhonon(omega_q,u_ph,q)
          
              do i_branch=1,3
                nq=1.d0/(exp(omega_q(i_branch)/kbt)-1.d0)
            
                call fnCalculateF1(f1,kc1,q,u_ph(:,i_branch))
                call fnCalculateF2(f2,kc2,q,u_ph(:,i_branch))
            
                MatElm=f1+f2
            
                ScatRate=2.d0*pi*(g0**2.d0)*(hb*nq)/(2.d0*Nu*mass*omega_q(i_branch))*norm2(t_vec)*dos*((abs(MatElm))**2.d0)
            
                !ScatRate=((abs(MatElm))**2.d0)
                !ScatRate=dos
                !ScatRate=(nq)/(omega_q(i_branch))*dos*((abs(MatElm))**2.d0)
  
            
                TotScatRate=TotScatRate+ScatRate
              enddo
            endif
          !endif
        enddo
      enddo
    enddo

    if (ik1 .le. 0) then
      write(fh4,10) Ek(mu1,ik1,1), TotScatRate
      write(fh6,11) ik1, TotScatRate
    endif
    
    if (ik1 .ge. 0) then
      write(fh5,10) Ek(mu1,ik1,1), TotScatRate
      write(fh6,11) ik1, TotScatRate
    endif
    
  enddo
  
  
  
10 FORMAT (E16.8,E16.8)
11 FORMAT (I5.4,E16.8)
  
  return
end