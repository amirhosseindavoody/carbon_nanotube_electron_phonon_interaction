subroutine fnExPhScatRate(ScatRate,iKcm1,iX1,iKcm2,iX2,i_branch,ratio1,ratio2)
  use comparams
  implicit none

  integer :: ikr1, ikr2
  integer :: iKcm1, iKcm2, iX1, iX2, i_branch
  real*8 :: ratio1, ratio2
  real*8 :: ScatRate
  real*8 :: nq, DOS
  real*8, dimension(2) :: k,q
  real*8, dimension(6) :: omega_tmp
  complex*16 :: MatElm, Psi2_tmp, f1, f2
  complex*16, dimension(:), allocatable :: Psi2
  complex*16, dimension(6,6) :: u_ph
  
  allocate(Psi2(ikr_low:ikr_high))
  
  MatElm=(0.d0,0.d0)
  
  Psi2=ratio2*Psi0_A2(iX2,iKcm2,:)+ratio1*Psi0_A2(iX2,iKcm2+1,:)
  Psi2=Psi2/sqrt(norm2(abs(Psi2)))
  
  q=2.d0*(dble(iKcm1)-ratio2*dble(iKcm2)-ratio1*dble(iKcm2+1))*dk*K2
  call fnGraphenePhonon(omega_tmp,u_ph,q)
  
  do  ikr1=(ikr_low+abs(iKcm2-iKcm1)),(ikr_high-abs(iKcm2-iKcm1)-1)
    ikr2=ikr1-iKcm1+iKcm2
    Psi2_tmp=ratio2*Psi2(ikr2)+ratio1*Psi2(ikr2+1)
    
    k=dble(min_sub)*K1+dble(ikr1-iKcm1+2*(iKcm2+ratio1))*dk*K2
    
    
    
    call fnCalculateF1(f1,k+q,q,u_ph(:,i_branch))
    call fnCalculateF2(f2,k,q,u_ph(:,i_branch))
    
    MatElm=MatElm+conjg(Psi0_A2(iX1,iKcm1,ikr1))*Psi2_tmp*(f1+f2)
  enddo
  
  nq=1.d0/(exp(omega_tmp(i_branch)/kbt)-1.d0)
  DOS=dk/abs(Ex0_A2(iX2,iKcm2)-Ex0_A2(iX2,iKcm2+1))
  !ScatRate=2.d0*pi*(g0**2.d0)*(hb*nq)/(2.d0*Nu*mass*omega_tmp(i_branch))*norm2(t_vec)*DOS*(abs(MatElm))**2.d0
  ScatRate=DOS
  
  return
end