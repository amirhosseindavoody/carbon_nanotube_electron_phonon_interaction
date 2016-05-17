!**********************************************************************************************************************
subroutine fnCalculatePhase(phase,k)
  use comparams, only: a1, a2, i1
  implicit none
  
  real*8, dimension(2) :: k
  complex*16 :: f_k
  complex*16 :: phase
  
  f_k=exp(i1*dot_product(k,(a1+a2)/3.d0))+exp(i1*dot_product(k,(a1-2.d0*a2)/3.d0))+exp(i1*dot_product(k,(a2-2.d0*a1)/3.d0))
  phase=f_k/abs(f_k)
  
  return
end

!**********************************************************************************************************************    
subroutine fnCalculateF1(f1,k,q,u_ph)
  use comparams
  implicit none
  
  real*8 :: norm
  real*8, dimension(2) :: k,q
  complex*16, dimension(3) :: e1, e2, e3
  complex*16, dimension(3) :: eAq, eBq
  complex*16 :: f1, phase
  complex*16, dimension(6) :: u_ph
  
  
  e1=(/a1(1)+a2(1), a1(2)+a2(2), 0 /)
  norm=dot_product(e1,e1)
  e1=e1/sqrt(norm)
  
  e2=(/a1(1)-2.d0*a2(1), a1(2)-2.d0*a2(2), 0 /)
  norm=dot_product(e2,e2)
  e2=e2/sqrt(norm)
  
  e3=(/a2(1)-2.d0*a1(1), a2(2)-2.d0*a1(2), 0 /)
  norm=dot_product(e3,e3)
  e3=e3/sqrt(norm)
  
  eAq=u_ph(1:3)
  norm=dot_product(eAq,eAq)
  eAq=eAq/sqrt(norm)
  
  eBq=u_ph(4:6)
  norm=dot_product(eBq,eBq)
  eBq=eBq/sqrt(norm)
  
  f1=(0.d0,0.d0)
  
  f1=f1+exp(i1*dot_product(k-q,(a1+a2)/3.d0))*dot_product(e1,eAq-eBq)
  f1=f1+exp(i1*dot_product(k-q,(a1-2.d0*a2)/3.d0))*dot_product(e2,eAq-eBq*exp(i1*dot_product(q,-a2)))
  f1=f1+exp(i1*dot_product(k-q,(a2-2.d0*a1)/3.d0))*dot_product(e3,eAq-eBq*exp(i1*dot_product(q,-a1)))
  
  call fnCalculatePhase(phase,k-q)
  f1=f1*conjg(phase)
  return
end

!**********************************************************************************************************************    
subroutine fnCalculateF2(f2,k,q,u_ph)
  use comparams
  implicit none
  
  real*8 :: norm
  real*8, dimension(2) :: k,q
  complex*16, dimension(3) :: e1, e2, e3
  complex*16, dimension(3) :: eAq, eBq
  complex*16 :: f2, phase
  complex*16, dimension(6) :: u_ph
  
  e1=(/a1(1)+a2(1), a1(2)+a2(2), 0 /)
  norm=dot_product(e1,e1)
  e1=e1/sqrt(norm)
  
  e2=(/a1(1)-2.d0*a2(1), a1(2)-2.d0*a2(2), 0 /)
  norm=dot_product(e2,e2)
  e2=e2/sqrt(norm)
  
  e3=(/a2(1)-2.d0*a1(1), a2(2)-2.d0*a1(2), 0 /)
  norm=dot_product(e3,e3)
  e3=e3/sqrt(norm)
  
  eAq=u_ph(1:3)
  norm=dot_product(eAq,eAq)
  eAq=eAq/sqrt(norm)
  
  eBq=u_ph(4:6)
  norm=dot_product(eBq,eBq)
  eBq=eBq/sqrt(norm)
  
  f2=(0.d0,0.d0)
  
  f2=f2+exp(-i1*dot_product(k+q,(a1+a2)/3.d0))*dot_product(e1,eAq-eBq)
  f2=f2+exp(-i1*dot_product(k+q,(a1-2.d0*a2)/3.d0))*dot_product(e2,eAq-eBq*exp(i1*dot_product(q,-a2)))
  f2=f2+exp(-i1*dot_product(k+q,(a2-2.d0*a1)/3.d0))*dot_product(e3,eAq-eBq*exp(i1*dot_product(q,-a1)))
  
  call fnCalculatePhase(phase,k+q)
  f2=f2*phase
  return
end