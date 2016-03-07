!**********************************************************************************************************************
! This subroutines interprets the input arguments of the simulation
!**********************************************************************************************************************
subroutine fnGeomProp()
  use comparams
  implicit none
  integer :: i,j,k
  real*8 :: cosTh, sinTh
  real*8, dimension(2,2) :: Rot
  real*8, dimension(2) :: tmp_vec
  
  ! unit vectors and reciprocal lattice vectors************************************************************************
  a1=(/ sqrt(3.d0)/2.d0*a_l , +1.d0/2.d0*a_l /)
  a2=(/ sqrt(3.d0)/2.d0*a_l , -1.d0/2.d0*a_l /)
  b1=(/ 1.d0/sqrt(3.d0)*2.d0*pi/a_l , +1.d0*2.d0*pi/a_l /)
  b2=(/ 1.d0/sqrt(3.d0)*2.d0*pi/a_l , -1.d0*2.d0*pi/a_l /)
  
  aCC_vec=1.d0/3.d0*(a1+a2)
  
  ! calculate chirality and translational vectors of CNT unit cell.****************************************************
  ch_vec=dble(n_ch)*a1+dble(m_ch)*a2

  len_ch=a_l*sqrt(dble(n_ch)**2+dble(m_ch)**2+dble(n_ch)*dble(m_ch))
  radius=len_ch/2.d0/pi
  
  call gcd(dR,2*n_ch+m_ch,2*m_ch+n_ch)
  
  t1=+(2.d0*dble(m_ch)+dble(n_ch))/dble(dR)
  t2=-(2.d0*dble(n_ch)+dble(m_ch))/dble(dR)
  
  t_vec=t1*a1+t2*a2
  
  Nu=2*(n_ch**2+m_ch**2+n_ch*m_ch)/dR
  
  ! rotate basis vectors so that ch_vec is along x-axis
  cosTh=ch_vec(1)/norm2(ch_vec)
  sinTh=ch_vec(2)/norm2(ch_vec)
  Rot=reshape((/ cosTh, -sinTh , sinTh, cosTh /), (/2,2/))
  ch_vec=matmul(Rot,ch_vec)
  t_vec=matmul(Rot,t_vec)
  a1=matmul(Rot,a1)
  a2=matmul(Rot,a2)
  b1=matmul(Rot,b1)
  b2=matmul(Rot,b2)
  aCC_vec=matmul(Rot,aCC_vec)
  
  ! calculate reciprocal lattice of CNT.*******************************************************************************
  dk=norm2(b1)/(dble(nkg)-1.d0)
  
  K1=(-dble(t2)*b1+dble(t1)*b2)/dble(Nu)
  K2=(dble(m_ch)*b1-dble(n_ch)*b2)/dble(Nu)
  K2=K2/norm2(K2)
  
  ! calculate coordinates of atoms in the unrolled CNT unit cell.******************************************************
  allocate(posA(Nu,2))
  allocate(posB(Nu,2))
    
  k=0
  do i=0,t1+n_ch
    do j=t2,m_ch
      if ((dble(t2)/dble(t1)*i .le. j) .and. (dble(m_ch)/dble(n_ch)*i .ge. j) .and. (dble(t2)/dble(t1)*(i-n_ch) .gt. (j-m_ch)) .and. (dble(m_ch)/dble(n_ch)*(i-t1) .lt. (j-t2))) then
        k=k+1
        posA(k,1)=dble(i)*a1(1)+dble(j)*a2(1)
        posA(k,2)=dble(i)*a1(2)+dble(j)*a2(2)
        posB(k,1)=posA(k,1)+aCC_vec(1)
        posB(k,2)=posA(k,2)+aCC_vec(2)
        
        if (posA(k,1) .gt. ch_vec(1))   posA(k,1)=posA(k,1)-ch_vec(1);
        if (posA(k,1) .lt. 0)           posA(k,1)=posA(k,1)+ch_vec(1);
        if (posA(k,2) .gt. t_vec(2))    posA(k,2)=posA(k,2)-t_vec(2);
        if (posA(k,2) .lt. 0)           posA(k,2)=posA(k,2)+t_vec(2);
        
        if (posB(k,1) .gt. ch_vec(1))   posB(k,1)=posB(k,1)-ch_vec(1);
        if (posB(k,1) .lt. 0)           posB(k,1)=posB(k,1)+ch_vec(1);
        if (posB(k,2) .gt. t_vec(2))    posB(k,2)=posB(k,2)-t_vec(2);
        if (posB(k,2) .lt. 0)           posB(k,2)=posB(k,2)+t_vec(2);
        
      endif
    enddo
  enddo
    
  if (k .ne. Nu) stop "*** Error in calculating atom positions ***"
  
  ! calculate distances between atoms in a warped CNT unit cell.*******************************************************
  allocate(posAA(Nu,2))
  allocate(posAB(Nu,2))
  allocate(posBA(Nu,2))
  allocate(posBB(Nu,2))
  
  do i=1,Nu
    posAA(i,:)=posA(i,:)-posA(1,:)
    posAB(i,:)=posA(i,:)-posB(1,:)
    posBA(i,:)=posB(i,:)-posA(1,:)
    posBB(i,:)=posB(i,:)-posB(1,:)
    if (posAA(i,1) .gt. ch_vec(1)/2.d0) posAA(i,1)=posAA(i,1)-ch_vec(1)
    if (posAB(i,1) .gt. ch_vec(1)/2.d0) posAB(i,1)=posAB(i,1)-ch_vec(1)
    if (posBA(i,1) .gt. ch_vec(1)/2.d0) posBA(i,1)=posBA(i,1)-ch_vec(1)
    if (posBB(i,1) .gt. ch_vec(1)/2.d0) posBB(i,1)=posBB(i,1)-ch_vec(1)
    
  end do
   
  return
end