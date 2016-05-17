!**********************************************************************************************************************
! This subroutine reads the dispersion of excitons.
!**********************************************************************************************************************
subroutine fnFindMode()
  use comparams
  implicit none
  
  integer :: iKcm1, iKcm2, iX1, iX2, i_branch
  real*8 :: tmpr1, tmpr2, tmpr3, tmpr4
  real*8 :: ratio1, ratio2
  real*8 :: ScatRate, TotScatRate
  real*8, dimension (:,:), allocatable:: omega_shifted
  
  
  do iKcm1=1,iKcm_max
    Ex0_A2(:,-iKcm1)=Ex0_A2(:,iKcm1)
  enddo
  
  iX1=1
  !iKcm1=iKcm_min+39-1
  do iKcm1=1,iKcm_max
    print *, 'Calculating Ex-Ph scattering rates -->', iKcm1
    TotScatRate=0
    do iX2=1,nX
      do iKcm2=iKcm_min,iKcm_max-1
        do i_branch=1,3
          !tmpr1=(Ex0_A2(iX2,iKcm2)-Ex0_A2(iX1,iKcm1))-omega(-2*(iKcm1-iKcm2),i_branch)
          !tmpr2=(Ex0_A2(iX2,iKcm2+1)-Ex0_A2(iX1,iKcm1))-omega(-2*(iKcm1-iKcm2-1),i_branch)
          tmpr1=(Ex0_A2(iX2,iKcm2)-Ex0_A2(iX1,iKcm1))
          tmpr2=(Ex0_A2(iX2,iKcm2+1)-Ex0_A2(iX1,iKcm1))
          if ((tmpr1*tmpr2) .le. 0.d0) then
            ratio1=abs(tmpr1)/(abs(tmpr1)+abs(tmpr2))
            ratio2=abs(tmpr2)/(abs(tmpr1)+abs(tmpr2))
            
            call fnExPhScatRate(ScatRate,iKcm1,iX1,iKcm2,iX2,i_branch,ratio1,ratio2)
            TotScatRate=TotScatRate+ScatRate
            
          endif
        enddo
      enddo
    enddo
    write(fh9,10) dble(iKcm1)/abs(dble(iKcm1))*Ex0_A2(iX1,iKcm1), TotScatRate
  enddo
  
  10 FORMAT (E16.8,E16.8)
  
  return
end