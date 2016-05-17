!*******************************************************************************
! This subroutines opens all the files that are used in the simulation
!*******************************************************************************
subroutine fnOpenFiles()
  use comparams
  use ifport
  implicit none
  
  character*100 :: dirname
  integer(4) :: istat
  logical(4) :: result
  

  write(dirname,"('CNT_',I2.2,'_',I2.2,'_',I4.4,'_',I4.4,'_',F3.1,'_',I1.1)") n_ch, m_ch, nkg, nr, E_th, i_sub
  
  istat=chdir(dirname)
  if (istat .ne. 0) then 
    print *, 'Directory did not change!!!'
    pause
    stop
  endif
  
  fh1=1
  open(unit=fh1,file='phonon_dispersion.dat',status="unknown")
  
  fh2=2
  open(unit=fh2,file='miscellaneous.dat',status="old", action='read')
  
  fh3=3
  open(unit=fh3,file='Ex_A1.dat',status="old", action='read')
  
  fh4=4
  open(unit=fh4,file='Ex0_A2.dat',status="old", action='read')
  
  fh5=5
  open(unit=fh5,file='Ex1_A2.dat',status="old", action='read')
  
  fh6=6
  open(unit=fh6,file='Psi_A1.dat',status="old", action='read')
  
  fh7=7
  open(unit=fh7,file='Psi0_A2.dat',status="old", action='read')
  
  fh8=8
  open(unit=fh8,file='Psi1_A2.dat',status="old", action='read')
  
  fh9=9
  open(unit=fh9,file='TotScatRate.dat',status="unknown")

  return
end
    
!*******************************************************************************
! This subroutines closes all the files that are used in the simulation
!*******************************************************************************
subroutine fnCloseFiles()
  use comparams
  implicit none
  
  close(fh1)
  close(fh2)
  close(fh3)
  close(fh4)
  close(fh5)
  close(fh6)
  close(fh7)
  close(fh8)
  close(fh9)
  
  return
end