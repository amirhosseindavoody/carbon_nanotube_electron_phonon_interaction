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
  
  !dirname='test'
  !istat=chdir(dirname)
  !if (istat .ne. 0) then 
  !  print *, 'Directory did not change!!!'
  !  pause
  !  stop
  !endif
  
  fh1=1
  open(unit=fh1,file='phonon_dispersion.dat',status="unknown")
  
  fh2=2
  open(unit=fh2,file='electron_conduction.dat',status="unknown")
  
  fh3=3
  open(unit=fh3,file='electron_valence.dat',status="unknown")
  
  fh4=4
  open(unit=fh4,file='ScatRateN.dat',status="unknown")
  
  fh5=5
  open(unit=fh5,file='ScatRateP.dat',status="unknown")

  fh6=6
  open(unit=fh6,file='ScatRate.dat',status="unknown")
  !
  !fh7=7
  !open(unit=fh7,file='Psi0_A2.dat',status="old", action='read')
  !
  !fh8=8
  !open(unit=fh8,file='Psi1_A2.dat',status="old", action='read')
  !
  !fh9=9
  !open(unit=fh9,file='TotScatRate.dat',status="unknown")

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
  !close(fh7)
  !close(fh8)
  !close(fh9)
  
  return
end