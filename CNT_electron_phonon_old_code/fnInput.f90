!*******************************************************************************
! This subroutines interprets the input arguments of the simulation
!*******************************************************************************
subroutine fnInput()
  use comparams
  implicit none
  
  integer i,count
  character*20 :: buffer
  
  ! set the simulation variables to default values
  n_ch=10
  m_ch=0
  nkg=501
  
  ! update the simulation variables according to input variables
  count=nargs()
  i=1
  
  do while (i .le. count-1)
    call getarg(i,buffer)
    if (buffer .eq. 'ch') then
      i=i+1
      call getarg(i,buffer)
      read(buffer,*) n_ch
      i=i+1
      call getarg(i,buffer)
      read(buffer,*) m_ch
    elseif (buffer .eq. 'nkg') then
      i=i+1
      call getarg(i,buffer)
      read(buffer,*) nkg
    else
        print *, "ERROR in input arguments!"
        pause
        stop
    end if
    i=i+1
  end do
  
  return
end