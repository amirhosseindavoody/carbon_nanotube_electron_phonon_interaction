!*******************************************************************************
! This subroutines interprets the input arguments of the simulation
!*******************************************************************************
subroutine fnInput()
  use comparams
  implicit none
  
  integer i,count
  character*20 :: buffer
  
  ! set the simulation variables to default values
  n_ch=17
  m_ch=0
  nkg=501
  nr=200
  E_th=1.5d0
  Kcm_max=1.5d9
  i_sub=1
  
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
    elseif (buffer .eq. 'nr') then
      i=i+1
      call getarg(i,buffer)
      read(buffer,*) nr
    elseif (buffer .eq. 'E_th') then
      i=i+1
      call getarg(i,buffer)
      read(buffer,*) E_th      
    elseif (buffer .eq. 'Kcm_max') then
      i=i+1
      call getarg(i,buffer)
      read(buffer,*) Kcm_max
      Kcm_max=Kcm_max*1.d9
    elseif (buffer .eq. 'i_sub') then
      i=i+1
      call getarg(i,buffer)
      read(buffer,*) i_sub 
    else
        print *, "ERROR in input arguments!"
        pause
        stop
    end if
    i=i+1
  end do
  
  return
end