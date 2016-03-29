!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Declaration of input parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module write_log_mod
	implicit none
	private
	public :: write_log
	character(len=1000), public :: log_input

contains
	!*******************************************************************************
	! This subroutines opens the log file add new log and closes the file
	!*******************************************************************************

	subroutine write_log(message)
		character(len=*) :: message
		logical :: flgexist
		integer :: logFile = 10

		inquire(file="log.dat",exist=flgexist)
		if (flgexist) then
			open(logFile, file="log.dat", status="old", position="append", action="write")
		else
			open(logFile, file="log.dat", status="new", action="write")
		end if
		write(logFile, *) message
		close(logFile)

		return
	end subroutine write_log

end module write_log_mod
