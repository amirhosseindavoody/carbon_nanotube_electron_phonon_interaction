!**************************************************************************************************************************
! this module calculates and saves the occupation of excitonic states in a certain CNT
!**************************************************************************************************************************

module occupation_mod
	implicit none
	private
	public  :: calculate_occupation_table
	
contains
	!**************************************************************************************************************************
	! calculate the partition function for a given carbon nanotube
	!**************************************************************************************************************************
	
	subroutine calculate_occupation_table(currcnt)
		use comparams, only : Temperature
		use cnt_class, only: cnt
		use physicalConstants, only : kb, eV

		type(cnt), intent(in) :: currcnt
		integer :: table_size, iX, iKcm, counter
		integer :: iKcm_max_fine, iKcm_min_fine
		real*8 :: partitionFunction
		real*8 , dimension(:,:), allocatable :: occupation_table
		real*8 :: occupation_Atype, occupation_Etype
		character(len=100) :: filename
		
		iKcm_min_fine = currcnt%dk_dkx_ratio * currcnt%iKcm_min
		iKcm_max_fine = currcnt%dk_dkx_ratio * currcnt%iKcm_max

		table_size = size(currcnt%Ex_t)
		allocate(occupation_table(table_size,3))

		occupation_table = 0.d0*occupation_table

		partitionFunction = 0.d0

		occupation_Atype = 0.d0
		do ix = 1,currcnt%nX_a
			do iKcm = iKcm_min_fine,iKcm_max_fine
				occupation_Atype = occupation_Atype + exp(-currcnt%Ex0_A2(ix,iKcm)/kb/Temperature)
			end do
		end do

		occupation_Etype = 0.d0
		do ix = 1,currcnt%nX_e
			do iKcm = iKcm_min_fine,iKcm_max_fine
				occupation_Etype = occupation_Etype + exp(-currcnt%Ex0_Ep(ix,iKcm)/kb/Temperature)
				occupation_Etype = occupation_Etype + exp(-currcnt%Ex0_Em(ix,iKcm)/kb/Temperature)
			end do
		end do
		write(*,*) "minimum energy of A2  excitons: ", minval(currcnt%Ex0_A2)/eV
		write(*,*) "minimum energy of E+ excitons: ", minval(currcnt%Ex0_Ep)/eV
		write(*,*) "minimum energy of E- excitons: ", minval(currcnt%Ex0_Em)/eV
		write(*,*) "A-type occupation percentage is: ", occupation_Atype/(occupation_Atype+occupation_Etype)
		write(*,*) "E-type occupation percentage is: ", occupation_Etype/(occupation_Atype+occupation_Etype)

! 		call exit()


		counter = 0
		do ix = 1,currcnt%nX_t
			do iKcm = iKcm_min_fine,iKcm_max_fine
				counter = counter + 1
				occupation_table(counter,1) = currcnt%Ex_t(ix,iKcm)
				occupation_table(counter,2) = exp(-currcnt%Ex_t(ix,iKcm)/kb/Temperature)
				partitionFunction = partitionFunction + occupation_table(counter,2)
			end do
		end do

		occupation_table(:,3) = partitionFunction

		write(filename,"( 'occupation-cnt(', I2.2, ',', I2.2, ')-i_sub(', I1.1, ')-Ckappa(', F3.1, ').dat' )") currcnt%n_ch, currcnt%m_ch, currcnt%i_sub, currcnt%Ckappa

		!write occupation_table to file
		open(unit=100,file=filename,status="unknown")
		do counter=1,table_size
			write(100,'(E16.8, E16.8)') occupation_table(counter,1)/eV, occupation_table(counter,2)/partitionFunction
		enddo
		close(100)

		return
		
	end subroutine calculate_occupation_table

	
			
end module occupation_mod