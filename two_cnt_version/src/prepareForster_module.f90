module prepareForster_module
	implicit none
	private
    public  :: calculateDOS,calculatePartitionFunction, saveDOS
    
contains
	!**************************************************************************************************************************
	! calculate the partition function for target or all exciton types of a given carbon nanotube
	!**************************************************************************************************************************
	
	subroutine calculatePartitionFunction(partition_function_type, currcnt, partitionFunction)
		use cnt_class, only: cnt
		use comparams, only : Temperature
		use physicalConstants, only : kb

		integer, intent(in) :: partition_function_type
		type(cnt), intent(in) :: currcnt
		real*8, intent(out) :: partitionFunction
		integer :: ix, iKcm
		real*8 :: min_energy, deltaE

		min_energy = minval(currcnt%Ex_t)
		deltaE = (-1.d0) * log(1.d-3) * kb * Temperature
        
		partitionFunction = 0.d0
        
        if (partition_function_type .eq. 0) then
			! partition function due to A-type excitons
			do ix = 1,currcnt%nX_a
				do iKcm = currcnt%iKcm_min_fine,currcnt%iKcm_max_fine
					if((currcnt%Ex_A1(ix,iKcm)-min_energy) .le. deltaE) then
						partitionFunction = partitionFunction + exp(-currcnt%Ex_A1(ix,iKcm)/kb/Temperature)    
					endif
					if((currcnt%Ex0_A2(ix,iKcm)-min_energy) .le. deltaE) then
						partitionFunction = partitionFunction + exp(-currcnt%Ex0_A2(ix,iKcm)/kb/Temperature)    
					endif
					if((currcnt%Ex1_A2(ix,iKcm)-min_energy) .le. deltaE) then
						partitionFunction = partitionFunction + exp(-currcnt%Ex1_A2(ix,iKcm)/kb/Temperature)    
					endif
				end do
			end do
			! partition function due to E-type excitons
			do ix = 1,currcnt%nX_e
				do iKcm = currcnt%iKcm_min_fine,currcnt%iKcm_max_fine
					if((currcnt%Ex0_Ep(ix,iKcm)-min_energy) .le. deltaE) then
						partitionFunction = partitionFunction + exp(-currcnt%Ex0_Ep(ix,iKcm)/kb/Temperature)    
					endif
					if((currcnt%Ex0_Em(ix,iKcm)-min_energy) .le. deltaE) then
						partitionFunction = partitionFunction + exp(-currcnt%Ex0_Em(ix,iKcm)/kb/Temperature)    
					endif
					if((currcnt%Ex1_Ep(ix,iKcm)-min_energy) .le. deltaE) then
						partitionFunction = partitionFunction + exp(-currcnt%Ex1_Ep(ix,iKcm)/kb/Temperature)    
					endif
					if((currcnt%Ex1_Em(ix,iKcm)-min_energy) .le. deltaE) then
						partitionFunction = partitionFunction + exp(-currcnt%Ex1_Em(ix,iKcm)/kb/Temperature)    
					endif
				end do
			end do

		elseif (partition_function_type .eq. 1) then
			do ix = 1,currcnt%nX_t
				do iKcm = currcnt%iKcm_min_fine,currcnt%iKcm_max_fine
					if((currcnt%Ex_t(ix,iKcm)-min_energy) .le. deltaE) then
						partitionFunction = partitionFunction + exp(-currcnt%Ex_t(ix,iKcm)/kb/Temperature)    
					endif
				end do
			end do
		else
			write(*,*) "Error in determining the type of partition function!"
		endif
        
	end subroutine calculatePartitionFunction

	!**************************************************************************************************************************
	! calculate the density of states at a given point
	!**************************************************************************************************************************
	
	subroutine calculateDOS(currcnt,iKcm,iX,dos)
		use cnt_class, only: cnt

		type(cnt), intent(in) :: currcnt
		integer, intent(in) :: iKcm,iX
		real*8, intent(out) :: dos
		
		if (iKcm .le. (currcnt%iKcm_min_fine+1)) then
			dos = 2.d0*currcnt%dkx/abs(-3.d0*currcnt%Ex_t(iX,iKcm)+4.d0*currcnt%Ex_t(iX,iKcm+1)-currcnt%Ex_t(iX,iKcm+2))
		else if(iKcm .ge. (currcnt%iKcm_max_fine-1)) then
			dos = 2.d0*currcnt%dkx/abs(3.d0*currcnt%Ex_t(iX,iKcm)-4.d0*currcnt%Ex_t(iX,iKcm-1)+currcnt%Ex_t(iX,iKcm-2))
		else if(iKcm == 0) then
			dos = 2.d0*currcnt%dkx/abs(3.d0*currcnt%Ex_t(iX,iKcm)-4.d0*currcnt%Ex_t(iX,iKcm-1)+currcnt%Ex_t(iX,iKcm-2))
		else
			dos = 12.d0*currcnt%dkx / abs(currcnt%Ex_t(iX,iKcm-2)-8.d0*currcnt%Ex_t(iX,iKcm-1)+8.d0*currcnt%Ex_t(iX,iKcm+1)-currcnt%Ex_t(iX,iKcm+2))
		end if
		return
	end subroutine calculateDOS

	!**************************************************************************************************************************
	! save total exciton density of states for a given cnt
	!**************************************************************************************************************************
	
	subroutine saveDOS(currcnt)
		use cnt_class, only: cnt

		type(cnt), intent(in) :: currcnt
		real*8 :: dos
		integer :: iKcm, iX
		character(len=100) :: filename

		write(filename,"( 'dos-cnt(', I2.2, ',', I2.2, ')-i_sub(', I1.1, ')-Ckappa(', F3.1, ').dat' )") currcnt%n_ch, currcnt%m_ch, currcnt%i_sub, currcnt%Ckappa

		!write currcnt dos
		open(unit=100,file=filename,status="unknown")
		do iKcm=currcnt%iKcm_min_fine,currcnt%iKcm_max_fine
			do iX=1,currcnt%nX_t
				call calculateDOS(currcnt,iKcm,iX,DOS)
				write(100,'(E16.8)', advance='no') dos
			enddo
			write(100,'(E16.8)')
		enddo
		close(100)

		return    
	end subroutine saveDOS
			
end module prepareForster_module