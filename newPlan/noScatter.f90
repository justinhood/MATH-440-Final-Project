program par2
  use Mod2
  use mpi
  implicit none

  integer :: i, j, n, grid_row, grid_col
  integer :: ierror, my_rank, num_cores, num_steps, last_core
  integer :: div, rem, first, last
  integer :: source
  integer :: master = 0
  integer :: tag = 100
  integer :: fail_code
  integer, allocatable, dimension(:) :: index_arr
  integer, dimension(MPI_STATUS_SIZE) :: mpi_status
  double precision :: start_time, end_time, x_scale, y_scale
  double precision :: xmin, xmax, ymin, ymax, t_step, t_final
  double precision, allocatable, dimension(:,:) :: master_grid,  concat_grid
  
  CALL MPI_Init(ierror)
  start_time = MPI_Wtime()
  CALL MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierror)
  CALL MPI_Comm_SIZE(MPI_COMM_WORLD, num_cores, ierror)

  last_core = num_cores-1
  open(unit = 20, file = 'newInput.txt')
  read(20,*) grid_row
  read(20,*) grid_col
  read(20,*) xmin
  read(20,*) ymin 
  read(20,*) xmax
  read(20,*) ymax
  read(20,*) t_final
  read(20, *) num_steps
  close(20)
  t_step = t_final/dble(num_steps)
  x_scale = (xmax-xmin)/(grid_row-1)
  y_scale = (ymax-ymin)/(grid_col-1)

  ! Initialize the grid
  if (my_rank == master) then
     allocate(master_grid(grid_row,grid_col))
     CALL initializeGrid(master_grid,grid_row,grid_col,x_scale,y_scale,xmin,ymin)
     ! Open up file
     open(unit = 21, file = "test.txt")
  end if

  CALL MPI_Barrier(MPI_COMM_WORLD, ierror)

  !allocate each matrix depending on the core 
  div = grid_col/num_cores
  rem = mod(grid_col, num_cores)
 
  !Start Allocation
  if(rem .ne. 0) then
	  !with remainder
	  if(my_rank .le. (rem-1)) then
		  if(my_rank == master) then
			  allocate(concat_grid(grid_row, div+2)) !1 from wrap and one from overlap
		  else
			  allocate(concat_grid(grid_row, div+3)) !1 from wrap and 2 from overlap
		  endif
	  else if(my_rank .ne. last_core) then
		  allocate(concat_grid(grid_row, div+2)) !0 from wrap and 2 from overlap
	  else
		  allocate(concat_grid(grid_row, div+1)) !0 from wrap and 1 from overlap
	  endif
  else
	  !without remainder
	  if(my_rank == master) then
		  allocate(concat_grid(grid_row, div+1)) !0 from wrap and 1 from overlap
	  else if(my_rank .ne. last_core) then
		  allocate(concat_grid(grid_row, div+2)) !0 from wrap and 2 from overlap
	  else
		  allocate(concat_grid(grid_row, div+1)) !0 from wrap and 1 from overlap
	  endif
  endif

  !Compute the indexes for sending
  allocate(index_arr(num_cores*2))
  first=1
  last=1
  if(my_rank == master) then
	  if(rem .ne. 0) then
		  last=div+2
		  index_arr(1)=first
		  index_arr(2)=last
		  first=last-1
		  if((rem-1) .ne. 0) then
			  do i=1, rem-1
				  last=first+div+2
				  index_arr(2*i+1)=first
				  index_arr(2*i+2)=last
				  first=last-1
			  enddo
		  endif
		  if((rem-1) .ne. (last_core-1)) then
			  do i=rem, last_core-1
				  last=first+div+1
				  index_arr(2*i+1)=first
				  index_arr(2*i+2)=last
				  first=last-1
			  enddo
		  endif
		  index_arr(2*last_core+1)=first
		  index_arr(2*last_core+2)=grid_col
	  else
		  last=div+1
		  index_arr(1)=first
		  index_arr(2)=last
		  first=last-1
		  do i=1, last_core-1
			  last=first+div+1
			  index_arr(2*i+1)=first
			  index_arr(2*i+2)=last
			  first=last-1
		  enddo
		  index_arr(2*last_core+1)=first
		  index_arr(2*last_core+2)=grid_col
	  endif
  endif
  call MPI_Bcast(index_arr, num_cores*2, MPI_INTEGER, master, MPI_COMM_WORLD, ierror)
  if (my_rank == master) then 
        do i = 1, grid_row
           write(21, *) (master_grid(i,j), j = 1, grid_col)
        end do
  end if
  write(21,*) ''



  !start the looping over time steps
  do n = 1, num_steps
	  if(my_rank==master) then
		  concat_grid=master_grid(:, index_arr(1):index_arr(2))
		  do i=1, last_core
			  CALL MPI_SEND(master_grid(:,index_arr(2*i+1):index_arr(2*i+2)), grid_row*(index_arr(2*i+2)-index_arr(2*i+1)+1),&
                                  MPI_DOUBLE_PRECISION, i, tag*i, MPI_COMM_WORLD, ierror)
		  enddo
	  else
		  CALL MPI_RECV(concat_grid(:,:), grid_row*(index_arr(2*my_rank+2)-index_arr(2*my_rank+1)+1), MPI_DOUBLE_PRECISION, &
			  master, tag*my_rank, MPI_COMM_WORLD, mpi_status, ierror)
	  endif

	  call doStep(concat_grid, grid_row, size(concat_grid,2), t_step, x_scale, y_scale)
	  
          if(my_rank .ne. master) then
                  if(my_rank .ne. last_core) then
                          call MPI_Send(concat_grid(:,2:(size(concat_grid,2)-1)), grid_row*(size(concat_grid,2)-2),&
                                 MPI_DOUBLE_PRECISION, master, tag*my_rank*500, MPI_COMM_WORLD, ierror)
                  else
                          call MPI_Send(concat_grid(:, 2:(size(concat_grid,2))), grid_row*(size(concat_grid,2)-1),&
                                 MPI_DOUBLE_PRECISION, master, tag*my_rank*500, MPI_COMM_WORLD, ierror)
                  endif
          else
                  master_grid(:,index_arr(1):(index_arr(2)-1))=concat_grid(:,1:(size(concat_grid,2)-1))
                  if(last_core .gt. 1) then
                        do i=1, last_core-1
                                call MPI_Recv(master_grid(:,(index_arr(i*2+1)+1):(index_arr(i*2+2)-1)), grid_row*(index_arr(2*i+2) &
                                        -index_arr(2*i+1)-1), MPI_DOUBLE_PRECISION, i, tag*i*500, mPI_COMM_WORLD, &
                                        mpi_status, ierror)
                        enddo
                  endif
                  call MPI_recv(master_grid(:, (index_arr(last_core*i+1)+1):grid_col), grid_row*(index_arr(2*last_core+2)&
                          -index_arr(2*last_core+1)), MPI_DOUBLE_PRECISION, last_core, tag*last_core*500, MPI_COMM_WORLD,&
                          mpi_status, ierror)
          endif
          if(my_rank==master) then
                  do i=1, grid_row
                          write(21,*) (master_grid(i,j), j=1, grid_col)
                  enddo
          endif
          write(21,*) ''
  end do

  ! Close file
  close(21)
 ! close(23)
  end_time=MPI_Wtime()
  open(unit=22, file='timer.txt')
  if(my_rank==master) then
        write(22,*) end_time-start_time
  endif
  CALL MPI_Finalize(ierror)  

end program par2

  
