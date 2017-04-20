program par2
  use initialMod
  use mpi
  implicit none

  integer :: i, j, n, grid_row, grid_col
  integer :: ierror, my_rank, num_cores, num_steps, last_core
  integer :: div, rem
  integer :: source
  integer :: master = 0
  integer :: tag = 100
  integer :: fail_code
  integer, dimension(MPI_STATUS_SIZE) :: mpi_status
  double precision :: start_time, end_time, x_scale, y_scale
  double precision :: xmin, xmax, ymin, ymax, t_step, t_final
  double precision, allocatable, dimension(:,:) :: master_grid, recv_grid, concat_grid, send_grid
  
  CALL MPI_Init(ierror)
  start_time = MPI_Wtime()
  CALL MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierror)
  CALL MPI_Comm_SIZE(MPI_COMM_WORLD, num_cores, ierror)

  last_core = num_cores-1
  open(unit = 20, file = 'input.txt')
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
  if (my_rank == last_core) then
     allocate(master_grid(grid_row,grid_col))
     CALL initializeGrid(master_grid,grid_row,grid_col,x_scale,y_scale,xmin,ymin)
     ! Open up file
     open(unit = 21, file = "test.txt")
  end if

  CALL MPI_Barrier(MPI_COMM_WORLD, ierror)

  !allocate each matrix depending on the core 
  div = grid_col/num_cores
  rem = div + mod(grid_col, num_cores)
  if (my_rank == master) then
     allocate(recv_grid(grid_row, div),send_grid(grid_row, div), concat_grid(grid_row,div+1))
  else if (my_rank .NE. num_cores-1) then
     allocate(recv_grid(grid_row, div),send_grid(grid_row, div),concat_grid(grid_row,div+2))
  else
     allocate(recv_grid(grid_row, rem),send_grid(grid_row, rem),concat_grid(grid_row,rem+1))
  end if

  !start the looping over time steps
  do n = 1, num_steps
     !write to file the current timestep
     if (my_rank == last_core) then 
        do i = 1, grid_row
           write(21, *) (master_grid(i,j), j = 1, grid_col)
        end do
     end if
     write(21,*) '' 

     !give all the cores the info
     CALL MPI_Scatter(master_grid, grid_row*div, MPI_DOUBLE_PRECISION, recv_grid, &
          grid_row*div, MPI_DOUBLE_PRECISION, last_core, MPI_COMM_WORLD, ierror)

     if (my_rank == last_core) then
        recv_grid = master_grid(:,grid_col-rem+1:grid_col)
       ! do i = 1, grid_row
          ! print *, my_rank, (recv_grid(i,j), j = 1, rem)
       ! end do

     else
       ! do i = 1, grid_row
         !  print *, my_rank, (recv_grid(i,j), j = 1, div)
      !  end do

        
     end if 

     !communicate and update grid
     if (my_rank == master) then
        
        CALL MPI_Send(recv_grid(:,div), grid_row, MPI_DOUBLE_PRECISION, my_rank+1, tag*2*(my_rank+1), &
             MPI_COMM_WORLD, ierror)
        concat_grid(:,1:div) = recv_grid
        CALL MPI_Recv(concat_grid(:,div+1), grid_row, MPI_DOUBLE_PRECISION, my_rank+1, tag*3*(my_rank+1), &
             MPI_COMM_WORLD, mpi_status, ierror)

        ! Do one step of the numerical method
        CALL doStep(concat_grid, grid_row, div+1, t_step, x_scale, y_scale)

        send_grid = concat_grid(:,1:div)

      !  do i = 1,grid_row
       !    print *, my_rank, (send_grid(i,j), j = 1, div)
      !  end do
      !  print *, ''

        
     else if (my_rank .NE. num_cores-1) then
        
        CALL MPI_Send(recv_grid(:,1), grid_row, MPI_DOUBLE_PRECISION, my_rank-1, tag*3*my_rank, &
             MPI_COMM_WORLD, ierror)
        CALL MPI_Send(recv_grid(:,div), grid_row, MPI_DOUBLE_PRECISION, my_rank+1, tag*2*(my_rank+1), &
             MPI_COMM_WORLD, ierror)
        
        concat_grid(:,2:1+div) = recv_grid
        
        CALL MPI_Recv(concat_grid(:,1), grid_row, MPI_DOUBLE_PRECISION, my_rank-1, tag*2*my_rank, &
             MPI_COMM_WORLD, mpi_status, ierror)
        CALL MPI_Recv(concat_grid(:,div+2), grid_row, MPI_DOUBLE_PRECISION, my_rank+1, tag*3*(my_rank+1), &
             MPI_COMM_WORLD, mpi_status, ierror)

        ! Do one step of the numerical method
        CALL doStep(concat_grid, grid_row, div+2, t_step, x_scale, y_scale)
        
        send_grid = concat_grid(:,2:div+1)
      !  do i = 1,grid_row
      !     print *, my_rank, (send_grid(i,j), j = 1, div)
      !  end do
      !  print *, ''
     else

        CALL MPI_Send(recv_grid(:,1), grid_row, MPI_DOUBLE_PRECISION, my_rank-1, tag*3*(my_rank), &
             MPI_COMM_WORLD, ierror)
        
        concat_grid(:,2:rem+1) = recv_grid
        
        CALL MPI_Recv(concat_grid(:,1), grid_row, MPI_DOUBLE_PRECISION, my_rank-1, tag*2*(my_rank), &
             MPI_COMM_WORLD, mpi_status, ierror)

        ! Do one step of the numerical method
        CALL doStep(concat_grid, grid_row, rem+1, t_step, x_scale, y_scale)
        
        send_grid = concat_grid(:,2:rem+1)

       ! do i = 1,grid_row
      !     print *, my_rank, (send_grid(i,j), j = 1, rem)
      !  end do
      !  print *, ''

      end if

   ! last core gathers instead of master
     CALL MPI_Gather(send_grid, grid_row*div, MPI_DOUBLE_PRECISION, master_grid, & 
         grid_row*div, MPI_DOUBLE_PRECISION, num_cores-1, MPI_COMM_WORLD, ierror)
     
     if (my_rank == last_core) then 
         master_grid(:,grid_col-rem+1:grid_col) = send_grid
     end if

    ! if (my_rank == last_core) then 
    !    open(unit = 23, file = 'test2.txt')
        
      !  do i = 1, grid_row
      !     write(23, *) (master_grid(i,j), j = 1, grid_col)
      !  end do

    ! end if

  end do

  ! Close file
  close(21)
 ! close(23)

  CALL MPI_Finalize(ierror)  

end program par2

  
