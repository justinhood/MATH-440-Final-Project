program par3
  use Mod2
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
  double precision :: start_time, end_time, x_scale, y_scale, alpha, beta
  double precision :: xmin, xmax, ymin, ymax, t_step, t_final
  double precision, allocatable, dimension(:,:) :: uold, u, uold_con, u_con, unew, old_grid, grid, new_grid, temp
  
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
  alpha = (0.5d0*t_step)/(x_scale**2)
  beta = (0.5d0*t_step)/(y_scale**2)

  ! Initialize the grid
  if (my_rank == last_core) then
     allocate(old_grid(grid_row,grid_col),grid(grid_row,grid_col),new_grid(grid_row,grid_col))
     CALL initializeGrid(old_grid,grid_row,grid_col,x_scale,y_scale,xmin,ymin)
     grid = old_grid
     new_grid=old_grid
     ! Open up file
     open(unit = 21, file = "implicit.txt")
  end if

  CALL MPI_Barrier(MPI_COMM_WORLD, ierror)

  !allocate each matrix depending on the core 
  div = grid_col/num_cores
  rem = div + mod(grid_col, num_cores)
  if (my_rank == master) then
     allocate(uold(grid_row, div),u(grid_row, div), uold_con(grid_row,div+1), u_con(grid_row,div+1), &
          temp(grid_row,div+1), unew(grid_row,div))
  else if (my_rank .NE. num_cores-1) then
     allocate(uold(grid_row, div),u(grid_row, div), uold_con(grid_row,div+2), u_con(grid_row,div+2), &
          temp(grid_row,div+2), unew(grid_row,div))
  else
     allocate(uold(grid_row, rem),u(grid_row, rem), uold_con(grid_row,rem+1), u_con(grid_row,rem+1), &
          temp(grid_row,rem+1), unew(grid_row,rem))
  end if

  !start the looping over time steps
  do n = 1, num_steps
     !write to file the current timestep
     if (my_rank == last_core) then 
        do i = 1, grid_row
           write(21, *) (grid(i,j), j = 1, grid_col)
        end do
     end if
     write(21,*) '' 

     !give all the cores the info
     CALL MPI_Scatter(old_grid, grid_row*div, MPI_DOUBLE_PRECISION, uold, &
          grid_row*div, MPI_DOUBLE_PRECISION, last_core, MPI_COMM_WORLD, ierror)
     CALL MPI_Scatter(grid, grid_row*div, MPI_DOUBLE_PRECISION, u, &
          grid_row*div, MPI_DOUBLE_PRECISION, last_core, MPI_COMM_WORLD, ierror)

     if (my_rank == last_core) then
        uold = old_grid(:,grid_col-rem+1:grid_col)
        u = grid(:,grid_col-rem+1:grid_col)
     else
    
     end if 

     !communicate and update grid
     if (my_rank == master) then
        
        ! send stuff for uold 
        CALL MPI_Send(uold(:,div), grid_row, MPI_DOUBLE_PRECISION, my_rank+1, tag*2*(my_rank+1), &
             MPI_COMM_WORLD, ierror)
        uold_con(:,1:div) = uold
        CALL MPI_Recv(uold_con(:,div+1), grid_row, MPI_DOUBLE_PRECISION, my_rank+1, tag*3*(my_rank+1), &
             MPI_COMM_WORLD, mpi_status, ierror)

        ! send stuff for u 
        CALL MPI_Send(u(:,div), grid_row, MPI_DOUBLE_PRECISION, my_rank+1, tag*4*(my_rank+1), &
             MPI_COMM_WORLD, ierror)
        u_con(:,1:div) = u
        CALL MPI_Recv(u_con(:,div+1), grid_row, MPI_DOUBLE_PRECISION, my_rank+1, tag*5*(my_rank+1), &
             MPI_COMM_WORLD, mpi_status, ierror)

        ! Do one step of the numerical method - unew is calculated and returned
        CALL doImplicitStep(uold_con, u_con, unew, x_scale, y_scale, t_step, alpha, beta, &
             grid_row, div+1, my_rank, temp)
        !do i = 1, grid_row
        !     write(*, *) (temp(i,j), j = 1, div+1)
        !end do
        !write(*,*) ''
        unew = temp(:,1:div)
       


     else if (my_rank .NE. num_cores-1) then

        ! send uold 
        CALL MPI_Send(uold(:,1), grid_row, MPI_DOUBLE_PRECISION, my_rank-1, tag*3*my_rank, &
             MPI_COMM_WORLD, ierror)
        CALL MPI_Send(uold(:,div), grid_row, MPI_DOUBLE_PRECISION, my_rank+1, tag*2*(my_rank+1), &
             MPI_COMM_WORLD, ierror)

        ! send u 
        CALL MPI_Send(u(:,1), grid_row, MPI_DOUBLE_PRECISION, my_rank-1, tag*5*my_rank, &
             MPI_COMM_WORLD, ierror)
        CALL MPI_Send(u(:,div), grid_row, MPI_DOUBLE_PRECISION, my_rank+1, tag*4*(my_rank+1), &
             MPI_COMM_WORLD, ierror)
        
        uold_con(:,2:1+div) = uold
        u_con(:,2:1+div) = u

        ! recv uold 
        CALL MPI_Recv(uold_con(:,1), grid_row, MPI_DOUBLE_PRECISION, my_rank-1, tag*2*my_rank, &
             MPI_COMM_WORLD, mpi_status, ierror)
        CALL MPI_Recv(uold_con(:,div+2), grid_row, MPI_DOUBLE_PRECISION, my_rank+1, tag*3*(my_rank+1), &
             MPI_COMM_WORLD, mpi_status, ierror)

        !recv u
        CALL MPI_Recv(u_con(:,1), grid_row, MPI_DOUBLE_PRECISION, my_rank-1, tag*4*my_rank, &
             MPI_COMM_WORLD, mpi_status, ierror)
        CALL MPI_Recv(u_con(:,div+2), grid_row, MPI_DOUBLE_PRECISION, my_rank+1, tag*5*(my_rank+1), &
             MPI_COMM_WORLD, mpi_status, ierror)
        

        ! Do one step of the numerical method
        CALL doImplicitStep(uold_con, u_con, unew, x_scale, y_scale, t_step, alpha, beta, &
             grid_row, div+2, my_rank, temp)
        !do i = 1, grid_row
        !     write(*, *) (temp(i,j), j = 1, div+2)
        !end do
        !write(*,*) ''
        unew = temp(:,2:1+div)
        
     else

        CALL MPI_Send(uold(:,1), grid_row, MPI_DOUBLE_PRECISION, my_rank-1, tag*3*(my_rank), &
             MPI_COMM_WORLD, ierror)
        CALL MPI_Send(u(:,1), grid_row, MPI_DOUBLE_PRECISION, my_rank-1, tag*5*(my_rank), &
             MPI_COMM_WORLD, ierror)
        
        uold_con(:,2:rem+1) = uold
        u_con(:,2:rem+1) = u
        
        CALL MPI_Recv(uold_con(:,1), grid_row, MPI_DOUBLE_PRECISION, my_rank-1, tag*2*(my_rank), &
             MPI_COMM_WORLD, mpi_status, ierror)
        CALL MPI_Recv(u_con(:,1), grid_row, MPI_DOUBLE_PRECISION, my_rank-1, tag*4*(my_rank), &
             MPI_COMM_WORLD, mpi_status, ierror)
        
        ! Do one step of the numerical method
        CALL doImplicitStep(uold_con, u_con, unew, x_scale, y_scale, t_step, alpha, beta, &
             grid_row, rem+1, my_rank, temp)
        !do i = 1, grid_row
        !     write(*, *) (temp(i,j), j = 1, rem+1)
        !end do
        !write(*,*) '***************************************************************'
        unew = temp(:,2:rem+1)
        
     end if
     
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
     
     ! last core gathers instead of master
     CALL MPI_Gather(unew, grid_row*div, MPI_DOUBLE_PRECISION, new_grid, & 
         grid_row*div, MPI_DOUBLE_PRECISION, num_cores-1, MPI_COMM_WORLD, ierror)
     
     if (my_rank == last_core) then 
         new_grid(:,grid_col-rem+1:grid_col) = unew
         old_grid = grid  
         grid = new_grid

     end if
      
       end do

  ! Close file
  close(21)
 ! close(23)
  end_time=MPI_Wtime()
  open(unit=22, file='timerimp.txt')
  if(my_rank==master) then
        write(22,*) end_time-start_time
  endif
  CALL MPI_Finalize(ierror)  

end program par3

  
