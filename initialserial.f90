program initialserial
  use initialMod
  implicit none

  integer :: i, j, n, grid_row, grid_col
  double precision, allocatable, dimension(:,:) :: master_grid, new_grid

  ! Initialize the grid 
  CALL initializeGrid(master_grid)

  ! Create temporary grid same size as master
  grid_row = size(master_grid,1)
  grid_col = size(master_grid,2)
  allocate(new_grid(grid_row,grid_col))

  ! Open up file
  open(unit = 20, file = "test.txt")

  do i = 1, grid_row
     write(20, *) (master_grid(i,j), j = 1, grid_col)
  end do
  write(20, *) ''
  
  do n = 1,10
  
     ! Do one step of the numerical method
     CALL doStep(master_grid, new_grid, grid_row, grid_col)


     ! Write the next time step grid to a file 
     do i = 1, grid_row
        write(20, *) (new_grid(i,j), j = 1, grid_col)
     end do
     
     write(20, *) ''
     
     master_grid = new_grid 

  end do

  ! Close file
  close(20)

     

end program initialserial

  
