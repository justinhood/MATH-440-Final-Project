program initialserial
  use initialMod
  implicit none

  integer :: i, j, n, grid_row, grid_col
  double precision, allocatable, dimension(:,:) :: master_grid, new_grid

  ! Initialize the grid 
  CALL initializeGrid(master_grid)

  ! Create temporary grid same size as master
  grid_row = size(master_grid,1)
  grid_col = size(master_grid,2
  allocate(new_grid(grid_row,grid_col))

  ! Do one step of the numerical method 
  CALL doStep(master_grid, new_grid, grid_row, grid_col)

  
  ! Write the grid to a file 
  open(unit = 20, file = "test.txt")
  do i = 1,grid_row
     write(20, *) (new_grid(i,j), j = 1,grid_col)
  end do
  close(20)

     

end program initialserial

  
