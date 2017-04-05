program initialserial
  use initialMod
  implicit none

  integer :: i, j, n
  double precision, allocatable, dimension(:,:) :: master_grid, new_grid


  CALL initializeGrid(master_grid)
  
  open(unit = 20, file = "test.txt")
  
  do i = 1,size(master_grid,1)
     write(20, *) (master_grid(i,j), j = 1,size(master_grid,2))
  end do

  close(20)

     
  

end program initialserial

  