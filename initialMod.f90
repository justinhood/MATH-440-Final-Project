module initialMod
contains
        subroutine initializeGrid(master_grid)
                implicit none
                double precision, allocatable, dimension(:,:), intent(inout) :: master_grid
                integer :: i, j
                allocate(master_grid(10,10))
                

                do i=1, size(master_grid, 1)
                        do j=1, size(master_grid, 2)
                                master_grid(i,j)=initialCondition(i,j)
                        enddo
                enddo
        end subroutine initializeGrid
       
        subroutine doStep(master_grid, new_grid, grid_row, grid_col)
                implicit none
                integer :: i,j, grid_row, grid_col
                double precision, dimension(grid_row, grid_col), intent(in) :: master_grid
                double precision, dimension(grid_row, grid_col), intent(inout) :: new_grid
                
                do i=1, grid_row
                        do j=1, grid_col
                                new_grid(i,j)=stepFunction(master_grid(i,j))
                        enddo
                enddo

        end subroutine doStep

        function initialCondition(x,y)
                implicit none
                integer, intent(in) :: x,y
                integer :: initialCondition

                initialCondition=x*y
        end function initialCondition

        function stepFunction(x)
                implicit none
                integer, intent(in):: x
                integer :: stepFunction
                stepFunction=x-1
        end function stepFunction


end module initialMod
